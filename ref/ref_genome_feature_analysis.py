import gzip
import logging
import shutil
from pathlib import Path
from typing import NamedTuple, Optional

import pysam

from contig_classification import ContigCategorizer
from contig_name_translation import ContigNameTranslator
from ref_util import assert_file_exists, get_nucleotides_from_fasta, get_nucleotides_from_string, \
    STANDARD_NUCLEOTIDES, SOFTMASKED_NUCLEOTIDES, UNKNOWN_NUCLEOTIDES


class ReferenceGenomeFeatureAnalysis(NamedTuple):
    has_unplaced_contigs: bool
    has_unlocalized_contigs: bool
    has_alts: bool
    has_decoys: bool
    has_patches: bool
    has_ebv: Optional[bool]
    has_rcrs: Optional[bool]
    uses_canonical_chrom_names: bool
    has_only_hardmasked_nucleotides_at_y_par1: Optional[bool]
    has_semi_ambiguous_iub_codes: bool
    has_softmasked_nucleotides: bool
    alts_are_padded: Optional[bool]


class ReferenceGenomeFeatureAnalyzer(object):
    Y_PAR1_TEST_REGION = (20000, 2640000)  # this region lies within the Y PAR1 for both GRCh37 and GRCh38

    @classmethod
    def do_analysis(
            cls, ref_genome_path: Path, rcrs_path: Optional[Path], contig_name_translator: ContigNameTranslator,
    ) -> ReferenceGenomeFeatureAnalysis:
        # Assert that files exist where they should
        assert_file_exists(ref_genome_path)
        if rcrs_path is not None:
            assert_file_exists(rcrs_path)
        with pysam.Fastafile(ref_genome_path) as genome_f:
            contig_names = list(genome_f.references)

        contig_categorizer = ContigCategorizer(contig_name_translator)
        categorized_contig_names = contig_categorizer.get_categorized_contig_names(contig_names)

        logging.debug(categorized_contig_names)

        if len(categorized_contig_names.autosomes) != 22:
            warn_msg = f"Did not find exactly 22 autosome contigs: autosomes={categorized_contig_names.autosomes}"
            logging.warning(warn_msg)
        if len(categorized_contig_names.x_contigs) != 1:
            warn_msg = f"Did not find exactly one X contig: x={categorized_contig_names.x_contigs}"
            logging.warning(warn_msg)

        nucleotides = get_nucleotides_from_fasta(ref_genome_path)
        logging.info(f"nucleotides: {sorted(nucleotides)}")

        has_unplaced_contigs = bool(categorized_contig_names.unplaced_contigs)
        has_unlocalized_contigs = bool(categorized_contig_names.unlocalized_contigs)
        has_alts = bool(categorized_contig_names.alt_contigs)
        has_decoys = bool(categorized_contig_names.decoys)
        has_patches = (
            bool(categorized_contig_names.fix_patch_contigs)
            or bool(categorized_contig_names.novel_patch_contigs)
        )
        has_ebv: Optional[bool]
        if len(categorized_contig_names.ebv_contigs) == 1:
            has_ebv = True
        elif len(categorized_contig_names.ebv_contigs) == 0:
            has_ebv = False
        else:
            logging.warning(f"Found more than one EBV contig: {categorized_contig_names.ebv_contigs}")
            has_ebv = None
        has_rcrs: Optional[bool]
        if rcrs_path is not None and len(categorized_contig_names.mitochondrial_contigs) == 1:
            has_rcrs = cls._mitochondrial_sequence_is_rcrs(
                ref_genome_path, rcrs_path, categorized_contig_names.mitochondrial_contigs[0],
            )
        elif rcrs_path is None:
            warn_msg = (
                f"rCRS argument not provided, so skipping comparison of mitochondrial sequence."
            )
            logging.warning(warn_msg)
            has_rcrs = None
        else:
            warn_msg = (
                f"Did not find exactly one mitochondrial contig: "
                f"mitochondrial={categorized_contig_names.mitochondrial_contigs}"
            )
            logging.warning(warn_msg)
            has_rcrs = None
        uses_canonical_chrom_names = all(
            contig_name_translator.is_canonical(contig_name)
            for contig_name in categorized_contig_names.get_contig_names()
        )
        has_only_hardmasked_nucleotides_at_y_par1: Optional[bool]
        if len(categorized_contig_names.y_contigs) == 1:
            y_test_nucleotides = get_nucleotides_from_string(
                cls._get_y_test_sequence(categorized_contig_names.y_contigs[0], ref_genome_path)
            )
            logging.info(f"nucleotides at y par1 test region: {sorted(y_test_nucleotides)}")
            has_only_hardmasked_nucleotides_at_y_par1 = not bool(y_test_nucleotides.difference(UNKNOWN_NUCLEOTIDES))
        else:
            warn_msg = (
                f"Did not find exactly one Y contig: "
                f"y={categorized_contig_names.y_contigs}"
            )
            logging.warning(warn_msg)
            has_only_hardmasked_nucleotides_at_y_par1 = None
        has_semi_ambiguous_iub_codes = bool(
            nucleotides.difference(STANDARD_NUCLEOTIDES).difference(SOFTMASKED_NUCLEOTIDES)
        )
        has_softmasked_nucleotides = bool(nucleotides.intersection(SOFTMASKED_NUCLEOTIDES))
        alts_are_padded: Optional[bool]
        if categorized_contig_names.alt_contigs:
            alt_is_definitely_padded_list = [
                cls._is_definitely_padded_with_n(contig, ref_genome_path)
                for contig in categorized_contig_names.alt_contigs
            ]
            if all(is_padded for is_padded in alt_is_definitely_padded_list):
                alts_are_padded = True
            elif all(not is_padded for is_padded in alt_is_definitely_padded_list):
                alts_are_padded = False
            else:
                logging.warning(f"Could not determine whether alts are padded with N's (or n's)")
                alts_are_padded = None
        else:
            logging.warning(f"There are no alts that could be padded with N's (or n's)")
            alts_are_padded = None
        analysis = ReferenceGenomeFeatureAnalysis(
            has_unplaced_contigs,
            has_unlocalized_contigs,
            has_alts,
            has_decoys,
            has_patches,
            has_ebv,
            has_rcrs,
            uses_canonical_chrom_names,
            has_only_hardmasked_nucleotides_at_y_par1,
            has_semi_ambiguous_iub_codes,
            has_softmasked_nucleotides,
            alts_are_padded,
        )
        return analysis

    @classmethod
    def _mitochondrial_sequence_is_rcrs(
            cls,
            ref_genome_path: Path,
            rcrs_path: Path,
            ref_mitochondrial_contig_name: str,
    ) -> bool:
        do_decompress = rcrs_path.suffix == ".gz"
        if do_decompress:
            uncompressed_rcrs_path = rcrs_path.with_suffix("")
            if uncompressed_rcrs_path.exists():
                error_msg = (
                    f"Cannot decompress file {rcrs_path} without overwriting existing file {uncompressed_rcrs_path}"
                )
                raise FileExistsError(error_msg)
            with gzip.open(rcrs_path, "rb") as source_f, open(uncompressed_rcrs_path, "wb") as target_f:
                shutil.copyfileobj(source_f, target_f)
        else:
            uncompressed_rcrs_path = rcrs_path
        with pysam.Fastafile(uncompressed_rcrs_path) as rcrs_f:
            if rcrs_f.nreferences != 1:
                raise ValueError(f"rCRS FASTA file has more than one contig")
            rcrs_genome = rcrs_f.fetch(rcrs_f.references[0])
        if do_decompress:
            uncompressed_rcrs_path.unlink()
            (uncompressed_rcrs_path.parent / f"{uncompressed_rcrs_path.name}.fai") .unlink()
        with pysam.Fastafile(ref_genome_path) as genome_f:
            mitochondrial_from_ref = genome_f.fetch(ref_mitochondrial_contig_name)
        return bool(rcrs_genome == mitochondrial_from_ref)

    @classmethod
    def _get_y_test_sequence(cls, y_contig_name: str, ref_genome_path: Path) -> str:
        with pysam.Fastafile(ref_genome_path) as genome_f:
            return str(genome_f.fetch(y_contig_name, cls.Y_PAR1_TEST_REGION[0], cls.Y_PAR1_TEST_REGION[1]))

    @classmethod
    def _is_definitely_padded_with_n(cls, contig_name: str, ref_genome_path: Path) -> bool:
        with pysam.Fastafile(ref_genome_path) as genome_f:
            first_1000_nucleotides = get_nucleotides_from_string(
                genome_f.fetch(contig_name, 0, 1000)
            )
            last_1000_nucleotides = get_nucleotides_from_string(
                genome_f.fetch(contig_name, start=genome_f.get_reference_length(contig_name) - 1000)
            )

        first_1000_nucleotides_all_unknown = first_1000_nucleotides.issubset(UNKNOWN_NUCLEOTIDES)
        last_1000_nucleotides_all_unknown = last_1000_nucleotides.issubset(UNKNOWN_NUCLEOTIDES)
        return first_1000_nucleotides_all_unknown and last_1000_nucleotides_all_unknown
