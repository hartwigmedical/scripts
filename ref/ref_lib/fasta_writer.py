import gzip
import hashlib
import logging
import shutil
from pathlib import Path
from typing import List

import pysam

from ref_lib.contig_classification import ContigCategorizer
from ref_lib.contig_name_translation import ContigNameTranslator
from ref_lib.contig_types import ContigType, ContigTypeDesirabilities, Assembly


class FastaWriter(object):
    FASTA_LINE_WRAP = 70
    FASTA_HEADER_SEPARATOR = "  "
    CONTIG_TYPE_TO_ROLE = {
        ContigType.AUTOSOME: "Chromosome",
        ContigType.X: "Chromosome",
        ContigType.Y: "Chromosome",
        ContigType.MITOCHONDRIAL: "Mitochondrion",
        ContigType.EBV: "decoy",
        ContigType.DECOY: "decoy",
        ContigType.UNLOCALIZED: "unlocalized",
        ContigType.UNPLACED: "unplaced",
        ContigType.ALT: "alt-scaffold",
        ContigType.FIX_PATCH: "fix-patch",
        ContigType.NOVEL_PATCH: "novel-patch",
    }

    @classmethod
    def combine_files(cls, sources: List[Path], target: Path) -> None:
        with open(target, "wb") as f_out:
            for source in sources:
                with open(source, "rb") as f_in:
                    shutil.copyfileobj(f_in, f_out)

    @classmethod
    def combine_compressed_files(cls, sources: List[Path], target: Path) -> None:
        with open(target, "wb") as f_out:
            for source in sources:
                with gzip.open(source, "rb") as f_in:
                    shutil.copyfileobj(f_in, f_out)

    @classmethod
    def write_hmf_ref_genome_fasta(
            cls,
            master_fasta: Path,
            target_fasta: Path,
            contig_categorizer: ContigCategorizer,
            contig_name_translator: ContigNameTranslator,
            contig_type_desirabilities: ContigTypeDesirabilities,
    ) -> None:
        cls._assert_master_fasta_contig_types_match_expected(
            master_fasta,
            contig_categorizer,
            contig_type_desirabilities,
        )

        with pysam.Fastafile(master_fasta) as master_f:
            contig_names = list(master_f.references)

            with open(target_fasta, "w") as out_f:
                for contig_name in contig_names:
                    logging.info(f"Handling {contig_name}")
                    contig_type = contig_categorizer.get_contig_type(contig_name)
                    if contig_type in contig_type_desirabilities.desired_contig_types:
                        logging.info(f"Include {contig_name} in output file")
                        sequence = master_f.fetch(contig_name).upper()
                        header = cls._get_header(
                            contig_name,
                            contig_name_translator.standardize(contig_name),
                            contig_type,
                            sequence,
                        )
                        logging.info(f"Header: {header}")
                        out_f.write(header + "\n")
                        for i in range(0, len(sequence), cls.FASTA_LINE_WRAP):
                            out_f.write(sequence[i:i + cls.FASTA_LINE_WRAP] + "\n")

    @classmethod
    def _assert_master_fasta_contig_types_match_expected(
            cls,
            master_fasta_path: Path,
            contig_categorizer: ContigCategorizer,
            contig_type_desirabilities: ContigTypeDesirabilities,
    ) -> None:
        with pysam.Fastafile(master_fasta_path) as master_f:
            seen_contig_types = {contig_categorizer.get_contig_type(name) for name in master_f.references}

        expected_contig_types = contig_type_desirabilities.get_expected_contig_types()
        if seen_contig_types != expected_contig_types:
            sorted_seen_contig_names = sorted(contig_type.name for contig_type in seen_contig_types)
            sorted_expected_contig_names = sorted(
                contig_type.name for contig_type in expected_contig_types
            )
            error_msg = (
                f"Seen contig types don't perfectly match the expected contig types: "
                f"seen={sorted_seen_contig_names}, expected={sorted_expected_contig_names}"
            )
            raise ValueError(error_msg)

    @classmethod
    def _get_header(cls, contig_name: str, standardized_contig_name: str, contig_type: ContigType, sequence: str) -> str:
        header_entries = [f">{standardized_contig_name}", f"AC:{contig_name}", f"LN:{len(sequence)}"]
        if contig_type == ContigType.UNLOCALIZED:
            region = standardized_contig_name.split("_")[0]
            header_entries.append(f"rg:{region}")
        elif contig_type in {ContigType.ALT, ContigType.FIX_PATCH, ContigType.NOVEL_PATCH}:
            # TODO: Add rg field entry for ALT and patch contigs, if we want these in ref genome.
            #       Would need the coordinates of where they align.
            #       See https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GRCh38_major_release_seqs_for_alignment_pipelines/README_analysis_sets.txt
            error_msg = (
                "Proper version of 'rg' field in FASTA header not implemented for ALT, FIX_PATCH and NOVEL_PATCH"
            )
            raise NotImplementedError(error_msg)
        header_entries.append(f"rl:{cls._get_contig_role(contig_type)}")
        header_entries.append(f"M5:{hashlib.md5(sequence.encode('utf-8')).hexdigest()}")
        assembly = contig_type.get_assembly()
        if assembly == Assembly.GRCH38:
            header_entries.append("AS:GRCh38")
        elif assembly == Assembly.HS38D1:
            header_entries.append("AS:hs38d1")
        elif assembly == Assembly.OTHER:
            pass
        else:
            raise ValueError(f"Unrecognized assembly: {assembly}")
        if contig_type == ContigType.EBV:
            header_entries.append("SP:Human_herpesvirus_4")
        if contig_type in {ContigType.MITOCHONDRIAL, ContigType.EBV}:
            header_entries.append("tp:circular")
        return cls.FASTA_HEADER_SEPARATOR.join(header_entries)

    @classmethod
    def _get_contig_role(cls, contig_type: ContigType) -> str:
        if contig_type in cls.CONTIG_TYPE_TO_ROLE.keys():
            return cls.CONTIG_TYPE_TO_ROLE[contig_type]
        else:
            raise ValueError(f"Encountered contig type without an assigned role: {contig_type}")
