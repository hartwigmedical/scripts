import argparse
import logging
import sys
from pathlib import Path
from typing import List, NamedTuple

import pysam

from contig_classification import ContigCategorizer
from contig_name_translation import ContigNameTranslator
from contig_types import ContigTypeDesirabilities
from fasta_writer import FastaWriter
from ref_genome_feature_analysis import ReferenceGenomeFeatureAnalyzer, ReferenceGenomeFeatureAnalysis
from ref_util import set_up_logging, assert_file_exists_in_bucket, \
    assert_file_does_not_exist, assert_dir_does_not_exist, download_file_over_https, decompress, get_text_from_bucket_file

# See gs://hmf-crunch-experiments/211005_david_DEV-2170_GRCh38-ref-genome-comparison/ for required files.

SCRIPT_NAME = "create_hmf_ref_genome_fasta_wip"

REFSEQ_FASTA_SOURCE = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.fna.gz"
DECOY_FASTA_SOURCE = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/786/075/GCA_000786075.2_hs38d1/GCA_000786075.2_hs38d1_genomic.fna.gz"
EBV_FASTA_SOURCE = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/402/265/GCF_002402265.1_ASM240226v1/GCF_002402265.1_ASM240226v1_genomic.fna.gz"

REFSEQ_FASTA_FILE_NAME = "REFSEQ.fasta"
DECOY_FASTA_FILE_NAME = "DECOY.fasta"
EBV_FASTA_FILE_NAME = "EBV.fasta"
MASTER_FASTA_FILE_NAME = "MASTER.fasta"
SOURCES_LIST_FILE_NAME = "sources.txt"


class Config(NamedTuple):
    contig_alias_bucket_path: str
    working_dir: Path
    output_path: Path

    def validate(self) -> None:
        assert_file_exists_in_bucket(self.contig_alias_bucket_path)
        assert_dir_does_not_exist(self.working_dir)
        assert_file_does_not_exist(self.output_path)

    def get_source_list_path(self) -> Path:
        return self.working_dir / SOURCES_LIST_FILE_NAME

    def get_refseq_fasta_path(self) -> Path:
        return self.working_dir / REFSEQ_FASTA_FILE_NAME

    def get_decoy_fasta_path(self) -> Path:
        return self.working_dir / DECOY_FASTA_FILE_NAME

    def get_ebv_fasta_path(self) -> Path:
        return self.working_dir / EBV_FASTA_FILE_NAME

    def get_master_fasta_path(self) -> Path:
        return self.working_dir / MASTER_FASTA_FILE_NAME

    def get_temp_output_path(self) -> Path:
        return self.output_path.parent / f"{self.output_path.name}.tmp"


def main(config: Config) -> None:
    set_up_logging()

    logging.info(f"Starting {SCRIPT_NAME}")

    config.validate()

    config.working_dir.mkdir(parents=True, exist_ok=True)

    with open(config.get_source_list_path(), "w") as f:
        logging.info(f"Writing sources to file: {config.get_source_list_path()}")
        text = (
            f"REFSEQ FASTA SOURCE: {REFSEQ_FASTA_SOURCE}\n"
            f"DECOY FASTA SOURCE: {DECOY_FASTA_SOURCE}\n"
            f"EBV FASTA SOURCE: {EBV_FASTA_SOURCE}"
        )
        f.write(text)

    logging.info(f"Creating contig name translator")
    contig_name_translator = ContigNameTranslator.from_contig_alias_text(
        get_text_from_bucket_file(config.contig_alias_bucket_path),
    )

    get_local_copy_fasta_file(REFSEQ_FASTA_SOURCE, config.get_refseq_fasta_path())
    get_local_copy_fasta_file(DECOY_FASTA_SOURCE, config.get_decoy_fasta_path())
    get_local_copy_fasta_file(EBV_FASTA_SOURCE, config.get_ebv_fasta_path())

    contig_categorizer = ContigCategorizer(contig_name_translator)
    contig_type_desirabilities = ContigTypeDesirabilities.create()

    logging.info(f"Combining FASTA files into a single file.")
    FastaWriter.combine_files(
        [config.get_refseq_fasta_path(), config.get_decoy_fasta_path(), config.get_ebv_fasta_path()],
        config.get_master_fasta_path(),
    )

    logging.info(f"Creating temp version of output FASTA file: {config.get_temp_output_path()}.")
    FastaWriter.write_hmf_ref_genome_fasta(
        config.get_master_fasta_path(),
        config.get_temp_output_path(),
        contig_categorizer,
        contig_name_translator,
        contig_type_desirabilities,
    )

    logging.info("Asserting that output is as expected")
    assert_created_ref_genome_matches_expectations(
        config, contig_categorizer, contig_name_translator, contig_type_desirabilities,
    )

    logging.info("Moving temp output file to real output file")
    config.get_temp_output_path().rename(config.output_path)

    logging.info("Clean up ")
    Path(f"{config.get_temp_output_path()}.fai").unlink()
    # TODO: Add option to exclude decoys
    # TODO: Add option to skip removing softmasks
    # TODO: Put used files and result in bucket? Or just do this manually afterwards.
    # TODO: Create requirements file to make this script properly reproducible with venv.
    #       Maybe add script to call venv and run script automatically. Maybe to create venv too, if needed.

    logging.info(f"Finished {SCRIPT_NAME}")


def assert_created_ref_genome_matches_expectations(
        config: Config,
        contig_categorizer: ContigCategorizer,
        contig_name_translator: ContigNameTranslator,
        contig_type_desirabilities: ContigTypeDesirabilities,
) -> None:
    with pysam.Fastafile(config.get_temp_output_path()) as temp_f:
        with pysam.Fastafile(config.get_master_fasta_path()) as master_f:
            contigs_expected_to_be_copied = {
                contig_name for contig_name in master_f.references
                if contig_categorizer.get_contig_type(contig_name) in contig_type_desirabilities.desired_contig_types
            }
            contigs_expected_in_output = {
                contig_name_translator.standardize(contig_name) for contig_name in contigs_expected_to_be_copied
            }
            if set(temp_f.references) != contigs_expected_in_output:
                error_msg = (
                    f"Contigs in output file not as expected: "
                    f"expected={contigs_expected_in_output}, actual={sorted(temp_f.references)}"
                )
                raise ValueError(error_msg)

            for contig_name in contigs_expected_to_be_copied:
                expected_sequence = master_f.fetch(contig_name).upper()
                actual_sequence = temp_f.fetch(contig_name_translator.standardize(contig_name))
                if actual_sequence != expected_sequence:
                    raise ValueError(f"Contig sequences are not identical: contig={contig_name}")
    feature_analysis = ReferenceGenomeFeatureAnalyzer.do_analysis(
        config.get_temp_output_path(), None, contig_name_translator,
    )
    expected_feature_analysis = ReferenceGenomeFeatureAnalysis(
        has_unplaced_contigs=True,
        has_unlocalized_contigs=True,
        has_alts=False,
        has_decoys=True,
        has_patches=False,
        has_ebv=True,
        has_rcrs=None,
        uses_canonical_chrom_names=True,
        has_only_hardmasked_nucleotides_at_y_par1=False,
        has_semi_ambiguous_iub_codes=True,
        has_softmasked_nucleotides=False,
        alts_are_padded=None,
    )
    if feature_analysis != expected_feature_analysis:
        error_msg = f"Not all features as expected: expected={expected_feature_analysis}, actual={feature_analysis}"
        raise ValueError(error_msg)


def get_local_copy_fasta_file(source: str, target: Path) -> None:
    compressed_target = target.parent / f"{target.name}.gz"

    logging.info(f"Downloading compressed file {source} to {compressed_target}.")
    download_file_over_https(source, compressed_target)

    logging.info(f"Decompressing downloaded file {compressed_target}.")
    decompress(compressed_target, target)


def parse_args(sys_args: List[str]) -> Config:
    parser = argparse.ArgumentParser(
        prog="create_hmf_ref_genome_fasta",
        description="Create FASTA file for GRCh38 HMF reference genome.",
    )
    parser.add_argument(
        "--contig_alias",
        "-c",
        type=str,
        required=True,
        help="Bucket path to TSV file with contig name translations. Source: create_contig_translation_file.py.",
    )
    parser.add_argument(
        "--working_dir",
        "-w",
        type=Path,
        required=True,
        help="Path to working directory. If directory does not exist, it is created.",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=Path,
        required=True,
        help="Path to output file.",
    )

    args = parser.parse_args(sys_args)

    return Config(args.contig_alias, args.working_dir, args.output)


if __name__ == "__main__":
    main(parse_args(sys.argv[1:]))
