import argparse
import logging
import sys
from pathlib import Path
from typing import List, NamedTuple, Optional

import pysam

from ref_lib.contig_name_translation import AliasToCanonicalContigNameTextWriter, ContigNameTranslator
from ref_lib.contig_classification import ContigCategorizer
from ref_lib.contig_types import ContigTypeDesirabilities
from ref_lib.fasta_writer import FastaWriter
from ref_lib.ref_genome_feature_analysis import ReferenceGenomeFeatureAnalyzer, ReferenceGenomeFeatureAnalysis
from ref_lib.ref_util import set_up_logging, assert_dir_does_not_exist, assert_bucket_dir_does_not_exist, \
    upload_directory_to_bucket, get_temp_path, make_temp_version_final, ALIAS_TO_CANONICAL_CONTIG_NAME_FILE_NAME, \
    MASTER_FASTA_FILE_NAME, SOURCE_FILES_DIR_NAME, combine_compressed_files
from ref_lib.source_files import SourceFile, SourceFileDownloader, SourceFileLocator

SCRIPT_NAME = "create_hmf_ref_genome_fasta"


class Config(NamedTuple):
    working_dir: Path
    output_fasta_name: str
    output_bucket_dir: Optional[str]
    reuse_existing_files: bool
    source_files_from_bucket_dir: Optional[str]

    def get_local_source_file_dir(self) -> Path:
        return self.working_dir / SOURCE_FILES_DIR_NAME

    def get_alias_to_canonical_contig_name_path(self) -> Path:
        return self.working_dir / ALIAS_TO_CANONICAL_CONTIG_NAME_FILE_NAME

    def get_local_uncompressed_master_fasta_path(self) -> Path:
        return self.working_dir / MASTER_FASTA_FILE_NAME

    def get_output_fasta_path(self) -> Path:
        return self.working_dir / self.output_fasta_name


def main(config: Config) -> None:
    set_up_logging()

    logging.info(f"Starting {SCRIPT_NAME}.")

    logging.info(f"Config values:\n{config}")

    # Sanity checks
    if not config.reuse_existing_files:
        assert_dir_does_not_exist(config.working_dir)
    if config.output_bucket_dir is not None:
        assert_bucket_dir_does_not_exist(config.output_bucket_dir)

    if not config.get_local_source_file_dir().exists():
        config.get_local_source_file_dir().mkdir(parents=True)

    logging.info("Downloading source files.")
    SourceFileDownloader.download_source_files(
        SourceFile.get_all(),
        config.get_local_source_file_dir(),
        config.source_files_from_bucket_dir,
    )

    logging.info(f"Creating {ALIAS_TO_CANONICAL_CONTIG_NAME_FILE_NAME} file.")
    if not config.get_alias_to_canonical_contig_name_path().exists():
        AliasToCanonicalContigNameTextWriter.create_contig_alias_file(
            config.get_alias_to_canonical_contig_name_path(),
            config.get_local_source_file_dir(),
        )
    else:
        logging.info(f"Skipping creation of {ALIAS_TO_CANONICAL_CONTIG_NAME_FILE_NAME} file. Already exists.")

    logging.info(f"Creating master FASTA file.")
    if not config.get_local_uncompressed_master_fasta_path().exists():
        create_master_fasta_file(config)
    else:
        logging.info(f"Skipping creation of master FASTA file. Already exists.")

    logging.info(f"Creating temp version of HMFref FASTA file.")
    with open(config.get_alias_to_canonical_contig_name_path(), "r") as f:
        contig_name_translator = ContigNameTranslator.from_contig_alias_text(f.read())
    contig_categorizer = ContigCategorizer(contig_name_translator)
    contig_type_desirabilities = ContigTypeDesirabilities.create()

    FastaWriter.write_hmf_ref_genome_fasta(
        config.get_local_uncompressed_master_fasta_path(),
        get_temp_path(config.get_output_fasta_path()),
        contig_categorizer,
        contig_name_translator,
        contig_type_desirabilities,
    )

    logging.info("Asserting that output is as expected.")
    assert_created_ref_genome_matches_expectations(
        config, contig_categorizer, contig_name_translator, contig_type_desirabilities,
    )
    logging.info("Output is as expected.")

    logging.info("Moving temp output file to final output location.")
    make_temp_version_final(config.get_output_fasta_path())
    Path(f"{get_temp_path(config.get_output_fasta_path())}.fai").rename(Path(f"{config.get_output_fasta_path()}.fai"))

    if config.output_bucket_dir is not None:
        logging.info("Upload results to bucket.")
        upload_directory_to_bucket(config.working_dir, config.output_bucket_dir)
    else:
        logging.info("Skip upload of results to bucket.")

    # TODO: Maybe add option to keep softmasks
    # TODO: Maybe add option to exclude decoys

    logging.info(f"Finished {SCRIPT_NAME}.")


def create_master_fasta_file(config: Config) -> None:
    combine_compressed_files(
        [
            SourceFileLocator().get_location(source_file, config.get_local_source_file_dir())
            for source_file in [SourceFile.REFSEQ_FASTA, SourceFile.DECOY_FASTA, SourceFile.EBV_FASTA]
        ],
        get_temp_path(config.get_local_uncompressed_master_fasta_path()),
    )
    make_temp_version_final(config.get_local_uncompressed_master_fasta_path())


def assert_created_ref_genome_matches_expectations(
        config: Config,
        contig_categorizer: ContigCategorizer,
        contig_name_translator: ContigNameTranslator,
        contig_type_desirabilities: ContigTypeDesirabilities,
) -> None:
    with pysam.Fastafile(get_temp_path(config.get_output_fasta_path())) as temp_f:
        with pysam.Fastafile(config.get_local_uncompressed_master_fasta_path()) as master_f:
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
        get_temp_path(config.get_output_fasta_path()),
        SourceFileLocator().get_location(SourceFile.RCRS_FASTA, config.get_local_source_file_dir()),
        contig_name_translator,
    )
    expected_feature_analysis = ReferenceGenomeFeatureAnalysis(
        has_unplaced_contigs=True,
        has_unlocalized_contigs=True,
        has_alts=False,
        has_decoys=True,
        has_patches=False,
        has_ebv=True,
        has_rcrs=True,
        uses_canonical_chrom_names=True,
        has_only_hardmasked_nucleotides_at_y_par1=False,
        has_semi_ambiguous_iub_codes=True,
        has_softmasked_nucleotides=False,
        alts_are_padded=None,
        has_ki270752=False,
    )
    if feature_analysis != expected_feature_analysis:
        error_msg = f"Not all features as expected: expected={expected_feature_analysis}, actual={feature_analysis}"
        raise ValueError(error_msg)


def parse_args(sys_args: List[str]) -> Config:
    parser = argparse.ArgumentParser(
        prog=f"{SCRIPT_NAME}",
        description="Create FASTA file for GRCh38 HMF reference genome.",
    )
    parser.add_argument(
        "--working_dir",
        "-w",
        type=Path,
        required=True,
        help="Path to local working directory. If directory does not exist, it is created.",
    )
    parser.add_argument(
        "--output_fasta_name",
        "-o",
        type=str,
        required=True,
        help="Name of created FASTA file.",
    )
    parser.add_argument(
        "--output_bucket_dir",
        "-b",
        type=str,
        default=None,
        help=(
            "Optional argument. Path to output bucket dir. Argument should be of the form 'gs://some/kind/of/path'. "
            "This script will not overwrite existing files in the bucket."
        ),
    )
    parser.add_argument(
        "--reuse_existing_files",
        "-u",
        help="Optional argument. Reuse local source files from a previous run of this script.",
        action="store_true",
    )
    parser.add_argument(
        "--source_files_from_bucket_dir",
        "-s",
        type=str,
        default=None,
        help=(
            "Optional argument. Get source files from this GCP bucket directory "
            "instead of downloading them from their original source."
        ),
    )
    parser.add_argument(
        "--venv_dir",
        "-v",
        type=Path,
        default=None,
        help="Optional argument. Directory to use for a Python venv. If not provided, no venv is used.",
    )  # Does nothing. Should be caught and used in shell script already.

    args = parser.parse_args(sys_args)

    if args.venv_dir is not None:
        error_msg = (
            f"Venv argument should be caught and used in the wrapping shell script, "
            f"not in the Python script itself: {args.venv_dir}"
        )
        raise SyntaxError(error_msg)

    config = Config(
        args.working_dir,
        args.output_fasta_name,
        args.output_bucket_dir,
        args.reuse_existing_files,
        args.source_files_from_bucket_dir,
    )
    return config


if __name__ == "__main__":
    main(parse_args(sys.argv[1:]))
