import argparse
import logging
import sys
from pathlib import Path
from typing import List, NamedTuple, Optional

from ref_lib.contig_name_translation import AliasToCanonicalContigNameTextWriter, ContigNameTranslator
from ref_lib.ref_genome_feature_analysis import ReferenceGenomeFeatureAnalyzer
from ref_lib.ref_util import set_up_logging, assert_dir_does_not_exist, ALIAS_TO_CANONICAL_CONTIG_NAME_FILE_NAME, \
    SOURCE_FILES_DIR_NAME, assert_file_exists
from ref_lib.source_files import SourceFile, SourceFileDownloader, SourceFileLocator

SCRIPT_NAME = "check_ref_genome_features"


class Config(NamedTuple):
    ref_genome_path: Path
    working_dir: Path
    reuse_existing_files: bool
    source_files_from_bucket_dir: Optional[str]

    def get_local_source_file_dir(self) -> Path:
        return self.working_dir / SOURCE_FILES_DIR_NAME

    def get_alias_to_canonical_contig_name_path(self) -> Path:
        return self.working_dir / ALIAS_TO_CANONICAL_CONTIG_NAME_FILE_NAME


def main(config: Config) -> None:
    set_up_logging()

    logging.info(f"Starting {SCRIPT_NAME}.")

    logging.info(f"Config values:\n{config}")

    # Sanity checks
    assert_file_exists(config.ref_genome_path)

    if not config.reuse_existing_files:
        assert_dir_does_not_exist(config.working_dir)

    if not config.get_local_source_file_dir().exists():
        config.get_local_source_file_dir().mkdir(parents=True)

    logging.info("Downloading source files.")
    SourceFileDownloader.download_source_files(
        SourceFile.get_required_for_contig_alias_file() + [SourceFile.RCRS_FASTA],
        config.get_local_source_file_dir(),
        config.source_files_from_bucket_dir,
        create_file_with_sources=False,
    )

    logging.info(f"Creating {ALIAS_TO_CANONICAL_CONTIG_NAME_FILE_NAME} file.")
    if not config.get_alias_to_canonical_contig_name_path().exists():
        AliasToCanonicalContigNameTextWriter.create_contig_alias_file(
            config.get_alias_to_canonical_contig_name_path(),
            config.get_local_source_file_dir(),
        )
    else:
        logging.info(f"Skipping creation of {ALIAS_TO_CANONICAL_CONTIG_NAME_FILE_NAME} file. Already exists.")

    logging.info(f"Checking features.")
    with open(config.get_alias_to_canonical_contig_name_path(), "r") as f:
        contig_name_translator = ContigNameTranslator.from_contig_alias_text(f.read())

    rcrs_path = SourceFileLocator().get_location(SourceFile.RCRS_FASTA, config.get_local_source_file_dir())
    analysis = ReferenceGenomeFeatureAnalyzer.do_analysis(
        config.ref_genome_path, rcrs_path, contig_name_translator,
    )

    logging.info(f"FEATURES GENOME:")

    logging.info(f"Unplaced contigs: {analysis.has_unplaced_contigs}")
    logging.info(f"Unlocalized contigs: {analysis.has_unlocalized_contigs}")
    logging.info(f"Alts: {analysis.has_alts}")
    logging.info(f"Decoys (hs38d1): {analysis.has_decoys}")
    logging.info(f"Patches: {analysis.has_patches}")
    logging.info(f"EBV: {analysis.has_ebv}")
    logging.info(f"rCRS mitochondrial sequence: {analysis.has_rcrs}")
    logging.info(f"Uses canonical contig names, so 'chr1' etc.: {analysis.uses_canonical_chrom_names}")
    logging.info(f"PAR hardmask (not fully accurate): {analysis.has_only_hardmasked_nucleotides_at_y_par1}")
    logging.info(f"Semi ambiguous IUB codes: {analysis.has_semi_ambiguous_iub_codes}")
    logging.info(f"Has softmasked nucleotides: {analysis.has_softmasked_nucleotides}")
    logging.info(f"Alts are padded with N: {analysis.alts_are_padded}")
    logging.info(f"PhiX: False?")
    logging.info(f"")
    logging.info(f"For easy copy-paste:")
    answers = [
        analysis.has_unplaced_contigs,
        analysis.has_unlocalized_contigs,
        analysis.has_alts,
        analysis.has_decoys,
        analysis.has_patches,
        analysis.has_ebv,
        analysis.has_rcrs,
        analysis.uses_canonical_chrom_names,
        analysis.has_only_hardmasked_nucleotides_at_y_par1,
        analysis.has_semi_ambiguous_iub_codes,
        analysis.has_softmasked_nucleotides,
        analysis.alts_are_padded,
    ]
    value_to_answer = {True: "Yes", False: "No", None: "?"}
    print("\n".join([value_to_answer[answer] for answer in answers]))
    logging.info(f"Finished {SCRIPT_NAME}.")


def parse_args(sys_args: List[str]) -> Config:
    parser = argparse.ArgumentParser(
        prog=f"{SCRIPT_NAME}",
        description=(
            "Check important features for ref genome FASTA file. Correctness is not guaranteed, especially for hg19."
        ),
    )
    parser.add_argument(
        "--ref_genome",
        "-i",
        type=Path,
        required=True,
        help="Path to FASTA file with reference genome to check.",
    )
    parser.add_argument(
        "--working_dir",
        "-w",
        type=Path,
        required=True,
        help="Path to local working directory. If directory does not exist, it is created.",
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

    args = parser.parse_args(sys_args)

    config = Config(
        args.ref_genome,
        args.working_dir,
        args.reuse_existing_files,
        args.source_files_from_bucket_dir,
    )
    return config


if __name__ == "__main__":
    main(parse_args(sys.argv[1:]))
