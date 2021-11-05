import argparse
import logging
import sys
from pathlib import Path
from typing import List, NamedTuple

from ref_lib.contig_name_translation import AliasToCanonicalContigNameTextWriter
from ref_lib.ref_util import set_up_logging, assert_file_does_not_exist, get_text_from_bucket_file

SCRIPT_NAME = "create_chrom_translation_file"
SOURCE_FILE_BUCKET = "hmf-crunch-experiments"
NON_DECOY_OLD_TRANSLATION_FILE_PATH = "211005_david_DEV-2170_GRCh38-ref-genome-comparison/GCA_000001405.15_GRCh38_assembly_report.txt"
# original source: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/GCA_000001405.15_GRCh38_assembly_report.txt
NON_DECOY_NEW_TRANSLATION_FILE_PATH = "211005_david_DEV-2170_GRCh38-ref-genome-comparison/GCA_000001405.28_GRCh38.p13_assembly_report.txt"
# original source: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.28_GRCh38.p13/GCA_000001405.28_GRCh38.p13_assembly_report.txt
DECOY_TRANSLATION_FILE_PATH = "211005_david_DEV-2170_GRCh38-ref-genome-comparison/GCA_000786075.2_hs38d1_assembly_report.txt"
# original source: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/786/075/GCA_000786075.2_hs38d1/GCA_000786075.2_hs38d1_assembly_report.txt


class Config(NamedTuple):
    output_path: Path

    def validate(self) -> None:
        assert_file_does_not_exist(self.output_path)


def main(config: Config) -> None:
    set_up_logging()

    logging.info(f"Starting {SCRIPT_NAME}")

    config.validate()

    assembly_reports_text = get_translation_text_from_bucket_files()

    contig_alias_text = AliasToCanonicalContigNameTextWriter.create_text_from_assembly_reports_text(
        assembly_reports_text,
    )

    with open(config.output_path, "w") as f:
        f.write(contig_alias_text)

    logging.info(f"Finished {SCRIPT_NAME}")


def get_translation_text_from_bucket_files() -> str:
    non_decoy_old_text = get_text_from_bucket_file(f"gs://{SOURCE_FILE_BUCKET}/{NON_DECOY_OLD_TRANSLATION_FILE_PATH}")
    non_decoy_new_text = get_text_from_bucket_file(f"gs://{SOURCE_FILE_BUCKET}/{NON_DECOY_NEW_TRANSLATION_FILE_PATH}")
    decoy_text = get_text_from_bucket_file(f"gs://{SOURCE_FILE_BUCKET}/{DECOY_TRANSLATION_FILE_PATH}")
    combined_text = f"{non_decoy_old_text}\n{non_decoy_new_text}\n{decoy_text}"
    return combined_text


def parse_args(sys_args: List[str]) -> Config:
    parser = argparse.ArgumentParser(
        prog=f"{SCRIPT_NAME}",
        description=(
            "Create translation file for ref genome fasta file comparison."
        ),
    )
    parser.add_argument(
        "--output",
        "-o",
        type=Path,
        required=True,
        help="Output file path.",
    )

    args = parser.parse_args(sys_args)

    config = Config(Path(args.output))
    return config


if __name__ == "__main__":
    main(parse_args(sys.argv[1:]))
