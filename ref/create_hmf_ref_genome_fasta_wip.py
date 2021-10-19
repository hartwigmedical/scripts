import argparse
import gzip
import logging
import re
import shutil
import sys
from pathlib import Path
from typing import List, NamedTuple

import requests

from ref_util import set_up_logging, ContigNameTranslator, assert_file_exists_in_bucket, get_blob, \
    assert_file_does_not_exist, assert_dir_does_not_exist


# See gs://hmf-crunch-experiments/211005_david_DEV-2170_GRCh38-ref-genome-comparison/ for required files.

SCRIPT_NAME = "create_hmf_ref_genome_fasta"

REFSEQ_FASTA_SOURCE = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.26_GRCh38/GCF_000001405.26_GRCh38_genomic.fna.gz"
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
        if not re.fullmatch(r"gs://.+", self.contig_alias_bucket_path):
            raise ValueError(f"Contig alias bucket path is not of the form 'gs://some/file/path'")
        assert_file_exists_in_bucket(self.contig_alias_bucket_path)
        # assert_dir_does_not_exist(self.working_dir)  # TODO: reenable
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


def main(config: Config) -> None:
    set_up_logging()

    logging.info(f"Starting {SCRIPT_NAME}")

    config.validate()

    config.working_dir.mkdir(parents=True)

    with open(config.get_source_list_path(), "w") as f:
        logging.info(f"Writing sources to file: {config.get_source_list_path()}")
        text = (
            f"REFSEQ FASTA SOURCE: {REFSEQ_FASTA_SOURCE}\n"
            f"DECOY FASTA SOURCE: {DECOY_FASTA_SOURCE}\n"
            f"EBV FASTA SOURCE: {EBV_FASTA_SOURCE}"
        )
        f.write(text)

    logging.info(f"Creating contig name translator")
    contig_name_translator = ContigNameTranslator.from_blob(get_blob(config.contig_alias_bucket_path))

    # get_local_copy_fasta_file(REFSEQ_FASTA_SOURCE, config.get_refseq_fasta_path())  # TODO: reenable
    # get_local_copy_fasta_file(DECOY_FASTA_SOURCE, config.get_decoy_fasta_path())  # TODO: reenable
    # get_local_copy_fasta_file(EBV_FASTA_SOURCE, config.get_ebv_fasta_path())  # TODO: reenable

    # logging.info(f"Combining FASTA files into a single file.")
    # combine_files(
    #     [config.get_refseq_fasta_path(), config.get_decoy_fasta_path(), config.get_ebv_fasta_path()],
    #     config.get_master_fasta_path(),  # TODO: reenable
    # )


    # TODO: Extract contig names from master file, and figure out which contigs we want
    # TODO: Put these contigs into output file, including proper header lines for each contig,
    #   and filtering out decoys and removing soft-masks if configured to do this
    #       TODO: Figure out what these header liens should look like

    logging.info(f"Finished {SCRIPT_NAME}")


def get_local_copy_fasta_file(source: str, target: Path) -> None:
    compressed_target = target.parent / f"{target.name}.gz"

    logging.info(f"Downloading compressed file {source} to {compressed_target}.")
    download_file(source, compressed_target)

    logging.info(f"Decompressing downloaded file {compressed_target}.")
    decompress(compressed_target, target)


def download_file(source: str, target: Path) -> None:
    response = requests.get(source, stream=True)
    response.raise_for_status()
    with open(target, 'wb') as f:
        for block in response.iter_content(1024):
            f.write(block)


def decompress(source: Path, target: Path) -> None:
    with gzip.open(source, "rb") as f_in:
        with open(target, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)


def combine_files(sources: List[Path], target: Path) -> None:
    with open(target, "wb") as f_out:
        for source in sources:
            with open(source, "rb") as f_in:
                shutil.copyfileobj(f_in, f_out)


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
