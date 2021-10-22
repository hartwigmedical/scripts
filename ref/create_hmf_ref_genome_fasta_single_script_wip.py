import argparse
import logging
import sys
from pathlib import Path
from typing import List, NamedTuple

from ref_util import set_up_logging, assert_dir_does_not_exist, assert_bucket_dir_does_not_exist, download_file

# See gs://hmf-crunch-experiments/211005_david_DEV-2170_GRCh38-ref-genome-comparison/ for required files.

SCRIPT_NAME = "create_hmf_ref_genome_fasta"

REFSEQ_FASTA_SOURCE = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.fna.gz"
DECOY_FASTA_SOURCE = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/786/075/GCA_000786075.2_hs38d1/GCA_000786075.2_hs38d1_genomic.fna.gz"
EBV_FASTA_SOURCE = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/402/265/GCF_002402265.1_ASM240226v1/GCF_002402265.1_ASM240226v1_genomic.fna.gz"
RCRS_FASTA_SOURCE = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_assembly_structure/non-nuclear/assembled_chromosomes/FASTA/chrMT.fna.gz"

REFSEQ_ASSEMBLY_REPORT_SOURCE = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_assembly_report.txt"
DECOY_ASSEMBLY_REPORT_SOURCE = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/786/075/GCA_000786075.2_hs38d1/GCA_000786075.2_hs38d1_assembly_report.txt"
EBV_ASSEMBLY_REPORT_SOURCE = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/402/265/GCF_002402265.1_ASM240226v1/GCF_002402265.1_ASM240226v1_assembly_report.txt"

SOURCES_LIST_FILE_NAME = "sources.txt"


class Config(NamedTuple):
    output_bucket_dir: str
    working_dir: Path

    def get_local_source_file_dir(self) -> Path:
        return self.working_dir / "source_files"

    def get_local_compressed_refseq_fasta_path(self) -> Path:
        return self.get_local_source_file_dir() / REFSEQ_FASTA_SOURCE.split("/")[-1]

    def get_local_compressed_decoy_fasta_path(self) -> Path:
        return self.get_local_source_file_dir() / DECOY_FASTA_SOURCE.split("/")[-1]

    def get_local_compressed_ebv_fasta_path(self) -> Path:
        return self.get_local_source_file_dir() / EBV_FASTA_SOURCE.split("/")[-1]

    def get_local_compressed_rcrs_fasta_path(self) -> Path:
        return self.get_local_source_file_dir() / "rCRS.fasta"

    def get_local_refseq_assembly_report_path(self) -> Path:
        return self.get_local_source_file_dir() / REFSEQ_ASSEMBLY_REPORT_SOURCE.split("/")[-1]

    def get_local_decoy_assembly_report_path(self) -> Path:
        return self.get_local_source_file_dir() / DECOY_ASSEMBLY_REPORT_SOURCE.split("/")[-1]

    def get_local_ebv_assembly_report_path(self) -> Path:
        return self.get_local_source_file_dir() / EBV_ASSEMBLY_REPORT_SOURCE.split("/")[-1]

    def get_source_list_path(self) -> Path:
        return self.get_local_source_file_dir() / SOURCES_LIST_FILE_NAME


def main(config: Config) -> None:
    set_up_logging()

    logging.info(f"Starting {SCRIPT_NAME}")

    # Sanity checks
    assert_bucket_dir_does_not_exist(config.output_bucket_dir)
    assert_dir_does_not_exist(config.working_dir)

    config.get_local_source_file_dir().mkdir(parents=True)
    logging.info("Downloading source files")
    download_file(REFSEQ_FASTA_SOURCE, config.get_local_compressed_refseq_fasta_path())
    download_file(DECOY_FASTA_SOURCE, config.get_local_compressed_decoy_fasta_path())
    download_file(EBV_FASTA_SOURCE, config.get_local_compressed_ebv_fasta_path())
    download_file(RCRS_FASTA_SOURCE, config.get_local_compressed_rcrs_fasta_path())
    download_file(REFSEQ_ASSEMBLY_REPORT_SOURCE, config.get_local_refseq_assembly_report_path())
    download_file(DECOY_ASSEMBLY_REPORT_SOURCE, config.get_local_decoy_assembly_report_path())
    download_file(EBV_ASSEMBLY_REPORT_SOURCE, config.get_local_ebv_assembly_report_path())

    with open(config.get_source_list_path(), "w") as f:
        logging.info(f"Writing sources to file: {config.get_source_list_path()}")
        lines = [
            f"REFSEQ_FASTA_SOURCE: {REFSEQ_FASTA_SOURCE}",
            f"DECOY_FASTA_SOURCE: {DECOY_FASTA_SOURCE}",
            f"EBV_FASTA_SOURCE: {EBV_FASTA_SOURCE}",
            f"RCRS_FASTA_SOURCE: {RCRS_FASTA_SOURCE}",
            f"REFSEQ_ASSEMBLY_REPORT_SOURCE: {REFSEQ_ASSEMBLY_REPORT_SOURCE}",
            f"DECOY_ASSEMBLY_REPORT_SOURCE: {DECOY_ASSEMBLY_REPORT_SOURCE}",
            f"EBV_ASSEMBLY_REPORT_SOURCE: {EBV_ASSEMBLY_REPORT_SOURCE}",
        ]
        f.write("\n".join(lines))

    logging.info(f"Finished {SCRIPT_NAME}")


def parse_args(sys_args: List[str]) -> Config:
    parser = argparse.ArgumentParser(
        prog=f"{SCRIPT_NAME}",
        description="Create FASTA file for GRCh38 HMF reference genome.",
    )
    parser.add_argument(
        "--output_bucket_dir",
        "-o",
        type=str,
        required=True,
        help="Path to output bucket dir. Argument should be of the form 'gs://some/kind/of/path'.",
    )
    parser.add_argument(
        "--working_dir",
        "-w",
        type=Path,
        required=True,
        help="Path to local working directory. If directory does not exist, it is created.",
    )

    args = parser.parse_args(sys_args)

    return Config(args.output_bucket_dir, args.working_dir)


if __name__ == "__main__":
    main(parse_args(sys.argv[1:]))
