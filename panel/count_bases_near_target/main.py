import argparse
import logging
import subprocess
import sys
from pathlib import Path
from typing import List, NamedTuple, Set

import pysam

WIDER_BED = "wider.bed"
SAMTOOLS_FILTERED_BAM = "samtools_filtered.bam"
PYTHON_FILTERED_BAM = "python_filtered.bam"
DEPTH_FILE = "relevant.depth"

MAX_TARGET_DISTANCE = 500


class Config(NamedTuple):
    bam_path: Path
    bed_path: Path
    samtools: Path
    working_directory: Path

    def validate(self) -> None:
        if self.bam_path.suffix != ".bam":
            raise ValueError(f"Bam argument is not a bam file: {self.bam_path}")
        if self.bed_path.suffix != ".bed":
            raise ValueError(f"Bed argument is not a bed file: {self.bed_path}")
        if not self.bam_path.is_file():
            raise ValueError(f"Bam file does not exist: {self.bam_path}")
        if not self.bed_path.is_file():
            raise ValueError(f"Bed file does not exist: {self.bed_path}")
        if not self.samtools.is_file():
            raise ValueError(f"Samtools path does not exist: {self.samtools}")

    @property
    def wider_bed_path(self) -> Path:
        return self.working_directory / WIDER_BED

    @property
    def samtools_filtered_bam_path(self) -> Path:
        return self.working_directory / SAMTOOLS_FILTERED_BAM

    @property
    def python_filtered_bam_path(self) -> Path:
        return self.working_directory / PYTHON_FILTERED_BAM

    @property
    def depth_path(self) -> Path:
        return self.working_directory / DEPTH_FILE


def main(config: Config) -> None:
    set_up_logging()
    config.validate()

    logging.info(f"Started count_bases_near_target for {config.bam_path}")

    config.working_directory.mkdir(parents=True, exist_ok=True)

    if not config.wider_bed_path.exists():
        logging.info(f"Creating wider bed file")
        config.samtools_filtered_bam_path.unlink()
        create_wider_bed(config)
        assert config.wider_bed_path.exists(), "Wider bed creation failed"
    else:
        logging.info(f"Wider bed file already exists")

    if not config.samtools_filtered_bam_path.exists():
        logging.info(f"Creating samtools-filtered bam file")
        config.python_filtered_bam_path.unlink()
        create_samtools_filtered_bam(config)
        assert config.samtools_filtered_bam_path.exists(), "Samtools filtering failed"
    else:
        logging.info(f"Samtools-filtered bam file already exists")

    if not config.python_filtered_bam_path.exists():
        logging.info(f"Creating python-filtered bam file")
        config.depth_path.unlink()
        create_python_filtered_bam(config)
        assert config.python_filtered_bam_path.exists(), "Python filtering failed"
    else:
        logging.info(f"Python-filtered bam file already exists")

    if not config.depth_path.exists():
        logging.info(f"Creating samtools depth file")
        create_depth_file(config)
        assert config.depth_path.exists(), "Samtools depth failed"
    else:
        logging.info(f"Samtools depth file already exists")

    logging.info(f"Calculating total number of useful bases close to target:")
    get_base_count(config)


def create_wider_bed(config: Config) -> None:
    with open(config.bed_path, "r") as input_f:
        with open(config.wider_bed_path, "w") as output_f:
            lines = input_f.read()
            new_lines = []
            for line in lines:
                chrom, start, end, section_id = line.split("\t")
                wider_start = max(start - MAX_TARGET_DISTANCE, 0)
                wider_end = end + MAX_TARGET_DISTANCE
                new_line = "\t".join([chrom, wider_start, wider_end, section_id])
                new_lines.append(new_line)
            output_f.write("\n".join(new_lines))


def create_samtools_filtered_bam(config: Config) -> None:
    cli_args = [
        config.samtools,
        "view",
        "--bam",
        "--with-header",
        "--target-file",
        config.wider_bed_path,
        config.bam_path,
        ">",
        config.samtools_filtered_bam_path
    ]
    subprocess.run(cli_args)


def create_python_filtered_bam(config: Config) -> None:
    relevant_read_names: Set[str] = set()
    with pysam.AlignmentFile(config.samtools_filtered_bam_path, "rb") as filter_f:
        for read in filter_f:
            relevant_read_names.add(read.reference_name)

    with pysam.AlignmentFile(config.bam_path, "rb") as input_f:
        with pysam.AlignmentFile(config.python_filtered_bam_path, "wb") as output_f:
            for read in input_f:
                if read.reference_name in relevant_read_names:
                    output_f.write(read)


def create_depth_file(config: Config) -> None:
    cli_args = [
        config.samtools,
        "depth",
        "-s",
        config.python_filtered_bam_path,
        ">",
        config.depth_path
    ]
    subprocess.run(cli_args)


def get_base_count(config: Config) -> None:
    cli_args = ["awk", '{sum+=$3;} END{printf "%.0f", sum;}', config.depth_path]
    subprocess.run(cli_args)


def set_up_logging() -> None:
    logging.basicConfig(
        format="%(asctime)s - [%(levelname)-8s] - %(message)s", level=logging.INFO, datefmt="%Y-%m-%d %H:%M:%S"
    )


def parse_args(sys_args: List[str]) -> Config:
    parser = argparse.ArgumentParser(
        prog="count_bases_near_target",
        description=(
            "Count bases of read pairs for which part of at least one read is close to a region in a target bed file."
        ),
    )
    parser.add_argument("--bam", "-i", type=str, required=True, help="Input bam.")
    parser.add_argument("--bed", "-b", type=str, required=True, help="Bed file of target.")
    parser.add_argument("--samtools", "-s", type=str, required=True, help="Samtools.")
    parser.add_argument("--working_dir", "-d", type=str, required=True, help="Working dir.")

    args = parser.parse_args(sys_args)

    config = Config(Path(args.bam), Path(args.bed), Path(args.samtools), Path(args.working_dir))
    return config


if __name__ == "__main__":
    main(parse_args(sys.argv[1:]))
