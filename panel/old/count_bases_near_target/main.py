import argparse
import logging
import multiprocessing
import subprocess
import sys
from pathlib import Path
from typing import List, NamedTuple, Set

import pysam

WIDER_BED = "wider.bed"
SAMTOOLS_FILTERED_BAM = "samtools_filtered.bam"
SAMTOOLS_FILTERED_BAM_INDEX = f"{SAMTOOLS_FILTERED_BAM}.bai"
PYTHON_FILTERED_BAM = "python_filtered.bam"
PYTHON_FILTERED_BAM_INDEX = f"{PYTHON_FILTERED_BAM}.bai"
DEPTH_FILE = "relevant.depth"
COUNT_FILE = "count.txt"

MAX_TARGET_DISTANCE = 500
THREAD_COUNT = multiprocessing.cpu_count()


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
    def samtools_filtered_bam_index_path(self) -> Path:
        return self.working_directory / SAMTOOLS_FILTERED_BAM_INDEX

    @property
    def python_filtered_bam_path(self) -> Path:
        return self.working_directory / PYTHON_FILTERED_BAM

    @property
    def python_filtered_bam_index_path(self) -> Path:
        return self.working_directory / PYTHON_FILTERED_BAM_INDEX

    @property
    def depth_path(self) -> Path:
        return self.working_directory / DEPTH_FILE

    @property
    def count_path(self) -> Path:
        return self.working_directory / COUNT_FILE


def main(config: Config) -> None:
    set_up_logging()
    config.validate()

    logging.info(f"Started count_bases_near_target for {config.bam_path}")

    config.working_directory.mkdir(parents=True, exist_ok=True)

    if not config.count_path.exists():
        if not config.wider_bed_path.exists():
            logging.info(f"Creating wider bed file")
            delete_if_exists(config.samtools_filtered_bam_path)
            create_wider_bed(config)
            assert config.wider_bed_path.exists(), "Wider bed creation failed"
        else:
            logging.info(f"Wider bed file already exists")

        if not config.samtools_filtered_bam_path.exists():
            logging.info(f"Creating samtools-filtered bam file")
            delete_if_exists(config.samtools_filtered_bam_index_path)
            create_samtools_filtered_bam(config)
            assert config.samtools_filtered_bam_path.exists(), "Samtools filtering failed"
        else:
            logging.info(f"Samtools-filtered bam file already exists")

        if not config.samtools_filtered_bam_index_path.exists():
            logging.info(f"Creating samtools-filtered bam index file")
            delete_if_exists(config.python_filtered_bam_path)
            create_bam_index(config.samtools_filtered_bam_path, config.samtools)
            assert config.samtools_filtered_bam_index_path.exists(), "Samtools filtered bam indexing failed"
        else:
            logging.info(f"Samtools-filtered bam file already exists")

        if not config.python_filtered_bam_path.exists():
            logging.info(f"Creating python-filtered bam file")
            delete_if_exists(config.python_filtered_bam_index_path)
            create_python_filtered_bam(config)
            assert config.python_filtered_bam_path.exists(), "Python filtering failed"
        else:
            logging.info(f"Python-filtered bam file already exists")

        if not config.python_filtered_bam_index_path.exists():
            logging.info(f"Creating python-filtered bam index file")
            delete_if_exists(config.depth_path)
            create_bam_index(config.python_filtered_bam_path, config.samtools)
            assert config.python_filtered_bam_index_path.exists(), "Python filtered bam indexing failed"
        else:
            logging.info(f"Python-filtered bam file already exists")

        if not config.depth_path.exists():
            logging.info(f"Creating samtools depth file")
            delete_if_exists(config.count_path)
            create_depth_file(config)
            assert config.depth_path.exists(), "Samtools depth failed"
        else:
            logging.info(f"Samtools depth file already exists")

        if not config.count_path.exists():
            logging.info(f"Creating count file")
            create_count_file(config)
            assert config.count_path.exists(), "Counting bases failed"
        else:
            logging.info(f"Count file already exists")

    if config.count_path.exists():
        print_base_count(config)

        delete_if_exists(config.wider_bed_path)
        delete_if_exists(config.samtools_filtered_bam_path)
        delete_if_exists(config.samtools_filtered_bam_index_path)
        delete_if_exists(config.python_filtered_bam_path)
        delete_if_exists(config.python_filtered_bam_index_path)
        delete_if_exists(config.depth_path)
    else:
        raise FileNotFoundError(f"Config count file {config.count_path} does not exist.")


def create_wider_bed(config: Config) -> None:
    with open(config.bed_path, "r") as input_f:
        with open(config.wider_bed_path, "w") as output_f:
            for line in input_f:
                chrom, start, end = line.split("\t")[0:3]
                wider_start = str(max(int(start) - MAX_TARGET_DISTANCE, 0))
                wider_end = str(int(end) + MAX_TARGET_DISTANCE)
                new_line = "\t".join([chrom, wider_start, wider_end]) + "\n"
                output_f.write(new_line)


def create_samtools_filtered_bam(config: Config) -> None:
    cli_args = [
        config.samtools,
        "view",
        "--bam",
        "--with-header",
        "--target-file",
        config.wider_bed_path,
        "--threads",
        str(THREAD_COUNT - 1),
        "--exclude-flags",
        "DUP",
        config.bam_path,
    ]
    with open(config.samtools_filtered_bam_path, "w") as f:
        subprocess.run(cli_args, stdout=f)


def create_python_filtered_bam(config: Config) -> None:
    relevant_read_names: Set[str] = set()
    with pysam.AlignmentFile(config.samtools_filtered_bam_path, "rb", threads=THREAD_COUNT) as filter_f:
        for read in filter_f:
            relevant_read_names.add(read.reference_name)

    with pysam.AlignmentFile(config.bam_path, "rb", threads=THREAD_COUNT) as input_f:
        with pysam.AlignmentFile(config.python_filtered_bam_path, "wb", template=input_f, threads=THREAD_COUNT) as output_f:
            for read in input_f:
                if read.reference_name in relevant_read_names:
                    output_f.write(read)


def create_depth_file(config: Config) -> None:
    cli_args = [
        config.samtools,
        "depth",
        "-s",
        config.python_filtered_bam_path,
    ]
    with open(config.depth_path, "w") as f:
        subprocess.run(cli_args, stdout=f)


def create_bam_index(bam_path: Path, samtools: Path) -> None:
    cli_args = [samtools, "index", "-@", str(THREAD_COUNT - 1), bam_path]
    subprocess.run(cli_args)


def create_count_file(config: Config) -> None:
    cli_args = ["awk", '{sum+=$3;} END{printf "%.0f", sum;}', config.depth_path]
    with open(config.count_path, "w") as f:
        subprocess.run(cli_args, stdout=f)


def print_base_count(config: Config) -> None:
    with open(config.count_path, "r") as f:
        logging.info(f"Total number of useful bases close to target: {f.read()}")


def set_up_logging() -> None:
    logging.basicConfig(
        format="%(asctime)s - [%(levelname)-8s] - %(message)s", level=logging.INFO, datefmt="%Y-%m-%d %H:%M:%S"
    )


def delete_if_exists(path: Path) -> None:
    if path.exists():
        path.unlink()


def parse_args(sys_args: List[str]) -> Config:
    parser = argparse.ArgumentParser(
        prog="count_bases_near_target",
        description=(
            "Count bases of read pairs for which part of at least one read is close to a region in a target bed file."
        ),
    )
    parser.add_argument("--bam", "-i", type=str, required=True, help="Input bam.")
    parser.add_argument("--bed", "-b", type=str, required=True, help="Bed file of target.")
    parser.add_argument("--samtools", "-s", type=str, required=True, help="Samtools. Version 1.13 or greater.")
    parser.add_argument("--working_dir", "-d", type=str, required=True, help="Working dir.")

    args = parser.parse_args(sys_args)

    config = Config(Path(args.bam), Path(args.bed), Path(args.samtools), Path(args.working_dir))
    return config


if __name__ == "__main__":
    main(parse_args(sys.argv[1:]))
