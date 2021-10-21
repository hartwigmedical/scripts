import concurrent.futures
import logging
from pathlib import Path
from typing import Set

import pysam
from google.cloud import storage

STANDARD_NUCLEOTIDES = {"A", "C", "G", "T", "N"}
SOFTMASKED_NUCLEOTIDES = {"a", "c", "g", "t", "n"}
UNKNOWN_NUCLEOTIDES = {"N", "n"}


def set_up_logging() -> None:
    logging.basicConfig(
        format="%(asctime)s - [%(levelname)-8s] - %(message)s", level=logging.INFO, datefmt="%Y-%m-%d %H:%M:%S"
    )


def assert_file_exists(path: Path) -> None:
    if not path.is_file():
        raise ValueError(f"File does not exist: {path}")


def assert_file_does_not_exist(path: Path) -> None:
    if path.is_file():
        raise ValueError(f"File exists: {path}")


def assert_dir_does_not_exist(path: Path) -> None:
    if path.is_dir():
        raise ValueError(f"Dir exists: {path}")


def delete_if_exists(path: Path) -> None:
    if path.exists():
        path.unlink()


def assert_file_exists_in_bucket(path: str) -> None:
    if not get_blob(path).exists():
        raise ValueError(f"File in bucket does not exist: {path}")


def get_blob(path: str) -> storage.Blob:
    bucket_name = path.split("/")[2]
    relative_path = "/".join(path.split("/")[3:])
    return storage.Client().get_bucket(bucket_name).get_blob(relative_path)


def get_nucleotides_from_fasta(fasta_path: Path) -> Set[str]:
    futures = []
    with concurrent.futures.ProcessPoolExecutor() as executor:
        with pysam.Fastafile(fasta_path) as genome_f:
            for contig_name in genome_f.references:
                contig = genome_f.fetch(contig_name)
                futures.append(executor.submit(get_nucleotides_from_string, contig))

    nucleotides = set()
    for future in futures:
        try:
            nucleotides = nucleotides.union(future.result())
        except Exception as exc:
            raise ValueError(exc)
    return nucleotides


def get_nucleotides_from_string(sequence: str) -> Set[str]:
    nucleotides = set()
    for nucleotide in sequence:
        nucleotides.add(nucleotide)
    return nucleotides
