import logging
from copy import deepcopy
from pathlib import Path
from typing import Dict

from google.cloud import storage


class ContigNameTranslator(object):
    """Standardizes names if it can. Returns argument as is if it cannot."""
    def __init__(self, contig_name_to_canonical_name: Dict[str, str]) -> None:
        self._contig_name_to_canonical_name = deepcopy(contig_name_to_canonical_name)

    @classmethod
    def from_blob(cls, blob: storage.Blob) -> "ContigNameTranslator":
        contig_name_to_canonical_name = {}
        for line in blob.download_as_text().split("\n"):
            split_line = line.split("\t")
            if len(split_line) != 2:
                raise ValueError(f"Incorrect length line: {line}")
            contig_name, canonical_name = split_line
            if contig_name in contig_name_to_canonical_name:
                raise ValueError(f"Encountered contig name multiple times: {contig_name}")
            contig_name_to_canonical_name[contig_name] = canonical_name
        return ContigNameTranslator(contig_name_to_canonical_name)

    def standardize(self, contig_name: str) -> str:
        if contig_name in self._contig_name_to_canonical_name:
            return self._contig_name_to_canonical_name[contig_name]
        else:
            raise ValueError(f"Could not standardize '{contig_name}'")

    def is_canonical(self, contig_name: str) -> bool:
        return contig_name in self._contig_name_to_canonical_name.values()


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