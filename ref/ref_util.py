import concurrent.futures
import gzip
import logging
import re
import shutil
from pathlib import Path
from typing import Set, Tuple

import pysam
import requests
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


def assert_dir_exists(path: Path) -> None:
    if not path.is_dir():
        raise ValueError(f"Dir does not exist: {path}")


def delete_if_exists(path: Path) -> None:
    if path.exists():
        path.unlink()


def assert_file_exists_in_bucket(bucket_path: str) -> None:
    if not re.fullmatch(r"gs://.+", bucket_path):
        raise ValueError(f"Path is not of the form 'gs://some/file/path': {bucket_path}")
    if not get_blob(bucket_path).exists():
        raise ValueError(f"File in bucket does not exist: {bucket_path}")


def assert_bucket_dir_does_not_exist(bucket_path: str) -> None:
    if not re.fullmatch(r"gs://.+", bucket_path):
        raise ValueError(f"Path is not of the form 'gs://some/file/path': {bucket_path}")
    bucket_name, relative_path = split_bucket_path(bucket_path)
    if list(storage.Client().get_bucket(bucket_name).list_blobs(prefix=relative_path)):
        raise ValueError(f"Bucket dir exists: {bucket_path}")


def upload_directory_to_bucket(source_dir: Path, bucket_dir: str) -> None:
    assert_dir_exists(source_dir)

    futures = []
    with concurrent.futures.ProcessPoolExecutor() as executor:
        for file_path in source_dir.glob("**/*"):
            if file_path.is_file():
                target_path = f"{bucket_dir}/{file_path.relative_to(source_dir)}"
                futures.append(executor.submit(upload_file_to_bucket, file_path, target_path))

    for future in concurrent.futures.as_completed(futures):
        try:
            future.result()
        except Exception as exc:
            raise ValueError(exc)


def upload_file_to_bucket(source_path: Path, bucket_path: str) -> None:
    logging.info(f"Uploading {source_path} to {bucket_path}")
    assert_file_exists(source_path)
    blob = get_blob(bucket_path)
    if not blob.exists():
        blob.upload_from_filename(str(source_path))
        logging.info(f"Finished upload of {source_path} to {bucket_path}")
    else:
        raise FileExistsError(f"Cannot upload file {source_path} since it would overwrite an existing file.")


def get_text_from_bucket_file(path: str) -> str:
    text: str = get_blob(path).download_as_text()
    return text


def get_blob(path: str) -> storage.Blob:
    bucket_name, relative_path = split_bucket_path(path)
    return storage.Client().get_bucket(bucket_name).blob(relative_path)


def split_bucket_path(path: str) -> Tuple[str, str]:
    bucket_name = path.split("/")[2]
    relative_path = "/".join(path.split("/")[3:])
    return bucket_name, relative_path


def get_nucleotides_from_fasta(fasta_path: Path) -> Set[str]:
    futures = []
    with concurrent.futures.ProcessPoolExecutor() as executor:
        with pysam.Fastafile(fasta_path) as genome_f:
            for contig_name in genome_f.references:
                contig = genome_f.fetch(contig_name)
                futures.append(executor.submit(get_nucleotides_from_string, contig))

    nucleotides: Set[str] = set()
    for future in concurrent.futures.as_completed(futures):
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


def download_file_over_https(source: str, target: Path) -> None:
    response = requests.get(source, stream=True)
    response.raise_for_status()
    with open(get_temp_path(target), 'wb') as f:
        for block in response.iter_content(1024):
            f.write(block)
    make_temp_version_final(target)


def decompress(source: Path, target: Path) -> None:
    with gzip.open(source, "rb") as f_in:
        with open(target, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)


def make_temp_version_final(path: Path) -> None:
    temp_path = get_temp_path(path)
    logging.info(f"Started moving temporary file {temp_path} to final location {path}")
    if not temp_path.exists():
        error_msg = (
            f"Cannot make temporary version of '{path}' permanent "
            f"since the temporary version '{temp_path}' doesn't exist."
        )
        raise ValueError(error_msg)
    temp_path.rename(path)
    logging.info(f"Finished moving temporary file {temp_path} to final location {path}")


def get_temp_path(path: Path) -> Path:
    return path.parent / f"{path.name}.tmp"
