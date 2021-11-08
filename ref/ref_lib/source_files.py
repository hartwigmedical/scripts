import logging
from enum import Enum, auto
from pathlib import Path
from typing import List, NamedTuple, overload, Any, Optional

from ref_lib.ref_util import make_temp_version_final, get_temp_path, download_bucket_file, download_file_over_https


class SourceFile(Enum):
    REFSEQ_FASTA = auto()
    DECOY_FASTA = auto()
    EBV_FASTA = auto()
    RCRS_FASTA = auto()
    REFSEQ_WITH_PATCHES_ASSEMBLY_REPORT = auto()
    REFSEQ_WITHOUT_PATCHES_ASSEMBLY_REPORT = auto()
    DECOY_ASSEMBLY_REPORT = auto()
    EBV_ASSEMBLY_REPORT = auto()

    @classmethod
    def get_all(cls) -> List["SourceFile"]:
        return [
            SourceFile.REFSEQ_FASTA,
            SourceFile.DECOY_FASTA,
            SourceFile.EBV_FASTA,
            SourceFile.RCRS_FASTA,
            SourceFile.REFSEQ_WITH_PATCHES_ASSEMBLY_REPORT,
            SourceFile.REFSEQ_WITHOUT_PATCHES_ASSEMBLY_REPORT,
            SourceFile.DECOY_ASSEMBLY_REPORT,
            SourceFile.EBV_ASSEMBLY_REPORT,
        ]

    @classmethod
    def get_required_for_contig_alias_file(cls) -> List["SourceFile"]:
        return [
            SourceFile.REFSEQ_WITH_PATCHES_ASSEMBLY_REPORT,
            SourceFile.REFSEQ_WITHOUT_PATCHES_ASSEMBLY_REPORT,
            SourceFile.DECOY_ASSEMBLY_REPORT,
            SourceFile.EBV_ASSEMBLY_REPORT,
        ]


class SourceFileLocator(object):
    SOURCE_FILE_TO_ORIGINAL_SOURCE = {
        SourceFile.REFSEQ_FASTA: "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.fna.gz",
        SourceFile.DECOY_FASTA: "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/786/075/GCA_000786075.2_hs38d1/GCA_000786075.2_hs38d1_genomic.fna.gz",
        SourceFile.EBV_FASTA: "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/402/265/GCF_002402265.1_ASM240226v1/GCF_002402265.1_ASM240226v1_genomic.fna.gz",
        SourceFile.RCRS_FASTA: "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_assembly_structure/non-nuclear/assembled_chromosomes/FASTA/chrMT.fna.gz",
        SourceFile.REFSEQ_WITH_PATCHES_ASSEMBLY_REPORT: "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_assembly_report.txt",
        SourceFile.REFSEQ_WITHOUT_PATCHES_ASSEMBLY_REPORT: "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.26_GRCh38/GCF_000001405.26_GRCh38_assembly_report.txt",
        SourceFile.DECOY_ASSEMBLY_REPORT: "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/786/075/GCA_000786075.2_hs38d1/GCA_000786075.2_hs38d1_assembly_report.txt",
        SourceFile.EBV_ASSEMBLY_REPORT: "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/402/265/GCF_002402265.1_ASM240226v1/GCF_002402265.1_ASM240226v1_assembly_report.txt",
    }

    # Overloading is for typing, it does not actually overload the operator
    @overload
    def get_location(self, source_file: SourceFile, base_dir: None = None) -> str:
        ...

    @overload
    def get_location(self, source_file: SourceFile, base_dir: str) -> str:
        ...

    @overload
    def get_location(self, source_file: SourceFile, base_dir: Path) -> Path:
        ...

    def get_location(self, source_file: Any, base_dir: Any = None) -> Any:
        original_source_location = self.SOURCE_FILE_TO_ORIGINAL_SOURCE[source_file]
        file_name = original_source_location.split("/")[-1]
        if base_dir is None:
            return original_source_location
        elif isinstance(base_dir, Path):
            return base_dir / file_name
        elif isinstance(base_dir, str):
            return f"{base_dir}/{file_name}"
        else:
            error_msg = "'get_location' only implemented for base_dir arguments that are None, Path or str."
            raise NotImplementedError(error_msg)


class DownloadJob(NamedTuple):
    source_file: SourceFile
    source: str
    target: Path


class SourceFileDownloader(object):
    SOURCES_LIST_FILE_NAME = "sources.txt"
    
    @classmethod
    def download_source_files(
            cls,
            source_files: List[SourceFile],
            target_dir: Path,
            bucket_dir: Optional[str] = None,
            create_file_with_sources: bool = True,
    ) -> None:
        logging.info(f"Starting download of source files: {[file.name for file in source_files]}")
        download_jobs = [
            DownloadJob(
                source_file,
                SourceFileLocator().get_location(source_file, bucket_dir),
                SourceFileLocator().get_location(source_file, target_dir)
            ) for source_file in source_files
        ]
        local_sources_list_file_path = target_dir / cls.SOURCES_LIST_FILE_NAME

        local_file_exists_list = [job.target.exists() for job in download_jobs]
        if create_file_with_sources:
            local_file_exists_list.append(local_sources_list_file_path.exists())

        if all(local_file_exists_list):
            logging.info("Skipping downloads. Source files already exist locally.")
            return
        elif any(local_file_exists_list):
            error_msg = (
                f"Some of the expected source files already exist locally and other do not. "
                f"Please delete the existing local source files so new version can be downloaded."
            )
            raise ValueError(error_msg)

        if not target_dir.is_dir():
            target_dir.mkdir(parents=True)

        if create_file_with_sources:
            if bucket_dir is None:
                logging.info(f"Writing original sources of files to a file: {local_sources_list_file_path}")
                cls._write_local_sources_list_file(download_jobs, local_sources_list_file_path)
            else:
                logging.info(f"Downloading file with original sources from bucket: {bucket_dir}")
                download_bucket_file(f"{bucket_dir}/{cls.SOURCES_LIST_FILE_NAME}", local_sources_list_file_path)

        logging.info("Starting downloads of source files")
        download_failed = False
        for job in download_jobs:
            logging.info(f"Start download of {job.source_file.name}")
            try:
                if bucket_dir is None:
                    logging.info(f"Download over https: {job.source}")
                    download_file_over_https(job.source, job.target)
                else:
                    logging.info(f"Download from bucket: {job.source}")
                    download_bucket_file(job.source, job.target)
            except Exception as exc:
                logging.error(
                    f"Download of {job.source_file.name} from {job.source} to {job.target} "
                    f"has generated an exception: {exc}"
                )
                download_failed = True
            else:
                logging.info(f"Finished download of {job.source_file.name}")
        if download_failed:
            raise ValueError("Download of at least one file has failed")
        else:
            logging.info(f"Finished downloads of source files: {[file.name for file in source_files]}")

    @classmethod
    def _write_local_sources_list_file(
            cls, download_jobs: List[DownloadJob], local_sources_list_file_path: Path,
    ) -> None:
        local_sources_list_text = "\n".join([f"{job.source_file.name}: {job.source}" for job in download_jobs])
        with open(get_temp_path(local_sources_list_file_path), "w") as f:
            f.write(local_sources_list_text)
        make_temp_version_final(local_sources_list_file_path)
