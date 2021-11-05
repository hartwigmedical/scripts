import argparse
import logging
import sys
from pathlib import Path
from typing import List, NamedTuple, Optional

import pysam

from contig_name_translation import AliasToCanonicalContigNameTextWriter, ContigNameTranslator
from contig_classification import ContigCategorizer
from contig_types import ContigTypeDesirabilities
from fasta_writer import FastaWriter
from ref_genome_feature_analysis import ReferenceGenomeFeatureAnalyzer, ReferenceGenomeFeatureAnalysis
from ref_util import set_up_logging, assert_dir_does_not_exist, assert_bucket_dir_does_not_exist, download_file_over_https, \
    upload_directory_to_bucket, get_temp_path, make_temp_version_final

# See gs://hmf-crunch-experiments/211005_david_DEV-2170_GRCh38-ref-genome-comparison/ for required files.

SCRIPT_NAME = "create_hmf_ref_genome_fasta"

REFSEQ_FASTA_SOURCE = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.fna.gz"
DECOY_FASTA_SOURCE = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/786/075/GCA_000786075.2_hs38d1/GCA_000786075.2_hs38d1_genomic.fna.gz"
EBV_FASTA_SOURCE = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/402/265/GCF_002402265.1_ASM240226v1/GCF_002402265.1_ASM240226v1_genomic.fna.gz"
RCRS_FASTA_SOURCE = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_assembly_structure/non-nuclear/assembled_chromosomes/FASTA/chrMT.fna.gz"

REFSEQ_WITH_PATCHES_ASSEMBLY_REPORT_SOURCE = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_assembly_report.txt"
REFSEQ_WITHOUT_PATCHES_ASSEMBLY_REPORT_SOURCE = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.26_GRCh38/GCF_000001405.26_GRCh38_assembly_report.txt"
DECOY_ASSEMBLY_REPORT_SOURCE = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/786/075/GCA_000786075.2_hs38d1/GCA_000786075.2_hs38d1_assembly_report.txt"
EBV_ASSEMBLY_REPORT_SOURCE = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/402/265/GCF_002402265.1_ASM240226v1/GCF_002402265.1_ASM240226v1_assembly_report.txt"

SOURCES_LIST_FILE_NAME = "sources.txt"

ALIAS_TO_CANONICAL_CONTIG_NAME_FILE_NAME = "alias_to_canonical_contig_name.tsv"
MASTER_FASTA_FILE_NAME = "master.fasta"


class Config(NamedTuple):
    working_dir: Path
    output_fasta_name: str
    output_bucket_dir: Optional[str]
    reuse_existing_files: bool

    def get_local_source_file_dir(self) -> Path:
        return self.working_dir / "source_files"

    def get_local_compressed_refseq_fasta_path(self) -> Path:
        return self.get_local_source_file_dir() / REFSEQ_FASTA_SOURCE.split("/")[-1]

    def get_local_compressed_decoy_fasta_path(self) -> Path:
        return self.get_local_source_file_dir() / DECOY_FASTA_SOURCE.split("/")[-1]

    def get_local_compressed_ebv_fasta_path(self) -> Path:
        return self.get_local_source_file_dir() / EBV_FASTA_SOURCE.split("/")[-1]

    def get_local_compressed_rcrs_fasta_path(self) -> Path:
        return self.get_local_source_file_dir() / RCRS_FASTA_SOURCE.split("/")[-1]

    def get_local_refseq_with_patches_assembly_report_path(self) -> Path:
        return self.get_local_source_file_dir() / REFSEQ_WITH_PATCHES_ASSEMBLY_REPORT_SOURCE.split("/")[-1]

    def get_local_refseq_without_patches_assembly_report_path(self) -> Path:
        return self.get_local_source_file_dir() / REFSEQ_WITHOUT_PATCHES_ASSEMBLY_REPORT_SOURCE.split("/")[-1]

    def get_local_decoy_assembly_report_path(self) -> Path:
        return self.get_local_source_file_dir() / DECOY_ASSEMBLY_REPORT_SOURCE.split("/")[-1]

    def get_local_ebv_assembly_report_path(self) -> Path:
        return self.get_local_source_file_dir() / EBV_ASSEMBLY_REPORT_SOURCE.split("/")[-1]

    def get_source_list_path(self) -> Path:
        return self.get_local_source_file_dir() / SOURCES_LIST_FILE_NAME

    def get_alias_to_canonical_contig_name_path(self) -> Path:
        return self.working_dir / ALIAS_TO_CANONICAL_CONTIG_NAME_FILE_NAME

    def get_local_uncompressed_master_fasta_path(self) -> Path:
        return self.working_dir / MASTER_FASTA_FILE_NAME

    def get_output_fasta_path(self) -> Path:
        return self.working_dir / self.output_fasta_name


def main(config: Config) -> None:
    set_up_logging()

    logging.info(f"Starting {SCRIPT_NAME}.")

    # Sanity checks
    if not config.reuse_existing_files:
        assert_dir_does_not_exist(config.working_dir)
    if config.output_bucket_dir is not None:
        assert_bucket_dir_does_not_exist(config.output_bucket_dir)

    if not config.get_local_source_file_dir().exists():
        config.get_local_source_file_dir().mkdir(parents=True)

    logging.info("Downloading source files.")
    download_source_files(config)

    logging.info(f"Creating {ALIAS_TO_CANONICAL_CONTIG_NAME_FILE_NAME} file.")
    if not config.get_alias_to_canonical_contig_name_path().exists():
        create_contig_alias_file(config)
    else:
        logging.info(f"Skipping creation of {ALIAS_TO_CANONICAL_CONTIG_NAME_FILE_NAME} file. Already exists.")

    logging.info(f"Creating master FASTA file.")
    if not config.get_local_uncompressed_master_fasta_path().exists():
        create_master_fasta_file(config)
    else:
        logging.info(f"Skipping creation of master FASTA file. Already exists.")

    logging.info(f"Creating temp version of HMFref FASTA file.")
    with open(config.get_alias_to_canonical_contig_name_path(), "r") as f:
        contig_name_translator = ContigNameTranslator.from_contig_alias_text(f.read())
    contig_categorizer = ContigCategorizer(contig_name_translator)
    contig_type_desirabilities = ContigTypeDesirabilities.create()

    FastaWriter.write_hmf_ref_genome_fasta(
        config.get_local_uncompressed_master_fasta_path(),
        get_temp_path(config.get_output_fasta_path()),
        contig_categorizer,
        contig_name_translator,
        contig_type_desirabilities,
    )

    logging.info("Asserting that output is as expected.")
    assert_created_ref_genome_matches_expectations(
        config, contig_categorizer, contig_name_translator, contig_type_desirabilities,
    )
    logging.info("Output is as expected.")

    logging.info("Moving temp output file to final output location.")
    make_temp_version_final(config.get_output_fasta_path())
    Path(f"{get_temp_path(config.get_output_fasta_path())}.fai").rename(Path(f"{config.get_output_fasta_path()}.fai"))

    if config.output_bucket_dir is not None:
        logging.info("Upload results to bucket.")
        upload_directory_to_bucket(config.working_dir, config.output_bucket_dir)
    else:
        logging.info("Skip upload of results to bucket.")

    # TODO: Make it possible to use sources from bucket dir
    # TODO: Remove unused scripts
    # TODO: Create requirements file to make this script properly reproducible with venv.
    #   Maybe add script to call venv and run script automatically. Maybe to create venv too, if needed.
    # TODO: Add option to exclude decoys
    # TODO: Add option to skip removing softmasks
    # TODO: Make check_ref_genome_features.py work in similar way to this script

    logging.info(f"Finished {SCRIPT_NAME}.")


def create_master_fasta_file(config: Config) -> None:
    FastaWriter.combine_compressed_files(
        [
            config.get_local_compressed_refseq_fasta_path(),
            config.get_local_compressed_decoy_fasta_path(),
            config.get_local_compressed_ebv_fasta_path(),
        ],
        get_temp_path(config.get_local_uncompressed_master_fasta_path()),
    )
    make_temp_version_final(config.get_local_uncompressed_master_fasta_path())


def create_contig_alias_file(config: Config) -> None:
    assembly_reports_text = "\n".join([
        get_text_from_file(config.get_local_refseq_with_patches_assembly_report_path()),
        get_text_from_file(config.get_local_refseq_without_patches_assembly_report_path()),
        get_text_from_file(config.get_local_decoy_assembly_report_path()),
        get_text_from_file(config.get_local_ebv_assembly_report_path()),
    ])
    contig_alias_text = AliasToCanonicalContigNameTextWriter.create_text_from_assembly_reports_text(
        assembly_reports_text
    )
    with open(get_temp_path(config.get_alias_to_canonical_contig_name_path()), "w") as f:
        f.write(contig_alias_text)
    make_temp_version_final(config.get_alias_to_canonical_contig_name_path())


def download_source_files(config: Config) -> None:
    download_name_source_target_triples = [
        ("REFSEQ_FASTA_SOURCE", REFSEQ_FASTA_SOURCE, config.get_local_compressed_refseq_fasta_path()),
        ("DECOY_FASTA_SOURCE", DECOY_FASTA_SOURCE, config.get_local_compressed_decoy_fasta_path()),
        ("EBV_FASTA_SOURCE", EBV_FASTA_SOURCE, config.get_local_compressed_ebv_fasta_path()),
        ("RCRS_FASTA_SOURCE", RCRS_FASTA_SOURCE, config.get_local_compressed_rcrs_fasta_path()),
        (
            "REFSEQ_WITH_PATCHES_ASSEMBLY_REPORT_SOURCE",
            REFSEQ_WITH_PATCHES_ASSEMBLY_REPORT_SOURCE,
            config.get_local_refseq_with_patches_assembly_report_path()
        ),
        (
            "REFSEQ_WITHOUT_PATCHES_ASSEMBLY_REPORT_SOURCE",
            REFSEQ_WITHOUT_PATCHES_ASSEMBLY_REPORT_SOURCE,
            config.get_local_refseq_without_patches_assembly_report_path()
        ),
        ("DECOY_ASSEMBLY_REPORT_SOURCE", DECOY_ASSEMBLY_REPORT_SOURCE, config.get_local_decoy_assembly_report_path()),
        ("EBV_ASSEMBLY_REPORT_SOURCE", EBV_ASSEMBLY_REPORT_SOURCE, config.get_local_ebv_assembly_report_path()),
    ]

    all_files_already_exist_locally = (
            all(triple[2].exists() for triple in download_name_source_target_triples)
            and config.get_source_list_path().exists()
    )
    some_files_already_exist_locally = (
            any(triple[2].exists() for triple in download_name_source_target_triples)
            or config.get_source_list_path().exists()
    )
    if all_files_already_exist_locally:
        logging.info("Skipping downloads. Source files already exist locally.")
        return
    elif some_files_already_exist_locally:
        error_msg = (
            f"Some of the expected source files already exist locally and other do not. "
            f"Please delete the existing local source files so new version can be downloaded."
        )
        raise ValueError(error_msg)

    logging.info(f"Writing sources to file: {config.get_source_list_path()}")
    lines = [f"{name}: {source}" for name, source, _ in download_name_source_target_triples]
    sources_text = "\n".join(lines)
    with open(get_temp_path(config.get_source_list_path()), "w") as f:
        f.write(sources_text)
    make_temp_version_final(config.get_source_list_path())

    logging.info("Starting downloads of source files")
    download_failed = False
    for name, source, target in download_name_source_target_triples:
        logging.info(f"Start download of {name}")
        try:
            download_file_over_https(source, target)
        except Exception as exc:
            logging.error(f"Download of {name} from {source} to {target} has generated an exception: {exc}")
            download_failed = True
        else:
            logging.info(f"Finished download of {name}")
    if download_failed:
        raise ValueError("Download of at least one file has failed")
    else:
        logging.info("Finished downloads of source files")


def assert_created_ref_genome_matches_expectations(
        config: Config,
        contig_categorizer: ContigCategorizer,
        contig_name_translator: ContigNameTranslator,
        contig_type_desirabilities: ContigTypeDesirabilities,
) -> None:
    with pysam.Fastafile(get_temp_path(config.get_output_fasta_path())) as temp_f:
        with pysam.Fastafile(config.get_local_uncompressed_master_fasta_path()) as master_f:
            contigs_expected_to_be_copied = {
                contig_name for contig_name in master_f.references
                if contig_categorizer.get_contig_type(contig_name) in contig_type_desirabilities.desired_contig_types
            }
            contigs_expected_in_output = {
                contig_name_translator.standardize(contig_name) for contig_name in contigs_expected_to_be_copied
            }
            if set(temp_f.references) != contigs_expected_in_output:
                error_msg = (
                    f"Contigs in output file not as expected: "
                    f"expected={contigs_expected_in_output}, actual={sorted(temp_f.references)}"
                )
                raise ValueError(error_msg)

            for contig_name in contigs_expected_to_be_copied:
                expected_sequence = master_f.fetch(contig_name).upper()
                actual_sequence = temp_f.fetch(contig_name_translator.standardize(contig_name))
                if actual_sequence != expected_sequence:
                    raise ValueError(f"Contig sequences are not identical: contig={contig_name}")
    feature_analysis = ReferenceGenomeFeatureAnalyzer.do_analysis(
        get_temp_path(config.get_output_fasta_path()),
        config.get_local_compressed_rcrs_fasta_path(),
        contig_name_translator,
    )
    expected_feature_analysis = ReferenceGenomeFeatureAnalysis(
        has_unplaced_contigs=True,
        has_unlocalized_contigs=True,
        has_alts=False,
        has_decoys=True,
        has_patches=False,
        has_ebv=True,
        has_rcrs=True,
        uses_canonical_chrom_names=True,
        has_only_hardmasked_nucleotides_at_y_par1=False,
        has_semi_ambiguous_iub_codes=True,
        has_softmasked_nucleotides=False,
        alts_are_padded=None,
    )
    if feature_analysis != expected_feature_analysis:
        error_msg = f"Not all features as expected: expected={expected_feature_analysis}, actual={feature_analysis}"
        raise ValueError(error_msg)


def get_text_from_file(path: Path) -> str:
    with open(path, "r") as f:
        return f.read().replace("\r", "")


def parse_args(sys_args: List[str]) -> Config:
    parser = argparse.ArgumentParser(
        prog=f"{SCRIPT_NAME}",
        description="Create FASTA file for GRCh38 HMF reference genome.",
    )
    parser.add_argument(
        "--working_dir",
        "-w",
        type=Path,
        required=True,
        help="Path to local working directory. If directory does not exist, it is created.",
    )
    parser.add_argument(
        "--output_fasta_name",
        "-o",
        type=str,
        required=True,
        help="Name of created FASTA file.",
    )
    parser.add_argument(
        "--output_bucket_dir",
        "-b",
        type=str,
        default=None,
        help=(
            "Optional argument. Path to output bucket dir. Argument should be of the form 'gs://some/kind/of/path'. "
            "This script will not overwrite existing files in the bucket."
        ),
    )
    parser.add_argument(
        "--reuse_existing_files",
        "-u",
        help=(
            "Optional argument. Reuse local source files from a previous run of this script."
        ),
        action="store_true"
    )

    args = parser.parse_args(sys_args)

    return Config(args.working_dir, args.output_fasta_name, args.output_bucket_dir, args.reuse_existing_files)


if __name__ == "__main__":
    main(parse_args(sys.argv[1:]))
