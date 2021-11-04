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
from ref_util import set_up_logging, assert_dir_does_not_exist, assert_bucket_dir_does_not_exist, download_file, \
    upload_directory_to_bucket

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
OUTPUT_FASTA_FILE_NAME = "output.fasta"


class Config(NamedTuple):
    working_dir: Path
    output_bucket_dir: Optional[str]

    def get_local_source_file_dir(self) -> Path:
        return self.working_dir / "source_files"

    def get_local_compressed_refseq_fasta_path(self) -> Path:
        return self.get_local_source_file_dir() / REFSEQ_FASTA_SOURCE.split("/")[-1]

    def get_local_compressed_decoy_fasta_path(self) -> Path:
        return self.get_local_source_file_dir() / DECOY_FASTA_SOURCE.split("/")[-1]

    def get_local_compressed_ebv_fasta_path(self) -> Path:
        return self.get_local_source_file_dir() / EBV_FASTA_SOURCE.split("/")[-1]

    def get_local_compressed_rcrs_fasta_path(self) -> Path:
        return self.get_local_source_file_dir() / "rCRS.fasta.gz"

    def get_local_uncompressed_rcrs_fasta_path(self) -> Path:
        return self.get_local_source_file_dir() / "rCRS.fasta"

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

    def get_temp_output_fasta_path(self) -> Path:
        return self.working_dir / f"{OUTPUT_FASTA_FILE_NAME}.tmp"

    def get_output_fasta_path(self) -> Path:
        return self.working_dir / f"{OUTPUT_FASTA_FILE_NAME}"


def main(config: Config) -> None:
    set_up_logging()

    logging.info(f"Starting {SCRIPT_NAME}.")

    # Sanity checks
    assert_dir_does_not_exist(config.working_dir)
    if config.output_bucket_dir is not None:
        assert_bucket_dir_does_not_exist(config.output_bucket_dir)

    config.get_local_source_file_dir().mkdir(parents=True)

    logging.info("Downloading source files.")
    download_source_files(config)

    logging.info(f"Creating {ALIAS_TO_CANONICAL_CONTIG_NAME_FILE_NAME} file.")
    assembly_reports_text = "\n".join([
        get_text_from_file(config.get_local_refseq_with_patches_assembly_report_path()),
        get_text_from_file(config.get_local_refseq_without_patches_assembly_report_path()),
        get_text_from_file(config.get_local_decoy_assembly_report_path()),
        get_text_from_file(config.get_local_ebv_assembly_report_path()),
    ])

    contig_alias_text = AliasToCanonicalContigNameTextWriter.create_text_from_assembly_reports_text(
        assembly_reports_text
    )

    with open(config.get_alias_to_canonical_contig_name_path(), "w") as f:
        f.write(contig_alias_text)

    logging.info(f"Creating master FASTA file.")
    FastaWriter.combine_compressed_files(
        [
            config.get_local_compressed_refseq_fasta_path(),
            config.get_local_compressed_decoy_fasta_path(),
            config.get_local_compressed_ebv_fasta_path(),
        ],
        config.get_local_uncompressed_master_fasta_path(),
    )

    logging.info(f"Creating temp version HMFref FASTA file.")
    contig_name_translator = ContigNameTranslator.from_contig_alias_text(contig_alias_text)
    contig_categorizer = ContigCategorizer(contig_name_translator)
    contig_type_desirabilities = ContigTypeDesirabilities.create()

    FastaWriter.write_hmf_ref_genome_fasta(
        config.get_local_uncompressed_master_fasta_path(),
        config.get_temp_output_fasta_path(),
        contig_categorizer,
        contig_name_translator,
        contig_type_desirabilities,
    )

    logging.info("Asserting that output is as expected.")
    assert_created_ref_genome_matches_expectations(
        config, contig_categorizer, contig_name_translator, contig_type_desirabilities,
    )
    logging.info("Output is as expected.")

    logging.info("Moving temp output file to real output file.")
    config.get_temp_output_fasta_path().rename(config.get_output_fasta_path())

    logging.info("Clean up.")
    Path(f"{config.get_temp_output_fasta_path()}.fai").unlink()

    if config.output_bucket_dir is not None:
        logging.info("Upload results to bucket.")
        upload_directory_to_bucket(config.working_dir, config.output_bucket_dir)
    else:
        logging.info("Skip upload of results to bucket.")

    # TODO: Do all downloading and writing etc. to a temp version of the file first,
    #   then change the name when it succeeds.
    #   Should the temp file contain a random substring or something to distinguish between runs? I think no.
    # TODO: Change output FASTA file name and include version number.
    #   Or maybe just make the final name an input argument.
    # TODO: Remove unused scripts
    # TODO: Create requirements file to make this script properly reproducible with venv.
    #   Maybe add script to call venv and run script automatically. Maybe to create venv too, if needed.
    # TODO: Make the script idempotent, with proper protections against errors partway through ruining things
    # TODO: Add option to exclude decoys
    # TODO: Add option to skip removing softmasks

    logging.info(f"Finished {SCRIPT_NAME}.")


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
    download_failed = False
    for name, source, target in download_name_source_target_triples:
        try:
            download_file(source, target)
        except Exception as exc:
            logging.error(f"Download of {name} from {source} to {target} has generated an exception: {exc}")
            download_failed = True
        else:
            logging.info(f"Downloaded {name}")
    if download_failed:
        raise ValueError("Download of at least one file has failed")

    with open(config.get_source_list_path(), "w") as f:
        logging.info(f"Writing sources to file: {config.get_source_list_path()}")
        lines = [f"{name}: {source}" for name, source, _ in download_name_source_target_triples]
        f.write("\n".join(lines))


def assert_created_ref_genome_matches_expectations(
        config: Config,
        contig_categorizer: ContigCategorizer,
        contig_name_translator: ContigNameTranslator,
        contig_type_desirabilities: ContigTypeDesirabilities,
) -> None:
    with pysam.Fastafile(config.get_temp_output_fasta_path()) as temp_f:
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
        config.get_temp_output_fasta_path(), config.get_local_compressed_rcrs_fasta_path(), contig_name_translator,
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
        "--output_bucket_dir",
        "-o",
        type=str,
        default=None,
        help=(
            "Optional argument. Path to output bucket dir. Argument should be of the form 'gs://some/kind/of/path'. "
            "This script will not overwrite existing files in the bucket."
        ),
    )

    args = parser.parse_args(sys_args)

    return Config(args.working_dir, args.output_bucket_dir)


if __name__ == "__main__":
    main(parse_args(sys.argv[1:]))
