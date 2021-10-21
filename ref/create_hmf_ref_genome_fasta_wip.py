import argparse
import gzip
import hashlib
import logging
import re
import shutil
import sys
from pathlib import Path
from typing import List, NamedTuple, FrozenSet, Set

import pysam
import requests

from ref_util import set_up_logging, ContigNameTranslator, assert_file_exists_in_bucket, get_blob, \
    assert_file_does_not_exist, assert_dir_does_not_exist, ContigCategorizer, ContigType, Assembly

# See gs://hmf-crunch-experiments/211005_david_DEV-2170_GRCh38-ref-genome-comparison/ for required files.

SCRIPT_NAME = "create_hmf_ref_genome_fasta"

REFSEQ_FASTA_SOURCE = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.fna.gz"
DECOY_FASTA_SOURCE = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/786/075/GCA_000786075.2_hs38d1/GCA_000786075.2_hs38d1_genomic.fna.gz"
EBV_FASTA_SOURCE = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/402/265/GCF_002402265.1_ASM240226v1/GCF_002402265.1_ASM240226v1_genomic.fna.gz"

REFSEQ_FASTA_FILE_NAME = "REFSEQ.fasta"
DECOY_FASTA_FILE_NAME = "DECOY.fasta"
EBV_FASTA_FILE_NAME = "EBV.fasta"
MASTER_FASTA_FILE_NAME = "MASTER.fasta"
SOURCES_LIST_FILE_NAME = "sources.txt"

FASTA_LINE_WRAP = 70
FASTA_HEADER_SEPARATOR = "  "

CONTIG_TYPE_TO_ROLE = {
    ContigType.AUTOSOME: "Chromosome",
    ContigType.X: "Chromosome",
    ContigType.Y: "Chromosome",
    ContigType.MITOCHONDRIAL: "Mitochondrion",
    ContigType.EBV: "decoy",
    ContigType.DECOY: "decoy",
    ContigType.UNLOCALIZED: "unlocalized",
    ContigType.UNPLACED: "unplaced",
    ContigType.ALT: "alt-scaffold",
    ContigType.FIX_PATCH: "fix-patch",
    ContigType.NOVEL_PATCH: "novel-patch",
}


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

    def get_temp_output_path(self) -> Path:
        return self.output_path.parent / f"{self.output_path.name}.tmp"


class ContigDesirabilityClassification(NamedTuple):
    desired_contig_types: FrozenSet[ContigType]
    undesired_contig_types: FrozenSet[ContigType]
    unexpected_contig_types: FrozenSet[ContigType]

    @classmethod
    def create(cls) -> "ContigDesirabilityClassification":
        desired_contig_types = frozenset({
            ContigType.AUTOSOME,
            ContigType.X,
            ContigType.Y,
            ContigType.MITOCHONDRIAL,
            ContigType.EBV,
            ContigType.DECOY,
            ContigType.UNLOCALIZED,
            ContigType.UNPLACED,
        })
        undesired_contig_types = frozenset({
            ContigType.ALT,
            ContigType.FIX_PATCH,
            ContigType.NOVEL_PATCH,
        })
        unexpected_contig_types = frozenset({
            ContigType.UNCATEGORIZABLE,
        })
        result = ContigDesirabilityClassification(desired_contig_types, undesired_contig_types, unexpected_contig_types)
        result.validate()
        return result

    def validate(self) -> None:
        handleable_contig_types = self.desired_contig_types.union(
            self.undesired_contig_types.union(self.unexpected_contig_types)
        )
        if set(contig_type for contig_type in ContigType) != handleable_contig_types:
            raise ValueError(f"Not all contig types can be handled: handleable_contig_types={handleable_contig_types}")
        contig_types_categories_overlap = (
                self.desired_contig_types.intersection(self.undesired_contig_types)
                or self.desired_contig_types.intersection(self.unexpected_contig_types)
                or self.undesired_contig_types.intersection(self.unexpected_contig_types)
        )
        if contig_types_categories_overlap:
            raise ValueError(f"Contig type categories overlap (so desired, undesired and unexpected are not disjoint)")

    def get_expected_contig_types(self) -> Set[ContigType]:
        return set(self.desired_contig_types.union(self.undesired_contig_types))


def main(config: Config) -> None:
    set_up_logging()

    logging.info(f"Starting {SCRIPT_NAME}")

    config.validate()

    config.working_dir.mkdir(parents=True, exist_ok=True)

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
    #
    # logging.info(f"Combining FASTA files into a single file.")
    # combine_files(
    #     [config.get_refseq_fasta_path(), config.get_decoy_fasta_path(), config.get_ebv_fasta_path()],
    #     config.get_master_fasta_path(),  # TODO: reenable
    # )

    logging.info("Categorize contigs in master FASTA file")
    contig_categorizer = ContigCategorizer(contig_name_translator)

    with pysam.Fastafile(config.get_master_fasta_path()) as master_f:
        contig_names = list(master_f.references)

    contig_name_to_type = {name: contig_categorizer.get_contig_type(name) for name in contig_names}

    classification = ContigDesirabilityClassification.create()

    assert_seen_contig_types_match_expected(set(contig_name_to_type.values()), classification)

    logging.info(f"Creating temp version of output FASTA file: {config.get_temp_output_path()}.")
    with pysam.Fastafile(config.get_master_fasta_path()) as master_f:
        with open(config.get_temp_output_path(), "w") as out_f:
            for contig_name in contig_names:
                logging.info(f"Handling {contig_name}")
                contig_type = contig_name_to_type[contig_name]
                if contig_type in classification.desired_contig_types:
                    logging.info(f"Include {contig_name} in output file")
                    sequence = master_f.fetch(contig_name).upper()
                    header = get_header(
                        contig_name,
                        contig_name_translator.standardize(contig_name),
                        contig_type,
                        sequence,
                    )
                    logging.info(f"Header: {header}")
                    out_f.write(header + "\n")
                    for i in range(0, len(sequence), FASTA_LINE_WRAP):
                        out_f.write(sequence[i:i + FASTA_LINE_WRAP] + "\n")

    logging.info("Asserting that output is as expected")
    with pysam.Fastafile(config.get_temp_output_path()) as temp_f:
        contigs_expected_to_be_copied = {
            contig_name for contig_name, contig_type in contig_name_to_type.items()
            if contig_type in classification.desired_contig_types
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
        with pysam.Fastafile(config.get_master_fasta_path()) as master_f:
            for contig_name in contigs_expected_to_be_copied:
                expected_sequence = master_f.fetch(contig_name).upper()
                actual_sequence = temp_f.fetch(contig_name_translator.standardize(contig_name))
                if actual_sequence != expected_sequence:
                    raise ValueError(f"Contig sequences are not identical: contig={contig_name}")

    logging.info("Moving temp output file to real output file")
    config.get_temp_output_path().rename(config.output_path)
    Path(f"{config.get_temp_output_path()}.fai").unlink(missing_ok=True)
    # TODO: Extract contig names from master file, and figure out which contigs we want
    # TODO: Put these contigs into output file, including proper header lines for each contig,
    #   and filtering out decoys and removing soft-masks if configured to do this
    #       A first include decoys and remove soft-masks. Can always add options to not do these things later
    #       TODO: Figure out what these header lines should look like
    # TODO: Add automatic checks to see whether the file is correct and as expected.
    # TODO: move tmp output file to actual output file location
    # TODO: Put used files and result in bucket? Or just do this manually, afterwards.
    # TODO: Create requirements file to make this script properly reproducible with venv.
    #       Maybe add script to call venv and run script automatically. Maybe to create venv too, if needed.

    logging.info(f"Finished {SCRIPT_NAME}")


def get_header(contig_name: str, standardized_contig_name: str, contig_type: ContigType, sequence: str) -> str:
    header_entries = [f">{standardized_contig_name}", f"AC:{contig_name}", f"LN:{len(sequence)}"]
    if contig_type == ContigType.UNLOCALIZED:
        # TODO: Add rg field entry for ALT and patch contigs, if we want these in ref genome.
        #       Would need the coordinates of where they align.
        #       See https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GRCh38_major_release_seqs_for_alignment_pipelines/README_analysis_sets.txt
        region = standardized_contig_name.split("_")[0]
        header_entries.append(f"rg:{region}")
    header_entries.append(f"rl:{get_contig_role(contig_type)}")
    header_entries.append(f"M5:{hashlib.md5(sequence.encode('utf-8')).hexdigest()}")
    assembly = contig_type.get_assembly()
    if assembly == Assembly.GRCH38:
        header_entries.append("AS:GRCh38")
    elif assembly == Assembly.HS38D1:
        header_entries.append("AS:hs38d1")
    elif assembly == Assembly.OTHER:
        pass
    else:
        raise ValueError(f"Unrecognized assembly: {assembly}")
    if contig_type == ContigType.EBV:
        header_entries.append("SP:Human_herpesvirus_4")
    if contig_type in {ContigType.MITOCHONDRIAL, ContigType.EBV}:
        header_entries.append("tp:circular")
    return FASTA_HEADER_SEPARATOR.join(header_entries)


def get_contig_role(contig_type: ContigType) -> str:
    if contig_type in CONTIG_TYPE_TO_ROLE.keys():
        return CONTIG_TYPE_TO_ROLE[contig_type]
    else:
        raise ValueError(f"Encountered contig type without an assigned role: {contig_type}")


def assert_seen_contig_types_match_expected(
    seen_contig_types: Set[ContigType], classification: ContigDesirabilityClassification
) -> None:
    expected_contig_types = classification.get_expected_contig_types()
    if seen_contig_types != expected_contig_types:
        sorted_seen_contig_names = sorted(contig_type.name for contig_type in seen_contig_types)
        sorted_expected_contig_names = sorted(
            contig_type.name for contig_type in expected_contig_types
        )
        error_msg = (
            f"Seen contig types don't perfectly match the expected contig types: "
            f"seen={sorted_seen_contig_names}, expected={sorted_expected_contig_names}"
        )
        raise ValueError(error_msg)


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
