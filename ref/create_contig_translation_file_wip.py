import argparse
import logging
import re
import sys
from pathlib import Path
from typing import List, NamedTuple, Dict, Optional, Set

from google.cloud import storage

SCRIPT_NAME = "create_chrom_translation_file"
SOURCE_FILE_BUCKET = "hmf-crunch-experiments"
NON_DECOY_TRANSLATION_FILE_PATH = "211005_david_DEV-2170_GRCh38-ref-genome-comparison/GCA_000001405.28_GRCh38.p13_assembly_report.txt"
# original source: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.28_GRCh38.p13/GCA_000001405.28_GRCh38.p13_assembly_report.txt
DECOY_TRANSLATION_FILE_PATH = "211005_david_DEV-2170_GRCh38-ref-genome-comparison/GCA_000786075.2_hs38d1_assembly_report.txt"
# original source: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/786/075/GCA_000786075.2_hs38d1/GCA_000786075.2_hs38d1_assembly_report.txt

ASSEMBLED_CONTIG_ROLE = "assembled-molecule"
UNLOCALIZED_CONTIG_ROLE = "unlocalized-scaffold"
UNPLACED_CONTIG_ROLE = "unplaced-scaffold"
FIX_PATCH_CONTIG_ROLE = "fix-patch"
NOVEL_PATCH_CONTIG_ROLE = "novel-patch"
ALT_CONTIG_ROLE = "alt-scaffold"

ASSIGNED_MOLECULE_TYPE_CHROMOSOME = "Chromosome"
ASSIGNED_MOLECULE_TYPE_UNKNOWN = "na"
ASSIGNED_MOLECULE_TYPE_MITOCHONDRION = "Mitochondrion"

PRIMARY_ASSEMBLY_UNIT = "Primary Assembly"
PATCHES_ASSEMBLY_UNIT = "PATCHES"
NON_NUCLEAR_UNIT = "non-nuclear"
ALTS_ASSEMBLY_UNITS = {f"ALT_REF_LOCI_{i}" for i in range(1, 36)}

UNKNOWN_CHROM = "chrUn"


class Config(NamedTuple):
    output_path: Path

    def validate(self) -> None:
        assert_file_does_not_exist(self.output_path)


class Contig(NamedTuple):
    sequence_name: str
    sequence_role: str
    assigned_molecule: str
    assigned_molecule_type: str
    genbank_accession_number: str
    refseq_accession_number: Optional[str]
    assembly_unit: str
    ucsc_style_name: Optional[str]

    @classmethod
    def from_line(cls, line: str) -> "Contig":
        split_line = line.split("\t")

        if len(split_line) != 10:
            raise ValueError(f"Incorrect line length: {line}")

        sequence_name = split_line[0]
        sequence_role = split_line[1]
        assigned_molecule = split_line[2]
        assigned_molecule_type = split_line[3]
        genbank_accession_number = split_line[4]
        relationship = split_line[5]
        refseq_accession_number = split_line[6]
        assembly_unit = split_line[7]
        # sequence_length = split_line[8]
        ucsc_style_name = split_line[9]

        if relationship == "<>":
            if genbank_accession_number == "na" or refseq_accession_number != "na":
                raise ValueError(
                    f"Can only handle '<>' relationship if RefSeq accession number equals 'na' "
                    f"and Genbank accession number doesn't equal 'na': "
                    f"genbank={genbank_accession_number}, refseq={refseq_accession_number}"
                )
            if refseq_accession_number == "na":
                refseq_accession_number = None
        elif relationship != "=":
            raise ValueError(f"Cannot handle relationship that is not '=' or '<>': {relationship}")

        if ucsc_style_name == "na":
            ucsc_style_name = None

        contig = Contig(
            sequence_name, sequence_role, assigned_molecule, assigned_molecule_type, genbank_accession_number,
            refseq_accession_number, assembly_unit, ucsc_style_name,
        )
        return contig

    def get_canonical_name(self) -> str:
        proper_chrom_name = self.get_proper_chrom_name()
        genbank_accession_without_dot = self.genbank_accession_number.replace(".", "v")
        if self.is_assembled():
            return proper_chrom_name
        elif self.is_unlocalized():
            return f"{proper_chrom_name}_{genbank_accession_without_dot}_random"
        elif self.is_unplaced():
            return f"{proper_chrom_name}_{genbank_accession_without_dot}"
        elif self.is_fix_patch():
            return f"{proper_chrom_name}_{genbank_accession_without_dot}_fix"
        elif self.is_novel_patch():
            return f"{proper_chrom_name}_{genbank_accession_without_dot}_novel"
        elif self.is_alt():
            return f"{proper_chrom_name}_{genbank_accession_without_dot}_alt"
        elif self.is_mitochondrion():
            return proper_chrom_name
        elif self.is_decoy():
            return f"{proper_chrom_name}_{genbank_accession_without_dot}_decoy"
        else:
            raise ValueError(f"Contig has no canonical name: {self}")

    def get_aliases(self) -> Set[str]:
        aliases = {
            self.get_canonical_name(),
            self.sequence_name,
            self.genbank_accession_number,
        }
        if self.refseq_accession_number is not None:
            aliases.add(self.refseq_accession_number)
        if self.ucsc_style_name is not None:
            aliases.add(self.ucsc_style_name)

        return aliases

    def validate(self) -> None:
        is_contig_type_list = [
            self.is_assembled(),
            self.is_unlocalized(),
            self.is_unplaced(),
            self.is_fix_patch(),
            self.is_novel_patch(),
            self.is_alt(),
            self.is_mitochondrion(),
            self.is_decoy(),
        ]
        true_count = len([is_type for is_type in is_contig_type_list if is_type])
        if true_count != 1:
            raise ValueError(f"Contig is not a known type: {self}")

    def is_assembled(self) -> bool:
        result = (
                self.sequence_role == ASSEMBLED_CONTIG_ROLE
                and self.get_proper_chrom_name() != UNKNOWN_CHROM
                and self.assigned_molecule_type == ASSIGNED_MOLECULE_TYPE_CHROMOSOME
                and self.assembly_unit == PRIMARY_ASSEMBLY_UNIT
        )
        return result

    def is_unlocalized(self) -> bool:
        result = (
                self.sequence_role == UNLOCALIZED_CONTIG_ROLE
                and self.get_proper_chrom_name() != UNKNOWN_CHROM
                and self.assigned_molecule_type == ASSIGNED_MOLECULE_TYPE_CHROMOSOME
                and self.assembly_unit == PRIMARY_ASSEMBLY_UNIT
        )
        return result

    def is_unplaced(self) -> bool:
        result = (
                self.sequence_role == UNPLACED_CONTIG_ROLE
                and self.get_proper_chrom_name() == UNKNOWN_CHROM
                and self.assigned_molecule_type == ASSIGNED_MOLECULE_TYPE_UNKNOWN
                and self.assembly_unit == PRIMARY_ASSEMBLY_UNIT
                and "decoy" not in self.sequence_name
        )
        return result

    def is_fix_patch(self) -> bool:
        result = (
                self.sequence_role == FIX_PATCH_CONTIG_ROLE
                and self.get_proper_chrom_name() != UNKNOWN_CHROM
                and self.assigned_molecule_type == ASSIGNED_MOLECULE_TYPE_CHROMOSOME
                and self.assembly_unit == PATCHES_ASSEMBLY_UNIT
        )
        return result

    def is_novel_patch(self) -> bool:
        result = (
                self.sequence_role == NOVEL_PATCH_CONTIG_ROLE
                and self.get_proper_chrom_name() != UNKNOWN_CHROM
                and self.assigned_molecule_type == ASSIGNED_MOLECULE_TYPE_CHROMOSOME
                and self.assembly_unit == PATCHES_ASSEMBLY_UNIT
        )
        return result

    def is_alt(self) -> bool:
        result = (
                self.sequence_role == ALT_CONTIG_ROLE
                and self.get_proper_chrom_name() != UNKNOWN_CHROM
                and self.assigned_molecule_type == ASSIGNED_MOLECULE_TYPE_CHROMOSOME
                and self.assembly_unit in ALTS_ASSEMBLY_UNITS
        )
        return result

    def is_mitochondrion(self) -> bool:
        result = (
                self.sequence_role == ASSEMBLED_CONTIG_ROLE
                and self.get_proper_chrom_name() != UNKNOWN_CHROM
                and self.assigned_molecule_type == ASSIGNED_MOLECULE_TYPE_MITOCHONDRION
                and self.assembly_unit == NON_NUCLEAR_UNIT
        )
        return result

    def is_decoy(self) -> bool:
        result = (
                self.sequence_role == UNPLACED_CONTIG_ROLE
                and self.get_proper_chrom_name() == UNKNOWN_CHROM
                and self.assigned_molecule_type == ASSIGNED_MOLECULE_TYPE_UNKNOWN
                and self.assembly_unit == PRIMARY_ASSEMBLY_UNIT
                and re.fullmatch(r"decoy\d{5}", self.sequence_name)
        )
        return result

    def get_proper_chrom_name(self) -> str:
        if self.assigned_molecule in {str(i) for i in range(1, 23)} or self.assigned_molecule in {"X", "Y"}:
            return f"chr{self.assigned_molecule}"
        elif self.assigned_molecule == "MT":
            return "chrM"
        elif self.assigned_molecule == "na":
            return UNKNOWN_CHROM
        else:
            raise ValueError(f"Do not recognize chrom name: {self.assigned_molecule}")


def main(config: Config) -> None:
    # We don't use the USCS-style names from the file because not all contigs have such a name in the file,
    # and because alt scaffolds ans novel patches share the '_alt' suffix in this file.
    # We instead use '_novel' for novel patches.
    set_up_logging()

    logging.info(f"Starting {SCRIPT_NAME}")

    config.validate()

    contigs = get_contigs()

    for contig in contigs:
        canonical_name_adjustment = contig.get_canonical_name() != contig.ucsc_style_name
        if canonical_name_adjustment and contig.ucsc_style_name is not None and not contig.is_novel_patch():
            error_msg = (
                f"Our canonical name does not match the given USCS-style name for no clear reason: {contig}"
            )
            raise ValueError(error_msg)

    canonical_contig_name_to_aliases = {contig.get_canonical_name(): contig.get_aliases() for contig in contigs}
    logging.info(canonical_contig_name_to_aliases)

    seen_aliases = set()
    for canonical_name, aliases in canonical_contig_name_to_aliases.items():
        if seen_aliases.intersection(aliases):
            error_msg = (
                f"Some alias(es) present for more than one contig: "
                f"contig={canonical_name}, overlap={seen_aliases.intersection(aliases)}"
            )
            raise ValueError(error_msg)
        else:
            seen_aliases = seen_aliases.union(aliases)

    total_alias_count = len(set.union(*canonical_contig_name_to_aliases.values()))
    sum_of_alias_counts = sum(len(aliases) for aliases in canonical_contig_name_to_aliases.values())
    if total_alias_count != sum_of_alias_counts:
        error_msg = f"At least one alias present for more than one canonical contig name: {total_alias_count}, {sum_of_alias_counts}"
        raise ValueError(error_msg)


def get_contigs() -> List[Contig]:
    text = get_translation_text_from_bucket_files()
    # skip empty lines and headers starting with "'"#"
    contigs = [
        Contig.from_line(line.replace("\r", "")) for line in text.split("\n") if line and line[0] != "#"
    ]
    return contigs


def get_translation_text_from_bucket_files() -> str:
    bucket = storage.Client().get_bucket(SOURCE_FILE_BUCKET)
    non_decoy_text = bucket.get_blob(NON_DECOY_TRANSLATION_FILE_PATH).download_as_text()
    decoy_text = bucket.get_blob(DECOY_TRANSLATION_FILE_PATH).download_as_text()
    combined_text = f"{non_decoy_text}\n{decoy_text}"
    return combined_text


def set_up_logging() -> None:
    logging.basicConfig(
        format="%(asctime)s - [%(levelname)-8s] - %(message)s", level=logging.INFO, datefmt="%Y-%m-%d %H:%M:%S"
    )


def assert_file_does_not_exist(path: Path) -> None:
    if path.is_file():
        raise ValueError(f"File exists: {path}")


def parse_args(sys_args: List[str]) -> Config:
    parser = argparse.ArgumentParser(
        prog=f"{SCRIPT_NAME}",
        description=(
            "Create translation file for ref genome fasta file comparison."
        ),
    )
    parser.add_argument(
        "--output",
        "-o",
        type=Path,
        required=True,
        help="Output file path.",
    )

    args = parser.parse_args(sys_args)

    config = Config(Path(args.output))
    return config


if __name__ == "__main__":
    main(parse_args(sys.argv[1:]))
