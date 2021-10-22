import logging
import re
from collections import defaultdict
from copy import deepcopy
from typing import Dict, NamedTuple, Optional, Set, DefaultDict

from contig_types import ContigType


class ContigNameTranslator(object):
    """Standardizes names if it can. Returns argument as is if it cannot."""
    def __init__(self, contig_name_to_canonical_name: Dict[str, str]) -> None:
        self._contig_name_to_canonical_name = deepcopy(contig_name_to_canonical_name)

    @classmethod
    def from_contig_alias_text(cls, text: str) -> "ContigNameTranslator":
        contig_name_to_canonical_name = {}
        for line in text.split("\n"):
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
HARDCODED_CANONICAL_NAME_TO_ALIASES = {"chrEBV": {"chrEBV", "EBV", "AJ507799.2", "NC_007605.1"}}


class ContigSummary(NamedTuple):
    sequence_name: str
    sequence_role: str
    assigned_molecule: str
    assigned_molecule_type: str
    genbank_accession_number: str
    refseq_accession_number: Optional[str]
    assembly_unit: str
    ucsc_style_name: Optional[str]

    @classmethod
    def from_line(cls, line: str) -> "ContigSummary":
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

        contig = ContigSummary(
            sequence_name, sequence_role, assigned_molecule, assigned_molecule_type, genbank_accession_number,
            refseq_accession_number, assembly_unit, ucsc_style_name,
        )
        return contig


class ContigSummaryInterpreter(object):
    @classmethod
    def get_canonical_name(cls, summary: ContigSummary) -> str:
        contig_type = cls.get_contig_type(summary)
        assigned_molecule_proper_name = cls.get_assigned_molecule_proper_name(summary.assigned_molecule)
        genbank_accession_without_dot = summary.genbank_accession_number.replace(".", "v")
        if contig_type == ContigType.AUTOSOME:
            canonical_name = assigned_molecule_proper_name
        elif contig_type == ContigType.X:
            canonical_name = assigned_molecule_proper_name
        elif contig_type == ContigType.Y:
            canonical_name = assigned_molecule_proper_name
        elif contig_type == ContigType.MITOCHONDRIAL:
            canonical_name = assigned_molecule_proper_name
        elif contig_type == ContigType.EBV:
            canonical_name = assigned_molecule_proper_name
        elif contig_type == ContigType.DECOY:
            canonical_name = f"{assigned_molecule_proper_name}_{genbank_accession_without_dot}_decoy"
        elif contig_type == ContigType.UNLOCALIZED:
            canonical_name = f"{assigned_molecule_proper_name}_{genbank_accession_without_dot}_random"
        elif contig_type == ContigType.UNPLACED:
            canonical_name = f"{assigned_molecule_proper_name}_{genbank_accession_without_dot}"
        elif contig_type == ContigType.ALT:
            canonical_name = f"{assigned_molecule_proper_name}_{genbank_accession_without_dot}_alt"
        elif contig_type == ContigType.FIX_PATCH:
            canonical_name = f"{assigned_molecule_proper_name}_{genbank_accession_without_dot}_fix"
        elif contig_type == ContigType.NOVEL_PATCH:
            canonical_name = f"{assigned_molecule_proper_name}_{genbank_accession_without_dot}_novel"
        else:
            raise NotImplementedError(f"Unrecognized contig type: {contig_type}")

        unexplained_canonical_name_mismatch = (
                canonical_name != summary.ucsc_style_name
                and summary.ucsc_style_name is not None
                and not cls.is_novel_patch(summary)
        )
        if unexplained_canonical_name_mismatch:
            error_msg = (
                f"Our canonical name does not match the given UCSC-style name for no clear reason: {summary}"
            )
            raise ValueError(error_msg)

        return canonical_name

    @classmethod
    def get_aliases(cls, summary: ContigSummary) -> Set[str]:
        aliases = {
            cls.get_canonical_name(summary),
            summary.sequence_name,
            summary.genbank_accession_number,
            f"CHR_{summary.sequence_name}"
        }
        if summary.refseq_accession_number is not None:
            aliases.add(summary.refseq_accession_number)
        if summary.ucsc_style_name is not None:
            aliases.add(summary.ucsc_style_name)

        return aliases

    @classmethod
    def get_contig_type(cls, summary: ContigSummary) -> ContigType:
        contig_type_to_is_type = {
            ContigType.AUTOSOME: cls.is_autosome(summary),
            ContigType.X: cls.is_x(summary),
            ContigType.Y: cls.is_y(summary),
            ContigType.MITOCHONDRIAL: cls.is_mitochondrion(summary),
            ContigType.EBV: cls.is_ebv(summary),
            ContigType.DECOY: cls.is_decoy(summary),
            ContigType.UNLOCALIZED: cls.is_unlocalized(summary),
            ContigType.UNPLACED: cls.is_unplaced(summary),
            ContigType.ALT: cls.is_alt(summary),
            ContigType.FIX_PATCH: cls.is_fix_patch(summary),
            ContigType.NOVEL_PATCH: cls.is_novel_patch(summary),
        }
        matching_contig_types = [contig_type for contig_type, is_type in contig_type_to_is_type.values() if is_type]
        if len(matching_contig_types) == 1:
            return matching_contig_types[0]
        else:
            error_msg = (
                f"Did not find a single matching contig type: "
                f"matching_types={matching_contig_types}, summary={summary}"
            )
            raise ValueError(error_msg)

    @classmethod
    def get_assigned_molecule_proper_name(cls, assigned_molecule: str) -> str:
        if assigned_molecule in {str(i) for i in range(1, 23)} or assigned_molecule in {"X", "Y"}:
            return f"chr{assigned_molecule}"
        elif assigned_molecule == "MT":
            return "chrM"
        elif assigned_molecule == "na":
            return "chrUn"
        else:
            raise ValueError(f"Do not recognize chrom name: {assigned_molecule}")

    @classmethod
    def is_autosome(cls, summary: ContigSummary) -> bool:
        result = (
                summary.sequence_role == ASSEMBLED_CONTIG_ROLE
                and summary.assigned_molecule in {str(i) for i in range(1, 23)}
                and summary.assigned_molecule_type == ASSIGNED_MOLECULE_TYPE_CHROMOSOME
                and summary.assembly_unit == PRIMARY_ASSEMBLY_UNIT
        )
        return result

    @classmethod
    def is_x(cls, summary: ContigSummary) -> bool:
        result = (
                summary.sequence_role == ASSEMBLED_CONTIG_ROLE
                and summary.assigned_molecule == "X"
                and summary.assigned_molecule_type == ASSIGNED_MOLECULE_TYPE_CHROMOSOME
                and summary.assembly_unit == PRIMARY_ASSEMBLY_UNIT
        )
        return result

    @classmethod
    def is_y(cls, summary: ContigSummary) -> bool:
        result = (
                summary.sequence_role == ASSEMBLED_CONTIG_ROLE
                and summary.assigned_molecule == "Y"
                and summary.assigned_molecule_type == ASSIGNED_MOLECULE_TYPE_CHROMOSOME
                and summary.assembly_unit == PRIMARY_ASSEMBLY_UNIT
        )
        return result

    @classmethod
    def is_ebv(cls, summary: ContigSummary) -> bool:
        result = (
                summary.sequence_role == ASSEMBLED_CONTIG_ROLE
                and summary.assigned_molecule == "na"
                and summary.assigned_molecule_type == ASSIGNED_MOLECULE_TYPE_CHROMOSOME
                and summary.assembly_unit == PRIMARY_ASSEMBLY_UNIT
        )
        return result

    @classmethod
    def is_unlocalized(cls, summary: ContigSummary) -> bool:
        result = (
                summary.sequence_role == UNLOCALIZED_CONTIG_ROLE
                and summary.assigned_molecule != "na"
                and summary.assigned_molecule_type == ASSIGNED_MOLECULE_TYPE_CHROMOSOME
                and summary.assembly_unit == PRIMARY_ASSEMBLY_UNIT
        )
        return result

    @classmethod
    def is_unplaced(cls, summary: ContigSummary) -> bool:
        result = (
                summary.sequence_role == UNPLACED_CONTIG_ROLE
                and summary.assigned_molecule == "na"
                and summary.assigned_molecule_type == ASSIGNED_MOLECULE_TYPE_UNKNOWN
                and summary.assembly_unit == PRIMARY_ASSEMBLY_UNIT
                and "decoy" not in summary.sequence_name
        )
        return result

    @classmethod
    def is_fix_patch(cls, summary: ContigSummary) -> bool:
        result = (
                summary.sequence_role == FIX_PATCH_CONTIG_ROLE
                and summary.assigned_molecule != "na"
                and summary.assigned_molecule_type == ASSIGNED_MOLECULE_TYPE_CHROMOSOME
                and summary.assembly_unit == PATCHES_ASSEMBLY_UNIT
        )
        return result

    @classmethod
    def is_novel_patch(cls, summary: ContigSummary) -> bool:
        result = (
                summary.sequence_role == NOVEL_PATCH_CONTIG_ROLE
                and summary.assigned_molecule != "na"
                and summary.assigned_molecule_type == ASSIGNED_MOLECULE_TYPE_CHROMOSOME
                and summary.assembly_unit == PATCHES_ASSEMBLY_UNIT
        )
        return result

    @classmethod
    def is_alt(cls, summary: ContigSummary) -> bool:
        result = (
                summary.sequence_role == ALT_CONTIG_ROLE
                and summary.assigned_molecule != "na"
                and summary.assigned_molecule_type == ASSIGNED_MOLECULE_TYPE_CHROMOSOME
                and summary.assembly_unit in ALTS_ASSEMBLY_UNITS
        )
        return result

    @classmethod
    def is_mitochondrion(cls, summary: ContigSummary) -> bool:
        result = (
                summary.sequence_role == ASSEMBLED_CONTIG_ROLE
                and summary.assigned_molecule == "MT"
                and summary.assigned_molecule_type == ASSIGNED_MOLECULE_TYPE_MITOCHONDRION
                and summary.assembly_unit == NON_NUCLEAR_UNIT
        )
        return result

    @classmethod
    def is_decoy(cls, summary: ContigSummary) -> bool:
        result = (
                summary.sequence_role == UNPLACED_CONTIG_ROLE
                and summary.assigned_molecule == "na"
                and summary.assigned_molecule_type == ASSIGNED_MOLECULE_TYPE_UNKNOWN
                and summary.assembly_unit == PRIMARY_ASSEMBLY_UNIT
                and re.fullmatch(r"decoy\d{5}", summary.sequence_name)
        )
        return result


def get_contig_alias_text_from_assembly_reports_text(assembly_reports_text: str) -> str:
    # We don't use the UCSC-style names from the file because not all contigs have such a name in the file,
    # and because alt scaffolds ans novel patches share the '_alt' suffix in this file.
    # We instead use '_novel' for novel patches.
    canonical_contig_name_to_aliases = get_canonical_name_to_aliases(assembly_reports_text)
    logging.debug(f"canonical_contig_name_to_aliases=\n{canonical_contig_name_to_aliases}")
    assert_no_contradictions(canonical_contig_name_to_aliases)
    contig_alias_text = get_text_for_contig_alias_file(canonical_contig_name_to_aliases)
    logging.debug(f"contig_alias_text=\n{contig_alias_text}")
    return contig_alias_text


def get_canonical_name_to_aliases(assembly_reports_text: str) -> Dict[str, Set[str]]:
    result: DefaultDict[str, Set[str]] = defaultdict(set)
    for line in assembly_reports_text.split("\n"):
        if not line or line[0] == "#":
            # skip empty lines and headers starting with "'"#"
            continue

        summary = ContigSummary.from_line(line.replace("\r", ""))

        canonical_name = ContigSummaryInterpreter.get_canonical_name(summary)
        aliases = ContigSummaryInterpreter.get_aliases(summary)

        result[canonical_name] = result[canonical_name].union(aliases)

    for canonical_name, aliases in HARDCODED_CANONICAL_NAME_TO_ALIASES.items():  
        # TODO: remove this when no longer needed
        result[canonical_name] = result[canonical_name].union(aliases)

    return dict(result)


def assert_no_contradictions(canonical_contig_name_to_aliases: Dict[str, Set[str]]) -> None:
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
        error_msg = (
            f"At least one alias present for more than one canonical contig name: "
            f"total_count={total_alias_count}, sum_of_counts={sum_of_alias_counts}"
        )
        raise ValueError(error_msg)


def get_text_for_contig_alias_file(canonical_contig_name_to_aliases: Dict[str, Set[str]]) -> str:
    alias_canonical_name_pairs = [
        (alias, canonical_name)
        for canonical_name, aliases in canonical_contig_name_to_aliases.items()
        for alias in sorted(aliases)
    ]
    contig_alias_text = "\n".join(f"{alias}\t{canonical_name}" for alias, canonical_name in alias_canonical_name_pairs)
    return contig_alias_text
