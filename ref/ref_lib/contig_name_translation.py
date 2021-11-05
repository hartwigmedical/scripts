import logging
import re
from collections import defaultdict
from copy import deepcopy
from typing import Dict, NamedTuple, Optional, Set, DefaultDict, List

from ref_lib.contig_types import ContigType


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
        refseq_accession_number: Optional[str] = split_line[6]
        assembly_unit = split_line[7]
        # sequence_length = split_line[8]
        ucsc_style_name: Optional[str] = split_line[9]

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
    ASSIGNED_MOLECULES_AUTOSOME = {str(i) for i in range(1, 23)}
    ASSIGNED_MOLECULE_X = "X"
    ASSIGNED_MOLECULE_Y = "Y"
    ASSIGNED_MOLECULE_MITOCHONDRION = "MT"
    ASSIGNED_MOLECULE_NA = "na"

    ASSIGNED_MOLECULE_TYPE_CHROMOSOME = "Chromosome"
    ASSIGNED_MOLECULE_TYPE_NA = "na"
    ASSIGNED_MOLECULE_TYPE_MITOCHONDRION = "Mitochondrion"
    ASSIGNED_MOLECULE_TYPE_SEGMENT = "Segment"

    ASSEMBLY_UNIT_PRIMARY = "Primary Assembly"
    ASSEMBLY_UNIT_PATCHES = "PATCHES"
    ASSEMBLY_UNIT_NON_NUCLEAR = "non-nuclear"
    ASSEMBLY_UNITS_ALT = {f"ALT_REF_LOCI_{i}" for i in range(1, 36)}

    CONTIG_ROLE_ASSEMBLED = "assembled-molecule"
    CONTIG_ROLE_UNLOCALIZED = "unlocalized-scaffold"
    CONTIG_ROLE_UNPLACED = "unplaced-scaffold"
    CONTIG_ROLE_FIX_PATCH = "fix-patch"
    CONTIG_ROLE_NOVEL_PATCH = "novel-patch"
    CONTIG_ROLE_ALT = "alt-scaffold"

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
            canonical_name = "chrEBV"
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
            f"CHR_{summary.sequence_name}",
            cls.get_canonical_name(summary)[3:],
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
        matching_contig_types = [contig_type for (contig_type, is_type) in contig_type_to_is_type.items() if is_type]
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
        if (assigned_molecule in cls.ASSIGNED_MOLECULES_AUTOSOME
                or assigned_molecule in {cls.ASSIGNED_MOLECULE_X, cls.ASSIGNED_MOLECULE_Y}):
            return f"chr{assigned_molecule}"
        elif assigned_molecule == cls.ASSIGNED_MOLECULE_MITOCHONDRION:
            return "chrM"
        elif assigned_molecule == cls.ASSIGNED_MOLECULE_NA:
            return "chrUn"
        else:
            raise ValueError(f"Do not recognize chrom name: {assigned_molecule}")

    @classmethod
    def is_autosome(cls, summary: ContigSummary) -> bool:
        result = (
                summary.sequence_role == cls.CONTIG_ROLE_ASSEMBLED
                and summary.assigned_molecule in cls.ASSIGNED_MOLECULES_AUTOSOME
                and summary.assigned_molecule_type == cls.ASSIGNED_MOLECULE_TYPE_CHROMOSOME
                and summary.assembly_unit == cls.ASSEMBLY_UNIT_PRIMARY
        )
        return result

    @classmethod
    def is_x(cls, summary: ContigSummary) -> bool:
        result = (
                summary.sequence_role == cls.CONTIG_ROLE_ASSEMBLED
                and summary.assigned_molecule == cls.ASSIGNED_MOLECULE_X
                and summary.assigned_molecule_type == cls.ASSIGNED_MOLECULE_TYPE_CHROMOSOME
                and summary.assembly_unit == cls.ASSEMBLY_UNIT_PRIMARY
        )
        return result

    @classmethod
    def is_y(cls, summary: ContigSummary) -> bool:
        result = (
                summary.sequence_role == cls.CONTIG_ROLE_ASSEMBLED
                and summary.assigned_molecule == cls.ASSIGNED_MOLECULE_Y
                and summary.assigned_molecule_type == cls.ASSIGNED_MOLECULE_TYPE_CHROMOSOME
                and summary.assembly_unit == cls.ASSEMBLY_UNIT_PRIMARY
        )
        return result

    @classmethod
    def is_ebv(cls, summary: ContigSummary) -> bool:
        result = (
                summary.sequence_role == cls.CONTIG_ROLE_ASSEMBLED
                and summary.assigned_molecule == cls.ASSIGNED_MOLECULE_NA
                and summary.assigned_molecule_type == cls.ASSIGNED_MOLECULE_TYPE_SEGMENT
                and summary.assembly_unit == cls.ASSEMBLY_UNIT_PRIMARY
        )
        return result

    @classmethod
    def is_unlocalized(cls, summary: ContigSummary) -> bool:
        result = (
                summary.sequence_role == cls.CONTIG_ROLE_UNLOCALIZED
                and summary.assigned_molecule != cls.ASSIGNED_MOLECULE_NA
                and summary.assigned_molecule_type == cls.ASSIGNED_MOLECULE_TYPE_CHROMOSOME
                and summary.assembly_unit == cls.ASSEMBLY_UNIT_PRIMARY
        )
        return result

    @classmethod
    def is_unplaced(cls, summary: ContigSummary) -> bool:
        result = (
                summary.sequence_role == cls.CONTIG_ROLE_UNPLACED
                and summary.assigned_molecule == cls.ASSIGNED_MOLECULE_NA
                and summary.assigned_molecule_type == cls.ASSIGNED_MOLECULE_TYPE_NA
                and summary.assembly_unit == cls.ASSEMBLY_UNIT_PRIMARY
                and "decoy" not in summary.sequence_name
        )
        return result

    @classmethod
    def is_fix_patch(cls, summary: ContigSummary) -> bool:
        result = (
                summary.sequence_role == cls.CONTIG_ROLE_FIX_PATCH
                and summary.assigned_molecule != cls.ASSIGNED_MOLECULE_NA
                and summary.assigned_molecule_type == cls.ASSIGNED_MOLECULE_TYPE_CHROMOSOME
                and summary.assembly_unit == cls.ASSEMBLY_UNIT_PATCHES
        )
        return result

    @classmethod
    def is_novel_patch(cls, summary: ContigSummary) -> bool:
        result = (
                summary.sequence_role == cls.CONTIG_ROLE_NOVEL_PATCH
                and summary.assigned_molecule != cls.ASSIGNED_MOLECULE_NA
                and summary.assigned_molecule_type == cls.ASSIGNED_MOLECULE_TYPE_CHROMOSOME
                and summary.assembly_unit == cls.ASSEMBLY_UNIT_PATCHES
        )
        return result

    @classmethod
    def is_alt(cls, summary: ContigSummary) -> bool:
        result = (
                summary.sequence_role == cls.CONTIG_ROLE_ALT
                and summary.assigned_molecule != cls.ASSIGNED_MOLECULE_NA
                and summary.assigned_molecule_type == cls.ASSIGNED_MOLECULE_TYPE_CHROMOSOME
                and summary.assembly_unit in cls.ASSEMBLY_UNITS_ALT
        )
        return result

    @classmethod
    def is_mitochondrion(cls, summary: ContigSummary) -> bool:
        result = (
                summary.sequence_role == cls.CONTIG_ROLE_ASSEMBLED
                and summary.assigned_molecule == cls.ASSIGNED_MOLECULE_MITOCHONDRION
                and summary.assigned_molecule_type == cls.ASSIGNED_MOLECULE_TYPE_MITOCHONDRION
                and summary.assembly_unit == cls.ASSEMBLY_UNIT_NON_NUCLEAR
        )
        return result

    @classmethod
    def is_decoy(cls, summary: ContigSummary) -> bool:
        result = (
                summary.sequence_role == cls.CONTIG_ROLE_UNPLACED
                and summary.assigned_molecule == cls.ASSIGNED_MOLECULE_NA
                and summary.assigned_molecule_type == cls.ASSIGNED_MOLECULE_TYPE_NA
                and summary.assembly_unit == cls.ASSEMBLY_UNIT_PRIMARY
                and bool(re.fullmatch(r"decoy\d{5}", summary.sequence_name))
        )
        return result


class AliasToCanonicalContigNameTextWriter(object):
    @classmethod
    def create_text_from_assembly_reports_text(cls, assembly_reports_text: str) -> str:
        # We don't use the UCSC-style names from the file because not all contigs have such a name in the file,
        # and because alt scaffolds ans novel patches share the '_alt' suffix in this file.
        # We instead use '_novel' for novel patches.
        contig_summaries = cls._get_contig_summaries(assembly_reports_text)
        canonical_contig_name_to_aliases = cls._get_canonical_name_to_aliases(contig_summaries)
        logging.debug(f"canonical_contig_name_to_aliases=\n{canonical_contig_name_to_aliases}")

        cls._assert_no_contradictions(canonical_contig_name_to_aliases)

        sorted_canonical_contig_names = cls._get_sorted_canonical_contig_names(contig_summaries)

        alias_canonical_name_pairs = [
            (alias, canonical_name)
            for canonical_name in sorted_canonical_contig_names
            for alias in sorted(canonical_contig_name_to_aliases[canonical_name])
        ]
        contig_alias_text = "\n".join(
            f"{alias}\t{canonical_name}" for alias, canonical_name in alias_canonical_name_pairs
        )

        logging.debug(f"contig_alias_text=\n{contig_alias_text}")

        return contig_alias_text

    @classmethod
    def _get_contig_summaries(cls, assembly_reports_text: str) -> List[ContigSummary]:
        summaries: List[ContigSummary] = []
        for line in assembly_reports_text.split("\n"):
            if not line or line[0] == "#":
                # skip empty lines and headers starting with "'"#"
                continue
            summaries.append(ContigSummary.from_line(line.replace("\r", "")))
        return summaries

    @classmethod
    def _get_canonical_name_to_aliases(cls, contig_summaries: List[ContigSummary]) -> Dict[str, Set[str]]:
        result: DefaultDict[str, Set[str]] = defaultdict(set)
        for summary in contig_summaries:
            canonical_name = ContigSummaryInterpreter.get_canonical_name(summary)
            aliases = ContigSummaryInterpreter.get_aliases(summary)

            result[canonical_name] = result[canonical_name].union(aliases)

        return dict(result)

    @classmethod
    def _assert_no_contradictions(cls, canonical_contig_name_to_aliases: Dict[str, Set[str]]) -> None:
        seen_aliases: Set[str] = set()
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

    @classmethod
    def _get_sorted_canonical_contig_names(cls, contig_summaries: List[ContigSummary]) -> List[str]:
        contig_type_canonical_name_pairs = {
            (ContigSummaryInterpreter.get_contig_type(summary), ContigSummaryInterpreter.get_canonical_name(summary))
            for summary in contig_summaries
        }
        sorted_contig_type_canonical_name_pairs = sorted(
            contig_type_canonical_name_pairs,
            key=lambda pair: (pair[0].value, cls._get_chromosome_number(pair[1]), pair[1])
        )
        return [pair[1] for pair in sorted_contig_type_canonical_name_pairs]

    @classmethod
    def _get_chromosome_number(cls, canonical_contig_name: str) -> int:
        matches = re.findall(r"^chr\d*", canonical_contig_name)
        if len(matches) != 1:
            raise ValueError("Canonical contig name does not start with chr.")
        match = matches[0]
        if match == "chr":
            return 1000
        else:
            return int(match[3:])
