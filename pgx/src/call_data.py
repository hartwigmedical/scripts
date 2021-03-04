from typing import NamedTuple, Tuple, Optional, Set

from base.gene_coordinate import GeneCoordinate
from base.util import get_covered_coordinates


HAPLOTYPE_HOMOZYGOUS_SUFFIX = "_HOM"
HAPLOTYPE_HETEROZYGOUS_SUFFIX = "_HET"


class Grch37Call(NamedTuple):
    start_coordinate: GeneCoordinate
    ref_allele: str
    alleles: Tuple[str, str]  # The order is (ref, alt) when there is one of each
    gene: str
    rs_ids: Tuple[str, ...]
    variant_annotation: str
    filter: str


class Grch37CallData(object):
    def __init__(self, calls: Tuple[Grch37Call, ...]) -> None:
        self.__calls = calls

    def __eq__(self, other: object) -> bool:
        return (
                isinstance(other, Grch37CallData)
                and self.__calls == other.__calls
        )

    def __repr__(self) -> str:
        return (
            f"Grch37CallData("
            f"calls={self.__calls!r}, "
            f")"
        )

    @property
    def calls(self) -> Tuple[Grch37Call, ...]:
        return self.__calls


class AnnotatedAllele(object):
    def __init__(self, allele: str, is_variant_vs_grch37: Optional[bool], is_variant_vs_grch38: Optional[bool]) -> None:
        self.__allele = allele
        self.__is_variant_vs_grch37 = is_variant_vs_grch37  # is None if unknown
        self.__is_variant_vs_grch38 = is_variant_vs_grch38  # is None if unknown

    def __eq__(self, other: object) -> bool:
        return (
                isinstance(other, AnnotatedAllele)
                and self.__allele == other.__allele
                and self.__is_variant_vs_grch37 == other.__is_variant_vs_grch37
                and self.__is_variant_vs_grch38 == other.__is_variant_vs_grch38
        )

    def __repr__(self) -> str:
        return (
            f"AnnotatedAllele2("
            f"allele={self.__allele!r}, "
            f"is_variant_vs_grch37={self.__is_variant_vs_grch37!r}, "
            f"is_variant_vs_grch38={self.__is_variant_vs_grch38!r}, "
            f")"
        )

    @property
    def allele(self) -> str:
        return self.__allele

    @property
    def is_variant_vs_grch37(self) -> bool:
        if self.__is_variant_vs_grch37 is None:
            raise ValueError("Cannot get is_variant_vs_grch37 for unannotated allele")
        return self.__is_variant_vs_grch37

    @property
    def is_variant_vs_grch38(self) -> bool:
        if self.__is_variant_vs_grch38 is None:
            raise ValueError("Cannot get is_variant_vs_grch38 for unannotated allele")
        return self.__is_variant_vs_grch38

    def is_annotated(self) -> bool:
        return self.__is_variant_vs_grch37 is not None and self.__is_variant_vs_grch38 is not None


class FullCall(NamedTuple):
    start_coordinate_grch37: GeneCoordinate
    reference_allele_grch37: str
    start_coordinate_grch38: Optional[GeneCoordinate]  # is None if unknown
    reference_allele_grch38: Optional[str]  # is None if unknown
    annotated_alleles: Tuple[AnnotatedAllele, AnnotatedAllele]  # The order is (grch37_ref, grch37_alt) when there is one of each
    gene: str
    rs_ids: Tuple[str, ...]
    variant_annotation: str
    filter: str

    def get_relevant_grch37_coordinates(self) -> Set[GeneCoordinate]:
        return get_covered_coordinates(self.start_coordinate_grch37, self.reference_allele_grch37)


class HaplotypeCall(object):
    def __init__(self, haplotype_name: str, count: int) -> None:
        if not 1 <= count <= 2:
            error_msg = f"Illegal haplotype count {count} for haplotype {haplotype_name}"
            raise SyntaxError(error_msg)

        self.__haplotype_name = haplotype_name
        self.__count = count

    def __eq__(self, other: object) -> bool:
        return (
                isinstance(other, HaplotypeCall)
                and self.__haplotype_name == other.__haplotype_name
                and self.__count == other.__count
        )

    def __hash__(self) -> int:
        return hash((self.__haplotype_name, self.__count))

    def __str__(self) -> str:
        if self.__count == 2:
            return self.haplotype_name + HAPLOTYPE_HOMOZYGOUS_SUFFIX
        else:
            return self.haplotype_name + HAPLOTYPE_HETEROZYGOUS_SUFFIX

    def __repr__(self) -> str:
        return (
            f"HaplotypeCall("
            f"{self.__haplotype_name!r}, "
            f"{self.__count!r}, "
            f")"
        )

    @property
    def haplotype_name(self) -> str:
        return self.__haplotype_name

    @property
    def count(self) -> int:
        return self.__count
