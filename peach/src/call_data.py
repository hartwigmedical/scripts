from typing import NamedTuple, Tuple, Optional, Set, FrozenSet

from base.filter import Filter
from base.gene_coordinate import GeneCoordinate
from base.util import get_covered_coordinates


class V37Call(NamedTuple):
    start_coordinate: GeneCoordinate
    ref_allele: str
    alleles: Tuple[str, str]  # The order is (ref, alt) when there is one of each
    gene: str
    rs_ids: Tuple[str, ...]
    variant_annotation: str
    filter: Filter


class V37CallData(NamedTuple):
    calls: FrozenSet[V37Call]


class AnnotatedAllele(object):
    def __init__(self, allele: str, is_variant_vs_v37: bool, is_variant_vs_v38: Optional[bool]) -> None:
        self.__allele = allele
        self.__is_variant_vs_v37 = is_variant_vs_v37
        self.__is_variant_vs_v38 = is_variant_vs_v38  # is None if unknown

    @classmethod
    def from_alleles(cls, allele: str, reference_allele_v37: str,
                     reference_allele_v38: Optional[str]) -> "AnnotatedAllele":
        is_v37_variant = (allele != reference_allele_v37)
        is_v38_variant = (
            (allele != reference_allele_v38) if reference_allele_v38 is not None else None
        )
        return AnnotatedAllele(allele, is_v37_variant, is_v38_variant)

    def __eq__(self, other: object) -> bool:  # pragma: no cover
        return (
                isinstance(other, AnnotatedAllele)
                and self.__allele == other.__allele
                and self.__is_variant_vs_v37 == other.__is_variant_vs_v37
                and self.__is_variant_vs_v38 == other.__is_variant_vs_v38
        )

    def __repr__(self) -> str:  # pragma: no cover
        return (
            f"AnnotatedAllele("
            f"allele={self.__allele!r}, "
            f"is_variant_vs_v37={self.__is_variant_vs_v37!r}, "
            f"is_variant_vs_v38={self.__is_variant_vs_v38!r}, "
            f")"
        )

    @property
    def allele(self) -> str:
        return self.__allele

    @property
    def is_variant_vs_v37(self) -> bool:
        return self.__is_variant_vs_v37

    @property
    def is_variant_vs_v38(self) -> bool:
        if self.__is_variant_vs_v38 is None:
            raise ValueError("Cannot get is_variant_vs_v38 if it is None")
        return self.__is_variant_vs_v38

    def is_annotated_vs_v38(self) -> bool:
        return self.__is_variant_vs_v38 is not None


class FullCall(NamedTuple):
    start_coordinate_v37: GeneCoordinate
    reference_allele_v37: str
    start_coordinate_v38: Optional[GeneCoordinate]  # is None if unknown
    reference_allele_v38: Optional[str]  # is None if unknown
    alleles: Tuple[str, str]
    gene: str
    rs_ids: Tuple[str, ...]
    variant_annotation_v37: str
    filter_v37: Filter
    variant_annotation_v38: str
    filter_v38: Filter

    def get_relevant_v37_coordinates(self) -> Set[GeneCoordinate]:
        return get_covered_coordinates(self.start_coordinate_v37, self.reference_allele_v37)

    def get_annotated_alleles(self) -> Tuple[AnnotatedAllele, AnnotatedAllele]:
        annotated_alleles = self.__annotate_allele(self.alleles[0]), self.__annotate_allele(self.alleles[1])
        self.__assert_alleles_in_expected_order(annotated_alleles)
        return annotated_alleles

    def __annotate_allele(self, allele: str) -> AnnotatedAllele:
        return AnnotatedAllele.from_alleles(allele, self.reference_allele_v37, self.reference_allele_v38)

    @classmethod
    def __assert_alleles_in_expected_order(cls, annotated_alleles: Tuple[AnnotatedAllele, AnnotatedAllele]) -> None:
        alleles_in_unexpected_order = (
                annotated_alleles[0].is_variant_vs_v37
                and not annotated_alleles[1].is_variant_vs_v37
        )
        if alleles_in_unexpected_order:
            error_msg = (f"Alleles are in unexpected order, alt before ref: "
                         f"alleles=({annotated_alleles[0]}, {annotated_alleles[1]})")
            raise ValueError(error_msg)


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

    def __repr__(self) -> str:  # pragma: no cover
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
