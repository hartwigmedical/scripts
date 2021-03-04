from typing import NamedTuple, Tuple, Optional

import pandas as pd

from base.gene_coordinate import GeneCoordinate
from dataframe_format import GRCH37_DATAFRAME_COLUMNS


class Grch37Call(NamedTuple):
    start_coordinate: GeneCoordinate
    ref_call: str
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
    start_coordinate_grch38: Optional[GeneCoordinate]  # is None if unknown
    annotated_alleles: Tuple[AnnotatedAllele, AnnotatedAllele]  # The order should be the same as it was in vcf file
    gene: str
    rs_ids: Tuple[str, ...]
    variant_annotation: str
    filter: str
