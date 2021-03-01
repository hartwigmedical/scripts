from typing import NamedTuple, Tuple, Optional

import pandas as pd

from base.gene_coordinate import GeneCoordinate
from dataframe_format import GRCH37_DATAFRAME_COLUMNS


class Grch37Call(NamedTuple):
    start_coordinate: GeneCoordinate
    alleles: Tuple[str, str]  # The order is (ref, alt) when there is one of each
    gene: str
    rs_ids: Tuple[str, ...]
    variant_annotation: str
    filter: str


class Grch37CallData(object):
    def __init__(self, calls: Tuple[Grch37Call, ...]) -> None:
        self.calls = calls

    def __eq__(self, other: object) -> bool:
        return (
                isinstance(other, Grch37CallData)
                and self.calls == other.calls
        )

    def __repr__(self) -> str:
        return (
            f"Grch37CallData("
            f"calls={self.calls!r}, "
            f")"
        )


class AnnotatedAllele(NamedTuple):
    allele: str
    is_variant_vs_grch37: Optional[bool]  # is None if unknown
    is_variant_vs_grch38: Optional[bool]  # is None if unknown


class FullCall(NamedTuple):
    start_coordinate_grch37: GeneCoordinate
    start_coordinate_grch38: Optional[GeneCoordinate]  # is None if unknown
    annotated_alleles: Tuple[AnnotatedAllele, AnnotatedAllele]  # The order should be the same as it was in vcf file
    gene: str
    rs_ids: Tuple[str, ...]
    variant_annotation: str
    filter: str
