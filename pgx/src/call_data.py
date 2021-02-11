from copy import deepcopy
from typing import List, NamedTuple, Tuple

import pandas as pd

from base.gene_coordinate import GeneCoordinate


class Grch37Call(NamedTuple):
    start_coordinate: GeneCoordinate
    alleles: Tuple[str, str]  # The order is (ref, alt) when there is one of each
    gene: str
    rs_ids: Tuple[str, ...]
    variant_annotation: str
    filter: str


class Grch37CallData(object):
    DATAFRAME_COLUMNS = ['position_GRCh37', 'ref_GRCh37', 'alt_GRCh37', 'rsid', 'variant_annotation', 'gene', 'filter']

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

    def get_data_frame(self) -> pd.DataFrame:
        data_frame = pd.DataFrame(columns=self.DATAFRAME_COLUMNS)
        for call in self.__calls:
            new_id = {
                'position_GRCh37': call.start_coordinate.get_position_string(),
                'ref_GRCh37': call.alleles[0],
                'alt_GRCh37': call.alleles[1],
                'rsid': ";".join(list(call.rs_ids)),
                'variant_annotation': call.variant_annotation,
                'gene': call.gene,
                'filter': call.filter,
            }
            data_frame = data_frame.append(new_id, ignore_index=True)

        return data_frame
