from collections import defaultdict
from typing import List, Tuple, TypeVar, Dict, Set

from base.gene_coordinate import GeneCoordinate

S = TypeVar("S")
T = TypeVar("T")


def get_key_to_multiple_values(key_value_pairs: List[Tuple[S, T]]) -> Dict[S, List[T]]:
    key_to_values = defaultdict(list)
    for key, value in key_value_pairs:
        key_to_values[key].append(value)
    key_to_overlapping_values = {key: values for key, values in key_to_values.items() if len(values) > 1}
    return key_to_overlapping_values


def get_covered_coordinates(start_coordinate: GeneCoordinate, allele: str) -> Set[GeneCoordinate]:
    return {GeneCoordinate(start_coordinate.chromosome, start_coordinate.position + i) for i in range(len(allele))}
