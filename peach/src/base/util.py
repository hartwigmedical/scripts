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


def replace_file_extension_of_path(path: str, new_file_extension: str) -> str:
    assert "." in path, f"Path does not contain point to signal file extension: path={path}"
    split_path = path.split(".")
    new_path = ".".join(split_path[0:-1]) + "." + new_file_extension
    return new_path


def get_chromosome_name_to_index() -> Dict[str, int]:
    chromosome_index_map = {str(i): i for i in range(1, 23)}
    chromosome_index_map["X"] = 23
    chromosome_index_map["Y"] = 24
    return chromosome_index_map
