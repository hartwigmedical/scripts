from copy import deepcopy
from typing import Set, Collection, Dict, List, FrozenSet

from json_alias import Json
from util import get_key_to_multiple_values
from variant import Variant, assert_no_overlap_variant_rs_ids


class Haplotype(object):
    def __init__(self, name: str, function: str, variants: FrozenSet[Variant]) -> None:
        assert_no_overlap_variant_rs_ids(variants, f"haplotype {name}")

        self.__name = name
        self.__function = function
        self.__variants = variants

    def __eq__(self, other: object) -> bool:
        return (
            isinstance(other, Haplotype)
            and self.__name == other.__name
            and self.__function == other.__function
            and self.__variants == other.__variants
        )

    def __repr__(self) -> str:
        return (
            f"Haplotype("
            f"name={self.__name}, "
            f"function={self.__function}, "
            f"variants={self.__variants}, "
            f")"
        )

    @property
    def name(self) -> str:
        return self.__name

    @property
    def function(self) -> str:
        return self.__function

    @property
    def variants(self) -> FrozenSet[Variant]:
        return self.__variants

    @classmethod
    def from_json(cls, data: Json) -> "Haplotype":
        name = data["alleleName"]
        function = data["function"]
        variants = frozenset({Variant.from_json(variant_json) for variant_json in data["alleleVariants"]})
        return Haplotype(name, function, variants)


def assert_no_overlap_haplotype_names(haplotypes: Collection[Haplotype], source_name: str) -> None:
    if names_of_haplotypes_overlap(haplotypes):
        name_to_multiple_haplotypes = get_name_to_multiple_haplotypes(haplotypes)
        raise ValueError(
            ("The {source_name} contains haplotypes with the same name but different summaries. "
             "Duplicates: {name_to_multiple_haplotypes}").format(
                source_name=source_name,
                name_to_multiple_haplotypes=name_to_multiple_haplotypes
            )
        )


def names_of_haplotypes_overlap(haplotypes: Collection[Haplotype]) -> bool:
    return len({haplotype.name for haplotype in haplotypes}) != len(haplotypes)


def get_name_to_multiple_haplotypes(haplotypes: Collection[Haplotype]) -> Dict[str, List[Haplotype]]:
    return get_key_to_multiple_values([(haplotype.name, haplotype) for haplotype in haplotypes])


def assert_no_overlap_haplotype_variant_combinations(haplotypes: Collection[Haplotype], source_name: str) -> None:
    if variant_combinations_of_haplotypes_overlap(haplotypes):
        variant_combination_to_multiple_haplotypes = get_variant_combination_to_multiple_haplotypes(haplotypes)
        raise ValueError(
            ("The {source_name} contains haplotypes with the same variaant combination but different names. "
             "Duplicates: {variant_combination_to_multiple_haplotypes}").format(
                source_name=source_name,
                variant_combination_to_multiple_haplotypes=variant_combination_to_multiple_haplotypes
            )
        )


def variant_combinations_of_haplotypes_overlap(haplotypes: Collection[Haplotype]) -> bool:
    return len({haplotype.variants for haplotype in haplotypes}) != len(haplotypes)


def get_variant_combination_to_multiple_haplotypes(haplotypes: Collection[Haplotype]
                                                   ) -> Dict[FrozenSet[Variant], List[Haplotype]]:
    return get_key_to_multiple_values([(haplotype.variants, haplotype) for haplotype in haplotypes])
