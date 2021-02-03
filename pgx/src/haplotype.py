from copy import deepcopy
from typing import List, Dict, Any

from json_alias import Json
from variant import Variant


class Haplotype(object):
    def __init__(self, name: str, function: str, variants: List[Variant]) -> None:
        # TODO: make sure variants cannot have the same rs id.
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
    def variants(self) -> List[Variant]:
        return deepcopy(self.__variants)

    @classmethod
    def from_json(cls, data: Json) -> "Haplotype":
        name = data["alleleName"]
        function = data["function"]
        variants = [Variant.from_json(variant_json) for variant_json in data["alleleVariants"]]
        return Haplotype(name, function, variants)
