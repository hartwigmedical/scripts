from typing import NamedTuple

from json_alias import Json


class GeneCoordinate(NamedTuple):
    chromosome: int
    position: int

    def get_position_string(self) -> str:
        return f"{self.chromosome}:{self.position}"

    def matches_position_string(self, position_string: str) -> bool:
        return self.get_position_string() == position_string
