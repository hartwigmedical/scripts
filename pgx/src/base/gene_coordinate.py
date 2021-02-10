from typing import NamedTuple


class GeneCoordinate(NamedTuple):
    chromosome: str
    position: int

    def get_position_string(self) -> str:
        return f"{self.chromosome}:{self.position}"
