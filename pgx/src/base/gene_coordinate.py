from typing import NamedTuple


class GeneCoordinate(NamedTuple):
    chromosome: str
    position: int

    def __str__(self) -> str:
        return f"{self.chromosome}:{self.position}"
