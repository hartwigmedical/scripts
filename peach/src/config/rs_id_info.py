from typing import NamedTuple, Set

from base.gene_coordinate import GeneCoordinate
from base.json_alias import Json
from base.util import get_covered_coordinates


class RsIdInfo(NamedTuple):
    rs_id: str
    reference_allele_v37: str
    reference_allele_v38: str
    start_coordinate_v37: GeneCoordinate
    start_coordinate_v38: GeneCoordinate

    @classmethod
    def from_json(cls, data: Json, chromosome: str) -> "RsIdInfo":
        rs_id = str(data['rsid'])
        reference_allele_v37 = str(data['referenceAlleleV37'])
        reference_allele_v38 = str(data['referenceAlleleV38'])
        start_coordinate_v37 = GeneCoordinate(chromosome, int(data['positionV37']))
        start_coordinate_v38 = GeneCoordinate(chromosome, int(data['positionV38']))
        info = RsIdInfo(
            rs_id,
            reference_allele_v37,
            reference_allele_v38,
            start_coordinate_v37,
            start_coordinate_v38,
        )
        return info

    def is_compatible(self, other: "RsIdInfo") -> bool:
        if self.rs_id == other.rs_id:
            return self == other
        else:
            return (
                not self.get_relevant_v37_coordinates().intersection(other.get_relevant_v37_coordinates())
                and not self.get_relevant_v38_coordinates().intersection(other.get_relevant_v38_coordinates())
            )

    def get_relevant_v37_coordinates(self) -> Set[GeneCoordinate]:
        return get_covered_coordinates(self.start_coordinate_v37, self.reference_allele_v37)

    def get_relevant_v38_coordinates(self) -> Set[GeneCoordinate]:
        return get_covered_coordinates(self.start_coordinate_v38, self.reference_allele_v38)
