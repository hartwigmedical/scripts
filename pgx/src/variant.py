from typing import NamedTuple, Dict

from json_alias import Json


class Variant(NamedTuple):
    rs_id: str
    alt_allele_grch37: str
    alt_allele_grch38: str

    @classmethod
    def from_json(cls, data: Json) -> "Variant":
        rs_id = data["rsid"]
        alt_allele_grch37 = data["altAllele"]
        alt_allele_grch38 = data["altAlleleGRCh38"]
        return Variant(rs_id, alt_allele_grch37, alt_allele_grch38)
