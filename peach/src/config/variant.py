from typing import NamedTuple, Dict, List, Collection

from base.json_alias import Json
from base.util import get_key_to_multiple_values


class Variant(NamedTuple):
    rs_id: str
    variant_allele: str

    @classmethod
    def from_json(cls, data: Json) -> "Variant":
        rs_id = str(data["rsid"])
        variant_allele = str(data["altAlleleV38"])
        return Variant(rs_id, variant_allele)


def assert_no_overlap_variant_rs_ids(variants: Collection[Variant], source_name: str) -> None:
    if rs_ids_of_variants_overlap(variants):
        rs_id_to_multiple_variants = get_rs_id_to_multiple_variants(variants)
        error_msg = (
            f"The {source_name} contains variants with the same rs_id but different summaries. "
            f"Duplicates: {rs_id_to_multiple_variants}"
        )
        raise ValueError(error_msg)


def rs_ids_of_variants_overlap(variants: Collection[Variant]) -> bool:
    return len({info.rs_id for info in variants}) != len(variants)


def get_rs_id_to_multiple_variants(variants: Collection[Variant]) -> Dict[str, List[Variant]]:
    return get_key_to_multiple_values([(info.rs_id, info) for info in variants])
