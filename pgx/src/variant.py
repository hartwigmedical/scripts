from typing import NamedTuple, Dict, List, Collection

from json_alias import Json
from util import get_key_to_multiple_values


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


def assert_no_overlap_variant_rs_ids(variants: Collection[Variant], source_name: str) -> None:
    if rs_ids_of_variants_overlap(variants):
        rs_id_to_multiple_variants = get_rs_id_to_multiple_variants(variants)
        raise ValueError(
            ("The {source_name} contains variants with the same rs_id but different summaries. "
             "Duplicates: {rs_id_to_multiple_variants}").format(
                source_name=source_name,
                rs_id_to_multiple_variants=rs_id_to_multiple_variants
            )
        )


def rs_ids_of_variants_overlap(variants: Collection[Variant]) -> bool:
    return len({info.rs_id for info in variants}) != len(variants)


def get_rs_id_to_multiple_variants(variants: Collection[Variant]) -> Dict[str, List[Variant]]:
    return get_key_to_multiple_values([(info.rs_id, info) for info in variants])
