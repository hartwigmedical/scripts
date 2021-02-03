from typing import NamedTuple, List, Dict, Collection

from json_alias import Json
from util import get_key_to_multiple_values


class DrugInfo(NamedTuple):
    name: str
    url_prescription_info: str

    @classmethod
    def from_json(cls, data: Json) -> "DrugInfo":
        return DrugInfo(data["name"], data["url_prescription_info"])


def assert_no_overlap_drug_names(drugs: Collection[DrugInfo], source_name: str) -> None:
    if drug_names_overlap(drugs):
        name_to_multiple_drug_infos = get_drug_name_to_multiple_infos(drugs)
        raise ValueError(
            ("The {source_name} contains drug summaries with the same drug name but different summaries. "
             "Duplicates: {name_to_multiple_drug_infos}").format(
                source_name=source_name,
                name_to_multiple_drug_infos=name_to_multiple_drug_infos
            )
        )


def drug_names_overlap(drug_infos: Collection[DrugInfo]) -> bool:
    return len({info.name for info in drug_infos}) != len(drug_infos)


def get_drug_name_to_multiple_infos(drug_infos: Collection[DrugInfo]) -> Dict[str, List[DrugInfo]]:
    return get_key_to_multiple_values([(info.name, info) for info in drug_infos])
