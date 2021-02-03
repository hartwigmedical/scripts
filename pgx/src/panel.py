from collections import defaultdict
from copy import deepcopy
from typing import NamedTuple, Dict, Any, List, Set

from util import get_key_to_multiple_values

Json = Any


class GeneCoordinate(NamedTuple):
    chromosome: int
    position: int

    @classmethod
    def from_json(cls, data: Json) -> "GeneCoordinate":
        chromosome = int(data['chromosome'])
        position = int(data['position'])
        return GeneCoordinate(chromosome, position)

    def get_position_string(self) -> str:
        return f"{self.chromosome}:{self.position}"

    def matches_position_string(self, position_string: str) -> bool:
        return self.get_position_string() == position_string


class RsIdInfo(NamedTuple):
    rs_id: str
    start_coordinate: GeneCoordinate

    @classmethod
    def from_json(cls, data: Json) -> "RsIdInfo":
        rs_id = data['rsid']
        coordinate = GeneCoordinate.from_json(data)
        return RsIdInfo(rs_id, coordinate)


class DrugInfo(NamedTuple):
    name: str
    url_prescription_info: str

    @classmethod
    def from_json(cls, data: Json) -> "DrugInfo":
        return DrugInfo(data["name"], data["url_prescription_info"])


class GeneInfo(object):
    def __init__(self, alleles: List[Dict[str, Any]], drugs: List[DrugInfo], gene: str, genome_build: str,
                 reference_allele: str, variants: List[Dict[str, str]]) -> None:
        if genome_build != "GRCh37":
            raise ValueError("Exiting, we only support GRCh37, not " + genome_build)

        if self.__drug_names_overlap(drugs):
            name_to_multiple_drug_infos = self.__get_drug_name_to_multiple_infos(drugs)
            raise ValueError(
                ("GeneInfo json contains drug summaries with the same drug name but different summaries. "
                 "Duplicates: {name_to_multiple_drug_infos}").format(
                    name_to_multiple_drug_infos=name_to_multiple_drug_infos
                )
            )

        self.__alleles = deepcopy(alleles)
        self.__drugs = deepcopy(drugs)
        self.__gene = gene
        self.__genome_build = genome_build
        self.__reference_allele = reference_allele
        self.__variants = deepcopy(variants)

    def __eq__(self, other: object) -> bool:
        return (
                isinstance(other, GeneInfo)
                and self.__alleles == other.__alleles
                and self.__drugs == other.__drugs
                and self.__gene == other.__gene
                and self.__genome_build == other.__genome_build
                and self.__reference_allele == other.__reference_allele
                and self.__variants == other.__variants
        )

    def __repr__(self) -> str:
        return (
            f"GeneInfo("
            f"alleles={self.__alleles}, "
            f"drugs={self.__drugs}, "
            f"gene={self.__gene}, "
            f"genome_build={self.__genome_build}, "
            f"reference_allele={self.__reference_allele}, "
            f"variants={self.__variants}, "
            f")"
        )

    @property
    def allelles(self) -> List[Dict[str, Any]]:
        return deepcopy(self.__alleles)

    @property
    def drugs(self) -> List[DrugInfo]:
        return deepcopy(self.__drugs)

    @property
    def gene(self) -> str:
        return self.__gene

    @property
    def genome_build(self) -> str:
        return self.__genome_build

    @property
    def reference_allele(self) -> str:
        return self.__reference_allele

    @property
    def variants(self) -> List[Dict[str, str]]:
        return deepcopy(self.__variants)

    @classmethod
    def from_json(cls, data: Json) -> "GeneInfo":
        alleles = data["alleles"]
        drugs = [DrugInfo.from_json(drug_json) for drug_json in data["drugs"]]
        gene = data['gene']
        genome_build = data["genomeBuild"]
        reference_allele = data["referenceAllele"]
        variants = data["variants"]
        return GeneInfo(alleles, drugs, gene, genome_build, reference_allele, variants)

    @staticmethod
    def __drug_names_overlap(drug_infos: List[DrugInfo]) -> bool:
        return len({info.name for info in drug_infos}) != len(drug_infos)

    @staticmethod
    def __get_drug_name_to_multiple_infos(drug_infos: List[DrugInfo]) -> Dict[str, List[DrugInfo]]:
        return get_key_to_multiple_values([(info.name, info) for info in drug_infos])


class Panel(object):
    def __init__(self, gene_infos: List[GeneInfo], rs_id_infos: Set[RsIdInfo]) -> None:
        if self.__rs_ids_overlap(rs_id_infos):
            rs_id_to_multiple_infos = self.__get_rs_id_to_multiple_infos(rs_id_infos)
            raise ValueError(
                ("Panel json contains rs id summaries with the same rs id but different positions. "
                 "Duplicates: {rs_id_to_multiple_infos}").format(
                    rs_id_to_multiple_infos=rs_id_to_multiple_infos
                )
            )
        if self.__gene_names_overlap(gene_infos):
            gene_name_to_multiple_infos = self.__get_gene_name_to_multiple_infos(gene_infos)
            raise ValueError(
                ("Panel json contains gene summaries with the same gene names but different contents. "
                 "Duplicates: {gene_name_to_multiple_infos}").format(
                    gene_name_to_multiple_infos=gene_name_to_multiple_infos
                )
            )
        self.__gene_infos = deepcopy(gene_infos)
        self.__rs_id_infos = deepcopy(rs_id_infos)

    def __eq__(self, other: object) -> bool:
        return (
                isinstance(other, Panel)
                and self.__gene_infos == other.__gene_infos
                and self.__rs_id_infos == other.__rs_id_infos
        )

    def __repr__(self) -> str:
        return (
            f"Panel("
            f"gene_infos={self.__gene_infos}, "
            f"rs_id_infos={self.__rs_id_infos}, "
            f")"
        )

    @classmethod
    def from_json(cls, data: Json) -> "Panel":
        gene_infos: List[GeneInfo] = []
        rs_id_infos = set()
        for gene_info_json in data['genes']:
            gene_info = GeneInfo.from_json(gene_info_json)
            gene_infos.append(gene_info)
            for variant in gene_info.variants:
                rs_id_infos.add(RsIdInfo.from_json(variant))

        return Panel(gene_infos, rs_id_infos)

    def get_gene_infos(self) -> List[GeneInfo]:
        return deepcopy(self.__gene_infos)

    def get_rs_id_infos(self) -> Set[RsIdInfo]:
        return deepcopy(self.__rs_id_infos)

    def contains_rs_id_with_position(self, position_string: str) -> bool:
        for info in self.__rs_id_infos:
            if info.start_coordinate.matches_position_string(position_string):
                return True
        return False

    def get_rs_id_with_position(self, position_string: str) -> str:
        matching_rs_ids = []
        for info in self.__rs_id_infos:
            if info.start_coordinate.matches_position_string(position_string):
                matching_rs_ids.append(info.rs_id)

        if matching_rs_ids and len(matching_rs_ids) == 1:
            return matching_rs_ids.pop()
        elif not matching_rs_ids:
            raise ValueError("No rs ids match position")
        else:
            raise ValueError("Multiple rs ids match position")

    def contains_rs_id(self, rs_id: str) -> bool:
        return rs_id in self.get_rs_ids()

    def get_rs_ids(self) -> Set[str]:
        return {info.rs_id for info in self.__rs_id_infos}

    def is_empty(self) -> bool:
        return not self.__gene_infos and not self.__rs_id_infos

    def get_genes(self) -> List[str]:
        return [info.gene for info in self.__gene_infos]

    @staticmethod
    def __rs_ids_overlap(infos: Set[RsIdInfo]) -> bool:
        return len({info.rs_id for info in infos}) != len(infos)

    @staticmethod
    def __get_rs_id_to_multiple_infos(infos: Set[RsIdInfo]) -> Dict[str, List[RsIdInfo]]:
        return get_key_to_multiple_values([(info.rs_id, info) for info in infos])

    @staticmethod
    def __gene_names_overlap(gene_infos: List[GeneInfo]) -> bool:
        return len({info.gene for info in gene_infos}) != len(gene_infos)

    @staticmethod
    def __get_gene_name_to_multiple_infos(gene_infos: List[GeneInfo]) -> Dict[str, List[GeneInfo]]:
        return get_key_to_multiple_values([(info.gene, info) for info in gene_infos])
