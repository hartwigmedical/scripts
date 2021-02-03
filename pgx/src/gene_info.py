from copy import deepcopy
from typing import List, Dict, Any, Set

from drug_info import DrugInfo, assert_no_overlap_drug_names
from json_alias import Json
from rs_id_info import RsIdInfo, assert_no_overlap_rs_ids
from util import get_key_to_multiple_values


class GeneInfo(object):
    def __init__(self, alleles: List[Dict[str, Any]], drugs: List[DrugInfo], gene: str, genome_build: str,
                 reference_allele: str, variants: Set[RsIdInfo]) -> None:
        if genome_build != "GRCh37":
            raise ValueError("Exiting, we only support GRCh37, not " + genome_build)

        assert_no_overlap_drug_names(drugs, f"GeneInfo json for {gene}")
        assert_no_overlap_rs_ids(variants, f"GeneInfo json for {gene}")

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
    def variants(self) -> Set[RsIdInfo]:
        return deepcopy(self.__variants)

    @classmethod
    def from_json(cls, data: Json) -> "GeneInfo":
        alleles = data["alleles"]
        drugs = [DrugInfo.from_json(drug_json) for drug_json in data["drugs"]]
        gene = data['gene']
        genome_build = data["genomeBuild"]
        reference_allele = data["referenceAllele"]
        variants = {RsIdInfo.from_json(variant) for variant in data["variants"]}
        return GeneInfo(alleles, drugs, gene, genome_build, reference_allele, variants)


def assert_no_overlap_gene_names(gene_infos: List[GeneInfo], source_name: str) -> None:
    if gene_names_overlap(gene_infos):
        gene_name_to_multiple_infos = get_gene_name_to_multiple_infos(gene_infos)
        raise ValueError(
            ("The {source_name} contains gene summaries with the same gene names but different contents. "
             "Duplicates: {gene_name_to_multiple_infos}").format(
                source_name=source_name,
                gene_name_to_multiple_infos=gene_name_to_multiple_infos
            )
        )


def gene_names_overlap(gene_infos: List[GeneInfo]) -> bool:
    return len({info.gene for info in gene_infos}) != len(gene_infos)


def get_gene_name_to_multiple_infos(gene_infos: List[GeneInfo]) -> Dict[str, List[GeneInfo]]:
    return get_key_to_multiple_values([(info.gene, info) for info in gene_infos])
