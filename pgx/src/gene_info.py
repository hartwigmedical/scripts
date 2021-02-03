from copy import deepcopy
from typing import List, Dict, Set

from drug_info import DrugInfo, assert_no_overlap_drug_names
from haplotype import Haplotype
from json_alias import Json
from rs_id_info import RsIdInfo, assert_no_overlap_rs_ids
from util import get_key_to_multiple_values


class GeneInfo(object):
    def __init__(self, haplotypes: List[Haplotype], drugs: List[DrugInfo], gene: str, genome_build: str,
                 reference_allele: str, rs_id_infos: Set[RsIdInfo]) -> None:
        # TODO: make sure that haplotypes cannot have the same name
        #  and maybe also not the exact same combination of variants
        if genome_build != "GRCh37":
            raise ValueError("Exiting, we only support GRCh37, not " + genome_build)

        assert_no_overlap_drug_names(drugs, f"GeneInfo json for {gene}")
        assert_no_overlap_rs_ids(rs_id_infos, f"GeneInfo json for {gene}")

        self.__haplotypes = deepcopy(haplotypes)
        self.__drugs = deepcopy(drugs)
        self.__gene = gene
        self.__genome_build = genome_build
        self.__reference_allele = reference_allele
        self.__rs_id_infos = deepcopy(rs_id_infos)

    def __eq__(self, other: object) -> bool:
        return (
            isinstance(other, GeneInfo)
            and self.__haplotypes == other.__haplotypes
            and self.__drugs == other.__drugs
            and self.__gene == other.__gene
            and self.__genome_build == other.__genome_build
            and self.__reference_allele == other.__reference_allele
            and self.__rs_id_infos == other.__rs_id_infos
        )

    def __repr__(self) -> str:
        return (
            f"GeneInfo("
            f"haplotypes={self.__haplotypes}, "
            f"drugs={self.__drugs}, "
            f"gene={self.__gene}, "
            f"genome_build={self.__genome_build}, "
            f"reference_allele={self.__reference_allele}, "
            f"rs_id_infos={self.__rs_id_infos}, "
            f")"
        )

    @property
    def haplotypes(self) -> List[Haplotype]:
        return deepcopy(self.__haplotypes)

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
    def rs_id_infos(self) -> Set[RsIdInfo]:
        return deepcopy(self.__rs_id_infos)

    @classmethod
    def from_json(cls, data: Json) -> "GeneInfo":
        alleles = [Haplotype.from_json(haplotype_json) for haplotype_json in data["alleles"]]
        drugs = [DrugInfo.from_json(drug_json) for drug_json in data["drugs"]]
        gene = data['gene']
        genome_build = data["genomeBuild"]
        reference_allele = data["referenceAllele"]
        rs_id_infos = {RsIdInfo.from_json(rs_id_info_json) for rs_id_info_json in data["variants"]}
        return GeneInfo(alleles, drugs, gene, genome_build, reference_allele, rs_id_infos)


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
