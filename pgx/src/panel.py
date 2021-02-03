from collections import defaultdict
from copy import deepcopy
from typing import NamedTuple, Dict, Any, List, Optional, Set


class SNP(NamedTuple):
    rs_id: str
    chromosome: int
    position: int

    @classmethod
    def from_json(cls, data: Dict[str, str]) -> "SNP":
        rs_id = data['rsid']
        chromosome = int(data['chromosome'])
        position = int(data['position'])
        return SNP(rs_id, chromosome, position)

    def get_position_string(self) -> str:
        return f"{self.chromosome}:{self.position}"

    def matches_position_string(self, position_string: str) -> bool:
        return self.get_position_string() == position_string


class GeneInfo(object):
    def __init__(self, alleles: List[Dict[str, Any]], drugs: List[Dict[str, str]], gene: str, genome_build: str,
                 reference_allele: str, variants: List[Dict[str, Any]]) -> None:
        if genome_build != "GRCh37":
            raise ValueError("Exiting, we only support GRCh37, not " + str(genome_build))

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
    def drugs(self) -> List[Dict[str, str]]:
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
    def variants(self) -> List[Dict[str, Any]]:
        return deepcopy(self.__variants)

    @classmethod
    def from_json(cls, data: Dict[str, Any]) -> "GeneInfo":
        alleles = deepcopy(data["alleles"])
        drugs = deepcopy(data["drugs"])
        gene = data['gene']
        genome_build = data["genomeBuild"]
        reference_allele = data["referenceAllele"]
        variants = deepcopy(data["variants"])
        return GeneInfo(alleles, drugs, gene, genome_build, reference_allele, variants)


class Panel(object):
    def __init__(self, gene_infos: List[GeneInfo], snps: Set[SNP]) -> None:
        if self.__rs_ids_overlap(snps):
            rs_id_to_duplicate_snps = self.__get_rs_id_to_duplicate_snps(snps)
            raise ValueError(
                ("Panel json contains snps with the same rs id but different positions. "
                 "Duplicates: {rs_id_to_duplicate_snps}").format(
                    rs_id_to_duplicate_snps=rs_id_to_duplicate_snps
                )
            )
        if self.__gene_names_overlap(gene_infos):
            gene_name_to_duplicate_infos = self.__get_gene_name_to_duplicate_infos(gene_infos)
            raise ValueError(
                ("Panel json contains gene summaries with the same gene names but different contents. "
                 "Duplicates: {gene_name_to_duplicate_infos}").format(
                    gene_name_to_duplicate_infos=gene_name_to_duplicate_infos
                )
            )
        self.__gene_infos = deepcopy(gene_infos)
        self.__snps = deepcopy(snps)

    def __eq__(self, other: object) -> bool:
        return (
                isinstance(other, Panel)
                and self.__gene_infos == other.__gene_infos
                and self.__snps == other.__snps
        )

    def __repr__(self) -> str:
        return (
            f"Panel("
            f"gene_infos={self.__gene_infos}, "
            f"snps={self.__snps}, "
            f")"
        )

    @classmethod
    def from_json(cls, data: Dict[str, Any]) -> "Panel":
        gene_infos: List[GeneInfo] = []
        snps = set()
        for gene_info_json in data['genes']:
            gene_info = GeneInfo.from_json(gene_info_json)
            gene_infos.append(gene_info)
            for variant in gene_info.variants:
                snps.add(SNP.from_json(variant))

        return Panel(gene_infos, snps)

    def get_gene_infos(self) -> List[GeneInfo]:
        return deepcopy(self.__gene_infos)

    def get_snps(self) -> Set[SNP]:
        return deepcopy(self.__snps)

    def contains_snp_with_position(self, position_string: str) -> bool:
        for snp in self.__snps:
            if snp.matches_position_string(position_string):
                return True
        return False

    def get_snp_with_position(self, position_string: str) -> SNP:
        matching_snps = []
        for snp in self.__snps:
            if snp.matches_position_string(position_string):
                matching_snps.append(snp)

        if matching_snps and len(matching_snps) == 1:
            return matching_snps.pop()
        elif not matching_snps:
            raise ValueError("No snps match position")
        else:
            raise ValueError("Multiple snps match position")

    def contains_rs_id(self, rs_id: str) -> bool:
        return rs_id in self.get_rs_ids()

    def get_rs_ids(self) -> Set[str]:
        return {snp.rs_id for snp in self.__snps}

    def is_empty(self) -> bool:
        return not self.__gene_infos and not self.__snps

    def get_genes(self) -> List[str]:
        return [info.gene for info in self.__gene_infos]

    @staticmethod
    def __rs_ids_overlap(snps: Set[SNP]) -> bool:
        return len({snp.rs_id for snp in snps}) != len(snps)

    @staticmethod
    def __get_rs_id_to_duplicate_snps(snps: Set[SNP]) -> Dict[str, List[SNP]]:
        rs_id_to_snps = defaultdict(list)
        for snp in snps:
            rs_id_to_snps[snp.rs_id].append(snp)
        rs_id_to_duplicate_snps = {rs_id: snps for rs_id, snps in rs_id_to_snps.items() if len(snps) > 1}
        return rs_id_to_duplicate_snps

    @staticmethod
    def __gene_names_overlap(gene_infos: List[GeneInfo]) -> bool:
        return len({info.gene for info in gene_infos}) != len(gene_infos)

    @staticmethod
    def __get_gene_name_to_duplicate_infos(gene_infos: List[GeneInfo]) -> Dict[str, List[GeneInfo]]:
        gene_name_to_infos = defaultdict(list)
        for info in gene_infos:
            gene_name_to_infos[info.gene].append(info)
        gene_name_to_duplicate_infos = {gene: info for gene, info in gene_name_to_infos.items() if len(info) > 1}
        return gene_name_to_duplicate_infos


