import itertools
from copy import deepcopy
from typing import List, Dict, Collection, FrozenSet

from base.json_alias import Json
from base.util import get_key_to_multiple_values
from config.drug_info import DrugInfo, assert_no_overlap_drug_names
from config.haplotype import Haplotype, assert_no_overlap_haplotype_names, \
    assert_no_overlap_haplotype_variant_combinations
from config.rs_id_info import RsIdInfo


class GeneInfo(object):
    def __init__(self, gene: str, chromosome: str, genome_build: str, reference_haplotype_name: str,
                 haplotypes: List[Haplotype], rs_id_infos: FrozenSet[RsIdInfo], drugs: List[DrugInfo],
                 rs_id_to_ref_seq_difference_annotation: Dict[str, str]) -> None:
        if genome_build != "GRCh37":
            raise ValueError("Exiting, we only support GRCh37, not " + genome_build)

        assert_no_overlap_haplotype_names(haplotypes, f"gene info for {gene}")
        assert_no_overlap_haplotype_variant_combinations(haplotypes, f"gene info for {gene}")
        assert_no_overlap_drug_names(drugs, f"GeneInfo json for {gene}")
        self.__assert_info_exists_for_all_rs_ids_in_haplotypes(haplotypes, rs_id_infos)
        self.__assert_rs_ids_with_ref_seq_differences_match_annotations(
            rs_id_infos, rs_id_to_ref_seq_difference_annotation
        )
        self.__assert_rs_id_infos_compatible(rs_id_infos)
        self.__assert_rs_id_infos_match_chrosome(rs_id_infos, chromosome)
        self.__assert_gene_locations_each_rs_id_info_agree_on_chromosome(rs_id_infos)

        self.__gene = gene
        self.__chromosome = chromosome
        self.__genome_build = genome_build
        self.__reference_haplotype_name = reference_haplotype_name
        self.__haplotypes = deepcopy(haplotypes)
        self.__rs_id_infos = rs_id_infos
        self.__drugs = deepcopy(drugs)
        self.__rs_id_to_ref_seq_difference_annotation = deepcopy(rs_id_to_ref_seq_difference_annotation)

    def __eq__(self, other: object) -> bool:
        return (
            isinstance(other, GeneInfo)
            and self.__gene == other.__gene
            and self.__chromosome == other.__chromosome
            and self.__genome_build == other.__genome_build
            and self.__reference_haplotype_name == other.__reference_haplotype_name
            and self.__haplotypes == other.__haplotypes
            and self.__rs_id_infos == other.__rs_id_infos
            and self.__drugs == other.__drugs
            and self.__rs_id_to_ref_seq_difference_annotation == other.__rs_id_to_ref_seq_difference_annotation
        )

    def __repr__(self) -> str:
        return (
            f"GeneInfo("
            f"gene={self.__gene!r}, "
            f"chromosome={self.__chromosome!r}, "
            f"genome_build={self.__genome_build!r}, "
            f"reference_haplotype_name={self.__reference_haplotype_name!r}, "
            f"haplotypes={self.__haplotypes!r}, "
            f"rs_id_infos={self.__rs_id_infos!r}, "
            f"drugs={self.__drugs!r}, "
            f"rs_id_to_ref_seq_difference_annotation={self.__rs_id_to_ref_seq_difference_annotation!r}, "
            f")"
        )

    @property
    def gene(self) -> str:
        return self.__gene

    @property
    def chromosome(self) -> str:
        return self.__chromosome

    @property
    def genome_build(self) -> str:
        return self.__genome_build

    @property
    def reference_haplotype_name(self) -> str:
        return self.__reference_haplotype_name

    @property
    def haplotypes(self) -> List[Haplotype]:
        return deepcopy(self.__haplotypes)

    @property
    def rs_id_infos(self) -> FrozenSet[RsIdInfo]:
        return self.__rs_id_infos

    @property
    def drugs(self) -> List[DrugInfo]:
        return deepcopy(self.__drugs)

    @classmethod
    def from_json(cls, data: Json) -> "GeneInfo":
        gene = str(data['gene'])
        chromosome = str(data['chromosome'])
        genome_build = str(data["genomeBuild"])
        reference_allele = str(data["referenceAllele"])
        rs_id_infos = frozenset({RsIdInfo.from_json(rs_id_info_json) for rs_id_info_json in data["variants"]})
        haplotypes = [Haplotype.from_json(haplotype_json) for haplotype_json in data["alleles"]]
        drugs = [DrugInfo.from_json(drug_json) for drug_json in data["drugs"]]
        rs_id_to_ref_seq_difference_annotation = {
            str(annotation_json["rsid"]): str(annotation_json["annotationGRCh38"])
            for annotation_json in data["refSeqDifferenceAnnotations"]
        }
        gene_info = GeneInfo(
            gene,
            chromosome,
            genome_build,
            reference_allele,
            haplotypes,
            rs_id_infos,
            drugs,
            rs_id_to_ref_seq_difference_annotation,
        )
        return gene_info

    def get_ref_sequence_difference_annotation(self, rs_id: str) -> str:
        return self.__rs_id_to_ref_seq_difference_annotation[rs_id]

    @staticmethod
    def __assert_info_exists_for_all_rs_ids_in_haplotypes(
            haplotypes: List[Haplotype], rs_id_infos: FrozenSet[RsIdInfo]) -> None:
        rs_ids_in_haplotypes = {variant.rs_id for haplotype in haplotypes for variant in haplotype.variants}
        rs_ids_with_info = {info.rs_id for info in rs_id_infos}
        if not rs_ids_in_haplotypes.issubset(rs_ids_with_info):
            rs_ids_without_info = rs_ids_in_haplotypes.difference(rs_ids_with_info)
            error_msg = f"No info available for some of the rs ids in the haplotypes. Rs ids: {rs_ids_without_info}"
            raise ValueError(error_msg)

    @staticmethod
    def __assert_rs_ids_with_ref_seq_differences_match_annotations(
            rs_id_infos: FrozenSet[RsIdInfo], rs_id_to_ref_seq_difference_annotation: Dict[str, str]) -> None:
        rs_id_infos_with_ref_seq_difference = {
            info for info in rs_id_infos if info.reference_allele_grch37 != info.reference_allele_grch38
        }
        rs_ids_from_infos = {info.rs_id for info in rs_id_infos_with_ref_seq_difference}
        rs_ids_from_annotation = set(rs_id_to_ref_seq_difference_annotation.keys())
        if rs_ids_from_infos != rs_ids_from_annotation:
            rs_ids_with_only_info = rs_ids_from_infos.difference(rs_ids_from_annotation)
            rs_ids_with_only_annotation = rs_ids_from_annotation.difference(rs_ids_from_infos)
            error_msg = (
                f"Rs ids with differences between GRCh37 and GRCh38 do not match "
                f"the rs ids with annotations for these differences. "
                f"Only info: {rs_ids_with_only_info} "
                f"Only annotation: {rs_ids_with_only_annotation}"
            )
            raise ValueError(error_msg)

    @staticmethod
    def __assert_rs_id_infos_compatible(rs_id_infos: FrozenSet[RsIdInfo]) -> None:
        for left, right in itertools.combinations(rs_id_infos, 2):
            if not left.is_compatible(right):
                error_msg = f"Incompatible rs id infos in gene info. left: {left}, right: {right}"
                raise ValueError(error_msg)

    @staticmethod
    def __assert_rs_id_infos_match_chrosome(rs_id_infos: FrozenSet[RsIdInfo], chromosome: str) -> None:
        for info in rs_id_infos:
            if info.start_coordinate_grch37.chromosome != chromosome:
                error_msg = (
                    f"Rs id and gene disagree on chromosome, "
                    f"'{info.start_coordinate_grch37.chromosome}' vs '{chromosome}'. "
                    f"Rs id info: {info}"
                )
                raise ValueError(error_msg)
            if info.start_coordinate_grch38.chromosome != chromosome:
                error_msg = (
                    f"Rs id and gene disagree on chromosome, "
                    f"'{info.start_coordinate_grch38.chromosome}' vs '{chromosome}'. "
                    f"Rs id info: {info}"
                )
                raise ValueError(error_msg)

    @staticmethod
    def __assert_gene_locations_each_rs_id_info_agree_on_chromosome(rs_id_infos: FrozenSet[RsIdInfo]) -> None:
        for rs_id_info in rs_id_infos:
            if rs_id_info.start_coordinate_grch37.chromosome != rs_id_info.start_coordinate_grch38.chromosome:
                error_msg = (
                    f"Panel only supports rs ids where the GRCh37 and GRCh38 positions are on the same chromosome."
                )
                raise ValueError(error_msg)


def assert_no_overlap_gene_names(gene_infos: Collection[GeneInfo], source_name: str) -> None:
    if gene_names_overlap(gene_infos):
        gene_name_to_multiple_infos = get_gene_name_to_multiple_infos(gene_infos)
        error_msg = (
            f"The {source_name} contains gene summaries with the same gene names but different contents. "
            f"Duplicates: {gene_name_to_multiple_infos}"
        )
        raise ValueError(error_msg)


def gene_names_overlap(gene_infos: Collection[GeneInfo]) -> bool:
    return len({info.gene for info in gene_infos}) != len(gene_infos)


def get_gene_name_to_multiple_infos(gene_infos: Collection[GeneInfo]) -> Dict[str, List[GeneInfo]]:
    return get_key_to_multiple_values([(info.gene, info) for info in gene_infos])
