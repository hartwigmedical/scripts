import itertools
from copy import deepcopy
from typing import List, Dict, Collection, FrozenSet

from base.constants import NORMAL_FUNCTION_STRING
from base.json_alias import Json
from base.util import get_key_to_multiple_values
from config.drug_info import DrugInfo, assert_no_overlap_drug_names
from config.haplotype import Haplotype, assert_no_overlap_haplotype_names, \
    assert_no_overlap_haplotype_variant_combinations
from config.rs_id_info import RsIdInfo


class GeneInfo(object):
    """This object is meant to be immutable"""
    def __init__(self, gene: str, chromosome: str, wild_type_haplotype_name: str,
                 haplotypes: FrozenSet[Haplotype], rs_id_infos: FrozenSet[RsIdInfo], drugs: FrozenSet[DrugInfo],
                 rs_id_to_ref_seq_difference_annotation_v38: Dict[str, str]) -> None:
        assert_no_overlap_haplotype_names(haplotypes, f"gene info for {gene}")
        assert_no_overlap_haplotype_variant_combinations(haplotypes, f"gene info for {gene}")
        assert_no_overlap_drug_names(drugs, f"GeneInfo json for {gene}")
        self.__assert_rs_id_infos_compatible(rs_id_infos)
        self.__assert_rs_id_infos_match_chromosome(rs_id_infos, chromosome)
        self.__assert_info_exists_for_all_rs_ids_in_haplotypes(haplotypes, rs_id_infos)
        self.__assert_variants_in_haplotypes_compatible_with_rs_id_infos(haplotypes, rs_id_infos)
        self.__assert_rs_ids_with_ref_seq_differences_match_annotations(
            rs_id_infos, rs_id_to_ref_seq_difference_annotation_v38
        )

        self.__gene = gene
        self.__chromosome = chromosome
        self.__wild_type_haplotype_name = wild_type_haplotype_name
        self.__haplotypes = haplotypes
        self.__rs_id_infos = rs_id_infos
        self.__drugs = drugs
        self.__rs_id_to_ref_seq_difference_annotation_v38 = deepcopy(rs_id_to_ref_seq_difference_annotation_v38)

    def __eq__(self, other: object) -> bool:
        return (
                isinstance(other, GeneInfo)
                and self.__gene == other.__gene
                and self.__chromosome == other.__chromosome
                and self.__wild_type_haplotype_name == other.__wild_type_haplotype_name
                and self.__haplotypes == other.__haplotypes
                and self.__rs_id_infos == other.__rs_id_infos
                and self.__drugs == other.__drugs
                and self.__rs_id_to_ref_seq_difference_annotation_v38 == other.__rs_id_to_ref_seq_difference_annotation_v38
        )

    def __hash__(self) -> int:
        return hash((
            self.__gene,
            self.__chromosome,
            self.__wild_type_haplotype_name,
            self.__haplotypes,
            self.__rs_id_infos,
            self.__drugs,
            tuple(sorted(self.__rs_id_to_ref_seq_difference_annotation_v38.items())),
        ))

    def __repr__(self) -> str:  # pragma: no cover
        return (
            f"GeneInfo("
            f"gene={self.__gene!r}, "
            f"chromosome={self.__chromosome!r}, "
            f"wild_type_haplotype_name={self.__wild_type_haplotype_name!r}, "
            f"haplotypes={self.__haplotypes!r}, "
            f"rs_id_infos={self.__rs_id_infos!r}, "
            f"drugs={self.__drugs!r}, "
            f"rs_id_to_ref_seq_difference_annotation={self.__rs_id_to_ref_seq_difference_annotation_v38!r}, "
            f")"
        )

    @property
    def gene(self) -> str:
        return self.__gene

    @property
    def chromosome(self) -> str:
        return self.__chromosome

    @property
    def wild_type_haplotype_name(self) -> str:
        return self.__wild_type_haplotype_name

    @property
    def haplotypes(self) -> FrozenSet[Haplotype]:
        return self.__haplotypes

    @property
    def rs_id_infos(self) -> FrozenSet[RsIdInfo]:
        return self.__rs_id_infos

    @property
    def drugs(self) -> FrozenSet[DrugInfo]:
        return self.__drugs

    @classmethod
    def from_json(cls, data: Json) -> "GeneInfo":
        gene = str(data['gene'])
        chromosome = str(data['chromosome'])
        wild_type_haplotype = str(data["wildTypeHaplotype"])
        rs_id_infos = frozenset({
            RsIdInfo.from_json(rs_id_info_json, chromosome) for rs_id_info_json in data["variants"]
        })
        haplotypes = frozenset({Haplotype.from_json(haplotype_json) for haplotype_json in data["haplotypes"]})
        drugs = frozenset({DrugInfo.from_json(drug_json) for drug_json in data["drugs"]})
        rs_id_to_ref_seq_difference_annotation_v38 = {
            str(annotation_json["rsid"]): str(annotation_json["annotationV38"])
            for annotation_json in data["refSeqDifferenceAnnotations"]
        }
        gene_info = GeneInfo(
            gene,
            chromosome,
            wild_type_haplotype,
            haplotypes,
            rs_id_infos,
            drugs,
            rs_id_to_ref_seq_difference_annotation_v38,
        )
        return gene_info

    def has_ref_sequence_difference_annotation(self, rs_id: str) -> bool:
        return rs_id in self.__rs_id_to_ref_seq_difference_annotation_v38.keys()

    def get_ref_sequence_difference_annotation_v38(self, rs_id: str) -> str:
        return self.__rs_id_to_ref_seq_difference_annotation_v38[rs_id]

    def get_haplotype_function(self, haplotype_name: str) -> str:
        if haplotype_name == self.__wild_type_haplotype_name:
            return NORMAL_FUNCTION_STRING
        else:
            return self.__get_haplotype(haplotype_name).function

    def __get_haplotype(self, haplotype_name: str) -> Haplotype:
        matching_haplotypes = [haplotype for haplotype in self.__haplotypes if haplotype.name == haplotype_name]
        assert len(matching_haplotypes) == 1, (
            f"No unique haplotype with name: {haplotype_name} for gene {self.gene}:\n"
            f"matching_haplotypes={matching_haplotypes}"
        )
        return matching_haplotypes[0]

    @staticmethod
    def __assert_info_exists_for_all_rs_ids_in_haplotypes(
            haplotypes: FrozenSet[Haplotype], rs_id_infos: FrozenSet[RsIdInfo]) -> None:
        rs_ids_in_haplotypes = {variant.rs_id for haplotype in haplotypes for variant in haplotype.variants}
        rs_ids_with_info = {info.rs_id for info in rs_id_infos}
        if not rs_ids_in_haplotypes.issubset(rs_ids_with_info):
            rs_ids_without_info = rs_ids_in_haplotypes.difference(rs_ids_with_info)
            error_msg = f"No info available for some of the rs ids in the haplotypes. Rs ids: {rs_ids_without_info}"
            raise ValueError(error_msg)

    @staticmethod
    def __assert_rs_ids_with_ref_seq_differences_match_annotations(
            rs_id_infos: FrozenSet[RsIdInfo], rs_id_to_ref_seq_difference_annotation: Dict[str, str]) -> None:
        rs_ids_from_infos = {
            info.rs_id for info in rs_id_infos if info.reference_allele_v37 != info.reference_allele_v38
        }
        rs_ids_from_annotation = set(rs_id_to_ref_seq_difference_annotation.keys())
        if rs_ids_from_infos != rs_ids_from_annotation:
            rs_ids_with_only_info = rs_ids_from_infos.difference(rs_ids_from_annotation)
            rs_ids_with_only_annotation = rs_ids_from_annotation.difference(rs_ids_from_infos)
            error_msg = (
                f"Rs ids with differences between v37 and v38 do not match "
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
    def __assert_rs_id_infos_match_chromosome(rs_id_infos: FrozenSet[RsIdInfo], chromosome: str) -> None:
        for info in rs_id_infos:
            if info.start_coordinate_v37.chromosome != chromosome:
                error_msg = (
                    f"Rs id and gene disagree on chromosome, "
                    f"'{info.start_coordinate_v37.chromosome}' vs '{chromosome}'. "
                    f"Rs id info: {info}"
                )
                raise ValueError(error_msg)
            if info.start_coordinate_v38.chromosome != chromosome:
                error_msg = (
                    f"Rs id and gene disagree on chromosome, "
                    f"'{info.start_coordinate_v38.chromosome}' vs '{chromosome}'. "
                    f"Rs id info: {info}"
                )
                raise ValueError(error_msg)

    @staticmethod
    def __assert_variants_in_haplotypes_compatible_with_rs_id_infos(
            haplotypes: FrozenSet[Haplotype], rs_id_infos: FrozenSet[RsIdInfo]) -> None:
        variants = {variant for haplotype in haplotypes for variant in haplotype.variants}
        for variant in variants:
            matching_rs_id_infos = [rs_id_info for rs_id_info in rs_id_infos if rs_id_info.rs_id == variant.rs_id]
            assert len(matching_rs_id_infos) == 1, (
                f"Unexpected number of rs id infos match rs id with variant from haplotype:\n"
                f"variant={variant}, matches={matching_rs_id_infos}"
            )
            if variant.variant_allele == matching_rs_id_infos[0].reference_allele_v38:
                error_msg = (f"Allele of variant matches reference allele from rs id info:\n"
                             f"variant={variant}, rs_id_info={matching_rs_id_infos[0]}")
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
