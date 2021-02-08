from copy import deepcopy
from typing import List, Set, FrozenSet, Optional, Dict, Any

from gene_info import GeneInfo, assert_no_overlap_gene_names
from json_alias import Json
from rs_id_info import RsIdInfo, assert_no_overlap_rs_ids


class Panel(object):
    def __init__(self, gene_infos: List[GeneInfo]) -> None:
        assert_no_overlap_gene_names(gene_infos, "panel json")
        self.__assert_genes_do_not_share_rs_ids(gene_infos)
        self.__assert_gene_locations_each_rs_id_info_agree_on_chromosome(gene_infos)

        self.__gene_infos = deepcopy(gene_infos)

    def __eq__(self, other: object) -> bool:
        return (
            isinstance(other, Panel)
            and self.__gene_infos == other.__gene_infos
        )

    def __repr__(self) -> str:
        return (
            f"Panel("
            f"gene_infos={self.__gene_infos!r}, "
            f")"
        )

    @classmethod
    def from_json(cls, data: Json) -> "Panel":
        gene_infos = [GeneInfo.from_json(gene_info_json) for gene_info_json in data['genes']]
        return Panel(gene_infos)

    def get_ref_seq_differences(self) -> List[Json]:
        results = []
        for gene_info in self.__gene_infos:
            for rs_id_info in gene_info.rs_id_infos:
                if rs_id_info.reference_allele_grch37 != rs_id_info.reference_allele_grch38:
                    ref_seq_difference = {
                        'rsid': rs_id_info.rs_id,
                        'gene': gene_info.gene,
                        'referenceAlleleGRCh38': rs_id_info.reference_allele_grch38,
                        'altAlleleGRCh38': rs_id_info.reference_allele_grch37,
                        'chromosome': rs_id_info.start_coordinate_grch38.chromosome,
                        'position': str(rs_id_info.start_coordinate_grch37.position),
                        'positionGRCh38': str(rs_id_info.start_coordinate_grch38.position),
                        'annotationGRCh38': gene_info.get_ref_sequence_difference_annotation(rs_id_info.rs_id)
                    }
                    results.append(ref_seq_difference)
        results.sort(key=lambda diff: (diff["gene"], diff["position"]))
        return results

    def get_gene_infos(self) -> List[GeneInfo]:
        return deepcopy(self.__gene_infos)

    def get_rs_id_infos(self) -> Set[RsIdInfo]:
        return set(self.__get_rs_id_infos_from_gene_infos(self.__gene_infos))

    def contains_rs_id_with_position(self, position_string: str) -> bool:
        for info in self.get_rs_id_infos():
            if info.start_coordinate_grch37.get_position_string() == position_string:
                return True
        return False

    def get_rs_id_with_position(self, position_string: str) -> str:
        matching_rs_ids = []
        for info in self.get_rs_id_infos():
            if info.start_coordinate_grch37.get_position_string() == position_string:
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
        return {info.rs_id for info in self.get_rs_id_infos()}

    def is_empty(self) -> bool:
        return not self.__gene_infos

    def get_genes(self) -> List[str]:
        return [info.gene for info in self.__gene_infos]

    @classmethod
    def __assert_genes_do_not_share_rs_ids(cls, gene_infos: List[GeneInfo]) -> None:
        assert_no_overlap_rs_ids(cls.__get_rs_id_infos_from_gene_infos(gene_infos), "panel json")

    @staticmethod
    def __assert_gene_locations_each_rs_id_info_agree_on_chromosome(gene_infos: List[GeneInfo]) -> None:
        for gene_info in gene_infos:
            for rs_id_info in gene_info.rs_id_infos:
                if rs_id_info.start_coordinate_grch37.chromosome != rs_id_info.start_coordinate_grch38.chromosome:
                    error_msg = (
                        f"Panel only supports rs ids where the GRCh37 and GRCh38 positions are on the same chromosome."
                    )
                    raise ValueError(error_msg)

    @staticmethod
    def __get_rs_id_infos_from_gene_infos(gene_infos: List[GeneInfo]) -> List[RsIdInfo]:
        return [rs_id_info for gene_info in gene_infos for rs_id_info in gene_info.rs_id_infos]

