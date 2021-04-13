import itertools
from typing import List, Set, Tuple, FrozenSet

from base.gene_coordinate import GeneCoordinate
from base.json_alias import Json
from config.gene_info import GeneInfo, assert_no_overlap_gene_names
from config.rs_id_info import RsIdInfo, assert_no_overlap_rs_ids


class Panel(object):
    def __init__(self, name: str, version: str, gene_infos: FrozenSet[GeneInfo]) -> None:
        assert_no_overlap_gene_names(gene_infos, "config json")
        self.__assert_all_rs_id_infos_compatible(gene_infos)

        self.__name = name
        self.__version = version
        self.__gene_infos = gene_infos

    def __eq__(self, other: object) -> bool:
        return (
            isinstance(other, Panel)
            and self.__name == other.__name
            and self.__version == other.__version
            and self.__gene_infos == other.__gene_infos
        )

    def __repr__(self) -> str:  # pragma: no cover
        return (
            f"Panel("
            f"name={self.__name!r}, "
            f"version={self.__version!r}, "
            f"gene_infos={self.__gene_infos!r}, "
            f")"
        )

    @classmethod
    def from_json(cls, data: Json) -> "Panel":
        name = str(data["panelName"])
        version = str(data["panelVersion"])
        gene_infos = frozenset({GeneInfo.from_json(gene_info_json) for gene_info_json in data['genes']})
        return Panel(name, version, gene_infos)

    @property
    def name(self) -> str:
        return self.__name

    @property
    def version(self) -> str:
        return self.__version

    def get_ref_seq_differences(self) -> List[Tuple[RsIdInfo, str, str]]:
        results = []
        for gene_info in self.__gene_infos:
            for rs_id_info in gene_info.rs_id_infos:
                if rs_id_info.reference_allele_v37 != rs_id_info.reference_allele_v38:
                    annotation = gene_info.get_ref_sequence_difference_annotation(rs_id_info.rs_id)
                    results.append((rs_id_info, gene_info.gene, annotation))
        results.sort(key=lambda diff: (diff[1], diff[2]))

        assert_no_overlap_rs_ids([diff[0] for diff in results], "get_ref_seq_differences")
        return results

    def has_ref_seq_difference_annotation(self, gene: str, rs_id: str) -> bool:
        return self.__get_gene_info(gene).has_ref_sequence_difference_annotation(rs_id)

    def get_ref_seq_difference_annotation(self, gene: str, rs_id: str) -> str:
        return self.__get_gene_info(gene).get_ref_sequence_difference_annotation(rs_id)

    def get_gene_infos(self) -> Set[GeneInfo]:
        return set(self.__gene_infos)

    def contains_rs_id_with_v37_coordinate(self, v37_coordinate: GeneCoordinate) -> bool:
        for info in self.__get_rs_id_infos():
            if info.start_coordinate_v37 == v37_coordinate:
                return True
        return False

    def contains_matching_rs_id_info(
            self, v37_coordinate: GeneCoordinate, v37_reference_allele: str) -> bool:
        for info in self.__get_rs_id_infos():
            if (info.start_coordinate_v37 == v37_coordinate
                    and info.reference_allele_v37 == v37_reference_allele):
                return True
        return False

    def get_matching_rs_id_info(
            self, v37_coordinate: GeneCoordinate, v37_reference_allele: str) -> RsIdInfo:
        matching_rs_id_infos = []
        for info in self.__get_rs_id_infos():
            if (info.start_coordinate_v37 == v37_coordinate
                    and info.reference_allele_v37 == v37_reference_allele):
                matching_rs_id_infos.append(info)

        if matching_rs_id_infos and len(matching_rs_id_infos) == 1:
            return matching_rs_id_infos.pop()
        elif not matching_rs_id_infos:
            raise ValueError("No rs id infos match position")
        else:
            raise ValueError("Multiple rs id infos match position")

    def contains_rs_id(self, rs_id: str) -> bool:
        return rs_id in self.__get_rs_ids()

    def is_empty(self) -> bool:
        return not self.__gene_infos

    def get_genes(self) -> Set[str]:
        return {info.gene for info in self.__gene_infos}

    def get_haplotype_function(self, gene: str, haplotype_name: str) -> str:
        gene_info = self.__get_gene_info(gene)
        return gene_info.get_haplotype_function(haplotype_name)

    def __get_gene_info(self, gene: str) -> GeneInfo:
        matching_gene_infos = [gene_info for gene_info in self.__gene_infos if gene_info.gene == gene]
        if len(matching_gene_infos) == 1:
            return matching_gene_infos[0]
        else:
            raise ValueError(f"Not exactly one matching gene info in panel: gene={gene}")

    def __get_rs_ids(self) -> Set[str]:
        return {info.rs_id for info in self.__get_rs_id_infos()}

    def __get_rs_id_infos(self) -> Set[RsIdInfo]:
        return {rs_id_info for gene_info in self.__gene_infos for rs_id_info in gene_info.rs_id_infos}

    @staticmethod
    def __assert_all_rs_id_infos_compatible(gene_infos: FrozenSet[GeneInfo]) -> None:
        for left_gene_info, right_gene_info in itertools.combinations(gene_infos, 2):
            for left_info in left_gene_info.rs_id_infos:
                for right_info in right_gene_info.rs_id_infos:
                    if not left_info.is_compatible(right_info):
                        error_msg = f"Incompatible rs id infos in config. left: {left_info}, right: {right_info}"
                        raise ValueError(error_msg)
