import itertools
from typing import Set, FrozenSet

from base.gene_coordinate import GeneCoordinate
from base.json_alias import Json
from call_data import V37Call
from config.gene_info import GeneInfo, assert_no_overlap_gene_names
from config.rs_id_info import RsIdInfo


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

    def has_ref_seq_difference_annotation(
            self, gene: str, v37_coordinate: GeneCoordinate, v37_reference_allele: str) -> bool:
        if not self.contains_rs_id_with_v37_coordinate_and_reference_allele(v37_coordinate, v37_reference_allele):
            return False
        rs_id = self.get_matching_rs_id_info(v37_coordinate, v37_reference_allele).rs_id
        return self.__get_gene_info(gene).has_ref_sequence_difference_annotation(rs_id)

    def get_ref_seq_difference_annotation(
            self, gene: str, v37_coordinate: GeneCoordinate, v37_reference_allele: str) -> str:
        rs_id = self.get_matching_rs_id_info(v37_coordinate, v37_reference_allele).rs_id
        return self.__get_gene_info(gene).get_ref_sequence_difference_annotation_v38(rs_id)

    def get_gene_infos(self) -> Set[GeneInfo]:
        return set(self.__gene_infos)

    def contains_rs_id_with_v37_coordinate(self, v37_coordinate: GeneCoordinate) -> bool:
        for info in self.__get_rs_id_infos():
            if info.start_coordinate_v37 == v37_coordinate:
                return True
        return False

    def contains_rs_id_with_v37_coordinate_and_reference_allele(
            self, v37_coordinate: GeneCoordinate, v37_reference_allele: str) -> bool:
        for info in self.__get_rs_id_infos():
            if (info.start_coordinate_v37 == v37_coordinate
                    and info.reference_allele_v37 == v37_reference_allele):
                return True
        return False

    def contains_rs_id_matching_v37_call(self, v37_call: V37Call) -> bool:
        if self.contains_rs_id_with_v37_coordinate_and_reference_allele(v37_call.start_coordinate, v37_call.ref_allele):
            return True
        elif any(self.contains_rs_id(rs_id) for rs_id in v37_call.rs_ids):
            error_msg = (
                f"Match v37 call with rs id info from panel on an rs id but not position:\n"
                f"rs ids: {v37_call.rs_ids}, input file position: {v37_call.start_coordinate}"
            )
            raise ValueError(error_msg)
        else:
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

    def get_id(self) -> str:
        return f"{self.__name}_v{self.__version}"

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
