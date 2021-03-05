from copy import deepcopy
from typing import Dict, Set, FrozenSet

import pandas as pd

from call_data import Grch37CallData, FullCall, HaplotypeCall
from config.panel import Panel
from grch37_call_translator import Grch37CallTranslator
from haplotype_caller import HaplotypeCaller


class PgxAnalysis(object):
    def __init__(self, all_full_calls: FrozenSet[FullCall],
                 gene_to_haplotype_calls: Dict[str, Set[HaplotypeCall]]) -> None:
        self.__all_full_calls = all_full_calls
        self.__gene_to_haplotype_calls = gene_to_haplotype_calls

    def __eq__(self, other: object) -> bool:
        return (
            isinstance(other, PgxAnalysis)
            and self.__all_full_calls == other.__all_full_calls
            and self.__gene_to_haplotype_calls == other.__gene_to_haplotype_calls
        )

    def __repr__(self) -> str:
        return (
            f"PgxAnalysis("
            f"all_full_calls={self.__all_full_calls!r}, "
            f"gene_to_haplotype_calls={self.__gene_to_haplotype_calls!r}, "
            f")"
        )

    def get_all_full_calls(self) -> FrozenSet[FullCall]:
        return self.__all_full_calls

    def get_gene_to_haplotype_calls(self) -> Dict[str, Set[HaplotypeCall]]:
        return deepcopy(self.__gene_to_haplotype_calls)


class PgxAnalyser(object):
    @classmethod
    def create_pgx_analysis(cls, call_data: Grch37CallData, panel: Panel) -> PgxAnalysis:
        full_calls_found_in_patient = Grch37CallTranslator.get_full_calls(call_data, panel)
        all_full_calls = cls.get_all_full_calls(full_calls_found_in_patient, panel)
        gene_to_haplotype_calls = HaplotypeCaller.get_gene_to_haplotypes_call(all_full_calls, panel)
        return PgxAnalysis(all_full_calls, gene_to_haplotype_calls)

    @classmethod
    def get_all_full_calls(
            cls, full_calls_found_in_patient: FrozenSet[FullCall], panel: Panel) -> FrozenSet[FullCall]:
        full_calls_not_found_in_patient = cls.__get_full_calls_not_found_in_patient(full_calls_found_in_patient, panel)
        all_full_calls = full_calls_found_in_patient.union(full_calls_not_found_in_patient)
        return all_full_calls

    @classmethod
    def __get_full_calls_not_found_in_patient(
            cls, full_calls: FrozenSet[FullCall], panel: Panel) -> FrozenSet[FullCall]:
        rs_ids_found_in_patient = {rs_id for full_call in full_calls for rs_id in full_call.rs_ids if rs_id != "."}
        grch37_coordinates_covered_by_found_calls = {
            coordinate for full_call in full_calls for coordinate in full_call.get_relevant_grch37_coordinates()
        }

        grch38_ref_full_calls = set()
        for gene_info in panel.get_gene_infos():
            for rs_id_info in gene_info.rs_id_infos:
                grch37_coordinates_partially_handled = bool(
                    rs_id_info.get_relevant_grch37_coordinates().intersection(
                        grch37_coordinates_covered_by_found_calls)
                )
                if rs_id_info.rs_id not in rs_ids_found_in_patient and not grch37_coordinates_partially_handled:
                    # Assuming REF/REF relative to GRCh38

                    # TODO: make strings into constants or similar
                    grch38_ref_full_call = FullCall(
                        rs_id_info.start_coordinate_grch37,
                        rs_id_info.reference_allele_grch37,
                        rs_id_info.start_coordinate_grch38,
                        rs_id_info.reference_allele_grch38,
                        (rs_id_info.reference_allele_grch38, rs_id_info.reference_allele_grch38),
                        gene_info.gene,
                        (rs_id_info.rs_id,),
                        "REF_CALL",
                        "NO_CALL"
                    )
                    grch38_ref_full_calls.add(grch38_ref_full_call)
        return frozenset(grch38_ref_full_calls)
