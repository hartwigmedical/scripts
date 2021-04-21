from copy import deepcopy
from typing import Dict, Set, FrozenSet

from base.constants import REF_CALL_ANNOTATION_STRING
from base.filter import Filter
from call_data import V37CallData, FullCall, HaplotypeCall
from config.panel import Panel
from v37_call_translator import V37CallTranslator
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

    def __repr__(self) -> str:  # pragma: no cover
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
    def create_pgx_analysis(cls, call_data: V37CallData, panel: Panel) -> PgxAnalysis:
        full_calls_found_in_patient_v37 = V37CallTranslator.get_full_calls(call_data, panel)
        all_full_calls = cls.get_all_full_calls(full_calls_found_in_patient_v37, panel)
        gene_to_haplotype_calls = HaplotypeCaller.get_gene_to_haplotypes_call(all_full_calls, panel)
        return PgxAnalysis(all_full_calls, gene_to_haplotype_calls)

    @classmethod
    def get_all_full_calls(
            cls, full_calls_found_in_patient_v37: FrozenSet[FullCall], panel: Panel) -> FrozenSet[FullCall]:
        full_calls_not_found_in_patient_v37 = cls.__get_full_calls_not_found_in_patient_v37(
            full_calls_found_in_patient_v37, panel
        )
        all_full_calls = full_calls_found_in_patient_v37.union(full_calls_not_found_in_patient_v37)
        return all_full_calls

    @classmethod
    def __get_full_calls_not_found_in_patient_v37(
            cls, full_calls: FrozenSet[FullCall], panel: Panel) -> FrozenSet[FullCall]:
        rs_ids_found_in_patient = {rs_id for full_call in full_calls for rs_id in full_call.rs_ids if rs_id != "."}
        v37_coordinates_covered_by_found_calls = {
            coordinate for full_call in full_calls for coordinate in full_call.get_relevant_v37_coordinates()
        }

        v38_ref_full_calls = set()
        for gene_info in panel.get_gene_infos():
            for rs_id_info in gene_info.rs_id_infos:
                v37_coordinates_partially_handled = bool(
                    rs_id_info.get_relevant_v37_coordinates().intersection(
                        v37_coordinates_covered_by_found_calls)
                )
                if rs_id_info.rs_id not in rs_ids_found_in_patient and not v37_coordinates_partially_handled:
                    # Assuming REF/REF relative to v37

                    if rs_id_info.reference_allele_v37 == rs_id_info.reference_allele_v38:
                        annotation_v38 = REF_CALL_ANNOTATION_STRING
                        filter_v38 = Filter.NO_CALL
                    else:
                        annotation_v38 = panel.get_ref_seq_difference_annotation(
                            gene_info.gene, rs_id_info.start_coordinate_v37, rs_id_info.reference_allele_v37)
                        filter_v38 = Filter.INFERRED_PASS

                    v38_ref_full_call = FullCall(
                        rs_id_info.start_coordinate_v37,
                        rs_id_info.reference_allele_v37,
                        rs_id_info.start_coordinate_v38,
                        rs_id_info.reference_allele_v38,
                        (rs_id_info.reference_allele_v37, rs_id_info.reference_allele_v37),
                        gene_info.gene,
                        (rs_id_info.rs_id,),
                        REF_CALL_ANNOTATION_STRING,
                        Filter.NO_CALL,
                        annotation_v38,
                        filter_v38,
                    )
                    v38_ref_full_calls.add(v38_ref_full_call)
        return frozenset(v38_ref_full_calls)
