from copy import deepcopy
from typing import Dict, Tuple, Set

import pandas as pd

from call_data import Grch37CallData, FullCall, HaplotypeCall, AnnotatedAllele
from config.panel import Panel
from grch37_call_translator import Grch37CallTranslator
from haplotype_caller import HaplotypeCaller


class PgxAnalysis(object):
    CALLS_DATAFRAME_COLUMNS = (
        'gene', 'position_GRCh37', 'ref_GRCh37', 'alt_GRCh37', 'position_GRCh38', 'ref_GRCh38', 'alt_GRCh38',
        'rsid', 'variant_annotation', 'filter'
    )

    def __init__(self, all_full_calls: Tuple[FullCall, ...],
                 gene_to_haplotype_calls: Dict[str, Set[HaplotypeCall]]) -> None:
        # TODO: maybe full calls in a frozenset instead of tuple, since maybe order shouldn't matter,
        #  and duplicates shouldn't be possible anyway
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

    def get_all_full_calls(self) -> Tuple[FullCall, ...]:
        return self.__all_full_calls

    def get_gene_to_haplotype_calls(self) -> Dict[str, Set[HaplotypeCall]]:
        return deepcopy(self.__gene_to_haplotype_calls)

    def get_panel_calls_df(self) -> pd.DataFrame:
        return self.__get_data_frame_from_full_calls(self.__all_full_calls)

    def __get_data_frame_from_full_calls(cls, full_calls: Tuple[FullCall, ...]) -> pd.DataFrame:
        # TODO: add ref alleles to data frame?
        # TODO: do sorting of data frame here
        data_frame = pd.DataFrame(columns=cls.CALLS_DATAFRAME_COLUMNS)
        for full_call in full_calls:
            annotated_alleles = full_call.get_annotated_alleles()
            grch37_alleles = ([
                annotated.allele for annotated in annotated_alleles if not annotated.is_variant_vs_grch37
            ] + [
                annotated.allele for annotated in annotated_alleles if annotated.is_variant_vs_grch37
            ])
            grch38_alleles = ([
                annotated.allele for annotated in annotated_alleles
                if annotated.is_annotated_vs_grch38() and not annotated.is_variant_vs_grch38
            ] + [
                annotated.allele for annotated in annotated_alleles
                if not annotated.is_annotated_vs_grch38() or annotated.is_variant_vs_grch38
            ])

            position_grch38 = (
                str(full_call.start_coordinate_grch38)
                if full_call.start_coordinate_grch38 is not None
                else "UNKNOWN"
            )

            new_id = {'position_GRCh37': str(full_call.start_coordinate_grch37),
                      'rsid': ";".join(list(full_call.rs_ids)),
                      'ref_GRCh37': grch37_alleles[0],
                      'alt_GRCh37': grch37_alleles[1],
                      'variant_annotation': full_call.variant_annotation,
                      'filter': full_call.filter,
                      'gene': full_call.gene,
                      'position_GRCh38': position_grch38,
                      'ref_GRCh38': grch38_alleles[0],
                      'alt_GRCh38': grch38_alleles[1]}
            data_frame = data_frame.append(new_id, ignore_index=True)

        if pd.isna(data_frame).any(axis=None):
            raise ValueError(f"This should never happen: Unhandled NaN values:\n{data_frame}")

        return data_frame


class PgxAnalyser(object):
    @classmethod
    def create_pgx_analysis(cls, call_data: Grch37CallData, panel: Panel) -> PgxAnalysis:
        full_calls_found_in_patient = Grch37CallTranslator.get_full_calls(call_data, panel)
        all_full_calls = cls.get_all_full_calls(full_calls_found_in_patient, panel)
        gene_to_haplotype_calls = HaplotypeCaller.get_gene_to_haplotypes_call(all_full_calls, panel)
        return PgxAnalysis(all_full_calls, gene_to_haplotype_calls)

    @classmethod
    def get_all_full_calls(
            cls, full_calls_found_in_patient: Tuple[FullCall, ...], panel: Panel) -> Tuple[FullCall, ...]:
        full_calls_not_found_in_patient = cls.__get_full_calls_not_found_in_patient(full_calls_found_in_patient, panel)
        all_full_calls = tuple(sorted(
            list(full_calls_found_in_patient) + list(full_calls_not_found_in_patient),
            key=lambda call: call.start_coordinate_grch37
        ))
        return all_full_calls

    @classmethod
    def __get_full_calls_not_found_in_patient(
            cls, full_calls: Tuple[FullCall, ...], panel: Panel) -> Tuple[FullCall, ...]:
        rs_ids_found_in_patient = {rs_id for full_call in full_calls for rs_id in full_call.rs_ids if rs_id != "."}
        grch37_coordinates_covered_by_found_calls = {
            coordinate for full_call in full_calls for coordinate in full_call.get_relevant_grch37_coordinates()
        }

        grch38_ref_full_calls = []
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
                    grch38_ref_full_calls.append(grch38_ref_full_call)
        return tuple(grch38_ref_full_calls)
