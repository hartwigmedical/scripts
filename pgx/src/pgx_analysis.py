from copy import deepcopy
from typing import Dict, Tuple, Set

import pandas as pd

from call_data import Grch37CallData, FullCall, HaplotypeCall
from config.panel import Panel
from grch37_call_translator import Grch37CallTranslator
from haplotype_caller import HaplotypeCaller


class PgxReport(object):
    def __init__(self, gene_to_haplotype_calls: Dict[str, Set[HaplotypeCall]], panel_calls_df: pd.DataFrame) -> None:
        self.__gene_to_haplotype_calls = gene_to_haplotype_calls
        self.__panel_calls_df = panel_calls_df

    def __eq__(self, other: object) -> bool:
        return (
            isinstance(other, PgxReport)
            and self.__gene_to_haplotype_calls == other.__gene_to_haplotype_calls
            and self.__panel_calls_df.equals(other.__panel_calls_df)
        )

    def __repr__(self) -> str:
        return (
            f"PgxAnalysis("
            f"gene_to_haplotype_calls={self.__gene_to_haplotype_calls!r}, "
            f"panel_calls_df=\n{self.__panel_calls_df!r}\n"
            f")"
        )

    def get_gene_to_haplotype_calls(self) -> Dict[str, Set[HaplotypeCall]]:
        return deepcopy(self.__gene_to_haplotype_calls)

    def get_panel_calls_df(self) -> pd.DataFrame:
        # TODO: add ref alleles to data frame?
        return deepcopy(self.__panel_calls_df)


class PgxAnalyser(object):
    CALLS_DATAFRAME_COLUMNS = (
        'gene', 'position_GRCh37', 'ref_GRCh37', 'alt_GRCh37', 'position_GRCh38', 'ref_GRCh38', 'alt_GRCh38',
        'rsid', 'variant_annotation', 'filter'
    )

    @classmethod
    def create_pgx_report(cls, call_data: Grch37CallData, panel: Panel) -> PgxReport:
        full_calls = Grch37CallTranslator.get_full_calls(call_data, panel)

        gene_to_haplotype_calls = HaplotypeCaller.get_gene_to_haplotypes_call(full_calls, panel)
        panel_calls_for_patient_df = cls.__get_panel_calls_for_patient_df(full_calls, panel)

        if pd.isna(panel_calls_for_patient_df).any(axis=None):
            raise ValueError(f"This should never happen: Unhandled NaN values:\n{panel_calls_for_patient_df}")

        return PgxReport(gene_to_haplotype_calls, panel_calls_for_patient_df)

    @classmethod
    def __get_panel_calls_for_patient_df(cls, full_calls: Tuple[FullCall, ...], panel: Panel) -> pd.DataFrame:
        panel_calls_found_in_patient_df = cls.__get_calls_found_in_patient_df(full_calls)
        panel_calls_not_found_in_patient_df = cls.__get_calls_not_found_in_patient_df(full_calls, panel)

        panel_calls_for_patient_df = pd.concat(
            [panel_calls_found_in_patient_df, panel_calls_not_found_in_patient_df], sort=True
        )
        panel_calls_for_patient_df = panel_calls_for_patient_df.sort_values(by='position_GRCh37').reset_index(drop=True)
        panel_calls_for_patient_df = panel_calls_for_patient_df[list(cls.CALLS_DATAFRAME_COLUMNS)]
        return panel_calls_for_patient_df

    @classmethod
    def __get_calls_found_in_patient_df(cls, full_calls: Tuple[FullCall, ...]) -> pd.DataFrame:
        # TODO: add ref alleles to data frame?
        panel_calls_found_in_patient_df = pd.DataFrame(columns=cls.CALLS_DATAFRAME_COLUMNS)
        for full_call in full_calls:
            grch37_alleles = [
                annotated.allele for annotated in full_call.annotated_alleles
                if annotated.is_annotated() and not annotated.is_variant_vs_grch37
            ] + [
                annotated.allele for annotated in full_call.annotated_alleles
                if not annotated.is_annotated() or annotated.is_variant_vs_grch37
            ]
            grch38_alleles = [
                annotated.allele for annotated in full_call.annotated_alleles
                if annotated.is_annotated() and not annotated.is_variant_vs_grch38
            ] + [
                annotated.allele for annotated in full_call.annotated_alleles
                if not annotated.is_annotated() or annotated.is_variant_vs_grch38
            ]

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
            panel_calls_found_in_patient_df = panel_calls_found_in_patient_df.append(new_id, ignore_index=True)
        return panel_calls_found_in_patient_df

    @classmethod
    def __get_calls_not_found_in_patient_df(cls, full_calls: Tuple[FullCall, ...], panel: Panel) -> pd.DataFrame:
        rs_ids_found_in_patient = {rs_id for full_call in full_calls for rs_id in full_call.rs_ids if rs_id != "."}
        grch37_coordinates_covered_by_found_calls = {
            coordinate for full_call in full_calls for coordinate in full_call.get_relevant_grch37_coordinates()
        }

        calls_not_found_in_patient_df = pd.DataFrame(columns=cls.CALLS_DATAFRAME_COLUMNS)
        for gene_info in panel.get_gene_infos():
            for rs_id_info in gene_info.rs_id_infos:
                grch37_coordinates_partially_handled = bool(rs_id_info.get_relevant_grch37_coordinates().intersection(
                    grch37_coordinates_covered_by_found_calls))
                if rs_id_info.rs_id not in rs_ids_found_in_patient and not grch37_coordinates_partially_handled:
                    # Assuming REF/REF relative to GRCh38
                    new_id = {
                        'position_GRCh37': str(rs_id_info.start_coordinate_grch37),
                        'rsid': rs_id_info.rs_id,
                        'ref_GRCh37': rs_id_info.reference_allele_grch38,
                        'alt_GRCh37': rs_id_info.reference_allele_grch38,
                        'variant_annotation': "REF_CALL",
                        'filter': "NO_CALL",
                        'gene': gene_info.gene,
                        'ref_GRCh38': rs_id_info.reference_allele_grch38,
                        'alt_GRCh38': rs_id_info.reference_allele_grch38,
                        'position_GRCh38': str(rs_id_info.start_coordinate_grch38)
                    }
                    calls_not_found_in_patient_df = calls_not_found_in_patient_df.append(new_id, ignore_index=True)
        return calls_not_found_in_patient_df
