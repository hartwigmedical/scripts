from typing import Dict, List, Tuple, Set

import pandas as pd

from call_data import Grch37CallData, FullCall
from config.panel import Panel
from dataframe_format import COMBINED_DATAFRAME_COLUMNS
from grch37_call_translator import Grch37CallTranslator
from haplotype_caller import HaplotypeCaller


class PgxAnalyser(object):
    @classmethod
    def create_pgx_analysis(cls, call_data: Grch37CallData, panel: Panel) -> Tuple[Dict[str, Set[str]], pd.DataFrame]:
        full_calls = Grch37CallTranslator.get_full_calls(call_data, panel)

        gene_to_haplotype_calls = HaplotypeCaller.get_gene_to_haplotypes_call(full_calls, panel)
        panel_calls_for_patient = cls.__get_panel_calls_for_patient(full_calls, panel)

        if pd.isna(panel_calls_for_patient).any(axis=None):
            raise ValueError(f"This should never happen: Unhandled NaN values:\n{panel_calls_for_patient}")

        return gene_to_haplotype_calls, panel_calls_for_patient

    @classmethod
    def __get_panel_calls_for_patient(cls, full_calls: List[FullCall], panel: Panel) -> pd.DataFrame:
        panel_calls_found_in_patient_df = cls.__get_dataframe_from_full_calls(full_calls)
        panel_calls_not_found_in_patient_df_df = cls.__get_calls_not_found_in_patient_df(panel_calls_found_in_patient_df, panel)

        panel_calls_for_patient = pd.concat(
            [panel_calls_found_in_patient_df, panel_calls_not_found_in_patient_df_df], sort=True
        )
        panel_calls_for_patient = panel_calls_for_patient.sort_values(by='position_GRCh37').reset_index(drop=True)
        panel_calls_for_patient = panel_calls_for_patient[list(COMBINED_DATAFRAME_COLUMNS)]
        return panel_calls_for_patient

    @classmethod
    def __get_dataframe_from_full_calls(cls, full_calls: List[FullCall]) -> pd.DataFrame:
        panel_calls_found_in_patient_df = pd.DataFrame(columns=COMBINED_DATAFRAME_COLUMNS)
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
                full_call.start_coordinate_grch38.get_position_string()
                if full_call.start_coordinate_grch38 is not None
                else "UNKNOWN"
            )

            new_id = {'position_GRCh37': full_call.start_coordinate_grch37.get_position_string(),
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
    def __get_calls_not_found_in_patient_df(cls, panel_calls_found_in_patient_df: pd.DataFrame, panel: Panel) -> pd.DataFrame:
        rs_ids_found_in_patient = set(panel_calls_found_in_patient_df.rsid.tolist())
        positions_found_in_patient = set(panel_calls_found_in_patient_df.position_GRCh37.tolist())
        calls_not_found_in_patient_df = pd.DataFrame(columns=COMBINED_DATAFRAME_COLUMNS)
        for gene_info in panel.get_gene_infos():
            for rs_id_info in gene_info.rs_id_infos:
                if (rs_id_info.rs_id not in rs_ids_found_in_patient and
                        rs_id_info.start_coordinate_grch37.get_position_string() not in positions_found_in_patient):
                    # TODO: check whether this properly takes variants of more than one base pair,
                    #       so MNV's etc., into account. The fact that only a single position is checked is suspicious.
                    # Assuming REF/REF relative to GRCh38
                    new_id = {
                        'position_GRCh37': rs_id_info.start_coordinate_grch37.get_position_string(),
                        'rsid': rs_id_info.rs_id,
                        'ref_GRCh37': rs_id_info.reference_allele_grch38,
                        'alt_GRCh37': rs_id_info.reference_allele_grch38,
                        'variant_annotation': "REF_CALL",
                        'filter': "NO_CALL",
                        'gene': gene_info.gene,
                        'ref_GRCh38': rs_id_info.reference_allele_grch38,
                        'alt_GRCh38': rs_id_info.reference_allele_grch38,
                        'position_GRCh38': rs_id_info.start_coordinate_grch38.get_position_string()
                    }
                    calls_not_found_in_patient_df = calls_not_found_in_patient_df.append(new_id, ignore_index=True)
        return calls_not_found_in_patient_df
