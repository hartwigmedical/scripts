from typing import FrozenSet

import pandas as pd

from base.gene_coordinate import GeneCoordinate
from call_data import FullCall, HaplotypeCall
from config.panel import Panel
from pgx_analysis import PgxAnalysis


class PgxReporter(object):
    CALLS_TSV_COLUMNS = (
        "gene", "position_GRCh37", "ref_GRCh37", "alt_GRCh37", "position_GRCh38", "ref_GRCh38", "alt_GRCh38",
        "rsid", "variant_annotation", "filter"
    )
    GENOTYPE_TSV_COLUMNS = (
        "gene", "haplotype", "function", "linked_drugs", "url_prescription_info", "panel_version", "repo_version")
    TSV_SEPARATOR = "\t"
    DRUG_SEPARATOR = ";"
    RS_ID_SEPARATOR = ";"

    @classmethod
    def get_calls_tsv_text(cls, pgx_analysis: PgxAnalysis) -> str:
        return str(cls.__get_panel_calls_df(pgx_analysis).to_csv(sep=cls.TSV_SEPARATOR, index=False))

    @classmethod
    def get_genotype_tsv_text(cls, pgx_analysis: PgxAnalysis, panel: Panel, panel_path: str, version: str) -> str:
        gene_to_drug_info = {}
        for gene_info in panel.get_gene_infos():
            sorted_drugs = sorted(
                [drug for drug in gene_info.drugs],
                key=lambda info: (info.name, info.url_prescription_info)
            )
            gene_to_drug_info[gene_info.gene] = (
                cls.DRUG_SEPARATOR.join([drug.name for drug in sorted_drugs]),
                cls.DRUG_SEPARATOR.join([drug.url_prescription_info for drug in sorted_drugs])
            )

        gene_to_haplotype_calls = pgx_analysis.get_gene_to_haplotype_calls()

        header = cls.TSV_SEPARATOR.join(cls.GENOTYPE_TSV_COLUMNS)
        lines = [header]
        for gene in sorted(list(gene_to_haplotype_calls.keys())):
            if gene_to_haplotype_calls[gene]:
                for haplotype_call in sorted(list(gene_to_haplotype_calls[gene]), key=lambda call: call.haplotype_name):
                    lines.append(cls.TSV_SEPARATOR.join([
                        gene,
                        PgxReporter.__get_haplotype_call_string(haplotype_call),
                        panel.get_haplotype_function(gene, haplotype_call.haplotype_name),
                        gene_to_drug_info[gene][0],
                        gene_to_drug_info[gene][1],
                        panel_path,
                        version,
                    ]))
            else:
                # TODO: make these strings into constants or enums?
                lines.append(cls.TSV_SEPARATOR.join([
                    gene,
                    "Unresolved Haplotype",
                    "Unknown Function",
                    gene_to_drug_info[gene][0],
                    gene_to_drug_info[gene][1],
                    panel_path,
                    version,
                ]))
        text = "\n".join(lines) + "\n"
        return text

    @classmethod
    def __get_panel_calls_df(cls, pgx_analysis: PgxAnalysis) -> pd.DataFrame:
        return cls.__get_calls_data_frame_from_full_calls(pgx_analysis.get_all_full_calls())

    @classmethod
    def __get_calls_data_frame_from_full_calls(cls, full_calls: FrozenSet[FullCall]) -> pd.DataFrame:
        # TODO: add ref alleles to data frame?
        data_frame = pd.DataFrame(columns=cls.CALLS_TSV_COLUMNS)
        for full_call in full_calls:
            annotated_alleles = full_call.get_annotated_alleles()

            grch37_ref_alleles = [annotated.allele for annotated in annotated_alleles if not annotated.is_variant_vs_grch37]
            grch37_variant_alleles = [annotated.allele for annotated in annotated_alleles if annotated.is_variant_vs_grch37]
            grch37_alleles = grch37_ref_alleles + grch37_variant_alleles

            grch38_ref_alleles = [
                annotated.allele for annotated in annotated_alleles
                if annotated.is_annotated_vs_grch38() and not annotated.is_variant_vs_grch38
            ]
            grch38_variant_alleles = [
                annotated.allele for annotated in annotated_alleles
                if not annotated.is_annotated_vs_grch38() or annotated.is_variant_vs_grch38
            ]
            grch38_alleles = (grch38_ref_alleles + grch38_variant_alleles)

            position_grch38 = (
                cls.__get_coordinate_string(full_call.start_coordinate_grch38)
                if full_call.start_coordinate_grch38 is not None
                else "UNKNOWN"
            )

            new_id = {
                "position_GRCh37": cls.__get_coordinate_string(full_call.start_coordinate_grch37),
                "rsid": cls.RS_ID_SEPARATOR.join(list(full_call.rs_ids)),
                "ref_GRCh37": grch37_alleles[0],
                "alt_GRCh37": grch37_alleles[1],
                "variant_annotation": full_call.variant_annotation,
                "filter": full_call.filter,
                "gene": full_call.gene,
                "position_GRCh38": position_grch38,
                "ref_GRCh38": grch38_alleles[0],
                "alt_GRCh38": grch38_alleles[1]
            }
            data_frame = data_frame.append(new_id, ignore_index=True)

        if pd.isna(data_frame).any(axis=None):
            raise ValueError(f"This should never happen: Unhandled NaN values:\n{data_frame}")

        data_frame = cls.__sort_call_data_frame(data_frame)

        return data_frame

    @classmethod
    def __sort_call_data_frame(cls, data_frame: pd.DataFrame) -> pd.DataFrame:
        data_frame["chromosome"] = (data_frame["position_GRCh37"]).str.split(pat=":", n=1, expand=True).iloc[:, 0]
        data_frame = data_frame.sort_values(by=["chromosome", "position_GRCh37"]).reset_index(drop=True)
        data_frame = data_frame.loc[:, cls.CALLS_TSV_COLUMNS]
        return data_frame

    @classmethod
    def __get_haplotype_call_string(cls, haplotype_call: HaplotypeCall) -> str:
        return str(haplotype_call)

    @classmethod
    def __get_coordinate_string(cls, coordinate: GeneCoordinate) -> str:
        return str(coordinate)
