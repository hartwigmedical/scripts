from typing import FrozenSet, Dict, Union

import pandas as pd

from base.constants import UNKNOWN_FUNCTION_STRING
from base.util import get_chromosome_name_to_index
from call_data import FullCall, HaplotypeCall
from config.panel import Panel
from pgx_analysis import PgxAnalysis


class GenotypeReporter(object):
    GENE_COLUMN_NAME = "gene"
    CHROMOSOME_COLUMN_NAME = "chromosome"
    POSITION_V37_COLUMN_NAME = "position_v37"
    POSITION_V38_COLUMN_NAME = "position_v38"
    REF_ALLELE_V37_COLUMN_NAME = "ref_v37"
    REF_ALLELE_V38_COLUMN_NAME = "ref_v38"
    FIRST_ALLELE_COLUMN_NAME = "allele1"
    SECOND_ALLELE_COLUMN_NAME = "allele2"
    RS_IDS_COLUMN_NAME = "rsid"
    ANNOTATION_V37_COLUMN_NAME = "variant_annotation_v37"
    ANNOTATION_V38_COLUMN_NAME = "variant_annotation_v38"
    FILTER_V37_COLUMN_NAME = "filter_v37"
    FILTER_V38_COLUMN_NAME = "filter_v38"
    PANEL_VERSION_COLUMN_NAME = "panel_version"
    TOOL_VERSION_COLUMN_NAME = "repo_version"
    NEW_CALLS_TSV_COLUMNS = (
        GENE_COLUMN_NAME, CHROMOSOME_COLUMN_NAME,
        POSITION_V37_COLUMN_NAME, POSITION_V38_COLUMN_NAME,
        REF_ALLELE_V37_COLUMN_NAME, REF_ALLELE_V38_COLUMN_NAME,
        FIRST_ALLELE_COLUMN_NAME, SECOND_ALLELE_COLUMN_NAME, RS_IDS_COLUMN_NAME,
        ANNOTATION_V37_COLUMN_NAME, FILTER_V37_COLUMN_NAME, ANNOTATION_V38_COLUMN_NAME, FILTER_V38_COLUMN_NAME,
        PANEL_VERSION_COLUMN_NAME, TOOL_VERSION_COLUMN_NAME,
    )
    CHROMOSOME_INDEX_NAME = "chromosome_index"

    UNKNOWN_POSITION_STRING = "UNKNOWN"
    UNKNOWN_REF_ALLELE_STRING = "UNKNOWN"

    TSV_SEPARATOR = "\t"
    RS_ID_SEPARATOR = ";"

    @classmethod
    def get_calls_tsv_text(cls, pgx_analysis: PgxAnalysis, panel_id: str, version: str) -> str:
        text = str(
            cls.__get_panel_calls_df(pgx_analysis, panel_id, version).to_csv(sep=cls.TSV_SEPARATOR, index=False)
        )
        return text

    @classmethod
    def __get_panel_calls_df(cls, pgx_analysis: PgxAnalysis, panel_id: str, version: str) -> pd.DataFrame:
        return cls.__get_calls_data_frame_from_full_calls(pgx_analysis.get_all_full_calls(), panel_id, version)

    @classmethod
    def __get_calls_data_frame_from_full_calls(
            cls, full_calls: FrozenSet[FullCall], panel_id: str, version: str) -> pd.DataFrame:
        data_frame = pd.DataFrame(columns=cls.NEW_CALLS_TSV_COLUMNS)
        for full_call in full_calls:
            assert (
                    full_call.start_coordinate_v38 is None
                    or full_call.start_coordinate_v37.chromosome == full_call.start_coordinate_v38.chromosome), \
                f"Chromosomes of two coordinates are different: call={full_call}"

            sorted_alleles = sorted(full_call.alleles)
            position_v38 = (
                str(full_call.start_coordinate_v38.position)
                if full_call.start_coordinate_v38 is not None
                else cls.UNKNOWN_POSITION_STRING
            )
            reference_allele_v38 = (
                full_call.reference_allele_v38
                if full_call.reference_allele_v38 is not None
                else cls.UNKNOWN_REF_ALLELE_STRING
            )

            new_id: Dict[str, Union[str, int]] = {
                cls.GENE_COLUMN_NAME: full_call.gene,
                cls.CHROMOSOME_COLUMN_NAME: full_call.start_coordinate_v37.chromosome,
                cls.POSITION_V37_COLUMN_NAME: full_call.start_coordinate_v37.position,
                cls.POSITION_V38_COLUMN_NAME: position_v38,
                cls.REF_ALLELE_V37_COLUMN_NAME: full_call.reference_allele_v37,
                cls.REF_ALLELE_V38_COLUMN_NAME: reference_allele_v38,
                cls.FIRST_ALLELE_COLUMN_NAME: sorted_alleles[0],
                cls.SECOND_ALLELE_COLUMN_NAME: sorted_alleles[1],
                cls.RS_IDS_COLUMN_NAME: cls.RS_ID_SEPARATOR.join(full_call.rs_ids),
                cls.ANNOTATION_V37_COLUMN_NAME: full_call.variant_annotation_v37,
                cls.FILTER_V37_COLUMN_NAME: full_call.filter_v37.name,
                cls.ANNOTATION_V38_COLUMN_NAME: full_call.variant_annotation_v38,
                cls.FILTER_V38_COLUMN_NAME: full_call.filter_v38.name,
                cls.PANEL_VERSION_COLUMN_NAME: panel_id,
                cls.TOOL_VERSION_COLUMN_NAME: version,
            }
            data_frame = data_frame.append(new_id, ignore_index=True)

        if pd.isna(data_frame).any(axis=None):
            raise ValueError(f"This should never happen: Unhandled NaN values:\n{data_frame}")

        data_frame = cls.__sort_call_data_frame(data_frame)

        return data_frame

    @classmethod
    def __sort_call_data_frame(cls, data_frame: pd.DataFrame) -> pd.DataFrame:
        data_frame.loc[:, cls.CHROMOSOME_INDEX_NAME] = (
            data_frame.loc[:, cls.CHROMOSOME_COLUMN_NAME].map(get_chromosome_name_to_index())
        )
        data_frame = data_frame.sort_values(
            by=[cls.CHROMOSOME_INDEX_NAME, cls.POSITION_V37_COLUMN_NAME, cls.REF_ALLELE_V37_COLUMN_NAME]
        ).reset_index(drop=True)
        data_frame = data_frame.loc[:, cls.NEW_CALLS_TSV_COLUMNS]
        return data_frame


class HaplotypeReporter(object):
    GENOTYPE_TSV_COLUMNS = (
        "gene", "haplotype", "function", "linked_drugs", "url_prescription_info", "panel_version", "repo_version")

    HAPLOTYPE_HOMOZYGOUS_SUFFIX = "_HOM"
    HAPLOTYPE_HETEROZYGOUS_SUFFIX = "_HET"

    UNRESOLVED_HAPLOTYPE_STRING = "Unresolved Haplotype"

    TSV_SEPARATOR = "\t"
    DRUG_SEPARATOR = ";"

    @classmethod
    def get_genotype_tsv_text(cls, pgx_analysis: PgxAnalysis, panel: Panel, panel_id: str, version: str) -> str:
        gene_to_haplotype_calls = pgx_analysis.get_gene_to_haplotype_calls()

        genes_in_analysis = set(gene_to_haplotype_calls.keys())
        assert genes_in_analysis == panel.get_genes(), (
            f"Gene lists inconsistent.\n"
            f"From analysis={sorted(list(genes_in_analysis))}\n"
            f"From panel={sorted(list(panel.get_genes()))}"
        )

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

        header = cls.TSV_SEPARATOR.join(cls.GENOTYPE_TSV_COLUMNS)
        lines = [header]
        for gene in sorted(gene_to_haplotype_calls.keys()):
            if gene_to_haplotype_calls[gene]:
                for haplotype_call in sorted(gene_to_haplotype_calls[gene], key=lambda call: call.haplotype_name):
                    lines.append(cls.TSV_SEPARATOR.join([
                        gene,
                        cls.__get_haplotype_call_string(haplotype_call),
                        panel.get_haplotype_function(gene, haplotype_call.haplotype_name),
                        gene_to_drug_info[gene][0],
                        gene_to_drug_info[gene][1],
                        panel_id,
                        version,
                    ]))
            else:
                lines.append(cls.TSV_SEPARATOR.join([
                    gene,
                    cls.UNRESOLVED_HAPLOTYPE_STRING,
                    UNKNOWN_FUNCTION_STRING,
                    gene_to_drug_info[gene][0],
                    gene_to_drug_info[gene][1],
                    panel_id,
                    version,
                ]))
        text = "\n".join(lines) + "\n"
        return text

    @classmethod
    def __get_haplotype_call_string(cls, haplotype_call: HaplotypeCall) -> str:
        if haplotype_call.count == 2:
            return haplotype_call.haplotype_name + cls.HAPLOTYPE_HOMOZYGOUS_SUFFIX
        elif haplotype_call.count == 1:
            return haplotype_call.haplotype_name + cls.HAPLOTYPE_HETEROZYGOUS_SUFFIX
        else:
            error_msg = f"Invalid haplotype count: haplotype call={haplotype_call}"
            raise ValueError(error_msg)
