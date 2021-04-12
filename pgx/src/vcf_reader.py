from typing import Dict, Any, Tuple, Optional

import allel

from base.constants import REF_CALL_ANNOTATION_STRING
from base.filter import Filter
from base.gene_coordinate import GeneCoordinate
from call_data import Grch37CallData, Grch37Call
from config.panel import Panel


class VcfReader(object):
    ALT_ALLELE_FIELD_NAME = "variants/ALT"
    ANNOTATION_FIELD_NAME = "variants/ANN"
    CHROMOSOME_FIELD_NAME = "variants/CHROM"
    FILTER_FIELD_NAME = "variants/FILTER"
    GENOTYPE_FIELD_NAME = "calldata/GT"
    POSITION_FIELD_NAME = "variants/POS"
    QUALITY_FIELD_NAME = "variants/QUAL"
    REF_ALLELE_FIELD_NAME = "variants/REF"
    RS_IDS_FIELD_NAME = "variants/ID"
    SAMPLE_FIELD_NAME = "samples"
    FIELD_NAMES = [
        SAMPLE_FIELD_NAME, GENOTYPE_FIELD_NAME, ALT_ALLELE_FIELD_NAME, CHROMOSOME_FIELD_NAME, FILTER_FIELD_NAME,
        RS_IDS_FIELD_NAME, POSITION_FIELD_NAME, QUALITY_FIELD_NAME, REF_ALLELE_FIELD_NAME, ANNOTATION_FIELD_NAME,
    ]

    RS_ID_SEPARATOR = ";"

    @classmethod
    def get_grch37_call_data(cls, filtered_vcf: str, panel: Panel) -> Grch37CallData:
        variants = cls.__get_variants_from_filtered_vcf(filtered_vcf)
        return cls.__get_call_data_from_variants(variants, panel)

    @classmethod
    def __get_variants_from_filtered_vcf(cls, filtered_vcf: str) -> Optional[Dict[str, Any]]:
        # variants is None precisely when filtered vcf file has no variants
        try:
            variants = allel.read_vcf(filtered_vcf, fields=cls.FIELD_NAMES, transformers=allel.ANNTransformer())
        except IOError:
            raise FileNotFoundError("File " + filtered_vcf + " not found or cannot be opened.")
        return variants

    @classmethod
    def __get_call_data_from_variants(cls, variants: Optional[Dict[str, Any]], panel: Panel) -> Grch37CallData:
        match_on_rsid = 0
        match_on_location = 0
        filtered_calls = set()
        if variants is not None:
            for i, rs_ids_string in enumerate(variants[cls.RS_IDS_FIELD_NAME]):
                chromosome = str(variants[cls.CHROMOSOME_FIELD_NAME][i])
                position = int(variants[cls.POSITION_FIELD_NAME][i])
                reference_allele = str(variants[cls.REF_ALLELE_FIELD_NAME][i])

                rs_ids = cls.__get_rs_ids_from_string(rs_ids_string)
                relevant_coordinates = cls.__get_relevant_coordinates(chromosome, position, reference_allele)

                rs_id_match_to_panel_exists = any(panel.contains_rs_id(rs_id) for rs_id in rs_ids)
                coordinate_match_to_panel_exists = any(
                    panel.contains_rs_id_with_grch37_coordinate(coord) for coord in relevant_coordinates
                )
                if rs_id_match_to_panel_exists or coordinate_match_to_panel_exists:
                    if rs_id_match_to_panel_exists:
                        match_on_rsid += 1
                    if coordinate_match_to_panel_exists:
                        match_on_location += 1
                    if variants[f"{cls.FILTER_FIELD_NAME}_PASS"][i]:
                        filter_type = Filter.PASS
                    else:
                        # Ignore all calls with filter != PASS
                        continue
                    alts = [str(allele) for allele in variants[cls.ALT_ALLELE_FIELD_NAME][i]]
                    variant_annotation = str(variants[f"{cls.ANNOTATION_FIELD_NAME}_HGVS_c"][i])
                    gene = str(variants[f"{cls.ANNOTATION_FIELD_NAME}_Gene_Name"][i])
                    genotype = variants[cls.GENOTYPE_FIELD_NAME][i][0].tolist()
                    if genotype == [0, 1]:
                        alleles = (reference_allele, alts[0])
                    elif genotype == [1, 1]:
                        alleles = (alts[0], alts[0])
                    elif genotype == [1, 2]:
                        alleles = (alts[0], alts[1])
                    elif genotype == [0, 0]:
                        alleles = (reference_allele, reference_allele)
                        variant_annotation = REF_CALL_ANNOTATION_STRING
                    else:
                        error_msg = f"Genotype not found: {genotype}"
                        raise ValueError(error_msg)

                    call = Grch37Call(
                        GeneCoordinate(chromosome, position),
                        reference_allele,
                        alleles,
                        gene,
                        rs_ids,
                        variant_annotation,
                        filter_type,
                    )
                    filtered_calls.add(call)

        print("[INFO] Matches on RS id: " + str(match_on_rsid))
        print("[INFO] Matches on location: " + str(match_on_location))

        return Grch37CallData(frozenset(filtered_calls))

    @classmethod
    def __get_rs_ids_from_string(cls, rs_ids_string: str) -> Tuple[str, ...]:
        if cls.RS_ID_SEPARATOR in rs_ids_string:
            return tuple(str(rs) for rs in rs_ids_string.split(cls.RS_ID_SEPARATOR) if rs.startswith("rs"))
        else:
            return (str(rs_ids_string),)

    @classmethod
    def __get_relevant_coordinates(cls, chromosome: str, position: int, ref_allele: str) -> Tuple[GeneCoordinate, ...]:
        return tuple(GeneCoordinate(chromosome, position + i) for i in range(len(ref_allele)))
