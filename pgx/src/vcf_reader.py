from typing import Dict, Any, Tuple

import allel

from base.gene_coordinate import GeneCoordinate
from call_data import Grch37CallData, Grch37Call
from config.panel import Panel


class VcfReader(object):
    @classmethod
    def get_grch37_call_data(cls, filtered_vcf: str, panel: Panel) -> Grch37CallData:
        variants = cls.__get_variants_from_filtered_vcf(filtered_vcf)
        return cls.__get_call_data_from_variants(variants, panel)

    @classmethod
    def __get_variants_from_filtered_vcf(cls, filtered_vcf: str) -> Dict[str, Any]:
        try:
            field_names = ['samples', 'calldata/GT', 'variants/ALT', 'variants/CHROM', 'variants/FILTER', 'variants/ID',
                           'variants/POS', 'variants/QUAL', 'variants/REF', 'variants/ANN']
            variants = allel.read_vcf(filtered_vcf, fields=field_names, transformers=allel.ANNTransformer())
        except IOError:
            raise FileNotFoundError("File " + filtered_vcf + " not found or cannot be opened.")
        return variants

    @classmethod
    def __get_call_data_from_variants(cls, variants: Dict[str, Any], panel: Panel) -> Grch37CallData:
        match_on_rsid = 0
        match_on_location = 0
        filtered_calls = set()
        for i, rs_ids_string in enumerate(variants['variants/ID']):
            chromosome = str(variants['variants/CHROM'][i])
            position = int(variants['variants/POS'][i])
            reference_allele = str(variants['variants/REF'][i])

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
                if variants['variants/FILTER_PASS'][i]:
                    filter_type = "PASS"
                else:
                    filter_type = "FILTERED"
                alts = [str(allele) for allele in variants['variants/ALT'][i]]
                variant_annotation = str(variants['variants/ANN_HGVS_c'][i])
                gene = str(variants['variants/ANN_Gene_Name'][i])
                genotype = variants['calldata/GT'][i][0].tolist()
                if genotype == [0, 1]:
                    alleles = (reference_allele, alts[0])
                elif genotype == [1, 1]:
                    alleles = (alts[0], alts[0])
                elif genotype == [1, 2]:
                    alleles = (alts[0], alts[1])
                elif genotype == [0, 0]:
                    alleles = (reference_allele, reference_allele)
                    variant_annotation = "REF_CALL"
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
        if ";" in rs_ids_string:
            return tuple(str(rs) for rs in rs_ids_string.split(";") if rs.startswith("rs"))
        else:
            return (str(rs_ids_string),)

    @classmethod
    def __get_relevant_coordinates(cls, chromosome: str, position: int, ref_allele: str) -> Tuple[GeneCoordinate, ...]:
        return tuple(GeneCoordinate(chromosome, position + i) for i in range(len(ref_allele)))
