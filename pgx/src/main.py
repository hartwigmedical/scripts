from typing import List, DefaultDict, Dict, Any

import allel
import argparse
import collections
import itertools
import json
import os
import pandas as pd
import subprocess
import sys
from shutil import copyfile

from panel import Panel


def main(vcf: str, sample_t_id: str, sample_r_id: str, version: str, panel_path: str, outputdir: str,
         recreate_bed: bool, vcftools: str, sourcedir: str):
    """ Run pharmacogenomics analysis on sample """
    print("\n[INFO] ## START PHARMACOGENOMICS ANALYSIS")

    # Check if output dir exists, create if it does not
    if not os.path.exists(outputdir):
        try:
            os.makedirs(outputdir)
        except FileExistsError:
            # Directory already exists
            pass

    panel = load_panel_from_json(panel_path)

    if panel.is_empty():
        sys.exit("[ERROR] No panel variants are given, no analysis is performed.")

    bed_file = get_bed_file(panel_path, recreate_bed, panel, sourcedir)
    filtered_vcf = get_filtered_vcf(vcf, bed_file, sample_r_id, sample_t_id, outputdir, vcftools)

    variants = get_variants_from_filtered_vcf(filtered_vcf)
    ids_found_in_patient = get_ids_found_in_patient_from_variants(variants, panel)

    ids_found_in_patient, results, severity, all_ids_in_panel, drug_info = \
        convert_results_into_haplotypes(ids_found_in_patient, panel, panel_path)

    out = outputdir + "/" + sample_t_id
    print_calls_to_file(out + "_calls.txt", all_ids_in_panel)
    print_haplotypes_to_file(out + "_genotype.txt", drug_info, panel_path, results, severity, version)
    # Also copy the bed-filtered VCF file for research purposes
    copyfile(filtered_vcf, out + "_PGx.vcf")

    # Clean up filtered_vcf
    if os.path.exists(filtered_vcf):
        if os.path.exists(filtered_vcf):
            os.remove(filtered_vcf)
            print("[INFO] " + filtered_vcf + " removed.")
        if os.path.exists(filtered_vcf.replace(".recode.vcf", ".log")):
            os.remove(filtered_vcf.replace(".recode.vcf", ".log"))
            print("[INFO] " + filtered_vcf.replace(".recode.vcf", ".log") + " removed.")

    # TODO: add genes CYP2D6, CYP3A4, CYP3A5

    print("[INFO] ## PHARMACOGENOMICS ANALYSIS FINISHED\n")


def convert_results_into_haplotypes(ids_found_in_patient: pd.DataFrame, panel: Panel, panel_path: str):
    rsid_to_gene_to_haplotype_variant: DefaultDict[str, Dict[str, Any]] = collections.defaultdict(dict)
    for gene_info in panel.gene_to_gene_info.values():
        for variant in gene_info.variants:
            rsid_to_gene_to_haplotype_variant[variant['rsid']][gene_info.gene] = variant

    results = {}
    severity = {}
    drug_info = {}

    # Process the exceptions
    ids_found_in_patient = process_exceptions(ids_found_in_patient, panel_path)
    # Fill in the GRCh38 gaps, they should be the same as the GRCh37 equivalent
    if not ids_found_in_patient.empty:
        for index, row in ids_found_in_patient.iterrows():
            if 'ref_GRCh38' in row:
                if pd.isna(row['ref_GRCh38']):
                    ids_found_in_patient.at[index, 'ref_GRCh38'] = row['ref_GRCh37']
                    ids_found_in_patient.at[index, 'alt_GRCh38'] = row['alt_GRCh37']
            else:
                ids_found_in_patient.at[index, 'ref_GRCh38'] = row['ref_GRCh37']
                ids_found_in_patient.at[index, 'alt_GRCh38'] = row['alt_GRCh37']

        for index, row in ids_found_in_patient.iterrows():
            # Fill in the rsid gaps, if annotation is missing on the location.
            if 'rsid' in row.index and row['rsid'] == '.':
                position_string = row['position_GRCh37']
                if panel.contains_snp_with_position(position_string):
                    rs_id = panel.get_snp_with_position(position_string).rs_id
                else:
                    rs_id = '.'
                ids_found_in_patient.at[index, 'rsid'] = rs_id

        for index, row in ids_found_in_patient.iterrows():
            matching_variants = rsid_to_gene_to_haplotype_variant[ids_found_in_patient.at[index, 'rsid']].values()
            grch38_locations = {
                str(variant['chromosome']) + ":" + str(variant['positionGRCh38']) for variant in matching_variants
            }
            if len(grch38_locations) == 1:
                grch38_location = grch38_locations.pop()
                if 'position_GRCh38' not in row.index or pd.isna(row['position_GRCh38']):
                    ids_found_in_patient.at[index, 'position_GRCh38'] = grch38_location
                elif row['position_GRCh38'] != grch38_location:
                    raise ValueError(
                        "[ERROR] Inconsistent GRCh38 locations for variants:\n"
                        "location from exceptions: " + row['position_GRCh38'] + "\n"
                        "location from variants: " + grch38_location + "\n"
                    )
            elif len(grch38_locations) > 1:
                matching_variants_string = ",".join([str(variant) for variant in matching_variants])
                raise ValueError("[ERROR] Inconsistent GRCh38 locations for variants:\n"
                                 "matching variants: " + matching_variants_string + "\n"
                                 "GRCh38 locations: " + ", ".join(grch38_locations) + "\n")

    # Generate a list of ids not found in patient
    ids_not_found_in_patient = pd.DataFrame(columns=['position_GRCh37', 'ref_GRCh37', 'alt_GRCh37', 'position_GRCh38',
                                                     'ref_GRCh38', 'alt_GRCh38', 'rsid', 'variant_annotation', 'gene',
                                                     'filter'])

    rs_ids_found_in_patient = set(ids_found_in_patient.rsid.tolist())
    positions_found_in_patient = set(ids_found_in_patient.position_GRCh37.tolist())

    for snp in panel.snps:
        if snp.rs_id not in rs_ids_found_in_patient and snp.get_position_string() not in positions_found_in_patient:
            new_id = {}
            for gene, variant in rsid_to_gene_to_haplotype_variant[snp.rs_id].items():
                new_id['position_GRCh37'] = snp.get_position_string()
                new_id['rsid'] = snp.rs_id
                new_id['ref_GRCh37'] = variant['referenceAlleleGRCh38']
                new_id['alt_GRCh37'] = variant['referenceAlleleGRCh38']  # Assuming REF/REF relative to GRCh38
                new_id['variant_annotation'] = "REF_CALL"
                new_id['filter'] = "NO_CALL"
                new_id['gene'] = gene
                new_id['ref_GRCh38'] = variant['referenceAlleleGRCh38']  # Again assuming REF/REF relative to GRCh38
                new_id['alt_GRCh38'] = variant['referenceAlleleGRCh38']
                new_id['position_GRCh38'] = str(variant['chromosome']) + ":" + str(variant['positionGRCh38'])
                ids_not_found_in_patient = ids_not_found_in_patient.append(new_id, ignore_index=True)

    # Now we want to process all the variants in terms of the alleles
    all_ids_in_panel = pd.concat([ids_found_in_patient, ids_not_found_in_patient], sort=True)
    all_ids_in_panel = all_ids_in_panel.sort_values(by='position_GRCh37').reset_index(drop=True)
    all_ids_in_panel = all_ids_in_panel[['gene', 'position_GRCh37', 'ref_GRCh37', 'alt_GRCh37', 'position_GRCh38',
                                         'ref_GRCh38', 'alt_GRCh38', 'rsid', 'variant_annotation', 'filter']]

    for gene_info in panel.gene_to_gene_info.values():
        print("[INFO] PROCESSING GENE " + gene_info.gene)
        ids_found_in_gene = all_ids_in_panel[all_ids_in_panel['gene'].str.contains(gene_info.gene)]
        perfect_match = False
        severity[gene_info.reference_allele] = "Normal Function"
        severity['Unresolved'] = "Unknown Function"
        drug_info[gene_info.gene] = [
            ";".join([x['name'] for x in gene_info.drugs]),
            ";".join(x['url_prescription_info'] for x in gene_info.drugs)
        ]

        # If all variants are assumed_ref, return reference allele
        if len(ids_found_in_gene.loc[ids_found_in_gene['variant_annotation'] == "REF_CALL"]) == len(ids_found_in_gene):
            print("[INFO] Found reference allele")
            results[gene_info.gene] = [gene_info.reference_allele + "_HOM"]
        else:
            results[gene_info.gene] = []
            haplotypes_matching = []
            for allele in gene_info.allelles:
                severity[allele['alleleName']] = allele['function']
                vars_found_in_gene = ids_found_in_gene.loc[ids_found_in_gene['variant_annotation'] != "REF_CALL"]
                variants_sample = list(zip(vars_found_in_gene.rsid.tolist(), vars_found_in_gene.alt_GRCh38.tolist()))
                if set(variants_sample) == set([(x['rsid'], x['altAlleleGRCh38']) for x in allele['alleleVariants']]):
                    perfect_match = True
                    print("[INFO] Found 1:1 match with allele " + allele['alleleName'])
                    # Now we want to see if we have hetrozygous or homozygous calls
                    allele_statuses = []
                    for index, row in vars_found_in_gene.iterrows():
                        if row['ref_GRCh38'] == row['alt_GRCh38']:
                            allele_statuses.append("HOM")
                        else:
                            allele_statuses.append("HET")
                    if all(x == allele_statuses[0] for x in allele_statuses):
                        allele_status = allele_statuses[0]
                    else:
                        allele_status = "HOMHET"
                    # Add to results
                    results[gene_info.gene].append(allele['alleleName'] + "_" + str(allele_status))
                    if allele_status == "HET":
                        # Assume if perfect match with HET, we are also looking at reference allele
                        results[gene_info.gene].append(gene_info.reference_allele + "_HET")
                    break
                else:
                    # print("Processing " + str(allele['alleleName']))
                    # print(set([(x['rsid'], x['altAlleleGRCh38']) for x in allele['alleleVariants']]))
                    # print(set(variants_sample))
                    if set([(x['rsid'], x['altAlleleGRCh38']) for x in allele['alleleVariants']]) <= \
                            set(variants_sample):
                        print("[INFO] A subset of rsids matches " + str(allele['alleleName']) + " in part")
                        haplotypes_matching.append(allele)

            if not perfect_match:
                if not haplotypes_matching:
                    print(
                        f"[WARN] No haplotype match found for {gene_info.gene}. "
                        f"Probable cause is that the variant is not in line with previously "
                        f"determined within defined haplotype."
                    )
                    results[gene_info.gene].append("Unresolved_Haplotype")
                else:
                    print("[INFO] Test all possible combinations of haplotypes to see if a perfect match can be found")
                    optimal_set: List[List[Any]] = []
                    for k in range(len(haplotypes_matching) + 1, 0, -1):  # TODO: shouldn't this order be reversed?
                        for subset in itertools.combinations(haplotypes_matching, k):
                            if perfect_match:
                                continue
                            # See if this combination results in a perfect match, otherwise store the score
                            allele_variants = [x['alleleVariants'] for x in subset]
                            rs_ids_subset = []
                            for x in allele_variants:
                                for var in x:
                                    rs_ids_subset.append((var['rsid'], var['altAlleleGRCh38']))
                            if compare_collection(variants_sample, rs_ids_subset):
                                print("[INFO] Perfect haplotype combination found!")
                                perfect_match = True
                                for allele in subset:
                                    allele_statuses = []
                                    rs_ids_in_allele = [x['rsid'] for x in allele['alleleVariants']]
                                    found_vars = ids_found_in_gene[ids_found_in_gene['rsid'].isin(rs_ids_in_allele)]
                                    for index, row in found_vars.iterrows():
                                        if row['ref_GRCh38'] == row['alt_GRCh38']:
                                            allele_statuses.append("HOM")
                                        else:
                                            allele_statuses.append("HET")
                                    if all(x == allele_statuses[0] for x in allele_statuses):
                                        allele_status = allele_statuses[0]
                                    else:
                                        allele_status = "HOMHET"
                                    results[gene_info.gene].append(allele['alleleName'] + "_" + str(allele_status))
                            else:
                                rs_matched = list(set(variants_sample) & set(rs_ids_subset))
                                rs_not_in_haplotype = list(set(variants_sample) - set(rs_ids_subset))
                                rs_not_found = list(set(rs_ids_subset) - set(variants_sample))
                                optimal_set.append([len(rs_matched), subset, rs_matched, rs_not_found,
                                                    rs_not_in_haplotype])
                    if not perfect_match:
                        # Get best scoring set and give as result
                        print("[INFO] No perfect match found. Start testing all possible combinations and choose best "
                              "one")
                        optimal_set.sort(key=lambda x: x[0], reverse=True)
                        for options in optimal_set:
                            print("[INFO] Option to be tested:")
                            print("[INFO]\t\t# Matches: " + str(options[0]))
                            print("[INFO]\t\tSubset: " + str(options[1]))
                            print("[INFO]\t\tVariants matched for sample and haplotype: " + str(options[2]))
                            print("[INFO]\t\tVariants in tested haplotype, but not in sample: " + str(options[3]))
                            print("[INFO]\t\tVariants in sample, not in tested haplotype: " + str(options[4]))
                        # Here we just pick the top option.
                        if optimal_set[0][0] >= 1:
                            subset = optimal_set[0][1]
                            # If we have a variant that is in sample or in haplotype that is not matched > undetermined
                            if len(optimal_set[0][3]) > 0 or len(optimal_set[0][4]) > 0:
                                results[gene_info.gene].append("Unresolved_Haplotype")
                            else:
                                for allele in subset:
                                    allele_statuses = []
                                    rs_ids_in_allele = [x['rsid'] for x in allele['alleleVariants']]
                                    found_vars = ids_found_in_gene[ids_found_in_gene['rsid'].isin(rs_ids_in_allele)]
                                    for index, row in found_vars.iterrows():
                                        if row['ref_GRCh38'] == row['alt_GRCh38']:
                                            allele_statuses.append("HOM")
                                        else:
                                            allele_statuses.append("HET")
                                    if all(x == allele_statuses[0] for x in allele_statuses):
                                        allele_status = allele_statuses[0]
                                    else:
                                        allele_status = "HOMHET"
                                    results[gene_info.gene].append(allele['alleleName'] + "_" + str(allele_status))
                        else:
                            sys.exit("[ERROR] No haplotype match was found. Exiting.")
        # If we only find one allele and it is HET, assume we're also dealing with reference allele
        if len(results[gene_info.gene]) == 1:
            if results[gene_info.gene][0].split("_")[-1] == "HET":
                results[gene_info.gene].append(gene_info.reference_allele + "_HET")

    return ids_found_in_patient, results, severity, all_ids_in_panel, drug_info


def process_exceptions(ids_found: pd.DataFrame, panel_path) -> pd.DataFrame:
    try:
        panel_exceptions_file = panel_path.replace(".json", "") + "_exceptions.json"
        with open(panel_exceptions_file, 'r+', encoding='utf-8') as json_file:
            exceptions = json.load(json_file)
            exceptions = exceptions['exceptions']
            for variant in exceptions:
                variant_location = str(variant['chromosome']) + ":" + str(variant['position'])
                variant_location_found = variant_location in ids_found['position_GRCh37'].tolist()
                rs_id_found = variant['rsid'] in ids_found['rsid'].tolist()
                if rs_id_found or variant_location_found:
                    if rs_id_found:
                        # get line and index from ids_found
                        found_var = ids_found[ids_found['rsid'].str.contains(variant['rsid'])]
                    else:
                        found_var = ids_found[ids_found['position_GRCh37'].str.contains(variant_location)]
                    # Delete found variant from results if the ref base and alt base are the correctedRefBase
                    if found_var['ref_GRCh37'].values == variant['referenceAlleleGRCh38'] and \
                            found_var['alt_GRCh37'].values == variant['referenceAlleleGRCh38']:
                        ids_found = ids_found.drop(found_var.index[0])
                    elif found_var['ref_GRCh37'].values == variant['referenceAlleleGRCh38'] and \
                            found_var['alt_GRCh37'].values == variant['altAlleleGRCh38']:
                        print("[WARN] What should we do? ref = corRef, alt = corAlt")
                    elif found_var['ref_GRCh37'].values == variant['altAlleleGRCh38'] and \
                            found_var['alt_GRCh37'].values == variant['referenceAlleleGRCh38']:
                        # Change variant_annotation and ref and alt base
                        ids_found.at[found_var.index[0], 'variant_annotation'] = variant['annotationGRCh38']
                        ids_found.at[found_var.index[0], 'ref_GRCh38'] = variant['referenceAlleleGRCh38']
                        ids_found.at[found_var.index[0], 'alt_GRCh38'] = variant['altAlleleGRCh38']
                        ids_found.at[found_var.index[0], 'position_GRCh38'] = (
                                str(variant['chromosome']) + ":" + str(variant['positionGRCh38'])
                        )
                    elif found_var['ref_GRCh37'].values == variant['altAlleleGRCh38'] and \
                            found_var['alt_GRCh37'].values == variant['altAlleleGRCh38']:
                        # Add variant_annotation and ref and alt base for hg38
                        ids_found.at[found_var.index[0], 'variant_annotation'] = variant['annotationGRCh38']
                        ids_found.at[found_var.index[0], 'ref_GRCh38'] = variant['altAlleleGRCh38']
                        ids_found.at[found_var.index[0], 'alt_GRCh38'] = variant['altAlleleGRCh38']
                        ids_found.at[found_var.index[0], 'position_GRCh38'] = (
                                str(variant['chromosome']) + ":" + str(variant['positionGRCh38'])
                        )
                    else:
                        print("[ERROR] Complete mismatch:")
                        print(found_var)
                        raise ValueError("[ERROR] Exceptions cannot be processed. Please check. Exiting.")
                else:
                    print("[INFO] Exception variant not found in this patient. This means ref/ref call, but should be "
                          "flipped. Add to table.")
                    new_id = {'position_GRCh37': str(variant['chromosome']) + ":" + str(variant['position']),
                              'rsid': variant['rsid'],
                              'ref_GRCh37': variant['altAlleleGRCh38'],
                              'alt_GRCh37': variant['altAlleleGRCh38'],
                              'variant_annotation': variant['annotationGRCh38'],
                              'filter': "INFERRED_REF_CALL",
                              'gene': variant['gene'],
                              'position_GRCh38': str(variant['chromosome']) + ":" + str(variant['positionGRCh38']),
                              'ref_GRCh38': variant['altAlleleGRCh38'],
                              'alt_GRCh38': variant['altAlleleGRCh38']}
                    ids_found = ids_found.append(new_id, ignore_index=True)
    except IOError:
        sys.exit("[ERROR] Exceptions file not found or cannot be opened.")

    return ids_found


def get_ids_found_in_patient_from_variants(variants, panel: Panel) -> pd.DataFrame:
    match_on_rsid = 0
    match_on_location = 0
    ids_found_in_patient = pd.DataFrame(columns=['position_GRCh37', 'ref_GRCh37', 'alt_GRCh37', 'rsid',
                                                 'variant_annotation', 'gene', 'filter'])
    for i, rs_number in enumerate(variants['variants/ID']):
        chr = variants['variants/CHROM'][i]
        pos = variants['variants/POS'][i]

        if ";" in rs_number:
            rs_id_filter = []
            cur_rs = rs_number.split(";")
            for rs in cur_rs:
                if rs.startswith("rs"):
                    rs_id_filter.append(rs)
        else:
            rs_id_filter = [rs_number]

        rs_id_filter_to_panel_match_exists = any(panel.contains_rs_id(rs_id) for rs_id in rs_id_filter)
        if rs_id_filter_to_panel_match_exists or panel.contains_snp_with_position(f"{chr}:{pos}"):
            if rs_id_filter_to_panel_match_exists:
                match_on_rsid += 1
            else:
                match_on_location += 1
            new_id = {}
            if variants['variants/FILTER_PASS'][i]:
                filter = "PASS"
            else:
                filter = "FILTERED"
            alt = variants['variants/ALT'][i]
            ref = variants['variants/REF'][i]

            # print(variants['calldata/GT'][i][0])
            new_id['variant_annotation'] = variants['variants/ANN_HGVS_c'][i]
            if variants['calldata/GT'][i][0].tolist() == [0, 1]:
                new_id['ref_GRCh37'] = ref
                new_id['alt_GRCh37'] = alt[0]
            elif variants['calldata/GT'][i][0].tolist() == [1, 1]:
                new_id['ref_GRCh37'] = alt[0]
                new_id['alt_GRCh37'] = alt[0]
            elif variants['calldata/GT'][i][0].tolist() == [1, 2]:
                new_id['ref_GRCh37'] = alt[0]
                new_id['alt_GRCh37'] = alt[1]
            elif variants['calldata/GT'][i][0].tolist() == [0, 0]:
                new_id['ref_GRCh37'] = ref
                new_id['alt_GRCh37'] = ref
                new_id['variant_annotation'] = "REF_CALL"
            else:
                print("[ERROR] Genotype not found: " + str(variants['calldata/GT'][i][0].tolist()))

            new_id['position_GRCh37'] = str(chr) + ":" + str(pos)
            new_id['rsid'] = ";".join(rs_id_filter)
            new_id['filter'] = filter
            new_id['gene'] = variants['variants/ANN_Gene_Name'][i]
            ids_found_in_patient = ids_found_in_patient.append(new_id, ignore_index=True)
    print("[INFO] Matches on RS id: " + str(match_on_rsid))
    print("[INFO] Matches on location: " + str(match_on_location))
    return ids_found_in_patient


def get_variants_from_filtered_vcf(filtered_vcf: str):
    try:
        variants = allel.read_vcf(filtered_vcf, fields=['samples', 'calldata/GT', 'variants/ALT', 'variants/CHROM',
                                                        'variants/FILTER', 'variants/ID', 'variants/POS',
                                                        'variants/QUAL', 'variants/REF', 'variants/ANN'],
                                  transformers=allel.ANNTransformer())
    except IOError:
        sys.exit("[ERROR] File " + filtered_vcf + " not found or cannot be opened.")
    return variants


def get_filtered_vcf(vcf: str, bed_file: str, sample_r_id: str, sample_t_id: str, outputdir: str, vcftools: str):
    filtered_vcf_prefix = outputdir + '/' + sample_t_id + '_PGx'
    filtered_vcf = filtered_vcf_prefix + '.recode.vcf'
    # Check if output vcf does not already exist
    if os.path.exists(filtered_vcf):
        raise IOError("Temporary VCF file " + filtered_vcf + " already exists. Exiting.")
    subprocess.run([vcftools, '--gzvcf', vcf, '--bed', bed_file, '--out', filtered_vcf_prefix,
                    '--indv', sample_r_id, '--recode', '--recode-INFO-all'])
    print("[INFO] Subprocess completed.")
    return filtered_vcf


def load_panel_from_json(panel: str) -> Panel:
    """ Load manually annotated JSON panel file """
    try:
        with open(panel, 'r+', encoding='utf-8') as json_file:
            data = json.load(json_file)
            return Panel.from_json(data)
    except IOError:
        sys.exit("[ERROR] File " + panel + " not found or cannot be opened.")


def get_bed_file(panel_path: str, recreate_bed: bool, panel: Panel, sourcedir: str):
    bed_file = replace_file_extension_of_path(panel_path, "bed")
    if recreate_bed:
        create_bed_file(panel.get_genes(), panel_path, sourcedir, bed_file)
    if not os.path.exists(bed_file):
        sys.exit(
            "[ERROR] Could not locate bed-file. "
            "Could it be that it should be (re)created? "
            "Retry running with --recreate_bed."
        )
    return bed_file


def create_bed_file(genes_in_panel: List[str], panel_path: str, sourcedir: str, bed_path: str) -> None:
    """ Generate bed file from gene panel and save as panel_path.bed """
    print("[INFO] Recreating bed-file...")
    header = 'track name="' + panel_path + '" description="Bed file generated from ' + panel_path + \
             ' with HMF_PGx main.py"\n'
    bed_regions = []  # chrom, start, end, gene
    covered = []
    transcripts = open(sourcedir + "/all_genes.37.tsv", 'r')
    for line in transcripts:
        split_line = line.rstrip().split("\t")
        if split_line[4] in genes_in_panel:
            if split_line[4] not in covered:
                bed_regions.append([split_line[0], split_line[1], split_line[2], split_line[4]])
                covered.append(split_line[4])
    if set(covered) != set(genes_in_panel):
        raise ValueError("[ERROR] Missing genes from the gene panel in the transcript list. Please check:\nCovered:\n"
                         + str(covered) + "\nOriginal gene panel:\n" + str(genes_in_panel))

    with open(bed_path, 'w') as bed:
        bed.write(header)
        for entry in bed_regions:
            bed.write("\t".join(entry) + "\n")

    print("[INFO] Created " + bed_path)


def print_calls_to_file(calls_file: str, all_ids_in_panel) -> None:
    all_ids_in_panel.to_csv(calls_file, sep='\t', index=False)


def print_haplotypes_to_file(genotype_file: str, drug_info, panel, results, severity, version):
    with open(genotype_file, 'w') as f:
        f.write("gene\thaplotype\tfunction\tlinked_drugs\turl_prescription_info\tpanel_version\trepo_version\n")
        for gene in results:
            for haplotype in results[gene]:
                f.write(
                    gene + "\t" +
                    haplotype + "\t" +
                    severity[haplotype.split("_")[0]] + "\t" +
                    drug_info[gene][0] + "\t" +
                    drug_info[gene][1] + "\t" +
                    panel + "\t" +
                    version + "\n"
                )


def parse_args(sys_args: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=('Run pharmacogenomics panel on germline VCF file. The pharmacogenomic annotations are done on '
                     'GRCh38, so in the output both reference genome output is given where possible.')
    )
    parser.add_argument('vcf', type=str, help='VCF file to use for pharmacogenomics analysis')
    parser.add_argument('sample_t_id', type=str, help='The sample ID of the tumor')
    parser.add_argument('sample_r_id', type=str, help='The sample ID of the normal')
    parser.add_argument('version', type=str, help='The version of the tool')
    parser.add_argument('outputdir', type=str, help='Directory to store output of pharmacogenomic analysis')
    parser.add_argument('panel', type=str, help='Json file with the panel variants')
    parser.add_argument('vcftools', type=str, default='vcftools', help="Path to vcftools > 0.1.14 if not in $PATH")
    parser.add_argument(
        '--recreate_bed', default=False, action='store_true',
        help='Recreate bed-file from JSON files. If false, the panel file with extension .bed is searched for.'
    )
    parser.add_argument('--sourcedir', type=str, default='data', help="Optional path to location of source files")
    return parser.parse_args(sys_args)


def replace_file_extension_of_path(path: str, new_file_extension: str):
    split_path = path.split(".")
    new_path = ".".join(split_path[0:-1]) + "." + new_file_extension
    return new_path


def compare_collection(a, b):
    if collections.Counter(a) == collections.Counter(b):
        return True
    else:
        return False


if __name__ == '__main__':
    args = parse_args(sys.argv[1:])
    main(args.vcf, args.sample_t_id, args.sample_r_id, args.version, args.panel,
         args.outputdir, args.recreate_bed, args.vcftools, args.sourcedir)
