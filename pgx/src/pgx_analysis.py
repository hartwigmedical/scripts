import collections
import itertools
import sys
from typing import DefaultDict, Dict, Any, List

import pandas as pd

from config.panel import Panel


def create_pgx_analysis(ids_found_in_patient: pd.DataFrame, panel: Panel):
    # Process the differences between GRCh37 and GRCh38
    ids_found_in_patient = process_differences_in_ref_sequence(ids_found_in_patient, panel)

    rsid_to_gene_to_rs_id_info: DefaultDict[str, Dict[str, Any]] = collections.defaultdict(dict)
    for gene_info in panel.get_gene_infos():
        for rs_id_info in gene_info.rs_id_infos:
            rsid_to_gene_to_rs_id_info[rs_id_info.rs_id][gene_info.gene] = rs_id_info

    results = {}
    severity = {}
    drug_info = {}

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
                if panel.contains_rs_id_with_position(position_string):
                    rs_id = panel.get_rs_id_with_position(position_string)
                else:
                    rs_id = '.'
                ids_found_in_patient.at[index, 'rsid'] = rs_id

        for index, row in ids_found_in_patient.iterrows():
            matching_rs_id_infos = rsid_to_gene_to_rs_id_info[ids_found_in_patient.at[index, 'rsid']].values()
            grch38_locations = {
                rs_id_info.start_coordinate_grch38.get_position_string() for rs_id_info in matching_rs_id_infos
            }
            if len(grch38_locations) == 1:
                grch38_location = grch38_locations.pop()
                if 'position_GRCh38' not in row.index or pd.isna(row['position_GRCh38']):
                    ids_found_in_patient.at[index, 'position_GRCh38'] = grch38_location
                elif row['position_GRCh38'] != grch38_location:
                    raise ValueError(
                        "[ERROR] Inconsistent GRCh38 locations for variants:\n"
                        "location from exceptions: " + row['position_GRCh38'] + "\n"
                        "location from rs id info: " + grch38_location + "\n"
                    )
            elif len(grch38_locations) > 1:
                matching_rs_id_infos_string = ",".join([str(rs_id_info) for rs_id_info in matching_rs_id_infos])
                raise ValueError("[ERROR] Inconsistent GRCh38 locations for variants:\n"
                                 "matching rs id infos: " + matching_rs_id_infos_string + "\n"
                                 "GRCh38 locations: " + ", ".join(grch38_locations) + "\n")

    # Generate a list of ids not found in patient
    ids_not_found_in_patient = pd.DataFrame(columns=['position_GRCh37', 'ref_GRCh37', 'alt_GRCh37', 'position_GRCh38',
                                                     'ref_GRCh38', 'alt_GRCh38', 'rsid', 'variant_annotation', 'gene',
                                                     'filter'])

    rs_ids_found_in_patient = set(ids_found_in_patient.rsid.tolist())
    positions_found_in_patient = set(ids_found_in_patient.position_GRCh37.tolist())

    for rs_id_info in panel.get_rs_id_infos():
        if (rs_id_info.rs_id not in rs_ids_found_in_patient and
                rs_id_info.start_coordinate_grch37.get_position_string() not in positions_found_in_patient):
            # TODO: check whether this properly takes variants of more than one base pair, so MNV's etc., into account.
            #       The fact that only a single position is checked is suspicious.
            new_id = {}
            for gene, rs_id_info in rsid_to_gene_to_rs_id_info[rs_id_info.rs_id].items():
                new_id['position_GRCh37'] = rs_id_info.start_coordinate_grch37.get_position_string()
                new_id['rsid'] = rs_id_info.rs_id
                new_id['ref_GRCh37'] = rs_id_info.reference_allele_grch38
                new_id['alt_GRCh37'] = rs_id_info.reference_allele_grch38  # Assuming REF/REF relative to GRCh38
                new_id['variant_annotation'] = "REF_CALL"
                new_id['filter'] = "NO_CALL"
                new_id['gene'] = gene
                new_id['ref_GRCh38'] = rs_id_info.reference_allele_grch38  # Again assuming REF/REF relative to GRCh38
                new_id['alt_GRCh38'] = rs_id_info.reference_allele_grch38
                new_id['position_GRCh38'] = rs_id_info.start_coordinate_grch38.get_position_string()
                ids_not_found_in_patient = ids_not_found_in_patient.append(new_id, ignore_index=True)

    # Now we want to process all the variants in terms of the alleles
    all_ids_in_panel = pd.concat([ids_found_in_patient, ids_not_found_in_patient], sort=True)
    all_ids_in_panel = all_ids_in_panel.sort_values(by='position_GRCh37').reset_index(drop=True)
    all_ids_in_panel = all_ids_in_panel[['gene', 'position_GRCh37', 'ref_GRCh37', 'alt_GRCh37', 'position_GRCh38',
                                         'ref_GRCh38', 'alt_GRCh38', 'rsid', 'variant_annotation', 'filter']]

    for gene_info in panel.get_gene_infos():
        print("[INFO] PROCESSING GENE " + gene_info.gene)
        ids_found_in_gene = all_ids_in_panel[all_ids_in_panel['gene'].str.contains(gene_info.gene)]
        perfect_match = False
        severity[gene_info.reference_haplotype_name] = "Normal Function"
        severity['Unresolved'] = "Unknown Function"
        drug_info[gene_info.gene] = [
            ";".join([drug.name for drug in gene_info.drugs]),
            ";".join([drug.url_prescription_info for drug in gene_info.drugs])
        ]

        # If all variants are assumed_ref, return reference haplotype
        if len(ids_found_in_gene.loc[ids_found_in_gene['variant_annotation'] == "REF_CALL"]) == len(ids_found_in_gene):
            print("[INFO] Found reference haplotype")
            results[gene_info.gene] = [gene_info.reference_haplotype_name + "_HOM"]
        else:
            results[gene_info.gene] = []
            haplotypes_matching = []
            for haplotype in gene_info.haplotypes:
                severity[haplotype.name] = haplotype.function
                vars_found_in_gene = ids_found_in_gene.loc[ids_found_in_gene['variant_annotation'] != "REF_CALL"]
                variants_sample = list(zip(vars_found_in_gene.rsid.tolist(), vars_found_in_gene.alt_GRCh38.tolist()))
                if set(variants_sample) == set([(x.rs_id, x.alt_allele_grch38) for x in haplotype.variants]):
                    perfect_match = True
                    print("[INFO] Found 1:1 match with haplotype " + haplotype.name)
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
                    results[gene_info.gene].append(haplotype.name + "_" + str(allele_status))
                    if allele_status == "HET":
                        # Assume if perfect match with HET, we are also looking at reference haplotype
                        results[gene_info.gene].append(gene_info.reference_haplotype_name + "_HET")
                    break
                else:
                    # print("Processing " + str(haplotype['alleleName']))
                    # print(set([(x['rsid'], x['altAlleleGRCh38']) for x in haplotype['alleleVariants']]))
                    # print(set(variants_sample))
                    if set([(x.rs_id, x.alt_allele_grch38) for x in haplotype.variants]) <= \
                            set(variants_sample):
                        print("[INFO] A subset of rsids matches " + str(haplotype.name) + " in part")
                        haplotypes_matching.append(haplotype)

            if not perfect_match:
                if not haplotypes_matching:
                    print(
                        f"[WARN] No haplotype match found for {gene_info.gene}. "
                        f"Probable cause is that the rs_id_info is not in line with previously "
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
                            allele_variants = [x.variants for x in subset]
                            rs_ids_subset = []
                            for x in allele_variants:
                                for var in x:
                                    rs_ids_subset.append((var.rs_id, var.alt_allele_grch38))
                            if compare_collection(variants_sample, rs_ids_subset):
                                print("[INFO] Perfect haplotype combination found!")
                                perfect_match = True
                                for haplotype in subset:
                                    allele_statuses = []
                                    rs_ids_in_allele = [x.rs_id for x in haplotype.variants]
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
                                    results[gene_info.gene].append(haplotype.name + "_" + str(allele_status))
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
                            # If we have a rs_id_info that is in sample or in haplotype that is not matched > undetermined
                            if len(optimal_set[0][3]) > 0 or len(optimal_set[0][4]) > 0:
                                results[gene_info.gene].append("Unresolved_Haplotype")
                            else:
                                for haplotype in subset:
                                    allele_statuses = []
                                    rs_ids_in_allele = [x.rs_id for x in haplotype.variants]
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
                                    results[gene_info.gene].append(haplotype.name + "_" + str(allele_status))
                        else:
                            sys.exit("[ERROR] No haplotype match was found. Exiting.")
        # If we only find one haplotype and it is HET, assume we're also dealing with reference haplotype
        if len(results[gene_info.gene]) == 1:
            if results[gene_info.gene][0].split("_")[-1] == "HET":
                results[gene_info.gene].append(gene_info.reference_haplotype_name + "_HET")

    return results, severity, all_ids_in_panel, drug_info


def process_differences_in_ref_sequence(ids_found: pd.DataFrame, panel: Panel) -> pd.DataFrame:
    for rs_id_info, gene, annotation in panel.get_ref_seq_differences():
        variant_location_grch37 = rs_id_info.start_coordinate_grch37.get_position_string()
        variant_location_grch38 = rs_id_info.start_coordinate_grch38.get_position_string()
        ref_allele_grch37 = rs_id_info.reference_allele_grch37
        ref_allele_grch38 = rs_id_info.reference_allele_grch38

        variant_location_found = variant_location_grch37 in ids_found['position_GRCh37'].tolist()
        rs_id_found = rs_id_info.rs_id in ids_found['rsid'].tolist()
        if rs_id_found or variant_location_found:
            if rs_id_found:
                # get line and index from ids_found
                found_var = ids_found[ids_found['rsid'].str.contains(rs_id_info.rs_id)]
            else:
                found_var = ids_found[ids_found['position_GRCh37'].str.contains(variant_location_grch37)]
            # Delete found variant from results if the ref base and alt base are the correctedRefBase

            found_ref_allele = found_var['ref_GRCh37'].values
            found_alt_allele = found_var['alt_GRCh37'].values
            if found_ref_allele == ref_allele_grch38 and found_alt_allele == ref_allele_grch38:
                ids_found = ids_found.drop(found_var.index[0])
            elif found_ref_allele == ref_allele_grch38 and found_alt_allele == ref_allele_grch37:
                raise NotImplementedError("What should we do? ref = corRef, alt = corAlt")
            elif found_ref_allele == ref_allele_grch37 and found_alt_allele == ref_allele_grch38:
                # Change variant_annotation and ref and alt base
                ids_found.at[found_var.index[0], 'variant_annotation'] = annotation
                ids_found.at[found_var.index[0], 'ref_GRCh38'] = ref_allele_grch38
                ids_found.at[found_var.index[0], 'alt_GRCh38'] = ref_allele_grch37
                ids_found.at[found_var.index[0], 'position_GRCh38'] = variant_location_grch38
            elif found_ref_allele == ref_allele_grch37 and found_alt_allele == ref_allele_grch37:
                # Add variant_annotation and ref and alt base for hg38
                ids_found.at[found_var.index[0], 'variant_annotation'] = annotation
                ids_found.at[found_var.index[0], 'ref_GRCh38'] = ref_allele_grch37
                ids_found.at[found_var.index[0], 'alt_GRCh38'] = ref_allele_grch37
                ids_found.at[found_var.index[0], 'position_GRCh38'] = variant_location_grch38
            else:
                print("[ERROR] Complete mismatch:")
                print(found_var)
                raise ValueError("[ERROR] Exceptions cannot be processed. Please check. Exiting.")
        else:
            print("[INFO] Exception variant not found in this patient. This means ref/ref call, but should be "
                  "flipped. Add to table.")
            new_id = {'position_GRCh37': variant_location_grch37,
                      'rsid': rs_id_info.rs_id,
                      'ref_GRCh37': ref_allele_grch37,
                      'alt_GRCh37': ref_allele_grch37,
                      'variant_annotation': annotation,
                      'filter': "INFERRED_REF_CALL",
                      'gene': gene,
                      'position_GRCh38': variant_location_grch38,
                      'ref_GRCh38': ref_allele_grch37,
                      'alt_GRCh38': ref_allele_grch37}
            ids_found = ids_found.append(new_id, ignore_index=True)

    return ids_found


def compare_collection(a, b) -> bool:
    if collections.Counter(a) == collections.Counter(b):
        return True
    else:
        return False
