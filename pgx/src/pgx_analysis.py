import collections
import itertools
import sys
from typing import DefaultDict, Dict, Any, List, Tuple, Optional

import pandas as pd

from call_data import Grch37CallData
from config.panel import Panel
from config.rs_id_info import RsIdInfo


def create_pgx_analysis(call_data: Grch37CallData, panel: Panel) -> Tuple[Dict[str, List[str]], pd.DataFrame]:
    ids_found_in_patient_df = process_differences_in_ref_sequence(call_data, panel)

    assert_compliance_to_panel_and_fill_in_missing_values(ids_found_in_patient_df, panel)

    if pd.isna(ids_found_in_patient_df).any(axis=None):
        raise ValueError(f"Unhandled NaN values:\n{ids_found_in_patient_df}")

    # Generate a list of ids not found in patient
    rs_ids_found_in_patient = set(ids_found_in_patient_df.rsid.tolist())
    positions_found_in_patient = set(ids_found_in_patient_df.position_GRCh37.tolist())

    rsid_to_gene_to_rs_id_info: DefaultDict[str, Dict[str, Any]] = collections.defaultdict(dict)
    for gene_info in panel.get_gene_infos():
        for rs_id_info in gene_info.rs_id_infos:
            rsid_to_gene_to_rs_id_info[rs_id_info.rs_id][gene_info.gene] = rs_id_info

    ids_not_found_in_patient = pd.DataFrame(columns=['position_GRCh37', 'ref_GRCh37', 'alt_GRCh37', 'position_GRCh38',
                                                     'ref_GRCh38', 'alt_GRCh38', 'rsid', 'variant_annotation', 'gene',
                                                     'filter'])
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
    all_ids_in_panel = pd.concat([ids_found_in_patient_df, ids_not_found_in_patient], sort=True)
    all_ids_in_panel = all_ids_in_panel.sort_values(by='position_GRCh37').reset_index(drop=True)
    all_ids_in_panel = all_ids_in_panel[['gene', 'position_GRCh37', 'ref_GRCh37', 'alt_GRCh37', 'position_GRCh38',
                                         'ref_GRCh38', 'alt_GRCh38', 'rsid', 'variant_annotation', 'filter']]

    results = {}

    for gene_info in panel.get_gene_infos():
        print("[INFO] PROCESSING GENE " + gene_info.gene)
        ids_found_in_gene = all_ids_in_panel.loc[all_ids_in_panel['gene'] == gene_info.gene]
        perfect_match = False

        # If all variants are assumed_ref, return reference haplotype
        if len(ids_found_in_gene.loc[ids_found_in_gene['variant_annotation'] == "REF_CALL"]) == len(ids_found_in_gene):
            print("[INFO] Found reference haplotype")
            results[gene_info.gene] = [gene_info.reference_haplotype_name + "_HOM"]
        else:
            results[gene_info.gene] = []
            haplotypes_matching = []
            for haplotype in gene_info.haplotypes:
                vars_found_in_gene = ids_found_in_gene.loc[ids_found_in_gene['variant_annotation'] != "REF_CALL"]
                variants_sample = list(zip(vars_found_in_gene.rsid.tolist(), vars_found_in_gene.alt_GRCh38.tolist()))
                if set(variants_sample) == set([(x.rs_id, x.variant_allele) for x in haplotype.variants]):
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
                    if set([(x.rs_id, x.variant_allele) for x in haplotype.variants]) <= \
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
                                    rs_ids_subset.append((var.rs_id, var.variant_allele))
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

    return results, all_ids_in_panel


def assert_compliance_to_panel_and_fill_in_missing_values(ids_found_in_patient_df: pd.DataFrame, panel: Panel) -> None:
    if not ids_found_in_patient_df.empty:
        for index, row in ids_found_in_patient_df.iterrows():
            # Fill in empty alleles GRCh38
            if 'ref_GRCh38' not in row or pd.isna(row['ref_GRCh38']):
                ids_found_in_patient_df.at[index, 'ref_GRCh38'] = row['ref_GRCh37']
                ids_found_in_patient_df.at[index, 'alt_GRCh38'] = row['alt_GRCh37']

            # Check against panel and fill in empty fields when possible
            if panel.contains_rs_id_with_position(row['position_GRCh37']):
                panel_rs_id_info = panel.get_rs_id_info_with_position(row['position_GRCh37'])

                if row['rsid'] == '.':
                    ids_found_in_patient_df.at[index, 'rsid'] = panel_rs_id_info.rs_id
                elif row['rsid'] != panel_rs_id_info.rs_id:
                    error_msg = (
                        f"Rs id from input file matches different rs id from panel based on GRCh37 position\n:"
                        f"position GRCh37={row['position_GRCh37']}, rs id input file={row['rsid']}, "
                        f"rs id from panel={panel_rs_id_info.rs_id}"
                    )
                    raise ValueError(error_msg)

                panel_grch38_position = panel_rs_id_info.start_coordinate_grch38.get_position_string()
                if 'position_GRCh38' not in row or pd.isna(row['position_GRCh38']):
                    ids_found_in_patient_df.at[index, 'position_GRCh38'] = panel_grch38_position
                elif row['position_GRCh38'] != panel_grch38_position:
                    error_msg = (
                        f"[ERROR] Inconsistent GRCh38 locations for variants:\n" 
                        f"location from exceptions: {row['position_GRCh38']}\n" 
                        f"location from rs id info: {panel_grch38_position}\n"
                    )
                    raise ValueError(error_msg)
            elif pd.isna(row['rsid']) or not panel.contains_rs_id(row['rsid']):
                # No matching known rs ids
                if 'position_GRCh38' not in row.index or pd.isna(row['position_GRCh38']):
                    ids_found_in_patient_df.at[index, 'position_GRCh38'] = "UNKNOWN"
            else:
                error_msg = (
                    f"[ERROR] Match rs id info from panel on rs id but not position:\n" 
                    f"rs id: {row['rsid']}, input file position: {row['position_GRCh37']}"
                )
                raise ValueError(error_msg)


def process_differences_in_ref_sequence(ids_found: Grch37CallData, panel: Panel) -> pd.DataFrame:
    ids_found_df = ids_found.get_data_frame()
    for rs_id_info, gene, annotation in panel.get_ref_seq_differences():
        variant_location_grch37 = rs_id_info.start_coordinate_grch37.get_position_string()
        variant_location_grch38 = rs_id_info.start_coordinate_grch38.get_position_string()
        ref_allele_grch37 = rs_id_info.reference_allele_grch37
        ref_allele_grch38 = rs_id_info.reference_allele_grch38

        found_var = get_matching_variant_in_patient(ids_found_df, rs_id_info)

        if found_var is not None:
            found_ref_allele = found_var['ref_GRCh37'].iat[0]
            found_alt_allele = found_var['alt_GRCh37'].iat[0]
            if found_ref_allele == ref_allele_grch38 and found_alt_allele == ref_allele_grch38:
                # Delete found variant from results if the ref base and alt base are the correctedRefBase
                ids_found_df = ids_found_df.drop(found_var.index[0])
            elif found_ref_allele == ref_allele_grch38 and found_alt_allele == ref_allele_grch37:
                raise NotImplementedError("What should we do? ref = corRef, alt = corAlt")
            elif found_ref_allele == ref_allele_grch37 and found_alt_allele == ref_allele_grch38:
                # Change variant_annotation and ref and alt base
                ids_found_df.at[found_var.index[0], 'variant_annotation'] = annotation
                ids_found_df.at[found_var.index[0], 'ref_GRCh38'] = ref_allele_grch38
                ids_found_df.at[found_var.index[0], 'alt_GRCh38'] = ref_allele_grch37
                ids_found_df.at[found_var.index[0], 'position_GRCh38'] = variant_location_grch38
            elif found_ref_allele == ref_allele_grch37 and found_alt_allele == ref_allele_grch37:
                # Add variant_annotation and ref and alt base for hg38
                ids_found_df.at[found_var.index[0], 'variant_annotation'] = annotation
                ids_found_df.at[found_var.index[0], 'ref_GRCh38'] = ref_allele_grch37
                ids_found_df.at[found_var.index[0], 'alt_GRCh38'] = ref_allele_grch37
                ids_found_df.at[found_var.index[0], 'position_GRCh38'] = variant_location_grch38
            else:
                # Not completely sure what should happen, so try something and print warning
                new_annotation = found_var["variant_annotation"].iat[0] + "?"
                ids_found_df.at[found_var.index[0], 'variant_annotation'] = new_annotation
                if found_alt_allele == ref_allele_grch38:
                    ids_found_df.at[found_var.index[0], 'ref_GRCh38'] = found_alt_allele
                    ids_found_df.at[found_var.index[0], 'alt_GRCh38'] = found_ref_allele
                else:
                    ids_found_df.at[found_var.index[0], 'ref_GRCh38'] = found_ref_allele
                    ids_found_df.at[found_var.index[0], 'alt_GRCh38'] = found_alt_allele
                ids_found_df.at[found_var.index[0], 'position_GRCh38'] = variant_location_grch38
                print(
                    f"[WARN] Unexpected allele in ref seq difference. Check whether annotation is correct: "
                    f"rs id info={rs_id_info}, found alleles={[found_ref_allele, found_alt_allele]}, "
                    f"annotation={new_annotation}"
                )
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
            ids_found_df = ids_found_df.append(new_id, ignore_index=True)

    return ids_found_df


def get_matching_variant_in_patient(ids_found_df: pd.DataFrame, rs_id_info: RsIdInfo) -> Optional[pd.DataFrame]:
    variant_location_grch37 = rs_id_info.start_coordinate_grch37.get_position_string()
    if rs_id_info.rs_id in ids_found_df['rsid'].tolist():
        found_var_by_rs_id = ids_found_df.loc[ids_found_df['rsid'] == rs_id_info.rs_id]
    else:
        found_var_by_rs_id = None

    if variant_location_grch37 in ids_found_df['position_GRCh37'].tolist():
        found_var_by_location = ids_found_df.loc[ids_found_df['position_GRCh37'] == variant_location_grch37]
    else:
        found_var_by_location = None

    # TODO: (maybe) make this smart enough to also handle MNV's etc.

    if found_var_by_rs_id is not None and found_var_by_location is not None:
        if found_var_by_rs_id.equals(found_var_by_location):
            found_var = found_var_by_rs_id
        else:
            error_msg = f"Rs id and position match with different variants: rs_id_info={rs_id_info}."
            raise ValueError(error_msg)
    elif found_var_by_rs_id is not None and found_var_by_location is None:
        # found_var = found_var_by_rs_id
        error_msg = (f"[WARN] Matched ref seq difference by rs id, but not position: "
                     f"variant_location_grch37={variant_location_grch37}, rs id={rs_id_info.rs_id}")
        raise ValueError(error_msg)
    elif found_var_by_rs_id is None and found_var_by_location is not None:
        found_var = found_var_by_location
        print(f"[WARN] Matched ref seq difference by position, but not rs id: "
              f"variant_location_grch37={variant_location_grch37}, rs id={rs_id_info.rs_id}")
    else:
        # found_var_by_rs_id is None and found_var_by_location is None
        found_var = None

    if found_var is not None and found_var.shape[0] != 1:
        error_msg = f"Multiple variants match with rs id info: rs_id_info={rs_id_info}"
        raise ValueError(error_msg)

    return found_var


def compare_collection(a: List[Tuple[str, str]], b: List[Tuple[str, str]]) -> bool:
    if collections.Counter(a) == collections.Counter(b):
        return True
    else:
        return False
