import collections
import itertools
import sys
from typing import Dict, Any, List, Tuple, Optional, Set

import pandas as pd

from base.gene_coordinate import GeneCoordinate
from call_data import Grch37CallData, FullCall, AnnotatedAllele, Grch37Call
from config.gene_info import GeneInfo
from config.panel import Panel
from config.rs_id_info import RsIdInfo
from config.variant import Variant
from dataframe_format import COMBINED_DATAFRAME_COLUMNS


def create_pgx_analysis(call_data: Grch37CallData, panel: Panel) -> Tuple[Dict[str, Set[str]], pd.DataFrame]:
    panel_calls_found_in_patient_df = get_panel_calls_found_in_patient_df(call_data, panel)

    if pd.isna(panel_calls_found_in_patient_df).any(axis=None):
        raise ValueError(f"This should never happen: Unhandled NaN values:\n{panel_calls_found_in_patient_df}")

    # Now we want to process all the variants in terms of the alleles
    results = {}
    for gene_info in panel.get_gene_infos():
        print("[INFO] PROCESSING GENE " + gene_info.gene)
        # TODO: is this safe? to filter by gene like this? what if the name is different?
        #       then again, how else can I determine what variants belong to what gene?
        #       Maybe get gene chromosome and positions from bed file and use them to filter.
        #       But then the format of the dataframe is annoying for filtering,
        #       because chrom and position are combined in a single column
        ids_found_in_gene = panel_calls_found_in_patient_df.loc[panel_calls_found_in_patient_df['gene'] == gene_info.gene]

        results[gene_info.gene] = get_haplotypes(ids_found_in_gene, gene_info)

    panel_calls_for_patient = get_panel_calls_for_patient(panel_calls_found_in_patient_df, panel)

    return results, panel_calls_for_patient


# def get_haplotypes_alternative(ids_found_in_gene: pd.DataFrame, gene_info: GeneInfo) -> Set[str]:
#     variant_to_count = {}
#     for row in ids_found_in_gene.iterrows():
#         # Note that row.rs_id can include multiple rs ids separated by ;
#         # TODO: make solution cleaner
#         variant_to_count[Variant(row["rsid"], Variant)]


def get_haplotypes(ids_found_in_gene: pd.DataFrame, gene_info: GeneInfo) -> Set[str]:
    perfect_match = False
    if len(ids_found_in_gene.loc[ids_found_in_gene['variant_annotation'] == "REF_CALL"]) == len(ids_found_in_gene):
        # If all variants are assumed_ref, return reference haplotype
        print("[INFO] Found reference haplotype")
        haplotypes = {gene_info.reference_haplotype_name + "_HOM"}
    else:
        haplotypes = set()
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
                haplotypes.add(haplotype.name + "_" + str(allele_status))
                if allele_status == "HET":
                    # Assume if perfect match with HET, we are also looking at reference haplotype
                    haplotypes.add(gene_info.reference_haplotype_name + "_HET")
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
                haplotypes = {"Unresolved_Haplotype"}
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
                                haplotypes.add(haplotype.name + "_" + str(allele_status))
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
                            haplotypes = {"Unresolved_Haplotype"}
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
                                haplotypes.add(haplotype.name + "_" + str(allele_status))
                    else:
                        sys.exit("[ERROR] No haplotype match was found. Exiting.")
    # If we only find one haplotype and it is HET, assume we're also dealing with reference haplotype
    if len(haplotypes) == 1:
        (single_haplotype,) = haplotypes
        if single_haplotype.split("_")[-1] == "HET":
            haplotypes.add(gene_info.reference_haplotype_name + "_HET")
    return haplotypes


def get_panel_calls_for_patient(panel_calls_found_in_patient_df: pd.DataFrame, panel: Panel) -> pd.DataFrame:
    panel_calls_not_found_in_patient_df_df = get_calls_not_found_in_patient_df(panel_calls_found_in_patient_df, panel)
    panel_calls_for_patient = pd.concat(
        [panel_calls_found_in_patient_df, panel_calls_not_found_in_patient_df_df], sort=True
    )
    panel_calls_for_patient = panel_calls_for_patient.sort_values(by='position_GRCh37').reset_index(drop=True)
    panel_calls_for_patient = panel_calls_for_patient[list(COMBINED_DATAFRAME_COLUMNS)]
    return panel_calls_for_patient


def get_calls_not_found_in_patient_df(panel_calls_found_in_patient_df: pd.DataFrame, panel: Panel) -> pd.DataFrame:
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


def get_panel_calls_found_in_patient_df(grch37_call_data: Grch37CallData, panel: Panel) -> pd.DataFrame:
    # panel_calls_found_in_patient_df = process_differences_in_ref_sequence(call_data, panel)

    handled_grch37_start_coordinates: Set[GeneCoordinate] = set()
    handled_grch37_rs_ids: Set[str] = set()
    full_calls_from_grch37_calls = []

    for grch37_call in grch37_call_data.calls:
        full_call = get_full_call_from_grch37_call(grch37_call, panel)

        for rs_id in full_call.rs_ids:
            if rs_id != ".":
                if rs_id in handled_grch37_rs_ids:
                    error_msg = (f"Call for rs id that has already been handled:\n"
                                 f"call={grch37_call}\n"
                                 f"handled_rs_ids={handled_grch37_rs_ids}")
                    raise ValueError(error_msg)
                handled_grch37_rs_ids.add(rs_id)

        if full_call.start_coordinate_grch37 in handled_grch37_start_coordinates:
            error_msg = (f"Call for position that has already been handled:\n"
                         f"call={grch37_call}\n"
                         f"handled_coords={handled_grch37_start_coordinates}")
            raise ValueError(error_msg)
        handled_grch37_start_coordinates.add(full_call.start_coordinate_grch37)
        full_calls_from_grch37_calls.append(full_call)

    # The "inferred ref" calls are calls for ref seq differences that do not correspond to grch37 calls.
    # This means that they were ref vs GRCh37 and therefore variant calls vs GRCh38
    inferred_ref_calls = get_inferred_ref_calls(panel, handled_grch37_start_coordinates, handled_grch37_rs_ids)

    full_calls = full_calls_from_grch37_calls + inferred_ref_calls

    # make pandas dataframe
    panel_calls_found_in_patient_df_alt = pd.DataFrame(columns=COMBINED_DATAFRAME_COLUMNS)
    for full_call in full_calls:
        grch37_alleles = [
            annotated.allele for annotated in full_call.annotated_alleles if not annotated.is_variant_vs_grch37
        ] + [
            annotated.allele for annotated in full_call.annotated_alleles if annotated.is_variant_vs_grch37
        ]
        grch38_alleles = [
            annotated.allele for annotated in full_call.annotated_alleles if not annotated.is_variant_vs_grch38
        ] + [
            annotated.allele for annotated in full_call.annotated_alleles if annotated.is_variant_vs_grch38
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
        panel_calls_found_in_patient_df_alt = panel_calls_found_in_patient_df_alt.append(new_id, ignore_index=True)

    return panel_calls_found_in_patient_df_alt


def get_inferred_ref_calls(panel: Panel, handled_grch37_start_coordinates: Set[GeneCoordinate],
                           handled_grch37_rs_ids: Set[str]) -> List[FullCall]:
    # The "inferred ref" calls are calls for ref seq differences that do not correspond to grch37 calls.
    # This means that they were ref vs GRCh37 and therefore variant calls vs GRCh38
    inferred_ref_calls = []
    for rs_id_info, gene, annotation in panel.get_ref_seq_differences():
        if rs_id_info.start_coordinate_grch37 not in handled_grch37_start_coordinates:
            # TODO: rename INFERRED_REF_CALL to INFERRED_CALL or something
            if rs_id_info.rs_id in handled_grch37_rs_ids:
                error_msg = (f"Have seen rs id of ref seq difference, but not location. "
                             f"Indicates mismatch between input file and panel."
                             f"rs_id_info={rs_id_info}")
                raise ValueError(error_msg)

            annotated_alleles = (
                AnnotatedAllele(rs_id_info.reference_allele_grch37, True, False),
                AnnotatedAllele(rs_id_info.reference_allele_grch37, True, False),
            )
            full_call = FullCall(
                rs_id_info.start_coordinate_grch37,
                rs_id_info.start_coordinate_grch38,
                annotated_alleles,
                gene,
                (rs_id_info.rs_id,),
                annotation,
                "INFERRED_REF_CALL",
            )
            inferred_ref_calls.append(full_call)
    return inferred_ref_calls


def get_full_call_from_grch37_call(grch37_call: Grch37Call, panel: Panel) -> FullCall:
    gene = grch37_call.gene
    assert_gene_in_panel(gene, panel)

    start_coordinate_grch37 = grch37_call.start_coordinate

    # determine annotated_alleles, start_coordinate_grch38 and rs_ids
    rs_ids: Tuple[str, ...]
    start_coordinate_grch38: Optional[GeneCoordinate]
    if panel.contains_rs_id_with_position(start_coordinate_grch37.get_position_string()):
        rs_id_info = panel.get_rs_id_info_with_position(start_coordinate_grch37.get_position_string())
        assert_rs_id_call_matches_info(grch37_call.rs_ids, (rs_id_info.rs_id,))

        annotated_alleles = get_annotated_alleles_from_rs_id_info(rs_id_info, grch37_call)
        assert_alleles_in_expected_order(annotated_alleles)

        start_coordinate_grch38 = rs_id_info.start_coordinate_grch38
        if grch37_call.rs_ids == (".",) and rs_id_info is not None:
            rs_ids = (rs_id_info.rs_id,)
        else:
            rs_ids = grch37_call.rs_ids
    elif any(panel.contains_rs_id(rs_id) for rs_id in grch37_call.rs_ids):
        error_msg = (
            f"[ERROR] Match rs id info from panel on an rs id but not position:\n"
            f"rs ids: {grch37_call.rs_ids}, input file position: {start_coordinate_grch37}"
        )
        raise ValueError(error_msg)
    else:
        # unknown variant
        annotated_alleles = (
            AnnotatedAllele(grch37_call.alleles[0], None, None),
            AnnotatedAllele(grch37_call.alleles[1], None, None)
        )
        start_coordinate_grch38 = None
        rs_ids = grch37_call.rs_ids

    # determine variant annotation and filter
    if len(rs_ids) > 1 and any(panel.has_ref_seq_difference_annotation(gene, rs_id) for rs_id in rs_ids):
        error_msg = f"One of multiple rs ids is of a ref seq difference, so not sure how to annotate {rs_ids}"
        raise ValueError(error_msg)
    elif len(rs_ids) == 1 and panel.has_ref_seq_difference_annotation(gene, rs_ids[0]):
        ref_call_due_to_ref_sequence_difference = all(
            annotated.is_variant_vs_grch37 and not annotated.is_variant_vs_grch38 for annotated in annotated_alleles
        )
        all_variants_ref_to_grch37_or_grch38 = all(
            not annotated.is_variant_vs_grch37 or not annotated.is_variant_vs_grch38 for annotated in annotated_alleles
        )

        if ref_call_due_to_ref_sequence_difference:
            variant_annotation = "REF_CALL"
            filter = "NO_CALL"
        elif all_variants_ref_to_grch37_or_grch38:
            variant_annotation = panel.get_ref_seq_difference_annotation(gene, rs_ids[0])
            filter = grch37_call.filter
        else:
            variant_annotation = grch37_call.variant_annotation + "?"
            filter = grch37_call.filter
            print(
                f"[WARN] Unexpected allele in ref seq difference location. Check whether annotation is correct: "
                f"found alleles=({annotated_alleles[0]}, {annotated_alleles[1]}), "
                f"annotation={variant_annotation}"
            )
    else:
        # no ref seq differences involved
        variant_annotation = grch37_call.variant_annotation
        filter = grch37_call.filter

    full_call = FullCall(
        start_coordinate_grch37, start_coordinate_grch38, annotated_alleles, gene, rs_ids, variant_annotation, filter,
    )
    return full_call


def get_annotated_alleles_from_rs_id_info(
        rs_id_info: RsIdInfo, grch37_call: Grch37Call) -> Tuple[AnnotatedAllele, AnnotatedAllele]:
    return (
        AnnotatedAllele(
            grch37_call.alleles[0],
            grch37_call.alleles[0] != rs_id_info.reference_allele_grch37,
            grch37_call.alleles[0] != rs_id_info.reference_allele_grch38,
            ),
        AnnotatedAllele(
            grch37_call.alleles[1],
            grch37_call.alleles[1] != rs_id_info.reference_allele_grch37,
            grch37_call.alleles[1] != rs_id_info.reference_allele_grch38,
            )
    )


def assert_rs_id_call_matches_info(rs_ids_call: Tuple[str, ...], rs_ids_info: Tuple[str, ...]) -> None:
    if rs_ids_call != (".",) and rs_ids_call != rs_ids_info:
        # TODO: make this more flexible, if necessary
        error_msg = (f"Given rs id does not match rs id from panel: "
                     f"from call={rs_ids_call}, from panel={rs_ids_info}")
        raise ValueError(error_msg)


def assert_gene_in_panel(gene: str, panel: Panel) -> None:
    if gene not in panel.get_genes():
        error_msg = f"Call for unknown gene:\ngene={gene}"
        raise ValueError(error_msg)


def assert_alleles_in_expected_order(annotated_alleles: Tuple[AnnotatedAllele, AnnotatedAllele]) -> None:
    if annotated_alleles[0].is_variant_vs_grch37 and not annotated_alleles[1].is_variant_vs_grch37:
        error_msg = (f"Alleles are in unexpected order, alt before ref: "
                     f"alleles=({annotated_alleles[0]},{annotated_alleles[0]})")
        raise ValueError(error_msg)


def compare_collection(a: List[Tuple[str, str]], b: List[Tuple[str, str]]) -> bool:
    if collections.Counter(a) == collections.Counter(b):
        return True
    else:
        return False
