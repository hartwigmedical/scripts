import collections
from copy import deepcopy
from typing import Dict, List, Tuple, Optional, Set, FrozenSet, DefaultDict

import pandas as pd

from base.gene_coordinate import GeneCoordinate
from call_data import Grch37CallData, FullCall, AnnotatedAllele, Grch37Call
from config.gene_info import GeneInfo
from config.haplotype import Haplotype
from config.panel import Panel
from config.rs_id_info import RsIdInfo
from config.variant import Variant
from dataframe_format import COMBINED_DATAFRAME_COLUMNS


def create_pgx_analysis(call_data: Grch37CallData, panel: Panel) -> Tuple[Dict[str, Set[str]], pd.DataFrame]:
    full_calls = get_full_calls(call_data, panel)

    # Now we want to process all the variants in terms of the alleles
    gene_to_haplotype_calls = {}
    for gene_info in panel.get_gene_infos():
        print("[INFO] PROCESSING GENE " + gene_info.gene)
        gene_to_haplotype_calls[gene_info.gene] = get_haplotypes_call(full_calls, gene_info)

    panel_calls_for_patient = get_panel_calls_for_patient(full_calls, panel)

    if pd.isna(panel_calls_for_patient).any(axis=None):
        raise ValueError(f"This should never happen: Unhandled NaN values:\n{panel_calls_for_patient}")

    return gene_to_haplotype_calls, panel_calls_for_patient


def get_haplotypes_call(full_calls: List[FullCall], gene_info: GeneInfo) -> Set[str]:
    full_calls_for_gene = [call for call in full_calls if call.gene == gene_info.gene]
    try:
        variant_to_count: DefaultDict[Variant, int] = collections.defaultdict(int)
        for call in full_calls_for_gene:
            assert_handleable_call(call)
            rs_id = call.rs_ids[0]
            for annotated_allele in call.annotated_alleles:
                if annotated_allele.is_variant_vs_grch38 is None:
                    error_msg = f"Unknown variant: allele={annotated_allele}"
                    raise ValueError(error_msg)
                if annotated_allele.is_variant_vs_grch38:
                    variant_to_count[Variant(rs_id, annotated_allele.allele)] += 1

        explaining_haplotype_combinations = get_explaining_haplotype_combinations(variant_to_count, gene_info.haplotypes)

        if not explaining_haplotype_combinations:
            error_msg = f"No explaining haplotype combinations"
            raise ValueError(error_msg)

        minimal_explaining_haplotype_combination = get_minimal_haplotype_combination(explaining_haplotype_combinations)

        haplotype_to_count: DefaultDict[str, int] = collections.defaultdict(int)
        for haplotype in minimal_explaining_haplotype_combination:
            haplotype_to_count[haplotype] += 1

        haplotype_calls = set()
        for haplotype, count in haplotype_to_count.items():
            if count == 1:
                haplotype_calls.add(haplotype + "_HET")
            elif count == 2:
                haplotype_calls.add(haplotype + "_HOM")
            else:
                error_msg = f"Impossible count for haplotype: haplotype={haplotype}, count={count}"
                raise ValueError(error_msg)

        called_haplotypes_count = sum(haplotype_to_count.values())

        if called_haplotypes_count == 0:
            haplotype_calls.add(gene_info.reference_haplotype_name + "_HOM")
        elif called_haplotypes_count == 1:
            haplotype_calls.add(gene_info.reference_haplotype_name + "_HET")

        return haplotype_calls

    except ValueError as e:
        print(f"[Error] Cannot resolve haplotype for gene {gene_info.gene}. Error: {e}")
        return {"Unresolved_Haplotype"}


def get_explaining_haplotype_combinations(
        variant_to_count: DefaultDict[Variant, int], haplotypes: FrozenSet[Haplotype]) -> Set[Tuple[str, ...]]:
    """
    Gets combinations of haplotypes that explain all variants in the stated amounts. Uses recursion.
    Always makes sure that the haplotypes in a haplotype combination are ordered alphabetically to
    ensure that each haplotype combination exists only once in the result set.
    """
    if any(count < 0 for count in variant_to_count.values()):
        return set()
    if all(count == 0 for count in variant_to_count.values()):
        return {tuple()}

    result_set = set()
    for haplotype in haplotypes:
        reduced_variant_to_count = deepcopy(variant_to_count)
        for variant in haplotype.variants:
            reduced_variant_to_count[variant] -= 1
            
        combinations_for_reduced_variant_set = get_explaining_haplotype_combinations(
            reduced_variant_to_count, haplotypes
        )
        for combination in combinations_for_reduced_variant_set:
            result_set.add(tuple(sorted(list(combination) + [haplotype.name])))

    return result_set


def get_minimal_haplotype_combination(explaining_haplotype_combinations: Set[Tuple[str, ...]]) -> Tuple[str, ...]:
    min_haplotype_count = min(len(combination) for combination in explaining_haplotype_combinations)
    minimal_explaining_haplotype_combinations = {
        combination for combination in explaining_haplotype_combinations if len(combination) == min_haplotype_count
    }
    if len(minimal_explaining_haplotype_combinations) > 1:
        error_msg = (f"No unique minimal explaining haplotype combination: "
                     f"options={minimal_explaining_haplotype_combinations}")
        raise ValueError(error_msg)
    minimal_explaining_haplotype_combination = minimal_explaining_haplotype_combinations.pop()
    return minimal_explaining_haplotype_combination


def get_panel_calls_for_patient(full_calls: List[FullCall], panel: Panel) -> pd.DataFrame:
    panel_calls_found_in_patient_df = get_dataframe_from_full_calls(full_calls)
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


def get_dataframe_from_full_calls(full_calls: List[FullCall]) -> pd.DataFrame:
    panel_calls_found_in_patient_df = pd.DataFrame(columns=COMBINED_DATAFRAME_COLUMNS)
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
        panel_calls_found_in_patient_df = panel_calls_found_in_patient_df.append(new_id, ignore_index=True)
    return panel_calls_found_in_patient_df


def get_full_calls(grch37_call_data: Grch37CallData, panel: Panel) -> List[FullCall]:
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
    return full_calls


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
                AnnotatedAllele(rs_id_info.reference_allele_grch37, False, True),
                AnnotatedAllele(rs_id_info.reference_allele_grch37, False, True),
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


def assert_handleable_call(call: FullCall) -> None:
    if len(call.rs_ids) > 1:
        error_msg = f"Call has more than one rs id: rs ids={call.rs_ids}, call={call}"
        raise ValueError(error_msg)
    if len(call.rs_ids) < 1:
        error_msg = f"Call has zero rs ids: call={call}"
        raise ValueError(error_msg)
    if call.rs_ids[0] == ".":
        error_msg = f"Call has unknown rs id: call={call}"
        raise ValueError(error_msg)


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
