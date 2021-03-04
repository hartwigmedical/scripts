from typing import List, Set, Tuple, Optional

from base.gene_coordinate import GeneCoordinate
from call_data import Grch37CallData, FullCall, AnnotatedAllele, Grch37Call
from config.panel import Panel
from config.rs_id_info import RsIdInfo


class Grch37CallTranslator(object):
    @classmethod
    def get_full_calls(cls, grch37_call_data: Grch37CallData, panel: Panel) -> Tuple[FullCall, ...]:
        handled_grch37_coordinates: Set[GeneCoordinate] = set()
        handled_grch37_rs_ids: Set[str] = set()
        full_calls_from_grch37_calls = []
        for grch37_call in grch37_call_data.calls:
            full_call = cls.__get_full_call_from_grch37_call(grch37_call, panel)

            for rs_id in full_call.rs_ids:
                if rs_id != ".":
                    if rs_id in handled_grch37_rs_ids:
                        error_msg = (f"Call for rs id that has already been handled:\n"
                                     f"call={grch37_call}\n"
                                     f"handled_rs_ids={handled_grch37_rs_ids}")
                        raise ValueError(error_msg)
                    handled_grch37_rs_ids.add(rs_id)

            relevant_grch37_coordinates = full_call.get_relevant_grch37_coordinates()
            if relevant_grch37_coordinates.intersection(handled_grch37_coordinates):
                error_msg = (f"Call involves at least one position that has already been handled:\n"
                             f"call={grch37_call}\n"
                             f"handled_coords={handled_grch37_coordinates}")
                raise ValueError(error_msg)
            handled_grch37_coordinates.update(relevant_grch37_coordinates)

            full_calls_from_grch37_calls.append(full_call)
        # The "inferred ref" calls are calls for ref seq differences that do not correspond to grch37 calls.
        # This means that they were ref vs GRCh37 and therefore variant calls vs GRCh38
        inferred_ref_calls = cls.__get_inferred_ref_calls(panel, handled_grch37_coordinates, handled_grch37_rs_ids)
        full_calls = full_calls_from_grch37_calls + inferred_ref_calls
        return tuple(full_calls)

    @classmethod
    def __get_full_call_from_grch37_call(cls, grch37_call: Grch37Call, panel: Panel) -> FullCall:
        gene = grch37_call.gene
        cls.__assert_gene_in_panel(gene, panel)

        start_coordinate_grch37 = grch37_call.start_coordinate
        reference_allele_grch37 = grch37_call.ref_allele

        # determine annotated_alleles, start_coordinate_grch38, reference_allele_grch38 and rs_ids
        start_coordinate_grch38: Optional[GeneCoordinate]
        reference_allele_grch38: Optional[str]
        rs_ids: Tuple[str, ...]
        if panel.contains_matching_rs_id_info(start_coordinate_grch37, reference_allele_grch37):
            rs_id_info = panel.get_matching_rs_id_info(start_coordinate_grch37, reference_allele_grch37)
            cls.__assert_rs_id_call_matches_info(grch37_call.rs_ids, (rs_id_info.rs_id,))

            annotated_alleles = cls.__get_annotated_alleles_from_rs_id_info(rs_id_info, grch37_call)
            cls.__assert_alleles_in_expected_order(annotated_alleles)

            start_coordinate_grch38 = rs_id_info.start_coordinate_grch38
            reference_allele_grch38 = rs_id_info.reference_allele_grch38
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
                AnnotatedAllele(grch37_call.alleles[1], None, None),
            )
            start_coordinate_grch38 = None
            reference_allele_grch38 = None
            rs_ids = grch37_call.rs_ids

        # determine variant annotation and filter
        if len(rs_ids) == 1 and panel.has_ref_seq_difference_annotation(gene, rs_ids[0]):
            ref_call_due_to_ref_sequence_difference = all(
                annotated.is_annotated() and annotated.is_variant_vs_grch37 and not annotated.is_variant_vs_grch38
                for annotated in annotated_alleles
            )
            all_variants_ref_to_grch37_or_grch38 = all(
                annotated.is_annotated() and not annotated.is_variant_vs_grch37 or not annotated.is_variant_vs_grch38
                for annotated in annotated_alleles
            )

            if ref_call_due_to_ref_sequence_difference:
                variant_annotation = "REF_CALL"
                filter_type = "NO_CALL"
            elif all_variants_ref_to_grch37_or_grch38:
                variant_annotation = panel.get_ref_seq_difference_annotation(gene, rs_ids[0])
                filter_type = grch37_call.filter
            else:
                variant_annotation = grch37_call.variant_annotation + "?"
                filter_type = grch37_call.filter
                print(
                    f"[WARN] Unexpected allele in ref seq difference location. Check whether annotation is correct: "
                    f"found alleles=({annotated_alleles[0]}, {annotated_alleles[1]}), "
                    f"annotation={variant_annotation}"
                )
        elif len(rs_ids) > 1 and any(panel.has_ref_seq_difference_annotation(gene, rs_id) for rs_id in rs_ids):
            error_msg = f"One of multiple rs ids is of a ref seq difference, so not sure how to annotate {rs_ids}"
            raise ValueError(error_msg)
        else:
            # no ref seq differences involved
            variant_annotation = grch37_call.variant_annotation
            filter_type = grch37_call.filter

        full_call = FullCall(
            start_coordinate_grch37, reference_allele_grch37, start_coordinate_grch38, reference_allele_grch38,
            annotated_alleles, gene, rs_ids, variant_annotation, filter_type,
        )
        return full_call

    @classmethod
    def __get_inferred_ref_calls(cls, panel: Panel, handled_grch37_coordinates: Set[GeneCoordinate],
                                 handled_grch37_rs_ids: Set[str]) -> List[FullCall]:
        # The "inferred ref" calls are calls for ref seq differences that do not correspond to grch37 calls.
        # This means that they were ref vs GRCh37 and therefore variant calls vs GRCh38
        inferred_ref_calls = []
        for rs_id_info, gene, annotation in panel.get_ref_seq_differences():
            if not rs_id_info.get_relevant_grch37_coordinates().intersection(handled_grch37_coordinates):
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
                    rs_id_info.reference_allele_grch37,
                    rs_id_info.start_coordinate_grch38,
                    rs_id_info.reference_allele_grch38,
                    annotated_alleles,
                    gene,
                    (rs_id_info.rs_id,),
                    annotation,
                    "INFERRED_REF_CALL",
                )
                inferred_ref_calls.append(full_call)
        return inferred_ref_calls

    @classmethod
    def __get_annotated_alleles_from_rs_id_info(
            cls, rs_id_info: RsIdInfo, grch37_call: Grch37Call) -> Tuple[AnnotatedAllele, AnnotatedAllele]:
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

    @classmethod
    def __assert_rs_id_call_matches_info(cls, rs_ids_call: Tuple[str, ...], rs_ids_info: Tuple[str, ...]) -> None:
        if rs_ids_call != (".",) and rs_ids_call != rs_ids_info:
            # TODO: make this more flexible, if necessary
            error_msg = (f"Given rs id does not match rs id from panel: "
                         f"from call={rs_ids_call}, from panel={rs_ids_info}")
            raise ValueError(error_msg)

    @classmethod
    def __assert_gene_in_panel(cls, gene: str, panel: Panel) -> None:
        if gene not in panel.get_genes():
            error_msg = f"Call for unknown gene:\ngene={gene}"
            raise ValueError(error_msg)

    @classmethod
    def __assert_alleles_in_expected_order(cls, annotated_alleles: Tuple[AnnotatedAllele, AnnotatedAllele]) -> None:
        alleles_in_unexpected_order = (
            annotated_alleles[0].is_annotated()
            and annotated_alleles[1].is_annotated()
            and annotated_alleles[0].is_variant_vs_grch37
            and not annotated_alleles[1].is_variant_vs_grch37
        )
        if alleles_in_unexpected_order:
            error_msg = (f"Alleles are in unexpected order, alt before ref: "
                         f"alleles=({annotated_alleles[0]},{annotated_alleles[1]})")
            raise ValueError(error_msg)
