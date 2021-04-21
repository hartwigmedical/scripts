from typing import Set, Tuple, Optional, FrozenSet

from base.constants import REF_CALL_ANNOTATION_STRING
from base.filter import Filter
from base.gene_coordinate import GeneCoordinate
from call_data import V37CallData, FullCall, AnnotatedAllele, V37Call
from config.panel import Panel


class V37CallTranslator(object):
    @classmethod
    def get_full_calls(cls, v37_call_data: V37CallData, panel: Panel) -> FrozenSet[FullCall]:
        handled_v37_coordinates: Set[GeneCoordinate] = set()
        handled_v37_rs_ids: Set[str] = set()
        full_calls_from_v37_calls = set()
        for v37_call in v37_call_data.calls:
            full_call = cls.__get_full_call_from_v37_call(v37_call, panel)

            for rs_id in full_call.rs_ids:
                if rs_id != ".":
                    if rs_id in handled_v37_rs_ids:
                        error_msg = (f"Call for rs id that has already been handled:\n"
                                     f"call={v37_call}\n"
                                     f"handled_rs_ids={handled_v37_rs_ids}")
                        raise ValueError(error_msg)
                    handled_v37_rs_ids.add(rs_id)

            relevant_v37_coordinates = full_call.get_relevant_v37_coordinates()
            if relevant_v37_coordinates.intersection(handled_v37_coordinates):
                warning_msg = (
                    f"[WARN] Call involves at least one position that has already been handled:\n"
                    f"call={v37_call}\n"
                    f"handled_coords={handled_v37_coordinates}"
                )
                print(warning_msg)
            handled_v37_coordinates.update(relevant_v37_coordinates)

            full_calls_from_v37_calls.add(full_call)

        return frozenset(full_calls_from_v37_calls)

    @classmethod
    def __get_full_call_from_v37_call(cls, v37_call: V37Call, panel: Panel) -> FullCall:
        cls.__assert_gene_in_panel(v37_call.gene, panel)

        # determine start_coordinate_v38, reference_allele_v38 and rs_ids
        start_coordinate_v38: Optional[GeneCoordinate]
        reference_allele_v38: Optional[str]
        rs_ids: Tuple[str, ...]
        if panel.contains_rs_id_matching_v37_call(v37_call):
            rs_id_info = panel.get_matching_rs_id_info(v37_call.start_coordinate, v37_call.ref_allele)
            cls.__assert_rs_id_call_matches_info(v37_call.rs_ids, (rs_id_info.rs_id,))

            start_coordinate_v38 = rs_id_info.start_coordinate_v38
            reference_allele_v38 = rs_id_info.reference_allele_v38
            if v37_call.rs_ids == (".",) and rs_id_info is not None:
                rs_ids = (rs_id_info.rs_id,)
            else:
                rs_ids = v37_call.rs_ids
        else:
            # unknown variant
            start_coordinate_v38 = None
            reference_allele_v38 = None
            rs_ids = v37_call.rs_ids

        annotated_alleles = (
            AnnotatedAllele.from_alleles(v37_call.alleles[0], v37_call.ref_allele, reference_allele_v38),
            AnnotatedAllele.from_alleles(v37_call.alleles[1], v37_call.ref_allele, reference_allele_v38),
        )

        # determine variant annotation v38 and filter v38
        if panel.has_ref_seq_difference_annotation(v37_call.gene, v37_call.start_coordinate, v37_call.ref_allele):
            v38_ref_call_due_to_ref_sequence_difference = all(
                annotated.is_variant_vs_v37
                and annotated.is_annotated_vs_v38()
                and not annotated.is_variant_vs_v38
                for annotated in annotated_alleles
            )
            all_variants_ref_to_v37_or_v38 = all(
                not annotated.is_variant_vs_v37
                or (annotated.is_annotated_vs_v38() and not annotated.is_variant_vs_v38)
                for annotated in annotated_alleles
            )

            if v38_ref_call_due_to_ref_sequence_difference:
                variant_annotation_v38 = REF_CALL_ANNOTATION_STRING
                filter_type_v38 = Filter.PASS
            elif all_variants_ref_to_v37_or_v38:
                variant_annotation_v38 = panel.get_ref_seq_difference_annotation(
                    v37_call.gene, v37_call.start_coordinate, v37_call.ref_allele)
                filter_type_v38 = Filter.PASS
            else:
                variant_annotation_v38 = v37_call.variant_annotation + "?"
                filter_type_v38 = Filter.UNKNOWN
                print(
                    f"[WARN] Unexpected allele in ref seq difference location. Check whether annotation is correct: "
                    f"found alleles=({annotated_alleles[0]}, {annotated_alleles[1]}), "
                    f"annotation={variant_annotation_v38}"
                )
        elif panel.contains_rs_id_matching_v37_call(v37_call):
            # known variant and no ref seq differences involved
            variant_annotation_v38 = v37_call.variant_annotation
            filter_type_v38 = Filter.PASS
        else:
            # unknown variant, no ref seq difference involved
            variant_annotation_v38 = v37_call.variant_annotation + "?"
            filter_type_v38 = Filter.UNKNOWN
            print(
                f"[WARN] Unknown variant. Check whether annotation is correct: "
                f"found alleles=({annotated_alleles[0]}, {annotated_alleles[1]}), "
                f"annotation={variant_annotation_v38}"
            )

        full_call = FullCall(
            v37_call.start_coordinate, v37_call.ref_allele, start_coordinate_v38, reference_allele_v38,
            v37_call.alleles, v37_call.gene, rs_ids,
            v37_call.variant_annotation, v37_call.filter, variant_annotation_v38, filter_type_v38,
        )
        return full_call

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
