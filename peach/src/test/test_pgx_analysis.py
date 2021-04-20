import unittest
from typing import Dict, Set

from base.filter import Filter
from base.gene_coordinate import GeneCoordinate
from call_data import V37CallData, V37Call, HaplotypeCall, FullCall
from config.drug_info import DrugInfo
from config.gene_info import GeneInfo
from config.haplotype import Haplotype
from config.panel import Panel
from config.rs_id_info import RsIdInfo
from config.variant import Variant
from pgx_analysis import PgxAnalyser, PgxAnalysis


class TestPgxAnalysis(unittest.TestCase):
    @classmethod
    def __get_wide_example_panel(cls) -> Panel:
        dpyd_two_a_variant = Variant("rs3918290", "T")
        dpyd_two_b_variant = Variant("rs1801159", "C")
        dpyd_three_variant = Variant("rs72549303", "TG")
        fake_variant = Variant("rs1212125", "C")
        fake2_variant = Variant("rs1212127", "C")

        dpyd_haplotypes = frozenset({
            Haplotype("*2A", "No Function", frozenset({dpyd_two_a_variant})),
            Haplotype("*2B", "No Function", frozenset({dpyd_two_a_variant, dpyd_two_b_variant})),
            Haplotype("*3", "Normal Function", frozenset({dpyd_three_variant})),
        })
        dpyd_rs_id_infos = frozenset({
            RsIdInfo("rs3918290", "C", "C", GeneCoordinate("1", 97915614), GeneCoordinate("1", 97450058)),
            RsIdInfo("rs72549309", "GATGA", "GATGA", GeneCoordinate("1", 98205966), GeneCoordinate("1", 97740410)),
            RsIdInfo("rs1801159", "T", "T", GeneCoordinate("1", 97981395), GeneCoordinate("1", 97515839)),
            RsIdInfo("rs72549303", "TG", "TC", GeneCoordinate("1", 97915621), GeneCoordinate("1", 97450065)),
        })
        dpyd_drugs = frozenset({
            DrugInfo("5-Fluorouracil", "https://www.pharmgkb.org/chemical/PA128406956/guidelineAnnotation/PA166104939"),
            DrugInfo("Capecitabine", "https://www.pharmgkb.org/chemical/PA448771/guidelineAnnotation/PA166104963"),
        })
        dpyd_rs_id_to_difference_annotations = {
            "rs72549303": "6744GA>CA",
        }

        fake_haplotypes = frozenset({
            Haplotype("*4A", "Reduced Function", frozenset({fake_variant})),
        })
        fake_rs_id_infos = frozenset({
            RsIdInfo("rs1212125", "T", "T", GeneCoordinate("5", 97915617), GeneCoordinate("5", 97450060)),
        })
        fake_drugs = frozenset({
            DrugInfo("Aspirin", "https://www.pharmgkb.org/some_other_url"),
        })
        fake_rs_id_to_difference_annotations: Dict[str, str] = {}

        fake2_haplotypes = frozenset({
            Haplotype("*4A", "Reduced Function", frozenset({fake2_variant})),
        })
        fake2_rs_id_infos = frozenset({
            RsIdInfo("rs1212127", "C", "T", GeneCoordinate("16", 97915617), GeneCoordinate("16", 97450060)),
        })
        fake2_drugs = frozenset({
            DrugInfo("Aspirin", "https://www.pharmgkb.org/some_other_url"),
        })
        fake2_rs_id_to_difference_annotations: Dict[str, str] = {"rs1212127": "1324T>C"}

        gene_infos = frozenset({
            GeneInfo("DPYD", "1", "*1", dpyd_haplotypes, dpyd_rs_id_infos,
                     dpyd_drugs, dpyd_rs_id_to_difference_annotations),
            GeneInfo("FAKE", "5", "*1", fake_haplotypes, fake_rs_id_infos,
                     fake_drugs, fake_rs_id_to_difference_annotations),
            GeneInfo("FAKE2", "16", "*1", fake2_haplotypes, fake2_rs_id_infos,
                     fake2_drugs, fake2_rs_id_to_difference_annotations),
        })
        name = "WideTestPanel"
        version = "1.0"
        return Panel(name, version, gene_infos)

    @classmethod
    def __get_narrow_example_panel(cls, included_haplotypes: Set[str]) -> Panel:
        dpyd_two_a_variant = Variant("rs3918290", "T")
        dpyd_five_variant = Variant("rs1801159", "C")
        dpyd_three_variant = Variant("rs72549303", "TG")
        dpyd_seven_variant = Variant("rs2938101", "AGT")

        possible_dpyd_haplotypes = [
            Haplotype("*2A", "No Function", frozenset({dpyd_two_a_variant})),
            Haplotype("*2B", "No Function", frozenset({dpyd_two_a_variant, dpyd_five_variant})),
            Haplotype("*3", "Normal Function", frozenset({dpyd_three_variant})),
            Haplotype("*5", "Normal Function", frozenset({dpyd_five_variant})),
            Haplotype("*7", "Normal Function", frozenset({dpyd_seven_variant})),
            Haplotype("*9", "Normal Function", frozenset({dpyd_two_a_variant, dpyd_three_variant})),
            Haplotype("*10", "Normal Function", frozenset({dpyd_three_variant, dpyd_five_variant})),
        ]
        possible_rs_id_infos = [
            RsIdInfo("rs3918290", "C", "C", GeneCoordinate("1", 97915614), GeneCoordinate("1", 97450058)),
            RsIdInfo("rs72549309", "GATGA", "GATGA", GeneCoordinate("1", 98205966), GeneCoordinate("1", 97740410)),
            RsIdInfo("rs1801159", "T", "T", GeneCoordinate("1", 97981395), GeneCoordinate("1", 97515839)),
            RsIdInfo("rs72549303", "TG", "TC", GeneCoordinate("1", 97915621), GeneCoordinate("1", 97450065)),
            RsIdInfo("rs2938101", "A", "A", GeneCoordinate("1", 97912838), GeneCoordinate("1", 97453984)),
        ]

        unknown_haplotypes = included_haplotypes.difference({haplotype.name for haplotype in possible_dpyd_haplotypes})
        if unknown_haplotypes:
            raise ValueError(f"Not all dpyd haplotype names recognized: unknown={unknown_haplotypes}")

        # use a set to cause the order of haplotypes to be random, to make sure this order will not matter.
        included_dpyd_haplotypes = frozenset(
            {haplotype for haplotype in possible_dpyd_haplotypes if haplotype.name in included_haplotypes}
        )
        included_dpyd_rs_ids = {"rs72549303"}.union(
            {variant.rs_id for haplotype in included_dpyd_haplotypes for variant in haplotype.variants}
        )
        included_dpyd_rs_id_infos = frozenset(
            {rs_id_info for rs_id_info in possible_rs_id_infos if rs_id_info.rs_id in included_dpyd_rs_ids}
        )
        dpyd_drugs = frozenset({
            DrugInfo("5-Fluorouracil", "https://www.pharmgkb.org/chemical/PA128406956/guidelineAnnotation/PA166104939"),
            DrugInfo("Capecitabine", "https://www.pharmgkb.org/chemical/PA448771/guidelineAnnotation/PA166104963"),
        })
        dpyd_rs_id_to_difference_annotations = {
            "rs72549303": "6744GA>CA",
        }

        gene_infos = frozenset({
            GeneInfo("DPYD", "1", "*1", included_dpyd_haplotypes, included_dpyd_rs_id_infos,
                     dpyd_drugs, dpyd_rs_id_to_difference_annotations),
        })
        name = "NarrowTestPanel"
        version = "1.0"
        return Panel(name, version, gene_infos)

    def test_empty(self) -> None:
        """No variants wrt v37"""
        panel = self.__get_wide_example_panel()
        ids_found_in_patient = V37CallData(frozenset())
        pgx_analysis = PgxAnalyser.create_pgx_analysis(ids_found_in_patient, panel)

        gene_to_haplotype_calls_expected = {
            "DPYD": {HaplotypeCall("*3", 2)}, "FAKE": {HaplotypeCall("*1", 2)}, "FAKE2": {HaplotypeCall("*4A", 2)}
        }

        all_full_calls_expected = frozenset({
            FullCall(
                GeneCoordinate("1", 97915614), "C", GeneCoordinate("1", 97450058), "C",
                ("C", "C"), "DPYD", ("rs3918290",), "REF_CALL", Filter.NO_CALL, "REF_CALL", Filter.NO_CALL,
            ),
            FullCall(
                GeneCoordinate("1", 97915621), "TG", GeneCoordinate("1", 97450065), "TC",
                ("TG", "TG"), "DPYD", ("rs72549303",), "REF_CALL", Filter.NO_CALL, "6744GA>CA", Filter.INFERRED_PASS,
            ),
            FullCall(
                GeneCoordinate("1", 97981395), "T", GeneCoordinate("1", 97515839), "T",
                ("T", "T"), "DPYD", ("rs1801159",), "REF_CALL", Filter.NO_CALL, "REF_CALL", Filter.NO_CALL,
            ),
            FullCall(
                GeneCoordinate("1", 98205966), "GATGA", GeneCoordinate("1", 97740410), "GATGA",
                ("GATGA", "GATGA"), "DPYD", ("rs72549309",), "REF_CALL", Filter.NO_CALL, "REF_CALL", Filter.NO_CALL,
            ),
            FullCall(
                GeneCoordinate("16", 97915617), "C", GeneCoordinate("16", 97450060), "T",
                ("C", "C"), "FAKE2", ("rs1212127",), "REF_CALL", Filter.NO_CALL, "1324T>C", Filter.INFERRED_PASS,
            ),
            FullCall(
                GeneCoordinate("5", 97915617), "T", GeneCoordinate("5", 97450060), "T",
                ("T", "T"), "FAKE", ("rs1212125",), "REF_CALL", Filter.NO_CALL, "REF_CALL", Filter.NO_CALL,
            ),
        })
        pgx_analysis_expected = PgxAnalysis(all_full_calls_expected, gene_to_haplotype_calls_expected)
        self.assertEqual(pgx_analysis_expected, pgx_analysis)

    def test_hom_ref(self) -> None:
        """All haplotypes are  *1_HOM"""
        panel = self.__get_wide_example_panel()
        ids_found_in_patient = V37CallData(frozenset({
            V37Call(
                GeneCoordinate("16", 97915617), "C", ("T", "T"), "FAKE2", ("rs1212127",), "1324C>T", Filter.PASS),
            V37Call(
                GeneCoordinate("1", 97915621), "TG", ("TC", "TC"), "DPYD", ("rs72549303",), "6744CA>GA", Filter.PASS),
            V37Call(
                GeneCoordinate("1", 97915614), "C", ("C", "C"), "DPYD", ("rs3918290",), "REF_CALL", Filter.PASS),
        }))
        pgx_analysis = PgxAnalyser.create_pgx_analysis(ids_found_in_patient, panel)

        gene_to_haplotype_calls_expected = {
            "DPYD": {HaplotypeCall("*1", 2)}, "FAKE": {HaplotypeCall("*1", 2)}, "FAKE2": {HaplotypeCall("*1", 2)}
        }
        all_full_calls_expected = frozenset({
            FullCall(
                GeneCoordinate("1", 97915614), "C", GeneCoordinate("1", 97450058), "C",
                ("C", "C"), "DPYD", ("rs3918290",), "REF_CALL", Filter.PASS, "REF_CALL", Filter.PASS,
            ),
            FullCall(
                GeneCoordinate("1", 97915621), "TG", GeneCoordinate("1", 97450065), "TC",
                ("TC", "TC"), "DPYD", ("rs72549303",), "6744CA>GA", Filter.PASS, "REF_CALL", Filter.PASS,
            ),
            FullCall(
                GeneCoordinate("1", 97981395), "T", GeneCoordinate("1", 97515839), "T",
                ("T", "T"), "DPYD", ("rs1801159",), "REF_CALL", Filter.NO_CALL, "REF_CALL", Filter.NO_CALL,
            ),
            FullCall(
                GeneCoordinate("1", 98205966), "GATGA", GeneCoordinate("1", 97740410), "GATGA",
                ("GATGA", "GATGA"), "DPYD", ("rs72549309",), "REF_CALL", Filter.NO_CALL, "REF_CALL", Filter.NO_CALL,
            ),
            FullCall(
                GeneCoordinate("16", 97915617), "C", GeneCoordinate("16", 97450060), "T",
                ("T", "T"), "FAKE2", ("rs1212127",), "1324C>T", Filter.PASS, "REF_CALL", Filter.PASS,
            ),
            FullCall(
                GeneCoordinate("5", 97915617), "T", GeneCoordinate("5", 97450060), "T",
                ("T", "T"), "FAKE", ("rs1212125",), "REF_CALL", Filter.NO_CALL, "REF_CALL", Filter.NO_CALL,
            ),
        })
        pgx_analysis_expected = PgxAnalysis(all_full_calls_expected, gene_to_haplotype_calls_expected)
        self.assertEqual(pgx_analysis_expected, pgx_analysis)

    def test_heterozygous(self) -> None:
        """All haplotypes are heterozygous. Both variant/ref and variant/variant"""
        panel = self.__get_wide_example_panel()
        ids_found_in_patient = V37CallData(frozenset({
            V37Call(
                GeneCoordinate("16", 97915617), "C", ("C", "T"), "FAKE2", ("rs1212127",), "1324C>T", Filter.PASS),
            V37Call(
                GeneCoordinate("5", 97915617), "T", ("T", "C"), "FAKE", ("rs1212125",), "1005T>C", Filter.PASS),
            V37Call(
                GeneCoordinate("1", 97915621), "TG", ("TG", "TC"), "DPYD", ("rs72549303",), "6744CA>GA", Filter.PASS),
            V37Call(
                GeneCoordinate("1", 97915614), "C", ("C", "T"), "DPYD", ("rs3918290",), "35G>A", Filter.PASS),
            V37Call(
                GeneCoordinate("1", 97981395), "T", ("T", "C"), "DPYD", ("rs1801159",), "674A>G", Filter.PASS),
        }))
        pgx_analysis = PgxAnalyser.create_pgx_analysis(ids_found_in_patient, panel)

        gene_to_haplotype_calls_expected = {
            "DPYD": {HaplotypeCall("*2B", 1), HaplotypeCall("*3", 1)},
            "FAKE": {HaplotypeCall("*4A", 1), HaplotypeCall("*1", 1)},
            "FAKE2": {HaplotypeCall("*4A", 1), HaplotypeCall("*1", 1)},
        }
        all_full_calls_expected = frozenset({
            FullCall(
                GeneCoordinate("1", 97915614), "C", GeneCoordinate("1", 97450058), "C",
                ("C", "T"), "DPYD", ("rs3918290",), "35G>A", Filter.PASS, "35G>A", Filter.PASS,
            ),
            FullCall(
                GeneCoordinate("1", 97915621), "TG", GeneCoordinate("1", 97450065), "TC",
                ("TG", "TC"), "DPYD", ("rs72549303",), "6744CA>GA", Filter.PASS, "6744GA>CA", Filter.PASS,
            ),
            FullCall(
                GeneCoordinate("1", 97981395), "T", GeneCoordinate("1", 97515839), "T",
                ("T", "C"), "DPYD", ("rs1801159",), "674A>G", Filter.PASS, "674A>G", Filter.PASS,
            ),
            FullCall(
                GeneCoordinate("1", 98205966), "GATGA", GeneCoordinate("1", 97740410), "GATGA",
                ("GATGA", "GATGA"), "DPYD", ("rs72549309",), "REF_CALL", Filter.NO_CALL, "REF_CALL", Filter.NO_CALL,
            ),
            FullCall(
                GeneCoordinate("16", 97915617), "C", GeneCoordinate("16", 97450060), "T",
                ("C", "T"), "FAKE2", ("rs1212127",), "1324C>T", Filter.PASS, "1324T>C", Filter.PASS,
            ),
            FullCall(
                GeneCoordinate("5", 97915617), "T", GeneCoordinate("5", 97450060), "T",
                ("T", "C"), "FAKE", ("rs1212125",), "1005T>C", Filter.PASS, "1005T>C", Filter.PASS,
            ),
        })
        pgx_analysis_expected = PgxAnalysis(all_full_calls_expected, gene_to_haplotype_calls_expected)
        self.assertEqual(pgx_analysis_expected, pgx_analysis)

    def test_ref_call_on_ref_seq_differences(self) -> None:
        """Explicit ref calls wrt v37 at differences between v37 and v38"""
        panel = self.__get_wide_example_panel()
        ids_found_in_patient = V37CallData(frozenset({
            V37Call(
                GeneCoordinate("16", 97915617), "C", ("C", "C"), "FAKE2", ("rs1212127",), "REF_CALL", Filter.PASS),
            V37Call(
                GeneCoordinate("1", 97915621), "TG", ("TG", "TG"), "DPYD", ("rs72549303",), "REF_CALL", Filter.PASS),
        }))
        pgx_analysis = PgxAnalyser.create_pgx_analysis(ids_found_in_patient, panel)

        gene_to_haplotype_calls_expected = {
            "DPYD": {HaplotypeCall("*3", 2)}, "FAKE": {HaplotypeCall("*1", 2)}, "FAKE2": {HaplotypeCall("*4A", 2)}
        }
        all_full_calls_expected = frozenset({
            FullCall(
                GeneCoordinate("1", 97915614), "C", GeneCoordinate("1", 97450058), "C",
                ("C", "C"), "DPYD", ("rs3918290",), "REF_CALL", Filter.NO_CALL, "REF_CALL", Filter.NO_CALL,
            ),
            FullCall(
                GeneCoordinate("1", 97915621), "TG", GeneCoordinate("1", 97450065), "TC",
                ("TG", "TG"), "DPYD", ("rs72549303",), "REF_CALL", Filter.PASS, "6744GA>CA", Filter.PASS,
            ),
            FullCall(
                GeneCoordinate("1", 97981395), "T", GeneCoordinate("1", 97515839), "T",
                ("T", "T"), "DPYD", ("rs1801159",), "REF_CALL", Filter.NO_CALL, "REF_CALL", Filter.NO_CALL,
            ),
            FullCall(
                GeneCoordinate("1", 98205966), "GATGA", GeneCoordinate("1", 97740410), "GATGA",
                ("GATGA", "GATGA"), "DPYD", ("rs72549309",), "REF_CALL", Filter.NO_CALL, "REF_CALL", Filter.NO_CALL,
            ),
            FullCall(
                GeneCoordinate("16", 97915617), "C", GeneCoordinate("16", 97450060), "T",
                ("C", "C"), "FAKE2", ("rs1212127",), "REF_CALL", Filter.PASS, "1324T>C", Filter.PASS,
            ),
            FullCall(
                GeneCoordinate("5", 97915617), "T", GeneCoordinate("5", 97450060), "T",
                ("T", "T"), "FAKE", ("rs1212125",), "REF_CALL", Filter.NO_CALL, "REF_CALL", Filter.NO_CALL,
            ),
        })
        pgx_analysis_expected = PgxAnalysis(all_full_calls_expected, gene_to_haplotype_calls_expected)
        self.assertEqual(pgx_analysis_expected, pgx_analysis)

    def test_only_position_match_on_ref_seq_differences(self) -> None:
        """At reference sequence differences: heterozygous between ref v37 and v38, and no rs_id provided"""
        panel = self.__get_wide_example_panel()
        ids_found_in_patient = V37CallData(frozenset({
            V37Call(GeneCoordinate("16", 97915617), "C", ("C", "T"), "FAKE2", (".",), "1324C>T", Filter.PASS),
            V37Call(GeneCoordinate("1", 97915621), "TG", ("TG", "TC"), "DPYD", (".",), "6744CA>GA", Filter.PASS),
        }))
        pgx_analysis = PgxAnalyser.create_pgx_analysis(ids_found_in_patient, panel)

        gene_to_haplotype_calls_expected = {
            "DPYD": {HaplotypeCall("*3", 1), HaplotypeCall("*1", 1)},
            "FAKE": {HaplotypeCall("*1", 2)},
            "FAKE2": {HaplotypeCall("*4A", 1), HaplotypeCall("*1", 1)},
        }
        all_full_calls_expected = frozenset({
            FullCall(
                GeneCoordinate("1", 97915614), "C", GeneCoordinate("1", 97450058), "C",
                ("C", "C"), "DPYD", ("rs3918290",), "REF_CALL", Filter.NO_CALL, "REF_CALL", Filter.NO_CALL,
            ),
            FullCall(
                GeneCoordinate("1", 97915621), "TG", GeneCoordinate("1", 97450065), "TC",
                ("TG", "TC"), "DPYD", ("rs72549303",), "6744CA>GA", Filter.PASS, "6744GA>CA", Filter.PASS,
            ),
            FullCall(
                GeneCoordinate("1", 97981395), "T", GeneCoordinate("1", 97515839), "T",
                ("T", "T"), "DPYD", ("rs1801159",), "REF_CALL", Filter.NO_CALL, "REF_CALL", Filter.NO_CALL,
            ),
            FullCall(
                GeneCoordinate("1", 98205966), "GATGA", GeneCoordinate("1", 97740410), "GATGA",
                ("GATGA", "GATGA"), "DPYD", ("rs72549309",), "REF_CALL", Filter.NO_CALL, "REF_CALL", Filter.NO_CALL,
            ),
            FullCall(
                GeneCoordinate("16", 97915617), "C", GeneCoordinate("16", 97450060), "T",
                ("C", "T"), "FAKE2", ("rs1212127",), "1324C>T", Filter.PASS, "1324T>C", Filter.PASS,
            ),
            FullCall(
                GeneCoordinate("5", 97915617), "T", GeneCoordinate("5", 97450060), "T",
                ("T", "T"), "FAKE", ("rs1212125",), "REF_CALL", Filter.NO_CALL, "REF_CALL", Filter.NO_CALL,
            ),
        })
        pgx_analysis_expected = PgxAnalysis(all_full_calls_expected, gene_to_haplotype_calls_expected)
        self.assertEqual(pgx_analysis_expected, pgx_analysis)

    def test_wrong_rs_id_on_ref_seq_differences(self) -> None:
        """
        At reference sequence differences: heterozygous between ref v37 and v38,
        and incorrect rs id provided
        """
        panel = self.__get_wide_example_panel()
        ids_found_in_patient = V37CallData(frozenset({
            V37Call(
                GeneCoordinate("16", 97915617), "C", ("C", "T"), "FAKE2", ("rs939535",), "1324C>T", Filter.PASS),
            V37Call(
                GeneCoordinate("1", 97915621), "TG", ("TG", "TC"), "DPYD", ("rs4020942",), "6744CA>GA", Filter.PASS),
        }))
        with self.assertRaises(ValueError):
            PgxAnalyser.create_pgx_analysis(ids_found_in_patient, panel)

    def test_wrong_position_on_ref_seq_differences(self) -> None:
        """Explicit ref calls wrt v37 at differences between v37 and v38, except positions are incorrect"""
        panel = self.__get_wide_example_panel()
        ids_found_in_patient = V37CallData(frozenset({
            V37Call(
                GeneCoordinate("16", 97915618), "C", ("C", "C"), "FAKE2", ("rs1212127",), "REF_CALL", Filter.PASS),
            V37Call(
                GeneCoordinate("1", 97915623), "TG", ("TG", "TG"), "DPYD", ("rs72549303",), "REF_CALL", Filter.PASS),
        }))
        with self.assertRaises(ValueError):
            PgxAnalyser.create_pgx_analysis(ids_found_in_patient, panel)

    def test_wrong_chromosome_on_ref_seq_differences(self) -> None:
        """Explicit ref calls wrt v37 at differences between v37 and v38, except chromosomes are incorrect"""
        panel = self.__get_wide_example_panel()
        ids_found_in_patient = V37CallData(frozenset({
            V37Call(
                GeneCoordinate("7", 97915617), "C", ("C", "C"), "FAKE2", ("rs1212127",), "REF_CALL", Filter.PASS),
            V37Call(
                GeneCoordinate("8", 97915621), "TG", ("TG", "TG"), "DPYD", ("rs72549303",), "REF_CALL", Filter.PASS),
        }))
        with self.assertRaises(ValueError):
            PgxAnalyser.create_pgx_analysis(ids_found_in_patient, panel)

    def test_position_of_other_variant_on_ref_seq_differences(self) -> None:
        """
        Explicit ref calls wrt v37 at differences between v37 and v38, except position is of other variant
        """
        panel = self.__get_wide_example_panel()
        ids_found_in_patient = V37CallData(frozenset({
            V37Call(
                GeneCoordinate("16", 97915617), "C", ("C", "C"), "FAKE2", ("rs1212127",), "REF_CALL", Filter.PASS),
            V37Call(
                GeneCoordinate("1", 98205966), "TG", ("TG", "TG"), "DPYD", ("rs72549303",), "REF_CALL", Filter.PASS),
        }))
        with self.assertRaises(ValueError):
            PgxAnalyser.create_pgx_analysis(ids_found_in_patient, panel)

    def test_single_different_allele_on_ref_seq_differences(self) -> None:
        """At reference sequence differences: single allele that is ref v37 or v38, other allele is neither"""
        panel = self.__get_wide_example_panel()
        ids_found_in_patient = V37CallData(frozenset({
            V37Call(
                GeneCoordinate("16", 97915617), "C", ("C", "A"),
                "FAKE2", ("rs1212127",), "1324C>A", Filter.PASS),  # with ref v37
            V37Call(
                GeneCoordinate("1", 97915621), "TG", ("AC", "TC"),
                "DPYD", ("rs72549303",), "6744CT>GT;6744CT>GC", Filter.PASS),  # with ref v38
        }))
        pgx_analysis = PgxAnalyser.create_pgx_analysis(ids_found_in_patient, panel)

        gene_to_haplotype_calls_expected: Dict[str, Set[HaplotypeCall]] = {
            "DPYD": set(), "FAKE": {HaplotypeCall("*1", 2)}, "FAKE2": set()
        }
        all_full_calls_expected = frozenset({
            FullCall(
                GeneCoordinate("1", 97915614), "C", GeneCoordinate("1", 97450058), "C",
                ("C", "C"), "DPYD", ("rs3918290",), "REF_CALL", Filter.NO_CALL, "REF_CALL", Filter.NO_CALL,
            ),
            FullCall(
                GeneCoordinate("1", 97915621), "TG", GeneCoordinate("1", 97450065), "TC",
                ("AC", "TC"), "DPYD", ("rs72549303",), "6744CT>GT;6744CT>GC", Filter.PASS, "6744CT>GT;6744CT>GC?", Filter.UNKNOWN,
            ),
            FullCall(
                GeneCoordinate("1", 97981395), "T", GeneCoordinate("1", 97515839), "T",
                ("T", "T"), "DPYD", ("rs1801159",), "REF_CALL", Filter.NO_CALL, "REF_CALL", Filter.NO_CALL,
            ),
            FullCall(
                GeneCoordinate("1", 98205966), "GATGA", GeneCoordinate("1", 97740410), "GATGA",
                ("GATGA", "GATGA"), "DPYD", ("rs72549309",), "REF_CALL", Filter.NO_CALL, "REF_CALL", Filter.NO_CALL,
            ),
            FullCall(
                GeneCoordinate("16", 97915617), "C", GeneCoordinate("16", 97450060), "T",
                ("C", "A"), "FAKE2", ("rs1212127",), "1324C>A", Filter.PASS, "1324C>A?", Filter.UNKNOWN,
            ),
            FullCall(
                GeneCoordinate("5", 97915617), "T", GeneCoordinate("5", 97450060), "T",
                ("T", "T"), "FAKE", ("rs1212125",), "REF_CALL", Filter.NO_CALL, "REF_CALL", Filter.NO_CALL,
            ),
        })
        pgx_analysis_expected = PgxAnalysis(all_full_calls_expected, gene_to_haplotype_calls_expected)
        self.assertEqual(pgx_analysis_expected, pgx_analysis)

    def test_double_different_allele_on_ref_seq_differences(self) -> None:
        """At reference sequence differences: both alleles not ref v37 or v38"""
        panel = self.__get_wide_example_panel()
        ids_found_in_patient = V37CallData(frozenset({
            V37Call(
                GeneCoordinate("16", 97915617), "C", ("A", "G"),
                "FAKE2", ("rs1212127",), "1324C>A;1324C>G", Filter.PASS),
            V37Call(
                GeneCoordinate("1", 97915621), "TG", ("AC", "AG"),
                "DPYD", ("rs72549303",), "6744CT>GT;6744CT>GC", Filter.PASS),
        }))
        pgx_analysis = PgxAnalyser.create_pgx_analysis(ids_found_in_patient, panel)

        gene_to_haplotype_calls_expected: Dict[str, Set[HaplotypeCall]] = {
            "DPYD": set(), "FAKE": {HaplotypeCall("*1", 2)}, "FAKE2": set()
        }
        all_full_calls_expected = frozenset({
            FullCall(
                GeneCoordinate("1", 97915614), "C", GeneCoordinate("1", 97450058), "C",
                ("C", "C"), "DPYD", ("rs3918290",), "REF_CALL", Filter.NO_CALL, "REF_CALL", Filter.NO_CALL,
            ),
            FullCall(
                GeneCoordinate("1", 97915621), "TG", GeneCoordinate("1", 97450065), "TC",
                ("AC", "AG"), "DPYD", ("rs72549303",), "6744CT>GT;6744CT>GC", Filter.PASS, "6744CT>GT;6744CT>GC?", Filter.UNKNOWN,
            ),
            FullCall(
                GeneCoordinate("1", 97981395), "T", GeneCoordinate("1", 97515839), "T",
                ("T", "T"), "DPYD", ("rs1801159",), "REF_CALL", Filter.NO_CALL, "REF_CALL", Filter.NO_CALL,
            ),
            FullCall(
                GeneCoordinate("1", 98205966), "GATGA", GeneCoordinate("1", 97740410), "GATGA",
                ("GATGA", "GATGA"), "DPYD", ("rs72549309",), "REF_CALL", Filter.NO_CALL, "REF_CALL", Filter.NO_CALL,
            ),
            FullCall(
                GeneCoordinate("16", 97915617), "C", GeneCoordinate("16", 97450060), "T",
                ("A", "G"), "FAKE2", ("rs1212127",), "1324C>A;1324C>G", Filter.PASS, "1324C>A;1324C>G?", Filter.UNKNOWN,
            ),
            FullCall(
                GeneCoordinate("5", 97915617), "T", GeneCoordinate("5", 97450060), "T",
                ("T", "T"), "FAKE", ("rs1212125",), "REF_CALL", Filter.NO_CALL, "REF_CALL", Filter.NO_CALL,
            ),
        })
        pgx_analysis_expected = PgxAnalysis(all_full_calls_expected, gene_to_haplotype_calls_expected)
        self.assertEqual(pgx_analysis_expected, pgx_analysis)

    def test_unknown_variants(self) -> None:
        """Variants that are completely unknown, including with unknown rs id"""
        panel = self.__get_wide_example_panel()
        ids_found_in_patient = V37CallData(frozenset({
            V37Call(
                GeneCoordinate("16", 39593405), "A", ("A", "G"),
                "FAKE2", ("rs1949223",), "384C>T", Filter.PASS),  # unknown
            V37Call(
                GeneCoordinate("5", 97915617), "T", ("T", "C"),
                "FAKE", ("rs1212125",), "1005T>C", Filter.PASS),
            V37Call(
                GeneCoordinate("1", 2488242), "AC", ("AC", "AG"),
                "DPYD", (".",), "9213CT>GT", Filter.PASS),  # unknown
        }))
        pgx_analysis = PgxAnalyser.create_pgx_analysis(ids_found_in_patient, panel)

        gene_to_haplotype_calls_expected: Dict[str, Set[HaplotypeCall]] = {
            "DPYD": set(), "FAKE": {HaplotypeCall("*4A", 1), HaplotypeCall("*1", 1)}, "FAKE2": set()}
        all_full_calls_expected = frozenset({
            FullCall(
                GeneCoordinate("1", 2488242), "AC", None, None,
                ("AC", "AG"), "DPYD", (".",), "9213CT>GT", Filter.PASS, "9213CT>GT?", Filter.UNKNOWN,
            ),
            FullCall(
                GeneCoordinate("1", 97915614), "C", GeneCoordinate("1", 97450058), "C",
                ("C", "C"), "DPYD", ("rs3918290",), "REF_CALL", Filter.NO_CALL, "REF_CALL", Filter.NO_CALL,
            ),
            FullCall(
                GeneCoordinate("1", 97915621), "TG", GeneCoordinate("1", 97450065), "TC",
                ("TG", "TG"), "DPYD", ("rs72549303",), "REF_CALL", Filter.NO_CALL, "6744GA>CA", Filter.INFERRED_PASS,
            ),
            FullCall(
                GeneCoordinate("1", 97981395), "T", GeneCoordinate("1", 97515839), "T",
                ("T", "T"), "DPYD", ("rs1801159",), "REF_CALL", Filter.NO_CALL, "REF_CALL", Filter.NO_CALL,
            ),
            FullCall(
                GeneCoordinate("1", 98205966), "GATGA", GeneCoordinate("1", 97740410), "GATGA",
                ("GATGA", "GATGA"), "DPYD", ("rs72549309",), "REF_CALL", Filter.NO_CALL, "REF_CALL", Filter.NO_CALL,
            ),
            FullCall(
                GeneCoordinate("16", 39593405), "A", None, None,
                ("A", "G"), "FAKE2", ("rs1949223",), "384C>T", Filter.PASS, "384C>T?", Filter.UNKNOWN,
            ),
            FullCall(
                GeneCoordinate("16", 97915617), "C", GeneCoordinate("16", 97450060), "T",
                ("C", "C"), "FAKE2", ("rs1212127",), "REF_CALL", Filter.NO_CALL, "1324T>C", Filter.INFERRED_PASS,
            ),
            FullCall(
                GeneCoordinate("5", 97915617), "T", GeneCoordinate("5", 97450060), "T",
                ("T", "C"), "FAKE", ("rs1212125",), "1005T>C", Filter.PASS, "1005T>C", Filter.PASS,
            ),
        })
        pgx_analysis_expected = PgxAnalysis(all_full_calls_expected, gene_to_haplotype_calls_expected)
        self.assertEqual(pgx_analysis_expected, pgx_analysis)

    def test_unknown_gene(self) -> None:
        """Variants that are of an unknown gene"""
        panel = self.__get_wide_example_panel()
        ids_found_in_patient = V37CallData(frozenset({
            V37Call(
                GeneCoordinate("5", 97915617), "T", ("T", "C"), "FAKE", ("rs1212125",), "1005T>C", Filter.PASS),
            V37Call(
                GeneCoordinate("3", 18473423), "T", ("T", "C"), ".", ("rs2492932",), "12T>C", Filter.PASS),  # unknown
        }))
        with self.assertRaises(ValueError):
            PgxAnalyser.create_pgx_analysis(ids_found_in_patient, panel)

    def test_known_variants_with_incorrect_or_missing_rs_id(self) -> None:
        """Known variants (not ref seq differences) with incorrect or unknown rs id"""
        panel = self.__get_wide_example_panel()
        ids_found_in_patient = V37CallData(frozenset({
            V37Call(
                GeneCoordinate("16", 97915617), "C", ("C", "T"),
                "FAKE2", ("rs1212127",), "1324C>T", Filter.PASS),
            V37Call(
                GeneCoordinate("5", 97915617), "T", ("T", "C"),
                "FAKE", ("rs27384",), "1005T>C", Filter.PASS),  # incorrect
            V37Call(
                GeneCoordinate("1", 97915621), "TG", ("TG", "TC"),
                "DPYD", ("rs72549303",), "6744CA>GA", Filter.PASS),
            V37Call(
                GeneCoordinate("1", 97915614), "C", ("C", "T"),
                "DPYD", (".",), "35G>A", Filter.PASS),  # missing
            V37Call(
                GeneCoordinate("1", 97981395), "T", ("T", "C"),
                "DPYD", ("rs1801159",), "674A>G", Filter.PASS),
        }))
        with self.assertRaises(ValueError):
            PgxAnalyser.create_pgx_analysis(ids_found_in_patient, panel)

    def test_known_variant_with_incorrect_position(self) -> None:
        """Known variants (not ref seq differences), with one incorrect position"""
        panel = self.__get_wide_example_panel()
        ids_found_in_patient = V37CallData(frozenset({
            V37Call(
                GeneCoordinate("16", 97915617), "C", ("C", "T"),
                "FAKE2", ("rs1212127",), "1324C>T", Filter.PASS),
            V37Call(
                GeneCoordinate("5", 97915617), "T", ("T", "C"),
                "FAKE", ("rs1212125",), "1005T>C", Filter.PASS),
            V37Call(
                GeneCoordinate("1", 97915621), "TG", ("TG", "TC"),
                "DPYD", ("rs72549303",), "6744CA>GA", Filter.PASS),
            V37Call(
                GeneCoordinate("1", 6778543), "C", ("C", "T"),
                "DPYD", ("rs3918290",), "35G>A", Filter.PASS),  # incorrect
            V37Call(
                GeneCoordinate("1", 97981395), "T", ("T", "C"),
                "DPYD", ("rs1801159",), "674A>G", Filter.PASS),
        }))
        with self.assertRaises(ValueError):
            PgxAnalyser.create_pgx_analysis(ids_found_in_patient, panel)

    def test_known_variant_with_incorrect_chromosome(self) -> None:
        """Known variants (not ref seq differences), with one incorrect chromosome"""
        panel = self.__get_wide_example_panel()
        ids_found_in_patient = V37CallData(frozenset({
            V37Call(
                GeneCoordinate("16", 97915617), "C", ("C", "T"),
                "FAKE2", ("rs1212127",), "1324C>T", Filter.PASS),
            V37Call(
                GeneCoordinate("5", 97915617), "T", ("T", "C"),
                "FAKE", ("rs1212125",), "1005T>C", Filter.PASS),
            V37Call(
                GeneCoordinate("1", 97915621), "TG", ("TG", "TC"),
                "DPYD", ("rs72549303",), "6744CA>GA", Filter.PASS),
            V37Call(
                GeneCoordinate("3", 97915614), "C", ("C", "T"),
                "DPYD", ("rs3918290",), "35G>A", Filter.PASS),  # incorrect
            V37Call(
                GeneCoordinate("1", 97981395), "T", ("T", "C"),
                "DPYD", ("rs1801159",), "674A>G", Filter.PASS),
        }))
        with self.assertRaises(ValueError):
            PgxAnalyser.create_pgx_analysis(ids_found_in_patient, panel)

    def test_known_variant_with_multiple_rs_ids_not_matching_panel(self) -> None:
        """Multiple rs ids when panel says there should be one"""
        panel = self.__get_wide_example_panel()
        ids_found_in_patient = V37CallData(frozenset({
            V37Call(
                GeneCoordinate("16", 97915617), "C", ("C", "T"),
                "FAKE2", ("rs1212127", "rs394832"), "1324C>T", Filter.PASS),  # incorrect
            V37Call(
                GeneCoordinate("5", 97915617), "T", ("T", "C"),
                "FAKE", ("rs1212125",), "1005T>C", Filter.PASS),
            V37Call(
                GeneCoordinate("1", 97915621), "TG", ("TG", "TC"),
                "DPYD", ("rs72549303",), "6744CA>GA", Filter.PASS),
            V37Call(
                GeneCoordinate("1", 97915614), "C", ("C", "T"),
                "DPYD", ("rs3918290", "rs202093"), "35G>A", Filter.PASS),  # incorrect
            V37Call(
                GeneCoordinate("1", 97981395), "T", ("T", "C"),
                "DPYD", ("rs1801159",), "674A>G", Filter.PASS),
        }))
        with self.assertRaises(ValueError):
            PgxAnalyser.create_pgx_analysis(ids_found_in_patient, panel)

    def test_more_than_two_haplotypes_present(self) -> None:
        """More than two haplotypes are present, specifically three"""
        panel = self.__get_narrow_example_panel({"*2A", "*3", "*7"})
        ids_found_in_patient = V37CallData(frozenset({
            V37Call(
                GeneCoordinate("1", 97915614), "C", ("C", "T"), "DPYD", ("rs3918290",), "9213C>T", Filter.PASS),
            V37Call(
                GeneCoordinate("1", 97912838), "A", ("AGT", "AGT"), "DPYD", ("rs2938101",), "293A>AGT", Filter.PASS),
        }))
        pgx_analysis = PgxAnalyser.create_pgx_analysis(ids_found_in_patient, panel)

        gene_to_haplotype_calls_expected = {
            "DPYD": {HaplotypeCall("*2A", 1), HaplotypeCall("*3", 2), HaplotypeCall("*7", 2)}
        }
        all_full_calls_expected = frozenset({
            FullCall(
                GeneCoordinate("1", 97912838), "A", GeneCoordinate("1", 97453984), "A",
                ("AGT", "AGT"), "DPYD", ("rs2938101",), "293A>AGT", Filter.PASS, "293A>AGT", Filter.PASS,
            ),
            FullCall(
                GeneCoordinate("1", 97915614), "C", GeneCoordinate("1", 97450058), "C",
                ("C", "T"), "DPYD", ("rs3918290",), "9213C>T", Filter.PASS, "9213C>T", Filter.PASS,
            ),
            FullCall(
                GeneCoordinate("1", 97915621), "TG", GeneCoordinate("1", 97450065), "TC",
                ("TG", "TG"), "DPYD", ("rs72549303",), "REF_CALL", Filter.NO_CALL, "6744GA>CA", Filter.INFERRED_PASS,
            ),
        })
        pgx_analysis_expected = PgxAnalysis(all_full_calls_expected, gene_to_haplotype_calls_expected)
        self.assertEqual(pgx_analysis_expected, pgx_analysis)

    def test_ambiguous_haplotype_with_clear_winner_homozygous(self) -> None:
        """Ambiguous homozygous haplotype, where the simpler possibility should be preferred"""
        panel = self.__get_narrow_example_panel({"*2A", "*5", "*2B"})
        ids_found_in_patient = V37CallData(frozenset({
            V37Call(
                GeneCoordinate("1", 97915614), "C", ("T", "T"), "DPYD", ("rs3918290",), "9213C>T", Filter.PASS),
            V37Call(
                GeneCoordinate("1", 97981395), "T", ("C", "C"), "DPYD", ("rs1801159",), "293T>C", Filter.PASS),
            V37Call(
                GeneCoordinate("1", 97915621), "TG", ("TC", "TC"), "DPYD", ("rs72549303",), "6744CA>GA", Filter.PASS),
        }))
        pgx_analysis = PgxAnalyser.create_pgx_analysis(ids_found_in_patient, panel)

        gene_to_haplotype_calls_expected = {"DPYD": {HaplotypeCall("*2B", 2)}}
        all_full_calls_expected = frozenset({
            FullCall(
                GeneCoordinate("1", 97915614), "C", GeneCoordinate("1", 97450058), "C",
                ("T", "T"), "DPYD", ("rs3918290",), "9213C>T", Filter.PASS, "9213C>T", Filter.PASS,
            ),
            FullCall(
                GeneCoordinate("1", 97915621), "TG", GeneCoordinate("1", 97450065), "TC",
                ("TC", "TC"), "DPYD", ("rs72549303",), "6744CA>GA", Filter.PASS, "REF_CALL", Filter.PASS,
            ),
            FullCall(
                GeneCoordinate("1", 97981395), "T", GeneCoordinate("1", 97515839), "T",
                ("C", "C"), "DPYD", ("rs1801159",), "293T>C", Filter.PASS, "293T>C", Filter.PASS,
            ),
        })
        pgx_analysis_expected = PgxAnalysis(all_full_calls_expected, gene_to_haplotype_calls_expected)
        self.assertEqual(pgx_analysis_expected, pgx_analysis)

    def test_ambiguous_haplotype_with_clear_winner_heterozygous(self) -> None:
        """Ambiguous heterozygous haplotype, where the simpler possibility should be preferred"""
        panel = self.__get_narrow_example_panel({"*2A", "*5", "*2B"})
        ids_found_in_patient = V37CallData(frozenset({
            V37Call(
                GeneCoordinate("1", 97915614), "C", ("C", "T"), "DPYD", ("rs3918290",), "9213C>T", Filter.PASS),
            V37Call(
                GeneCoordinate("1", 97981395), "T", ("T", "C"), "DPYD", ("rs1801159",), "293T>C", Filter.PASS),
            V37Call(
                GeneCoordinate("1", 97915621), "TG", ("TC", "TC"), "DPYD", ("rs72549303",), "6744CA>GA", Filter.PASS),
        }))
        pgx_analysis = PgxAnalyser.create_pgx_analysis(ids_found_in_patient, panel)

        gene_to_haplotype_calls_expected = {"DPYD": {HaplotypeCall("*2B", 1), HaplotypeCall("*1", 1)}}
        all_full_calls_expected = frozenset({
            FullCall(
                GeneCoordinate("1", 97915614), "C", GeneCoordinate("1", 97450058), "C",
                ("C", "T"), "DPYD", ("rs3918290",), "9213C>T", Filter.PASS, "9213C>T", Filter.PASS,
            ),
            FullCall(
                GeneCoordinate("1", 97915621), "TG", GeneCoordinate("1", 97450065), "TC",
                ("TC", "TC"), "DPYD", ("rs72549303",), "6744CA>GA", Filter.PASS, "REF_CALL", Filter.PASS,
            ),
            FullCall(
                GeneCoordinate("1", 97981395), "T", GeneCoordinate("1", 97515839), "T",
                ("T", "C"), "DPYD", ("rs1801159",), "293T>C", Filter.PASS, "293T>C", Filter.PASS,
            ),
        })
        pgx_analysis_expected = PgxAnalysis(all_full_calls_expected, gene_to_haplotype_calls_expected)
        self.assertEqual(pgx_analysis_expected, pgx_analysis)

    def test_ambiguous_haplotype_with_clear_winner_mix(self) -> None:
        """Ambiguous mix of homozygous and heterozygous haplotype, where the simplest possibility should be preferred"""
        panel = self.__get_narrow_example_panel({"*2A", "*5", "*2B"})
        ids_found_in_patient = V37CallData(frozenset({
            V37Call(
                GeneCoordinate("1", 97915614), "C", ("T", "T"), "DPYD", ("rs3918290",), "9213C>T", Filter.PASS),
            V37Call(
                GeneCoordinate("1", 97981395), "T", ("T", "C"), "DPYD", ("rs1801159",), "293T>C", Filter.PASS),
            V37Call(
                GeneCoordinate("1", 97915621), "TG", ("TC", "TC"), "DPYD", ("rs72549303",), "6744CA>GA", Filter.PASS),
        }))
        pgx_analysis = PgxAnalyser.create_pgx_analysis(ids_found_in_patient, panel)

        gene_to_haplotype_calls_expected = {"DPYD": {HaplotypeCall("*2B", 1), HaplotypeCall("*2A", 1)}}
        all_full_calls_expected = frozenset({
            FullCall(
                GeneCoordinate("1", 97915614), "C", GeneCoordinate("1", 97450058), "C",
                ("T", "T"), "DPYD", ("rs3918290",), "9213C>T", Filter.PASS, "9213C>T", Filter.PASS,
            ),
            FullCall(
                GeneCoordinate("1", 97915621), "TG", GeneCoordinate("1", 97450065), "TC",
                ("TC", "TC"), "DPYD", ("rs72549303",), "6744CA>GA", Filter.PASS, "REF_CALL", Filter.PASS,
            ),
            FullCall(
                GeneCoordinate("1", 97981395), "T", GeneCoordinate("1", 97515839), "T",
                ("T", "C"), "DPYD", ("rs1801159",), "293T>C", Filter.PASS, "293T>C", Filter.PASS,
            ),
        })
        pgx_analysis_expected = PgxAnalysis(all_full_calls_expected, gene_to_haplotype_calls_expected)
        self.assertEqual(pgx_analysis_expected, pgx_analysis)

    def test_complicated_ambiguous_haplotype_with_a_clear_winner(self) -> None:
        """Ambiguous set of haplotypes, with no haplotype that is simplest"""
        panel = self.__get_narrow_example_panel({"*2A", "*5", "*2B", "*3", "*9", "*10", "*7"})
        ids_found_in_patient = V37CallData(frozenset({
            V37Call(
                GeneCoordinate("1", 97915614), "C", ("T", "T"), "DPYD", ("rs3918290",), "9213C>T", Filter.PASS),
            V37Call(
                GeneCoordinate("1", 97981395), "T", ("T", "C"), "DPYD", ("rs1801159",), "293T>C", Filter.PASS),
            V37Call(
                GeneCoordinate("1", 97915621), "TG", ("TG", "TC"), "DPYD", ("rs72549303",), "6744CA>GA", Filter.PASS),
            V37Call(
                GeneCoordinate("1", 97912838), "A", ("AGT", "AGT"), "DPYD", ("rs2938101",), "301A>AGT", Filter.PASS),
        }))
        pgx_analysis = PgxAnalyser.create_pgx_analysis(ids_found_in_patient, panel)

        gene_to_haplotype_calls_expected = {
            "DPYD": {HaplotypeCall("*2B", 1), HaplotypeCall("*9", 1), HaplotypeCall("*7", 2)}
        }
        all_full_calls_expected = frozenset({
            FullCall(
                GeneCoordinate("1", 97912838), "A", GeneCoordinate("1", 97453984), "A",
                ("AGT", "AGT"), "DPYD", ("rs2938101",), "301A>AGT", Filter.PASS, "301A>AGT", Filter.PASS,
            ),
            FullCall(
                GeneCoordinate("1", 97915614), "C", GeneCoordinate("1", 97450058), "C",
                ("T", "T"), "DPYD", ("rs3918290",), "9213C>T", Filter.PASS, "9213C>T", Filter.PASS,
            ),
            FullCall(
                GeneCoordinate("1", 97915621), "TG", GeneCoordinate("1", 97450065), "TC",
                ("TG", "TC"), "DPYD", ("rs72549303",), "6744CA>GA", Filter.PASS, "6744GA>CA", Filter.PASS,
            ),
            FullCall(
                GeneCoordinate("1", 97981395), "T", GeneCoordinate("1", 97515839), "T",
                ("T", "C"), "DPYD", ("rs1801159",), "293T>C", Filter.PASS, "293T>C", Filter.PASS,
            ),
        })
        pgx_analysis_expected = PgxAnalysis(all_full_calls_expected, gene_to_haplotype_calls_expected)
        self.assertEqual(pgx_analysis_expected, pgx_analysis)

    def test_ambiguous_homozygous_haplotype_with_a_less_clear_winner(self) -> None:
        """Ambiguous set of homozygous haplotypes, with winning haplotype that might not be ideal"""
        panel = self.__get_narrow_example_panel({"*2A", "*5", "*2B", "*3", "*9", "*10"})
        ids_found_in_patient = V37CallData(frozenset({
            V37Call(
                GeneCoordinate("1", 97915614), "C", ("T", "T"), "DPYD", ("rs3918290",), "9213C>T", Filter.PASS),
            V37Call(
                GeneCoordinate("1", 97981395), "T", ("C", "C"), "DPYD", ("rs1801159",), "293T>C", Filter.PASS),
            V37Call(
                GeneCoordinate("1", 97915621), "TG", ("TG", "TG"), "DPYD", ("rs72549303",), "REF_CALL", Filter.PASS),
        }))
        pgx_analysis = PgxAnalyser.create_pgx_analysis(ids_found_in_patient, panel)

        gene_to_haplotype_calls_expected = {
            "DPYD": {HaplotypeCall("*10", 1), HaplotypeCall("*2B", 1), HaplotypeCall("*9", 1)}
        }
        all_full_calls_expected = frozenset({
            FullCall(
                GeneCoordinate("1", 97915614), "C", GeneCoordinate("1", 97450058), "C",
                ("T", "T"), "DPYD", ("rs3918290",), "9213C>T", Filter.PASS, "9213C>T", Filter.PASS,
            ),
            FullCall(
                GeneCoordinate("1", 97915621), "TG", GeneCoordinate("1", 97450065), "TC",
                ("TG", "TG"), "DPYD", ("rs72549303",), "REF_CALL", Filter.PASS, "6744GA>CA", Filter.PASS,
            ),
            FullCall(
                GeneCoordinate("1", 97981395), "T", GeneCoordinate("1", 97515839), "T",
                ("C", "C"), "DPYD", ("rs1801159",), "293T>C", Filter.PASS, "293T>C", Filter.PASS,
            ),
        })
        pgx_analysis_expected = PgxAnalysis(all_full_calls_expected, gene_to_haplotype_calls_expected)
        self.assertEqual(pgx_analysis_expected, pgx_analysis)

    def test_ambiguous_heterozygous_haplotype_without_a_clear_winner(self) -> None:
        """
        Ambiguous set of heterozygous haplotypes, with no haplotype that is simplest.
        Could  be {*2A_HET,*10_HET}, {*5_HET,*9_HET} or {*3_HET,*2B_HET}
        """
        panel = self.__get_narrow_example_panel({"*2A", "*5", "*2B", "*3", "*9", "*10"})
        ids_found_in_patient = V37CallData(frozenset({
            V37Call(
                GeneCoordinate("1", 97915614), "C", ("C", "T"), "DPYD", ("rs3918290",), "9213C>T", Filter.PASS),
            V37Call(
                GeneCoordinate("1", 97981395), "T", ("T", "C"), "DPYD", ("rs1801159",), "293T>C", Filter.PASS),
            V37Call(
                GeneCoordinate("1", 97915621), "TG", ("TG", "TC"), "DPYD", ("rs72549303",), "6744CA>GA", Filter.PASS),
        }))
        pgx_analysis = PgxAnalyser.create_pgx_analysis(ids_found_in_patient, panel)

        gene_to_haplotype_calls_expected: Dict[str, Set[HaplotypeCall]] = {"DPYD": set()}
        all_full_calls_expected = frozenset({
            FullCall(
                GeneCoordinate("1", 97915614), "C", GeneCoordinate("1", 97450058), "C",
                ("C", "T"), "DPYD", ("rs3918290",), "9213C>T", Filter.PASS, "9213C>T", Filter.PASS,
            ),
            FullCall(
                GeneCoordinate("1", 97915621), "TG", GeneCoordinate("1", 97450065), "TC",
                ("TG", "TC"), "DPYD", ("rs72549303",), "6744CA>GA", Filter.PASS, "6744GA>CA", Filter.PASS,
            ),
            FullCall(
                GeneCoordinate("1", 97981395), "T", GeneCoordinate("1", 97515839), "T",
                ("T", "C"), "DPYD", ("rs1801159",), "293T>C", Filter.PASS, "293T>C", Filter.PASS,
            ),
        })
        pgx_analysis_expected = PgxAnalysis(all_full_calls_expected, gene_to_haplotype_calls_expected)
        self.assertEqual(pgx_analysis_expected, pgx_analysis)

    def test_unresolved_haplotype_because_of_unexpected_base(self) -> None:
        """No haplotype call because of unexpected base at known variant location"""
        panel = self.__get_narrow_example_panel({"*2A", "*5", "*2B"})
        ids_found_in_patient = V37CallData(frozenset({
            V37Call(
                GeneCoordinate("1", 97915614), "C", ("T", "A"), "DPYD", ("rs3918290",), "9213C>T", Filter.PASS),
            V37Call(
                GeneCoordinate("1", 97981395), "T", ("C", "C"), "DPYD", ("rs1801159",), "293T>C", Filter.PASS),
            V37Call(
                GeneCoordinate("1", 97915621), "TG", ("TC", "TC"), "DPYD", ("rs72549303",), "6744CA>GA", Filter.PASS),
        }))
        pgx_analysis = PgxAnalyser.create_pgx_analysis(ids_found_in_patient, panel)

        gene_to_haplotype_calls_expected: Dict[str, Set[HaplotypeCall]] = {"DPYD": set()}
        all_full_calls_expected = frozenset({
            FullCall(
                GeneCoordinate("1", 97915614), "C", GeneCoordinate("1", 97450058), "C",
                ("T", "A"), "DPYD", ("rs3918290",), "9213C>T", Filter.PASS, "9213C>T", Filter.PASS,
            ),
            FullCall(
                GeneCoordinate("1", 97915621), "TG", GeneCoordinate("1", 97450065), "TC",
                ("TC", "TC"), "DPYD", ("rs72549303",), "6744CA>GA", Filter.PASS, "REF_CALL", Filter.PASS,
            ),
            FullCall(
                GeneCoordinate("1", 97981395), "T", GeneCoordinate("1", 97515839), "T",
                ("C", "C"), "DPYD", ("rs1801159",), "293T>C", Filter.PASS, "293T>C", Filter.PASS,
            ),
        })
        pgx_analysis_expected = PgxAnalysis(all_full_calls_expected, gene_to_haplotype_calls_expected)
        self.assertEqual(pgx_analysis_expected, pgx_analysis)

    def test_unresolved_haplotype_because_of_only_half_of_haplotype(self) -> None:
        """No haplotype call because of missing half of haplotype. One half is present twice, other once"""
        panel = self.__get_narrow_example_panel({"*2B"})
        ids_found_in_patient = V37CallData(frozenset({
            V37Call(
                GeneCoordinate("1", 97915614), "C", ("T", "T"), "DPYD", ("rs3918290",), "9213C>T", Filter.PASS),
            V37Call(
                GeneCoordinate("1", 97981395), "T", ("T", "C"), "DPYD", ("rs1801159",), "293T>C", Filter.PASS),
            V37Call(
                GeneCoordinate("1", 97915621), "TG", ("TC", "TC"), "DPYD", ("rs72549303",), "6744CA>GA", Filter.PASS),
        }))
        pgx_analysis = PgxAnalyser.create_pgx_analysis(ids_found_in_patient, panel)

        gene_to_haplotype_calls_expected: Dict[str, Set[HaplotypeCall]] = {"DPYD": set()}
        all_full_calls_expected = frozenset({
            FullCall(
                GeneCoordinate("1", 97915614), "C", GeneCoordinate("1", 97450058), "C",
                ("T", "T"), "DPYD", ("rs3918290",), "9213C>T", Filter.PASS, "9213C>T", Filter.PASS,
            ),
            FullCall(
                GeneCoordinate("1", 97915621), "TG", GeneCoordinate("1", 97450065), "TC",
                ("TC", "TC"), "DPYD", ("rs72549303",), "6744CA>GA", Filter.PASS, "REF_CALL", Filter.PASS,
            ),
            FullCall(
                GeneCoordinate("1", 97981395), "T", GeneCoordinate("1", 97515839), "T",
                ("T", "C"), "DPYD", ("rs1801159",), "293T>C", Filter.PASS, "293T>C", Filter.PASS,
            ),
        })
        pgx_analysis_expected = PgxAnalysis(all_full_calls_expected, gene_to_haplotype_calls_expected)
        self.assertEqual(pgx_analysis_expected, pgx_analysis)

    def test_unresolved_haplotype_because_mnv_covers_snv_starting_early(self) -> None:
        """Unresolved haplotype because MNV covers where SNV was expected, where MNV starts before this location"""
        panel = self.__get_narrow_example_panel({"*2A", "*5", "*2B"})
        ids_found_in_patient = V37CallData(frozenset({
            V37Call(
                GeneCoordinate("1", 97915613), "GC", ("CT", "CT"), "DPYD", (".",), "9212GC>CT", Filter.PASS),  # unknown
            V37Call(
                GeneCoordinate("1", 97981395), "T", ("C", "C"), "DPYD", ("rs1801159",), "293T>C", Filter.PASS),
            V37Call(
                GeneCoordinate("1", 97915621), "TG", ("TC", "TC"), "DPYD", ("rs72549303",), "6744CA>GA", Filter.PASS),
        }))
        pgx_analysis = PgxAnalyser.create_pgx_analysis(ids_found_in_patient, panel)

        gene_to_haplotype_calls_expected: Dict[str, Set[HaplotypeCall]] = {"DPYD": set()}
        all_full_calls_expected = frozenset({
            FullCall(
                GeneCoordinate("1", 97915613), "GC", None, None,
                ("CT", "CT"), "DPYD", (".",), "9212GC>CT", Filter.PASS, "9212GC>CT?", Filter.UNKNOWN,
            ),
            FullCall(
                GeneCoordinate("1", 97915621), "TG", GeneCoordinate("1", 97450065), "TC",
                ("TC", "TC"), "DPYD", ("rs72549303",), "6744CA>GA", Filter.PASS, "REF_CALL", Filter.PASS,
            ),
            FullCall(
                GeneCoordinate("1", 97981395), "T", GeneCoordinate("1", 97515839), "T",
                ("C", "C"), "DPYD", ("rs1801159",), "293T>C", Filter.PASS, "293T>C", Filter.PASS,
            ),
        })
        pgx_analysis_expected = PgxAnalysis(all_full_calls_expected, gene_to_haplotype_calls_expected)
        self.assertEqual(pgx_analysis_expected, pgx_analysis)

    def test_unresolved_haplotype_because_mnv_covers_snv_starting_there(self) -> None:
        """Unresolved haplotype because MNV covers where SNV was expected, where MNV starts at this location"""
        panel = self.__get_narrow_example_panel({"*2A", "*5", "*2B"})
        ids_found_in_patient = V37CallData(frozenset({
            V37Call(
                GeneCoordinate("1", 97915614), "CG", ("TC", "TC"), "DPYD", (".",), "9212CG>TC", Filter.PASS),  # unknown
            V37Call(
                GeneCoordinate("1", 97981395), "T", ("C", "C"), "DPYD", ("rs1801159",), "293T>C", Filter.PASS),
            V37Call(
                GeneCoordinate("1", 97915621), "TG", ("TC", "TC"), "DPYD", ("rs72549303",), "6744CA>GA", Filter.PASS),
        }))
        pgx_analysis = PgxAnalyser.create_pgx_analysis(ids_found_in_patient, panel)

        gene_to_haplotype_calls_expected: Dict[str, Set[HaplotypeCall]] = {"DPYD": set()}
        all_full_calls_expected = frozenset({
            FullCall(
                GeneCoordinate("1", 97915614), "CG", None, None,
                ("TC", "TC"), "DPYD", (".",), "9212CG>TC", Filter.PASS, "9212CG>TC?", Filter.UNKNOWN,
            ),
            FullCall(
                GeneCoordinate("1", 97915621), "TG", GeneCoordinate("1", 97450065), "TC",
                ("TC", "TC"), "DPYD", ("rs72549303",), "6744CA>GA", Filter.PASS, "REF_CALL", Filter.PASS,
            ),
            FullCall(
                GeneCoordinate("1", 97981395), "T", GeneCoordinate("1", 97515839), "T",
                ("C", "C"), "DPYD", ("rs1801159",), "293T>C", Filter.PASS, "293T>C", Filter.PASS,
            ),
        })
        pgx_analysis_expected = PgxAnalysis(all_full_calls_expected, gene_to_haplotype_calls_expected)
        self.assertEqual(pgx_analysis_expected, pgx_analysis)

    def test_overlapping_mnv_and_snv(self) -> None:
        """Overlapping SNV and unrecognized MNV do not cause errors"""
        panel = self.__get_narrow_example_panel({"*2A", "*5", "*2B"})
        ids_found_in_patient = V37CallData(frozenset({
            V37Call(
                GeneCoordinate("1", 97915614), "C", ("T", "T"), "DPYD", ("rs3918290",), "9212C>T", Filter.PASS),
            V37Call(
                GeneCoordinate("1", 97915614), "CG", ("CG", "TC"), "DPYD", (".",), "9212CG>TC", Filter.PASS),  # unknown
            V37Call(
                GeneCoordinate("1", 97981395), "T", ("C", "C"), "DPYD", ("rs1801159",), "293T>C", Filter.PASS),
            V37Call(
                GeneCoordinate("1", 97915621), "TG", ("TC", "TC"), "DPYD", ("rs72549303",), "6744CA>GA", Filter.PASS),
        }))
        pgx_analysis = PgxAnalyser.create_pgx_analysis(ids_found_in_patient, panel)

        gene_to_haplotype_calls_expected: Dict[str, Set[HaplotypeCall]] = {"DPYD": set()}
        all_full_calls_expected = frozenset({
            FullCall(
                GeneCoordinate("1", 97915614), "C", GeneCoordinate("1", 97450058), "C",
                ("T", "T"), "DPYD", ("rs3918290",), "9212C>T", Filter.PASS, "9212C>T", Filter.PASS,
            ),
            FullCall(
                GeneCoordinate("1", 97915614), "CG", None, None,
                ("CG", "TC"), "DPYD", (".",), "9212CG>TC", Filter.PASS, "9212CG>TC?", Filter.UNKNOWN,
            ),
            FullCall(
                GeneCoordinate("1", 97915621), "TG", GeneCoordinate("1", 97450065), "TC",
                ("TC", "TC"), "DPYD", ("rs72549303",), "6744CA>GA", Filter.PASS, "REF_CALL", Filter.PASS,
            ),
            FullCall(
                GeneCoordinate("1", 97981395), "T", GeneCoordinate("1", 97515839), "T",
                ("C", "C"), "DPYD", ("rs1801159",), "293T>C", Filter.PASS, "293T>C", Filter.PASS,
            ),
        })
        pgx_analysis_expected = PgxAnalysis(all_full_calls_expected, gene_to_haplotype_calls_expected)
        self.assertEqual(pgx_analysis_expected, pgx_analysis)

    def test_overlapping_mnv_and_snv_at_ref_seq_difference(self) -> None:
        """Overlapping SNV and unrecognized MNV at ref seq difference do not cause error"""
        panel = self.__get_narrow_example_panel({"*2A", "*5", "*2B"})
        ids_found_in_patient = V37CallData(frozenset({
            V37Call(
                GeneCoordinate("1", 97915614), "C", ("T", "T"), "DPYD", ("rs3918290",), "9212C>T", Filter.PASS),
            V37Call(
                GeneCoordinate("1", 97981395), "T", ("C", "C"), "DPYD", ("rs1801159",), "293T>C", Filter.PASS),
            V37Call(
                GeneCoordinate("1", 97915621), "TG", ("TG", "AC"), "DPYD", ("rs72549303",), "6744CA>GT", Filter.PASS),
            V37Call(
                GeneCoordinate("1", 97915622), "G", ("C", "C"), "DPYD", (".",), "6744C>G", Filter.PASS),  # unknown
        }))
        pgx_analysis = PgxAnalyser.create_pgx_analysis(ids_found_in_patient, panel)

        gene_to_haplotype_calls_expected: Dict[str, Set[HaplotypeCall]] = {"DPYD": set()}
        all_full_calls_expected = frozenset({
            FullCall(
                GeneCoordinate("1", 97915614), "C", GeneCoordinate("1", 97450058), "C",
                ("T", "T"), "DPYD", ("rs3918290",), "9212C>T", Filter.PASS, "9212C>T", Filter.PASS,
            ),
            FullCall(
                GeneCoordinate("1", 97915621), "TG", GeneCoordinate("1", 97450065), "TC",
                ("TG", "AC"), "DPYD", ("rs72549303",), "6744CA>GT", Filter.PASS, '6744CA>GT?', Filter.UNKNOWN,
            ),
            FullCall(
                GeneCoordinate("1", 97915622), "G", None, None,
                ("C", "C"), "DPYD", (".",), "6744C>G", Filter.PASS, "6744C>G?", Filter.UNKNOWN,
            ),
            FullCall(
                GeneCoordinate("1", 97981395), "T", GeneCoordinate("1", 97515839), "T",
                ("C", "C"), "DPYD", ("rs1801159",), "293T>C", Filter.PASS, "293T>C", Filter.PASS,
            ),
        })
        pgx_analysis_expected = PgxAnalysis(all_full_calls_expected, gene_to_haplotype_calls_expected)
        self.assertEqual(pgx_analysis_expected, pgx_analysis)


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
