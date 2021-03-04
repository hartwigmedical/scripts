import unittest
from typing import Dict, Set

import pandas as pd

from base.gene_coordinate import GeneCoordinate
from call_data import Grch37CallData, Grch37Call
from config.drug_info import DrugInfo
from config.gene_info import GeneInfo
from config.haplotype import Haplotype
from config.panel import Panel
from config.rs_id_info import RsIdInfo
from config.variant import Variant
from pgx_analysis import PgxAnalyser

ALL_IDS_IN_PANEL_COLUMNS = [
    'gene', 'position_GRCh37', 'ref_GRCh37', 'alt_GRCh37', 'position_GRCh38', 'ref_GRCh38', 'alt_GRCh38',
    'rsid', 'variant_annotation', 'filter'
]


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
            GeneInfo("DPYD", "1", "GRCh37", "*1", dpyd_haplotypes, dpyd_rs_id_infos,
                     dpyd_drugs, dpyd_rs_id_to_difference_annotations),
            GeneInfo("FAKE", "5", "GRCh37", "*1", fake_haplotypes, fake_rs_id_infos,
                     fake_drugs, fake_rs_id_to_difference_annotations),
            GeneInfo("FAKE2", "16", "GRCh37", "*1", fake2_haplotypes, fake2_rs_id_infos,
                     fake2_drugs, fake2_rs_id_to_difference_annotations),
        })
        return Panel(gene_infos)

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
            GeneInfo("DPYD", "1", "GRCh37", "*1", included_dpyd_haplotypes, included_dpyd_rs_id_infos,
                     dpyd_drugs, dpyd_rs_id_to_difference_annotations),
        })
        return Panel(gene_infos)

    def test_empty(self) -> None:
        """No variants wrt GRCh37"""
        panel = self.__get_wide_example_panel()
        ids_found_in_patient = Grch37CallData(tuple())
        results, panel_calls_for_patient = PgxAnalyser.create_pgx_analysis(ids_found_in_patient, panel)

        panel_calls_for_patient_expected = pd.DataFrame(
            [
                ("FAKE2", "16:97915617", "C", "C", "16:97450060", "C", "C", "rs1212127", "1324T>C", "INFERRED_REF_CALL"),
                ("DPYD", "1:97915614", "C", "C", "1:97450058", "C", "C", "rs3918290", "REF_CALL", "NO_CALL"),
                ("DPYD", "1:97915621", "TG", "TG", "1:97450065", "TG", "TG", "rs72549303", "6744GA>CA", "INFERRED_REF_CALL"),
                ("DPYD", "1:97981395", "T", "T", "1:97515839", "T", "T", "rs1801159", "REF_CALL", "NO_CALL"),
                ("DPYD", "1:98205966", "GATGA", "GATGA", "1:97740410", "GATGA", "GATGA", "rs72549309", "REF_CALL", "NO_CALL"),
                ("FAKE", "5:97915617", "T", "T", "5:97450060", "T", "T", "rs1212125", "REF_CALL", "NO_CALL"),
            ], columns=ALL_IDS_IN_PANEL_COLUMNS
        )
        pd.testing.assert_frame_equal(panel_calls_for_patient_expected, panel_calls_for_patient)

        results_expected = {"DPYD": {"*3_HOM"}, "FAKE": {"*1_HOM"}, "FAKE2": {"*4A_HOM"}}
        self.assertEqual(results_expected, results)

    def test_hom_ref(self) -> None:
        """All haplotypes are  *1_HOM"""
        panel = self.__get_wide_example_panel()
        ids_found_in_patient = Grch37CallData((
            Grch37Call(GeneCoordinate("16", 97915617), "C", ("T", "T"), "FAKE2", ("rs1212127",), "1324C>T", "PASS"),
            Grch37Call(GeneCoordinate("1", 97915621), "TG", ("TC", "TC"), "DPYD", ("rs72549303",), "6744CA>GA", "PASS"),
            Grch37Call(GeneCoordinate("1", 97915614), "C", ("C", "C"), "DPYD", ("rs3918290",), "REF_CALL", "PASS"),
        ))
        results, panel_calls_for_patient = PgxAnalyser.create_pgx_analysis(ids_found_in_patient, panel)

        panel_calls_for_patient_expected = pd.DataFrame(
            [
                ("FAKE2", "16:97915617", "T", "T", "16:97450060", "T", "T", "rs1212127", "REF_CALL", "NO_CALL"),
                ("DPYD", "1:97915614", "C", "C", "1:97450058", "C", "C", "rs3918290", "REF_CALL", "PASS"),
                ("DPYD", "1:97915621", "TC", "TC", "1:97450065", "TC", "TC", "rs72549303", "REF_CALL", "NO_CALL"),
                ("DPYD", "1:97981395", "T", "T", "1:97515839", "T", "T", "rs1801159", "REF_CALL", "NO_CALL"),
                ("DPYD", "1:98205966", "GATGA", "GATGA", "1:97740410", "GATGA", "GATGA", "rs72549309", "REF_CALL", "NO_CALL"),
                ("FAKE", "5:97915617", "T", "T", "5:97450060", "T", "T", "rs1212125", "REF_CALL", "NO_CALL"),
            ], columns=ALL_IDS_IN_PANEL_COLUMNS
        )
        pd.testing.assert_frame_equal(panel_calls_for_patient_expected, panel_calls_for_patient)

        results_expected = {"DPYD": {"*1_HOM"}, "FAKE": {"*1_HOM"}, "FAKE2": {"*1_HOM"}}
        self.assertEqual(results_expected, results)

    def test_heterozygous(self) -> None:
        """All haplotypes are heterozygous. Both variant/ref and variant/variant"""
        panel = self.__get_wide_example_panel()
        ids_found_in_patient = Grch37CallData((
            Grch37Call(GeneCoordinate("16", 97915617), "C", ("C", "T"), "FAKE2", ("rs1212127",), "1324C>T", "PASS"),
            Grch37Call(GeneCoordinate("5", 97915617), "T", ("T", "C"), "FAKE", ("rs1212125",), "1005T>C", "PASS"),
            Grch37Call(GeneCoordinate("1", 97915621), "TG", ("TG", "TC"), "DPYD", ("rs72549303",), "6744CA>GA", "PASS"),
            Grch37Call(GeneCoordinate("1", 97915614), "C", ("C", "T"), "DPYD", ("rs3918290",), "35G>A", "PASS"),
            Grch37Call(GeneCoordinate("1", 97981395), "T", ("T", "C"), "DPYD", ("rs1801159",), "674A>G", "PASS"),
        ))
        results, panel_calls_for_patient = PgxAnalyser.create_pgx_analysis(ids_found_in_patient, panel)

        panel_calls_for_patient_expected = pd.DataFrame(
            [
                ("FAKE2", "16:97915617", "C", "T", "16:97450060", "T", "C", "rs1212127", "1324T>C", "PASS"),
                ("DPYD", "1:97915614", "C", "T", "1:97450058", "C", "T", "rs3918290", "35G>A", "PASS"),
                ("DPYD", "1:97915621", "TG", "TC", "1:97450065", "TC", "TG", "rs72549303", "6744GA>CA", "PASS"),
                ("DPYD", "1:97981395", "T", "C", "1:97515839", "T", "C", "rs1801159", "674A>G", "PASS"),
                ("DPYD", "1:98205966", "GATGA", "GATGA", "1:97740410", "GATGA", "GATGA", "rs72549309", "REF_CALL", "NO_CALL"),
                ("FAKE", "5:97915617", "T", "C", "5:97450060", "T", "C", "rs1212125", "1005T>C", "PASS"),
            ], columns=ALL_IDS_IN_PANEL_COLUMNS
        )
        pd.testing.assert_frame_equal(panel_calls_for_patient_expected, panel_calls_for_patient)

        results_expected = {"DPYD": {"*2B_HET", "*3_HET"}, "FAKE": {"*4A_HET", "*1_HET"}, "FAKE2": {"*4A_HET", "*1_HET"}}
        self.assertEqual(results_expected, results)

    def test_ref_call_on_ref_seq_differences(self) -> None:
        """Explicit ref calls wrt GRCh37 at differences between GRCh37 and GRCh38"""
        panel = self.__get_wide_example_panel()
        ids_found_in_patient = Grch37CallData((
            Grch37Call(GeneCoordinate("16", 97915617), "C", ("C", "C"), "FAKE2", ("rs1212127",), "REF_CALL", "PASS"),
            Grch37Call(GeneCoordinate("1", 97915621), "TG", ("TG", "TG"), "DPYD", ("rs72549303",), "REF_CALL", "PASS"),
        ))
        results, panel_calls_for_patient = PgxAnalyser.create_pgx_analysis(ids_found_in_patient, panel)

        panel_calls_for_patient_expected = pd.DataFrame(
            [
                ("FAKE2", "16:97915617", "C", "C", "16:97450060", "C", "C", "rs1212127", "1324T>C", "PASS"),
                ("DPYD", "1:97915614", "C", "C", "1:97450058", "C", "C", "rs3918290", "REF_CALL", "NO_CALL"),
                ("DPYD", "1:97915621", "TG", "TG", "1:97450065", "TG", "TG", "rs72549303", "6744GA>CA", "PASS"),
                ("DPYD", "1:97981395", "T", "T", "1:97515839", "T", "T", "rs1801159", "REF_CALL", "NO_CALL"),
                ("DPYD", "1:98205966", "GATGA", "GATGA", "1:97740410", "GATGA", "GATGA", "rs72549309", "REF_CALL", "NO_CALL"),
                ("FAKE", "5:97915617", "T", "T", "5:97450060", "T", "T", "rs1212125", "REF_CALL", "NO_CALL"),
            ], columns=ALL_IDS_IN_PANEL_COLUMNS
        )
        pd.testing.assert_frame_equal(panel_calls_for_patient_expected, panel_calls_for_patient)

        results_expected = {"DPYD": {"*3_HOM"}, "FAKE": {"*1_HOM"}, "FAKE2": {"*4A_HOM"}}
        self.assertEqual(results_expected, results)

    def test_only_position_match_on_ref_seq_differences(self) -> None:
        """At reference sequence differences: heterozygous between ref GRCh37 and GRCh38, and no rs_id provided"""
        panel = self.__get_wide_example_panel()
        ids_found_in_patient = Grch37CallData((
            Grch37Call(GeneCoordinate("16", 97915617), "C", ("C", "T"), "FAKE2", (".",), "1324C>T", "PASS"),
            Grch37Call(GeneCoordinate("1", 97915621), "TG", ("TG", "TC"), "DPYD", (".",), "6744CA>GA", "PASS"),
        ))
        results, panel_calls_for_patient = PgxAnalyser.create_pgx_analysis(ids_found_in_patient, panel)

        panel_calls_for_patient_expected = pd.DataFrame(
            [
                ("FAKE2", "16:97915617", "C", "T", "16:97450060", "T", "C", "rs1212127", "1324T>C", "PASS"),
                ("DPYD", "1:97915614", "C", "C", "1:97450058", "C", "C", "rs3918290", "REF_CALL", "NO_CALL"),
                ("DPYD", "1:97915621", "TG", "TC", "1:97450065", "TC", "TG", "rs72549303", "6744GA>CA", "PASS"),
                ("DPYD", "1:97981395", "T", "T", "1:97515839", "T", "T", "rs1801159", "REF_CALL", "NO_CALL"),
                ("DPYD", "1:98205966", "GATGA", "GATGA", "1:97740410", "GATGA", "GATGA", "rs72549309", "REF_CALL", "NO_CALL"),
                ("FAKE", "5:97915617", "T", "T", "5:97450060", "T", "T", "rs1212125", "REF_CALL", "NO_CALL"),
            ], columns=ALL_IDS_IN_PANEL_COLUMNS
        )
        pd.testing.assert_frame_equal(panel_calls_for_patient_expected, panel_calls_for_patient)

        results_expected = {"DPYD": {'*3_HET', '*1_HET'}, "FAKE": {"*1_HOM"}, "FAKE2": {'*4A_HET', '*1_HET'}}
        self.assertEqual(results_expected, results)

    def test_wrong_rs_id_on_ref_seq_differences(self) -> None:
        """
        At reference sequence differences: heterozygous between ref GRCh37 and GRCh38,
        and incorrect rs id provided
        """
        panel = self.__get_wide_example_panel()
        ids_found_in_patient = Grch37CallData((
            Grch37Call(GeneCoordinate("16", 97915617), "C", ("C", "T"), "FAKE2", ("rs939535",), "1324C>T", "PASS"),
            Grch37Call(GeneCoordinate("1", 97915621), "TG", ("TG", "TC"), "DPYD", ("rs4020942",), "6744CA>GA", "PASS"),
        ))
        with self.assertRaises(ValueError):
            PgxAnalyser.create_pgx_analysis(ids_found_in_patient, panel)

    def test_wrong_position_on_ref_seq_differences(self) -> None:
        """Explicit ref calls wrt GRCh37 at differences between GRCh37 and GRCh38, except positions are incorrect"""
        panel = self.__get_wide_example_panel()
        ids_found_in_patient = Grch37CallData((
            Grch37Call(GeneCoordinate("16", 97915618), "C", ("C", "C"), "FAKE2", ("rs1212127",), "REF_CALL", "PASS"),
            Grch37Call(GeneCoordinate("1", 97915623), "TG", ("TG", "TG"), "DPYD", ("rs72549303",), "REF_CALL", "PASS"),
        ))
        with self.assertRaises(ValueError):
            PgxAnalyser.create_pgx_analysis(ids_found_in_patient, panel)

    def test_wrong_chromosome_on_ref_seq_differences(self) -> None:
        """Explicit ref calls wrt GRCh37 at differences between GRCh37 and GRCh38, except chromosomes are incorrect"""
        panel = self.__get_wide_example_panel()
        ids_found_in_patient = Grch37CallData((
            Grch37Call(GeneCoordinate("7", 97915617), "C", ("C", "C"), "FAKE2", ("rs1212127",), "REF_CALL", "PASS"),
            Grch37Call(GeneCoordinate("8", 97915621), "TG", ("TG", "TG"), "DPYD", ("rs72549303",), "REF_CALL", "PASS"),
        ))
        with self.assertRaises(ValueError):
            PgxAnalyser.create_pgx_analysis(ids_found_in_patient, panel)

    def test_position_of_other_variant_on_ref_seq_differences(self) -> None:
        """
        Explicit ref calls wrt GRCh37 at differences between GRCh37 and GRCh38, except position is of other variant
        """
        panel = self.__get_wide_example_panel()
        ids_found_in_patient = Grch37CallData((
            Grch37Call(GeneCoordinate("16", 97915617), "C", ("C", "C"), "FAKE2", ("rs1212127",), "REF_CALL", "PASS"),
            Grch37Call(GeneCoordinate("1", 98205966), "TG", ("TG", "TG"), "DPYD", ("rs72549303",), "REF_CALL", "PASS"),
        ))
        with self.assertRaises(ValueError):
            PgxAnalyser.create_pgx_analysis(ids_found_in_patient, panel)

    def test_single_different_allele_on_ref_seq_differences(self) -> None:
        """At reference sequence differences: single allele that is ref GRCh37 or GRCh38, other allele is neither"""
        panel = self.__get_wide_example_panel()
        ids_found_in_patient = Grch37CallData((
            Grch37Call(GeneCoordinate("16", 97915617), "C", ("C", "A"), "FAKE2", ("rs1212127",), "1324C>A", "PASS"),  # with ref GRCh37
            Grch37Call(GeneCoordinate("1", 97915621), "TG", ("AC", "TC"), "DPYD", ("rs72549303",), "6744CT>GT;6744CT>GC", "PASS"),  # with ref GRCh38
        ))
        results, panel_calls_for_patient = PgxAnalyser.create_pgx_analysis(ids_found_in_patient, panel)

        panel_calls_for_patient_expected = pd.DataFrame(
            [
                ("FAKE2", "16:97915617", "C", "A", "16:97450060", "C", "A", "rs1212127", "1324C>A?", "PASS"),
                ("DPYD", "1:97915614", "C", "C", "1:97450058", "C", "C", "rs3918290", "REF_CALL", "NO_CALL"),
                ("DPYD", "1:97915621", "AC", "TC", "1:97450065", "TC", "AC", "rs72549303", "6744CT>GT;6744CT>GC?", "PASS"),
                ("DPYD", "1:97981395", "T", "T", "1:97515839", "T", "T", "rs1801159", "REF_CALL", "NO_CALL"),
                ("DPYD", "1:98205966", "GATGA", "GATGA", "1:97740410", "GATGA", "GATGA", "rs72549309", "REF_CALL", "NO_CALL"),
                ("FAKE", "5:97915617", "T", "T", "5:97450060", "T", "T", "rs1212125", "REF_CALL", "NO_CALL"),
            ], columns=ALL_IDS_IN_PANEL_COLUMNS
        )
        pd.testing.assert_frame_equal(panel_calls_for_patient_expected, panel_calls_for_patient)

        results_expected = {"DPYD": {'Unresolved_Haplotype'}, "FAKE": {"*1_HOM"}, "FAKE2": {'Unresolved_Haplotype'}}
        self.assertEqual(results_expected, results)

    def test_double_different_allele_on_ref_seq_differences(self) -> None:
        """At reference sequence differences: both alleles not ref GRCh37 or GRCh38"""
        panel = self.__get_wide_example_panel()
        ids_found_in_patient = Grch37CallData((
            Grch37Call(GeneCoordinate("16", 97915617), "C", ("A", "G"), "FAKE2", ("rs1212127",), "1324C>A;1324C>G", "PASS"),
            Grch37Call(GeneCoordinate("1", 97915621), "TG", ("AC", "AG"), "DPYD", ("rs72549303",), "6744CT>GT;6744CT>GC", "PASS"),
        ))
        results, panel_calls_for_patient = PgxAnalyser.create_pgx_analysis(ids_found_in_patient, panel)

        panel_calls_for_patient_expected = pd.DataFrame(
            [
                ("FAKE2", "16:97915617", "A", "G", "16:97450060", "A", "G", "rs1212127", "1324C>A;1324C>G?", "PASS"),
                ("DPYD", "1:97915614", "C", "C", "1:97450058", "C", "C", "rs3918290", "REF_CALL", "NO_CALL"),
                ("DPYD", "1:97915621", "AC", "AG", "1:97450065", "AC", "AG", "rs72549303", "6744CT>GT;6744CT>GC?", "PASS"),
                ("DPYD", "1:97981395", "T", "T", "1:97515839", "T", "T", "rs1801159", "REF_CALL", "NO_CALL"),
                ("DPYD", "1:98205966", "GATGA", "GATGA", "1:97740410", "GATGA", "GATGA", "rs72549309", "REF_CALL", "NO_CALL"),
                ("FAKE", "5:97915617", "T", "T", "5:97450060", "T", "T", "rs1212125", "REF_CALL", "NO_CALL"),
            ], columns=ALL_IDS_IN_PANEL_COLUMNS
        )
        pd.testing.assert_frame_equal(panel_calls_for_patient_expected, panel_calls_for_patient)

        results_expected = {"DPYD": {'Unresolved_Haplotype'}, "FAKE": {"*1_HOM"}, "FAKE2": {'Unresolved_Haplotype'}}
        self.assertEqual(results_expected, results)

    def test_unknown_variants(self) -> None:
        """Variants that are completely unknown, including with unknown rs id"""
        panel = self.__get_wide_example_panel()
        ids_found_in_patient = Grch37CallData((
            Grch37Call(GeneCoordinate("16", 39593405), "A", ("A", "G"), "FAKE2", ("rs1949223",), "384C>T", "PASS"),  # unknown
            Grch37Call(GeneCoordinate("5", 97915617), "T", ("T", "C"), "FAKE", ("rs1212125",), "1005T>C", "PASS"),
            Grch37Call(GeneCoordinate("1", 2488242), "AC", ("AC", "AG"), "DPYD", (".",), "9213CT>GT", "PASS"),  # unknown
        ))
        results, panel_calls_for_patient = PgxAnalyser.create_pgx_analysis(ids_found_in_patient, panel)

        panel_calls_for_patient_expected = pd.DataFrame(
            [
                ("FAKE2", "16:39593405", "A", "G", "UNKNOWN", "A", "G", "rs1949223", "384C>T", "PASS"),
                ("FAKE2", "16:97915617", "C", "C", "16:97450060", "C", "C", "rs1212127", "1324T>C", "INFERRED_REF_CALL"),
                ("DPYD", "1:2488242", "AC", "AG", "UNKNOWN", "AC", "AG", ".", "9213CT>GT", "PASS"),
                ("DPYD", "1:97915614", "C", "C", "1:97450058", "C", "C", "rs3918290", "REF_CALL", "NO_CALL"),
                ("DPYD", "1:97915621", "TG", "TG", "1:97450065", "TG", "TG", "rs72549303", "6744GA>CA", "INFERRED_REF_CALL"),
                ("DPYD", "1:97981395", "T", "T", "1:97515839", "T", "T", "rs1801159", "REF_CALL", "NO_CALL"),
                ("DPYD", "1:98205966", "GATGA", "GATGA", "1:97740410", "GATGA", "GATGA", "rs72549309", "REF_CALL", "NO_CALL"),
                ("FAKE", "5:97915617", "T", "C", "5:97450060", "T", "C", "rs1212125", "1005T>C", "PASS"),
            ], columns=ALL_IDS_IN_PANEL_COLUMNS
        )
        pd.testing.assert_frame_equal(panel_calls_for_patient_expected, panel_calls_for_patient)

        results_expected = {
            "DPYD": {'Unresolved_Haplotype'},
            "FAKE": {"*4A_HET", "*1_HET"},
            "FAKE2": {'Unresolved_Haplotype'},
        }
        self.assertEqual(results_expected, results)

    def test_unknown_gene(self) -> None:
        """Variants that are of an unknown gene"""
        panel = self.__get_wide_example_panel()
        ids_found_in_patient = Grch37CallData((
            Grch37Call(GeneCoordinate("5", 97915617), "T", ("T", "C"), "FAKE", ("rs1212125",), "1005T>C", "PASS"),
            Grch37Call(GeneCoordinate("3", 18473423), "T", ("T", "C"), ".", ("rs2492932",), "12T>C", "PASS"),  # unknown
        ))
        with self.assertRaises(ValueError):
            PgxAnalyser.create_pgx_analysis(ids_found_in_patient, panel)

    def test_known_variants_with_incorrect_or_missing_rs_id(self) -> None:
        """Known variants (not ref seq differences) with incorrect or unknown rs id"""
        panel = self.__get_wide_example_panel()
        ids_found_in_patient = Grch37CallData((
            Grch37Call(GeneCoordinate("16", 97915617), "C", ("C", "T"), "FAKE2", ("rs1212127",), "1324C>T", "PASS"),
            Grch37Call(GeneCoordinate("5", 97915617), "T", ("T", "C"), "FAKE", ("rs27384",), "1005T>C", "PASS"),  # incorrect
            Grch37Call(GeneCoordinate("1", 97915621), "TG", ("TG", "TC"), "DPYD", ("rs72549303",), "6744CA>GA", "PASS"),
            Grch37Call(GeneCoordinate("1", 97915614), "C", ("C", "T"), "DPYD", (".",), "35G>A", "PASS"),  # missing
            Grch37Call(GeneCoordinate("1", 97981395), "T", ("T", "C"), "DPYD", ("rs1801159",), "674A>G", "PASS"),
        ))
        with self.assertRaises(ValueError):
            PgxAnalyser.create_pgx_analysis(ids_found_in_patient, panel)

    def test_known_variant_with_incorrect_position(self) -> None:
        """Known variants (not ref seq differences), with one incorrect position"""
        panel = self.__get_wide_example_panel()
        ids_found_in_patient = Grch37CallData((
            Grch37Call(GeneCoordinate("16", 97915617), "C", ("C", "T"), "FAKE2", ("rs1212127",), "1324C>T", "PASS"),
            Grch37Call(GeneCoordinate("5", 97915617), "T", ("T", "C"), "FAKE", ("rs1212125",), "1005T>C", "PASS"),
            Grch37Call(GeneCoordinate("1", 97915621), "TG", ("TG", "TC"), "DPYD", ("rs72549303",), "6744CA>GA", "PASS"),
            Grch37Call(GeneCoordinate("1", 6778543), "C", ("C", "T"), "DPYD", ("rs3918290",), "35G>A", "PASS"),  # incorrect
            Grch37Call(GeneCoordinate("1", 97981395), "T", ("T", "C"), "DPYD", ("rs1801159",), "674A>G", "PASS"),
        ))
        with self.assertRaises(ValueError):
            PgxAnalyser.create_pgx_analysis(ids_found_in_patient, panel)

    def test_known_variant_with_incorrect_chromosome(self) -> None:
        """Known variants (not ref seq differences), with one incorrect chromosome"""
        panel = self.__get_wide_example_panel()
        ids_found_in_patient = Grch37CallData((
            Grch37Call(GeneCoordinate("16", 97915617), "C", ("C", "T"), "FAKE2", ("rs1212127",), "1324C>T", "PASS"),
            Grch37Call(GeneCoordinate("5", 97915617), "T", ("T", "C"), "FAKE", ("rs1212125",), "1005T>C", "PASS"),
            Grch37Call(GeneCoordinate("1", 97915621), "TG", ("TG", "TC"), "DPYD", ("rs72549303",), "6744CA>GA", "PASS"),
            Grch37Call(GeneCoordinate("3", 97915614), "C", ("C", "T"), "DPYD", ("rs3918290",), "35G>A", "PASS"),  # incorrect
            Grch37Call(GeneCoordinate("1", 97981395), "T", ("T", "C"), "DPYD", ("rs1801159",), "674A>G", "PASS"),
        ))
        with self.assertRaises(ValueError):
            PgxAnalyser.create_pgx_analysis(ids_found_in_patient, panel)

    def test_known_variant_with_multiple_rs_ids_not_matching_panel(self) -> None:
        """Multiple rs ids when panel says there should be one"""
        panel = self.__get_wide_example_panel()
        ids_found_in_patient = Grch37CallData((
            Grch37Call(GeneCoordinate("16", 97915617), "C", ("C", "T"), "FAKE2", ("rs1212127", "rs394832"), "1324C>T", "PASS"),  # incorrect
            Grch37Call(GeneCoordinate("5", 97915617), "T", ("T", "C"), "FAKE", ("rs1212125",), "1005T>C", "PASS"),
            Grch37Call(GeneCoordinate("1", 97915621), "TG", ("TG", "TC"), "DPYD", ("rs72549303",), "6744CA>GA", "PASS"),
            Grch37Call(GeneCoordinate("1", 97915614), "C", ("C", "T"), "DPYD", ("rs3918290", "rs202093"), "35G>A", "PASS"),  # incorrect
            Grch37Call(GeneCoordinate("1", 97981395), "T", ("T", "C"), "DPYD", ("rs1801159",), "674A>G", "PASS"),
        ))
        with self.assertRaises(ValueError):
            PgxAnalyser.create_pgx_analysis(ids_found_in_patient, panel)

    def test_more_than_two_haplotypes_present(self) -> None:
        """More than two haplotypes are present, specifically three"""
        panel = self.__get_narrow_example_panel({"*2A", "*3", "*7"})
        ids_found_in_patient = Grch37CallData((
            Grch37Call(GeneCoordinate("1", 97915614), "C", ("C", "T"), "DPYD", ("rs3918290",), "9213C>T", "PASS"),
            Grch37Call(GeneCoordinate("1", 97912838), "A", ("AGT", "AGT"), "DPYD", ("rs2938101",), "293A>AGT", "PASS"),
        ))
        results, panel_calls_for_patient = PgxAnalyser.create_pgx_analysis(ids_found_in_patient, panel)

        panel_calls_for_patient_expected = pd.DataFrame(
            [
                ("DPYD", "1:97912838", "AGT", "AGT", "1:97453984", "AGT", "AGT", "rs2938101", "293A>AGT", "PASS"),
                ("DPYD", "1:97915614", "C", "T", "1:97450058", "C", "T", "rs3918290", "9213C>T", "PASS"),
                ("DPYD", "1:97915621", "TG", "TG", "1:97450065", "TG", "TG", "rs72549303", "6744GA>CA", "INFERRED_REF_CALL"),
            ], columns=ALL_IDS_IN_PANEL_COLUMNS
        )
        pd.testing.assert_frame_equal(panel_calls_for_patient_expected, panel_calls_for_patient)

        results_expected = {
            "DPYD": {'*2A_HET', '*3_HOM', '*7_HOM'},
        }
        self.assertEqual(results_expected, results)

    def test_ambiguous_haplotype_with_clear_winner_homozygous(self) -> None:
        """Ambiguous homozygous haplotype, where the simpler possibility should be preferred"""
        panel = self.__get_narrow_example_panel({"*2A", "*5", "*2B"})
        ids_found_in_patient = Grch37CallData((
            Grch37Call(GeneCoordinate("1", 97915614), "C", ("T", "T"), "DPYD", ("rs3918290",), "9213C>T", "PASS"),
            Grch37Call(GeneCoordinate("1", 97981395), "T", ("C", "C"), "DPYD", ("rs1801159",), "293T>C", "PASS"),
            Grch37Call(GeneCoordinate("1", 97915621), "TG", ("TC", "TC"), "DPYD", ("rs72549303",), "6744GA>CA", "PASS"),
        ))
        results, panel_calls_for_patient = PgxAnalyser.create_pgx_analysis(ids_found_in_patient, panel)

        panel_calls_for_patient_expected = pd.DataFrame(
            [
                ("DPYD", "1:97915614", "T", "T", "1:97450058", "T", "T", "rs3918290", "9213C>T", "PASS"),
                ("DPYD", "1:97915621", "TC", "TC", "1:97450065", "TC", "TC", "rs72549303", "REF_CALL", "NO_CALL"),
                ("DPYD", "1:97981395", "C", "C", "1:97515839", "C", "C", "rs1801159", "293T>C", "PASS"),
            ], columns=ALL_IDS_IN_PANEL_COLUMNS
        )
        pd.testing.assert_frame_equal(panel_calls_for_patient_expected, panel_calls_for_patient)

        results_expected = {
            "DPYD": {'*2B_HOM'},
        }
        self.assertEqual(results_expected, results)

    def test_ambiguous_haplotype_with_clear_winner_heterozygous(self) -> None:
        """Ambiguous heterozygous haplotype, where the simpler possibility should be preferred"""
        panel = self.__get_narrow_example_panel({"*2A", "*5", "*2B"})
        ids_found_in_patient = Grch37CallData((
            Grch37Call(GeneCoordinate("1", 97915614), "C", ("C", "T"), "DPYD", ("rs3918290",), "9213C>T", "PASS"),
            Grch37Call(GeneCoordinate("1", 97981395), "T", ("T", "C"), "DPYD", ("rs1801159",), "293T>C", "PASS"),
            Grch37Call(GeneCoordinate("1", 97915621), "TG", ("TC", "TC"), "DPYD", ("rs72549303",), "6744GA>CA", "PASS"),
        ))
        results, panel_calls_for_patient = PgxAnalyser.create_pgx_analysis(ids_found_in_patient, panel)

        panel_calls_for_patient_expected = pd.DataFrame(
            [
                ("DPYD", "1:97915614", "C", "T", "1:97450058", "C", "T", "rs3918290", "9213C>T", "PASS"),
                ("DPYD", "1:97915621", "TC", "TC", "1:97450065", "TC", "TC", "rs72549303", "REF_CALL", "NO_CALL"),
                ("DPYD", "1:97981395", "T", "C", "1:97515839", "T", "C", "rs1801159", "293T>C", "PASS"),
            ], columns=ALL_IDS_IN_PANEL_COLUMNS
        )
        pd.testing.assert_frame_equal(panel_calls_for_patient_expected, panel_calls_for_patient)

        results_expected = {
            "DPYD": {'*2B_HET', '*1_HET'},
        }
        self.assertEqual(results_expected, results)

    def test_ambiguous_haplotype_with_clear_winner_mix(self) -> None:
        """Ambiguous mix of homozygous and heterozygous haplotype, where the simplest possibility should be preferred"""
        panel = self.__get_narrow_example_panel({"*2A", "*5", "*2B"})
        ids_found_in_patient = Grch37CallData((
            Grch37Call(GeneCoordinate("1", 97915614), "C", ("T", "T"), "DPYD", ("rs3918290",), "9213C>T", "PASS"),
            Grch37Call(GeneCoordinate("1", 97981395), "T", ("T", "C"), "DPYD", ("rs1801159",), "293T>C", "PASS"),
            Grch37Call(GeneCoordinate("1", 97915621), "TG", ("TC", "TC"), "DPYD", ("rs72549303",), "6744GA>CA", "PASS"),
        ))
        results, panel_calls_for_patient = PgxAnalyser.create_pgx_analysis(ids_found_in_patient, panel)

        panel_calls_for_patient_expected = pd.DataFrame(
            [
                ("DPYD", "1:97915614", "T", "T", "1:97450058", "T", "T", "rs3918290", "9213C>T", "PASS"),
                ("DPYD", "1:97915621", "TC", "TC", "1:97450065", "TC", "TC", "rs72549303", "REF_CALL", "NO_CALL"),
                ("DPYD", "1:97981395", "T", "C", "1:97515839", "T", "C", "rs1801159", "293T>C", "PASS"),
            ], columns=ALL_IDS_IN_PANEL_COLUMNS
        )
        pd.testing.assert_frame_equal(panel_calls_for_patient_expected, panel_calls_for_patient)

        results_expected = {
            "DPYD": {'*2B_HET', '*2A_HET'},
        }
        self.assertEqual(results_expected, results)

    def test_complicated_ambiguous_haplotype_with_a_clear_winner(self) -> None:
        """Ambiguous set of haplotypes, with no haplotype that is simplest"""
        panel = self.__get_narrow_example_panel({"*2A", "*5", "*2B", "*3", "*9", "*10", "*7"})
        ids_found_in_patient = Grch37CallData((
            Grch37Call(GeneCoordinate("1", 97915614), "C", ("T", "T"), "DPYD", ("rs3918290",), "9213C>T", "PASS"),
            Grch37Call(GeneCoordinate("1", 97981395), "T", ("T", "C"), "DPYD", ("rs1801159",), "293T>C", "PASS"),
            Grch37Call(GeneCoordinate("1", 97915621), "TG", ("TG", "TC"), "DPYD", ("rs72549303",), "6744CA>GA", "PASS"),
            Grch37Call(GeneCoordinate("1", 97912838), "A", ("AGT", "AGT"), "DPYD", ("rs2938101",), "301A>AGT", "PASS"),
        ))
        results, panel_calls_for_patient = PgxAnalyser.create_pgx_analysis(ids_found_in_patient, panel)

        panel_calls_for_patient_expected = pd.DataFrame(
            [
                ("DPYD", "1:97912838", "AGT", "AGT", "1:97453984", "AGT", "AGT", "rs2938101", "301A>AGT", "PASS"),
                ("DPYD", "1:97915614", "T", "T", "1:97450058", "T", "T", "rs3918290", "9213C>T", "PASS"),
                ("DPYD", "1:97915621", "TG", "TC", "1:97450065", "TC", "TG", "rs72549303", "6744GA>CA", "PASS"),
                ("DPYD", "1:97981395", "T", "C", "1:97515839", "T", "C", "rs1801159", "293T>C", "PASS"),
            ], columns=ALL_IDS_IN_PANEL_COLUMNS
        )
        pd.testing.assert_frame_equal(panel_calls_for_patient_expected, panel_calls_for_patient)

        results_expected = {
            "DPYD": {'*2B_HET', '*9_HET', "*7_HOM"},
        }
        self.assertEqual(results_expected, results)

    def test_ambiguous_homozygous_haplotype_with_a_less_clear_winner(self) -> None:
        """Ambiguous set of homozygous haplotypes, with winning haplotype that might not be ideal"""
        panel = self.__get_narrow_example_panel({"*2A", "*5", "*2B", "*3", "*9", "*10"})
        ids_found_in_patient = Grch37CallData((
            Grch37Call(GeneCoordinate("1", 97915614), "C", ("T", "T"), "DPYD", ("rs3918290",), "9213C>T", "PASS"),
            Grch37Call(GeneCoordinate("1", 97981395), "T", ("C", "C"), "DPYD", ("rs1801159",), "293T>C", "PASS"),
            Grch37Call(GeneCoordinate("1", 97915621), "TG", ("TG", "TG"), "DPYD", ("rs72549303",), "6744CA>GA", "PASS"),
        ))
        results, panel_calls_for_patient = PgxAnalyser.create_pgx_analysis(ids_found_in_patient, panel)

        panel_calls_for_patient_expected = pd.DataFrame(
            [
                ("DPYD", "1:97915614", "T", "T", "1:97450058", "T", "T", "rs3918290", "9213C>T", "PASS"),
                ("DPYD", "1:97915621", "TG", "TG", "1:97450065", "TG", "TG", "rs72549303", "6744GA>CA", "PASS"),
                ("DPYD", "1:97981395", "C", "C", "1:97515839", "C", "C", "rs1801159", "293T>C", "PASS"),
            ], columns=ALL_IDS_IN_PANEL_COLUMNS
        )
        pd.testing.assert_frame_equal(panel_calls_for_patient_expected, panel_calls_for_patient)

        results_expected = {
            "DPYD": {'*10_HET', '*2B_HET', '*9_HET'},
        }
        self.assertEqual(results_expected, results)

    def test_ambiguous_heterozygous_haplotype_without_a_clear_winner(self) -> None:
        """
        Ambiguous set of heterozygous haplotypes, with no haplotype that is simplest.
        Could  be {*2A_HET,*10_HET}, {*5_HET,*9_HET} or {*3_HET,*2B_HET}
        """
        panel = self.__get_narrow_example_panel({"*2A", "*5", "*2B", "*3", "*9", "*10"})
        ids_found_in_patient = Grch37CallData((
            Grch37Call(GeneCoordinate("1", 97915614), "C", ("C", "T"), "DPYD", ("rs3918290",), "9213C>T", "PASS"),
            Grch37Call(GeneCoordinate("1", 97981395), "T", ("T", "C"), "DPYD", ("rs1801159",), "293T>C", "PASS"),
            Grch37Call(GeneCoordinate("1", 97915621), "TG", ("TG", "TC"), "DPYD", ("rs72549303",), "6744CA>GA", "PASS"),
        ))
        results, panel_calls_for_patient = PgxAnalyser.create_pgx_analysis(ids_found_in_patient, panel)

        panel_calls_for_patient_expected = pd.DataFrame(
            [
                ("DPYD", "1:97915614", "C", "T", "1:97450058", "C", "T", "rs3918290", "9213C>T", "PASS"),
                ("DPYD", "1:97915621", "TG", "TC", "1:97450065", "TC", "TG", "rs72549303", "6744GA>CA", "PASS"),
                ("DPYD", "1:97981395", "T", "C", "1:97515839", "T", "C", "rs1801159", "293T>C", "PASS"),
            ], columns=ALL_IDS_IN_PANEL_COLUMNS
        )
        pd.testing.assert_frame_equal(panel_calls_for_patient_expected, panel_calls_for_patient)

        results_expected = {
            "DPYD": {'Unresolved_Haplotype'},
        }
        self.assertEqual(results_expected, results)

    def test_unresolved_haplotype_because_of_unexpected_base(self) -> None:
        """No haplotype call because of unexpected base at known variant location"""
        panel = self.__get_narrow_example_panel({"*2A", "*5", "*2B"})
        ids_found_in_patient = Grch37CallData((
            Grch37Call(GeneCoordinate("1", 97915614), "C", ("T", "A"), "DPYD", ("rs3918290",), "9213C>T", "PASS"),
            Grch37Call(GeneCoordinate("1", 97981395), "T", ("C", "C"), "DPYD", ("rs1801159",), "293T>C", "PASS"),
            Grch37Call(GeneCoordinate("1", 97915621), "TG", ("TC", "TC"), "DPYD", ("rs72549303",), "6744GA>CA", "PASS"),
        ))
        results, panel_calls_for_patient = PgxAnalyser.create_pgx_analysis(ids_found_in_patient, panel)

        panel_calls_for_patient_expected = pd.DataFrame(
            [
                ("DPYD", "1:97915614", "T", "A", "1:97450058", "T", "A", "rs3918290", "9213C>T", "PASS"),
                ("DPYD", "1:97915621", "TC", "TC", "1:97450065", "TC", "TC", "rs72549303", "REF_CALL", "NO_CALL"),
                ("DPYD", "1:97981395", "C", "C", "1:97515839", "C", "C", "rs1801159", "293T>C", "PASS"),
            ], columns=ALL_IDS_IN_PANEL_COLUMNS
        )
        pd.testing.assert_frame_equal(panel_calls_for_patient_expected, panel_calls_for_patient)

        results_expected = {
            "DPYD": {'Unresolved_Haplotype'},
        }
        self.assertEqual(results_expected, results)

    def test_unresolved_haplotype_because_of_only_half_of_haplotype(self) -> None:
        """No haplotype call because of missing half of haplotype. One half is present twice, other once"""
        panel = self.__get_narrow_example_panel({"*2B"})
        ids_found_in_patient = Grch37CallData((
            Grch37Call(GeneCoordinate("1", 97915614), "C", ("T", "T"), "DPYD", ("rs3918290",), "9213C>T", "PASS"),
            Grch37Call(GeneCoordinate("1", 97981395), "T", ("T", "C"), "DPYD", ("rs1801159",), "293T>C", "PASS"),
            Grch37Call(GeneCoordinate("1", 97915621), "TG", ("TC", "TC"), "DPYD", ("rs72549303",), "6744GA>CA", "PASS"),
        ))
        results, panel_calls_for_patient = PgxAnalyser.create_pgx_analysis(ids_found_in_patient, panel)

        panel_calls_for_patient_expected = pd.DataFrame(
            [
                ("DPYD", "1:97915614", "T", "T", "1:97450058", "T", "T", "rs3918290", "9213C>T", "PASS"),
                ("DPYD", "1:97915621", "TC", "TC", "1:97450065", "TC", "TC", "rs72549303", "REF_CALL", "NO_CALL"),
                ("DPYD", "1:97981395", "T", "C", "1:97515839", "T", "C", "rs1801159", "293T>C", "PASS"),
            ], columns=ALL_IDS_IN_PANEL_COLUMNS
        )
        pd.testing.assert_frame_equal(panel_calls_for_patient_expected, panel_calls_for_patient)

        results_expected = {
            "DPYD": {'Unresolved_Haplotype'},
        }
        self.assertEqual(results_expected, results)

    @unittest.skip("WIP")
    def test_unresolved_haplotype_because_mnv_covers_snv_starting_early(self) -> None:
        """Unresolved haplotype because MNV covers where SNV was expected, where MNV starts before this location"""
        panel = self.__get_narrow_example_panel({"*2A", "*5", "*2B"})
        ids_found_in_patient = Grch37CallData((
            Grch37Call(GeneCoordinate("1", 97915613), "GC", ("CT", "CT"), "DPYD", (".",), "9212GC>CT", "PASS"),
            Grch37Call(GeneCoordinate("1", 97981395), "T", ("C", "C"), "DPYD", ("rs1801159",), "293T>C", "PASS"),
            Grch37Call(GeneCoordinate("1", 97915621), "TG", ("TC", "TC"), "DPYD", ("rs72549303",), "6744GA>CA", "PASS"),
        ))
        results, panel_calls_for_patient = PgxAnalyser.create_pgx_analysis(ids_found_in_patient, panel)

        panel_calls_for_patient_expected = pd.DataFrame(
            [
                ("DPYD", "1:97915613", "CT", "CT", "UNKNOWN", "CT", "CT", ".", "9212GC>CT", "PASS"),
                ("DPYD", "1:97915621", "TC", "TC", "1:97450065", "TC", "TC", "rs72549303", "REF_CALL", "NO_CALL"),
                ("DPYD", "1:97981395", "C", "C", "1:97515839", "C", "C", "rs1801159", "293T>C", "PASS"),
            ], columns=ALL_IDS_IN_PANEL_COLUMNS
        )
        pd.testing.assert_frame_equal(panel_calls_for_patient_expected, panel_calls_for_patient)

        results_expected = {
            "DPYD": {'Unresolved_Haplotype'},
        }
        self.assertEqual(results_expected, results)
        self.fail("WIP")

    @unittest.skip("WIP")
    def test_unresolved_haplotype_because_mnv_covers_snv_starting_there(self) -> None:
        """Unresolved haplotype because MNV covers where SNV was expected, where MNV starts at this location"""
        panel = self.__get_narrow_example_panel({"*2A", "*5", "*2B"})
        ids_found_in_patient = Grch37CallData((
            Grch37Call(GeneCoordinate("1", 97915614), "C", ("TC", "TC"), "DPYD", (".",), "9212CG>TC", "PASS"),
            Grch37Call(GeneCoordinate("1", 97981395), "T", ("C", "C"), "DPYD", ("rs1801159",), "293T>C", "PASS"),
            Grch37Call(GeneCoordinate("1", 97915621), "TG", ("TC", "TC"), "DPYD", ("rs72549303",), "6744GA>CA", "PASS"),
        ))
        results, panel_calls_for_patient = PgxAnalyser.create_pgx_analysis(ids_found_in_patient, panel)

        panel_calls_for_patient_expected = pd.DataFrame(
            [
                ("DPYD", "1:97915614", "TC", "TC", "1:97450058", "TC", "TC", ".", "9212CG>TC", "PASS"),  # TODO: maybe GRCh38 should be unknown
                ("DPYD", "1:97915621", "TC", "TC", "1:97450065", "TC", "TC", "rs72549303", "REF_CALL", "NO_CALL"),
                ("DPYD", "1:97981395", "C", "C", "1:97515839", "C", "C", "rs1801159", "293T>C", "PASS"),
            ], columns=ALL_IDS_IN_PANEL_COLUMNS
        )
        pd.testing.assert_frame_equal(panel_calls_for_patient_expected, panel_calls_for_patient)

        results_expected = {
            "DPYD": {'Unresolved_Haplotype'},
        }
        self.assertEqual(results_expected, results)
        self.fail("WIP")

    @unittest.skip("WIP")
    def test_ambiguous_call(self) -> None:
        # TODO:
        #   Make sure that MNV's (mostly) work with variants output, but not with haplotype calling!!! (before SAGE GERMLINE!)
        #   Make sure ref seq differences can only be single alleles? should I? for now, probably yes.
        #   Bunch of errors
        #   SAGE GERMLINE as input
        #   (Maybe do "MNV instead of SNV" stuff after SAGE Germline stuff)
        #   MNV covering ref seq difference. Many possibilities for what should happen...
        #   MNV covering somewhere else. No haplotype?
        #   Currently filtering MNV from source file in or out? Which should it be?
        self.fail("WIP")


if __name__ == '__main__':
    unittest.main()
