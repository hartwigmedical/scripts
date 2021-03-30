import unittest
from typing import Dict

from base.gene_coordinate import GeneCoordinate
from config.drug_info import DrugInfo
from config.gene_info import GeneInfo
from config.haplotype import Haplotype
from config.panel import Panel
from config.rs_id_info import RsIdInfo
from config.variant import Variant
from main import load_panel
from test_resources.test_resource import get_panel_test_resource


class TestLoadConfig(unittest.TestCase):
    def test_load_panel(self) -> None:
        panel_path = get_panel_test_resource()
        panel = load_panel(str(panel_path))

        dpyd_two_a_variant = Variant("rs3918290", "T")
        dpyd_two_b_variant = Variant("rs1801159", "C")
        dpyd_three_variant = Variant("rs72549303", "TG")
        fake_variant = Variant("rs1212125", "C")
        fake2_variant = Variant("rs1212127", "C")

        dpyd_haplotypes_expected = frozenset({
            Haplotype("*2A", "No Function", frozenset({dpyd_two_a_variant})),
            Haplotype("*2B", "No Function", frozenset({dpyd_two_a_variant, dpyd_two_b_variant})),
            Haplotype("*3", "Normal Function", frozenset({dpyd_three_variant})),
        })
        dpyd_rs_id_infos_expected = frozenset({
            RsIdInfo("rs3918290", "C", "C", GeneCoordinate("1", 97915614), GeneCoordinate("1", 97450058)),
            RsIdInfo("rs72549309", "GATGA", "GATGA", GeneCoordinate("1", 98205966), GeneCoordinate("1", 97740410)),
            RsIdInfo("rs1801159", "T", "T", GeneCoordinate("1", 97981395), GeneCoordinate("1", 97515839)),
            RsIdInfo("rs72549303", "TG", "TC", GeneCoordinate("1", 97915621), GeneCoordinate("1", 97450065)),
            RsIdInfo("rs1801265", "G", "A", GeneCoordinate("1", 98348885), GeneCoordinate("1", 97883329)),
        })
        dpyd_drugs_expected = frozenset({
            DrugInfo("5-Fluorouracil", "https://www.pharmgkb.org/chemical/PA128406956/guidelineAnnotation/PA166104939"),
            DrugInfo("Capecitabine", "https://www.pharmgkb.org/chemical/PA448771/guidelineAnnotation/PA166104963"),
        })
        dpyd_rs_id_to_difference_annotations = {
            "rs72549303": "6744GA>CA",
            "rs1801265": "85T>C",
        }
        fake_haplotypes_expected = frozenset({
            Haplotype("*4A", "Reduced Function", frozenset({fake_variant})),
        })
        fake_rs_id_infos_expected = frozenset({
            RsIdInfo("rs1212125", "T", "T", GeneCoordinate("5", 97915617), GeneCoordinate("5", 97450060)),
        })
        fake_drugs_expected = frozenset({
            DrugInfo("Aspirin", "https://www.pharmgkb.org/some_other_url"),
        })
        fake_rs_id_to_difference_annotations: Dict[str, str] = {}

        fake2_haplotypes_expected = frozenset({
            Haplotype("*4A", "Reduced Function", frozenset({fake2_variant})),
        })
        fake2_rs_id_infos_expected = frozenset({
            RsIdInfo("rs1212127", "C", "T", GeneCoordinate("16", 97915617), GeneCoordinate("16", 97450060)),
        })
        fake2_drugs_expected = frozenset({
            DrugInfo("Aspirin", "https://www.pharmgkb.org/some_other_url"),
        })
        fake2_rs_id_to_difference_annotations: Dict[str, str] = {"rs1212127": "1324T>C"}

        gene_infos_expected = frozenset({
            GeneInfo("DPYD", "1", "*1", dpyd_haplotypes_expected, dpyd_rs_id_infos_expected,
                     dpyd_drugs_expected, dpyd_rs_id_to_difference_annotations),
            GeneInfo("FAKE", "5", "*1", fake_haplotypes_expected, fake_rs_id_infos_expected,
                     fake_drugs_expected, fake_rs_id_to_difference_annotations),
            GeneInfo("FAKE2", "16", "*1", fake2_haplotypes_expected, fake2_rs_id_infos_expected,
                     fake2_drugs_expected, fake2_rs_id_to_difference_annotations),
        })
        name_expected = "fake_panel"
        version_expected = "0.2"
        panel_expected = Panel(name_expected, version_expected, gene_infos_expected)

        self.assertEqual(panel_expected, panel)

    def test_load_ref_sequence_differences_from_panel(self) -> None:
        panel_path = get_panel_test_resource()
        panel = load_panel(str(panel_path))

        ref_seq_differences = panel.get_ref_seq_differences()

        ref_seq_differences_expected = [
            (
                RsIdInfo('rs72549303', 'TG', 'TC', GeneCoordinate("1", 97915621), GeneCoordinate("1", 97450065)),
                'DPYD',
                '6744GA>CA'
            ),
            (
                RsIdInfo('rs1801265', 'G', 'A', GeneCoordinate("1", 98348885), GeneCoordinate("1", 97883329)),
                'DPYD',
                '85T>C'
            ),
            (
                RsIdInfo('rs1212127', 'C', 'T', GeneCoordinate("16", 97915617), GeneCoordinate("16", 97450060)),
                'FAKE2',
                '1324T>C'
            ),
        ]

        self.assertEqual(ref_seq_differences_expected, ref_seq_differences)


if __name__ == '__main__':
    unittest.main()
