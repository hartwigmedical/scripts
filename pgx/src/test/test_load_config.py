import unittest
from typing import Dict

from drug_info import DrugInfo
from gene_coordinate import GeneCoordinate
from gene_info import GeneInfo
from haplotype import Haplotype
from main import load_panel, load_ref_sequence_differences
from panel import Panel
from rs_id_info import RsIdInfo
from test_resources.test_resource import get_panel_test_resource, get_ref_sequence_differences_test_resource
from variant import Variant


class TestLoadConfig(unittest.TestCase):
    def test_load_panel(self) -> None:
        panel_path = get_panel_test_resource()
        panel = load_panel(str(panel_path))

        dpyd_two_a_variant = Variant("rs3918290", "T", "T")
        dpyd_two_b_variant = Variant("rs1801159", "C", "C")
        dpyd_three_variant = Variant("rs72549303", "T", "A")
        fake_variant = Variant("rs1212125", "C", "C")

        dpyd_haplotypes_expected = [
            Haplotype("*2A", "No Function", frozenset({dpyd_two_a_variant})),
            Haplotype("*2B", "No Function", frozenset({dpyd_two_a_variant, dpyd_two_b_variant})),
            Haplotype("*3", "Normal Function", frozenset({dpyd_three_variant})),
        ]
        dpyd_rs_id_infos_expected = frozenset({
            RsIdInfo("rs3918290", "C", "C", GeneCoordinate(1, 97915614), GeneCoordinate(1, 97450058)),
            RsIdInfo("rs72549309", "GATGA", "GATGA", GeneCoordinate(1, 98205966), GeneCoordinate(1, 97740410)),
            RsIdInfo("rs1801159", "T", "T", GeneCoordinate(1, 97981395), GeneCoordinate(1, 97515839)),
            RsIdInfo("rs72549303", "TG", "TC", GeneCoordinate(1, 97915621), GeneCoordinate(1, 97450065)),
        })
        dpyd_drugs_expected = [
            DrugInfo("5-Fluorouracil", "https://www.pharmgkb.org/chemical/PA128406956/guidelineAnnotation/PA166104939"),
            DrugInfo("Capecitabine", "https://www.pharmgkb.org/chemical/PA448771/guidelineAnnotation/PA166104963"),
        ]
        dpyd_rs_id_to_difference_annotations = {"rs72549303": "6744GA>CA"}
        fake_haplotypes_expected = [
            Haplotype("*4A", "Reduced Function", frozenset({fake_variant})),
        ]
        fake_rs_id_infos_expected = frozenset({
            RsIdInfo("rs1212125", "T", "T", GeneCoordinate(1, 97915617), GeneCoordinate(1, 97450060)),
        })
        fake_drugs_expected = [
            DrugInfo("Aspirin", "https://www.pharmgkb.org/some_other_url")
        ]
        fake_rs_id_to_difference_annotations: Dict[str, str] = {}

        gene_infos_expected = [
            GeneInfo("DPYD", "GRCh37", "*1", dpyd_haplotypes_expected, dpyd_rs_id_infos_expected, dpyd_drugs_expected, dpyd_rs_id_to_difference_annotations),
            GeneInfo("FAKE", "GRCh37", "*1", fake_haplotypes_expected, fake_rs_id_infos_expected, fake_drugs_expected, fake_rs_id_to_difference_annotations),
        ]
        panel_expected = Panel(gene_infos_expected)

        self.assertEqual(panel_expected, panel)

    def test_load_ref_sequence_differences(self) -> None:
        ref_seq_diff_path = get_ref_sequence_differences_test_resource()
        ref_seq_differences = load_ref_sequence_differences(str(ref_seq_diff_path))

        ref_seq_differences_expected = [
            {
                'rsid': 'rs1801265',
                'gene': 'DPYD',
                'referenceAlleleGRCh38': 'A',
                'altAlleleGRCh38': 'G',
                'chromosome': 1,
                'position': '98348885',
                'positionGRCh38': '97883329',
                'annotationGRCh38': '85T>C'
            },
            {
                'rsid': 'rs1801270',
                'gene': 'FAKE',
                'referenceAlleleGRCh38': 'T',
                'altAlleleGRCh38': 'C',
                'chromosome': 5,
                'position': '98348890',
                'positionGRCh38': '97883345',
                'annotationGRCh38': '236T>C'
            }
        ]

        self.assertEqual(ref_seq_differences_expected, ref_seq_differences)


if __name__ == '__main__':
    unittest.main()
