import unittest
from typing import Dict

import pandas as pd

from base.gene_coordinate import GeneCoordinate
from call_data import Grch37CallData
from config.drug_info import DrugInfo
from config.gene_info import GeneInfo
from config.haplotype import Haplotype
from config.panel import Panel
from config.rs_id_info import RsIdInfo
from config.variant import Variant
from pgx_analysis import create_pgx_analysis


class TestPgxAnalysis(unittest.TestCase):
    @classmethod
    def __get_example_panel(cls) -> Panel:
        dpyd_two_a_variant = Variant("rs3918290", "T")
        dpyd_two_b_variant = Variant("rs1801159", "C")
        dpyd_three_variant = Variant("rs72549303", "TG")
        fake_variant = Variant("rs1212125", "C")
        fake2_variant = Variant("rs1212127", "C")

        dpyd_haplotypes = (
            Haplotype("*2A", "No Function", frozenset({dpyd_two_a_variant})),
            Haplotype("*2B", "No Function", frozenset({dpyd_two_a_variant, dpyd_two_b_variant})),
            Haplotype("*3", "Normal Function", frozenset({dpyd_three_variant})),
        )
        dpyd_rs_id_infos = frozenset({
            RsIdInfo("rs3918290", "C", "C", GeneCoordinate("1", 97915614), GeneCoordinate("1", 97450058)),
            RsIdInfo("rs72549309", "GATGA", "GATGA", GeneCoordinate("1", 98205966), GeneCoordinate("1", 97740410)),
            RsIdInfo("rs1801159", "T", "T", GeneCoordinate("1", 97981395), GeneCoordinate("1", 97515839)),
            RsIdInfo("rs72549303", "TG", "TC", GeneCoordinate("1", 97915621), GeneCoordinate("1", 97450065)),
        })
        dpyd_drugs = (
            DrugInfo("5-Fluorouracil", "https://www.pharmgkb.org/chemical/PA128406956/guidelineAnnotation/PA166104939"),
            DrugInfo("Capecitabine", "https://www.pharmgkb.org/chemical/PA448771/guidelineAnnotation/PA166104963"),
        )
        dpyd_rs_id_to_difference_annotations = {
            "rs72549303": "6744GA>CA",
        }
        
        fake_haplotypes = (
            Haplotype("*4A", "Reduced Function", frozenset({fake_variant})),
        )
        fake_rs_id_infos = frozenset({
            RsIdInfo("rs1212125", "T", "T", GeneCoordinate("5", 97915617), GeneCoordinate("5", 97450060)),
        })
        fake_drugs = (
            DrugInfo("Aspirin", "https://www.pharmgkb.org/some_other_url"),
        )
        fake_rs_id_to_difference_annotations: Dict[str, str] = {}

        fake2_haplotypes = (
            Haplotype("*4A", "Reduced Function", frozenset({fake2_variant})),
        )
        fake2_rs_id_infos = frozenset({
            RsIdInfo("rs1212127", "C", "T", GeneCoordinate("16", 97915617), GeneCoordinate("16", 97450060)),
        })
        fake2_drugs = (
            DrugInfo("Aspirin", "https://www.pharmgkb.org/some_other_url"),
        )
        fake2_rs_id_to_difference_annotations: Dict[str, str] = {"rs1212127": "1324T>C"}

        gene_infos = (
            GeneInfo("DPYD", "1", "GRCh37", "*1", dpyd_haplotypes, dpyd_rs_id_infos,
                     dpyd_drugs, dpyd_rs_id_to_difference_annotations),
            GeneInfo("FAKE", "5", "GRCh37", "*1", fake_haplotypes, fake_rs_id_infos,
                     fake_drugs, fake_rs_id_to_difference_annotations),
            GeneInfo("FAKE2", "16", "GRCh37", "*1", fake2_haplotypes, fake2_rs_id_infos,
                     fake2_drugs, fake2_rs_id_to_difference_annotations),
        )
        return Panel(gene_infos)

    def test_empty(self) -> None:
        panel = self.__get_example_panel()
        ids_found_in_patient = Grch37CallData(tuple())
        results, severity, all_ids_in_panel, drug_info = create_pgx_analysis(ids_found_in_patient, panel)

        results_expected = {"DPYD": ["*3_HOM"], "FAKE": ["*1_HOM"], "FAKE2": ["*4A_HOM"]}
        self.assertEqual(results_expected, results)

        severity_expected = {
            '*1': 'Normal Function',
            'Unresolved': 'Unknown Function',
            '*2A': 'No Function',
            '*2B': 'No Function',
            '*3': 'Normal Function',
            '*4A': 'Reduced Function'
        }
        self.assertEqual(severity_expected, severity)

        columns_expected = [
            'gene', 'position_GRCh37', 'ref_GRCh37', 'alt_GRCh37', 'position_GRCh38', 'ref_GRCh38', 'alt_GRCh38',
            'rsid', 'variant_annotation', 'filter'
        ]
        all_ids_in_panel_expected = pd.DataFrame(
            [
                ("FAKE2", "16:97915617", "C", "C", "16:97450060", "C", "C", "rs1212127", "1324T>C", "INFERRED_REF_CALL"),
                ("DPYD", "1:97915614", "C", "C", "1:97450058", "C", "C", "rs3918290", "REF_CALL", "NO_CALL"),
                ("DPYD", "1:97915621", "TG", "TG", "1:97450065", "TG", "TG", "rs72549303", "6744GA>CA", "INFERRED_REF_CALL"),
                ("DPYD", "1:97981395", "T", "T", "1:97515839", "T", "T", "rs1801159", "REF_CALL", "NO_CALL"),
                ("DPYD", "1:98205966", "GATGA", "GATGA", "1:97740410", "GATGA", "GATGA", "rs72549309", "REF_CALL", "NO_CALL"),
                ("FAKE", "5:97915617", "T", "T", "5:97450060", "T", "T", "rs1212125", "REF_CALL", "NO_CALL"),
            ], columns=columns_expected
        )
        pd.testing.assert_frame_equal(all_ids_in_panel_expected, all_ids_in_panel)

        drug_info_expected = {
            'DPYD': [
                '5-Fluorouracil;Capecitabine',
                'https://www.pharmgkb.org/chemical/PA128406956/guidelineAnnotation/PA166104939;https://www.pharmgkb.org/chemical/PA448771/guidelineAnnotation/PA166104963'
            ],
            'FAKE': [
                'Aspirin',
                'https://www.pharmgkb.org/some_other_url'
            ],
            'FAKE2': [
                'Aspirin',
                'https://www.pharmgkb.org/some_other_url'
            ]
        }
        self.assertEqual(drug_info_expected, drug_info)


if __name__ == '__main__':
    unittest.main()
