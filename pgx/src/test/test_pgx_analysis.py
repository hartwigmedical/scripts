import unittest
from typing import Dict

import pandas as pd

from base.gene_coordinate import GeneCoordinate
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
        dpyd_two_a_variant = Variant("rs3918290", "T", "T")
        dpyd_two_b_variant = Variant("rs1801159", "C", "C")
        dpyd_three_variant = Variant("rs72549303", "T", "A")
        fake_variant = Variant("rs1212125", "C", "C")
        fake2_variant = Variant("rs1212127", "T", "C")

        dpyd_haplotypes = [
            Haplotype("*2A", "No Function", frozenset({dpyd_two_a_variant})),
            Haplotype("*2B", "No Function", frozenset({dpyd_two_a_variant, dpyd_two_b_variant})),
            Haplotype("*3", "Normal Function", frozenset({dpyd_three_variant})),
        ]
        dpyd_rs_id_infos = frozenset({
            RsIdInfo("rs3918290", "C", "C", GeneCoordinate("1", 97915614), GeneCoordinate("1", 97450058)),
            RsIdInfo("rs72549309", "GATGA", "GATGA", GeneCoordinate("1", 98205966), GeneCoordinate("1", 97740410)),
            RsIdInfo("rs1801159", "T", "T", GeneCoordinate("1", 97981395), GeneCoordinate("1", 97515839)),
            RsIdInfo("rs72549303", "TG", "TC", GeneCoordinate("1", 97915621), GeneCoordinate("1", 97450065)),
            RsIdInfo("rs1801265", "G", "A", GeneCoordinate("1", 98348885), GeneCoordinate("1", 97883329)),
        })
        dpyd_drugs = [
            DrugInfo("5-Fluorouracil", "https://www.pharmgkb.org/chemical/PA128406956/guidelineAnnotation/PA166104939"),
            DrugInfo("Capecitabine", "https://www.pharmgkb.org/chemical/PA448771/guidelineAnnotation/PA166104963"),
        ]
        dpyd_rs_id_to_difference_annotations = {
            "rs72549303": "6744GA>CA",
            "rs1801265": "85T>C",
        }
        fake_haplotypes = [
            Haplotype("*4A", "Reduced Function", frozenset({fake_variant})),
        ]
        fake_rs_id_infos = frozenset({
            RsIdInfo("rs1212125", "T", "T", GeneCoordinate("5", 97915617), GeneCoordinate("5", 97450060)),
        })
        fake_drugs = [
            DrugInfo("Aspirin", "https://www.pharmgkb.org/some_other_url")
        ]
        fake_rs_id_to_difference_annotations: Dict[str, str] = {}

        fake2_haplotypes = [
            Haplotype("*4A", "Reduced Function", frozenset({fake2_variant})),
        ]
        fake2_rs_id_infos = frozenset({
            RsIdInfo("rs1212127", "C", "T", GeneCoordinate("16", 97915617), GeneCoordinate("16", 97450060)),
        })
        fake2_drugs = [
            DrugInfo("Aspirin", "https://www.pharmgkb.org/some_other_url")
        ]
        fake2_rs_id_to_difference_annotations: Dict[str, str] = {"rs1212127": "1324T>C"}

        gene_infos = [
            GeneInfo("DPYD", "1", "GRCh37", "*1", dpyd_haplotypes, dpyd_rs_id_infos,
                     dpyd_drugs, dpyd_rs_id_to_difference_annotations),
            GeneInfo("FAKE", "5", "GRCh37", "*1", fake_haplotypes, fake_rs_id_infos,
                     fake_drugs, fake_rs_id_to_difference_annotations),
            GeneInfo("FAKE2", "16", "GRCh37", "*1", fake2_haplotypes, fake2_rs_id_infos,
                     fake2_drugs, fake2_rs_id_to_difference_annotations),
        ]
        return Panel(gene_infos)

    @unittest.skip("WIP")
    def test_empty(self) -> None:
        panel = self.__get_example_panel()
        ids_found_in_patient = pd.DataFrame()
        results, severity, all_ids_in_panel, drug_info = create_pgx_analysis(ids_found_in_patient, panel)
        print(results)
        print(severity)
        print(all_ids_in_panel)
        print(drug_info)
        self.fail("WIP")


if __name__ == '__main__':
    unittest.main()
