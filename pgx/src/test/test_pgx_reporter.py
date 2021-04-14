import unittest
from typing import Dict, Set

from base.filter import Filter
from base.gene_coordinate import GeneCoordinate
from call_data import FullCall, HaplotypeCall
from config.drug_info import DrugInfo
from config.gene_info import GeneInfo
from config.haplotype import Haplotype
from config.panel import Panel
from config.rs_id_info import RsIdInfo
from config.variant import Variant
from pgx_analysis import PgxAnalysis
from pgx_reporter import GenotypeReporter, HaplotypeReporter


class TestPgxReporter(unittest.TestCase):
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
        version = "1.1"
        return Panel(name, version, gene_infos)

    def test_genotype_reporter_empty(self) -> None:
        pgx_analysis = PgxAnalysis(frozenset(), {})
        panel_id = "Panel_v0.2"
        version = "V1"
        result = GenotypeReporter.get_calls_tsv_text(pgx_analysis, panel_id, version)

        result_expected = (
            "gene\tchromosome\tposition_v37\tposition_v38\tref_v37\tref_v38\t"
            "allele1\tallele2\trsid\tvariant_annotation_v37\tfilter_v37\tvariant_annotation_v38\tfilter_v38\tpanel_version\trepo_version\n"
        )
        self.assertEqual(result_expected, result)

    def test_genotype_reporter_non_empty(self) -> None:
        all_full_calls = frozenset({
            FullCall(GeneCoordinate("1", 5), "A", GeneCoordinate("1", 25), "A", ("G", "C"),
                     "DPYD", (".",), "25A>C;25A>G", Filter.PASS, "25A>C;25A>G", Filter.PASS),
            FullCall(GeneCoordinate("1", 15), "C", None, None, ("C", "CAG"),
                     "DPYD", ("rs536",), "35A>C;35A>G", Filter.PASS, "35A>C;35A>G?", Filter.UNKNOWN),
            FullCall(GeneCoordinate("X", 15), "TT", GeneCoordinate("X", 40), "AA", ("TT", "TT"),
                     "GENE", ("rs23",), "REF_CALL", Filter.NO_CALL, "627AA>TT", Filter.PASS),
            FullCall(GeneCoordinate("2", 154663), "T", GeneCoordinate("2", 40565464), "T", ("T", "T"),
                     "BRAF", ("rs154", "rs8839"), "REF_CALL", Filter.NO_CALL, "REF_CALL", Filter.NO_CALL),
            FullCall(GeneCoordinate("15", 24113), "A", GeneCoordinate("15", 684633), "T", ("T", "T"),
                     ".", ("rs462", "rs9820", "rs536"), "29482A>T", Filter.PASS, "REF_CALL", Filter.NO_CALL),
        })
        pgx_analysis = PgxAnalysis(all_full_calls, {})
        panel_id = "Panel_v0.2"
        version = "V1"
        result = GenotypeReporter.get_calls_tsv_text(pgx_analysis, panel_id, version)

        result_expected = (
            "gene\tchromosome\tposition_v37\tposition_v38\tref_v37\tref_v38\tallele1\tallele2\trsid\tvariant_annotation_v37\tfilter_v37\tvariant_annotation_v38\tfilter_v38\tpanel_version\trepo_version\n"
            "DPYD\t1\t5\t25\tA\tA\tC\tG\t.\t25A>C;25A>G\tPASS\t25A>C;25A>G\tPASS\tPanel_v0.2\tV1\n"
            "DPYD\t1\t15\tUNKNOWN\tC\tUNKNOWN\tC\tCAG\trs536\t35A>C;35A>G\tPASS\t35A>C;35A>G?\tUNKNOWN\tPanel_v0.2\tV1\n"
            "BRAF\t2\t154663\t40565464\tT\tT\tT\tT\trs154;rs8839\tREF_CALL\tNO_CALL\tREF_CALL\tNO_CALL\tPanel_v0.2\tV1\n"
            ".\t15\t24113\t684633\tA\tT\tT\tT\trs462;rs9820;rs536\t29482A>T\tPASS\tREF_CALL\tNO_CALL\tPanel_v0.2\tV1\n"
            "GENE\tX\t15\t40\tTT\tAA\tTT\tTT\trs23\tREF_CALL\tNO_CALL\t627AA>TT\tPASS\tPanel_v0.2\tV1\n"
        )
        self.assertEqual(result_expected, result)

    def test_haplotype_reporter_empty(self) -> None:
        pgx_analysis = PgxAnalysis(frozenset(), {})
        panel = Panel("EmptyPanel", "0.3", frozenset())
        panel_id = "Panel_v0.2"
        version = "V1"
        result = HaplotypeReporter.get_genotype_tsv_text(pgx_analysis, panel, panel_id, version)

        result_expected = "gene\thaplotype\tfunction\tlinked_drugs\turl_prescription_info\tpanel_version\trepo_version\n"
        self.assertEqual(result_expected, result)

    def test_haplotype_reporter_non_empty(self) -> None:
        empty_haplotype_set: Set[HaplotypeCall] = set()
        gene_to_haplotype_calls = {
            "DPYD": {HaplotypeCall("*2B", 2), HaplotypeCall("*3", 1), HaplotypeCall("*2A", 1)},
            "FAKE": empty_haplotype_set,
            "FAKE2": {HaplotypeCall("*1", 1), HaplotypeCall("*4A", 1)},
        }
        pgx_analysis = PgxAnalysis(frozenset(), gene_to_haplotype_calls)
        panel = self.__get_wide_example_panel()
        panel_id = "Panel_v0.2"
        version = "V1"
        result = HaplotypeReporter.get_genotype_tsv_text(pgx_analysis, panel, panel_id, version)

        result_expected = (
            "gene\thaplotype\tfunction\tlinked_drugs\turl_prescription_info\tpanel_version\trepo_version\n"
            "DPYD\t*2A_HET\tNo Function\t5-Fluorouracil;Capecitabine\thttps://www.pharmgkb.org/chemical/PA128406956/guidelineAnnotation/PA166104939;https://www.pharmgkb.org/chemical/PA448771/guidelineAnnotation/PA166104963\tPanel_v0.2\tV1\n"
            "DPYD\t*2B_HOM\tNo Function\t5-Fluorouracil;Capecitabine\thttps://www.pharmgkb.org/chemical/PA128406956/guidelineAnnotation/PA166104939;https://www.pharmgkb.org/chemical/PA448771/guidelineAnnotation/PA166104963\tPanel_v0.2\tV1\n"
            "DPYD\t*3_HET\tNormal Function\t5-Fluorouracil;Capecitabine\thttps://www.pharmgkb.org/chemical/PA128406956/guidelineAnnotation/PA166104939;https://www.pharmgkb.org/chemical/PA448771/guidelineAnnotation/PA166104963\tPanel_v0.2\tV1\n"
            "FAKE\tUnresolved Haplotype\tUnknown Function\tAspirin\thttps://www.pharmgkb.org/some_other_url\tPanel_v0.2\tV1\n"
            "FAKE2\t*1_HET\tNormal Function\tAspirin\thttps://www.pharmgkb.org/some_other_url\tPanel_v0.2\tV1\n"
            "FAKE2\t*4A_HET\tReduced Function\tAspirin\thttps://www.pharmgkb.org/some_other_url\tPanel_v0.2\tV1\n"
        )
        self.assertEqual(result_expected, result)


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
