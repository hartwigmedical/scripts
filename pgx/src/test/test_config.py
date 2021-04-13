import unittest
from typing import Dict, FrozenSet

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
        """Load panel from json"""
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
            DrugInfo("5-Fluorouracil", "https://www.source_url.org/5-Fluorouracil"),
            DrugInfo("Capecitabine", "https://www.source_url.org/Capecitabine"),
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
            DrugInfo("Aspirin", "https://www.source_url.org/Aspirin"),
        })
        fake_rs_id_to_difference_annotations: Dict[str, str] = {}

        fake2_haplotypes_expected = frozenset({
            Haplotype("*4A", "Reduced Function", frozenset({fake2_variant})),
        })
        fake2_rs_id_infos_expected = frozenset({
            RsIdInfo("rs1212127", "C", "T", GeneCoordinate("16", 97915617), GeneCoordinate("16", 97450060)),
        })
        fake2_drugs_expected = frozenset({
            DrugInfo("Aspirin", "https://www.source_url.org/Aspirin"),
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

    def test_panel_name_from_json(self) -> None:
        panel_path = get_panel_test_resource()
        panel = load_panel(str(panel_path))
        expected_panel_id = "fake_panel_v0.2"
        self.assertEqual(expected_panel_id, panel.get_id())

    def test_gene_info_with_overlapping_haplotype_names(self) -> None:
        """Error when haplotype name used multiple times for gene"""
        gene = "FAKE"
        chromosome = "X"
        reference_haplotype_name = "*1"
        drugs: FrozenSet[DrugInfo] = frozenset()
        rs_id_to_ref_seq_difference_annotation: Dict[str, str] = dict()

        variant1 = Variant("rs94982", "A")
        variant2 = Variant("rs394934", "T")

        rs_id_infos = frozenset([
            RsIdInfo(variant1.rs_id, "C", "C",
                     GeneCoordinate(chromosome, 4994545), GeneCoordinate(chromosome, 2993823)),
            RsIdInfo(variant2.rs_id, "G", "G",
                     GeneCoordinate(chromosome, 3993842), GeneCoordinate(chromosome, 2949923)),
        ])

        haplotypes1 = frozenset([
            Haplotype("*2", "No Function", frozenset([variant1])),
            Haplotype("*4", "Partial Function", frozenset([variant1, variant2])),
        ])
        haplotypes2 = frozenset([
            Haplotype("*2", "No Function", frozenset([variant1])),
            Haplotype("*2", "Partial Function", frozenset([variant1, variant2])),
        ])

        GeneInfo(gene, chromosome, reference_haplotype_name, haplotypes1, rs_id_infos,
                 drugs, rs_id_to_ref_seq_difference_annotation)
        with self.assertRaises(ValueError):
            GeneInfo(gene, chromosome, reference_haplotype_name, haplotypes2, rs_id_infos,
                     drugs, rs_id_to_ref_seq_difference_annotation)

    def test_gene_info_with_overlapping_haplotype_variants(self) -> None:
        """Error when different haplotypes have the exact same variant combination"""
        gene = "FAKE"
        chromosome = "X"
        reference_haplotype_name = "*1"
        drugs: FrozenSet[DrugInfo] = frozenset()
        rs_id_to_ref_seq_difference_annotation: Dict[str, str] = dict()

        variant1 = Variant("rs94982", "A")
        variant2 = Variant("rs394934", "T")
        variant3 = Variant("rs495825", "C")

        rs_id_infos = frozenset([
            RsIdInfo(variant1.rs_id, "C", "C",
                     GeneCoordinate(chromosome, 4994545), GeneCoordinate(chromosome, 2993823)),
            RsIdInfo(variant2.rs_id, "G", "G",
                     GeneCoordinate(chromosome, 3993842), GeneCoordinate(chromosome, 2949923)),
            RsIdInfo(variant3.rs_id, "A", "A",
                     GeneCoordinate(chromosome, 293923), GeneCoordinate(chromosome, 138812)),
        ])

        haplotypes1 = frozenset([
            Haplotype("*2", "No Function", frozenset([variant1, variant2])),
            Haplotype("*4", "Partial Function", frozenset([variant1, variant3])),
        ])
        haplotypes2 = frozenset([
            Haplotype("*2", "No Function", frozenset([variant1, variant2])),
            Haplotype("*4", "Partial Function", frozenset([variant1, variant2])),
        ])

        GeneInfo(gene, chromosome, reference_haplotype_name, haplotypes1, rs_id_infos,
                 drugs, rs_id_to_ref_seq_difference_annotation)
        with self.assertRaises(ValueError):
            GeneInfo(gene, chromosome, reference_haplotype_name, haplotypes2, rs_id_infos,
                     drugs, rs_id_to_ref_seq_difference_annotation)

    def test_gene_info_with_overlapping_drug_names(self) -> None:
        """Error when drug name used multiple times for gene"""
        gene = "FAKE"
        chromosome = "X"
        reference_haplotype_name = "*1"
        haplotypes: FrozenSet[Haplotype] = frozenset()
        rs_id_infos: FrozenSet[RsIdInfo] = frozenset()
        rs_id_to_ref_seq_difference_annotation: Dict[str, str] = dict()

        druginfo1 = DrugInfo("Paracetamol", "google.com")
        druginfo2 = DrugInfo("Paracetamol", "wikipedia.com")
        drugs = frozenset([druginfo1])
        overlapping_drugs = frozenset([druginfo1, druginfo2])

        GeneInfo(gene, chromosome, reference_haplotype_name, haplotypes, rs_id_infos,
                 drugs, rs_id_to_ref_seq_difference_annotation)
        with self.assertRaises(ValueError):
            GeneInfo(gene, chromosome, reference_haplotype_name, haplotypes, rs_id_infos,
                     overlapping_drugs, rs_id_to_ref_seq_difference_annotation)

    def test_gene_info_with_overlapping_rs_id_infos(self) -> None:
        """Error when gene info has rs id infos for which the relevant coordinates overlap"""
        gene = "FAKE"
        chromosome = "X"
        reference_haplotype_name = "*1"
        haplotypes: FrozenSet[Haplotype] = frozenset()
        drugs: FrozenSet[DrugInfo] = frozenset()
        rs_id_to_ref_seq_difference_annotation: Dict[str, str] = dict()

        rs_id_info1 = RsIdInfo("rs294924", "A", "A",
                               GeneCoordinate(chromosome, 499593), GeneCoordinate(chromosome, 399483))
        rs_id_info2 = RsIdInfo("rs294927", "AA", "AA",
                               GeneCoordinate(chromosome, 499592), GeneCoordinate(chromosome, 399482))
        single_rs_id_info = frozenset([rs_id_info1])
        overlapping_rs_id_infos = frozenset([rs_id_info1, rs_id_info2])

        GeneInfo(gene, chromosome, reference_haplotype_name, haplotypes, single_rs_id_info,
                 drugs, rs_id_to_ref_seq_difference_annotation)
        with self.assertRaises(ValueError):
            GeneInfo(gene, chromosome, reference_haplotype_name, haplotypes, overlapping_rs_id_infos,
                     drugs, rs_id_to_ref_seq_difference_annotation)

    def test_gene_info_with_rs_id_infos_for_different_chromosome(self) -> None:
        """Error when gene info has rs id infos for which the relevant coordinates have a different chromosome"""
        gene = "FAKE"
        chromosome = "X"
        reference_haplotype_name = "*1"
        haplotypes: FrozenSet[Haplotype] = frozenset()
        drugs: FrozenSet[DrugInfo] = frozenset()
        rs_id_to_ref_seq_difference_annotation: Dict[str, str] = dict()

        other_chromosome = "1"

        rs_id_info1 = RsIdInfo("rs294924", "A", "A",
                               GeneCoordinate(chromosome, 499593), GeneCoordinate(chromosome, 399483))
        rs_id_info2 = RsIdInfo("rs294924", "A", "A",
                               GeneCoordinate(other_chromosome, 499593), GeneCoordinate(chromosome, 399483))
        rs_id_info3 = RsIdInfo("rs294924", "A", "A",
                               GeneCoordinate(chromosome, 499593), GeneCoordinate(other_chromosome, 399483))
        rs_id_info4 = RsIdInfo("rs294924", "A", "A",
                               GeneCoordinate(other_chromosome, 499593), GeneCoordinate(other_chromosome, 399483))

        GeneInfo(gene, chromosome, reference_haplotype_name, haplotypes, frozenset([rs_id_info1]),
                 drugs, rs_id_to_ref_seq_difference_annotation)
        with self.assertRaises(ValueError):
            GeneInfo(gene, chromosome, reference_haplotype_name, haplotypes, frozenset([rs_id_info2]),
                     drugs, rs_id_to_ref_seq_difference_annotation)
        with self.assertRaises(ValueError):
            GeneInfo(gene, chromosome, reference_haplotype_name, haplotypes, frozenset([rs_id_info3]),
                     drugs, rs_id_to_ref_seq_difference_annotation)
        with self.assertRaises(ValueError):
            GeneInfo(gene, chromosome, reference_haplotype_name, haplotypes, frozenset([rs_id_info4]),
                     drugs, rs_id_to_ref_seq_difference_annotation)

    def test_gene_info_with_rs_id_in_haplotype_without_info(self) -> None:
        """Error when haplotype uses rs id for which there is no RsIdInfo object"""
        gene = "FAKE"
        chromosome = "X"
        reference_haplotype_name = "*1"
        rs_id_infos: FrozenSet[RsIdInfo] = frozenset()
        drugs: FrozenSet[DrugInfo] = frozenset()
        rs_id_to_ref_seq_difference_annotation: Dict[str, str] = dict()

        empty_haplotypes: FrozenSet[Haplotype] = frozenset()
        non_empty_haplotypes = frozenset([
            Haplotype("*2", "No Function", frozenset([Variant("rs238423", "A")]))
        ])

        GeneInfo(gene, chromosome, reference_haplotype_name, empty_haplotypes, rs_id_infos,
                 drugs, rs_id_to_ref_seq_difference_annotation)
        with self.assertRaises(ValueError):
            GeneInfo(gene, chromosome, reference_haplotype_name, non_empty_haplotypes, rs_id_infos,
                     drugs, rs_id_to_ref_seq_difference_annotation)

    def test_gene_info_with_variant_that_is_ref_v38(self) -> None:
        """Error when gene info has variant where variant allele is the ref allele"""
        gene = "FAKE"
        chromosome = "X"
        reference_haplotype_name = "*1"
        drugs: FrozenSet[DrugInfo] = frozenset()

        rs_id = "rs294924"

        rs_id_to_ref_seq_difference_annotation = {rs_id: "399483C>A"}
        haplotypes = frozenset([Haplotype("*3", "No Function", frozenset([Variant("rs294924", "G")]))])

        rs_id_infos1 = frozenset([RsIdInfo(rs_id, "A", "C",
                                           GeneCoordinate(chromosome, 499593), GeneCoordinate(chromosome, 399483))])
        rs_id_infos2 = frozenset([RsIdInfo(rs_id, "A", "G",
                                           GeneCoordinate(chromosome, 499593), GeneCoordinate(chromosome, 399483))])

        GeneInfo(gene, chromosome, reference_haplotype_name, haplotypes, rs_id_infos1,
                 drugs, rs_id_to_ref_seq_difference_annotation)
        with self.assertRaises(ValueError):
            GeneInfo(gene, chromosome, reference_haplotype_name, haplotypes, rs_id_infos2,
                     drugs, rs_id_to_ref_seq_difference_annotation)

    def test_gene_info_with_ref_seq_difference_without_annotation(self) -> None:
        """Error when a ref seq difference does not have an annotation"""
        gene = "FAKE"
        chromosome = "X"
        reference_haplotype_name = "*1"
        haplotypes: FrozenSet[Haplotype] = frozenset()
        drugs: FrozenSet[DrugInfo] = frozenset()
        rs_id_to_ref_seq_difference_annotation: Dict[str, str] = dict()

        empty_rs_id_infos: FrozenSet[RsIdInfo] = frozenset()
        non_empty_rs_id_infos = frozenset([
            RsIdInfo("rs294924", "A", "C", GeneCoordinate(chromosome, 499593), GeneCoordinate(chromosome, 399483))
        ])

        GeneInfo(gene, chromosome, reference_haplotype_name, haplotypes, empty_rs_id_infos,
                 drugs, rs_id_to_ref_seq_difference_annotation)
        with self.assertRaises(ValueError):
            GeneInfo(gene, chromosome, reference_haplotype_name, haplotypes, non_empty_rs_id_infos,
                     drugs, rs_id_to_ref_seq_difference_annotation)

    def test_gene_info_with_extra_ref_seq_difference_annotation(self) -> None:
        """Error when a ref seq difference annotation does not match known a known ref seq difference"""
        gene = "FAKE"
        chromosome = "X"
        reference_haplotype_name = "*1"
        haplotypes: FrozenSet[Haplotype] = frozenset()
        rs_id_infos: FrozenSet[RsIdInfo] = frozenset()
        drugs: FrozenSet[DrugInfo] = frozenset()

        empty_rs_id_to_ref_seq_difference_annotation: Dict[str, str] = dict()
        non_empty_rs_id_to_ref_seq_difference_annotation = {"rs2493023": "3445A>T"}

        GeneInfo(gene, chromosome, reference_haplotype_name, haplotypes, rs_id_infos,
                 drugs, empty_rs_id_to_ref_seq_difference_annotation)
        with self.assertRaises(ValueError):
            GeneInfo(gene, chromosome, reference_haplotype_name, haplotypes, rs_id_infos,
                     drugs, non_empty_rs_id_to_ref_seq_difference_annotation)

    def test_haplotype_with_repeat_rs_ids(self) -> None:
        """Error when haplotype has multiple variants with the same rs id"""
        name = "*5"
        function = "No Function"

        non_clashing_variants = frozenset([Variant("rs88293", "A"), Variant("rs39492", "T")])
        clashing_variants = frozenset([Variant("rs88293", "A"), Variant("rs88293", "T")])

        Haplotype(name, function, non_clashing_variants)
        with self.assertRaises(ValueError):
            Haplotype(name, function, clashing_variants)

    def test_empty_haplotype(self) -> None:
        """Error when haplotype has no variants"""
        name = "*5"
        function = "No Function"

        variants = frozenset([Variant("rs88293", "A"), Variant("rs39492", "T")])
        empty_variants: FrozenSet[Variant] = frozenset()

        Haplotype(name, function, variants)
        with self.assertRaises(ValueError):
            Haplotype(name, function, empty_variants)

    def test_panel_with_repeat_gene_names(self) -> None:
        """Error when panel has info for gene multiple times"""
        name = "FakePanel"
        version = "1.0"

        gene1 = "FAKE"
        chromosome1 = "X"
        gene2 = "OTHER"
        chromosome2 = "15"

        reference_haplotype_name = "*1"
        haplotypes: FrozenSet[Haplotype] = frozenset()
        rs_id_infos: FrozenSet[RsIdInfo] = frozenset()
        drugs: FrozenSet[DrugInfo] = frozenset()
        rs_id_to_ref_seq_difference_annotation: Dict[str, str] = dict()

        gene_info1 = GeneInfo(gene1, chromosome1, reference_haplotype_name, haplotypes,
                              rs_id_infos, drugs, rs_id_to_ref_seq_difference_annotation)
        gene_info2 = GeneInfo(gene2, chromosome2, reference_haplotype_name, haplotypes,
                              rs_id_infos, drugs, rs_id_to_ref_seq_difference_annotation)
        gene_info3 = GeneInfo(gene1, chromosome2, reference_haplotype_name, haplotypes,
                              rs_id_infos, drugs, rs_id_to_ref_seq_difference_annotation)
        gene_infos_without_repeat = frozenset([gene_info1, gene_info2])
        gene_infos_with_repeat = frozenset([gene_info1, gene_info3])

        Panel(name, version, gene_infos_without_repeat)
        with self.assertRaises(ValueError):
            Panel(name, version, gene_infos_with_repeat)

    def test_panel_with_overlapping_rs_id_infos_for_different_genes(self) -> None:
        """Error when panel has overlapping rs id infos for different genes, but not when they are exactly the same"""
        name = "FakePanel"
        version = "1.0"

        gene1 = "FAKE"
        gene2 = "OTHER"

        chromosome = "X"
        reference_haplotype_name = "*1"
        haplotypes: FrozenSet[Haplotype] = frozenset()
        drugs: FrozenSet[DrugInfo] = frozenset()
        rs_id_to_ref_seq_difference_annotation: Dict[str, str] = dict()

        rs_id_info1 = RsIdInfo("rs294924", "AT", "AT",
                               GeneCoordinate(chromosome, 499593), GeneCoordinate(chromosome, 399483))
        rs_id_info2 = RsIdInfo("rs3949923", "C", "C",
                               GeneCoordinate(chromosome, 293993), GeneCoordinate(chromosome, 1388323))
        rs_id_info3 = RsIdInfo("rs12993", "GG", "GG",
                               GeneCoordinate(chromosome, 499592), GeneCoordinate(chromosome, 399482))

        rs_id_infos1 = frozenset([rs_id_info1])
        rs_id_infos2 = frozenset([rs_id_info1, rs_id_info2])
        rs_id_infos3 = frozenset([rs_id_info3])
        gene_info1 = GeneInfo(gene1, chromosome, reference_haplotype_name, haplotypes,
                              rs_id_infos1, drugs, rs_id_to_ref_seq_difference_annotation)
        gene_info2 = GeneInfo(gene2, chromosome, reference_haplotype_name, haplotypes,
                              rs_id_infos2, drugs, rs_id_to_ref_seq_difference_annotation)
        gene_info3 = GeneInfo(gene2, chromosome, reference_haplotype_name, haplotypes,
                              rs_id_infos3, drugs, rs_id_to_ref_seq_difference_annotation)

        Panel(name, version, frozenset([gene_info1, gene_info2]))
        with self.assertRaises(ValueError):
            Panel(name, version, frozenset([gene_info1, gene_info3]))


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
