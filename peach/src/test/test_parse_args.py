import unittest
from argparse import Namespace
from copy import deepcopy
from io import StringIO
from unittest.mock import patch

from main import parse_args


class TestParseArgs(unittest.TestCase):
    def test_parse_args_representative_example(self) -> None:
        arguments = ["vcf_file", "tumor_sample_id", "ref_sample_id", "script_version", "output_directory",
                     "panel_location", "--transcript_tsv", "path_to_transcripts_tsv",
                     "vcftools_location"]
        actual_namespace = parse_args(arguments)

        expected_namespace = Namespace(
            vcf="vcf_file", sample_t_id="tumor_sample_id", sample_r_id="ref_sample_id", version="script_version",
            outputdir="output_directory", panel="panel_location", recreate_bed=False,
            vcftools="vcftools_location", transcript_tsv="path_to_transcripts_tsv"
        )
        self.assertEqual(expected_namespace, actual_namespace)

    def test_parse_args_defaults(self) -> None:
        arguments = ["vcf_file", "tumor_sample_id", "ref_sample_id", "script_version", "output_directory",
                     "panel_location", "vcftools_location"]
        actual_namespace = parse_args(arguments)

        expected_namespace = Namespace(
            vcf="vcf_file", sample_t_id="tumor_sample_id", sample_r_id="ref_sample_id", version="script_version",
            outputdir="output_directory", panel="panel_location", recreate_bed=False,
            vcftools="vcftools_location", transcript_tsv=None
        )
        self.assertEqual(expected_namespace, actual_namespace)

    def test_parse_args_no_defaults(self) -> None:
        arguments = ["vcf_file", "tumor_sample_id", "ref_sample_id", "script_version", "output_directory",
                     "panel_location", "--transcript_tsv", "path_to_transcripts_tsv",
                     "--recreate_bed", "vcftools_location"]
        actual_namespace = parse_args(arguments)

        expected_namespace = Namespace(
            vcf="vcf_file", sample_t_id="tumor_sample_id", sample_r_id="ref_sample_id", version="script_version",
            outputdir="output_directory", panel="panel_location", recreate_bed=True,
            vcftools="vcftools_location", transcript_tsv="path_to_transcripts_tsv"
        )
        self.assertEqual(expected_namespace, actual_namespace)

    def test_parse_args_missing_required_arguments(self) -> None:
        minimum_arguments = [
            "vcf_file", "tumor_sample_id", "ref_sample_id", "script_version",
            "output_directory", "panel_location", "vcftools_location"
        ]
        parse_args(minimum_arguments)

        too_few_arguments = deepcopy(minimum_arguments[:-1])
        with patch("sys.stderr", new_callable=StringIO):  # silence the help message in this test
            with self.assertRaises(SystemExit):
                parse_args(too_few_arguments)


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
