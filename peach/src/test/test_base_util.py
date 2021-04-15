import unittest
from typing import List, Tuple, Dict, Set

from base.gene_coordinate import GeneCoordinate
from base.util import get_key_to_multiple_values, get_covered_coordinates, replace_file_extension_of_path, \
    get_chromosome_name_to_index


class TestBaseUtil(unittest.TestCase):
    def test_get_key_to_multiple_values_empty(self) -> None:
        key_value_pairs: List[Tuple[str, int]] = []
        result = get_key_to_multiple_values(key_value_pairs)
        result_expected: Dict[str, List[int]] = {}
        self.assertEqual(result_expected, result)

    def test_get_key_to_multiple_values_non_empty(self) -> None:
        key_value_pairs = [("A", 1), ("A", 2), ("A", 3), ("B", 1), ("B", 5), ("C", 1), ("D", 9)]
        result = get_key_to_multiple_values(key_value_pairs)

        result_expected = {
            "A": [1, 2, 3],
            "B": [1, 5],
        }

        self.assertEqual(result_expected, result)

    def test_get_covered_coordinates_empty(self) -> None:
        start_coordinate = GeneCoordinate("X", 17)
        result = get_covered_coordinates(start_coordinate, "")

        result_expected: Set[GeneCoordinate] = set()
        self.assertEqual(result_expected, result)

    def test_get_covered_coordinates_single(self) -> None:
        start_coordinate = GeneCoordinate("X", 17)
        result = get_covered_coordinates(start_coordinate, "A")

        result_expected = {start_coordinate}
        self.assertEqual(result_expected, result)

    def test_get_covered_coordinates_multiple(self) -> None:
        start_coordinate = GeneCoordinate("X", 17)
        result = get_covered_coordinates(start_coordinate, "AGTA")

        result_expected = {
            GeneCoordinate("X", 17),
            GeneCoordinate("X", 18),
            GeneCoordinate("X", 19),
            GeneCoordinate("X", 20),
        }
        self.assertEqual(result_expected, result)

    def test_replace_file_extension_of_path_single_point(self) -> None:
        result = replace_file_extension_of_path("something.ext", "test")
        result_expected = "something.test"
        self.assertEqual(result_expected, result)

    def test_replace_file_extension_of_path_multiple_points(self) -> None:
        result = replace_file_extension_of_path("this.is.something.ext", "test")
        result_expected = "this.is.something.test"
        self.assertEqual(result_expected, result)

    def test_replace_file_extension_of_path_no_point(self) -> None:
        with self.assertRaises(AssertionError):
            replace_file_extension_of_path("something", "test")

    def test_get_chromosome_name_to_index(self) -> None:
        chromosome_name_to_index = get_chromosome_name_to_index()
        self.assertEqual(24, len(chromosome_name_to_index))

        for chrom in range(1, 22):
            with self.subTest(chrom=chrom):
                self.assertEqual(chrom, chromosome_name_to_index[str(chrom)])

        self.assertEqual(23, chromosome_name_to_index["X"])
        self.assertEqual(24, chromosome_name_to_index["Y"])


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
