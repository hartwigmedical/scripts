import unittest

from main import load_panel_from_json
from panel import Panel
from test_resources.test_resource import get_panel_test_resource


class TestLoadConfig(unittest.TestCase):
    @unittest.skip("WIP")
    def test_load_panel(self) -> None:
        panel_path = get_panel_test_resource()
        panel = load_panel_from_json(str(panel_path))

        # gene_infos_expected = []
        # rs_id_infos_expected
        # panel_expected = Panel(gene_infos_expected, rs_id_infos_expected)
        self.fail("WIP")


if __name__ == '__main__':
    unittest.main()
