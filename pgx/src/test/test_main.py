import unittest


class TestMain(unittest.TestCase):
    def test_trivial(self):
        self.assertEqual(5, 5)

    def test_trivial_false(self):
        self.assertEqual(5, 1)

    def test_trivial_boolean(self):
        self.assertEqual(False, False)


if __name__ == '__main__':
    unittest.main()
