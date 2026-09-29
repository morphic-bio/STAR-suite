import unittest

from compare_sam_records import compare


class SamRecordsTests(unittest.TestCase):
    def test_identical(self):
        self.assertEqual(compare(b"a\nb\n", b"a\nb\n")["status"], "identical")

    def test_order_only(self):
        self.assertEqual(compare(b"a\nb\na\n", b"b\na\na\n")["status"], "order_only")

    def test_duplicate_not_discarded(self):
        self.assertEqual(compare(b"a\nb\na\n", b"b\na\n")["status"], "different")

    def test_field_change_fails(self):
        self.assertEqual(compare(b"a\t1\n", b"a\t2\n")["status"], "different")

    def test_empty_matches_only_empty(self):
        self.assertEqual(compare(b"", b"")["status"], "identical")
        self.assertEqual(compare(b"a\n", b"")["status"], "different")


if __name__ == "__main__":
    unittest.main()
