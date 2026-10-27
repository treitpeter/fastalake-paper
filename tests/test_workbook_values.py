"""Keep Excel roundoff tolerance confined to a declared result column."""
import tempfile
import unittest
from pathlib import Path

import openpyxl

from tables.STab6.panels import compare_xlsx_values


class WorkbookComparisonTests(unittest.TestCase):
    def test_only_the_declared_ratio_can_differ_by_roundoff(self):
        with tempfile.TemporaryDirectory() as directory:
            paths = [Path(directory) / name for name in ("expected.xlsx", "actual.xlsx")]
            def write(path, count, ratio):
                book = openpyxl.Workbook()
                sheet = book.active
                sheet.title = "Acquisition_seed_rows"
                sheet.append(["accepted", "ratio_to_exact"])
                sheet.append([count, ratio])
                book.save(path)
                book.close()
            write(paths[0], 178, 0.9368421052631579)
            write(paths[1], 178, 0.936842105263158)
            allowed = {"Acquisition_seed_rows": {"ratio_to_exact"}}
            self.assertFalse(compare_xlsx_values(paths[1], paths[0])["identical"])
            self.assertTrue(compare_xlsx_values(paths[1], paths[0], allowed)["identical"])
            write(paths[1], 179, 0.936842105263158)
            self.assertFalse(compare_xlsx_values(paths[1], paths[0], allowed)["identical"])
            write(paths[1], 178, 0.94)
            self.assertFalse(compare_xlsx_values(paths[1], paths[0], allowed)["identical"])
            with self.assertRaises(ValueError):
                compare_xlsx_values(paths[1], paths[0], {"Acquisition_seed_rows": {"typo"}})


if __name__ == "__main__":
    unittest.main()
