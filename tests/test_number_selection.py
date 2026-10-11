"""A mistyped unit or missing package must not produce a successful empty check."""

import contextlib
import io
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from tools import check_numbers


class NumberSelectionTests(unittest.TestCase):
    def test_unknown_unit_fails(self):
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            status = check_numbers.main(["--unit", "not_a_paper_unit"])
        self.assertEqual(status, 1)
        self.assertIn("no NUMBERS.tsv files matched", output.getvalue())

    def test_missing_tables_fail(self):
        with tempfile.TemporaryDirectory() as folder:
            with patch.object(check_numbers, "REPO", Path(folder)):
                with contextlib.redirect_stdout(io.StringIO()):
                    self.assertEqual(check_numbers.main([]), 1)
