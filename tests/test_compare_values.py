"""Regressions for public CI failures and meaningful numerical discrepancies."""

import unittest

from tools.compare_values import same_numeric_value, same_prefixed_numeric_value


class NumericComparisonTests(unittest.TestCase):
    def test_public_ci_last_bit_difference(self):
        self.assertTrue(same_numeric_value("0.1561093348358916", "0.15610933483589162"))

    def test_same_printed_rounding_does_not_hide_a_changed_result(self):
        self.assertFalse(same_numeric_value("0.8951563068271853", "0.8922003883978662"))
        self.assertFalse(same_numeric_value("0.156109334", "0.156109335"))

    def test_small_probabilities_and_zero_remain_distinct(self):
        self.assertFalse(same_numeric_value("1e-20", "2e-20"))
        self.assertFalse(same_numeric_value(0, 1e-16))

    def test_non_finite_values_fail(self):
        for value in ("nan", "inf", "-inf"):
            self.assertFalse(same_numeric_value(value, value))

    def test_result_labels_are_exact(self):
        self.assertTrue(same_numeric_value("NA", "NA"))
        self.assertFalse(same_numeric_value("NA", "not measured"))

    def test_mantel_cpu_roundoff_preserves_the_label(self):
        self.assertTrue(same_prefixed_numeric_value(
            "Mantel rho 0.46298900555207534", "Mantel rho 0.4629890055520753", "Mantel rho "))
        self.assertFalse(same_prefixed_numeric_value(
            "Pearson r 0.4629890055520753", "Mantel rho 0.4629890055520753", "Mantel rho "))

    def test_mantel_changed_result_and_invalid_suffix_fail(self):
        for value in ("0.463", "nan", "inf", "NA", "0.4629890055520753 (n=77)"):
            self.assertFalse(same_prefixed_numeric_value(
                "Mantel rho " + value, "Mantel rho 0.4629890055520753", "Mantel rho "))


if __name__ == "__main__":
    unittest.main()
