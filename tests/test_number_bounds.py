"""Printed bounds retain their endpoints and require a finite measurement."""

import unittest

from tools.check_numbers import agrees_any


class NumberBoundTests(unittest.TestCase):
    def test_inclusive_endpoints(self):
        for printed in ("P ≤ 0.01", "P >= 0.01", "P <= 0.01", "P ≥ 0.01"):
            with self.subTest(printed=printed):
                self.assertTrue(agrees_any(printed, "0.01"))

    def test_strict_endpoints(self):
        for printed in ("P < 0.01", "P > 0.01"):
            with self.subTest(printed=printed):
                self.assertFalse(agrees_any(printed, "0.01"))

    def test_bounds_do_not_round_away_a_violation(self):
        for printed, measured in (("≤ 0.01", "0.010001"), ("≥ 0.01", "0.009999")):
            with self.subTest(printed=printed):
                self.assertFalse(agrees_any(printed, measured))

    def test_recomputed_value_must_be_a_measurement(self):
        for printed, measured in (("0.01", "< 0.01"), ("0.01", "> 0.01"),
                                  ("< 0.02", "> 0.01"), ("> 0.01", "< 0.02")):
            with self.subTest(printed=printed, measured=measured):
                self.assertFalse(agrees_any(printed, measured))

    def test_nonfinite_values_fail_without_crashing(self):
        for invalid in ("NaN", "sNaN", "Infinity", "-Infinity", "1e9999", "1e-9999"):
            for printed, measured in (("1.0", invalid), (invalid, invalid)):
                with self.subTest(printed=printed, measured=measured):
                    self.assertFalse(agrees_any(printed, measured))
