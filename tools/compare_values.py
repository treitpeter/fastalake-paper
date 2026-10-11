"""Compare recomputed numeric fields without requiring identical float formatting.

Use only for numeric result fields. Identity, source, scope, count, seed and
method fields must retain their explicit exact comparisons. Relative tolerance
1e-12 accommodates last-bit floating-point differences across supported
environments; zero has no absolute tolerance and non-finite numbers fail.
"""

import math


def same_numeric_value(actual, expected):
    """Accept finite equal numeric values within 1e-12 relative error.

    Non-numeric result labels, such as ``NA``, are compared exactly. This does
    not round results to manuscript precision or relax scientific thresholds.
    """
    try:
        left, right = float(actual), float(expected)
    except (TypeError, ValueError):
        return str(actual) == str(expected)
    return (
        math.isfinite(left)
        and math.isfinite(right)
        and math.isclose(left, right, rel_tol=1e-12, abs_tol=0.0)
    )
