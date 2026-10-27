"""Reproduce updated v42 claims from their saved analysis tables.

Allocation ranges summarize the checked-annotation analyses, not the older
unverified annotation population. Rescoring and null-control summaries start
from all acquisition/perturbation or acquisition/permutation outcomes. This
does not rerun searches or the underlying allocation and nested model fits.
"""
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
SD = REPO / "source_data" / "v42_updates"
RULES = ("unique", "complete", "unique_complete")


def numbers():
    result = {}
    scope = json.loads((SD / "crossfit_scope.json").read_text())
    result["cache_capscan"] = int(re.search(r"pooled over (\d+) acquisitions", scope["cache_roster"]).group(1))
    assert scope["n_capscan"] == 295 and result["cache_capscan"] == 347
    a = pd.read_csv(SD / "allocation_medians.tsv", sep="\t")
    a = a[a.profile.eq("other_unresolved")]
    eq = a[a.variant.eq("equal")].allocated_fraction
    rs = a[a.variant.isin(RULES)].allocated_fraction
    result.update(allocation_restricted_min=100 * rs.min(), allocation_restricted_max=100 * rs.max(),
                  allocation_equal_min=100 * eq.min(), allocation_equal_max=100 * eq.max())

    own = pd.read_csv(SD / "allocation_own_stool.tsv", sep="\t")
    own = own[own.variant.isin(("equal", *RULES))]
    assert (own.own_closer == own.participants).all() and (own.participants == 15).all()
    result.update(allocation_own_n=int(own.participants.min()), allocation_own_p=float(own.p.max()))

    m = pd.read_csv(SD / "allocation_predict.tsv", sep="\t")
    host = m.layer_a.eq("host_metaP") | m.layer_b.eq("host_metaP")
    eq, rs = m.variant.eq("equal"), m.variant.isin(RULES)
    for prefix, mask in [("other_restricted", rs & ~host), ("other_equal", eq & ~host),
                         ("host_restricted", rs & host), ("host_equal", eq & host)]:
        result[prefix + "_min"] = float(m.loc[mask, "rho"].min())
        result[prefix + "_max"] = float(m.loc[mask, "rho"].max())
    result.update(other_p_max=float(m.loc[(eq | rs) & ~host, "p"].max()),
                  host_restricted_p_max=float(m.loc[rs & host, "p"].max()))

    h = pd.read_csv(SD / "allocation_hstool.tsv", sep="\t")
    h = h[h.profile.eq("other_unresolved") & h.variant.isin(("equal", *RULES))]
    tech = h[h.design.eq("TechnicalReplicates")].r_per_pair
    inter = h[h.design.eq("Interindividual")]
    restricted = inter[inter.variant.isin(RULES)].r_per_pair
    result.update(technical_min=float(tech.min()), technical_max=float(tech.max()),
                  inter_equal=float(inter.loc[inter.variant.eq("equal"), "r_per_pair"].iloc[0]),
                  inter_restricted_min=float(restricted.min()), inter_restricted_max=float(restricted.max()))

    folds = pd.read_csv(SD / "fold_fractions.tsv", sep="\t")
    assert len(folds) == 591 and not folds.duplicated(["cohort", "sample"]).any()
    pct = folds[[f"fold{i}_pct" for i in range(10)]]
    assert np.allclose(pct.sum(axis=1), 100, atol=.001)
    result.update(fold_acquisitions=len(folds), fold_min=float(pct.min().min()), fold_max=float(pct.max().max()))

    null = pd.read_csv(SD / "label_permutations.tsv", sep="\t")
    assert null.seed.nunique() == 1000 and null.filename.nunique() == 58
    assert not null.duplicated(["seed", "filename"]).any() and len(null) == 58000
    per_seed = null.groupby("seed").file_plain.max()
    result.update(null_permutations=int(null.seed.nunique()), null_max=int(per_seed.max()),
                  null_median_max=float(per_seed.median()), null_plus1_max=int(null.file_plus1.max()))

    d = pd.read_csv(SD / "numerical_per_acquisition.tsv", sep="\t")
    assert len(d) == 591 and not d.duplicated(["cohort", "sample"]).any()
    # Recalculate gains from the accepted FL/MG counts for each implementation.
    gain = lambda method, length: 100 * (d[f"FL_{method}_{length}"] / d[f"MG_{method}_{length}"].replace(0, np.nan) - 1)
    base = gain("P0", 5)
    diffs = pd.DataFrame({i: (gain(f"PERT{i}", 5) - base).abs() for i in range(1, 21)})
    changed = diffs.max(axis=1) > .5
    for cohort in ("PREDICT", "CAPSCAN"):
        result["ulp_changed_" + cohort] = int((changed & d.cohort.eq(cohort)).sum())
    max_spread = 0.
    for cohort in ("PREDICT", "CAPSCAN"):
        mask = d.cohort.eq(cohort)
        for length in (5, 9):
            medians = [gain(f"PERT{i}", length)[mask].median() for i in range(1, 21)]
            max_spread = max(max_spread, max(medians) - min(medians))
    result.update(ulp_perturbations=20, ulp_max_median_spread=float(max_spread), ulp_acquisitions=len(d))

    eds = pd.read_csv(SD / "eds_validation.tsv", sep="\t")
    perms = pd.read_csv(SD / "eds_permutations.tsv", sep="\t")
    assert set(eds.seed) == {20260915, 20260916, 20260917}
    assert (eds.permutations == 199).all() and (eds.n_boot == 2000).all()
    result.update(eds_seeds=int(eds.seed.nunique()), eds_permutations=int(eds.permutations.min()),
                  eds_bootstraps=int(eds.n_boot.min()), eds_permutation_rows=len(perms))
    return {k: v.item() if hasattr(v, "item") else v for k, v in result.items()}


if __name__ == "__main__":
    print(json.dumps(numbers(), indent=2, allow_nan=False))
