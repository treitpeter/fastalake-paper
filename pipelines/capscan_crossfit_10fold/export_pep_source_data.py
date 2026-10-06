"""Write the cross-fitted PEP sensitivity tables into source_data/supplementary_methods_6_7/.

Added for the paper repository (not part of the original run). Reads the outputs of
cf_analyze.py from FASTALAKE_OUT_DIR and, for the per-acquisition table, replaces the raw
participant identifier with the pseudonym already used in
source_data/supplementary_methods_6_7/crossfit_per_acquisition.tsv (joined on the
acquisition name) and blanks the PREDICT visit code in `group`, as in that table.
Values are not changed.

Usage: python export_pep_source_data.py <paper repository root>
"""
import os
import sys
from pathlib import Path

import pandas as pd


def main(repo):
    """Copy the two summaries and write the pseudonymised per-acquisition table."""
    src = Path(os.environ["FASTALAKE_OUT_DIR"])
    sd = Path(repo) / "source_data" / "supplementary_methods_6_7"
    for name, out in [("CF_PEP_SUMMARY.tsv", "crossfit_pep_summary.tsv"),
                      ("CF_PEP_SUMMARY_BY_REGION.tsv", "crossfit_pep_summary_by_region.tsv")]:
        (sd / out).write_bytes((src / name).read_bytes())
    a = pd.read_csv(src / "CF_PEP_PER_ACQUISITION.tsv", sep="\t", dtype=str, keep_default_na=False)
    key = pd.read_csv(sd / "crossfit_per_acquisition.tsv", sep="\t", dtype=str, keep_default_na=False)
    key = key.set_index("sample")["subject_id"]
    assert set(a["sample"]) == set(key.index) and len(a) == len(key) == 578
    a["subject_id"] = a["sample"].map(key)
    a.loc[a["cohort"] == "PREDICT", "group"] = ""
    a.to_csv(sd / "crossfit_pep_per_acquisition.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main(sys.argv[1])
