"""Copy the tables Supplementary Figure 6 needs from the lab archive into source_data/SFig06/.

This script documents how every file in ``source_data/SFig06/`` was made. It is not run
by the notebook and it cannot run outside the lab, because it reads the original
producer tables. Set ``FASTALAKE_LAB_ROOT`` to the lab project directory (``LAB_ROOT``
in PROVENANCE.md) to rerun it.

Transformations (and nothing else):

* ``split_half_acquisitions_726.tsv``: ``ACQUISITION_REVIEW.tsv`` of the spectrum-split
  experiment with the 347-acquisition CapScan re-run (726 acquisitions x 2 length views),
  with the columns ``participant`` and ``specimen`` removed. For PREDICT these columns
  hold raw patient and specimen identifiers (CONVENTIONS.md section 8). No panel or
  number of this figure uses them. (The earlier 710-acquisition table of the v28f figure
  is no longer used.)
* ``same_spectrum_benchmark_0_predict.tsv`` and ``same_spectrum_benchmark_1_capscan.tsv``:
  Source Data package sheets ``Figure4/benchmark_0_per_acquisition.tsv`` and
  ``Figure4/benchmark_1_per_acquisition.tsv``, with the column ``specimen_id`` removed
  (it carries a pseudonymised tube identifier that no panel uses). For PREDICT the column
  ``group`` is removed as well: it holds the study visit code (S = scheduled, A =
  admission, R = readmission visit), a clinical event next to a participant id, and no
  panel or number uses it (CONVENTIONS.md section 8). For CapScan ``group`` is the
  sampling region and is kept. ``subject_id`` is the package's pseudonymised participant
  id and is kept.
* ``split_half_structure.tsv``: per cohort, the number of acquisitions, distinct
  participants and distinct specimens at 9-50 residues, counted from the
  ``participant`` and ``specimen`` columns of the 726 table before they are removed, and
  the number of acquisitions whose identifier is empty. Aggregate counts only; no
  identifier is written. Empty identifiers are not counted as a participant or specimen.

Every input sha256 is asserted, so a changed archive file stops the script.
"""

import hashlib
import os
from pathlib import Path

import pandas as pd

LAB = Path(os.environ["FASTALAKE_LAB_ROOT"])
OUT = Path(__file__).resolve().parents[2] / "source_data" / "SFig06"

FINAL = ("review/v15_integration/jess_integration_2026-09-29/analysis_followup/"
         "manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/")
PKG = FINAL + "figure_sd/source_data/package/FastaLake_Source_Data/"

#: (input relative to LAB_ROOT, sha256, output name, columns removed)
JOBS = [
    (FINAL + "public347/spectrum_split347/collect/ACQUISITION_REVIEW.tsv",
     "e0f257046686fb7e551bf2455539ab2f917eec88e5d4deec9339401bb0340e79",
     "split_half_acquisitions_726.tsv", ["participant", "specimen"]),
    (PKG + "Figure4/benchmark_0_per_acquisition.tsv",
     "fdd5383ae384b7ec59202820a3e501ab2d80a32fd8e98318d5747d96ca667808",
     "same_spectrum_benchmark_0_predict.tsv", ["specimen_id", "group"]),
    (PKG + "Figure4/benchmark_1_per_acquisition.tsv",
     "65ac87033c11d58c20246a6a70ac747109ce901c10dc952974262c8051ff8e31",
     "same_spectrum_benchmark_1_capscan.tsv", ["specimen_id"]),
]


def sha256(path):
    """Return the hex sha256 of a file."""
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def structure(df):
    """Write ``split_half_structure.tsv``: acquisitions, participants and specimens per cohort.

    Parameters
    ----------
    df : pandas.DataFrame
        The 726 spectrum-split table with its ``participant`` and ``specimen`` columns,
        read as text.
    """
    g = df[df.min_length == "9"]
    rows = []
    for co, d in g.groupby("cohort", sort=False):
        rows.append(dict(cohort=co, n_acquisitions=len(d),
                         n_participants=int(d.participant[d.participant != ""].nunique()),
                         n_specimens=int(d.specimen[d.specimen != ""].nunique()),
                         n_acquisitions_without_participant=int((d.participant == "").sum()),
                         n_acquisitions_without_specimen=int((d.specimen == "").sum())))
    out = OUT / "split_half_structure.tsv"
    pd.DataFrame(rows).to_csv(out, sep="\t", index=False, lineterminator="\n")
    print(out.name, len(rows), sha256(out))


def main():
    """Copy, drop the listed columns, write tab-separated without index; print hashes."""
    OUT.mkdir(parents=True, exist_ok=True)
    for rel, digest, name, drop in JOBS:
        src = LAB / rel
        assert sha256(src) == digest, f"{name}: input changed"
        # keep_default_na=False keeps empty fields empty and strings such as 'NA' unchanged
        df = pd.read_csv(src, sep="\t", dtype=str, keep_default_na=False)
        if name == "split_half_acquisitions_726.tsv":
            structure(df)
        df = df.drop(columns=drop)
        df.to_csv(OUT / name, sep="\t", index=False, lineterminator="\n")
        print(name, len(df), sha256(OUT / name))


if __name__ == "__main__":
    main()
