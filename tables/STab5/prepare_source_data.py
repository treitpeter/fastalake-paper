"""Write source_data/STab5/ from the original analysis outputs (run once by the authors).

This script is not part of the reviewer path; the notebook reads only
``source_data/STab5/``. It documents exactly how each shipped file was made from the
original files under ``LAB_ROOT`` (the lab project directory, given by the environment
variable ``FASTALAKE_LAB_ROOT``). It also needs two controlled files (the PREDICT sample
roster with raw patient identifiers and the package pseudonym key), read from
``FASTALAKE_CONTROLLED_DIR`` or, if unset, from their LAB_ROOT locations.

Files written
-------------
Byte-identical copies of the analysis outputs (no transformation):

- ``crossomics_same_function_correlations.tsv``  revised allocation, 30,581 tests
- ``same_function_summary.tsv``                  nine-family summary
- ``eligibility.tsv``                            every candidate function and layer pair
- ``previous_allocation_same_function_correlations.tsv``  previous allocation, 17,914 tests

Patient-level inputs of the correlation step, pseudonymised (one file per level and
layer, long format, non-zero values only, exact float64 values):

- ``profiles_<level>_<layer>.tsv.gz`` with columns ``patient``, ``feature``, ``value``
- ``profile_features.tsv``: every feature column of those matrices (``level``, ``layer``,
  ``feature``), including features that are zero in all 77 patients
- ``patients.tsv``: the 77 pseudonymised patients

The profile rows are exactly the matrices that ``individual_functions.py`` (refit mode)
indexes with ``.loc[first.index, terms]``: the 77 PREDICT patients with complete
DNA/RNA/protein/metabolite data, restricted to features that occur in at least one
layer-pair vocabulary of that level. Protein profiles are first collapsed from
acquisitions to patients by the original ``collapse`` function (mean over the patient's
acquisitions). Raw patient identifiers are replaced by the Source Data package pseudonyms
(``PREDICT_Pnnn``) and rows are sorted by pseudonym, so that row order does not reveal the
key.

Usage (on Slurm)::

    FASTALAKE_LAB_ROOT=... python tables/STab5/prepare_source_data.py
"""

import hashlib
import json
import os
import shutil
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

UNIT = Path(__file__).resolve().parent
REPO = UNIT.parents[1]
SD = REPO / "source_data" / UNIT.name
LAB = Path(os.environ["FASTALAKE_LAB_ROOT"])
CONTROLLED = os.environ.get("FASTALAKE_CONTROLLED_DIR")

BR = LAB / "review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit"
V50 = BR / "claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50"
REFIT = V50 / "B/individual_refit"
PROFILES = V50 / "C/profiles/PREDICT_noclust"
OLD = LAB / "analyses/enterotypes_2026-09-07"
KEYDIR = BR / "claude_review_2026-10-03/final/figure_sd/source_data/internal_only"

COPIES = {
    "crossomics_same_function_correlations.tsv": REFIT / "crossomics_same_function_correlations.tsv",
    "same_function_summary.tsv": REFIT / "same_function_summary.tsv",
    "eligibility.tsv": REFIT / "eligibility.tsv",
    "previous_allocation_same_function_correlations.tsv": BR / "individual_baseline/crossomics_same_function_correlations.tsv",
}
LEVELS = ["OG", "KO", "Pathway"]
LAYERS = ["metaG", "metaT", "microbial_metaP"]


def sha(p):
    h = hashlib.sha256()
    with Path(p).open("rb") as f:
        for b in iter(lambda: f.read(8 << 20), b""):
            h.update(b)
    return h.hexdigest()


def rd(p, **kw):
    """Read a TSV exactly as the original analysis did (pandas default float parser)."""
    return pd.read_csv(p, sep="\t", **kw)


def frame(folder, key):
    """Original ``common.frame``: dense samples x features frame from a sparse npz."""
    a = json.loads((folder / (key + "_axes.json")).read_text())
    v = sparse.load_npz(folder / (key + ".npz"))
    return pd.DataFrame(v.toarray(), index=a["samples"], columns=a["features"])


def collapse(d):
    """Original ``individual_functions.collapse``: mean over acquisitions of each patient."""
    ids = d.index.to_series().str.extract(r"(MicrobPred_\d+)", expand=False)
    assert ids.notna().all()
    return d.groupby(ids).mean()


def controlled(name, fallback):
    p = Path(CONTROLLED) / name if CONTROLLED else fallback
    if not p.is_file():
        sys.exit(f"needs controlled input {name}")
    return p


def main():
    SD.mkdir(parents=True, exist_ok=True)
    receipt = {"copies": {}, "profiles": {}}
    for name, src in COPIES.items():
        shutil.copyfile(src, SD / name)
        assert sha(src) == sha(SD / name)
        receipt["copies"][name] = {"origin": "LAB_ROOT/" + str(src.relative_to(LAB)), "sha256": sha(src)}

    first_path = controlled("PREDICT_fullstack_primary_samples.tsv", OLD / "PREDICT_fullstack_primary_samples.tsv")
    key_path = controlled("PREDICT_PARTICIPANT_KEY.tsv", KEYDIR / "PREDICT_PARTICIPANT_KEY.tsv")
    first = rd(first_path, index_col=0)
    assert len(first) == first.idpatient.nunique() == 77
    key = rd(key_path, dtype=str).set_index("study_code").package_code
    idp = first.idpatient.astype(int).astype(str)
    pseudo = idp.map(key)
    assert pseudo.notna().all() and pseudo.nunique() == 77
    receipt["inputs"] = {"primary_samples_sha256": sha(first_path), "participant_key_sha256": sha(key_path)}

    vocab = []
    pd.DataFrame({"patient": sorted(pseudo)}).to_csv(SD / "patients.tsv", sep="\t", index=False)
    for lev in LEVELS:
        fs = dict(
            metaG=rd(OLD / ("metaG_OG_TPM.tsv" if lev == "OG" else "PREDICT_metaG_" + lev + ".tsv"), index_col=0),
            metaT=rd(OLD / ("metaT_OG_TPM.tsv" if lev == "OG" else "PREDICT_metaT_" + lev + ".tsv"), index_col=0),
            microbial_metaP=collapse(frame(PROFILES, lev + "_other_unresolved")),
        )
        used = {k: set() for k in fs}
        for a, b in combinations(fs, 2):
            terms = set(fs[a]) & set(fs[b])
            used[a] |= terms
            used[b] |= terms
        for lay, d in fs.items():
            cols = sorted(used[lay])
            vocab += [(lev, lay, c) for c in cols]
            m = d.loc[first.index, cols]
            assert np.isfinite(m.to_numpy()).all(), (lev, lay)
            m.index = pseudo.loc[m.index].to_numpy()
            long = m.stack()
            long = long[long != 0]
            long.index.names = ["patient", "feature"]
            long = long.rename("value").reset_index().sort_values(["patient", "feature"], kind="stable")
            out = SD / f"profiles_{lev}_{lay}.tsv.gz"
            long.to_csv(out, sep="\t", index=False, compression={"method": "gzip", "mtime": 0})
            receipt["profiles"][out.name] = {"rows": int(len(long)), "features": len(cols), "patients": 77,
                                             "sha256": sha(out), "bytes": out.stat().st_size}
            print(lev, lay, len(cols), len(long), out.stat().st_size, flush=True)
    pd.DataFrame(vocab, columns=["level", "layer", "feature"]).to_csv(SD / "profile_features.tsv", sep="\t", index=False)
    (Path(os.environ.get("STAB5_WORK", ".")) / "PREPARE_RECEIPT.json").write_text(json.dumps(receipt, indent=1) + "\n")


if __name__ == "__main__":
    main()
