"""Copy the inputs of Supplementary Table 2 into ``source_data/STab2/``.

Run once by the unit author, not by the notebook. It reads the original producer
outputs under the lab project directory named by the environment variable
``FASTALAKE_LAB_ROOT`` (written ``LAB_ROOT`` in PROVENANCE.md) and writes plain copies:

* the four per-cohort tables of the full-route control experiment
  (``LAB_ROOT/canonical/publication/analyses/full_route_entrapment_2026-09-06/tables/<cohort>/``),
  concatenated in cohort order CAMPI_MOCK, HSTOOL, PREDICT, CAPSCAN with one header,
  exactly as ``build_package.py`` did with ``pandas.concat``; text is concatenated line
  by line so no value is re-formatted;
* ``control_discoveries.tsv`` per cohort (byte copies; the concatenation exceeds 20 MB);
* the completed fusion tables and the CAMPI independent-selection tables that
  ``package_edit3.py`` merged into the workbook as ``fusion_*`` and ``indep_*`` sheets
  (byte copies of the archived inputs of ``package_streamlining_sources.py`` and
  ``package_v10_sources.py``);
* the healthy-stool pooling table behind the Methods sentence at v28d line 743
  (byte copy of ``FIGURATOR/Figure3/review_2026-09-07/entrapment_pooling/per_acquisition.tsv``,
  identical to Source Data package ``Figure3/per_acquisition.tsv``).

Every copy is recorded with its origin sha256 in ``work/tables/STab2/prepare_receipt.json``.
"""
import hashlib
import json
import os
import shutil
from pathlib import Path

UNIT = Path(__file__).resolve().parent
REPO = UNIT.parents[1]
SD = REPO / "source_data" / UNIT.name
LAB = Path(os.environ["FASTALAKE_LAB_ROOT"])
WORK = Path(os.environ["FASTALAKE_WORK_DIR"])

COHORTS = ["CAMPI_MOCK", "HSTOOL", "PREDICT", "CAPSCAN"]
FULL_ROUTE = "canonical/publication/analyses/full_route_entrapment_2026-09-06/tables"
ARCHIVE = "../BiorXiv/99_archive/2026-09-11_before_weekend_organisation/review_2026-09-08_v10_AP/source_data"
POOLING = "FIGURATOR/Figure3/review_2026-09-07/entrapment_pooling/per_acquisition.tsv"


def sha(p):
    """Return the sha256 hex digest of file ``p``."""
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def main():
    """Write every source_data file and a receipt of origins and hashes."""
    SD.mkdir(parents=True, exist_ok=True)
    receipt = []
    for name in ["per_acquisition", "stage_control_counts", "baseline_comparison", "summary"]:
        out = SD / ("quantile_summary.tsv" if name == "summary" else name + ".tsv")
        origins = []
        with open(out, "w", newline="") as dst:
            for i, c in enumerate(COHORTS):
                src = LAB / FULL_ROUTE / c / (name + ".tsv")
                lines = src.read_text().splitlines(keepends=True)
                dst.writelines(lines if i == 0 else lines[1:])
                origins.append({"origin": f"LAB_ROOT/{FULL_ROUTE}/{c}/{name}.tsv", "sha256": sha(src)})
        receipt.append({"file": out.name, "sha256": sha(out), "from": origins, "how": "line concatenation, one header"})
    for c in COHORTS:
        src = LAB / FULL_ROUTE / c / "control_discoveries.tsv"
        out = SD / f"control_discoveries_{c}.tsv"
        shutil.copyfile(src, out)
        receipt.append({"file": out.name, "sha256": sha(out),
                        "from": [{"origin": f"LAB_ROOT/{FULL_ROUTE}/{c}/control_discoveries.tsv", "sha256": sha(src)}],
                        "how": "byte copy"})
    copies = {
        "fusion_cohort_q001.tsv": "completed_fusion/cohort_q001.tsv",
        "fusion_cohort_threshold_summary.tsv": "completed_fusion/cohort_threshold_summary.tsv",
        "fusion_matched_seed_summary.tsv": "completed_fusion/matched_seed_summary.tsv",
        "fusion_short_5_to_8_summary.tsv": "completed_fusion/short_5_to_8_summary.tsv",
        "fusion_acquisition_outliers_q001.tsv": "completed_fusion/acquisition_outliers_q001.tsv",
        "fusion_completeness.tsv": "completed_fusion/completeness.tsv",
        "fusion_all_acquisitions.tsv": "completed_fusion/snapshot/per_acquisition.tsv",
        "fusion_validation.json": "completed_fusion/VALIDATION.json",
        "indep_ordinary_counts.tsv": "campi_independent_selection/ordinary_counts.tsv",
        "indep_fusion_all_thresholds.tsv": "campi_independent_selection/fusion_diagnostics.tsv",
        "indep_validation.json": "campi_independent_selection/VALIDATION.json",
    }
    for dst, rel in copies.items():
        src = LAB / ARCHIVE / rel
        shutil.copyfile(src, SD / dst)
        receipt.append({"file": dst, "sha256": sha(SD / dst),
                        "from": [{"origin": "LAB_ROOT/" + ARCHIVE + "/" + rel, "sha256": sha(src)}],
                        "how": "byte copy"})
    src = LAB / POOLING
    shutil.copyfile(src, SD / "hstool_pooling_per_acquisition.tsv")
    receipt.append({"file": "hstool_pooling_per_acquisition.tsv", "sha256": sha(SD / "hstool_pooling_per_acquisition.tsv"),
                    "from": [{"origin": "LAB_ROOT/" + POOLING, "sha256": sha(src)}], "how": "byte copy"})
    (WORK / "prepare_receipt.json").write_text(json.dumps(receipt, indent=2))
    for r in receipt:
        print(r["file"], r["sha256"])


if __name__ == "__main__":
    main()
