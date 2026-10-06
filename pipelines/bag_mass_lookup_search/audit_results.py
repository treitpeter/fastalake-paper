"""Independent table reconciliation; preserves unfinished/failed acquisitions."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent


def main():
    output = HERE / "audit"
    output.mkdir(exist_ok=True)
    manifest = json.loads((HERE / "MANIFEST.json").read_text())
    mechanism = pd.read_csv(HERE.parent / "results/PEPTIDE_MECHANISMS.tsv", sep="\t")
    summaries, locals_, receipts, sources = [], [], [], []
    metrics = ["legacy_ics_raw_z1_20ppm", "sage_backbone_coverage", "sage_unique_by_coverage"]
    keys = ["cohort", "sample", "arm", "filename", "scannr"]
    for index, case in enumerate(manifest["cases"]):
        root = HERE / "runs" / str(index)
        if not (root / "VALIDATION.json").exists():
            continue
        receipt = json.loads((root / "VALIDATION.json").read_text())
        assert receipt["status"] == "PASS"
        psm = pd.read_csv(root / "PSM_COVERAGE.tsv", sep="\t")
        pairs = pd.read_csv(root / "PAIRED_EVIDENCE.tsv", sep="\t")
        assert not psm.duplicated(keys).any()
        assert not pairs.duplicated(["cohort", "sample", "contrast", "change", "filename", "scannr"]).any()
        assert len(psm) == receipt["scored_psms"] and len(pairs) == receipt["paired_observations"]
        for metric in metrics:
            assert psm[metric].notna().all() and psm[metric].between(0,1).all()
        assert psm.sage_backbone_coverage.gt(0).all(), "Missing fragment annotations cannot be silently zero"
        assert (psm.sage_unique_by_coverage * 2*(psm.length-1) <= psm.matched_peaks + 1e-6).all()
        for metric, denominator in (("legacy_ics_raw_z1_20ppm", 2*(psm.length-1)), ("sage_backbone_coverage", psm.length-1)):
            values = psm[metric] * denominator
            assert np.allclose(values, np.round(values))
        own = pairs.merge(psm, left_on=["cohort","sample","own_arm","filename","scannr"], right_on=keys, validate="many_to_one")
        assert len(own) == len(pairs) and (own.own_peptide == own.peptide).all()
        assert (own.label == 1).all() and (own.peptide_q <= .01).all()
        for metric in metrics:
            assert np.allclose(own["own_"+metric], own[metric])
        for _, row in pairs[pairs.other_peptide.notna()].iterrows():
            for metric in metrics + ["hyperscore"]:
                assert abs(row["delta_"+metric] - (row["own_"+metric]-row["other_"+metric])) < 1e-7
        same = pairs[pairs.own_peptide == pairs.other_peptide]
        assert np.allclose(same.delta_legacy_ics_raw_z1_20ppm, 0)
        diagnostics = pairs[pairs.pairwise_diagnostic_site_ions.notna()]
        assert (diagnostics.pairwise_diagnostic_site_ions <= diagnostics.pairwise_site_ions).all()
        assert diagnostics.pairwise_diagnostic_site_ions.ge(0).all()
        expected = mechanism[(mechanism.cohort == case["cohort"]) & (mechanism["sample"] == case["sample"])]
        for (contrast, change), group in pairs.groupby(["contrast", "change"]):
            n = int(expected[(expected.contrast == contrast) & (expected.change == change)].own_accepted_spectra.sum())
            assert len(group) == n, (case, contrast, change, len(group), n)
            reported = group[group.other_peptide.notna()]
            summaries.append(dict(cohort=case["cohort"], sample=case["sample"], contrast=contrast, change=change,
                spectrum_events=len(group), comparator_reported=len(reported), comparator_missing=int(group.other_peptide.isna().sum()),
                decoy_comparators=int((group.relation == "decoy_comparator").sum()),
                identical_peptidoform=int((group.own_peptide == group.other_peptide).sum()),
                median_delta_hyperscore=float(reported.delta_hyperscore.median()),
                median_delta_legacy_ics=float(reported.delta_legacy_ics_raw_z1_20ppm.median()),
                median_delta_backbone=float(reported.delta_sage_backbone_coverage.median())))
        for (contrast, change, relationship), group in diagnostics.groupby(["contrast", "change", "relation"]):
            locals_.append(dict(cohort=case["cohort"], sample=case["sample"], contrast=contrast, change=change, relation=relationship,
                spectrum_events=len(group), zero_diagnostic=int(group.pairwise_diagnostic_site_ions.eq(0).sum()),
                positive_diagnostic=int(group.pairwise_diagnostic_site_ions.gt(0).sum()),
                median_diagnostic=float(group.pairwise_diagnostic_site_ions.median())))
        receipts.append(dict(index=index, cohort=case["cohort"], sample=case["sample"], psms=len(psm), pairs=len(pairs)))
        for name in ("PSM_COVERAGE.tsv", "PAIRED_EVIDENCE.tsv", "VALIDATION.json", "REPLAY_COMPARISON.json"):
            path = root / name
            sources.append(dict(path=str(path), sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    pd.DataFrame(summaries).to_csv(output / "COMPARISON_SUMMARY.tsv", sep="\t", index=False)
    pd.DataFrame(locals_).to_csv(output / "LOCAL_DIAGNOSTIC_SUMMARY.tsv", sep="\t", index=False)
    result = dict(status="PASS_COMPLETED_CASES", cases=receipts, sources=sources, utc=datetime.now(timezone.utc).isoformat(),
                  checks=["Per-contrast spectrum counts reconcile to frozen gain/loss table", "Own assignment and coverage joins reproduce PSM table",
                          "No missing or out-of-range coverage; positive Sage annotated coverage", "Discrete ion/cleavage denominators",
                          "Deduplicated Sage ions do not exceed reported matched peaks", "Identical peptidoforms have identical raw-spectrum ICS",
                          "Comparator missingness retained; deltas reconcile; diagnostic ions bounded by site-spanning ions"])
    (output / "VALIDATION.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(dict(status=result["status"], cases=len(receipts))))


if __name__ == "__main__":
    main()
