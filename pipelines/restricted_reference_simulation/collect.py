"""Collect all planned units, preserving mixture-level replication and uncertainty."""

import json
import math
import sys

import numpy as np
import pandas as pd
from scipy.stats import t

from common import ROOT, guard, location, save, sha, table, units


def summarize(values, rng):
    values = np.asarray(values, dtype=float)
    n = len(values)
    mean = float(values.mean())
    if n < 2:
        return dict(
            n_units=n, mean=mean, bootstrap_low=None, bootstrap_high=None, t_low=None, t_high=None
        )
    draws = values[rng.integers(n, size=(20000, n))].mean(axis=1)
    half = float(t.ppf(0.975, n - 1) * values.std(ddof=1) / math.sqrt(n))
    return dict(
        n_units=n,
        mean=mean,
        bootstrap_low=float(np.quantile(draws, 0.025)),
        bootstrap_high=float(np.quantile(draws, 0.975)),
        t_low=mean - half,
        t_high=mean + half,
    )


def main(phase):
    frozen = guard()
    output = ROOT / "runs" / phase
    destination = output / "collected"
    destination.mkdir(exist_ok=False)
    primary, diagnostic, incremental, audits = [], [], [], {}
    for index, unit in enumerate(units(phase)):
        _, base = location(phase, index)
        audit = json.loads((base / "evaluation/AUDIT.json").read_text())
        for filename, expected in audit["outputs"].items():
            assert sha(base / "evaluation" / filename) == expected
        for filename, expected in audit["input_hashes"].items():
            assert sha(filename) == expected
        assert (
            audit["independent_recount_pass"]
            and audit["donor_exclusion_pass"]
            and audit["matched_truth_pass"]
        )
        audits[unit["id"]] = audit
        for name, target in [
            ("primary", primary),
            ("entrapment_vs_truth", diagnostic),
            ("incremental", incremental),
        ]:
            frame = pd.read_csv(base / "evaluation" / f"{name}.tsv", sep="\t")
            frame["unit_id"], frame["regime"] = unit["id"], unit["regime"]
            target.append(frame)
    primary = pd.concat(primary, ignore_index=True)
    diagnostic = pd.concat(diagnostic, ignore_index=True)
    incremental = pd.concat(incremental, ignore_index=True)
    for name, frame in [
        ("primary", primary),
        ("entrapment_vs_truth", diagnostic),
        ("incremental", incremental),
    ]:
        frame.to_csv(destination / f"{name}.tsv", sep="\t", index=False)
    rng = np.random.default_rng(2026091051)
    estimates = []
    ordinary = primary[primary.arm.isin(["fixed", "self", "independent"])]
    for (regime, arm), frame in ordinary.groupby(["regime", "arm"], sort=True):
        estimates.append(
            dict(
                regime=regime,
                endpoint=arm,
                **summarize(frame.realised_FDP, rng),
                total_discoveries=int(frame.discoveries.sum()),
                total_false=int(frame.false_discoveries.sum()),
                pooled_FDP=float(frame.false_discoveries.sum() / max(frame.discoveries.sum(), 1)),
            )
        )
    paired = ordinary.pivot(index=["regime", "unit_id"], columns="arm", values="realised_FDP")
    for regime, frame in paired.groupby(level="regime", sort=True):
        estimates.append(
            dict(
                regime=regime,
                endpoint="self_minus_independent",
                **summarize(frame.self - frame.independent, rng),
                total_discoveries=None,
                total_false=None,
                pooled_FDP=None,
            )
        )
    table(destination / "uncertainty.tsv", estimates)
    # Averages of control repeats are formed before mixtures are resampled.
    diagnostic["arm_family"] = diagnostic.arm.str.replace(r"_[123]$", "", regex=True)
    nested = diagnostic.groupby(
        ["regime", "unit_id", "arm_family", "q_multiplier", "q_threshold", "min_residues"],
        sort=True,
    )[["paired_upper_diagnostic_with_unpaired_worst_case", "actual_fdp"]].agg(
        lambda values: values.mean(skipna=False)
    )
    nested.to_csv(destination / "control_seed_means.tsv", sep="\t")
    bias = []
    for key, frame in nested.reset_index().groupby(
        ["regime", "arm_family", "q_multiplier", "q_threshold", "min_residues"], sort=True
    ):
        difference = frame.paired_upper_diagnostic_with_unpaired_worst_case - frame.actual_fdp
        if difference.isna().any():
            result = dict(
                n_units=len(frame),
                mean=None,
                bootstrap_low=None,
                bootstrap_high=None,
                t_low=None,
                t_high=None,
            )
        else:
            result = summarize(difference, rng)
        bias.append(
            dict(zip(["regime", "arm_family", "q_multiplier", "q_threshold", "min_residues"], key))
            | result
        )
    table(destination / "entrapment_bias.tsv", bias)
    ready = all(a["ordinary_outputs_assessable"] for a in audits.values())
    report = [
        "# Reference benchmark results",
        "",
        f"Phase: **{phase}**. Completed {len(audits)} independent mixtures, 12 searches each.",
        "",
        "The directly counted errors below describe this simulated reference domain.",
        "Peptide q ≤ 0.01; unmodified I/L-canonical sequences, lengths 5–50.",
        "",
        "| Unit | Regime | Arm | False / accepted | Observed FDP |",
        "|---|---|---|---:|---:|",
    ]
    for row in ordinary.itertuples():
        report.append(
            f"| {row.unit_id} | {row.regime} | {row.arm} | "
            f"{row.false_discoveries}/{row.discoveries} | {row.realised_FDP:.3%} |"
        )
    report.extend(
        [
            "",
            "Read `uncertainty.tsv` for unweighted unit means and paired comparisons;",
            "`entrapment_vs_truth.tsv` retains every control seed, threshold and q perturbation.",
            "`incremental.tsv` records the truth of extra and shared discoveries.",
            "",
            "The development gate tests integrity and nondegenerate output. It does not",
            "require favorable error rates. Simulation remains limited by the shared",
            "spectral-predictor/training context and restricted reference. Protein groups,",
            "cohort aggregation and measured-spectrum calibration are not established.",
            "",
        ]
    )
    (destination / "REPORT.md").write_text("\n".join(report))
    if phase == "development":
        save(
            output / "DEVELOPMENT_GATE.json",
            dict(
                ready_for_reserved_validation=ready,
                criteria="all planned units/searches, hashes, truth/recount/input-exclusion, nonempty ordinary output and decoy winners",
                outcome_direction_used_for_gate=False,
                units={k: a["ordinary_outputs_assessable"] for k, a in audits.items()},
                frozen_sha256=sha(ROOT / "FROZEN.json"),
                source=frozen["source"],
                report_sha256=sha(destination / "REPORT.md"),
            ),
        )
    save(
        destination / "COMPLETE.json",
        dict(
            phase=phase,
            n_units=len(audits),
            integrity_pass=True,
            all_ordinary_outputs_assessable=ready,
            whole_procedure_fdr_control_established=False,
            outputs={p.name: sha(p) for p in destination.iterdir() if p.is_file()},
        ),
    )
    print("\n".join(report), flush=True)


if __name__ == "__main__":
    main(sys.argv[1])
