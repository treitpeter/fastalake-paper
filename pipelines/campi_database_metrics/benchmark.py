#!/usr/bin/env python3
"""Recalculate CAMPI coverage on accepted archived PSMs; no new identification search.
Ion-series coverage collapses charge: observed b/y series positions / 2(L-1).
Backbone coverage: cleavage positions supported by b OR y / (L-1).
Sage v0.14.6 scoring and spectrum preparation are reproduced and checked against matched_peaks.
"""

import csv
import gzip
import hashlib
import json
import os
import re
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
from numba import njit
from pyteomics import mzml

HERE = Path(os.environ["FIG2_WORK"])
CONFIG = json.loads((HERE / "RUN_INPUTS.json").read_text())
ARMS = CONFIG["arms"]
SAMPLES = CONFIG["samples"]
MASS = np.array(
    [
        71.03711,
        0,
        103.00919,
        115.02694,
        129.04259,
        147.0684,
        57.02146,
        137.05891,
        113.08406,
        0,
        128.09496,
        113.08406,
        131.0405,
        114.04293,
        237.14774,
        97.05276,
        128.05858,
        156.1011,
        87.03203,
        101.04768,
        150.95363,
        99.06841,
        186.07932,
        0,
        163.06332,
        0,
    ],
    dtype=np.float32,
)
PROTON = np.float32(1.0072764)
TOK = re.compile(r"([A-Z])(?:\[([+-]?[0-9.]+)\])?")


def parse(s):
    t = TOK.findall(s)
    assert "".join(a + ("[" + m + "]" if m else "") for a, m in t) == s, s
    return "".join(a for a, m in t), np.array(
        [np.float32(MASS[ord(a) - 65] + np.float32(m or 0)) for a, m in t], dtype=np.float32
    )


def canon(s):
    return re.sub(r"\[[^\]]*\]", "", s).replace("I", "L")


@njit(cache=True)
def match(masses, ints, residues, total, charge):
    n = len(residues)
    b = np.float32(0)
    y = np.float32(total)
    bmask = np.uint64(0)
    ymask = np.uint64(0)
    matched = 0
    eligible = 0
    for i in range(n - 1):
        b = np.float32(b + residues[i])
        y = np.float32(y - residues[i])
        for series in range(2):
            neutral = b if series == 0 else y
            for z in range(1, max(2, min(charge, 3))):
                center = np.float32(neutral / np.float32(z))
                delta = np.float32(np.float32(center * np.float32(20)) / np.float32(1000000))
                lo = np.float32(center - delta)
                hi = np.float32(center + delta)
                if center + PROTON >= 200 and center + PROTON <= 2000:
                    eligible += 1
                index = np.searchsorted(masses, lo)
                found = False
                while index < len(masses) and masses[index] <= hi:
                    if ints[index] >= 0:
                        found = True
                    index += 1
                if found:
                    matched += 1
                    if series == 0:
                        bmask |= np.uint64(1) << np.uint64(i)
                    else:
                        ymask |= np.uint64(1) << np.uint64(i)
    return bmask, ymask, matched, eligible


def load_fasta(p):
    out = {}
    acc = None
    buf = []
    with Path(p).open() as f:
        for line in f:
            if line.startswith(">"):
                if acc is not None:
                    out[acc] = "".join(buf).replace("I", "L")
                acc = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if acc is not None:
        out[acc] = "".join(buf).replace("I", "L")
    return out


def protein_metrics(evidence, db):
    out = []
    for acc, peps in sorted(evidence.items()):
        seq = db.get(acc)
        assert seq is not None, acc
        mask = np.zeros(len(seq), dtype=bool)
        first = np.zeros(len(seq), dtype=bool)
        for p in peps:
            start = seq.find(p)
            assert start >= 0, (acc, p)
            first[start : start + len(p)] = True
            while start >= 0:
                mask[start : start + len(p)] = True
                start = seq.find(p, start + 1)
        out.append(
            {
                "accession": acc,
                "protein_length": len(seq),
                "n_peptides": len(peps),
                "covered_residues_all_occurrences": int(mask.sum()),
                "coverage_pct": float(100 * mask.mean()),
                "coverage_first_occurrence_pct": float(100 * first.mean()),
                "peptides": ";".join(sorted(peps)),
            }
        )
    return out


def run(sample, limit=0):
    start = time.time()
    out = HERE / ("benchmark_pilot" if limit else "benchmark") / sample
    out.mkdir(parents=True, exist_ok=True)
    man = CONFIG["manifest"]
    wanted = defaultdict(lambda: defaultdict(list))
    evidence = {k: defaultdict(set) for k in ARMS}
    counts = {}
    configs = {}
    prov = []
    for k, arm in ARMS.items():
        row = next(r for r in man if r["arm"] == arm["arm"] and r["sample"] == sample)
        d = Path(row["sagedir"])
        cfg = json.loads((d / "config.json").read_text())
        configs[k] = cfg
        assert (
            cfg["database"]["ion_kinds"] == ["b", "y"]
            and not cfg["deisotope"]
            and cfg["max_peaks"] == 150
            and cfg["fragment_tol"] == {"ppm": [-20, 20]}
            and cfg["max_fragment_charge"] == 2
        )
        assert (
            cfg["database"]["fragment_min_mz"] == 200 and cfg["database"]["fragment_max_mz"] == 2000
        )
        n = 0
        accepted_peptides = set()
        matched_peaks = []
        matched_intensity = []
        with (d / "results.sage.tsv").open() as f:
            for r in csv.DictReader(f, delimiter="\t"):
                if r["label"] != "1" or float(r["peptide_q"]) > 0.01:
                    continue
                assert r["rank"] == "1"
                bare, m = parse(r["peptide"])
                assert len(bare) == int(r["peptide_len"]) and 5 <= len(bare) <= 50
                n += 1
                accepted_peptides.add(canon(r["peptide"]))
                matched_peaks.append(int(r["matched_peaks"]))
                matched_intensity.append(float(r["matched_intensity_pct"]))
                key = (r["peptide"], int(r["charge"]))
                wanted[r["scannr"]][key].append(
                    (
                        k,
                        int(r["matched_peaks"]),
                        float(r["matched_intensity_pct"]),
                        np.float32(r["calcmass"]),
                        r["psm_id"],
                    )
                )
                if float(r["protein_q"]) <= 0.01:
                    for a in r["proteins"].split(";"):
                        if a and not a.startswith(("rev_", "DECOY_")):
                            evidence[k][a].add(canon(r["peptide"]))
        counts[k] = n
        prov.append(
            {
                "strategy": k,
                "sage": str(d / "results.sage.tsv"),
                "sha256": hashlib.sha256((d / "results.sage.tsv").read_bytes()).hexdigest(),
                "fasta": cfg["database"]["fasta"],
            }
        )
        pr = protein_metrics(evidence[k], load_fasta(cfg["database"]["fasta"]))
        with gzip.open(out / f"{k}_protein_coverage.tsv.gz", "wt") as f:
            w = csv.DictWriter(f, fieldnames=list(pr[0]), delimiter="\t")
            w.writeheader()
            w.writerows(pr)
        counts[k] = {
            "peptides": len(accepted_peptides),
            "peptides_per_protein_mean": float(np.mean([r["n_peptides"] for r in pr])),
            "peptides_per_protein_median": float(np.median([r["n_peptides"] for r in pr])),
            "matched_peaks_median": float(np.median(matched_peaks)),
            "matched_intensity_pct_median": float(np.median(matched_intensity)),
            "psms": n,
            "proteins": len(pr),
            "protein_coverage_pct_median": float(np.median([r["coverage_pct"] for r in pr])),
            "protein_coverage_first_occurrence_pct_median": float(
                np.median([r["coverage_first_occurrence_pct"] for r in pr])
            ),
            "protein_coverage_pct_mean": float(np.mean([r["coverage_pct"] for r in pr])),
            "protein_ge50pct_coverage_pct": 100
            * sum(r["coverage_pct"] >= 50 for r in pr)
            / len(pr),
            "protein_ge2peptide_pct": 100 * sum(r["n_peptides"] >= 2 for r in pr) / len(pr),
        }
    mzpath = configs["FL"]["mzml_paths"][0]
    assert all(c["mzml_paths"] == [mzpath] for c in configs.values())
    print(sample, "reading", mzpath, "wanted spectra", len(wanted), flush=True)
    stats = {k: defaultdict(list) for k in ARMS}
    seen = set()
    n = 0
    totalmatches = 0
    mismatches = []
    handles = {k: gzip.open(out / f"{k}_fragment_coverage.tsv.gz", "wt") for k in ARMS}
    writers = {k: csv.writer(f, delimiter="\t") for k, f in handles.items()}
    head = [
        "scan",
        "peptide",
        "charge",
        "psm_id",
        "length",
        "matched_b_mask",
        "matched_y_mask",
        "series_coverage_pct",
        "backbone_coverage_pct",
        "matched_ion_charge_hypotheses",
        "eligible_ion_charge_hypotheses",
        "eligible_ion_coverage_pct",
        "sage_matched_peaks",
        "shared_all_strategies",
    ]
    for w in writers.values():
        w.writerow(head)
    with mzml.MzML(mzpath, use_index=False, iterative=True, decode_binary=True) as reader:
        for sp in reader:
            sid = sp["id"]
            if sid not in wanted:
                continue
            mz = np.asarray(sp["m/z array"], dtype=np.float32)
            it = np.asarray(sp["intensity array"], dtype=np.float32)
            keep = (mz >= 200) & (mz <= 2000)
            mz = mz[keep]
            it = it[keep]
            ix = np.argsort(it, kind="stable")[-150:]
            mz = mz[ix]
            it = it[ix]
            ix = np.argsort(mz)
            masses = np.asarray(mz[ix] - PROTON, dtype=np.float32)
            ints = it[ix]
            for (peptide, charge), records in wanted[sid].items():
                bare, res = parse(peptide)
                assert all(t[3] == records[0][3] for t in records)
                bm, ym, matched, eligible = match(masses, ints, res, records[0][3], charge)
                b = int(bm)
                y = int(ym)
                series = 100 * (b.bit_count() + y.bit_count()) / (2 * (len(bare) - 1))
                backbone = 100 * (b | y).bit_count() / (len(bare) - 1)
                eligiblecov = 100 * matched / eligible if eligible else float("nan")
                assert 0 <= series <= 100 and 0 <= backbone <= 100 and 0 <= eligiblecov <= 100
                all4 = len({r[0] for r in records}) == len(ARMS)
                for k, expected, intpct, total, pid in records:
                    totalmatches += 1
                    if matched != expected:
                        if len(mismatches) < 1000:
                            mismatches.append([sid, peptide, charge, k, expected, matched])
                        stats[k]["mismatch"].append(1)
                    writers[k].writerow(
                        [
                            sid,
                            peptide,
                            charge,
                            pid,
                            len(bare),
                            b,
                            y,
                            series,
                            backbone,
                            matched,
                            eligible,
                            eligiblecov,
                            expected,
                            int(all4),
                        ]
                    )
                    for key, value in [
                        ("series", series),
                        ("backbone", backbone),
                        ("eligible", eligiblecov),
                        ("length", len(bare)),
                    ]:
                        stats[k][key].append(value)
                    if all4:
                        stats[k]["shared_series"].append(series)
            seen.add(sid)
            n += 1
            if n % 10000 == 0:
                print(sample, n, "spectra", round(time.time() - start, 1), "seconds", flush=True)
            if limit and n >= limit:
                break
    for f in handles.values():
        f.close()
    if not limit:
        assert seen == set(wanted), (len(seen), len(wanted))
    for k in ARMS:
        if not limit:
            assert len(stats[k]["series"]) == counts[k]["psms"]
        counts[k].update(
            {
                key + "_coverage_pct_median": float(np.median(stats[k][key]))
                for key in ["series", "backbone", "eligible"]
            }
        )
        counts[k].update(
            {
                "fragment_rows": len(stats[k]["series"]),
                "matched_peaks_mismatches": len(stats[k]["mismatch"]),
                "shared_all_psms": len(stats[k]["shared_series"]),
                "shared_all_series_coverage_pct_median": float(
                    np.median(stats[k]["shared_series"])
                ),
                "peptide_length_median": float(np.median(stats[k]["length"])),
            }
        )
    result = {
        "sample": sample,
        "limit": limit,
        "seconds": time.time() - start,
        "mzml": mzpath,
        "n_spectra_scored": n,
        "n_spectra_expected": len(wanted),
        "arms": counts,
        "mismatch_examples": mismatches,
        "provenance": prov,
    }
    (out / "summary.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2)[:5000], flush=True)
    if not limit:
        assert all(v["matched_peaks_mismatches"] == 0 for v in counts.values()), (
            "Fragment counts disagree with Sage"
        )
    if not limit:
        (out / "COMPLETE").write_text("All requested spectra found; all accepted PSMs measured.\n")


if __name__ == "__main__":
    arg = sys.argv[1]
    sample = SAMPLES[int(arg)] if arg.isdigit() else arg
    run(sample, int(sys.argv[2]) if len(sys.argv) > 2 else 0)
