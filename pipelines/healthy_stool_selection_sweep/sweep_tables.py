#!/usr/bin/env python3
"""Per-sample and per-arm tables for the HSTOOL length / score / k-mer sweeps.

Metric definitions are those of analyses/minlen10_vs_9_2026-09-03/src/compare_cohort.py
(read_sage): PSMs = target rows with spectrum_q <= 0.01; decoy fraction over all rows at
spectrum_q <= 0.01; peptides = distinct mod-stripped sequences with peptide_q <= 0.01;
proteins = distinct non-decoy accessions with protein_q <= 0.01. Razor and stage-3 numbers
come from the chain's own JSON. Arms whose stage 1 pooled 58 samples are directly comparable;
`ml9_pool97` (noclust_2026-08-07) pooled 97 and is reported for reference only.

    sweep_tables.py OUT_DIR [--threads N]
"""
import sys, os, json, csv, re, statistics as st, argparse
from multiprocessing import Pool
csv.field_size_limit(1 << 30)
CANON = "LAB_ROOT/canonical"
ARMS_ROOT = f"{CANON}/publication_data/arms"
LIST58 = f"{CANON}/publication/analyses/fig1_revision_2026-09-02/cohorts/hstool_sample_list_ruled58.tsv"
MOD_RE = re.compile(r"\[[^\]]*\]")

def arms():
    out = []
    for L in (6, 7, 8, 11, 12):
        out.append((f"ml{L}", "length", L, 0.30, f"{ARMS_ROOT}/peptide_length_sweep_runs/minlen{L}/cohorts/HSTOOL"))
    out.append(("ml10", "length", 10, 0.30, f"{ARMS_ROOT}/10peptide_length_runs/cohorts/HSTOOL"))
    out.append(("ml9_pool97", "reference", 9, 0.30, f"{ARMS_ROOT}/noclust_2026-08-07/cohorts/HSTOOL"))
    for P in ("0.05", "0.10", "0.15", "0.20", "0.30", "0.40", "0.50", "0.60", "0.70", "0.80", "0.90", "1.00"):
        out.append((f"tp{P}", "score", 9, float(P), f"{ARMS_ROOT}/score_sweep_runs/top{P}/cohorts/HSTOOL"))
    for P in ("1.00", "0.30"):
        out.append((f"tp{P}_kmer", "kmer", 9, float(P), f"{ARMS_ROOT}/score_sweep_runs/top{P}_kmer/cohorts/HSTOOL"))
    return out

def samples():
    with open(LIST58) as fh:
        r = csv.reader(fh, delimiter="\t"); next(r)
        return [row[0] for row in r if row and row[0].strip()]

def jload(p):
    with open(p) as fh: return json.load(fh)

def read_sage(path):
    psms = decoy = 0; peps = set(); prots = set(); sc = []
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            lab = int(r["label"]); sq = float(r["spectrum_q"])
            if sq <= 0.01:
                if lab == 1: psms += 1
                else: decoy += 1
            if lab != 1: continue
            if float(r["peptide_q"]) <= 0.01: peps.add(MOD_RE.sub("", r["peptide"]))
            if sq <= 0.01: sc.append(float(r["scored_candidates"]))
            if float(r["protein_q"]) <= 0.01:
                for a in r["proteins"].split(";"):
                    if a and not a.startswith("rev_"): prots.add(a)
    return dict(psms_q01=psms, decoy_frac_q01=(decoy / (psms + decoy)) if psms + decoy else float("nan"),
                peptides_q01=len(peps), proteins_q01=len(prots),
                scored_candidates_median=st.median(sc) if sc else float("nan"))

def one(task):
    arm, kind, L, P, root, sample = task
    row = dict(arm=arm, kind=kind, min_length=L, top_percent=P, sample=sample)
    try:
        s3 = jload(f"{root}/3_persample/{sample}/sample_lake.stats.json")
        row.update(s3_seeds=s3["total_peptides"], s3_match_rate=s3["match_rate"], s3_proteins=s3["proteins_matched"])
        rz = jload(f"{root}/4_parsimony/razor_hash-acc/{sample}/{sample}_razor_stats.json")
        row.update(razor_selected=rz["selected_proteins"], razor_unique_peptides=rz["unique_peptides"],
                   razor_mapped_peptides=rz["total_peptides_mapped"])
        row.update(read_sage(f"{root}/5_sage/razor_hash-acc/{sample}/results.sage.tsv"))
        row["status"] = "ok"
    except Exception as e:
        row["status"] = f"missing: {type(e).__name__}: {e}"
    return row

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("out"); ap.add_argument("--threads", type=int, default=8)
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
    sm = samples(); tasks = []; cohort = []
    for arm, kind, L, P, root in arms():
        if not os.path.isdir(root):
            print("skip (no root):", arm, root); continue
        try:
            s1 = jload(f"{root}/1_evidence/evidence_lake.stats.json")
            cohort.append(dict(arm=arm, kind=kind, min_length=L, top_percent=P, s1_pooled_peptides=s1["total_peptides"],
                               s1_evidence_proteins=s1["total_output_proteins"], s1_score_filter=s1.get("score_filter"),
                               n_samples_pooled=97 if arm == "ml9_pool97" else 58))
        except Exception as e:
            cohort.append(dict(arm=arm, kind=kind, min_length=L, top_percent=P, s1_pooled_peptides="", s1_evidence_proteins="",
                               s1_score_filter=f"missing {e}", n_samples_pooled=""))
        tasks += [(arm, kind, L, P, root, s) for s in sm]
    with Pool(a.threads) as pool: rows = pool.map(one, tasks, chunksize=4)
    cols = ["arm", "kind", "min_length", "top_percent", "sample", "s3_seeds", "s3_match_rate", "s3_proteins", "razor_selected",
            "razor_unique_peptides", "razor_mapped_peptides", "psms_q01", "decoy_frac_q01", "peptides_q01", "proteins_q01",
            "scored_candidates_median", "status"]
    with open(f"{a.out}/per_sample_by_arm.tsv", "w", newline="") as g:
        w = csv.DictWriter(g, fieldnames=cols, delimiter="\t"); w.writeheader()
        for r in rows: w.writerow({k: r.get(k, "") for k in cols})
    with open(f"{a.out}/cohort_by_arm.tsv", "w", newline="") as g:
        w = csv.DictWriter(g, fieldnames=list(cohort[0].keys()), delimiter="\t"); w.writeheader(); w.writerows(cohort)
    metrics = cols[5:16]
    with open(f"{a.out}/arm_summary.tsv", "w", newline="") as g:
        w = csv.writer(g, delimiter="\t")
        w.writerow(["arm", "kind", "min_length", "top_percent", "n_ok", "metric", "median", "q1", "q3", "min", "max"])
        for arm, kind, L, P, root in arms():
            ok = [r for r in rows if r["arm"] == arm and r["status"] == "ok"]
            if not ok: continue
            for m in metrics:
                v = sorted(float(r[m]) for r in ok if r.get(m, "") != "" and r[m] == r[m])
                if not v: continue
                q = st.quantiles(v, n=4) if len(v) >= 2 else [v[0], v[0], v[0]]
                w.writerow([arm, kind, L, P, len(ok), m, st.median(v), q[0], q[2], v[0], v[-1]])
    bad = [r for r in rows if r["status"] != "ok"]
    print(f"rows {len(rows)}, ok {len(rows) - len(bad)}, missing {len(bad)}")
    for r in bad[:10]: print("  ", r["arm"], r["sample"], r["status"][:120])

if __name__ == "__main__":
    main()
