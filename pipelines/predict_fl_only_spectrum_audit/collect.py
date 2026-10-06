#!/usr/bin/env python
"""Pool the 296 per-acquisition audits into scoring_counts_public_vs_specimen.tsv and SUMMARY.json."""
import csv, hashlib, json, os, sys, time
A = os.environ["FASTALAKE_OUT_DIR"].rstrip("/") + "/"
CODE = os.path.dirname(os.path.abspath(__file__))
CATS = ["A_absent", "D_displaced", "B_not_cleavable", "C_top_fails_q", "E_no_target"]
def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(8 << 20), b""): h.update(b)
    return h.hexdigest()
m = json.load(open(A + "inputs/manifest.json"))
samples = [t["sample"] for t in m["tasks"]]; assert len(samples) == 296
rows = []; missing = []; inputs = {}
for s in samples:
    p = A + "out/%s.json" % s
    if not os.path.exists(p): missing.append(s); continue
    j = json.load(open(p)); assert j["sample"] == s
    assert sum(j["counts"].values()) == j["n_fl_only"] == j["n_fl_spectra_q01"] - j["n_shared_spectra_q01"]
    assert j["fl_max_rank"] == 5 and j["mg_max_rank"] == 5, (s, j["fl_max_rank"], j["mg_max_rank"])
    rows.append(dict(sample=s, total=j["n_fl_only"], **{k: j["counts"][k] for k in CATS},
                     n_fl_spectra_q01=j["n_fl_spectra_q01"], n_mg_spectra_q01=j["n_mg_spectra_q01"], n_shared_spectra_q01=j["n_shared_spectra_q01"]))
    inputs[s] = dict(FL_results_sha256=j["inputs"]["FL_results"]["sha256"], MG_results_sha256=j["inputs"]["MG_results"]["sha256"],
                     MG_fasta_sha256=j["comparator_fasta"]["sha256"], MG_fasta_bytes=j["comparator_fasta"]["bytes"], audit_json_sha256=sha(p))
if missing:
    print("INCOMPLETE: %d acquisitions missing" % len(missing), missing[:10]); sys.exit(1)
with open(A + "scoring_counts_public_vs_specimen.tsv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t"); w.writeheader(); w.writerows(rows)
tot = {k: sum(r[k] for r in rows) for k in CATS}; N = sum(tot.values())
assert N == sum(r["total"] for r in rows)
summary = dict(status="PASS", analysis=m["analysis"], created=time.strftime("%Y-%m-%dT%H:%M:%S%z"),
    arms=dict(FL="public-only FastaLake, abl_nometag_noclust 4_parsimony/razor_hash-acc/<acq>/parsimony.fasta",
              MG="specimen metagenome (sha256 headers) + human Swiss-Prot, metag_hashed_human/PREDICT/_db/<acq>.metag_hashed_plus_human.fasta"),
    search="Sage 0.14.6 (sha256 %s), R8A_loss_psm5 configuration: production settings, lfq false, report_psms 5" % m["sage"]["sha256"],
    definitions="R8A_loss_psm5.sbatch: best TARGET row per spectrum by hyperscore among reported ranks; FL-only = FL peptide_q <= 0.01 and MG absent or peptide_q > 0.01; residue letters only, exact I/L; R3_digest tryptic screen (KR, not before P, 2 missed cleavages, 5-50 aa, 500-5000 Da, unmodified mass)",
    n_acquisitions=len(rows), n_spectra_fl_only=N, counts=tot,
    percent={k: round(100.0 * tot[k] / N, 3) for k in CATS},
    n_fl_spectra_q01=sum(r["n_fl_spectra_q01"] for r in rows), n_mg_spectra_q01=sum(r["n_mg_spectra_q01"] for r in rows),
    n_shared_spectra_q01=sum(r["n_shared_spectra_q01"] for r in rows),
    roster_source=m["roster_source"], original_pipeline=m["original_pipeline"], digest_module=dict(path=m["digest_module"], sha256=sha(m["digest_module"])),
    code=dict(audit_py_sha256=sha(os.path.join(CODE, "audit.py")), task_sh_sha256=sha(os.path.join(CODE, "task.sh")), collect_py_sha256=sha(__file__), manifest_sha256=sha(A + "inputs/manifest.json")),
    input_hashes=inputs, table=A + "scoring_counts_public_vs_specimen.tsv", table_sha256=sha(A + "scoring_counts_public_vs_specimen.tsv"),
    job_id=os.environ.get("SLURM_JOB_ID"))
json.dump(summary, open(A + "SUMMARY.json", "w"), indent=1)
print("A32 COLLECT PASS n=%d spectra=%d" % (len(rows), N), {k: summary["percent"][k] for k in CATS})
