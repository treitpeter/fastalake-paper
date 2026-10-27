"""Run one M09 search. Usage: run_task.py <index> [--pilot]. Writes a receipt JSON; exits non-zero on any failed check."""
import sys, os, json, csv, hashlib, subprocess, time, re, socket, copy
RUN = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SAGE = "LAB_ROOT/MM_review_analysis/A13_4_7_predict_unspecific/bin/sage"
SAGE_SHA = "d3820543a31bfa2a556e04f91719204ce0ac6f34dcca6670314a13df87064625"
csv.field_size_limit(10**9)

def sha(p, bs=1 << 24):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(bs), b""): h.update(b)
    return h.hexdigest()

idx = int(sys.argv[1]); pilot = "--pilot" in sys.argv
rows = list(csv.DictReader(open(f"{RUN}/tasks.tsv"), delimiter="\t"))
t = rows[idx]; assert int(t["index"]) == idx
cfg = json.load(open(t["config"]))
src = json.load(open(t["source_config"]))
assert sha(t["source_config"]) == t["source_sha256"]
out = f"{RUN}/pilot/{idx:04d}" if pilot else t["out_dir"]
if os.path.exists(out + "/results.sage.tsv"):
    sys.exit(f"refusing to overwrite existing output {out}")
os.makedirs(out, exist_ok=True)
cfg["output_directory"] = out
run_cfg = out + "/config.json"; json.dump(cfg, open(run_cfg, "w"), indent=2)
rec = dict(index=idx, cohort=t["cohort"], arm=t["arm"], sample=t["sample"], pilot=pilot, host=socket.gethostname(),
           slurm_job=os.environ.get("SLURM_JOB_ID"), array_job=os.environ.get("SLURM_ARRAY_JOB_ID"),
           array_task=os.environ.get("SLURM_ARRAY_TASK_ID"), out_dir=out, config=run_cfg)
rec["sage_sha256"] = sha(SAGE); assert rec["sage_sha256"] == SAGE_SHA
rec["sage_version"] = subprocess.run([SAGE, "--version"], capture_output=True, text=True).stdout.strip()
rec["fasta"] = cfg["database"]["fasta"]; rec["fasta_sha256"] = sha(rec["fasta"])
rec["mzml"] = cfg["mzml_paths"][0]; rec["mzml_bytes"] = os.path.getsize(rec["mzml"]); rec["mzml_sha256"] = sha(rec["mzml"])
t0 = time.time(); rec["started"] = time.strftime("%Y-%m-%dT%H:%M:%S%z")
p = subprocess.run([SAGE, run_cfg, "--write-pin"], stdout=open(out + "/sage.stdout", "w"), stderr=open(out + "/sage.stderr", "w"))
rec["sage_exit"] = p.returncode; rec["elapsed_s"] = round(time.time() - t0, 1)
checks = {}
try:
    res = json.load(open(out + "/results.json"))
    checks["results_json_min_len_7"] = res["database"]["enzyme"]["min_len"] == 7
    a = copy.deepcopy(res); b = copy.deepcopy(src)
    for d in (a, b):
        for k in ("version", "output_paths", "output_directory"): d.pop(k, None)
        d["database"]["enzyme"].pop("min_len")
    checks["resolved_settings_equal_source_except_min_len"] = a == b
    checks["version_0.14.6"] = res["version"] == "0.14.6"
    n = acc = short = 0; lens = {}
    for r in csv.DictReader(open(out + "/results.sage.tsv"), delimiter="\t"):
        n += 1
        if r["rank"] == "1" and r["label"] == "1" and float(r["peptide_q"]) <= 0.01:
            acc += 1; L = len(re.sub(r"\[[^\]]*\]|[^A-Z]", "", r["peptide"]))
            lens[L] = lens.get(L, 0) + 1; short += L < 7
    rec.update(psm_rows=n, accepted_psms=acc, accepted_length_hist=lens)
    checks["accepted_psms_positive"] = acc > 0
    checks["no_accepted_peptide_below_7"] = short == 0
    checks["pin_written"] = os.path.getsize(out + "/results.sage.pin") > 0
    checks["lfq_as_source"] = os.path.exists(out + "/lfq.tsv") == bool(src["quant"]["lfq"])
except Exception as e:
    checks["exception"] = repr(e)
rec["checks"] = checks
rec["pass"] = p.returncode == 0 and all(v is True for v in checks.values())
rec["finished"] = time.strftime("%Y-%m-%dT%H:%M:%S%z")
rd = f"{RUN}/pilot" if pilot else f"{RUN}/receipts"
json.dump(rec, open(f"{rd}/receipt_{idx:04d}.json", "w"), indent=1)
print(json.dumps(rec)[:2000]); sys.exit(0 if rec["pass"] else 1)
