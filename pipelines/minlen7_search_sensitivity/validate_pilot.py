"""Gate: every pilot receipt passes. Also records baseline comparison (not a gate) and Slurm MaxRSS."""
import json, os, sys, csv, glob, re, subprocess
RUN = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
csv.field_size_limit(10**9)
idx = [int(x) for x in open(f"{RUN}/pilot/pilot_indices.txt").read().split()]
tasks = list(csv.DictReader(open(f"{RUN}/tasks.tsv"), delimiter="\t"))
out = {"tasks": {}, "missing": []}
def peps(path, minlen):
    s = set()
    for r in csv.DictReader(open(path), delimiter="\t"):
        if r["rank"] == "1" and r["label"] == "1" and float(r["peptide_q"]) <= 0.01:
            p = re.sub(r"\[[^\]]*\]|[^A-Z]", "", r["peptide"]).replace("I", "L")
            if len(p) >= minlen: s.add(p)
    return s
for i in idx:
    f = f"{RUN}/pilot/receipt_{i:04d}.json"
    if not os.path.exists(f): out["missing"].append(i); continue
    rec = json.load(open(f)); t = tasks[i]
    fresh = peps(f"{RUN}/pilot/{i:04d}/results.sage.tsv", 7)
    nested = peps(t["baseline_dir"] + "/results.sage.tsv", 7)
    base = peps(t["baseline_dir"] + "/results.sage.tsv", 1)
    out["tasks"][i] = dict(route=f"{t['cohort']}/{t['arm']}", pass_=rec["pass"], checks=rec["checks"], elapsed_s=rec["elapsed_s"],
                           baseline_5_50=len(base), nested_ge7=len(nested), fresh_ge7=len(fresh),
                           shared=len(fresh & nested), fresh_only=len(fresh - nested), nested_only=len(nested - fresh),
                           slurm_job=rec["slurm_job"])
jobs = {v["slurm_job"] for v in out["tasks"].values()}
out["sacct"] = subprocess.run(["sacct", "-j", ",".join(j for j in jobs if j), "-o", "JobID,State,Elapsed,MaxRSS,ReqMem", "-P", "-n"],
                              capture_output=True, text=True).stdout
out["pass"] = not out["missing"] and all(v["pass_"] for v in out["tasks"].values())
json.dump(out, open(f"{RUN}/pilot/PILOT_VALIDATION.json", "w"), indent=1, default=str)
print(json.dumps(out, indent=1, default=str)[:4000]); sys.exit(0 if out["pass"] else 1)
