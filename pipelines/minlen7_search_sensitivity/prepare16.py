"""min_len=7 FL_public configs for the 16 participant-15 A08_347 searches, by the M09 prepare.py transformation
(copy source results.json, drop version/output_paths, min_len 5->7, set output_directory; assert nothing else differs).
  --control : apply the same transformation to the 331 production A08 search configs with M09's output paths and
              require byte-equal JSON objects to the M09 FL_public configs (exit 1 otherwise).
  --build   : write configs/ and tasks.tsv (16 rows, index 0-15) for A08_347 indices 331-346.
"""
import json, os, sys, copy, csv, hashlib
F = "LAB_ROOT"
RUN = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
M09 = F + "/00_PACKAGE_2026-09-17/cleanup_2026-09-22/v10_mm_philip_comment_audit_2026-09-22/nine_point_followup/manuscript_final_pass_2026-09-23/mm_v11_revision_2026-09-23/reappraisal_2026-09-24/main_methods_pass_2026-09-24/review_decisions/M09_short_peptides/run_2026-09-24_minlen7"
A08N = os.path.join(os.path.dirname(RUN), "A08_347")
def sha(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()
def transform(srcp, od):
    src = json.load(open(srcp))
    assert src["version"] == "0.14.6" and src["database"]["enzyme"]["min_len"] == 5 and src["database"]["enzyme"]["max_len"] == 50
    new = copy.deepcopy(src); del new["version"], new["output_paths"]
    new["database"]["enzyme"]["min_len"] = 7; new["output_directory"] = od
    a, b = copy.deepcopy(src), copy.deepcopy(new)
    for d in (a, b):
        d.pop("version", None); d.pop("output_paths", None); d.pop("output_directory", None); d["database"]["enzyme"].pop("min_len")
    assert a == b
    return src, new
mode = sys.argv[1]
if mode == "--control":
    rows = [r for r in csv.DictReader(open(M09 + "/tasks.tsv"), delimiter="\t") if r["cohort"] == "CAPSCAN" and r["arm"] == "FL_public"]
    assert len(rows) == 331
    bad = [r["sample"] for r in rows if transform(r["source_config"], r["out_dir"])[1] != json.load(open(r["config"]))]
    print(json.dumps(dict(control="M09 transformation reproduces 331 CapScan FL_public configs", n=331, differing=len(bad), status="PASS" if not bad else "FAIL")))
    sys.exit(1 if bad else 0)
assert mode == "--build"
m = json.load(open(A08N + "/inputs/manifest.json")); new16 = [r["sample"] for r in m["samples"][331:]]; assert len(new16) == 16
tasks = []
for i, s in enumerate(new16):
    srcp = f"{A08N}/results/{s}/search/results.json"
    od = f"{RUN}/searches/{s}"; assert not os.path.exists(od), od
    src, new = transform(srcp, od)
    fasta = src["database"]["fasta"]; assert fasta == f"{A08N}/results/{s}/parsimony/{s}_razor.fasta", fasta
    cfg = f"{RUN}/configs/{i:04d}_CAPSCAN_FL_public_{s}.json"; json.dump(new, open(cfg, "w"), indent=2)
    tasks.append(dict(index=i, cohort="CAPSCAN", arm="FL_public", sample=s, lfq=src["quant"]["lfq"], fasta_mb=round(os.path.getsize(fasta) / 1e6, 2),
                      mem_gb=4, config=cfg, source_config=srcp, source_sha256=sha(srcp), baseline_dir=os.path.dirname(srcp), out_dir=od))
with open(f"{RUN}/tasks.tsv", "w") as fh:
    cols = list(tasks[0]); fh.write("\t".join(cols) + "\n")
    for t in tasks: fh.write("\t".join(str(t[c]) for c in cols) + "\n")
print(json.dumps(dict(built=len(tasks))))
