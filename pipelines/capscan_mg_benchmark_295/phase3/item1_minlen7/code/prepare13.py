"""min_len=7 MG configs for the 13 mg295 participant-15 acquisitions (plus the 1311 control), by the M09
prepare.py transformation (copy source results.json, drop version/output_paths, min_len 5->7, set
output_directory; assert nothing else differs). Same transform as public347/minlen7_public347/code/prepare16.py.
  --control : apply the transformation to the 282 M09 CapScan MG source configs with M09's output paths and
              require byte-equal JSON objects to the archived M09 MG configs (exit 1 otherwise).
  --build   : write configs/ and tasks.tsv. index 0 = control 20210513_..._1311 (source: mg295/mg_full/out,
              the mg295 re-search of the production pair); indices 1-13 = mg295/inputs/ACQ13.tsv order.
              Source config for every task = mg295/mg_full/out/<s>/results.json (the min_len-5 MG search).
"""
import json, os, sys, copy, csv, hashlib
F = "LAB_ROOT"
RUN = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MG295 = os.path.dirname(os.path.dirname(RUN))
M09 = F + "/00_PACKAGE_2026-09-17/cleanup_2026-09-22/v10_mm_philip_comment_audit_2026-09-22/nine_point_followup/manuscript_final_pass_2026-09-23/mm_v11_revision_2026-09-23/reappraisal_2026-09-24/main_methods_pass_2026-09-24/review_decisions/M09_short_peptides/run_2026-09-24_minlen7"
CONTROL = "20210513_EXPL0_PeTr_SA_T760_1311"
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
def diffkeys(a, b, p=""):
    out = []
    if isinstance(a, dict) and isinstance(b, dict):
        for k in sorted(set(a) | set(b)): out += diffkeys(a.get(k), b.get(k), p + "/" + k)
    elif a != b: out.append(p)
    return out
mode = sys.argv[1]
if mode == "--control":
    rows = [r for r in csv.DictReader(open(M09 + "/tasks.tsv"), delimiter="\t") if r["cohort"] == "CAPSCAN" and r["arm"] == "MG"]
    assert len(rows) == 282
    bad = [r["sample"] for r in rows if transform(r["source_config"], r["out_dir"])[1] != json.load(open(r["config"]))]
    c = [r for r in rows if r["sample"] == CONTROL][0]
    _, mine = transform(f"{MG295}/mg_full/out/{CONTROL}/results.json", f"{RUN}/searches/{CONTROL}")
    d = diffkeys(json.load(open(c["config"])), mine)
    res = dict(control="M09 transformation reproduces the 282 CapScan MG configs", n=282, differing=len(bad),
               control_1311_config_vs_M09_differs_only_at=d,
               status="PASS" if not bad and set(d) <= {"/database/fasta", "/mzml_paths", "/output_directory"} else "FAIL")
    print(json.dumps(res)); sys.exit(0 if res["status"] == "PASS" else 1)
assert mode == "--build"
acq = [r["sample"] for r in csv.DictReader(open(MG295 + "/inputs/ACQ13.tsv"), delimiter="\t")]; assert len(acq) == 13
mem = int(os.environ.get("MEM_GB", "16"))
tasks = []
for i, s in enumerate([CONTROL] + acq):
    srcp = f"{MG295}/mg_full/out/{s}/results.json"
    od = f"{RUN}/searches/{s}"; assert not os.path.exists(od + "/results.sage.tsv"), od
    src, new = transform(srcp, od)
    fasta = src["database"]["fasta"]
    assert fasta.endswith(f"/{s}.metag_hashed_plus_human.fasta"), fasta
    cfg = f"{RUN}/configs/{i:04d}_CAPSCAN_MG_{s}.json"; json.dump(new, open(cfg, "w"), indent=2)
    tasks.append(dict(index=i, cohort="CAPSCAN", arm="MG", sample=s, lfq=src["quant"]["lfq"], fasta_mb=round(os.path.getsize(fasta) / 1e6, 2),
                      mem_gb=mem, config=cfg, source_config=srcp, source_sha256=sha(srcp), baseline_dir=os.path.dirname(srcp), out_dir=od))
with open(f"{RUN}/tasks.tsv", "w") as fh:
    cols = list(tasks[0]); fh.write("\t".join(cols) + "\n")
    for t in tasks: fh.write("\t".join(str(t[c]) for c in cols) + "\n")
print(json.dumps(dict(built=len(tasks))))
