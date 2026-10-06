"""M09 grouping: run the EXECUTED production grouping (study.py sha256 5f69c457…) unchanged on fresh min_len-7 searches.
Usage: group_fresh.py <task_id>. Tasks mirror the 10 Sep manuscript_cohorts manifest (same acquisition rosters).
Builds a symlink search root (one dir per acquisition -> fresh output dir; PREDICT's five repaired acquisitions -> LFQ-on addendum outputs),
requires a PASS receipt for every acquisition, writes to grouping/out/<task>/ (new directory only)."""
import sys, os, json, csv, hashlib, importlib.util, collections
F = "LAB_ROOT"
RUN = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BASE = F + "/simplified_grouping_results/manuscript_cohorts_2026-09-10_2200"
STUDY = BASE + "/code/fasta_lake/study.py"
assert hashlib.sha256(open(STUDY, "rb").read()).hexdigest().startswith("5f69c457"), "study.py changed"
spec = importlib.util.spec_from_file_location("study", STUDY); study = importlib.util.module_from_spec(spec); spec.loader.exec_module(study)
man = json.load(open(BASE + "/manifest.json"))
TASKMAP = {"HSTOOL_noclust": ("HSTOOL", "FL"), "PREDICT_noclust": ("PREDICT", "FL_cohortMG"), "CAPSCAN_noclust": ("CAPSCAN", "FL_cohortMG"),
           "CAMPI_FL_UNIPROT": ("CAMPI", "FL"), "CAMPI_FL_REF": ("CAMPI", "FL_R"), "CAMPI_FL_REF_MG": ("CAMPI", "FL_R_M"),
           "CAMPI_REF": ("CAMPI", "REF"), "CAMPI_MG": ("CAMPI", "MG")}
tid = sys.argv[1]; cohort, arm = TASKMAP[tid]
t = [x for x in man["tasks"] if x["id"] == tid][0]
names = open(t["samples_file"]).read().split()
tasks = list(csv.DictReader(open(RUN + "/tasks.tsv"), delimiter="\t"))
idx = {(r["cohort"], r["arm"], r["sample"]): r for r in tasks}
out = f"{RUN}/grouping/out/{tid}"; root = f"{RUN}/grouping/roots/{tid}"
assert not os.path.exists(out), "output exists"
os.makedirs(root, exist_ok=False)
used = collections.Counter()
for n in names:
    r = idx.get((cohort, arm + "_lfq", n)) if cohort == "PREDICT" else None
    r = r or idx[(cohort, arm, n)]
    rec = json.load(open(f"{RUN}/receipts/receipt_{int(r['index']):04d}.json")); assert rec["pass"] is True, n
    assert os.path.exists(r["out_dir"] + "/lfq.tsv"), ("no lfq", n)
    os.symlink(r["out_dir"], f"{root}/{n}"); used[r["arm"]] += 1
summary = study.group_searches(root, out, config_name="results.json", acquisition_names=names)
base_summary = json.load(open(f"{BASE}/cohorts/{tid}/direct/summary.json"))
keys = ["acquisitions", "inferred_study_groups", "quantified_study_groups", "pooled_canonical_peptides", "lfq_features", "searched_universe", "features_representative_outside_local_fasta"]
comp = {k: dict(baseline=base_summary.get(k), fresh7=summary.get(k)) for k in keys}
json.dump(dict(task=tid, arms_used=used, comparison=comp, baseline_summary=f"{BASE}/cohorts/{tid}/direct/summary.json",
               note="baseline = 10 Sep manuscript grouping (for PREDICT this run still holds the five damaged acquisitions; compare with the repaired regrouping as well)"),
          open(f"{RUN}/grouping/COMPARE_{tid}.json", "w"), indent=1)
print(json.dumps(comp, indent=1))
