"""A2 preflight. usage: preflight.py <searches|grouping331|grouping347> [--negative]
--negative plants one defect in memory (searches: a changed precursor tolerance in one config; grouping331: one receipt
marked failed; grouping347: one roster acquisition removed from the task table). The check must then exit 1."""
import sys, os, json, csv, copy, hashlib, importlib.util
F = "LAB_ROOT"
A = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RUN, RUN331 = A + "/run", A + "/run331"
M09 = F + "/00_PACKAGE_2026-09-17/cleanup_2026-09-22/v10_mm_philip_comment_audit_2026-09-22/nine_point_followup/manuscript_final_pass_2026-09-23/mm_v11_revision_2026-09-23/reappraisal_2026-09-24/main_methods_pass_2026-09-24/review_decisions/M09_short_peptides/run_2026-09-24_minlen7"
C2G = F + "/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig5/grouping"
STUDY = F + "/simplified_grouping_results/manuscript_cohorts_2026-09-10_2200/code/fasta_lake/study.py"
SAGE = F + "/MM_review_analysis/A13_4_7_predict_unspecific/bin/sage"
csv.field_size_limit(10**9)
mode = sys.argv[1]; neg = "--negative" in sys.argv
def sha(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()
def norm(c):
    c = copy.deepcopy(c)
    for k in ("output_directory", "mzml_paths", "version", "output_paths"): c.pop(k, None)
    c["database"].pop("fasta"); return c
def fail(msg): print("FAIL:", msg); sys.exit(1)
res = {}
if mode == "searches":
    rows = list(csv.DictReader(open(RUN + "/tasks.tsv"), delimiter="\t"))
    new = [r for r in rows if r["orig_index"] == "new"]
    if len(rows) != 347 or len(new) != 16: fail("task table size")
    ref = norm(json.load(open(rows[0]["config"])))            # an original M09 CapScan FL_cohortMG config
    if sha(SAGE) != "d3820543a31bfa2a556e04f91719204ce0ac6f34dcca6670314a13df87064625": fail("sage sha")
    for k, r in enumerate(new):
        c = json.load(open(r["config"]))
        if neg and k == 0: c["precursor_tol"] = {"ppm": [-10.0, 10.0]}
        if c["database"]["enzyme"]["min_len"] != 7: fail("min_len " + r["sample"])
        if norm(c) != ref: fail("settings differ from the original 331 configs: " + r["sample"])
        if sha(r["source_config"]) != r["source_sha256"]: fail("source sha " + r["sample"])
        if not (os.path.exists(c["database"]["fasta"]) and os.path.exists(c["mzml_paths"][0])): fail("inputs " + r["sample"])
        if c["output_directory"] != r["out_dir"] or os.path.exists(r["out_dir"] + "/results.sage.tsv"): fail("out_dir " + r["sample"])
    res = dict(new_searches=16, settings_equal_original=True)
elif mode in ("grouping331", "grouping347"):
    run = RUN331 if mode == "grouping331" else RUN
    if not sha(STUDY).startswith("5f69c457"): fail("study.py sha")
    spec = importlib.util.spec_from_file_location("study", STUDY); m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    assert hasattr(m, "group_searches")
    base = F + "/simplified_grouping_results/manuscript_cohorts_2026-09-10_2200" if mode == "grouping331" else C2G
    man = json.load(open(base + "/manifest.json")); t = [x for x in man["tasks"] if x["id"] == "CAPSCAN_noclust"][0]
    names = open(t["samples_file"]).read().split()
    if not os.path.exists(f"{base}/cohorts/CAPSCAN_noclust/direct/summary.json"): fail("baseline summary missing")
    rows = list(csv.DictReader(open(run + "/tasks.tsv"), delimiter="\t"))
    idx = {(r["cohort"], r["arm"], r["sample"]): r for r in rows}
    if neg and mode == "grouping347": idx.pop(("CAPSCAN", "FL_cohortMG", names[5]))
    pending = []
    for k, n in enumerate(names):
        r = idx.get(("CAPSCAN", "FL_cohortMG", n))
        if r is None: fail("no task for " + n)
        rp = f"{run}/receipts/receipt_{int(r['index']):04d}.json"
        if mode == "grouping347" and r.get("orig_index") == "new":
            pending.append(n); continue                      # produced by the search array (afterok)
        rec = json.load(open(rp))
        if neg and mode == "grouping331" and k == 7: rec["pass"] = False
        if rec["pass"] is not True: fail("receipt not PASS " + n)
        if not os.path.exists(r["out_dir"] + "/lfq.tsv"): fail("no lfq " + n)
    if os.path.exists(run + "/grouping/out/CAPSCAN_noclust"): fail("output exists")
    if mode == "grouping347":
        roster = {l.split("\t")[0] for l in open(C2G + "/../../roster/fixed_roster_347.tsv").read().splitlines()[1:]}
        if set(names) != roster or len(names) != 347: fail("manifest roster != fixed_roster_347")
        if len(pending) != 16: fail("pending new searches != 16")
    res = dict(acquisitions=len(names), receipts_checked=len(names) - len(pending), pending_new=len(pending))
else:
    fail("mode")
print("PASS", mode, json.dumps(res))
