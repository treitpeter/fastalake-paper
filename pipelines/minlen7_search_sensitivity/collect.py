"""A2 collector: writes RESULTS.md (old 331 vs new 347 for every printed number of the min-length-7 sensitivity).
usage: collect.py            real run (after the 16 searches and both groupings)
       collect.py --selftest  checks the collector against the ORIGINAL M09 outputs (repro path := original); writes nothing
       collect.py --selftest --negative  same, with one planted wrong reproduction value; must exit 1"""
import sys, os, json, glob, hashlib
import pandas as pd
F = "LAB_ROOT"
A = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = F + "/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03"
M09 = F + "/00_PACKAGE_2026-09-17/cleanup_2026-09-22/v10_mm_philip_comment_audit_2026-09-22/nine_point_followup/manuscript_final_pass_2026-09-23/mm_v11_revision_2026-09-23/reappraisal_2026-09-24/main_methods_pass_2026-09-24/review_decisions/M09_short_peptides/run_2026-09-24_minlen7"
XLSX = ROOT + "/final/figure_sd/source_data/package/FastaLake_Source_Data/Supplementary_Figure_04/LENGTH_SENSITIVITY_source_data.xlsx"
self_test = "--selftest" in sys.argv; neg = "--negative" in sys.argv
def sha(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()
fails = []
def check(ok, msg):
    if not ok: fails.append(msg)

# 1. packaged Source Data (the published values)
x = pd.read_excel(XLSX, sheet_name=None)
ab, c, e, g = x["a_b_benchmark_per_acquisition"], x["c_per_arm_change"], x["e_fusion"], x["g_grouping"]
med = lambda co, pop, col: ab[(ab.cohort == co) & (ab.population == pop)][col].median()
pk = g.set_index("task").loc["CAPSCAN_noclust"]

# 2. 331 reproduction of the grouping with the byte-identical copied driver
orig_cmp = json.load(open(M09 + "/grouping/COMPARE_CAPSCAN_noclust.json"))
rep_dir = (M09 if self_test else A + "/run331") + "/grouping"
rep_cmp = json.load(open(rep_dir + "/COMPARE_CAPSCAN_noclust.json"))
if neg: rep_cmp["comparison"]["inferred_study_groups"]["fresh7"] += 1
check(rep_cmp["comparison"] == orig_cmp["comparison"], "331 reproduction COMPARE differs from original M09")
check(rep_cmp["comparison"]["inferred_study_groups"]["fresh7"] == pk["fresh7"] and rep_cmp["comparison"]["inferred_study_groups"]["baseline"] == pk["baseline"],
      "331 reproduction differs from packaged Source Data sheet g")
file_eq = {}
for fn in ("study_groups.tsv", "study_group_matrix.tsv", "peptide_to_study_group.tsv", "study_group_long.tsv", "acquisition_summary.tsv"):
    a, b = f"{M09}/grouping/out/CAPSCAN_noclust/{fn}", f"{rep_dir}/out/CAPSCAN_noclust/{fn}"
    file_eq[fn] = sha(a) == sha(b); check(file_eq[fn], "331 reproduction file differs: " + fn)

if self_test:
    print("SELFTEST", "FAIL" if fails else "PASS", fails); sys.exit(1 if fails else 0)

# 3. the 16 new searches and the 347 grouping
new_rec = [json.load(open(p)) for p in sorted(glob.glob(A + "/run/receipts/receipt_03[3-4][0-9].json")) if 331 <= int(p[-9:-5]) <= 346]
check(len(new_rec) == 16 and all(r["pass"] is True for r in new_rec), "new search receipts not 16 PASS")
orig_receipts = sorted(glob.glob(M09 + "/receipts/receipt_*.json"))
orig_pass = sum(json.load(open(p))["pass"] is True for p in orig_receipts)
check(len(orig_receipts) == 1990 and orig_pass == 1990, "original 1,990 receipts")
cmp347 = json.load(open(A + "/run/grouping/COMPARE_CAPSCAN_noclust.json"))["comparison"]
check(cmp347["acquisitions"] == {"baseline": 347, "fresh7": 347}, "347 acquisitions")
b7, f7 = cmp347["inferred_study_groups"]["baseline"], cmp347["inferred_study_groups"]["fresh7"]
chg347 = 100 * (f7 / b7 - 1)
status = "FAIL" if fails else "PASS"

L = []
L.append("# A2: minimum-length-7 sensitivity (Supplementary Fig. 4), 331 vs 347\n")
L.append(f"Status: {status}{(' ' + '; '.join(fails)) if fails else ''}\n")
L.append("Only the CapScan FastaLake production arm (FL_cohortMG, assembly-inclusive) gains the 16 participant-15 searches. "
         "The public-sequence arm (FL_public, A08) has no FASTA for runs 1799-1814 and the metagenome arm (MG) has no matched database, "
         "so the matched 282-acquisition benchmark and every FL_public/MG value cannot change.\n")
L.append("## Reproduction (331)\n")
L.append(f"- Copied group_fresh.py (byte-identical) on the original 331 fresh-7 searches: COMPARE equal to the original M09 run: {rep_cmp['comparison'] == orig_cmp['comparison']}; "
         f"output tables identical (sha256): {file_eq}.")
L.append(f"- Packaged Source Data sheet g, CAPSCAN_noclust: baseline {int(pk['baseline'])}, fresh7 {int(pk['fresh7'])}, change {pk['change_pct']:.4f}%; reproduced: "
         f"{rep_cmp['comparison']['inferred_study_groups']['baseline']} / {rep_cmp['comparison']['inferred_study_groups']['fresh7']}.")
L.append(f"- 16 new searches: {sum(r['pass'] is True for r in new_rec)}/16 receipts PASS (min_len 7 in resolved results.json, settings equal to source except min_len, Sage 0.14.6, no accepted peptide <7, LFQ as source).\n")
L.append("## Printed numbers\n")
L.append("| manuscript location | printed (331) | 347 | basis |")
L.append("|---|---|---|---|")
L.append(f"| Discussion, Methods, Supp. Fig. 4 caption: production searches repeated | 1,990 (1,985 standard + 5 LFQ-enabled PREDICT) | 2,006 (2,001 standard + 5) | 1,990 original receipts PASS ({orig_pass}) + 16 new CapScan FL_cohortMG PASS |")
L.append(f"| Discussion: median same-spectrum gain, PREDICT, 5-50 -> fresh >=7 | 18.8% -> 14.5% | unchanged by construction | 296 matched PREDICT; recomputed from package: {med('PREDICT','baseline5','gain_pct'):.2f} -> {med('PREDICT','fresh7','gain_pct'):.2f} |")
L.append(f"| Discussion: same, CapScan | 51.7% -> 45.0% | unchanged by construction | 282 matched FL_public vs MG; no MG or A08 for participant 15; package: {med('CAPSCAN','baseline5','gain_pct'):.2f} -> {med('CAPSCAN','fresh7','gain_pct'):.2f} |")
L.append(f"| Discussion: MG identities recovered, fresh >=7 | 85.4% and 84.0% | unchanged by construction | package: {med('PREDICT','fresh7','comparator_identity_recovery_pct'):.2f} / {med('CAPSCAN','fresh7','comparator_identity_recovery_pct'):.2f} |")
L.append(f"| Discussion: PREDICT gain at >=9 | 9.4% | unchanged (PREDICT) | package: {med('PREDICT','nested9','gain_pct'):.2f} |")
fr = e.set_index(["seed", "population"]).pooled_diagnostic_pct
L.append(f"| Discussion: CAMPI fusion diagnostic, >=7 -> >=9 | 2.4-2.6% -> 1.6% (5-50: 2.59-2.79; nested >=9 1.55-1.56) | unchanged (CAMPI) | package fresh7 {fr.xs('fresh7', level=1).min():.2f}-{fr.xs('fresh7', level=1).max():.2f}, nested9 {fr.xs('nested9', level=1).min():.2f}-{fr.xs('nested9', level=1).max():.2f} |")
cs = c.set_index("arm").loc["CapScan FL"]
L.append(f"| Supp. Fig. 4c, CapScan FL bar (arm FL_public, n=331) | filtered {cs['filtered_to_7']:.2f}%, new search {cs['new_search_7']:.2f}% | unchanged: no A08 FASTA for participant 15 | FL_public cannot gain the 16 |")
L.append(f"| Supp. Fig. 4a,b panel n | 296 PREDICT, 282 CapScan | unchanged | matched rosters |")
L.append(f"| Supp. Fig. 4g, CapScan study groups, fresh >=7 vs 5-50 | CapScan (331): {int(pk['baseline']):,} -> {int(pk['fresh7']):,} ({pk['change_pct']:+.2f}%, bar label +0.3) | CapScan (347): {b7:,} -> {f7:,} ({chg347:+.2f}%, bar label {chg347:+.1f}) | baseline 347 = C2 grouping (Sept-10 driver, roster 347); fresh7 347 = copied group_fresh.py with C2 manifest |")
L.append(f"| Discussion: protein-group counts essentially unchanged | yes (+0.3% CapScan) | {'yes' if abs(chg347) < 1 else 'CHECK'} ({chg347:+.2f}%) | as above |")
L.append("")
L.append("Other CapScan grouping quantities (not printed), 347 baseline -> fresh7: " + "; ".join(f"{k} {v['baseline']} -> {v['fresh7']}" for k, v in cmp347.items()))
L.append("\nUnchanged and not involving CapScan: Supp. Fig. 4d (CAMPI members), 4e (CAMPI fusion), 4f (healthy stool), 4g healthy/PREDICT/CAMPI bars.")
open(A + "/RESULTS.md", "w").write("\n".join(L) + "\n")
print(status); sys.exit(1 if fails else 0)
