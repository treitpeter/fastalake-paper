"""Addendum 1: five LFQ-on repaired PREDICT FL_cohortMG searches (indices 1985-1989), appended; existing rows untouched."""
import json, os, copy, math, hashlib, glob, csv
F = "LAB_ROOT"
RUN = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = F + "/canonical/publication_data/arms/M09_minlen7_2026-09-24"
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()
rows = list(csv.DictReader(open(f"{RUN}/tasks.tsv"), delimiter="\t"))
assert len(rows) == 1985, len(rows)
cols = list(rows[0])
lfqoff = [r for r in rows if r["lfq"] == "False"]
assert len(lfqoff) == 5
new = []
for i, r in enumerate(sorted(lfqoff, key=lambda r: r["sample"])):
    idx = 1985 + i; a = r["sample"]
    src_path = glob.glob(f"{F}/MM_review_analysis/recovery_2026-09-15/quantification_recovery/{a}/results/{a}__tryptic_lfq/search/results.json")
    assert len(src_path) == 1; src_path = src_path[0]
    src = json.load(open(src_path)); off = json.load(open(r["source_config"]))
    a1, b1 = copy.deepcopy(src), copy.deepcopy(off)
    for d in (a1, b1):
        d.pop("output_paths"); d["quant"].pop("lfq")
    assert a1 == b1, a                       # identical to the inventory config except LFQ
    assert src["quant"]["lfq"] is True and src["database"]["enzyme"]["min_len"] == 5
    cfg = copy.deepcopy(src); del cfg["version"], cfg["output_paths"]
    cfg["database"]["enzyme"]["min_len"] = 7
    od = f"{OUT}/PREDICT/FL_cohortMG_lfq/{a}"; assert not os.path.exists(od)
    cfg["output_directory"] = od
    cp = f"{RUN}/configs/{idx:04d}_PREDICT_FL_cohortMG_lfq_{a}.json"; json.dump(cfg, open(cp, "w"), indent=2)
    fmb = os.path.getsize(src["database"]["fasta"]) / 1e6
    new.append(dict(index=idx, cohort="PREDICT", arm="FL_cohortMG_lfq", sample=a, lfq=True, fasta_mb=round(fmb, 2),
                    mem_gb=max(16, math.ceil(fmb * 0.361 * 1.2 + 6)), config=cp, source_config=src_path,
                    source_sha256=sha(src_path), baseline_dir=os.path.dirname(src_path), out_dir=od))
with open(f"{RUN}/tasks.tsv", "a") as fh:
    for t in new: fh.write("\t".join(str(t[c]) for c in cols) + "\n")
open(f"{RUN}/addendum1_indices.txt", "w").write("\n".join(str(t["index"]) for t in new) + "\n")
json.dump(dict(created=__import__("datetime").datetime.now().isoformat(), reason="the five repaired PREDICT FL_cohortMG inventory configs have LFQ off; the repaired grouping (biological_grouping_recovery) used LFQ-on quantification_recovery searches of the same FASTA/mzML", tasks=new, prepare_sha256=sha(os.path.abspath(__file__))), open(f"{RUN}/ADDENDUM1.json", "w"), indent=1)
print(len(new), [t["index"] for t in new])
