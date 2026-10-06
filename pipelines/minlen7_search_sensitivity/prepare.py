"""Build isolated min_len=7 Sage configs for M09 wave 1. Writes configs/, tasks.tsv, PREPARE.json."""
import json, os, copy, math, hashlib, collections, datetime
F = "LAB_ROOT"
RUN = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MAN = F + "/00_PACKAGE_2026-09-17/cleanup_2026-09-22/v10_mm_philip_comment_audit_2026-09-22/nine_point_followup/manuscript_final_pass_2026-09-23/mm_v11_revision_2026-09-23/reappraisal_2026-09-24/main_methods_pass_2026-09-24/reviews/design_selection_search/CONFIG_MANIFEST.json"
ROSTER = F + "/simplified_grouping_results/manuscript_cohorts_2026-09-10_2200/rosters/HSTOOL.txt"
OUT = F + "/canonical/publication_data/arms/M09_minlen7_2026-09-24"
NEW_MIN = 7

def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()

rows = json.load(open(MAN))
roster = set(open(ROSTER).read().split())
assert len(rows) == 2024 and len(roster) == 58
kept = [r for r in rows if not (r["cohort"] == "HSTOOL" and r["sample"] not in roster)]
assert len(kept) == 1985, len(kept)
kept.sort(key=lambda r: (r["cohort"], r["arm"], r["sample"]))
tasks = []
for i, r in enumerate(kept):
    assert sha(r["path"]) == r["sha256"], r["path"]          # source config unchanged since audit
    src = json.load(open(r["path"]))
    assert src["version"] == "0.14.6"
    assert src["database"]["enzyme"]["min_len"] == 5 and src["database"]["enzyme"]["max_len"] == 50
    new = copy.deepcopy(src)
    del new["version"], new["output_paths"]
    new["database"]["enzyme"]["min_len"] = NEW_MIN
    od = f"{OUT}/{r['cohort']}/{r['arm']}/{r['sample']}"
    assert not os.path.exists(od), od                       # never write into an existing folder
    new["output_directory"] = od
    # assert the only differences are the declared ones
    a = copy.deepcopy(src); b = copy.deepcopy(new)
    for d in (a, b):
        d.pop("version", None); d.pop("output_paths", None); d.pop("output_directory", None)
        d["database"]["enzyme"].pop("min_len")
    assert a == b
    fasta = src["database"]["fasta"]; mz = src["mzml_paths"]
    assert len(mz) == 1 and os.path.exists(fasta) and os.path.exists(mz[0])
    fmb = os.path.getsize(fasta) / 1e6
    mem = max(16, math.ceil(fmb * 0.361 * 1.2 + 6))          # measured RAM law + 20%, 16G floor as production
    base_dir = os.path.dirname(r["path"])
    assert os.path.exists(base_dir + "/results.sage.tsv")
    cfg = f"{RUN}/configs/{i:04d}_{r['cohort']}_{r['arm']}_{r['sample']}.json"
    json.dump(new, open(cfg, "w"), indent=2)
    tasks.append(dict(index=i, cohort=r["cohort"], arm=r["arm"], sample=r["sample"], lfq=src["quant"]["lfq"],
                      fasta_mb=round(fmb, 2), mem_gb=mem, config=cfg, source_config=r["path"],
                      source_sha256=r["sha256"], baseline_dir=base_dir, out_dir=od))
with open(f"{RUN}/tasks.tsv", "w") as fh:
    cols = list(tasks[0]); fh.write("\t".join(cols) + "\n")
    for t in tasks: fh.write("\t".join(str(t[c]) for c in cols) + "\n")
summ = collections.Counter((t["cohort"], t["arm"]) for t in tasks)
json.dump(dict(created=datetime.datetime.now().isoformat(), n_tasks=len(tasks), new_min_len=NEW_MIN,
               output_root=OUT, source_manifest=MAN, source_manifest_sha256=sha(MAN), roster=ROSTER,
               roster_sha256=sha(ROSTER), per_route={f"{k[0]}/{k[1]}": v for k, v in sorted(summ.items())},
               lfq_off=[t["index"] for t in tasks if not t["lfq"]], max_mem_gb=max(t["mem_gb"] for t in tasks),
               prepare_sha256=sha(os.path.abspath(__file__))), open(f"{RUN}/PREPARE.json", "w"), indent=1)
print(len(tasks), dict(summ), "max mem", max(t["mem_gb"] for t in tasks))
