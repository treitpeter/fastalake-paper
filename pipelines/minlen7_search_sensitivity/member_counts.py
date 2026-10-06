"""CAMPI member-specific peptide counts for baseline 5-50, nested >=7, fresh >=7, nested >=9 and fresh7_nested9.
Membership rule identical to A14 frozen/search.py (exact ' OS=<member> ' header match; non-members contribute no label;
specific = exactly one member; I/L-normalised exact containment). Pooled = union over 15 acquisitions; per-acquisition counts incl. zeros."""
import json, csv, os, gzip, collections, ahocorasick
RUN = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
A14 = "LAB_ROOT/MM_review_analysis/A14_4_3_campi_gap"
m = json.load(open(A14 + "/inputs/manifest.json"))
members = list(csv.DictReader(open(m["members"]["path"]), delimiter="\t"))
tasks = [r for r in csv.DictReader(open(RUN + "/tasks.tsv"), delimiter="\t") if r["cohort"] == "CAMPI"]
POPS = ["baseline5", "nested7", "fresh7", "nested9", "fresh7_nested9"]
sets = {}
for t in tasks:
    p = f"{RUN}/wave2/peptide_tasks/{int(t['index']):04d}"
    assert json.load(open(p + "/COMPLETE.json"))["status"] == "PASS"
    for pop in POPS:
        sets[(t["arm"], t["sample"], pop)] = set(gzip.open(f"{p}/{pop}.txt.gz", "rt").read().split())
union = set().union(*sets.values())
A = ahocorasick.Automaton()
for q in union: A.add_word(q, q)
A.make_automaton()
def fasta(p):
    h, s = None, []
    for line in open(p):
        line = line.rstrip("\n")
        if line.startswith(">"):
            if h: yield h, "".join(s)
            h, s = line[1:], []
        else: s.append(line.strip())
    if h: yield h, "".join(s)
mem = collections.defaultdict(set)
for h, s in fasta(m["reference_fasta"]["path"]):
    lab = [r["member"] for r in members if " OS=" + r["campi_os_string"] + " " in h + " "]
    if lab:
        for _, q in A.iter(s.replace("I", "L")): mem[q].add(lab[0])
names = [r["member"] for r in members]
rows = []
for arm in sorted({t["arm"] for t in tasks}):
    samples = sorted({t["sample"] for t in tasks if t["arm"] == arm})
    for pop in POPS:
        pooled = set().union(*(sets[(arm, s, pop)] for s in samples))
        for name in names:
            per = [sum(1 for q in sets[(arm, s, pop)] if mem[q] == {name}) for s in samples]
            rows.append(dict(arm=arm, population=pop, member=name,
                             pooled_specific=sum(1 for q in pooled if mem[q] == {name}),
                             acquisitions_ge1=sum(x > 0 for x in per), median_per_acq=sorted(per)[len(per)//2], per_acquisition=";".join(map(str, per))))
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "MEMBER_COUNTS_BY_LENGTH.tsv"), "w", newline="") as f:
    w = csv.DictWriter(f, list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)
for r in rows:
    if r["member"] in ("Lactobacillus plantarum", "Clostridium butyricum", "Bifidobacterium longum") and r["arm"] in ("FL", "REF", "MG"):
        print(r["arm"], r["population"], r["member"], r["pooled_specific"], r["acquisitions_ge1"], r["median_per_acq"])
