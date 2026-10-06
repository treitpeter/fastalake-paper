#!/usr/bin/env python3
"""Peptide-centric equal-split projection, every variant the revised framework needs.

Locked rule (2026-09-19). For every accepted quantified LFQ feature, take its canonical
peptide, find the compatible proteins in that cohort's study search space by exact
substring containment on I/L-normalised sequence, collect the distinct annotations those
proteins carry at the level, divide the feature intensity equally among them, and retain
features with no qualifying annotation explicitly as FunctionallyUnassigned.

Vocabularies are copied verbatim from
MM_review_analysis/A19_4_6_biology_robustness/frozen/common.py:annotation_terms.

Protein sets:
  all              every compatible protein
  microbial_single_og   compatible proteins with exactly one bacterial/archaeal OG token.
                   Peptide-level analogue of frozen common.aggregate(microbial_domain=True),
                   which kept groups with len(annotation_terms(OG,'OG')) == 1.
  microbial_broad  compatible proteins with an MGYG/GMGC/GMSC accession or at least one
                   bacterial/archaeal OG token, excluding '_HUMAN'. Peptide-level analogue
                   of the inline rule in frozen derive_capscan.py and derive_clinical.py.

Species labels stay representative-based: a feature's allocated share carries the
representative taxon of the study group the September direct assignment gave that feature.
'_HUMAN' study groups label as Host and 'Unassigned' as 'No species annotation', exactly
as frozen derive_capscan.py does. No LCA and no new taxonomic inference.

Outputs per cohort:
  <C>_<level>_<set>.tsv.gz      function x sample allocation, one file per level and set
  <C>_KO_species.npz            KO x species x sample under the simplified (direct)
                                grouping's study-group representative taxa
  <C>_KO_species_original.npz   the same allocation under the original grouping's taxa,
                                reconstructed from the run's own changed_feature_assignments
  <C>_variant_stats.json        conservation and unassigned accounting for every output
"""
import ahocorasick, gzip, json, os, re, sys
from collections import defaultdict
import numpy as np

BASE = "LAB_ROOT"
GROUPING = os.environ.get("FL_GROUPING_COHORTS",
    BASE + "/simplified_grouping_results/manuscript_cohorts_2026-09-10_2200/cohorts")
COH, OUT = sys.argv[1], sys.argv[2]
RUN = f"{GROUPING}/{COH}_noclust"
EGG = f"{BASE}/canonical/publication_data/arms/noclust_2026-08-07/cohorts/{COH}/11_eggnog/{COH}_eggnog_accession_keyed.tsv"
norm = lambda s: s.replace("I", "L")
say = lambda *a: print(*a, flush=True)

LEVELS = ['KO', 'OG', 'Pathway']
COL = {'KO': 'KEGG_ko', 'OG': 'OG', 'Pathway': 'KEGG_Pathway'}
SETS = ['all', 'microbial_single_og', 'microbial_broad']
UNASSIGNED = 'FunctionallyUnassigned'


def terms(value, level):                       # verbatim from frozen/common.py
    if level == 'KO':
        return sorted(set(re.findall(r'K\d{5}', value)))
    if level == 'OG':
        return sorted(set(re.findall(r'[^;,]+@(?:2\|Bacteria|2157\|Archaea)', value)))
    if level == 'Pathway':
        return sorted(t for t in set(re.findall(r'map\d{5}', value)) if not t.startswith('map011'))
    raise ValueError(level)


# ---- annotations, per accession, with the two microbial-domain memberships -------------
ann = {}
n_single_og = n_broad = 0
with open(EGG) as fh:
    ix = {c: i for i, c in enumerate(next(fh).rstrip("\n").split("\t"))}
    for line in fh:
        r = line.rstrip("\n").split("\t")
        acc = r[ix["accession"]]
        raw = {}
        for L in LEVELS:
            v = r[ix[COL[L]]] if ix[COL[L]] < len(r) else ""
            raw[L] = terms(v, L) if v and v != "-" else []
        og = raw['OG']
        single = len(og) == 1
        broad = (acc.startswith(('MGYG', 'GMGC', 'GMSC')) or bool(og)) and '_HUMAN' not in acc
        n_single_og += single
        n_broad += broad
        rec = {}
        for L in LEVELS:
            if raw[L]:
                rec[L] = frozenset(raw[L])
        if rec or single or broad:
            ann[acc] = (rec, single, broad)
say(f"{COH}: annotated accessions {len(ann):,}  single_og {n_single_og:,}  broad {n_broad:,}")

# ---- the study search space: every FASTA the direct assignment searched ----------------
assign = set()
with open(f"{RUN}/direct/peptide_to_study_group.tsv") as fh:
    next(fh)
    for line in fh:
        p, _ = line.rstrip("\n").split("\t")
        assign.add(norm(p))
fastas = set()
for e in json.load(open(f"{RUN}/direct/input_manifest.json")):
    if not e["path"].endswith("results.json"):
        continue
    def walk(o):
        if isinstance(o, dict):
            for k, v in o.items():
                if isinstance(v, str) and v.endswith(".fasta"):
                    fastas.add(v)
                else:
                    walk(v)
        elif isinstance(o, list):
            for x in o:
                walk(x)
    walk(json.load(open(e["path"])))
seq = {}
for f in sorted(fastas):
    acc, buf = None, []
    with open(f) as fh:
        for line in fh:
            if line[0] == ">":
                if acc is not None and acc not in seq:
                    seq[acc] = norm("".join(buf))
                acc = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if acc is not None and acc not in seq:
        seq[acc] = norm("".join(buf))
say(f"  universe {len(seq):,} proteins, {len(assign):,} assigned peptides")

# ---- containment scan: peptide -> annotation sets, per level and per protein set -------
A = ahocorasick.Automaton()
for p in assign:
    A.add_word(p, p)
A.make_automaton()
pep = {(L, S): defaultdict(set) for L in LEVELS for S in SETS}
for i, (acc, s) in enumerate(seq.items()):
    rec = ann.get(acc)
    if not rec:
        continue
    levels, single, broad = rec
    if not levels:
        continue
    hits = {p for _, p in A.iter(s)}
    if not hits:
        continue
    member = ['all'] + (['microbial_single_og'] if single else []) + (['microbial_broad'] if broad else [])
    for L, vals in levels.items():
        for S in member:
            d = pep[L, S]
            for p in hits:
                d[p] |= vals
    if i % 100000 == 0:
        say(f"    scanned {i:,}/{len(seq):,}")
say("  scan done")

# ---- study-group representative taxon, for the species dimension ----------------------
taxon = {}
with gzip.open(f"{RUN}/reaggregated/representative_annotations.tsv.gz", "rt") as fh:
    ix = {c: i for i, c in enumerate(next(fh).rstrip("\n").split("\t"))}
    for line in fh:
        r = line.rstrip("\n").split("\t")
        taxon[r[ix["study_group"]]] = r[ix["representative_taxon"]]


def species_of(group):
    if '_HUMAN' in group:
        return 'Host'
    t = taxon.get(group, 'Unassigned')
    return 'No species annotation' if t == 'Unassigned' else t


# ---- the original grouping's assignment for the same features -------------------------
# The run records only the features whose group changed, so the original group is the
# recorded original_group where present and the direct group everywhere else.
changed = {}
intern = {}
with gzip.open(f"{RUN}/reaggregated/changed_feature_assignments.tsv.gz", "rt") as fh:
    cx = {c: i for i, c in enumerate(next(fh).rstrip("\n").split("\t"))}
    for line in fh:
        r = line.rstrip("\n").split("\t")
        g = r[cx["original_group"]]
        changed[r[cx["sample"]] + "\t" + r[cx["input_line"]]] = intern.setdefault(g, g)
say(f"  original-rule reassignments {len(changed):,}")


# ---- allocation over the accepted quantified features ---------------------------------
prof = {(L, S): defaultdict(lambda: defaultdict(float)) for L in LEVELS for S in SETS}
joint = defaultdict(lambda: defaultdict(float))       # (KO, species) -> sample -> value
joint_original = defaultdict(lambda: defaultdict(float))
samples = set()
total = 0.0
rows = 0
relabelled = 0.0
unknown_peptides = set()
with gzip.open(f"{RUN}/direct/feature_assignments.tsv.gz", "rt") as fh:
    ix = {c: i for i, c in enumerate(next(fh).rstrip("\n").split("\t"))}
    for line in fh:
        r = line.rstrip("\n").split("\t")
        s = r[ix["sample"]]
        p = norm(r[ix["canonical_peptide"]])
        v = float(r[ix["intensity"]])
        g = r[ix["study_group"]]
        samples.add(s)
        total += v
        rows += 1
        if p not in assign:
            unknown_peptides.add(p)
        for L in LEVELS:
            for S in SETS:
                t = pep[L, S].get(p)
                d = prof[L, S]
                if t:
                    w = v / len(t)
                    for x in t:
                        d[x][s] += w
                else:
                    d[UNASSIGNED][s] += v
        t = pep['KO', 'all'].get(p)
        sp = species_of(g)
        go = changed.get(s + "\t" + r[ix["input_line"]], g)
        spo = sp if go == g else species_of(go)
        if spo != sp:
            relabelled += v
        if t:
            w = v / len(t)
            for k in t:
                joint[k, sp][s] += w
                joint_original[k, spo][s] += w
        else:
            joint[UNASSIGNED, sp][s] += v
            joint_original[UNASSIGNED, spo][s] += v
say(f"  features {rows:,} samples {len(samples)} total {total:.6e}")
# The automaton is built from the direct peptide table. A canonical peptide present in the
# feature ledger but absent there would be absorbed into the unassigned row without
# comment, inflating the framework's own headline accounting quantity.
assert not unknown_peptides, (f"{len(unknown_peptides):,} canonical peptides in the feature "
                              f"ledger are absent from the direct peptide table, e.g. "
                              f"{sorted(unknown_peptides)[:3]}")
say(f"  intensity relabelled between grouping rules {100*relabelled/total:.2f}%")

os.makedirs(OUT, exist_ok=True)
names = sorted(samples)
stats = {'cohort': COH, 'samples': len(names), 'features_read': rows,
         'total_intensity': total, 'variants': {}}
for L in LEVELS:
    for S in SETS:
        d = prof[L, S]
        path = f"{OUT}/{COH}_{L}_{S}.tsv.gz"
        alloc = 0.0
        with gzip.open(path, "wt") as w:
            w.write("feature\t" + "\t".join(names) + "\n")
            for feat in sorted(d):
                row = d[feat]
                alloc += sum(row.values())
                w.write(feat + "\t" + "\t".join(f"{row.get(s, 0.0):.10g}" for s in names) + "\n")
        un = sum(d[UNASSIGNED].values())
        stats['variants'][f"{L}|{S}"] = dict(features=len(d), allocated=alloc,
                                             conservation_ratio=alloc / total,
                                             unassigned_intensity=un,
                                             unassigned_pct=100 * un / total)
        say(f"  {L:8s} {S:20s} {len(d):7,} rows  conservation {alloc/total:.10f}  unassigned {100*un/total:6.2f}%")
        assert abs(alloc / total - 1.0) < 1e-9, f"{L}/{S} does not conserve intensity"

sindex = {s: i for i, s in enumerate(names)}


def save_joint(table, name, label):
    keys = sorted(table)
    say(f"  {label} pairs {len(keys):,}  array {len(keys)*len(names)*8/2**30:.2f} GiB")
    values = np.zeros((len(keys), len(names)), dtype=float)
    for i, k in enumerate(keys):
        for s, v in table[k].items():
            values[i, sindex[s]] = v
    alloc = float(values.sum())
    un = float(sum(v for (k, sp), row in table.items() if k == UNASSIGNED
                   for v in row.values()))
    np.savez_compressed(f"{OUT}/{COH}_{name}.npz",
                        ko=np.array([k for k, _ in keys], dtype=str),
                        species=np.array([sp for _, sp in keys], dtype=str),
                        values=values, samples=np.array(names, dtype=str))
    stats[name] = dict(pairs=len(keys), allocated=alloc, conservation_ratio=alloc / total,
                       unassigned_intensity=un, unassigned_pct=100 * un / total,
                       distinct_species=len({sp for _, sp in keys}),
                       distinct_KO=len({k for k, _ in keys} - {UNASSIGNED}))
    say(f"  {label} conservation {alloc/total:.10f}  unassigned {100*un/total:.2f}%  "
        f"species {stats[name]['distinct_species']:,}")
    assert abs(alloc / total - 1.0) < 1e-9, f"{label} does not conserve intensity"
    del values


stats['species_relabelled_pct_between_grouping_rules'] = 100 * relabelled / total
save_joint(joint, 'KO_species', 'joint KO x species (simplified)')
save_joint(joint_original, 'KO_species_original', 'joint KO x species (original)')
json.dump(stats, open(f"{OUT}/{COH}_variant_stats.json", "w"), indent=2)
say("BUILD DONE")
