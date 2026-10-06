#!/usr/bin/env python3
"""The primary species dimension: carrier-supported, under both grouping rules.

Measured 2026-09-19: 13.8% of species-labelled allocation in CapScan, 14.2% in PREDICT and
10.6% in healthy stool carried the representative taxon of the feature's own study group
while every compatible protein that carried the KO resolved to some other species. A
further 9.9 to 12.7% carried a label no compatible carrier could speak to at all.

This builder keeps the representative label only where a compatible carrier of that KO
resolves to it, and routes everything else to 'No species annotation', the unresolved
category the panels already exclude from their resolved fraction. Absent evidence and
contradicting evidence are therefore treated alike, which is conservative rather than
correct, and the split between them is reported in audit/measurements/.

No LCA is introduced. Labels remain representative taxa. Intensity conserves.

Both grouping rules are emitted from one pass, because a comparison presented as differing
only in the grouping rule has to use the same annotation rule on both sides. Species
resolution uses the published rule: one UHGG species from an MGYG genome named in the
accession, accepted only when exactly one such species is found.
"""
import ahocorasick, csv, gzip, json, os, re, sys
from collections import defaultdict
import numpy as np

BASE = "LAB_ROOT"
ARM = f"{BASE}/canonical/publication_data/arms/noclust_2026-08-07"
GENOMES = "LAB_ROOT/.._Microbpred/metalab_magnify_databases/human_gut/genomes-all_metadata.tsv"
COH, OUT = sys.argv[1], sys.argv[2]
GROUPING = os.environ.get("FL_GROUPING_COHORTS",
    BASE + "/simplified_grouping_results/manuscript_cohorts_2026-09-10_2200/cohorts")
RUN = f"{GROUPING}/{COH}_noclust"
EGG = f"{ARM}/cohorts/{COH}/11_eggnog/{COH}_eggnog_accession_keyed.tsv"
norm = lambda s: s.replace("I", "L")
say = lambda *a: print(*a, flush=True)
csv.field_size_limit(1 << 30)
MGYG = re.compile(r'MGYG\d+')
UNASSIGNED = 'FunctionallyUnassigned'
UNRESOLVED = 'No species annotation'


def kos(value):                                # verbatim from frozen common.annotation_terms
    return sorted(set(re.findall(r'K\d{5}', value)))


g2s = {}
with open(GENOMES) as fh:
    r = csv.reader(fh, delimiter='\t')
    h = next(r)
    iG, iRep, iLin = h.index('Genome'), h.index('Species_rep'), h.index('Lineage')
    rows, rep2lin = [], {}
    for row in r:
        rows.append((row[iG], row[iRep], row[iLin]))
        if row[iLin]:
            rep2lin[row[iRep]] = row[iLin]
    for gg, rep, lin in rows:
        for t in (lin or rep2lin.get(rep, '')).split(';'):
            if t.startswith('s__') and len(t) > 3:
                g2s[gg] = t[3:]
sp_id = {}


def species_id_of_accession(acc):
    names = {g2s[g] for g in MGYG.findall(acc) if g in g2s}
    if len(names) != 1:
        return None
    return sp_id.setdefault(next(iter(names)), len(sp_id))


ann = {}
with open(EGG) as fh:
    ix = {c: i for i, c in enumerate(next(fh).rstrip("\n").split("\t"))}
    for line in fh:
        r = line.rstrip("\n").split("\t")
        v = r[ix['KEGG_ko']] if ix['KEGG_ko'] < len(r) else ''
        t = kos(v) if v and v != '-' else []
        if t:
            ann[r[ix['accession']]] = t
say(f"{COH}: accessions with a KO {len(ann):,}")

assign = set()
with open(f"{RUN}/direct/peptide_to_study_group.tsv") as fh:
    next(fh)
    for line in fh:
        assign.add(norm(line.split("\t")[0]))
fastas = set()
for e in json.load(open(f"{RUN}/direct/input_manifest.json")):
    if not e['path'].endswith('results.json'):
        continue
    def walk(o):
        if isinstance(o, dict):
            for k, v in o.items():
                if isinstance(v, str) and v.endswith('.fasta'):
                    fastas.add(v)
                else:
                    walk(v)
        elif isinstance(o, list):
            for x in o:
                walk(x)
    walk(json.load(open(e['path'])))
seq = {}
for f in sorted(fastas):
    acc, buf = None, []
    with open(f) as fh:
        for line in fh:
            if line[0] == '>':
                if acc is not None and acc not in seq:
                    seq[acc] = norm(''.join(buf))
                acc = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if acc is not None and acc not in seq:
        seq[acc] = norm(''.join(buf))
say(f"universe {len(seq):,} proteins, {len(assign):,} assigned peptides")

A = ahocorasick.Automaton()
for p in assign:
    A.add_word(p, p)
A.make_automaton()
carriers = defaultdict(dict)          # peptide -> {ko: set of carrier species ids}
for i, (acc, s) in enumerate(seq.items()):
    t = ann.get(acc)
    if not t:
        continue
    hits = {p for _, p in A.iter(s)}
    if not hits:
        continue
    sid = species_id_of_accession(acc)
    for p in hits:
        d = carriers[p]
        for k in t:
            e = d.get(k)
            if e is None:
                e = d[k] = set()
            if sid is not None:
                e.add(sid)
    if i % 100000 == 0:
        say(f"  scanned {i:,}/{len(seq):,}")
say(f"scan done; peptides with a KO {len(carriers):,}")

taxon = {}
with gzip.open(f"{RUN}/reaggregated/representative_annotations.tsv.gz", 'rt') as fh:
    ix = {c: i for i, c in enumerate(next(fh).rstrip("\n").split("\t"))}
    for line in fh:
        r = line.rstrip("\n").split("\t")
        taxon[r[ix['study_group']]] = r[ix['representative_taxon']]


def label_of(group):
    if '_HUMAN' in group:
        return 'Host'
    t = taxon.get(group, 'Unassigned')
    return UNRESOLVED if t == 'Unassigned' else t


changed = {}
intern = {}
with gzip.open(f"{RUN}/reaggregated/changed_feature_assignments.tsv.gz", 'rt') as fh:
    cx = {c: i for i, c in enumerate(next(fh).rstrip("\n").split("\t"))}
    for line in fh:
        r = line.rstrip("\n").split("\t")
        g = r[cx['original_group']]
        changed[r[cx['sample']] + "\t" + r[cx['input_line']]] = intern.setdefault(g, g)
say(f"original-rule reassignments {len(changed):,}")

joint = {'simplified': defaultdict(lambda: defaultdict(float)),
         'original': defaultdict(lambda: defaultdict(float))}
tally = {'simplified': defaultdict(float), 'original': defaultdict(float)}
sample_ix = {}
tot = 0.0
rows = 0
unknown = set()
with gzip.open(f"{RUN}/direct/feature_assignments.tsv.gz", 'rt') as fh:
    ix = {c: i for i, c in enumerate(next(fh).rstrip("\n").split("\t"))}
    for line in fh:
        r = line.rstrip("\n").split("\t")
        v = float(r[ix['intensity']])
        g = r[ix['study_group']]
        p = norm(r[ix['canonical_peptide']])
        si = sample_ix.setdefault(r[ix['sample']], len(sample_ix))
        rows += 1
        tot += v
        if p not in assign:
            unknown.add(p)
        d = carriers.get(p)
        go = changed.get(r[ix['sample']] + "\t" + r[ix['input_line']], g)
        for rule, grp in (('simplified', g), ('original', go)):
            lab = label_of(grp)
            if not d:
                joint[rule][UNASSIGNED, lab][si] += v
                tally[rule]['unassigned'] += v
                continue
            w = v / len(d)
            for k, carrier_species in d.items():
                if lab in ('Host', UNRESOLVED):
                    tag, kind = lab, 'label is not a species claim'
                elif carrier_species and sp_id.get(lab) in carrier_species:
                    tag, kind = lab, 'supported'
                elif carrier_species:
                    tag, kind = UNRESOLVED, 'contradicted, withheld'
                else:
                    tag, kind = UNRESOLVED, 'undeterminable, withheld'
                joint[rule][k, tag][si] += w
                tally[rule][kind] += w
say(f"features {rows:,} total {tot:.6e}")
assert not unknown, (f"{len(unknown):,} canonical peptides in the feature ledger are absent "
                     f"from the direct peptide table, e.g. {sorted(unknown)[:3]}")

os.makedirs(OUT, exist_ok=True)
names = sorted(sample_ix, key=lambda s: sample_ix[s])
stats = {'cohort': COH, 'features': rows, 'total_intensity': tot, 'rules': {}}
for rule, stem in (('simplified', 'KO_species'), ('original', 'KO_species_original')):
    table = joint[rule]
    keys = sorted(table)
    values = np.zeros((len(keys), len(names)))
    for i, k in enumerate(keys):
        for si, x in table[k].items():
            values[i, si] = x
    alloc = float(values.sum())
    np.savez_compressed(f"{OUT}/{COH}_{stem}.npz",
                        ko=np.array([k for k, _ in keys], dtype=str),
                        species=np.array([sp for _, sp in keys], dtype=str),
                        values=values, samples=np.array(names, dtype=str))
    withheld = tally[rule]['contradicted, withheld'] + tally[rule]['undeterminable, withheld']
    stats['rules'][rule] = dict(
        pairs=len(keys), conservation_ratio=alloc / tot,
        distinct_species=len({sp for _, sp in keys} - {UNRESOLVED, 'Host'}),
        distinct_KO=len({k for k, _ in keys} - {UNASSIGNED}),
        allocation_pct={k: 100 * x / tot for k, x in tally[rule].items()},
        withheld_pct_of_cohort_intensity=100 * withheld / tot)
    say(f"  {rule:11s} {len(keys):7,} pairs  conservation {alloc/tot:.10f}  "
        f"withheld {100*withheld/tot:.2f}% of cohort intensity")
    for k, x in sorted(tally[rule].items()):
        say(f"      {k:32s} {100*x/tot:6.2f}%")
    assert abs(alloc / tot - 1.0) < 1e-9, f'{rule} does not conserve intensity'
    del values
json.dump(stats, open(f"{OUT}/{COH}_gated_joint_stats.json", 'w'), indent=2)
say("GATED JOINT DONE")
