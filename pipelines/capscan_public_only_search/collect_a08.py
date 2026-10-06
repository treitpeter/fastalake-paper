#!/usr/bin/env python3
"""A08 public arm, 331 (production) vs 347 (OUT/A08_347): identity proof + matched-benchmark numbers.

1. Every production sample (331) re-run in A08_347: razor FASTA sha, counts, MG comparisons and every peptide-set file identical?
2. 16 participant-15 acquisitions: PASS, FL-only counts (no MG).
3. Matched benchmark (282, FL_public vs specimen MG), recomputed by one code path from the per-acquisition peptide sets of
   (a) production A08 and (b) A08_347. (a) is gated against the packaged Source Data (gain_by_region.tsv and
   ordinary_all_pairs.tsv q_only 9-50); a gate failure stops the collector (exit 1).
4. Cohort-level FL-only medians over all public-arm acquisitions (331 vs 347).
Writes OUT/collect/a08_*.tsv and a08_summary.json.
"""
import csv, hashlib, json, pathlib, statistics as st, sys
OUT = pathlib.Path(sys.argv[1]).resolve()
NEW = OUT / 'A08_347'
SELFTEST = '--selftest' in sys.argv  # results base = production; checks the gate and code path before the run exists
PROD = pathlib.Path('LAB_ROOT/MM_review_analysis/A08_4_5_capscan_public_only')
SD = pathlib.Path('LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/figure_sd/source_data/package/FastaLake_Source_Data')
ROSTER = pathlib.Path('LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/roster/fixed_roster_347.tsv')
sys.path.insert(0, str(NEW / 'frozen'))
RES = PROD if SELFTEST else NEW
from search_helpers import recount  # frozen, byte-identical copy

def sha(p): return hashlib.sha256(open(p, 'rb').read()).hexdigest()
def rset(p): return set(l.strip() for l in open(p) if l.strip())
def tsv(p, rows):
    with open(p, 'w') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(rows)

m = json.load(open(NEW / 'inputs/manifest.json'))
samples = [r['sample'] for r in m['samples']]
prod_samples = samples[:331]; new16 = samples[331:]
assert len(samples) == 347 and len(new16) == 16
ev = json.load(open(RES / 'evidence/COMPLETE.json')); assert ev['status'] == 'PASS' and ev['all_selected_sequences_publicly_provenanced']
evp = json.load(open(PROD / 'evidence/COMPLETE.json'))

# 1. identity of the 331
ident = []
for s in prod_samples:
    a, b = PROD / 'results' / s, RES / 'results' / s
    ca, cb = json.load(open(a / 'COMPLETE.json')), json.load(open(b / 'COMPLETE.json'))
    assert cb['status'] == 'PASS'
    row = dict(sample=s, fasta_identical=ca['selected_fasta']['sha256'] == cb['selected_fasta']['sha256'],
               proteins_old=ca['selected_proteins'], proteins_new=cb['selected_proteins'],
               peptides_old=ca['counts']['peptide_q01'], peptides_new=cb['counts']['peptide_q01'],
               counts_identical=ca['counts'] == cb['counts'], comparisons_identical=ca['comparisons'] == cb['comparisons'])
    row['sets_identical'] = all(sha(a / f.name) == sha(b / f.name) for f in a.glob('*.txt'))
    row['all_identical'] = row['fasta_identical'] and row['counts_identical'] and row['comparisons_identical'] and row['sets_identical']
    ident.append(row)
tsv(OUT / ('collect/a08_identity_331_SELFTEST.tsv' if SELFTEST else 'collect/a08_identity_331.tsv'), ident)
n_ident = sum(r['all_identical'] for r in ident)

# 2. the 16 new
roster = {r['sample_id']: r for r in csv.DictReader(open(ROSTER), delimiter='\t')}
new_rows = []
for s in ([] if SELFTEST else new16):
    c = json.load(open(NEW / 'results' / s / 'COMPLETE.json')); assert c['status'] == 'PASS'
    new_rows.append(dict(sample=s, region=roster[s]['region'], participant=roster[s]['subject_id'], selected_proteins=c['selected_proteins'],
                         peptides_5_50=c['counts']['peptide_q01'], peptides_9_50=c['counts']['length9_50'],
                         noclust_baseline_recovery=c['comparisons']['cohort_metagenomes_available']['comparator_recovery'],
                         specimen_metagenome='none'))
tsv(OUT / 'collect/a08_participant15_new16.tsv', new_rows) if new_rows else None

# 3. matched benchmark from per-acquisition sets
labels = {r['sample']: r['group'] for r in csv.DictReader(open(SD / 'Figure4/benchmark_regions/paired_benchmark_rows.tsv'), delimiter='\t') if r['cohort'] == 'CAPSCAN'}
if '--negative' in sys.argv: labels.pop(sorted(labels)[0])  # negative control: one matched acquisition missing -> must fail
assert len(labels) == 282, ('matched roster incomplete', len(labels))
mg = {r['sample']: r['metagenome']['path'] for r in m['samples'] if r['metagenome']}
assert set(mg) == set(labels)
mgsets = {}
for s, p in mg.items():
    _, sets, _ = recount(p); mgsets[s] = (sets['peptide_q01'], sets['length9_50'])

def bench(base):
    rows = []
    for s in sorted(mg):
        d = base / 'results' / s
        fl5, fl9 = rset(d / 'peptide_q01.txt'), rset(d / 'length9_50.txt'); m5, m9 = mgsets[s]
        rows.append(dict(sample=s, group=labels[s], n_fl_5=len(fl5), n_mg_5=len(m5), rec_5=len(fl5 & m5) / len(m5),
                         n_fl_9=len(fl9), n_mg_9=len(m9), rec_9=len(fl9 & m9) / len(m9) if m9 else float('nan')))
    return rows
def summ(rows):
    out = {}
    strata = {'all': lambda g: True, 'capsules pooled': lambda g: g != 'stool'}
    strata.update({g: (lambda g0: lambda g: g == g0)(g) for g in ['duodenum', 'jejunum', 'ileum', 'colon', 'stool']})
    for k, f in strata.items():
        x = [r for r in rows if f(r['group'])]
        out[k] = dict(n=len(x),
                      gain_5_50=100 * st.median(r['n_fl_5'] / r['n_mg_5'] for r in x) - 100,
                      gain_9_50=100 * st.median(r['n_fl_9'] / r['n_mg_9'] for r in x if r['n_mg_9']) - 100,
                      recovery_5_50=100 * st.median(r['rec_5'] for r in x),
                      fl_more_5_50=sum(r['n_fl_5'] > r['n_mg_5'] for r in x),
                      fl_more_9_50=sum(r['n_fl_9'] > r['n_mg_9'] for r in x),
                      median_fl_5=st.median(r['n_fl_5'] for r in x), median_mg_5=st.median(r['n_mg_5'] for r in x))
    return out
old_rows, new_rows_b = bench(PROD), bench(RES)
old, new = summ(old_rows), summ(new_rows_b)
None if SELFTEST else (tsv(OUT / 'collect/a08_benchmark282_old.tsv', old_rows), tsv(OUT / 'collect/a08_benchmark282_new.tsv', new_rows_b))

# gate: old recompute reproduces packaged Source Data
gate = {}
for r in csv.DictReader(open(SD / 'Figure4/benchmark_regions/gain_by_region.tsv'), delimiter='\t'):
    if r['cohort'] == 'CAPSCAN':
        gate[r['subset']] = (abs(old[r['subset']]['gain_5_50'] - float(r['median_gain_percent'])) < 1e-6,
                             abs(old[r['subset']]['recovery_5_50'] - 100 * float(r['median_identity_recovery'])) < 1e-6)
gate['9_50_stool_13.95'] = (abs(old['stool']['gain_9_50'] - 13.9524) < 1e-3,)
reg9 = [old[g]['gain_9_50'] for g in ['duodenum', 'jejunum', 'ileum', 'colon']]
gate['9_50_regions_min_46.63_max_64.4'] = (abs(min(reg9) - 46.6275) < 1e-3, round(max(reg9), 1) == 64.4)
gate['225_of_226'] = (old['capsules pooled']['fl_more_9_50'] == 225,)
gate['56_of_56'] = (old['stool']['fl_more_9_50'] == 56,)
gate['pooled_yield_56pct'] = (abs(100*old['capsules pooled']['median_fl_5']/old['stool']['median_fl_5'] - 55.8759) < 1e-3, abs(100*old['capsules pooled']['median_mg_5']/old['stool']['median_mg_5'] - 35.9) < 0.1)
gate_ok = all(all(v) for v in gate.values())

# 4. cohort FL-only medians
def cohort(base, ss):
    c = [json.load(open(base / 'results' / s / 'COMPLETE.json')) for s in ss]
    return dict(n=len(c), median_peptides_5_50=st.median(x['counts']['peptide_q01'] for x in c),
                median_selected_proteins=st.median(x['selected_proteins'] for x in c))
summary = dict(evidence=dict(old_selected_unique_sequences=evp['selected_unique_sequences'], new_selected_unique_sequences=ev['selected_unique_sequences'],
                             old_unique_peptides=evp['stats']['total_peptides'], new_unique_peptides=ev['stats']['total_peptides'],
                             new_all_public=ev['all_selected_sequences_publicly_provenanced'], new_source_counts=ev['source_counts']),
               identity_331=dict(identical=n_ident, of=331), new16=len(new_rows),
               benchmark282=dict(old=old, new=new, identical=old == new), gate_old_vs_source_data=dict(ok=gate_ok, checks={k: list(v) for k, v in gate.items()}),
               cohort_public_arm=dict(old=cohort(PROD, prod_samples), new=cohort(RES, prod_samples if SELFTEST else samples)))
json.dump(summary, open(OUT / ('collect/a08_summary_SELFTEST.json' if SELFTEST else 'collect/a08_summary.json'), 'w'), indent=1)
print(json.dumps(dict(identity=n_ident, gate_ok=gate_ok, bench_identical=old == new)))
sys.exit(0 if gate_ok else 1)
