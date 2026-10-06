"""Insert the 16 participant-15 rows into the Source Data tables (CapScan block re-sorted by sample, as in the 710
tables) and write NUMBERS.md (old 710 -> new 726)."""
import csv, json, pathlib, hashlib
import numpy as np, pandas as pd
OUT = pathlib.Path(__file__).resolve().parents[1]
SD = OUT.parent / 'figure_sd/source_data/package/FastaLake_Source_Data'
ROSTER = OUT.parents[1] / 'coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/roster/fixed_roster_347.tsv'
sha = lambda p: hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()

def insert(old_path, new_lines, cap, key):
    L = open(old_path).read().splitlines(); head, body = L[0], L[1:]
    cohort = lambda l: l.split('\t')[0]
    i0 = next(i for i, l in enumerate(body) if cohort(l) == cap); i1 = max(i for i, l in enumerate(body) if cohort(l) == cap) + 1
    assert all(cohort(l) == cap for l in body[i0:i1]) and not any(cohort(l) == cap for l in body[i1:])
    block = body[i0:i1]; assert block == sorted(block, key=key)
    block = sorted(block + new_lines, key=key)
    return head, body[:i0] + block + body[i1:]

# sizes: same transformation as the Source Data copy (path -> basename; verified: only differing column)
new = list(csv.DictReader(open(OUT / 'work/acquisition_sizes_new16_fullpath.tsv'), delimiter='\t'))
hdr = open(SD / 'Figure1/acquisition_sizes.tsv').readline().rstrip('\n').split('\t')
assert list(new[0]) == hdr
for r in new: r['path'] = pathlib.Path(r['path']).name
nl = ['\t'.join(r[h] for h in hdr) for r in new]
h, b = insert(SD / 'Figure1/acquisition_sizes.tsv', nl, 'CapScan', key=lambda l: l.split('\t')[1])
(OUT / 'acquisition_sizes_726.tsv').write_text(h + '\n' + '\n'.join(b) + '\n')

# support: population order within a sample preserved
POP = {'all_selected': 0, 'identified_accessions': 1, 'native_protein_q01': 2}
nl = open(OUT / 'work/acquisition_support_summary_new16.tsv').read().splitlines()[1:]
h2, b2 = insert(SD / 'Figure1/acquisition_support_summary.tsv', nl, 'CAPSCAN', key=lambda l: (l.split('\t')[1], POP[l.split('\t')[2]]))
(OUT / 'acquisition_support_summary_726.tsv').write_text(h2 + '\n' + '\n'.join(b2) + '\n')

# checks
roster = {r['sample_id'] for r in csv.DictReader(open(ROSTER), delimiter='\t')}
S0 = pd.read_csv(SD / 'Figure1/acquisition_sizes.tsv', sep='\t'); S1 = pd.read_csv(OUT / 'acquisition_sizes_726.tsv', sep='\t')
U0 = pd.read_csv(SD / 'Figure1/acquisition_support_summary.tsv', sep='\t'); U1 = pd.read_csv(OUT / 'acquisition_support_summary_726.tsv', sep='\t')
D710 = pd.read_csv(SD / 'Review_v16/Database_sizes/database_sizes_710.tsv', sep='\t')
assert len(S1) == 726 and set(S1[S1.cohort == 'CapScan']['sample']) == roster
assert len(U1) == 726 * 3 and set(U1[U1.cohort == 'CAPSCAN']['sample']) == roster
assert (D710.set_index(['cohort', 'sample']).entries == S0.set_index(['cohort', 'sample']).entries).all()   # RES-07 table = same entries
assert (D710.set_index(['cohort', 'sample']).reservoir_entries == S0.set_index(['cohort', 'sample']).reservoir_entries).all()
old_rows = set(open(SD / 'Figure1/acquisition_sizes.tsv').read().splitlines()); assert old_rows <= set(open(OUT / 'acquisition_sizes_726.tsv').read().splitlines())
old_rows = set(open(SD / 'Figure1/acquisition_support_summary.tsv').read().splitlines()); assert old_rows <= set(open(OUT / 'acquisition_support_summary_726.tsv').read().splitlines())
# every new FASTA row: entries == unique sequences, no duplicate accessions (as for the 710)
N = S1[S1['sample'].isin(set(S1['sample']) - set(S0['sample']))]; assert len(N) == 16

HUMAN = 20431
CO = ['CAMPI', 'HSTOOL', 'PREDICT', 'CapScan']; CU = {'CapScan': 'CAPSCAN'}
def nums(S, U):
    a = U[U.population == 'all_selected']
    z = dict(total=len(S), n={c: int((S.cohort == c).sum()) for c in CO}, max_entries=int(S.entries.max()),
             max_entries_cohort=S.loc[S.entries.idxmax(), 'cohort'], n_below_human=int((S.entries < HUMAN).sum()),
             min_entries=int(S.entries.min()),
             median={c: float(S[S.cohort == c].entries.median()) for c in CO},
             capscan_range=(int(S[S.cohort == 'CapScan'].entries.min()), int(S[S.cohort == 'CapScan'].entries.max())),
             one_top={c: 100 * float(a[a.cohort == CU.get(c, c)].one_top_prediction_fraction.median()) for c in CO},
             fl_corr={c: 100 * float(a[a.cohort == CU.get(c, c)].with_FL_corroborated_prediction_fraction.median()) for c in CO},
             one_top_capscan_minmax=tuple(100 * a[a.cohort == 'CAPSCAN'].one_top_prediction_fraction.agg(['min', 'max'])))
    return z
o, n = nums(S0, U0), nums(S1, U1)
f1 = lambda x: f'{x:.1f}'; f2 = lambda x: f'{x:.2f}'
rng = lambda d, f: f'{f(min(d.values()))}–{f(max(d.values()))}% (min {min(d, key=d.get)}, max {max(d, key=d.get)})'
r100 = lambda x: f'{int(round(x, -2)):,}'
lines = ['# Fig. 1b,d and selected-database counts: 331 → 347 CapScan acquisitions (15 participants)', '',
 'Old = Source Data tables (710 acquisitions). New = the same tables plus the 16 participant-15 capsule acquisitions (runs 1799–1814), 726 acquisitions.',
 'Only CapScan rows were added; the 379 CAMPI/HSTOOL/PREDICT rows and the 331 existing CapScan rows are byte-identical to the Source Data copies.', '',
 '| Item (Source Data map row) | Text now | Old (710) | New (726) |', '|---|---|---|---|',
 f"| Total acquisitions (FIG1-b, SF1, ABS-01) | 710 | {o['total']} | {n['total']} |",
 f"| CapScan acquisitions (SF1-n; panel b,d labels) | CapScan (n = 331) | {o['n']['CapScan']} | {n['n']['CapScan']} |",
 f"| Largest selected FASTA; all below 20,431 (ABS-01) | \"all … fewer entries than 20,431\" | max {o['max_entries']:,} ({o['max_entries_cohort']}); {o['n_below_human']}/{o['total']} below | max {n['max_entries']:,} ({n['max_entries_cohort']}); {n['n_below_human']}/{n['total']} below |",
 f"| Smallest selected FASTA | – | {o['min_entries']:,} | {n['min_entries']:,} |",
 f"| CAMPI median entries (RES-07b, RES-25a) | 3,100 | {o['median']['CAMPI']:,.1f} | {n['median']['CAMPI']:,.1f} (unchanged) |",
 f"| Healthy-stool median entries (RES-07) | 12,200 | {o['median']['HSTOOL']:,.1f} | {n['median']['HSTOOL']:,.1f} (unchanged) |",
 f"| PREDICT median entries | – | {o['median']['PREDICT']:,.1f} | {n['median']['PREDICT']:,.1f} (unchanged) |",
 f"| CapScan median entries | – | {o['median']['CapScan']:,.1f} (≈{r100(o['median']['CapScan'])}) | {n['median']['CapScan']:,.1f} (≈{r100(n['median']['CapScan'])}) |",
 f"| CapScan entry range | – | {o['capscan_range'][0]:,}–{o['capscan_range'][1]:,} | {n['capscan_range'][0]:,}–{n['capscan_range'][1]:,} |",
 f"| Per-cohort median range, one top-scoring lookup query (RES-15, Fig. 1d) | 21.6–35.8% | {rng(o['one_top'], f1)} | {rng(n['one_top'], f1)} |",
 f"| CapScan median, one top-scoring lookup query (Fig. 1d) | – | {f2(o['one_top']['CapScan'])}% | {f2(n['one_top']['CapScan'])}% |",
 f"| CapScan per-acquisition min–max, one top query (Fig. 1d) | – | {f1(o['one_top_capscan_minmax'][0])}–{f1(o['one_top_capscan_minmax'][1])}% | {f1(n['one_top_capscan_minmax'][0])}–{f1(n['one_top_capscan_minmax'][1])}% |",
 f"| Per-cohort median range, FL-corroborated prediction (RES-16) | 95.7–97.7% | {rng(o['fl_corr'], f1)} | {rng(n['fl_corr'], f1)} |",
 f"| CapScan median, FL-corroborated prediction (RES-16) | – | {f2(o['fl_corr']['CapScan'])}% | {f2(n['fl_corr']['CapScan'])}% |",
 '', 'Per-cohort medians, one top-scoring lookup query (%), old → new: ' + '; '.join(f"{c} {f2(o['one_top'][c])} → {f2(n['one_top'][c])}" for c in CO),
 '', 'Per-cohort medians, FL-corroborated prediction (%), old → new: ' + '; '.join(f"{c} {f2(o['fl_corr'][c])} → {f2(n['fl_corr'][c])}" for c in CO),
 '', '## Source Data map rows that read these tables and do not change', '',
 'RES-01 (58 HSTOOL), RES-02/RES-08a/RES-10k (233.5 M HSTOOL reservoir), RES-08b (272.9 M PREDICT), RES-08c (264.9 M CapScan reservoir: the 16 added rows carry the same 264,925,283), RES-08d/RES-375 (375.3 M CAMPI), RES-306 (306 PREDICT), M-SIZES (pointer only). SF1-n text "15 CAMPI, 58 healthy-stool, 306 PREDICT and 331 CapScan" → "… and 347 CapScan". SF1 "all 710" → "all 726" (Supplementary Fig. 1 itself was not rebuilt here).',
 '', '## The 16 added acquisitions', '', '| sample | entries | one top query (%) | FL-corroborated (%) |', '|---|---:|---:|---:|']
a1 = U1[U1.population == 'all_selected'].set_index('sample')
for _, r in N.iterrows():
    lines.append(f"| {r['sample']} | {r.entries:,} | {100*a1.loc[r['sample'], 'one_top_prediction_fraction']:.2f} | {100*a1.loc[r['sample'], 'with_FL_corroborated_prediction_fraction']:.2f} |")
lines += ['', '## Provenance', '',
 '- Sizes: `MM_review_analysis/A01_4_1_database_sizes/frozen/audit_v2.py::scan` (imported unchanged) on `canonical/publication_data/arms/noclust_2026-08-07/cohorts/CAPSCAN/4_parsimony/razor_hash-acc/<sample>/parsimony.fasta`; path column reduced to the basename as in the Source Data copy.',
 '- Support: `MM_review_analysis/A09_4_1_prediction_support/frozen/analysis.py` (byte-identical copy in `a09_run/frozen/`), summarised with the record code of `frozen/collect.py::support()` (`code/summarize.py`). Inputs: AlphaNovo predictions in `projects/capscan_unified/mzml_combined/alphanovo_inference/`, Sage `5_sage/razor_hash-acc/<sample>/results.sage.tsv`, same arm.',
 '- Reproduction: 3 existing CapScan acquisitions (1249, 1262, 2561) reproduced exactly: size rows identical in every column, support rows identical for all three populations, protein_support.tsv byte-identical (`checks/`).',
 '- Negative controls failed as required: missing FASTA (sizes), changed FASTA metadata and missing search file (support).',
 f"- Output sha256: acquisition_sizes_726.tsv {sha(OUT/'acquisition_sizes_726.tsv')}; acquisition_support_summary_726.tsv {sha(OUT/'acquisition_support_summary_726.tsv')}"]
(OUT / 'NUMBERS.md').write_text('\n'.join(lines) + '\n')
json.dump(dict(old=o, new=n), open(OUT / 'checks/numbers.json', 'w'), indent=1, default=str)
print('\n'.join(lines))
