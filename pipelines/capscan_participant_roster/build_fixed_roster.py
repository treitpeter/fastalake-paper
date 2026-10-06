"""C2 (KC1177/JAG118): rebuild the CapScan sample->participant/region map and fix the identifier join.

Historical rule (CAPSCAN_NATMICRO/STOOL_PROVENANCE_RESOLVED.md): the map was "built by joining matrix
columns (run number) to updated_sample_overview_table.tsv" on the overview *Pellet* column, region
Capsule Type1..4 -> duodenum/jejunum/ileum/colon, location Stool -> stool; unmatched capsule runs were
labelled capsule_repeat_unassigned and excluded, giving conservative_331.tsv.

Defect: for participant 15 the overview Pellet field holds the workbook *supernatant* IDs (1961-1976),
while MS file names carry the *pellet* IDs (1799-1814). The direct join therefore misses all 16
participant-15 capsules.

Fix (only change): when the direct Pellet join finds no overview row, bridge MS pellet ID ->
workbook (Capscan_Main_PeTr.xlsx) Pellet -> workbook Main (device aliquot) -> overview Main.
The bridge is accepted only if it is unique at every step and Subject and Type agree between
workbook and overview. Region rule and inclusion rule are unchanged; labels of the 331 are unchanged.
"""
from pathlib import Path
import hashlib, json, re, sys
import pandas as pd

OUT = Path(sys.argv[1]); OUT.mkdir(parents=True, exist_ok=True)
P = Path('LAB_ROOT/..')
OVERVIEW = P/'PhD/Dari/metadata/updated_sample_overview_table.tsv'
WORKBOOK = P/'PhD/Dari/metadata/Capscan_Main_PeTr.xlsx'
HIST_MAP = P/'PhD/Dari/CAPSCAN_NATMICRO/precompute/authoritative_sample_map.tsv'
FROZEN = P/'FastaLake/canonical/publication/cohorts/CAPSCAN/conservative_331.tsv'
SEARCHED = P/'FastaLake/canonical/publication_data/arms/noclust_2026-08-07/cohorts/CAPSCAN/5_sage/razor_hash-acc'
REGION = {'Type1': 'duodenum', 'Type2': 'jejunum', 'Type3': 'ileum', 'Type4': 'colon'}
ORD = {'duodenum': 1, 'jejunum': 2, 'ileum': 3, 'colon': 4, 'stool': 5}

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def num(s): return pd.to_numeric(s, errors='coerce')

ov = pd.read_csv(OVERVIEW, sep='\t')
wb = pd.read_excel(WORKBOOK)
hist = pd.read_csv(HIST_MAP, sep='\t')
frozen = pd.read_csv(FROZEN, sep='\t')
ov_pellet, ov_main = num(ov.Pellet), num(ov.Main)
wb_pellet, wb_main = num(wb.Pellet), num(wb.Main)

def label(row):
    loc = str(row['Putative location']).strip().lower()
    if loc == 'stool': return 'stool'
    if loc == 'capsule' and row.Type in REGION: return REGION[row.Type]
    return None

def consistent(hits):
    # the overview repeats an aliquot once per DNA sequencing run; accept only if all repeats agree
    labs = {(label(r), float(r.Subject), r.Type) for _, r in hits.iterrows()}
    return hits.iloc[0] if len(hits) and len(labs) == 1 else None

def historical(run):
    hits = ov[ov_pellet.eq(run)]; r = consistent(hits)
    if r is None: return None
    reg = label(r)
    return None if reg is None else dict(region=reg, subject=int(r.Subject), how='direct_overview_Pellet', overview_rows=len(hits))

def bridged(run):
    w = wb[wb_pellet.eq(run)]
    if len(w) != 1: return None, f'workbook_pellet_rows={len(w)}'
    w = w.iloc[0]; main = num(pd.Series([w.Main])).iloc[0]
    hits = ov[ov_main.eq(main)]; r = consistent(hits)
    if r is None: return None, f'overview_main_rows={len(hits)}_or_inconsistent'
    reg = label(r)
    wtype = f"Type{int(w.Type)}" if pd.notna(num(pd.Series([w.Type])).iloc[0]) else str(w.Type)
    ok = reg is not None and float(r.Subject) == float(w.Subject) and r.Type == wtype
    ev = dict(region=reg, subject=int(r.Subject), how='bridge_workbook_Pellet_to_Main_to_overview_Main',
              overview_rows=len(hits), workbook_Main=int(main), workbook_Super=int(w.Super),
              overview_Pellet=int(r.Pellet), workbook_Subject=int(w.Subject), workbook_Set=str(w.Set),
              workbook_Type=wtype, overview_Type=r.Type, subject_and_type_agree=bool(ok))
    return (ev if ok else None), ev

rows = []
for sid, rn in zip(hist.sample_id, hist.run_number):
    # historical key = run number as parsed by the original map (its run_number column; 'None' for 1341_x)
    run = int(rn) if str(rn) not in ('None', 'nan') else None
    h = historical(run) if run is not None else None
    rows.append(dict(sample_id=sid, run_number=run, historical=h))

# Control 1: the historical direct join reproduces the authoritative map region/subject for every sample
mism = []
for r, (_, a) in zip(rows, hist.iterrows()):
    h = r['historical']
    exp = None if a.region == 'capsule_repeat_unassigned' else (a.region, int(a.subject))
    got = None if h is None else (h['region'], h['subject'])
    if exp != got: mism.append((r['sample_id'], exp, got))
assert not mism, mism[:10]
old = {r['sample_id']: r['historical'] for r in rows if r['historical']}
assert len(old) == 331
fz = {s: (reg, int(sub)) for s, reg, sub in zip(frozen.sample_id, frozen.region, frozen.subject_id)}
assert fz == {s: (h['region'], h['subject']) for s, h in old.items()}, 'historical join does not reproduce conservative_331'

# Fix: bridge only the rows the direct join could not place
new, evidence = {}, []
for r in rows:
    if r['historical'] or r['run_number'] is None: continue
    ev, info = bridged(r['run_number'])
    evidence.append(dict(sample_id=r['sample_id'], run_number=r['run_number'], resolved=ev is not None,
                         **(info if isinstance(info, dict) else dict(reason=info))))
    if ev: new[r['sample_id']] = ev
ev = pd.DataFrame(evidence); ev.to_csv(OUT/'join_fix_evidence.tsv', sep='\t', index=False)
assert set(v['subject'] for v in new.values()) == {15}, new
assert len(new) == 16 and all((SEARCHED/s).is_dir() for s in new)
# negative control: the bridge must not relabel any of the 331 (they already have direct matches)
for s, h in old.items():
    w = bridged(int(re.search(r'_T760_(\d+)', s).group(1)))[0]
    if w is not None: assert (w['region'], w['subject']) == (h['region'], h['subject']), s

fixed = [dict(sample_id=s, region=v['region'], region_ordinal=ORD[v['region']], subject_id=v['subject']) for s, v in {**old, **new}.items()]
fx = pd.DataFrame(fixed); fx = fx.set_index('sample_id').loc[[s for s in hist.sample_id if s in old or s in new]].reset_index()
fx.to_csv(OUT/'fixed_roster_347.tsv', sep='\t', index=False)
pd.DataFrame([dict(sample_id=s, region=v['region'], region_ordinal=ORD[v['region']], subject_id=v['subject'], **{k: v[k] for k in v if k not in ('region', 'subject')}) for s, v in new.items()]).to_csv(OUT/'participant15_added.tsv', sep='\t', index=False)
summary = dict(status='PASS', historical_join_reproduces_authoritative_map=True, historical_join_reproduces_conservative_331=True,
               frozen_rows=331, added_rows=len(new), fixed_rows=len(fx), participants=int(fx.subject_id.nunique()),
               added_by_region=pd.Series([v['region'] for v in new.values()]).value_counts().to_dict(),
               still_unassigned=int((~ev.resolved).sum()), still_unassigned_samples=ev.loc[~ev.resolved, 'sample_id'].tolist(),
               fixed_region_counts=fx.region.value_counts().to_dict(),
               inputs={str(p): sha(p) for p in [OVERVIEW, WORKBOOK, HIST_MAP, FROZEN]}, script_sha256=sha(__file__))
(OUT/'ROSTER_FIX.json').write_text(json.dumps(summary, indent=2) + '\n')
print(json.dumps(summary, indent=2))
