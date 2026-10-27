"""Item 1 collector: Supp. Fig. 4a,b / Supp. Methods 7 minimum-length-7 matched CapScan benchmark at 295.

read() and the matched-row formula are verbatim copies of M09 wave2/code/collect_peptides.py read() and
summarize_peptides.py (those scripts run on import, so they cannot be imported). The text statistics are
then taken by the release-repo producers themselves: stats.py:stat_same_spectrum and panels.py:stat_ab_matched.

Modes
  --reproduce282 : recompute every CapScan (and PREDICT) FL_public/MG row of M09 matched_comparisons.tsv from the
                   archived results.sage.tsv files; require equality; require the producers to give the committed
                   282 values (44.9721893631811 / 84.01297879503241 CapScan; 14.455684375777889 / 85.39185872471847
                   PREDICT). Negative control: one perturbed reference count must be detected.
  --final        : search control (1311 re-search equals the archived M09 MG search, sets at 5/7/9), consistency gates
                   for the 13 (FL_public rows equal public347 per_acquisition_new16.tsv; MG baseline counts equal
                   mg295/mg_full COUNTS.json), then 295-row tables and 295 statistics.
"""
from pathlib import Path
import sys, json, csv, statistics as st, importlib.util, shutil, copy
from multiprocessing import Pool
F = "LAB_ROOT"
sys.path.insert(0, F + '/MM_review_analysis/A10_4_3_4_5_identity_recovery/frozen')
from peptide import canonical
csv.field_size_limit(1 << 30)
RUN = Path(__file__).resolve().parent.parent
MG295 = RUN.parent.parent
FINAL = MG295.parent
REPO = MG295 / 'phase2/repo'
M09 = Path(F + "/00_PACKAGE_2026-09-17/cleanup_2026-09-22/v10_mm_philip_comment_audit_2026-09-22/nine_point_followup/manuscript_final_pass_2026-09-23/mm_v11_revision_2026-09-23/reappraisal_2026-09-24/main_methods_pass_2026-09-24/review_decisions/M09_short_peptides/run_2026-09-24_minlen7")
P347 = FINAL / 'public347'
POPS = ['baseline5', 'nested7', 'nested9', 'fresh7', 'fresh7_nested9']
CONTROL = "20210513_EXPL0_PeTr_SA_T760_1311"

def read(p):
 sets={5:set(),7:set(),9:set()};proteins=set();hist={};n=0
 with Path(p).open() as f:
  for r in csv.DictReader(f,delimiter='\t'):
   n+=1
   if r['rank']!='1' or r['label']!='1' or float(r['peptide_q'])>.01:continue
   seq=canonical(r['peptide'])
   if not 5<=len(seq)<=50:continue
   for low in sets:
    if len(seq)>=low:sets[low].add(seq)
   if float(r['protein_q'])<=.01:proteins.update(x for x in r['proteins'].split(';') if x and not x.startswith('rev_'))
 for seq in sets[5]:hist[len(seq)]=hist.get(len(seq),0)+1
 return sets,dict(psm_rows=n,native_identified_accessions=len(proteins),canonical_length_hist=hist)

def pops(base_tsv, fresh_tsv):
    base, _ = read(base_tsv); new, _ = read(fresh_tsv); assert new[5] == new[7]
    return {'baseline5': base[5], 'nested7': base[7], 'nested9': base[9], 'fresh7': new[7], 'fresh7_nested9': new[9]}

def matched(co, a, b, name, A, B):
    out = []
    for pop in POPS:
        aa, bb = A[pop], B[pop]
        out.append(dict(cohort=co, arm=a, comparator=b, sample=name, population=pop, arm_peptides=len(aa), comparator_peptides=len(bb), shared=len(aa & bb),
                        gain_pct=100*(len(aa)/len(bb)-1) if bb else None, comparator_identity_recovery_pct=100*len(aa & bb)/len(bb) if bb else None))
    return out

def wide(co, arm, sample, P):
    r = dict(cohort=co, arm=arm, sample=sample, baseline5=float(len(P['baseline5'])), fresh7=float(len(P['fresh7'])), fresh7_nested9=float(len(P['fresh7_nested9'])),
             nested7=float(len(P['nested7'])), nested9=float(len(P['nested9'])), gained_fresh_nested7=len(P['fresh7'] - P['nested7']), lost_fresh_nested7=len(P['nested7'] - P['fresh7']))
    r['fresh7_vs_base5'] = 100 * (r['fresh7'] / r['baseline5'] - 1); r['fresh7_vs_nested7'] = 100 * (r['fresh7'] / r['nested7'] - 1)
    r['nested7_vs_base5'] = 100 * (r['nested7'] / r['baseline5'] - 1)
    return r

def load_mod(name, path):
    spec = importlib.util.spec_from_file_location(name, path); m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m

def stats_values(sd_supp, sd_sf4):
    """295/282 values from the repo producers themselves."""
    S = load_mod('stats67', REPO / 'text_statistics/supplementary_methods_6_7/stats.py')
    ss = S.stat_same_spectrum(sd_supp)
    import pandas as pd
    sys.path.insert(0, str(REPO))
    try:
        Pn = load_mod('panels_sf4', REPO / 'figures/SFig04/panels.py')
        mc = pd.read_csv(Path(sd_sf4) / 'benchmark_matched_per_acquisition.tsv', sep='\t')
        d = dict(mc=mc)
        ab = {f'{c}_{p}': Pn.stat_ab_matched(d, c, p) for c in ('PREDICT', 'CAPSCAN') for p in ('fresh7', 'nested7', 'baseline5', 'nested9')}
        g58 = {c: Pn.stat_gain_5_8(d, c) for c in ('PREDICT', 'CAPSCAN')}
    except Exception as e:
        ab = dict(error=repr(e)); g58 = None
    return dict(stats_py=ss, panels_ab=ab, panels_gain_5_8=g58)

def _m09_pair(args):
    co, name, tF, tM = args
    A = pops(tF['baseline_dir'] + '/results.sage.tsv', tF['out_dir'] + '/results.sage.tsv')
    B = pops(tM['baseline_dir'] + '/results.sage.tsv', tM['out_dir'] + '/results.sage.tsv')
    return matched(co, 'FL_public', 'MG', name, A, B)

def same(a, b):
    for k in ('arm_peptides', 'comparator_peptides', 'shared'):
        if int(a[k]) != int(b[k]): return False
    for k in ('gain_pct', 'comparator_identity_recovery_pct'):
        if abs(float(a[k]) - float(b[k])) > 1e-9: return False
    return True

def reproduce282():
    ref = [r for r in csv.DictReader(open(M09 / 'wave2/peptide_summary/matched_comparisons.tsv'), delimiter='\t')
           if r['arm'] == 'FL_public' and r['comparator'] == 'MG' and r['cohort'] in ('PREDICT', 'CAPSCAN')]
    T = {(r['cohort'], r['arm'], r['sample']): r for r in csv.DictReader(open(M09 / 'tasks.tsv'), delimiter='\t')}
    names = sorted({(r['cohort'], r['sample']) for r in ref})
    assert sum(c == 'CAPSCAN' for c, _ in names) == 282 and sum(c == 'PREDICT' for c, _ in names) == 296
    with Pool(16) as p: rows = [x for g in p.map(_m09_pair, [(c, n, T[c, 'FL_public', n], T[c, 'MG', n]) for c, n in names], chunksize=4) for x in g]
    R = {(r['cohort'], r['sample'], r['population']): r for r in ref}
    mism = [(r['cohort'], r['sample'], r['population']) for r in rows if not same(r, R[r['cohort'], r['sample'], r['population']])]
    # producer on the committed repo source tables (282)
    v = stats_values(REPO / 'source_data/supplementary_methods_6_7', REPO / 'source_data/SFig04')
    # producer on tables rebuilt from my recomputed rows (282)
    tmp = RUN / 'sd282_rebuilt'; build_sd(tmp, rows, [])
    v2 = stats_values(tmp / 'supplementary_methods_6_7', tmp / 'SFig04')
    want = dict(CAPSCAN_gain=44.9721893631811, CAPSCAN_rec=84.01297879503241, PREDICT_gain=14.455684375777889, PREDICT_rec=85.39185872471847)
    def pick(v):
        s = v['stats_py']; return dict(CAPSCAN_gain=s['ml7_CAPSCAN_gain'], CAPSCAN_rec=s['ml7_CAPSCAN_recovery'], PREDICT_gain=s['ml7_PREDICT_gain'], PREDICT_rec=s['ml7_PREDICT_recovery'])
    ok_committed = all(abs(pick(v)[k] - want[k]) < 1e-12 for k in want)
    ok_rebuilt = all(abs(pick(v2)[k] - want[k]) < 1e-12 for k in want)
    # negative control: perturb one reference count; the same comparison must flag it
    bad = copy.deepcopy(R); k0 = next(k for k in bad if k[0] == 'CAPSCAN' and k[2] == 'fresh7'); bad[k0]['shared'] = str(int(bad[k0]['shared']) + 1)
    neg_detected = sum(not same(r, bad[r['cohort'], r['sample'], r['population']]) for r in rows) == 1
    res = dict(rows_recomputed=len(rows), mismatches=len(mism), mismatch_examples=mism[:5], committed_tables=pick(v), rebuilt_tables=pick(v2), expected=want,
               panels_ab_committed=v['panels_ab'], panels_gain_5_8_committed=v['panels_gain_5_8'], panels_ab_rebuilt=v2['panels_ab'],
               negative_control_perturbed_one_count_detected=neg_detected,
               status='PASS' if (not mism and ok_committed and ok_rebuilt and neg_detected and len(rows) == 578 * 5) else 'FAIL')
    (RUN / 'CONTROL_282.json').write_text(json.dumps(res, indent=1, default=str) + '\n'); print(json.dumps(res, default=str)[:3000])
    if res['status'] != 'PASS': sys.exit(1)

def build_sd(dest, m09rows, newrows):
    """Copy the two committed source tables; replace their FL_public/MG rows by m09rows only when given rows equal
    them (checked in reproduce282); append newrows. Columns and populations as the committed tables."""
    dest = Path(dest)
    for unit, fn, popset in (('supplementary_methods_6_7', 'minlen7_benchmark_per_acquisition.tsv', ('baseline5', 'fresh7', 'nested9')),
                             ('SFig04', 'benchmark_matched_per_acquisition.tsv', tuple(POPS))):
        (dest / unit).mkdir(parents=True, exist_ok=True)
        src = REPO / 'source_data' / unit / fn
        lines = src.read_text().splitlines(); hdr = lines[0].split('\t')
        if m09rows:
            order = {(r['cohort'], r['sample'], r['population']): r for r in m09rows}
            body = []
            for l in lines[1:]:
                c = dict(zip(hdr, l.split('\t')))
                r = order[c['cohort'], c['sample'], c['population']]
                body.append('\t'.join(str(r[h]) for h in hdr))
        else:
            body = lines[1:]
        body += ['\t'.join(str(r[h]) for h in hdr) for r in newrows if r['population'] in popset]
        (dest / unit / fn).write_text('\n'.join([lines[0]] + body) + '\n')
    s67 = REPO / 'source_data/supplementary_methods_6_7/confidence_ordinary_all_pairs.tsv'
    shutil.copy(s67, dest / 'supplementary_methods_6_7' / s67.name)

def final():
    tasks = list(csv.DictReader(open(RUN / 'tasks.tsv'), delimiter='\t')); assert len(tasks) == 14
    for t in tasks:
        rec = json.load(open(RUN / 'receipts' / ('receipt_%04d.json' % int(t['index'])))); assert rec['pass'], t['sample']
    # (1) search control
    t0 = tasks[0]; assert t0['sample'] == CONTROL
    m = [r for r in csv.DictReader(open(M09 / 'tasks.tsv'), delimiter='\t') if r['cohort'] == 'CAPSCAN' and r['arm'] == 'MG' and r['sample'] == CONTROL][0]
    mine, _ = read(t0['out_dir'] + '/results.sage.tsv'); arch, _ = read(m['out_dir'] + '/results.sage.tsv')
    ctl = {str(k): dict(rerun=len(mine[k]), archived=len(arch[k]), identical=mine[k] == arch[k]) for k in (5, 7, 9)}
    ctl_ok = all(v['identical'] for v in ctl.values())
    # (2) the 13
    acq = [t for t in tasks[1:]]
    new16 = {r['sample']: r for r in csv.DictReader(open(P347 / 'minlen7_public347/per_acquisition_new16.tsv'), delimiter='\t')}
    rows, widesMG, gates = [], [], []
    for t in acq:
        s = t['sample']
        A = pops(f"{P347}/A08_347/results/{s}/search/results.sage.tsv", f"{P347}/minlen7_public347/searches/{s}/results.sage.tsv")
        B = pops(t['baseline_dir'] + '/results.sage.tsv', t['out_dir'] + '/results.sage.tsv')
        wa = wide('CAPSCAN', 'FL_public', s, A); n16 = new16[s]
        g1 = all(float(n16[k]) == float(wa[k]) for k in ('baseline5', 'fresh7', 'fresh7_nested9', 'nested7', 'nested9', 'gained_fresh_nested7', 'lost_fresh_nested7'))
        cj = json.load(open(t['baseline_dir'] + '/COUNTS.json'))
        g2 = cj['peptides_5_50'] == len(B['baseline5']) and cj['peptides_9_50'] == len(B['nested9'])
        gates.append(dict(sample=s, FL_public_equals_public347_new16=g1, MG_baseline_equals_mg_full_COUNTS=g2))
        rows += matched('CAPSCAN', 'FL_public', 'MG', s, A, B); widesMG.append(wide('CAPSCAN', 'MG', s, B))
    gates_ok = all(g['FL_public_equals_public347_new16'] and g['MG_baseline_equals_mg_full_COUNTS'] for g in gates)
    hdr = list(rows[0]);
    with open(RUN / 'matched_comparisons_new13.tsv', 'w') as f:
        w = csv.DictWriter(f, fieldnames=hdr, delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(rows)
    with open(RUN / 'per_acquisition_MG_new13.tsv', 'w') as f:
        w = csv.DictWriter(f, fieldnames=list(widesMG[0]), delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(widesMG)
    sd = RUN / 'sd295'; build_sd(sd, None, rows)
    v295 = stats_values(sd / 'supplementary_methods_6_7', sd / 'SFig04')
    v282 = stats_values(REPO / 'source_data/supplementary_methods_6_7', REPO / 'source_data/SFig04')
    # panel c MG median (diagnostic; not printed in text)
    W = [r for r in csv.DictReader(open(M09 / 'wave2/peptide_summary/per_acquisition_wide_claude.tsv'), delimiter='\t') if r['cohort'] == 'CAPSCAN' and r['arm'] == 'MG']
    pc = dict(n282=len(W), filtered_to_7_282=st.median(float(r['nested7_vs_base5']) for r in W), new_search_7_282=st.median(float(r['fresh7_vs_base5']) for r in W),
              filtered_to_7_295=st.median([float(r['nested7_vs_base5']) for r in W] + [r['nested7_vs_base5'] for r in widesMG]),
              new_search_7_295=st.median([float(r['fresh7_vs_base5']) for r in W] + [r['fresh7_vs_base5'] for r in widesMG]))
    res = dict(search_control_1311=ctl, search_control_status='PASS' if ctl_ok else 'FAIL', gates_13=gates, gates_13_status='PASS' if gates_ok else 'FAIL',
               v282=v282, v295=v295, panel_c_capscan_MG=pc,
               search_count=dict(before=2022, standard_before=2017, after=2022 + 13, standard_after=2017 + 13))
    (RUN / 'RESULTS_295.json').write_text(json.dumps(res, indent=1, default=str) + '\n'); print(json.dumps(res, default=str)[:4000])
    if not (ctl_ok and gates_ok): sys.exit(1)

if __name__ == '__main__':
    {'--reproduce282': reproduce282, '--final': final}[sys.argv[1]]()
