"""CR-003 follow-up: rebuild Supplementary Tables 3 and 4 for one variant and extract every printed number.

Inputs are the shadow trees in ../shadow/<variant> (iteration-1 analysis outputs plus this window's figure builds).
Variant 'replay' is the control: its tables must equal the frozen tables cell for cell and its numbers must equal
the printed manuscript values. Nothing outside this window's directory is written.
"""
from pathlib import Path
import sys, json, hashlib
import numpy as np, pandas as pd, openpyxl

O = Path(__file__).resolve().parents[1]
V = sys.argv[1]
S = O/'shadow'/V; H = S/'B'
BR = Path('/f[s]/gpfs41/lv12/fileset02/pool/pool-mann-projects/Peter/FastaLake/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit')
PR = BR/'passage_review_2026-10-03'
T = O/'tables'/V; T.mkdir(parents=True, exist_ok=True)
N = O/'old_vs_new'/'per_variant'; N.mkdir(parents=True, exist_ok=True)
rd = lambda p, **k: pd.read_csv(p, sep='\t', **k)
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
rows = []
def rec(item, location, quantity, value, printed=None, source=''):
    rows.append(dict(variant=V, item=item, location=location, quantity=quantity,
                     value=value if not isinstance(value, (np.floating, np.integer)) else value.item(),
                     printed=printed, source=source))

# ---------------- Supplementary Table 3 (EDS display and focused update) ----------------
features = rd(H/'eds_refit/input_data/heatmap_features.tsv'); columns = rd(H/'eds_refit/input_data/heatmap_columns.tsv')
raw = np.load(H/'eds_refit/input_data/heatmap_values.npz')['raw']; patient = np.load(H/'eds_refit/input_data/heatmap_85_values.npz')
st3 = T/'Supplementary_Table_3.xlsx'
with pd.ExcelWriter(st3, engine='openpyxl') as w:   # verbatim view logic of build_review_packet.py
    pd.DataFrame(raw, index=features.feature_id, columns=columns.sample_id).to_excel(w, sheet_name='All20_features_259_observations')
    features.to_excel(w, sheet_name='Fixed_feature_key', index=False); columns.to_excel(w, sheet_name='Observation_metadata', index=False)
    rd(H/'eds_refit/input_data/heatmap_85_patient_order.tsv').to_excel(w, sheet_name='85_patient_order', index=False)
    pd.DataFrame(patient['raw'], index=features.feature_id, columns=patient['samples']).to_excel(w, sheet_name='85_patient_raw')
    pd.DataFrame(patient['z'], index=features.feature_id, columns=patient['samples']).to_excel(w, sheet_name='85_patient_z')
    pd.DataFrame({'Note': ['Fixed earlier selection; no reselection. HV366 remains in all20 tables and has no positive signal in the patient display.',
        'Human rows are human-assigned peptide signal within retained reporting groups; pathway rows exclude human-assigned signal and retain unresolved origin.',
        'Non-detections are not imputed. This descriptive display is not independent EDS validation.']}).to_excel(w, sheet_name='Read_me', index=False)
book = openpyxl.load_workbook(st3)
notes = book.create_sheet('Current_EDS_update', 0)
for line in ['3 October 2026: Fig. 4h now displays three post hoc metabolite–enzyme pairs in 85 patients.',
             'Supplementary Fig. 10a enlarges that display; 10b retains the earlier 19-feature overview.',
             'All previous sheets and their values are preserved; the earlier Read_me describes that original selection.',
             'Focused_raw contains the six displayed measurements and patient order; non-detections remain missing in the plots.',
             'Focused_pairs reports saved q values from complete screen families, not a correction restricted to the three examples.',
             'Focused_adjustments contains all 192 descriptive model estimates; no new discovery p values or intervals were computed.',
             'The full screen tables, original display files and analysis code accompany EDS_focused_20261003_SOURCE_DATA.zip.']:
    notes.append([line])
ES = S/'eds/EDS_approved/source_data'
for name, file in [('Focused_raw', 'focused_raw_measurements.tsv'), ('Focused_features', 'focused_features.tsv'),
                   ('Focused_pairs', 'selected_pair_summary.tsv'), ('Focused_adjustments', 'separate_adjustments_all_models.tsv')]:
    df = rd(ES/file); sh = book.create_sheet(name); sh.append(list(df.columns)); sh.freeze_panes = 'A2'
    for row in df.itertuples(index=False, name=None): sh.append([None if pd.isna(v) else v for v in row])
book.save(st3)

# ---------------- Supplementary Table 4 (individual functions) ----------------
ref = rd(H/'individual_refit/crossomics_same_function_correlations.tsv')
base = rd(BR/'individual_baseline/crossomics_same_function_correlations.tsv')   # previous allocation, frozen
keys = ['level', 'layer_a', 'layer_b', 'feature']
m = ref.merge(base[keys+['n_patients', 'rho']], on=keys, how='outer', suffixes=('_refit', '_baseline'), indicator='eligibility', validate='one_to_one')
ver = []
for (l, a, b), g in m.groupby(['level', 'layer_a', 'layer_b'], sort=True):
    both = g[g.eligibility == 'both']; new = g[g.eligibility == 'left_only']
    ver.append(dict(level=l, layer_a=a, layer_b=b, n_common=len(both), n_new=len(new), n_lost=int((g.eligibility == 'right_only').sum()),
                    median_refit_common=both.rho_refit.median(), median_baseline_common=both.rho_baseline.median(),
                    median_refit_new=new.rho_refit.median() if len(new) else np.nan))
ver = pd.DataFrame(ver); ver.to_csv(T/'common_feature_verification.tsv', sep='\t', index=False)
m.to_csv(T/'eligibility_matched_correlations.tsv.gz', sep='\t', index=False)
ntests = len(ref)
st4 = T/'Supplementary_Table_4.xlsx'
with pd.ExcelWriter(st4, engine='openpyxl') as w:
    for file, sheet in [('crossomics_same_function_correlations.tsv', f'All_{ntests}_tests'), ('same_function_summary.tsv', 'Nine_family_summary'), ('eligibility.tsv', 'All_eligibility')]:
        rd(H/'individual_refit'/file).to_excel(w, sheet_name=sheet, index=False)
book = openpyxl.load_workbook(st4)
if V == 'replay':   # control: reproduce the frozen workbook exactly, including its note sheet
    lines = ['3 October 2026: Peter approved a common-eligibility sensitivity and descriptive Results ending.',
             'Every earlier sheet/value remains unchanged, including all 30,632 original tests and complete correction families.',
             'Eligibility_common matches functional representation, layer pair and function identifier across previous/revised allocations.',
             'No profiles were renormalized and no correlations or p/q values were refitted for this subset.',
             'Full per-function membership and baseline/refit correlations are in Crossomics_common_eligibility_SOURCE_DATA.zip.',
             'The native annotation verification (CR-003) remains a separate dependency; this workbook uses the current frozen refit.']
    notes = book.create_sheet('Current_eligibility_update', 0)
else:               # neutral Read_me (Codex go/no-go item 3); no review-internal text
    lines = [f'Individual-function Spearman correlations between microbial DNA, RNA and protein layers across PREDICT patients ({ntests:,} tests).',
             'Native UHGG eggNOG annotations were used only after the protein aligned to the annotation’s seed orthologue with at least half of the reported bitscore.',
             'Benjamini–Hochberg correction was applied separately within each of the nine representation-by-layer-pair families.',
             'Eligibility_common compares functions eligible under the previous and revised allocations, matched by representation, layer pair and function identifier.',
             'No profiles were renormalized and no correlations or p/q values were refitted for the common-eligibility subset.']
    notes = book.create_sheet('Read_me', 0)
for line in lines: notes.append([line])
sh = book.create_sheet('Eligibility_common'); sh.append(list(ver.columns)); sh.freeze_panes = 'A2'
for row in ver.itertuples(index=False, name=None): sh.append([None if pd.isna(v) else v for v in row])
book.save(st4)

# ---------------- control: replay tables equal the frozen tables ----------------
def sheets(p):
    b = openpyxl.load_workbook(p, data_only=True); return {s.title: [tuple(r) for r in s.values] for s in b}
def same(a, b):
    if a.keys() != b.keys(): return False, f'sheets {list(a)} vs {list(b)}'
    for k in a:
        if len(a[k]) != len(b[k]): return False, f'{k} rows {len(a[k])} vs {len(b[k])}'
        for r1, r2 in zip(a[k], b[k]):
            for x, y in zip(r1, r2):
                if x == y or (isinstance(x, float) and isinstance(y, float) and (np.isclose(x, y, rtol=1e-12, atol=0) or (np.isnan(x) and np.isnan(y)))): continue
                return False, f'{k}: {r1[:3]} vs {r2[:3]} ({x!r} vs {y!r})'
    return True, 'equal'
frozen3 = PR/'EDS_approved/tables/Supplementary_Table3_EDS_updated.xlsx'
frozen4 = PR/'crossomics_approval/Supplementary_Table4_common_eligibility.xlsx'
cmp = {'ST3_vs_frozen': same(sheets(st3), sheets(frozen3)), 'ST4_vs_frozen': same(sheets(st4), sheets(frozen4))}
frozen_ver = rd(PR/'crossomics_approval/common_feature_verification.tsv')
cmp['common_feature_verification_vs_frozen'] = (bool(np.allclose(frozen_ver.select_dtypes('number').to_numpy(float), ver.select_dtypes('number').to_numpy(float), equal_nan=True)) if frozen_ver.shape == ver.shape else False, str(ver.shape))
if V == 'replay':
    assert all(v[0] for v in cmp.values()), cmp

# ---------------- printed numbers ----------------
sp = rd(ES/'selected_pair_summary.tsv')
for r in sp.itertuples():
    rec('XML71', 'Results XML71', f'Spearman rho {r.metabolite}-{r.KO}', r.rho, f'{r.rho:.2f}'.replace('-', '−'), 'EDS_approved/source_data/selected_pair_summary.tsv')
    rec('XML71', 'Results XML71', f'n {r.metabolite}-{r.KO}', r.n, str(r.n), 'selected_pair_summary.tsv')
    rec('XML71', 'Results XML71 (none pass after joint adjustment)', f'joint-adjusted full-family q {r.metabolite}-{r.KO}', r.joint_adjusted_full_family_q, 'none <0.05' if r.joint_adjusted_full_family_q >= .05 else 'PASSES', 'selected_pair_summary.tsv')
    rec('Supp Table 3 Focused_pairs', 'Supplementary Table 3', f'full-family q {r.metabolite}-{r.KO}', r.full_family_q, f'{r.full_family_q:.2g}', 'selected_pair_summary.tsv')
ff = rd(ES/'focused_features.tsv')
for r in ff.itertuples():
    rec('Supp Table 3 Focused_features', 'Supplementary Table 3', f'EDS rho {r.feature}', r.EDS_rho, f'{r.EDS_rho:.2f}', 'focused_features.tsv')
    rec('Supp Table 3 Focused_features', 'Supplementary Table 3', f'EDS q {r.feature}', r.EDS_full_family_q, f'{r.EDS_full_family_q:.2g}', 'focused_features.tsv')
mods = rd(ES/'separate_adjustments_all_models.tsv')
rec('XML379', 'Methods XML379', 'descriptive adjustment estimates', len(mods), str(len(mods)), 'separate_adjustments_all_models.tsv')
# candidate filter (XML379: "The candidate set comprised 21 native-signal pairs")
aps = rd(H/'eds_refit/augmented_pair_screen.tsv'); eds = rd(H/'eds_refit/all_EDS_associations.tsv')
e = eds[(eds.normalization == 'native') & (eds.model == 'unadjusted')]
met = e[e.layer == 'Metabolite'].set_index('feature'); ko = e[e.layer == 'KO'].set_index('feature')
c = aps[(aps.layer == 'KO') & (aps.normalization == 'native') & (aps.model == 'unadjusted') & (aps.q_augmented < .05) & (aps.rho < 0)].copy()
c = c[c.metabolite.isin(met.index) & c.protein_feature.isin(ko.index)]
c = c[(met.loc[c.metabolite, 'rho'].to_numpy() > 0) & (met.loc[c.metabolite, 'q'].to_numpy() < .05) & (ko.loc[c.protein_feature, 'rho'].to_numpy() < 0) & (ko.loc[c.protein_feature, 'q'].to_numpy() < .05)]
c.to_csv(T/'candidate_inverse_pairs.tsv', sep='\t', index=False)
rec('XML379', 'Methods XML379', 'candidate native-signal inverse pairs (q<0.05 in complete unadjusted families)', len(c), str(len(c)), 'eds_refit/augmented_pair_screen.tsv + all_EDS_associations.tsv')
for k, m3, ko3 in [('K01812', 'HILIC_FT23003', 'K01812'), ('K01667', 'HILIC_FT24070', 'K01667'), ('K01442', 'C18_F2423', 'K01442')]:
    rec('XML379', 'Methods XML379', f'displayed pair {m3}-{ko3} in candidate set', bool(((c.metabolite == m3) & (c.protein_feature == ko3)).any()), '', 'candidate_inverse_pairs.tsv')
# complete screens (Source Data / Supplementary Table 3 text)
rec('EDS screens', 'Source Data EDS / ST3', 'significant EDS associations (q<0.05, all layers, normalizations, models)', int((eds.q < .05).sum()), '', 'eds_refit/all_EDS_associations.tsv')
rec('EDS screens', 'Source Data EDS / ST3', 'significant augmented pair rows (q_augmented<0.05)', int((aps.q_augmented < .05).sum()), '', 'eds_refit/augmented_pair_screen.tsv')
rec('EDS screens', 'Source Data EDS / ST3', 'significant pair rows (q<0.05)', int((aps.q < .05).sum()), '', 'eds_refit/augmented_pair_screen.tsv')
# XML72
nat = rd(H/'cross_refit/native_correlations.tsv'); fun = rd(H/'cross_refit/functional_correlations.tsv')
for a in ['DNA', 'RNA', 'Metabolites']:
    r = nat[(nat.cohort == 'PREDICT') & (nat.analysis == 'sensitivity_all_features') & (((nat.layer_a == a) & (nat.layer_b == 'Non-host protein')) | ((nat.layer_b == a) & (nat.layer_a == 'Non-host protein')))].iloc[0]
    rec('XML72', 'Results XML72', f'native distance rho Non-host protein-{a}', r.rho, f'{r.rho:.2f}', 'cross_refit/native_correlations.tsv')
r = fun[(fun.cohort == 'PREDICT') & (fun.scope == 'displayed77') & (fun.level == 'Pathway') & (fun.prevalence == 0) & (fun.layer_a == 'metaT') & (fun.layer_b == 'microbial_metaP')].iloc[0]
rec('XML72', 'Results XML72', 'pathway RNA-protein distance rho', r.rho, f'{r.rho:.2f}', 'cross_refit/functional_correlations.tsv')
rec('XML72', 'Results XML72', 'pathway RNA-protein q', r.q, f'{r.q:.3f}', 'cross_refit/functional_correlations.tsv')
# XML73 and ST4
sm = rd(H/'individual_refit/same_function_summary.tsv')
prot = sm[sm.layer_b == 'microbial_metaP']; dr = sm[(sm.layer_a == 'metaG') & (sm.layer_b == 'metaT')]
rec('XML73', 'Results XML73', 'DNA/RNA-protein median rho, min', prot.median_rho.min(), f'{prot.median_rho.min():.2f}', 'individual_refit/same_function_summary.tsv')
rec('XML73', 'Results XML73', 'DNA/RNA-protein median rho, max', prot.median_rho.max(), f'{prot.median_rho.max():.2f}', 'individual_refit/same_function_summary.tsv')
rec('XML73', 'Results XML73', 'DNA-RNA median rho, min', dr.median_rho.min(), f'{dr.median_rho.min():.2f}', 'same_function_summary.tsv')
rec('XML73', 'Results XML73', 'DNA-RNA median rho, max', dr.median_rho.max(), f'{dr.median_rho.max():.2f}', 'same_function_summary.tsv')
vp = ver[ver.layer_b == 'microbial_metaP']
rec('XML73', 'Results XML73', 'common-eligible protein median rho, min', vp.median_refit_common.min(), f'{vp.median_refit_common.min():.2f}', 'tables/<v>/common_feature_verification.tsv')
rec('XML73', 'Results XML73', 'common-eligible protein median rho, max', vp.median_refit_common.max(), f'{vp.median_refit_common.max():.2f}', 'tables/<v>/common_feature_verification.tsv')
rec('ST4', 'Supplementary Table 4', 'eligible tests', ntests, f'{ntests:,}', 'individual_refit/crossomics_same_function_correlations.tsv')
rec('ST4', 'Supplementary Table 4', 'tests sheet name', f'All_{ntests}_tests', f'All_{ntests}_tests', 'Supplementary_Table_4.xlsx')
for r in sm.itertuples():
    rec('ST4 Nine_family_summary', 'Supplementary Table 4', f'{r.level} {r.layer_a}-{r.layer_b} features', r.features, str(r.features), 'same_function_summary.tsv')
    rec('ST4 Nine_family_summary', 'Supplementary Table 4', f'{r.level} {r.layer_a}-{r.layer_b} median rho', r.median_rho, f'{r.median_rho:.3f}', 'same_function_summary.tsv')
    rec('ST4 Nine_family_summary', 'Supplementary Table 4', f'{r.level} {r.layer_a}-{r.layer_b} significant', r.significant, str(r.significant), 'same_function_summary.tsv')
for r in ver.itertuples():
    for f in ['n_common', 'n_new', 'n_lost']:
        rec('ST4 Eligibility_common', 'Supplementary Table 4', f'{r.level} {r.layer_a}-{r.layer_b} {f}', getattr(r, f), str(getattr(r, f)), 'common_feature_verification.tsv')
    for f in ['median_refit_common', 'median_baseline_common', 'median_refit_new']:
        v = getattr(r, f); rec('ST4 Eligibility_common', 'Supplementary Table 4', f'{r.level} {r.layer_a}-{r.layer_b} {f}', v, '' if pd.isna(v) else f'{v:.3f}', 'common_feature_verification.tsv')
# Figure 3f labels
fv = json.loads((S/'fig3/FIGURE_VALIDATION.json').read_text())
for s in fv['summary_values']:
    rec('Fig 3f', f"Figure 3f {s['space']} {s['tier']} {s['cls']}", 'median r label', s['median_r'], f"{s['median_r']:.2f}", 'hstool_precision_refit/pairwise_native_features.tsv')
    rec('Fig 3f', f"Figure 3f {s['space']} {s['tier']} {s['cls']}", 'median CV label (%)', s['median_cv'], f"{s['median_cv']:.0f}", 'hstool_precision_refit/pairwise_native_features.tsv')
    rec('Fig 3f', f"Figure 3f {s['space']} {s['tier']} {s['cls']}", 'median jointly quantified features label', s['median_features'], f"{s['median_features']:,.0f}", 'hstool_precision_refit/pairwise_native_features.tsv')
# Figure 4i and Supplementary Figure 11 a,b labels
pv = rd(H/'figures/all_plotted_crossomics_values.tsv')
for r in pv.itertuples():
    where = {'panel_i': 'Figure 4i', 'Supplementary_crossomics_PREDICT': 'Supplementary Fig. 11a', 'Supplementary_crossomics_CAPSCAN': 'Supplementary Fig. 11b (and Fig. 5g refit values, D2)'}[r.panel]
    rec(where.split(' (')[0], where, f'{r.cohort} {r.level} {r.layer_a}-{r.layer_b} cell label', r.rho, r.label, 'cross_refit native/functional_correlations.tsv')
sv = rd(H/'figures/supplementary_plotted_values.tsv')
rec('Supp Fig 11c-f', 'Supplementary Fig. 11c-f', 'sha256 of plotted native values', hashlib.sha256(sv.to_csv(sep='\t', index=False).encode()).hexdigest(), '', 'supplementary_plotted_values.tsv')
# Figure 4h / Supplementary Figure 10 (no printed numbers; record displayed values)
fz = np.load(ES/'focused_values.npz'); lz = np.load(ES/'legacy_85_values.npz')
for i, f in enumerate(fz['features']):
    rec('Fig 4h/SF10a values', 'Figure 4h, Supplementary Fig. 10a', f'positive patients {f}', int(np.isfinite(fz['z'][i]).sum()), '', 'focused_values.npz')
lf = rd(ES/'legacy_features.tsv')
for i, f in enumerate(lf.feature_id if 'feature_id' in lf else lf.iloc[:, 0]):
    rec('SF10b values', 'Supplementary Fig. 10b', f'positive patients {f}', int(np.isfinite(lz['z'][i]).sum()), '', 'legacy_85_values.npz')
pd.DataFrame(rows).to_csv(N/f'numbers_{V}.tsv', sep='\t', index=False)
(N/f'TABLES_{V}.json').write_text(json.dumps(dict(variant=V, controls={k: list(v) for k, v in cmp.items()}, st3=str(st3), st3_sha256=sha(st3), st4=str(st4), st4_sha256=sha(st4), numbers=len(rows)), indent=1)+'\n')
print(json.dumps({k: list(v) for k, v in cmp.items()}), len(rows), 'numbers', flush=True)
