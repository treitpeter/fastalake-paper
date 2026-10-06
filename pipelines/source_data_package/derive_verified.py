"""Round 3: decide which package files are replaced by their verified-row versions, derive the
files that need a format step, and list the additions for v20.

Usage: derive_verified.py WORK_DIR [--negative-control]
Inputs (read only): round-2 SELECTION.tsv and MANIFEST_build.tsv, work/INDEX_variants.tsv
(index_variants.py), the cr003 and iteration-1 replay/verified50 trees, the joint-search,
Supplementary Fig. 10 and Fig. 5g outputs of round 2, and the CapScan raw-header audit.
Every replacement is accepted only if the replay counterpart equals the shipped source
(bytes, decompressed bytes for .gz, or cell values for .xlsx) -- the positive control.
--negative-control swaps in the verified file as the 'replay' side, which must fail.
Writes WORK_DIR/derived/{OVERRIDES.tsv, ADDITIONS.tsv, CONTROLS.tsv, files/*}."""
import csv, gzip, hashlib, io, json, os, sys, zipfile
import numpy as np, pandas as pd, scipy.sparse as sp
sys.path.insert(0, os.path.dirname(__file__))
import selection as S

C = S.B + 'claude_review_2026-10-03/'
R2 = C + 'iteration2/source_data/work/'
ROOTS = {
    'downstream': C + 'iteration1/fix_eggnog_native/downstream/',
    'profiles': C + 'iteration1/fix_eggnog_native/profiles/',
    'shadow': C + 'iteration2/cr003_followup/shadow/',
    'tables': C + 'iteration2/cr003_followup/tables/',
    'nav': C + 'iteration2/cr003_followup/source_data/Native_annotation_verification/',
}
NAV_SUB = {'replay': 'previous', 'verified50': 'corrected'}
PREF = ['nav', 'downstream', 'profiles', 'shadow', 'tables']


def vpath(tree, variant, rel):
    if tree == 'nav':
        return ROOTS['nav'] + NAV_SUB[variant] + '/' + rel
    return ROOTS[tree] + variant + '/' + rel


def sha(b):
    return hashlib.sha256(b).hexdigest()


def src_bytes(row):
    return __import__('build_package').read_member(row['source'], row['member'])


def xlsx_cells(data, sheets=None):
    import openpyxl
    wb = openpyxl.load_workbook(io.BytesIO(data), read_only=True)
    out = {}
    for ws in wb.worksheets:
        if sheets and ws.title not in sheets:
            continue
        out[ws.title] = [tuple(r) for r in ws.iter_rows(values_only=True)]
    return out


def cells_close(a, b, rtol=1e-12):
    """Workbooks equal sheet by sheet; numbers within relative tolerance 1e-12 (last-digit float differences)."""
    if set(a) != set(b):
        return False
    for k in a:
        if len(a[k]) != len(b[k]):
            return False
        for x, y in zip(a[k], b[k]):
            if len(x) != len(y):
                return False
            for u, v in zip(x, y):
                if isinstance(u, (int, float)) and isinstance(v, (int, float)) and not isinstance(u, bool):
                    if not np.isclose(u, v, rtol=rtol, atol=0):
                        return False
                elif u != v:
                    return False
    return True


def main():
    work = sys.argv[1]; neg = '--negative-control' in sys.argv
    D = os.path.join(work, 'derived'); F = os.path.join(D, 'files'); os.makedirs(F, exist_ok=True)
    sel = {r['package_path']: r for r in csv.DictReader(open(R2 + 'SELECTION.tsv'), delimiter='\t')}
    man = {r['package_path']: r for r in csv.DictReader(open(R2 + 'full/MANIFEST_build.tsv'), delimiter='\t')}
    ix = pd.read_csv(os.path.join(work, 'INDEX_variants.tsv'), sep='\t')
    rep = ix[ix.variant == 'replay']; ver = ix[ix.variant == 'verified50']
    vsha = {(r.tree, r.relpath): r.sha256 for r in ver.itertuples()}
    over, controls = [], []
    rp, vp = ('verified50', 'replay') if neg else ('replay', 'verified50')

    def ctrl(path, kind, ok, detail):
        controls.append([path, kind, 'PASS' if ok else 'FAIL', detail])
        return ok

    # 1. byte-identical replay counterparts (positive control: replay bytes == shipped source bytes)
    for pp, m in man.items():
        hits = rep[rep.sha256 == m['origin_sha256']]
        if hits.empty:
            continue
        hits = sorted(hits.itertuples(), key=lambda h: PREF.index(h.tree))
        h = hits[0]
        replay_file = vpath(h.tree, rp, h.relpath); ver_file = vpath(h.tree, vp, h.relpath)
        ok = ctrl(pp, 'replay bytes equal shipped source', sha(open(replay_file, 'rb').read()) == m['origin_sha256'], f'{h.tree}:{h.relpath}')
        if not ok:
            continue
        vb = open(ver_file, 'rb').read()
        if sha(vb) != m['origin_sha256']:
            over.append([pp, ver_file, 'verified-row version (threshold 0.5); replay equals the shipped file byte for byte'])
    # 2. gz files: compare decompressed content
    for pp, m in man.items():
        if not pp.endswith('.gz') or any(o[0] == pp for o in over):
            continue
        name = pp.split('/')[-1]
        cand = rep[rep.relpath.str.endswith(name)]
        if cand.empty:
            continue
        orig = gzip.decompress(src_bytes(sel[pp]))
        for h in sorted(cand.itertuples(), key=lambda h: PREF.index(h.tree)):
            rf = vpath(h.tree, rp, h.relpath)
            if gzip.decompress(open(rf, 'rb').read()) == orig:
                ctrl(pp, 'replay decompressed equals shipped decompressed', True, f'{h.tree}:{h.relpath}')
                vf = vpath(h.tree, vp, h.relpath)
                if gzip.decompress(open(vf, 'rb').read()) != orig:
                    over.append([pp, vf, 'verified-row version (threshold 0.5); replay equals the shipped file after decompression'])
                break
        else:
            ctrl(pp, 'replay decompressed equals shipped decompressed', False, 'no decompressed match')
    # 3. HSTOOL profile exports rebuilt from the verified sparse profiles
    for name in ['KO_all', 'KO_other_unresolved']:
        pp = f'HSTOOL_KO_refit/revised_profiles/{name}.tsv.gz'
        def export(variant):
            base = ROOTS['profiles'] + variant + '/HSTOOL_noclust/'
            z = np.load(base + name + '.npz'); ax = json.load(open(base + name + '_axes.json'))
            m = sp.csr_matrix((z['data'], z['indices'], z['indptr']), shape=tuple(z['shape']))
            return pd.DataFrame(m.toarray().T, index=pd.Index(ax['features'], name='term'), columns=ax['samples']).to_csv(sep='\t')
        orig = gzip.decompress(src_bytes(sel[pp])).decode()
        if ctrl(pp, 'replay profile export equals shipped decompressed text', export(rp) == orig, 'profiles/HSTOOL_noclust/' + name + '.npz'):
            buf = io.BytesIO()
            with gzip.GzipFile(fileobj=buf, mode='wb', mtime=0, filename='') as g:
                g.write(export(vp).encode())
            out = os.path.join(F, f'HSTOOL_{name}.tsv.gz'); open(out, 'wb').write(buf.getvalue())
            over.append([pp, out, 'rebuilt from the verified-row sparse profile (same export as the shipped file; replay export equals it exactly)'])
    # 4. Supplementary Tables 3 and 4, Table4 eligibility rows (cell-level positive control)
    for pp, rel, sheets in [('Supplementary_Tables/Supplementary_Table_3.xlsx', 'Supplementary_Table_3.xlsx', None),
                            ('Supplementary_Tables/Supplementary_Table_4.xlsx', 'Supplementary_Table_4.xlsx', ['All_30632_tests', 'Nine_family_summary', 'All_eligibility', 'Eligibility_common'])]:
        a = xlsx_cells(src_bytes(sel[pp]), sheets); b = xlsx_cells(open(ROOTS['tables'] + rp + '/' + rel, 'rb').read(), sheets)
        same = cells_close(a, b)
        if ctrl(pp, 'replay workbook equals shipped workbook cell for cell, numbers within rtol 1e-12' + (' (data sheets)' if sheets else ''), same, rel):
            over.append([pp, ROOTS['tables'] + vp + '/' + rel, 'verified-row workbook'])
    pp = 'Supplementary_Tables/Table4_support/eligibility_matched_correlations.tsv'
    rg = gzip.decompress(open(ROOTS['tables'] + rp + '/eligibility_matched_correlations.tsv.gz', 'rb').read())
    if ctrl(pp, 'replay rows equal shipped rows', rg == src_bytes(sel[pp]), 'tables/eligibility_matched_correlations.tsv.gz'):
        out = os.path.join(F, 'eligibility_matched_correlations.tsv')
        open(out, 'wb').write(gzip.decompress(open(ROOTS['tables'] + vp + '/eligibility_matched_correlations.tsv.gz', 'rb').read()))
        over.append([pp, out, 'verified-row eligibility comparison'])
    # 5. derived threshold-change tables stay header-only under both variants
    for v in (rp, vp):
        f = pd.read_csv(ROOTS['downstream'] + v + '/B/cross_refit/functional_correlations.tsv', sep='\t')
        k = ['cohort', 'scope', 'level', 'layer_a', 'layer_b']
        a = f[f.prevalence == 0].merge(f[f.prevalence == 0.1], on=k, suffixes=('_all', '_p10'))
        n1 = int(((a.q_all < 0.05) != (a.q_p10 < 0.05)).sum())
        mt = f[f.scope.str.contains('matched')].merge(f[f.scope.str.contains('mixed')], on=['cohort', 'level', 'layer_a', 'layer_b', 'prevalence'], suffixes=('_m', '_x'))
        n2 = int(((mt.q_m < 0.05) != (mt.q_x < 0.05)).sum())
        ctrl(f'Biology_refit/PREVALENCE_THRESHOLD_CHANGES.tsv;POPULATION_THRESHOLD_CHANGES.tsv ({v})', 'classification changes recomputed', n1 == 0 and n2 == 0, f'prevalence {n1} of {len(a)}, population {n2} of {len(mt)}')
    # 5b. EDS selection-review tables (no saved builder): rebuilt from the screen and EDS association tables.
    #     Rule found by reconstruction: all native unadjusted KO links, joined with each feature's native unadjusted
    #     EDS association; inverse pairs = rho < 0, q_augmented < 0.05, both EDS q < 0.05, opposite EDS signs, sorted by rho.
    #     Metabolite-level columns (assay annotation, frozen-heatmap membership) do not depend on the protein annotation.
    sr = 'EDS_focused_20261003/selection_review/'
    shipL = pd.read_csv(io.BytesIO(src_bytes(sel[sr + 'all_native_unadjusted_KO_links_with_EDS.tsv'])), sep='\t')
    shipI_bytes = src_bytes(sel[sr + 'inverse_pairs_with_opposing_EDS.tsv']); shipL_bytes = src_bytes(sel[sr + 'all_native_unadjusted_KO_links_with_EDS.tsv'])
    legacy = set(pd.read_csv(io.BytesIO(src_bytes(sel['EDS_focused_20261003/source_data/legacy_features.tsv'])), sep='\t').feature_id)
    ann = shipL.drop_duplicates('metabolite').set_index('metabolite').assay_annotation
    def selection_review(variant):
        base = ROOTS['downstream'] + variant + '/B/eds_refit/'
        a = pd.read_csv(base + 'augmented_pair_screen.tsv', sep='\t'); e = pd.read_csv(base + 'all_EDS_associations.tsv', sep='\t')
        s_ = a[(a.layer == 'KO') & (a.normalization == 'native') & (a.model == 'unadjusted')].reset_index(drop=True)
        eu = e[(e.normalization == 'native') & (e.model == 'unadjusted')]
        met = eu[eu.layer == 'Metabolite'].set_index('feature'); ko = eu[eu.layer == 'KO'].set_index('feature')
        r = s_.copy()
        r['assay_annotation'] = r.metabolite.map(ann).fillna(r['name'])
        r['metabolite_in_frozen_heatmap'] = r.metabolite.isin(legacy); r['KO_in_frozen_heatmap'] = r.protein_feature.isin(legacy)
        r['metabolite_EDS_rho'] = r.metabolite.map(met.rho); r['metabolite_EDS_q'] = r.metabolite.map(met.q)
        r['protein_feature_EDS_rho'] = r.protein_feature.map(ko.rho); r['protein_feature_EDS_q'] = r.protein_feature.map(ko.q)
        inv = r[(r.rho < 0) & (r.q_augmented < 0.05) & (r.metabolite_EDS_q < 0.05) & (r.protein_feature_EDS_q < 0.05)
                & (np.sign(r.metabolite_EDS_rho) == -np.sign(r.protein_feature_EDS_rho))].sort_values('rho').reset_index(drop=True)
        return r, inv
    rL, rI = selection_review(rp)
    tL = rL.to_csv(sep='\t', index=False).encode(); tI = rI.to_csv(sep='\t', index=False).encode()
    okL = ctrl(sr + 'all_native_unadjusted_KO_links_with_EDS.tsv', 'replay rebuild equals shipped bytes', tL == shipL_bytes, f'{len(rL)} rows')
    okI = ctrl(sr + 'inverse_pairs_with_opposing_EDS.tsv', 'replay rebuild equals shipped bytes', tI == shipI_bytes, f'{len(rI)} rows')
    if okL and okI:
        vL, vI = selection_review(vp)
        for nm, df in (('all_native_unadjusted_KO_links_with_EDS.tsv', vL), ('inverse_pairs_with_opposing_EDS.tsv', vI)):
            out = os.path.join(F, nm); df.to_csv(out, sep='\t', index=False)
            over.append([sr + nm, out, f'rebuilt from the verified-row screen and EDS association tables ({len(df)} rows); replay rebuild equals the shipped file byte for byte'])
    # resolve duplicates: one override per path
    seen = {}
    for o in over:
        seen.setdefault(o[0], o)
    over = list(seen.values())
    with open(os.path.join(D, 'OVERRIDES.tsv'), 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t', lineterminator='\n'); w.writerow(['package_path', 'source_file', 'note']); w.writerows(sorted(over))
    # 6. additions
    J = C + 'iteration2/joint_search/'
    add = [
        ('Supplementary_Figure_07/per_acquisition_confidence/identity_per_acquisition.tsv', J + 'out/compare/identity_per_acquisition.tsv', 'SF7a;Results joint search'),
        ('Supplementary_Figure_07/per_acquisition_confidence/identity_summary.tsv', J + 'out/compare/identity_summary.tsv', 'SF7a;Results joint search'),
        ('Supplementary_Figure_07/per_acquisition_confidence/completeness_all.tsv', J + 'out/compare/completeness_all.tsv', 'Results joint search'),
        ('Supplementary_Figure_07/per_acquisition_confidence/microbial_group_pairs.tsv', J + 'out/compare/microbial_group_pairs.tsv', 'SF7c,d'),
        ('Supplementary_Figure_07/per_acquisition_confidence/microbial_group_pair_summary.tsv', J + 'out/compare/microbial_group_pair_summary.tsv', 'SF7c,d;Results joint search'),
        ('Supplementary_Figure_07/Supplementary_Figure7_plotted_values.tsv', J + 'figure/Supplementary_Figure7_plotted_values.tsv', 'SF7'),
        ('Supplementary_Figure_07/Supplementary_Figure7a_per_acquisition_column.tsv', J + 'figure/Supplementary_Figure7a_new_column_rows.tsv', 'SF7a'),
        ('Supplementary_Figure_07/Supplementary_Figure7cd_pair_rows.tsv', J + 'figure/Supplementary_Figure7cd_pair_rows.tsv', 'SF7c,d'),
        ('Supplementary_Figure_10/Supplementary_Figure10_plotted_values.tsv', C + 'iteration2/consolidate/suppfig10_merged/out/SupplementaryFigure10_plotted_values.tsv', 'SF10'),
        ('Figure5/Figure5g_plotted_values.tsv', C + 'iteration2/figures/figure5/Figure5g_plotted_values.tsv', 'Fig5g'),
        ('Annotation_verification/native_row_verification.tsv.gz', ROOTS['nav'] + 'verification/native_row_verification.tsv.gz', 'Methods annotation check (9,368,708 rows)'),
        ('Annotation_verification/verification_summary.json', ROOTS['nav'] + 'verification/VERIFICATION_SUMMARY.json', 'Methods annotation check counts'),
        ('Annotation_verification/independent_mapper_concordance.json', ROOTS['nav'] + 'verification/INDEPENDENT_MAPPER_CONCORDANCE.json', 'annotation check against eggNOG-mapper'),
        ('Annotation_verification/previous_vs_verified_all_numbers.tsv', ROOTS['nav'] + 'comparison/OLD_VS_NEW_ALL_NUMBERS.tsv', 'values before and after the check'),
        ('Annotation_verification/previous_vs_verified_printed_numbers.tsv', ROOTS['nav'] + 'comparison/OLD_VS_NEW_PRINTED_CHANGES.tsv', 'printed values before and after the check'),
        ('Annotation_verification/cohort_profile_change_summary.tsv', ROOTS['nav'] + 'comparison/COHORT_PROFILE_CHANGE_SUMMARY.tsv', 'functional signal reallocated per acquisition'),
        ('Annotation_verification/heatmap_displayed_value_changes.tsv', ROOTS['nav'] + 'comparison/HEATMAP_DISPLAYED_VALUE_CHANGES.tsv', 'Fig4h/SF10 displayed changes'),
        ('Annotation_verification/selected_KO_cohort_totals.tsv', ROOTS['nav'] + 'comparison/selected_KO_cohort_totals.tsv', 'selected KO totals'),
    ]
    # authored: CapScan specimen-number discrepancies (R2-014) and joint-search arm definitions
    au = pd.read_csv(S.J + 'cohort_reconciliation/raw_identity_audit/all_RAW_header_names.tsv', sep='\t')
    x = au[au.header_specimen_number_changed == True].copy()
    num = lambda s: s.split('_T760_')[1].split('_')[0]
    disc = pd.DataFrame({'acquisition': x['sample'], 'in_analysis_roster': x['current_331'].map({True: 'yes', False: 'no'}),
                         'roster_participant': x['subject'].astype('Int64'),
                         'specimen_number_in_file_name': x['sample'].map(num),
                         'specimen_number_in_instrument_record': x['header_original_names'].map(num),
                         'status': 'unresolved: analysis labels follow the frozen roster'})
    ctrl('Supplementary_Data_1/CapScan_specimen_number_discrepancies.tsv', 'count equals the printed sixteen and fourteen',
         len(disc) == 16 and int((disc.in_analysis_roster == 'yes').sum()) == 14 and len(au) == 368, f'{len(disc)} acquisitions, {(disc.in_analysis_roster == "yes").sum()} in roster, {len(au)} audited')
    p = os.path.join(F, 'CapScan_specimen_number_discrepancies.tsv'); disc.to_csv(p, sep='\t', index=False)
    add.append(('Supplementary_Data_1/CapScan_specimen_number_discrepancies.tsv', p, 'Results XML64;SM1 (sixteen acquisitions, fourteen in the roster)'))
    arms = pd.DataFrame([
        ('original', 'acquisition-specific FastaLake searches; peptide q estimated by Sage within each acquisition'),
        ('original_fileq', 'acquisition-specific searches re-scored with the per-acquisition estimator used for joint_fileq (sensitivity)'),
        ('union_individual', 'separate search of each acquisition against the common 153,267-protein FASTA'),
        ('joint_corrected', 'joint search of all 58 acquisitions (corrected Sage build); study-level peptide q estimated once across the study'),
        ('joint_sageq_repro', 'joint_corrected re-extracted from the same PSMs with study-level q (reproduces joint_corrected; positive control)'),
        ('joint_fileq', 'joint search PSMs with peptide q re-estimated within each acquisition (per-acquisition confidence, the primary joint result)'),
    ], columns=['arm', 'definition'])
    p = os.path.join(F, 'arm_definitions.tsv'); arms.to_csv(p, sep='\t', index=False)
    add.append(('Supplementary_Figure_07/per_acquisition_confidence/arm_definitions.tsv', p, 'SF7 arm definitions'))
    for pp, src, use in add:
        if not os.path.exists(src):
            ctrl(pp, 'addition source exists', False, src)
    with open(os.path.join(D, 'ADDITIONS.tsv'), 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t', lineterminator='\n'); w.writerow(['package_path', 'source_file', 'use']); w.writerows(add)
    with open(os.path.join(D, 'CONTROLS.tsv'), 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t', lineterminator='\n'); w.writerow(['item', 'check', 'result', 'detail']); w.writerows(controls)
    nfail = sum(1 for c in controls if c[1] != 'FAIL' and c[2] == 'FAIL')
    print(json.dumps({'overrides': len(over), 'additions': len(add), 'controls': len(controls), 'control_failures': nfail, 'negative_control_mode': neg}))
    sys.exit(2 if nfail else 0)


if __name__ == '__main__':
    main()
