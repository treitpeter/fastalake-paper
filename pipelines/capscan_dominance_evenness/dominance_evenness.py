"""C4 (JAG55): per-acquisition species dominance (top-1 share) and Pielou evenness,
small-intestinal fluid (regions 1-4) versus stool, CapScan frozen 331 roster.

Metric definitions are copied from the pre-refit diagnostic
  analysis_followup/dominance_evenness_2026-09-30/analyze.py  (metrics(), named-species rule BAD,
  species summed across KOs; paired() exact sign-flip enumeration)
and from the annotation-repair arm
  annotation_repair/manuscript_benchmark_repaired/arm.py cap_species() (single-consistent KO rows,
  species-rank microbial lineage).
Only change: metrics are evaluated per acquisition, then summarised as participant medians.

Sources
  a33       A33 CAPSCAN__simplified Figure5/function_species.npz: the species annotation behind v26
            Fig. 5e/5f, Supplementary Fig. 12 and the 29.4% sentence (C2 fig5e_f control trace).
  repaired  annotation_repair evidence_repaired/full/CAPSCAN_noclust KO_evidence_signal.npz
            (single-consistent KO x species-rank microbial taxonomy).

Sanity checks (must pass before any result is written)
  S1 analytical metric cases;
  S2 legacy aggregation (total-LFQ normalised, specimen mean, participant-region mean) on the a33
     inputs reproduces the pre-refit community_metrics.tsv (simplified, all_named) and the
     20.8% / 11.8% and 0.589 / 0.653 summary values;
  S3 the same aggregation on the repaired inputs reproduces CAPSCAN_species_evenness.tsv top1/pielou.
"""
from pathlib import Path
import sys, json, gzip, hashlib
from collections import defaultdict
import numpy as np, pandas as pd
from scipy import sparse
from scipy.stats import wilcoxon

OUT = Path(__file__).resolve().parent
ROOT = Path('LAB_ROOT')
A33 = ROOT / 'MM_review_analysis/A33_4_6_peptide_centric_downstream_2026-09-19/results/CAPSCAN__simplified/source_data/Figure5'
JI = ROOT / 'review/v15_integration/jess_integration_2026-09-29/analysis_followup'
PRE = JI / 'dominance_evenness_2026-09-30'
REP = JI / 'github_workflow_reconciliation_2026-09-30/resumed_audit/annotation_repair'
EVD = REP / 'evidence_repaired/full/CAPSCAN_noclust'
REP_ARM = REP / 'manuscript_benchmark_repaired/arms/CAPSCAN_noclust'
INPUTS = {
    'design_331': A33 / 'design.tsv',
    'a33_function_species': A33 / 'function_species.npz',
    'a33_acquisitions': A33 / 'acquisition_composition_local_candidate.tsv',
    'repaired_KO_evidence_signal': EVD / 'KO_evidence_signal.npz',
    'repaired_KO_quantity_ids': EVD / 'KO_quantity_ids.json',
    'repaired_quantity_definitions': EVD / 'quantity_definitions.jsonl.gz',
    'repaired_sample_totals': EVD / 'sample_totals.npz',
    'prerefit_community_metrics': PRE / 'full/community_metrics.tsv',
    'repaired_species_evenness': REP_ARM / 'CAPSCAN_species_evenness.tsv',
}
REG = ['duodenum', 'jejunum', 'ileum', 'colon', 'stool']
FLUID = REG[:4]
BAD = ['Host', 'Mixed host/microbial', 'Conflicting species', 'No species annotation']  # pre-refit analyze.py


def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(1 << 23), b''):
            h.update(b)
    return h.hexdigest()


def metrics(v):
    """Copied from pre-refit analyze.py metrics(); v = species x columns."""
    v = np.asarray(v, float); tot = v.sum(0)
    p = np.divide(v, tot[None, :], out=np.zeros_like(v), where=tot[None, :] > 0)
    S = (p > 0).sum(0); logs = np.zeros_like(p); np.log(p, out=logs, where=p > 0)
    entropy = -(p * logs).sum(0)
    even = np.divide(entropy, np.log(np.maximum(S, 1)), out=np.full(len(tot), np.nan), where=S > 1)
    top1 = np.where(tot > 0, np.sort(p, axis=0)[-1], np.nan)
    return dict(top1=top1, pielou=even, richness=S.astype(float), resolved_mass=tot)


def signflip(v):
    """Exact two-sided paired sign-flip on the sum (pre-refit paired())."""
    v = np.asarray(v, float); v = v[np.isfinite(v)]; n = len(v)
    signs = (((np.arange(2 ** n, dtype=np.uint32)[:, None] >> np.arange(n)) & 1).astype(np.int8) * 2 - 1)
    return float(np.mean(np.abs(signs @ v) >= abs(v.sum()) - 1e-12))


def aggregate(v, meta, fields):
    g = meta.groupby(fields, sort=True).indices; keys = list(g)
    return np.stack([v[:, g[k]].mean(1) for k in keys], 1), pd.DataFrame(keys, columns=fields)


def legacy_cells(taxa_by_acq, meta, total):
    """Pre-refit aggregation: /total LFQ, specimen mean, participant-region mean (14 x 5)."""
    rel = taxa_by_acq / total[None, :]
    sv, sm = aggregate(rel, meta, ['participant', 'region', 'specimen_id'])
    cv, cm = aggregate(sv, sm, ['participant', 'region'])
    idx = pd.MultiIndex.from_frame(cm).get_indexer(pd.MultiIndex.from_product([range(1, 15), REG]))
    assert (idx >= 0).all()
    return cv[:, idx], cm.iloc[idx].reset_index(drop=True)


def load_a33(design):
    z = np.load(INPUTS['a33_function_species'], allow_pickle=True)
    vals, sp, samples = z['values'], z['species'], z['samples']
    named = ~np.isin(sp, BAD); nonhost = ~np.isin(sp, ['Host', 'Mixed host/microbial'])
    names, inv = np.unique(sp[named], return_inverse=True)
    taxa = np.zeros((len(names), vals.shape[1])); np.add.at(taxa, inv, vals[named])
    acq = pd.read_csv(INPUTS['a33_acquisitions'], sep='\t').set_index('sample_id').loc[samples]
    total = acq.total_LFQ.to_numpy(float)
    nonhost_mass = vals[nonhost].sum(0)
    return taxa, list(samples), total, nonhost_mass, names


def load_repaired():
    z = np.load(INPUTS['repaired_sample_totals']); samples = [str(s) for s in z['samples']]
    total = np.asarray(z['intensity'], float)
    mat = sparse.load_npz(INPUTS['repaired_KO_evidence_signal']).tocsr()
    ids = json.loads(Path(INPUTS['repaired_KO_quantity_ids']).read_text()); want = set(ids)
    defs = {}
    with gzip.open(INPUTS['repaired_quantity_definitions'], 'rt') as f:
        for line in f:
            r = json.loads(line)
            if r['quantity_id'] in want:
                defs[r['quantity_id']] = dict(function=r['function'], taxonomy=r['taxonomy'])
    rows = defaultdict(list); single = []
    for i, k in enumerate(ids):
        r = defs[k]
        if r['function']['state'] != 'single_consistent':
            continue
        single.append(i); t = r['taxonomy']; lin = t['supported_lineage']
        if t['rank'] == 'species' and lin and lin[0] in ('d__Bacteria', 'd__Archaea'):
            rows[lin[-1][3:]].append(i)
    names = np.array(sorted(rows))
    taxa = np.vstack([np.asarray(mat[rows[s]].sum(0)).ravel() for s in names])
    single_mass = np.asarray(mat[single].sum(0)).ravel()
    return taxa, samples, total, single_mass, names


def per_acquisition(source, taxa, samples, denom, design):
    m = metrics(taxa)
    d = design.set_index('sample_id').loc[samples].reset_index()
    out = d[['sample_id', 'participant', 'region', 'specimen_id']].copy()
    out['source'] = source
    out['compartment'] = np.where(out.region == 'stool', 'stool', 'fluid')
    for k in ['top1', 'pielou', 'richness', 'resolved_mass']:
        out[k] = m[k]
    out['species_resolved_fraction'] = np.divide(m['resolved_mass'], denom, out=np.full(len(denom), np.nan), where=denom > 0)
    # Roster acquisitions absent from the source matrix are kept, with no species signal.
    missing = sorted(set(design.sample_id) - set(samples))
    if missing:
        add = design[design.sample_id.isin(missing)][['sample_id', 'participant', 'region', 'specimen_id']].copy()
        add['source'] = source; add['compartment'] = np.where(add.region == 'stool', 'stool', 'fluid')
        for k in ['top1', 'pielou', 'species_resolved_fraction']:
            add[k] = np.nan
        add['richness'] = 0.0; add['resolved_mass'] = 0.0
        out = pd.concat([out, add], ignore_index=True)
    return out


def summarise(acq):
    rows = []; part_rows = []
    for (source, metric), _ in [((s, m), None) for s in acq.source.unique() for m in ['top1', 'pielou']]:
        a = acq[(acq.source == source)]
        good = a[np.isfinite(a[metric])]
        pm = good.groupby(['participant', 'compartment'])[metric].median().unstack()
        cnt = good.groupby(['participant', 'compartment'])[metric].size().unstack()
        both = pm.dropna(subset=['fluid', 'stool'])
        diff = (both.fluid - both.stool).to_numpy()
        w = wilcoxon(both.fluid, both.stool, alternative='two-sided', method='exact')
        for p_, r in pm.iterrows():
            part_rows.append(dict(source=source, metric=metric, participant=p_,
                                  fluid_median=r.get('fluid', np.nan), stool_median=r.get('stool', np.nan),
                                  fluid_acquisitions=int(cnt.loc[p_].get('fluid', 0) if pd.notna(cnt.loc[p_].get('fluid', np.nan)) else 0),
                                  stool_acquisitions=int(cnt.loc[p_].get('stool', 0) if pd.notna(cnt.loc[p_].get('stool', np.nan)) else 0)))
        rows.append(dict(
            source=source, metric=metric,
            acq_n_fluid=int((good.compartment == 'fluid').sum()), acq_n_stool=int((good.compartment == 'stool').sum()),
            acq_n_undefined=int(len(a) - len(good)),
            acq_median_fluid=float(good[good.compartment == 'fluid'][metric].median()),
            acq_median_stool=float(good[good.compartment == 'stool'][metric].median()),
            paired_participants=int(len(both)),
            participant_median_fluid=float(both.fluid.median()), participant_median_stool=float(both.stool.median()),
            median_paired_difference=float(np.median(diff)), mean_paired_difference=float(diff.mean()),
            participants_fluid_gt_stool=int((diff > 0).sum()),
            wilcoxon_W=float(w.statistic), wilcoxon_p_exact=float(w.pvalue),
            signflip_p_exact=signflip(diff), signflip_patterns=2 ** len(diff)))
    return pd.DataFrame(rows), pd.DataFrame(part_rows)


def main():
    hashes = {k: sha(p) for k, p in INPUTS.items()}
    design = pd.read_csv(INPUTS['design_331'], sep='\t')
    assert len(design) == 331 and design.sample_id.is_unique
    checks = {}

    # S1 analytical cases
    m = metrics(np.array([[1., 1.], [1., 0.], [1., 0.], [1., 0.]]))
    assert np.allclose(m['top1'], [.25, 1]) and np.isclose(m['pielou'][0], 1) and np.isnan(m['pielou'][1])
    assert signflip(np.ones(14)) == 2 / 2 ** 14
    checks['S1_analytical'] = 'PASS'

    # A33 source
    taxa, samples, total, nonhost, _ = load_a33(design)
    assert len(samples) == 330
    meta = design.set_index('sample_id').loc[samples].reset_index()
    # S2 reproduce pre-refit participant-region metrics and the 20.8/11.8, 0.589/0.653 values
    cv, cm = legacy_cells(taxa, meta, total); lm = metrics(cv)
    pre = pd.read_csv(INPUTS['prerefit_community_metrics'], sep='\t')
    pre = pre[(pre.rule == 'simplified') & (pre.scope == 'all_named')].set_index(['participant', 'region']).loc[pd.MultiIndex.from_frame(cm)]
    for k in ['top1', 'pielou']:
        assert np.allclose(lm[k], pre[k].to_numpy(), rtol=1e-9, atol=1e-12, equal_nan=True), k
    w = {k: pd.Series(lm[k], index=pd.MultiIndex.from_frame(cm)).unstack() for k in ['top1', 'pielou']}
    s2 = {k: dict(capsule_mean=float(w[k][FLUID].mean(1).mean()), stool=float(w[k]['stool'].mean())) for k in w}
    assert round(100 * s2['top1']['capsule_mean'], 1) == 20.8 and round(100 * s2['top1']['stool'], 1) == 11.8, s2
    assert round(s2['pielou']['capsule_mean'], 3) == 0.589 and round(s2['pielou']['stool'], 3) == 0.653, s2
    checks['S2_prerefit_reproduced'] = dict(status='PASS', cells=70, values=s2,
                                            note='participant means of (mean of 4 region cells) vs stool cell')
    a33 = per_acquisition('a33', taxa, samples, nonhost, design)

    # Repaired source
    rtaxa, rsamples, rtotal, rsingle, _ = load_repaired()
    pos = rtotal > 0
    rmeta = design.set_index('sample_id').loc[np.array(rsamples)[pos]].reset_index()
    cv, cm = legacy_cells(rtaxa[:, pos], rmeta, rtotal[pos]); lm = metrics(cv)
    ref = pd.read_csv(INPUTS['repaired_species_evenness'], sep='\t').set_index(['participant', 'region']).loc[pd.MultiIndex.from_frame(cm)]
    for k in ['top1', 'pielou']:
        assert np.allclose(lm[k], ref[k].to_numpy(), rtol=1e-8, atol=1e-12, equal_nan=True), k
    checks['S3_repaired_arm_reproduced'] = dict(status='PASS', cells=70, positive_acquisitions=int(pos.sum()), acquisitions=len(rsamples))
    rep = per_acquisition('repaired', rtaxa, rsamples, rsingle, design)

    acq = pd.concat([a33, rep], ignore_index=True)
    assert (acq.groupby('source').size() == 331).all()
    summary, part = summarise(acq)
    acq.to_csv(OUT / 'per_acquisition_metrics.tsv', sep='\t', index=False)
    part.to_csv(OUT / 'participant_medians.tsv', sep='\t', index=False)
    summary.to_csv(OUT / 'summary_tests.tsv', sep='\t', index=False)
    roster = design.groupby(['participant', 'region']).size().unstack(fill_value=0)
    (OUT / 'RUN.json').write_text(json.dumps(dict(
        status='PASS', checks=checks, input_sha256=hashes, code_sha256=sha(__file__),
        roster_acquisitions=331, roster_participants=int(design.participant.nunique()),
        participants_with_fluid_and_stool=int(((roster[FLUID].sum(1) > 0) & (roster['stool'] > 0)).sum()),
        definitions=dict(
            a33='named species = function_species rows whose species label is not Host / Mixed host/microbial / Conflicting species / No species annotation, summed over KOs; resolved fraction denominator = non-host KO allocation',
            repaired='species = single-consistent KO evidence rows with species-rank Bacteria/Archaea supported lineage, summed over KOs; resolved fraction denominator = all single-consistent KO signal',
            top1='largest species share of species-resolved signal in the acquisition',
            pielou='Shannon entropy / ln(observed species); undefined for <2 species',
            test='participant median per compartment; paired over participants with both; Wilcoxon signed-rank exact and exact sign-flip on the sum of differences (2^n patterns)')),
        indent=2, default=float) + '\n')
    print(summary.to_string(index=False))
    print(json.dumps(checks, indent=1, default=float))


if __name__ == '__main__':
    main()
