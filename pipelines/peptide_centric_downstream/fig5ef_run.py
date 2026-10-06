"""Fig. 5e/5f and Supp. Fig. 12c,d (A33 peptide-centric lineage) on a given roster.

Executes, for each grouping rule, exactly the blocks of the A33/A19 producers that the v25 numbers
come from, with the frozen modules imported unmodified:
  * revised/derive_capscan_pc.py   -> meta + function_species.npz construction (block copied verbatim)
  * frozen summarize_biology.joint_summary (imported, unmodified)  -> contributor tables, SF12 numbers
  * frozen derive_capscan_extras.py fixed 24-function mosaic block (copied verbatim; range(1,16) unchanged)
then the A33 frozen collector.py capsule/stool matching (copied verbatim, minus its PREDICT
statistical-endpoint alignment and receipt checks, which do not feed these numbers) and
render/make_leader_tables.py (copied; its `assert len(t)==14` replaced by an assertion against the
expected participant count passed on the command line).
Not run: host composition, target tests, BSH, geometry (other branches own those numbers); the
local-candidate Human/Microbial/Mixed/Unknown columns are therefore not written to meta.
usage: fig5ef_run.py <gated_joint_dir> <grouping_cohorts_dir> <design.tsv> <outdir> <expected_participants>
"""
import importlib.util, json, os, sys
from collections import defaultdict, Counter
from itertools import product
from pathlib import Path
import numpy as np
import pandas as pd

A19F = Path('LAB_ROOT/MM_review_analysis/A19_4_6_biology_robustness/frozen')
FIG = Path('LAB_ROOT/FIGURATOR')
sys.path.insert(0, str(A19F))
GATED, GROUPS, DESIGN, OUTDIR, NPART = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]), Path(sys.argv[4]), int(sys.argv[5])
OUTDIR.mkdir(parents=True, exist_ok=False)
# frozen basis.py (imported by summarize_biology) needs these and creates scratch dirs under MM_BIO_OUT
os.environ.update(MM_BIO_OUT=str(OUTDIR/'mm_bio_scratch'), MM_BIO_ROOT='LAB_ROOT')
import common
common.RUN = GROUPS.parent
UNASSIGNED = 'FunctionallyUnassigned'
REGIONS = common.REGIONS


def write(frame, path, index=False):
    common.write(frame, path, index)


def figure5(rule):
    common.METHOD = rule
    out = OUTDIR / f'CAPSCAN__{rule}' / 'source_data' / 'Figure5'; out.mkdir(parents=True)
    # ---- derive_capscan_pc.main(), meta and function_species blocks ----
    matrix, annot, summary = common.load('CAPSCAN_noclust')
    design = pd.read_csv(DESIGN, sep='\t', dtype={'participant': str})
    total_all = matrix.sum()
    samples = total_all[total_all > 0].index
    matrix = matrix[samples]
    meta = design.set_index('sample_id').loc[samples].rename_axis('sample_id').reset_index()
    meta['total_LFQ'] = total_all.loc[samples].to_numpy()
    write(meta, out / 'acquisition_composition_local_candidate.tsv')
    species_phyla, _ = common.taxonomy()
    (out / 'species_phyla.json').write_text(json.dumps(species_phyla, indent=2))
    joint = defaultdict(lambda: np.zeros(len(samples)))
    stem = 'KO_species' if rule == 'simplified' else 'KO_species_original'
    with np.load(GATED / f'CAPSCAN_{stem}.npz', allow_pickle=False) as zz:
        z = dict(ko=zz['ko'], species=zz['species'], values=zz['values'], samples=zz['samples'])
    order = pd.Index(z['samples']).get_indexer(np.asarray(samples))
    assert (order >= 0).all(), 'the peptide-centric matrix is missing acquisitions'
    for ko, sp, row in zip(z['ko'], z['species'], z['values']):
        vv = row[order]
        if ko == UNASSIGNED:
            continue
        joint[ko, sp] += vv
    keys = sorted(joint)
    values = np.array([joint[k] for k in keys])
    np.savez_compressed(out / 'function_species.npz', ko=np.array([k[0] for k in keys]), species=np.array([k[1] for k in keys]),
                        values=values, samples=np.array(samples, dtype=str))
    del z, joint, values
    # ---- frozen summarize_biology.joint_summary, unmodified ----
    spec = importlib.util.spec_from_file_location('petr3_capscan_stats', A19F / 'summarize_biology.py')
    analysis = importlib.util.module_from_spec(spec); spec.loader.exec_module(analysis)
    analysis.D = out
    analysis.save = lambda frame, name: write(frame, out / (name + '.tsv'))
    contributors = analysis.joint_summary('function_species.npz', meta, '')
    # ---- derive_capscan_extras.main(), fixed 24-function mosaic block (verbatim) ----
    aggregate = analysis.aggregate
    meta = pd.read_csv(out / 'acquisition_composition_local_candidate.tsv', sep='\t')
    z = np.load(out / 'function_species.npz')
    meta = meta.set_index('sample_id').loc[z['samples']].reset_index()
    ko, species = z['ko'], z['species']
    values = z['values'] / meta.total_LFQ.to_numpy()
    nonhost = ~np.isin(species, ['Host', 'Mixed host/microbial'])
    resolved = nonhost & ~np.isin(species, ['No species annotation', 'Conflicting species'])
    vs, ms = aggregate(values, meta, ['participant', 'region', 'specimen_id'])
    vr, mr = aggregate(vs, ms, ['participant', 'region'])
    vp, mp = aggregate(vr, mr, ['participant'])
    oldrows = pd.read_csv(FIG / 'Figure5/Fig5/panel_who_in_you_sourcedata.tsv', sep='\t')
    selected = list(oldrows.ko)
    assert len(selected) == len(set(selected)) == 24
    phyla = json.loads((out / 'species_phyla.json').read_text())
    abundance = {k: float(vp[(ko == k) & nonhost].sum(axis=0).mean()) for k in selected}
    selected.sort(key=lambda k: (-abundance[k], k))
    rows, regional = [], []
    def cell(k, v):
        ix = np.flatnonzero((ko == k) & resolved)
        total = v[(ko == k) & nonhost].sum()
        rr = v[ix].sum()
        if rr <= 0:
            return dict(species='', phylum='', resolved_fraction=0 if total > 0 else np.nan, top_share=np.nan, status='No species-resolved signal')
        lead = sorted(ix, key=lambda i: (-v[i], species[i]))[0]
        return dict(species=str(species[lead]), phylum=phyla.get(str(species[lead]), 'Other phyla'), resolved_fraction=float(rr / total), top_share=float(v[lead] / rr), status='Resolved')
    lookup = {(r.participant, r.region): i for i, r in mr.iterrows()}
    for rank, k in enumerate(selected, 1):
        for j, r in mp.iterrows():
            rows.append(dict(row=rank, KO=k, participant=int(r.participant), function_share_total=vp[(ko == k) & nonhost, j].sum(), **cell(k, vp[:, j])))
        for pid, region in product(range(1, 16), REGIONS):
            j = lookup.get((pid, region))
            result = cell(k, vr[:, j]) if j is not None else dict(species='', phylum='', resolved_fraction=np.nan, top_share=np.nan, status='Not sampled')
            regional.append(dict(row=rank, KO=k, participant=pid, region=region, **result))
    write(pd.DataFrame(rows), out / 'mosaic_donors.tsv')
    write(pd.DataFrame(regional), out / 'mosaic_donor_regions.tsv')
    write(pd.DataFrame([dict(row=i+1, KO=k, mean_share_percent=100*abundance[k]) for i, k in enumerate(selected)]), out / 'mosaic_abundance.tsv')
    (out / 'SUMMARY.json').write_text(json.dumps(dict(status='PASS', rule=rule, n_acquisitions=len(samples),
        n_participants=int(meta.participant.nunique()), contributors=contributors), indent=2))
    return out


def read(p): return pd.read_csv(p, sep='\t')


# ---- A33 frozen/collector.py, region_state and matched comparisons (verbatim) ----
def collect():
    def aggregate(v, d, keys):
        groups = d.groupby(keys, sort=True).indices
        a = np.stack([v[:, ix].mean(axis=1) for ix in groups.values()], axis=1)
        return a, pd.DataFrame(list(groups), columns=keys)
    def region_state(folder, eligible):
        d = folder/'source_data/Figure5'; z = np.load(d/'function_species.npz')
        ac = read(d/'acquisition_composition_local_candidate.tsv').set_index('sample_id').loc[z['samples']].reset_index()
        v = z['values']/ac.total_LFQ.to_numpy(); ko = z['ko']; species = z['species']
        vs, sm = aggregate(v, ac, ['participant', 'region', 'specimen_id']); vr, cm = aggregate(vs, sm, ['participant', 'region'])
        nonhost = ~np.isin(species, ['Host', 'Mixed host/microbial'])
        resolved = nonhost & ~np.isin(species, ['No species annotation', 'Conflicting species', 'Unassigned'])
        old = read(d/'contributor_region_cells.tsv').set_index(['KO', 'participant', 'region'])
        rows = []
        for k in sorted(eligible):
            ix = np.flatnonzero((ko == k) & resolved); allix = np.flatnonzero((ko == k) & nonhost)
            totals = vr[allix].sum(axis=0); rr = vr[ix].sum(axis=0)
            for j, c in cm.iterrows():
                key = (k, c.participant, c.region); leader = ''; gap = np.nan; fraction = rr[j]/totals[j] if totals[j] > 0 else np.nan
                if rr[j] > 0:
                    order = sorted(ix, key=lambda i: (-vr[i, j], species[i])); leader = str(species[order[0]])
                    gap = (vr[order[0], j]-(vr[order[1], j] if len(order) > 1 else 0))/rr[j]
                    prior = old.loc[key]; assert prior.leading_species == leader
                    assert np.isclose(prior.gap, gap, rtol=1e-10, atol=1e-12)
                    assert np.isclose(prior.species_resolved_fraction, fraction, rtol=1e-10, atol=1e-12)
                else: assert key not in old.index
                rows.append(dict(KO=k, participant=int(c.participant), region=c.region,
                    nonhost_function_fraction=totals[j], resolved_fraction_total=rr[j], species_resolved_fraction=fraction,
                    leading_species=leader, leader_gap=gap, strict=bool(rr[j] > 0 and fraction >= .2 and gap >= .05),
                    state='resolved' if rr[j] > 0 else ('unresolved_positive_signal' if totals[j] > 0 else 'not_detected')))
        return pd.DataFrame(rows)
    out = OUTDIR/'results'; out.mkdir()
    folders = {method: OUTDIR/('CAPSCAN__'+method) for method in ['original', 'simplified']}
    elig = {method: set(read(p/'source_data/Figure5/contributor_cells.tsv').KO) for method, p in folders.items()}
    common_ = elig['original'] & elig['simplified']; union = elig['original'] | elig['simplified']
    assert common_
    write(pd.DataFrame([dict(KO=k, original_eligible=k in elig['original'], simplified_eligible=k in elig['simplified'], common_eligible=k in common_) for k in sorted(union)]), out/'contributor_eligibility.tsv')
    states = {method: region_state(p, union) for method, p in folders.items()}
    paired = states['original'].merge(states['simplified'], on=['KO', 'participant', 'region'], suffixes=('_original', '_simplified'), validate='1:1', how='outer', indicator=True)
    assert paired['_merge'].eq('both').all(); paired = paired.drop(columns='_merge')
    paired['common_eligible_function'] = paired.KO.isin(common_)
    paired['both_resolved'] = paired.state_original.eq('resolved') & paired.state_simplified.eq('resolved')
    paired['leader_changed'] = np.where(paired.both_resolved, paired.leading_species_original.ne(paired.leading_species_simplified), np.nan)
    paired['both_strict'] = paired.strict_original & paired.strict_simplified
    write(paired, out/'contributor_method_cells.tsv')
    comparisons = []; summaries = []
    for scope, cut in [('common_functions', paired.common_eligible_function), ('eligible_union_sensitivity', pd.Series(True, index=paired.index))]:
        pp = paired[cut].copy()
        stool = pp[pp.region.eq('stool')].drop(columns='region')
        cp = pp[pp.region.ne('stool')].merge(stool, on=['KO', 'participant'], suffixes=('_capsule', '_stool'), validate='m:1')
        for threshold, ok in [('resolved', cp.both_resolved_capsule & cp.both_resolved_stool), ('strict_both_methods', cp.both_strict_capsule & cp.both_strict_stool)]:
            cc = cp[ok].copy()
            for method in ['original', 'simplified']:
                cc['agreement_'+method] = cc['leading_species_'+method+'_capsule'].eq(cc['leading_species_'+method+'_stool'])
            dd = cc.groupby(['participant', 'region'])[['agreement_original', 'agreement_simplified']].mean().reset_index()
            dd['cells'] = cc.groupby(['participant', 'region']).size().to_numpy(); dd['scope'] = scope; dd['threshold'] = threshold
            comparisons.append(dd)
            for method in ['original', 'simplified']:
                dv = dd.groupby('participant')['agreement_'+method].mean()
                summaries.append(dict(endpoint='capsule_stool_agreement_'+method, scope=scope, threshold=threshold, cells=len(cc), participants=len(dv), participant_mean=float(dv.mean()), pooled_cell_mean=float(cc['agreement_'+method].mean())))
    write(pd.concat(comparisons, ignore_index=True), out/'matched_capsule_stool_donors.tsv')
    write(pd.DataFrame(summaries), out/'matched_leader_summary.tsv')
    # ---- render/make_leader_tables.py (copied; participant-count assertion parameterised) ----
    src = read(out/'matched_capsule_stool_donors.tsv'); tab = OUTDIR/'render/tables'; tab.mkdir(parents=True)
    leaders = {}
    for threshold, name in [('resolved', 'participant_leaders_resolved.tsv'), ('strict_both_methods', 'participant_leaders_strict_both_methods.tsv')]:
        s = src[src.scope.eq('common_functions') & src.threshold.eq(threshold)]
        t = s.groupby('participant')[['agreement_original', 'agreement_simplified']].mean().reset_index()
        assert len(t) == NPART, (threshold, len(t))
        t.to_csv(tab/name, sep='\t', index=False)
        leaders[threshold] = dict(participants=len(t), original=float(t.agreement_original.mean()), simplified=float(t.agreement_simplified.mean()))
    return leaders


def mosaic_numbers(path):
    m = read(path)
    st = m[m.region.eq('stool')][['KO', 'participant', 'species', 'status']].rename(columns={'species': 's', 'status': 'ss'})
    c = m[m.region.ne('stool')].merge(st, on=['KO', 'participant'])
    c = c[c.participant.isin(c[c.status.ne('Not sampled')].participant.unique())]
    both = c[c.status.eq('Resolved') & c.ss.eq('Resolved')]
    ev = []
    for (k, p), x in c.groupby(['KO', 'participant']):
        if x.ss.iloc[0] != 'Resolved': continue
        r = x[x.status.eq('Resolved')]
        if len(r): ev.append(bool((r.species == r.s).all()))
    return dict(capsule_participants=int(c.participant.nunique()), same_leader=int((both.species == both.s).sum()), comparisons=len(c),
                both_resolved=len(both), persistence=int(sum(ev)), persistence_profiles=len(ev))


res = {rule: figure5(rule) for rule in ['simplified', 'original']}
leaders = collect()
s = json.loads((res['simplified']/'SUMMARY.json').read_text())
o = json.loads((res['original']/'SUMMARY.json').read_text())
numbers = dict(fig5f_leaders=leaders, fig5e_mosaic_direct=mosaic_numbers(res['simplified']/'mosaic_donor_regions.tsv'),
               fig5e_mosaic_original=mosaic_numbers(res['original']/'mosaic_donor_regions.tsv'),
               sf12_direct=dict(eligible_functions=s['contributors']['eligible_functions'], eligible_cells=s['contributors']['eligible_cells'],
                                median_resolved_fraction=s['contributors']['median_resolved_fraction'], participants=s['n_participants'],
                                acquisitions=s['n_acquisitions']),
               sf12_original=dict(eligible_functions=o['contributors']['eligible_functions'], eligible_cells=o['contributors']['eligible_cells'],
                                  median_resolved_fraction=o['contributors']['median_resolved_fraction']))
(OUTDIR/'NUMBERS.json').write_text(json.dumps(numbers, indent=2))
print(json.dumps(numbers, indent=2))
print('FIG5EF PASS', flush=True)
