"""CR-016: PREDICT antibiotic endpoints with rifaximin for hepatic encephalopathy (HE)
handled three ways. Re-derives antibiotic_state from the metadata flags (same rule as
FIGURATOR/Figure4/v5_2026-09-07/derive_clinical.py l.21-27), asserts it equals the
state used by C/clinical_refit, then reruns the paired inference of
C/clinical_refit.py::derive (same bootstrap/sign-flip sizes, seeds and BH family).

Scenarios
  A_current          rule as published (visit_rifaximin_for_he ignored)
  B_HE_as_exposure   visit_rifaximin_for_he added to the core flags (1 -> On; Off needs it = 0)
  C_HE_patients_excl patients with any visit_rifaximin_for_he = 1 acquisition removed
  D_HE_off_rows_excl only acquisitions labelled Off with visit_rifaximin_for_he = 1 set Unknown

Outputs contain counts and statistics only, no patient or sample identifiers.
Usage: rifaximin_sensitivity.py OUTDIR [--negative-control]
"""
from pathlib import Path
import sys, json, hashlib, os
from datetime import datetime, timezone
import numpy as np, pandas as pd
from scipy.stats import rankdata

ROOT = Path(__import__('os').environ['FASTALAKE_DATA_DIR'])
C = ROOT/'review/v15_integration/integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01'
META = Path(os.environ['FASTALAKE_CONTROLLED_DIR'])/'PREDICT_metadata.tsv'
AC = Path(os.environ['FASTALAKE_CONTROLLED_DIR'])/'acquisition_accounting.tsv'
STATS = C/'clinical_refit/clinical_paired_statistics.tsv'
EDS_ORDER = C/'biology_refit/eds_refit/input_data/heatmap_85_patient_order.tsv'
EXPECT = {'meta': 'fc5dbf07df2fd25c7def3e9d583b2745c9f237b12a6a9c6e78ecefb249e90238'}
CLASSES = ['Human', 'Microbial_species', 'Microbial_above_species', 'Other_annotated_origins', 'Unresolved_origin']

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def bh(p):
    p = np.asarray(p, float); o = np.argsort(p); r = p[o]*len(p)/np.arange(1, len(p)+1)
    q = np.empty_like(p); q[o] = np.minimum(1, np.minimum.accumulate(r[::-1])[::-1]); return q

def derive_state(d, extra_core=()):
    core = ['visit_antibiotic_prophylaxis_for_sbp', 'visit_antibiotic_prophylaxis_for_variceal_bleeding',
            'visit_antibiotic_prophylaxis_for_intervention', 'visit_other_type_of_antibiotic',
            'visit_antibiotic_for_infection_nonprophylactic', 'visit_rifaximin', 'visit_norfloxacin',
            'visit_cephalosporines', 'visit_penicillines'] + list(extra_core)
    flags = core + [c for c in d if c.startswith('visit_antibiotic_treatment') and c != 'visit_antibiotic_treatment_other_specify']
    assert set(d[flags].stack()) <= {'', '0', '1'}
    pos = d[flags].eq('1').any(axis=1)
    return pd.Series(np.where(pos, 'On', np.where(d[core].eq('0').all(axis=1), 'Off', 'Unknown')), index=d.index)

def paired_inference(ac, state):
    """Transcription of clinical_refit.py::derive paired block, without the n=42 assertion."""
    a = ac.assign(antibiotic_state=state.to_numpy())
    st = a[a.antibiotic_state.ne('Unknown')].groupby(['patient_id', 'antibiotic_state'])
    paired = st[CLASSES].mean().join(st.pg_parsimony.median()).join(st.size().rename('n_acquisitions')).reset_index()
    el = paired.groupby('patient_id').antibiotic_state.nunique()
    paired = paired[paired.patient_id.isin(el[el.eq(2)].index)]
    n = paired.patient_id.nunique(); assert len(paired) == 2*n
    rows = []
    for idx, (metric, scale) in enumerate([('Human', 100), ('pg_parsimony', 1)]):
        wide = paired.pivot(index='patient_id', columns='antibiotic_state', values=metric)[['Off', 'On']]*scale
        delta = (wide.On-wide.Off).to_numpy()
        rng = np.random.default_rng(20260907)
        boots = delta[rng.integers(0, len(delta), (50000, len(delta)))].mean(axis=1)
        lo, hi = np.quantile(boots, [.025, .975])
        rng = np.random.default_rng(20260907+idx); count = 0
        for _ in range(100):
            null = ((rng.integers(0, 2, (2000, len(delta)))*2-1)*delta).mean(axis=1)
            count += int((abs(null) >= abs(delta.mean())-1e-12).sum())
        rows.append(dict(metric=metric, n_patients=n, mean_off=wide.Off.mean(), mean_on=wide.On.mean(),
                         mean_on_minus_off=delta.mean(), median_on_minus_off=float(np.median(delta)),
                         n_patients_increase=int((delta > 0).sum()), CI95_low=lo, CI95_high=hi,
                         p_signflip_MC=(count+1)/200001))
    s = pd.DataFrame(rows); s['q_two_clinical_endpoints'] = bh(s.p_signflip_MC)
    comp = paired.groupby('antibiotic_state')[CLASSES].mean().reset_index()  # Fig. 4f analogue, mean over patients
    acq = paired.groupby('antibiotic_state').n_acquisitions.sum().to_dict()
    return s, comp, paired, acq

def resid_spearman(x, y, Z):
    Z = np.column_stack([np.ones(len(x))]+Z); r = np.column_stack([rankdata(x), rankdata(y)])
    e = r-Z @ np.linalg.lstsq(Z, r, rcond=None)[0]; return float(np.corrcoef(e.T)[0, 1])

def main(out, negative_control=False):
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    assert sha(META) == EXPECT['meta']
    d = pd.read_csv(META, sep='\t', dtype=str, keep_default_na=False).set_index('sample_id')
    ac = pd.read_csv(AC, sep='\t', dtype={'patient_id': str})
    assert len(ac) == ac.sample_id.nunique() == 306
    d = d.loc[ac.sample_id]
    he = d.visit_rifaximin_for_he
    if negative_control:
        # Planted control: flip every HE flag. If the pipeline is wired correctly the
        # B scenario must differ from the HE-as-exposure result and the A gate must still pass.
        he = he.map({'0': '1', '1': '0', '': ''})
        d = d.assign(visit_rifaximin_for_he=he)
    sA = derive_state(d)
    gate = (sA.to_numpy() == ac.antibiotic_state.to_numpy())
    assert gate.all(), f'{(~gate).sum()} state mismatches against clinical_refit'
    sB = derive_state(d, extra_core=['visit_rifaximin_for_he'])
    he1 = he.eq('1').to_numpy()
    pid = ac.patient_id.to_numpy()
    he_pat = set(pid[he1])
    sC = pd.Series(np.where(np.isin(pid, list(he_pat)), 'Unknown', sA), index=d.index)
    sD = pd.Series(np.where(he1 & (sA.to_numpy() == 'Off'), 'Unknown', sA), index=d.index)
    scen = {'A_current': sA, 'B_HE_as_exposure': sB, 'C_HE_patients_excluded': sC, 'D_HE_off_rows_excluded': sD}

    counts = dict(
        acquisitions=306,
        HE_flag_values={k: int(v) for k, v in he.value_counts().items()},
        HE_by_current_state={f'{h}|{s}': int(v) for (h, s), v in pd.crosstab(he, sA).stack().items()},
        HE_rows_patients=len(he_pat),
        HE_off_rows=int((he1 & (sA.to_numpy() == 'Off')).sum()),
        HE_off_row_patients=int(len(set(pid[he1 & (sA.to_numpy() == 'Off')]))),
        state_changes_A_to_B={f'{a}->{b}': int(v) for (a, b), v in pd.crosstab(sA, sB).stack().items() if v and a != b},
        HE_with_core_rifaximin=int((he1 & d.visit_rifaximin.eq('1').to_numpy()).sum()),
        HE_with_treatment_rifaximin=int((he1 & d.visit_antibiotic_treatment_rifaximin.eq('1').to_numpy()).sum()),
    )
    allstats, allcomp, membership = [], [], {}
    for name, s in scen.items():
        st, comp, paired, acq = paired_inference(ac, s)
        st.insert(0, 'scenario', name); comp.insert(0, 'scenario', name)
        st['acq_off'] = acq.get('Off'); st['acq_on'] = acq.get('On')
        allstats.append(st); allcomp.append(comp); membership[name] = set(paired.patient_id)
    A = membership['A_current']
    for name, m in membership.items():
        counts[f'paired_{name}'] = dict(n=len(m), dropped_vs_A=len(A-m), added_vs_A=len(m-A),
                                        HE_patients_in_pairs=len(m & he_pat))
    counts['HE_patients_in_A_pairs'] = len(A & he_pat)
    counts['A_pairs_with_HE_off_rows'] = len(A & set(pid[he1 & (sA.to_numpy() == 'Off')]))
    stats = pd.concat(allstats, ignore_index=True); comp = pd.concat(allcomp, ignore_index=True)

    # Baseline gate against the frozen refit statistics (scenario A must reproduce exactly).
    ref = pd.read_csv(STATS, sep='\t')
    a = stats[stats.scenario.eq('A_current')].reset_index(drop=True)
    for c in ['n_patients', 'mean_off', 'mean_on', 'mean_on_minus_off', 'CI95_low', 'CI95_high', 'p_signflip_MC', 'q_two_clinical_endpoints']:
        assert np.allclose(a[c], ref[c], rtol=1e-10, atol=1e-12), c
    counts['gate_A_reproduces_clinical_refit'] = True

    # EDS 85-patient covariate (descriptive): how many of the 85 displayed samples change state,
    # and the antibiotic-only adjusted Spearman rho for the three approved pairs under A and B.
    o = pd.read_csv(EDS_ORDER, sep='\t', dtype={'patient_id': str})
    o_he = d.loc[o.sample_id, 'visit_rifaximin_for_he'].eq('1').to_numpy() if set(o.sample_id) <= set(d.index) else None
    eds = dict(n=len(o))
    if o_he is not None:
        oA = sA.loc[o.sample_id].to_numpy(); oB = sB.loc[o.sample_id].to_numpy()
        assert (oA == o.antibiotic_state.to_numpy()).all()
        eds.update(HE_flag=int(o_he.sum()), changes_A_to_B=int((oA != oB).sum()),
                   state_counts_A={k: int(v) for k, v in pd.Series(oA).value_counts().items()},
                   state_counts_B={k: int(v) for k, v in pd.Series(oB).value_counts().items()})
    counts['EDS_85'] = eds
    if negative_control:
        counts['NEGATIVE_CONTROL'] = 'HE flags inverted; results are NOT for reporting'
    stats.to_csv(out/'paired_endpoints_by_scenario.tsv', sep='\t', index=False)
    comp.to_csv(out/'paired_composition_by_scenario.tsv', sep='\t', index=False)
    (out/'counts.json').write_text(json.dumps(counts, indent=1)+'\n')
    (out/'COMPLETE.json').write_text(json.dumps(dict(status='PASS', negative_control=negative_control,
        utc=datetime.now(timezone.utc).isoformat(), job_id=os.environ.get('SLURM_JOB_ID'),
        script_sha256=sha(__file__), inputs={str(p): sha(p) for p in [META, AC, STATS, EDS_ORDER]},
        outputs={p.name: sha(p) for p in out.iterdir() if p.is_file() and p.name != 'COMPLETE.json'}), indent=1)+'\n')
    print(stats.to_string()); print(json.dumps(counts, indent=1))

if __name__ == '__main__':
    main(sys.argv[1], '--negative-control' in sys.argv)
