"""C4 (JAG55) rerun with participant 15: identical method, 347-acquisition roster.

Imports the unmodified ../C4_dominance_evenness/dominance_evenness.py (sha recorded) and reuses
load_a33 / per_acquisition / summarise / metrics / signflip unchanged. Only inputs and roster change:
  design   C2_participant15/fig5e_f/design_347.tsv (331 + 16 participant-15 rows)
  species  C2_participant15/fig5e_f/run347/CAPSCAN__simplified/source_data/Figure5/function_species.npz
           (+ acquisition_composition_local_candidate.tsv), the run347 output behind Fig. 5e/f (n=15).
The 'repaired' annotation source has no 347 equivalent and is not computed here.

Positive control P1: the same imported functions on the old A33 inputs/331 roster must reproduce
19.8% / 10.9% and 0.630 / 0.673 and equal the a33 rows of ../C4_dominance_evenness/summary_tests.tsv.
"""
from pathlib import Path
import importlib.util, json
import numpy as np, pandas as pd

OUT = Path(__file__).resolve().parent
C = OUT.parent
ORIG = C / 'C4_dominance_evenness/dominance_evenness.py'
spec = importlib.util.spec_from_file_location('de', ORIG); de = importlib.util.module_from_spec(spec); spec.loader.exec_module(de)
W = C / 'C2_participant15/fig5e_f'
F347 = W / 'run347/CAPSCAN__simplified/source_data/Figure5'
NEW = {'design_347': W / 'design_347.tsv', 'a33_function_species': F347 / 'function_species.npz',
       'a33_acquisitions': F347 / 'acquisition_composition_local_candidate.tsv'}
OLD = {k: de.INPUTS[k] for k in ['design_331', 'a33_function_species', 'a33_acquisitions']}


def run(design, fs, acqf):
    de.INPUTS['a33_function_species'] = fs; de.INPUTS['a33_acquisitions'] = acqf
    taxa, samples, total, nonhost, _ = de.load_a33(design)
    acq = de.per_acquisition('a33', taxa, samples, nonhost, design)
    assert len(acq) == len(design)
    s, p = de.summarise(acq)
    return acq, s, p, len(samples)


def main():
    hashes = {k: de.sha(v) for k, v in {**OLD, **NEW}.items()}
    hashes['orig_script'] = de.sha(ORIG); hashes['this_script'] = de.sha(__file__)
    # P1 positive control
    d331 = pd.read_csv(OLD['design_331'], sep='\t'); assert len(d331) == 331
    _, s0, _, n0 = run(d331, OLD['a33_function_species'], OLD['a33_acquisitions'])
    s0 = s0.set_index('metric')
    assert int(s0.loc['top1', 'paired_participants']) == 14
    assert round(100 * s0.loc['top1', 'participant_median_fluid'], 1) == 19.8 and round(100 * s0.loc['top1', 'participant_median_stool'], 1) == 10.9
    assert round(s0.loc['pielou', 'participant_median_fluid'], 3) == 0.630 and round(s0.loc['pielou', 'participant_median_stool'], 3) == 0.673
    ref = pd.read_csv(C / 'C4_dominance_evenness/summary_tests.tsv', sep='\t'); ref = ref[ref.source == 'a33'].set_index('metric')
    num = [c for c in ref.columns if c not in ('source',)]
    assert np.allclose(s0.loc[ref.index, num].astype(float), ref[num].astype(float), rtol=1e-12, atol=0), 'P1 mismatch'
    p1 = dict(status='PASS', matrix_samples=n0, top1=[float(s0.loc['top1', 'participant_median_fluid']), float(s0.loc['top1', 'participant_median_stool'])],
              pielou=[float(s0.loc['pielou', 'participant_median_fluid']), float(s0.loc['pielou', 'participant_median_stool'])])
    # 347
    d = pd.read_csv(NEW['design_347'], sep='\t')
    assert len(d) == 347 and d.sample_id.is_unique and d.participant.nunique() == 15
    acq, s, p, n = run(d, NEW['a33_function_species'], NEW['a33_acquisitions'])
    assert (s.paired_participants == 15).all(), s.paired_participants
    acq.to_csv(OUT / 'per_acquisition_metrics.tsv', sep='\t', index=False)
    p.to_csv(OUT / 'participant_medians.tsv', sep='\t', index=False)
    s.to_csv(OUT / 'summary_tests.tsv', sep='\t', index=False)
    (OUT / 'RUN.json').write_text(json.dumps(dict(status='PASS', P1_positive_control=p1, matrix_samples_347=n,
        roster_acquisitions=347, undefined=acq[~np.isfinite(acq.top1)].sample_id.tolist(), input_sha256=hashes), indent=2, default=float) + '\n')
    print(s.to_string(index=False)); print(json.dumps(p1, default=float))


if __name__ == '__main__':
    main()
