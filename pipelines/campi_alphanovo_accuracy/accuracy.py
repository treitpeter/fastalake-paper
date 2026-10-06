#!/usr/bin/env python3
import os
"""AlphaNovo top-ranked prediction accuracy on CAMPI (15 acquisitions).

Denominator: spectra (scans) identified by the reference-database Sage search
(campi_ref_solo; target, rank 1, peptide_q <= 0.01), reference peptide length in window.
Join: prediction spec_idx + offset == Sage scannr 'scan=N'; offset taken from the
A03 same-spectrum audit (validated by precursor-mass agreement) and re-validated here.
Match: prediction 'peptide_prediction_detokenized_unmodified' reduced to A-Z, I->L,
equal to canonical(reference) (modifications removed, I->L).
Production filter: lib.load_queries (score > 0, 9-50 residues, complexity, top 30% by score).
Negative control: prediction sequences permuted across spectra within acquisition.
"""
import csv, json, random, re, sys, pathlib
sys.path.insert(0, os.path.join(os.environ['PIPELINES_DIR'], 'campi_near_miss_audit'))
from lib import clean, load_queries, COL
from peptide import canonical

A03 = pathlib.Path(os.environ['FASTALAKE_DATA_DIR'], 'MM_review_analysis/A03_4_2_near_misses')
OUT = pathlib.Path(__file__).resolve().parent / 'results'
PROTON = 1.007276466
WINDOWS = {'9-50': (9, 50), '5-50': (5, 50)}


def read_ref(path):
    refs = {}
    with open(path) as f:
        for r in csv.DictReader(f, delimiter='\t'):
            if r['label'] != '1' or int(r['rank']) != 1 or not float(r['peptide_q']) <= .01:
                continue
            scan = int(re.search(r'(?:^| )scan=(\d+)(?:$| )', r['scannr']).group(1))
            assert scan not in refs, 'duplicate accepted scan'
            refs[scan] = (canonical(r['peptide']), float(r['expmass']))
    return refs


def run(sample, pred_path, ref_path, offset, seed=20261005):
    refs = read_ref(ref_path)
    with open(pred_path) as f:
        rows = list(csv.DictReader(f))
    assert {r['beam_rank'] for r in rows} == {'0'}, 'expected top-ranked beam only'
    _, filt, kept = load_queries(pred_path)
    kept_rows = {i for i, _, _ in kept}
    pred = {}   # scan -> (seq_IL, retained, mass)
    for i, r in enumerate(rows, 1):
        scan = int(r['spec_idx']) + offset
        assert scan not in pred, 'duplicate prediction for spectrum'
        m = float(r['prec_mz']) * float(r['prec_charge']) - float(r['prec_charge']) * PROTON
        pred[scan] = (clean(r.get(COL) or '').replace('I', 'L'), i in kept_rows, m)
    # re-validate join by precursor mass (20 ppm)
    paired = [s for s in refs if s in pred]
    ok = sum(abs(pred[s][2] - refs[s][1]) / refs[s][1] * 1e6 <= 20 for s in paired)
    mass_agree = ok / len(paired)
    assert mass_agree >= .95, (sample, mass_agree)
    # negative control: permute predicted sequences across spectra
    scans = sorted(pred)
    seqs = [pred[s][0] for s in scans]
    rng = random.Random(seed + sum(map(ord, sample)))
    perm = seqs[:]; rng.shuffle(perm)
    shuffled = dict(zip(scans, perm))
    out = []
    for w, (lo, hi) in WINDOWS.items():
        den = [s for s, (p, _) in refs.items() if lo <= len(p) <= hi]
        has = [s for s in den if s in pred]
        has_seq = [s for s in has if pred[s][0]]
        exact = [s for s in has if pred[s][0] == refs[s][0]]
        ret = [s for s in has if pred[s][1]]
        ret_exact = [s for s in ret if pred[s][0] == refs[s][0]]
        shuf_exact = [s for s in has if shuffled[s] == refs[s][0]]
        out.append(dict(sample=sample, window=w, n_reference_spectra=len(den),
                        n_with_prediction=len(has), n_with_nonempty_prediction=len(has_seq),
                        n_exact_all=len(exact), n_retained_top30=len(ret), n_exact_retained=len(ret_exact),
                        n_exact_shuffled=len(shuf_exact),
                        frac_with_prediction=len(has) / len(den),
                        frac_exact_all=len(exact) / len(den),
                        frac_exact_among_predicted=len(exact) / len(has),
                        frac_exact_retained_of_reference=len(ret_exact) / len(den),
                        frac_exact_among_retained=len(ret_exact) / len(ret) if ret else float('nan'),
                        frac_exact_shuffled=len(shuf_exact) / len(den),
                        join_offset=offset, join_mass_agree_20ppm=mass_agree,
                        prediction_rows=len(rows), production_retained_rows=filt['retained_rows']))
    return out


def main():
    m = json.loads((A03 / 'inputs/manifest.json').read_text())
    rows = []
    pilot = len(sys.argv) > 1 and sys.argv[1] == 'pilot'
    for s in m['samples'][:1] if pilot else m['samples']:
        name = s['sample']
        offset = json.loads((A03 / 'results' / name / 'COMPLETE.json').read_text())['scan_offset']
        rows += run(name, s['predictions']['path'], s['ref_search']['path'], offset)
        print(name, 'done', flush=True)
    if pilot:
        print(json.dumps(rows, indent=1)); return
    assert len({r['sample'] for r in rows}) == 15
    OUT.mkdir(exist_ok=True)
    with open(OUT / 'per_acquisition.tsv', 'w') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        w.writeheader(); w.writerows(rows)
    import statistics as st
    pooled = []
    for win in WINDOWS:
        rs = [r for r in rows if r['window'] == win]
        S = lambda k: sum(r[k] for r in rs)
        d = dict(window=win, n_acquisitions=len(rs))
        for k in ['n_reference_spectra', 'n_with_prediction', 'n_with_nonempty_prediction', 'n_exact_all',
                  'n_retained_top30', 'n_exact_retained', 'n_exact_shuffled']:
            d[k] = S(k)
        d['pooled_frac_with_prediction'] = d['n_with_prediction'] / d['n_reference_spectra']
        d['pooled_frac_exact_all'] = d['n_exact_all'] / d['n_reference_spectra']
        d['pooled_frac_exact_among_predicted'] = d['n_exact_all'] / d['n_with_prediction']
        d['pooled_frac_exact_retained_of_reference'] = d['n_exact_retained'] / d['n_reference_spectra']
        d['pooled_frac_exact_among_retained'] = d['n_exact_retained'] / d['n_retained_top30']
        d['pooled_frac_exact_shuffled'] = d['n_exact_shuffled'] / d['n_reference_spectra']
        for k in ['frac_with_prediction', 'frac_exact_all', 'frac_exact_retained_of_reference',
                  'frac_exact_among_retained', 'frac_exact_shuffled']:
            v = [r[k] for r in rs]
            d['median_' + k] = st.median(v); d['min_' + k] = min(v); d['max_' + k] = max(v)
        pooled.append(d)
    with open(OUT / 'pooled.tsv', 'w') as f:
        w = csv.DictWriter(f, fieldnames=list(pooled[0]), delimiter='\t', lineterminator='\n')
        w.writeheader(); w.writerows(pooled)
    for d in pooled:
        assert d['pooled_frac_exact_shuffled'] < 0.01, 'negative control failed'
    print(json.dumps(pooled, indent=1))


if __name__ == '__main__':
    main()
