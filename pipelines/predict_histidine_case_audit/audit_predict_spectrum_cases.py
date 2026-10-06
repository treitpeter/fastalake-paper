#!/usr/bin/env python3
"""Independent selection, spectrum and source-context audit of the four corrected PREDICT cases."""
import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import re
from pathlib import Path

import numpy as np
from pyteomics import mzml

MM = Path(os.environ.get('MM_REVIEW_ROOT', Path(__file__).resolve().parents[1])).resolve()
SOURCE = MM / 'recovery_2026-09-15/cohort_recovery/A25/results_revision2/PREDICT'
OUT = MM / 'end_of_day_2026-09-15/predict_spectrum_review'
AA = dict(zip('ACDEFGHIKLMNPQRSTVWY', [71.03711,103.00919,115.02694,129.04259,147.0684,
    57.02146,137.05891,113.08406,128.09496,113.08406,131.0405,114.04293,97.05276,
    128.05858,156.1011,87.03203,101.04768,99.06841,186.07932,163.06332]))
PROTON = 1.0072764
csv.field_size_limit(64 * 1024 * 1024)


def load(path): return json.loads(Path(path).read_text())


def rows(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for x in iter(lambda: f.read(8*1024*1024), b''): h.update(x)
    return h.hexdigest()


def verify(r):
    p = Path(r['path']); s = p.stat()
    assert s.st_size == r['bytes']
    if 'mtime_ns' in r: assert s.st_mtime_ns == r['mtime_ns']
    if 'sha256' in r: assert sha(p) == r['sha256'], str(p)


def canonical(s): return re.sub(r'\[[^\]]*\]', '', s).replace('I', 'L')


def masses(peptide):
    parts = re.findall(r'([A-Z])(?:\[([+-]?[\d.]+)\])?', peptide)
    assert ''.join(a + ('['+v+']' if v else '') for a,v in parts) == peptide
    return np.array([AA[a] + float(v or 0) for a,v in parts], dtype=np.float64)


def hits(peptide, total, charge, mz, intensity):
    res = masses(peptide)
    prefix = np.cumsum(res)[:-1]
    found = []
    for offset, b in enumerate(prefix):
        for series, neutral, ordinal in [('b', b, offset+1), ('y', total-b, len(res)-offset-1)]:
            for z in range(1, max(2, min(charge, 3))):
                # Float64 direct mass formula, independent of the producer's float32 helper.
                expected = neutral / z + PROTON
                ppm = (mz - expected) / (expected - PROTON) * 1e6
                candidates = np.flatnonzero(np.abs(ppm) <= 20)
                if len(candidates):
                    ix = int(candidates[np.argmax(intensity[candidates])])
                    found.append(dict(series=series, ordinal=ordinal, charge=z, peak_index=ix,
                                      observed_mz=float(mz[ix]), ppm=float(ppm[ix]), cleavage=offset))
    return found


def fasta(path):
    header, seq = None, []
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                if header is not None: yield header, ''.join(seq)
                header, seq = line[1:].strip().split()[0], []
            else: seq.append(line.strip())
    if header is not None: yield header, ''.join(seq)


def write_table(path, records):
    with path.open('w') as f:
        w = csv.DictWriter(f, fieldnames=list(records[0]), delimiter='\t'); w.writeheader(); w.writerows(records)


def main():
    parser = argparse.ArgumentParser(); parser.add_argument('--phase', choices=['pilot','full'], required=True)
    args = parser.parse_args()
    if args.phase == 'full':
        pilot = load(OUT / 'pilot/COMPLETE.json')
        assert pilot['status'] == 'PASS' and pilot['cases'] == 2
        assert pilot['worker_sha256'] == sha(__file__)
    out = OUT / args.phase; out.mkdir(parents=True, exist_ok=False)
    receipt = load(SOURCE/'COMPLETE.json'); assert receipt['status']=='PASS' and receipt['acquisitions']==306 and receipt['cases']==4
    for r in receipt['outputs'] + receipt['source_receipts']: verify(r)
    manifest = load(MM/'recovery_2026-09-15/cohort_recovery/A25/inputs/manifest_revision2.json')
    tasks = {t['sample']:t for t in manifest['tasks'] if t['cohort']=='PREDICT'}; assert len(tasks)==306
    counts = rows(SOURCE/'complete_candidate_counts.tsv')
    assert len(counts)==4*306 and {r['sample'] for r in counts}==set(tasks)
    all_cases = rows(SOURCE/'selected_cases.tsv'); assert len(all_cases)==4
    strata = rows(SOURCE/'all_strata.tsv'); assert len(strata)==4 and all(s['status']=='SELECTED' for s in strata)
    eligible = rows(SOURCE/'eligible_case_selection_psms.tsv.gz')
    calls = rows(SOURCE/'same_scan_all_search_calls.tsv')
    result, plots = [], []
    selected = [r for r in all_cases if args.phase=='full' or r['case_id'] in ['HIS_002','HIS_004']]
    for row in selected:
        cid, sample = row['case_id'], row['sample']; t=tasks[sample]; base=SOURCE/cid
        prov=load(base/'PROVENANCE.json'); assert prov['selected_psm']['case_id']==cid
        verify(prov['raw']); verify(prov['fasta']); verify(prov['unspecific_config'])
        assert int(row['rank'])==1 and int(row['label'])==1 and float(row['peptide_q'])<=.01 and float(row['spectrum_q'])<=.01
        pep=canonical(row['peptide']); assert 9<=len(pep)<=50
        baseline=set(Path(t['baseline_peptides']['path']).read_text().splitlines()); verify(t['baseline_peptides']); assert pep not in baseline
        choices=sorted((int(c['candidates']),c['sample']) for c in counts if int(c['candidates'])>0 and all(c[k]==row[k] for k in ['cohort','compartment','origin','side']))
        assert choices[len(choices)//2][1]==sample and int(row['acquisition_rank'])==len(choices)//2
        group=[e for e in eligible if e['case_id']==cid]
        group.sort(key=lambda e:(float(e['matched_intensity_pct']),e['canonical_peptide'],e['scannr']))
        assert len(group)==int(row['eligible_peptides']) and len({g['canonical_peptide'] for g in group})==len(group)
        chosen=group[len(group)//2]; assert chosen['canonical_peptide']==pep and chosen['scannr']==row['scannr']
        loci={}; parents=[]
        for header, sequence in fasta(prov['fasta']['path']):
            seq=sequence.replace('I','L'); start=seq.find(pep)
            while start>=0:
                end=start+len(pep); pos=start if row['side']=='N' else end
                loci[(seq,start,end)]=pos; parents.append(header)
                assert 0<pos<len(seq) and seq[pos-1]=='H', (cid,header,start,end)
                start=seq.find(pep,start+1)
        assert len(loci)==1, (cid,len(loci))
        with mzml.MzML(prov['raw']['path'], use_index=True, decode_binary=True) as reader:
            sp=reader.get_by_id(row['scannr'])
        assert sp['id']==row['scannr'] and sp['ms level']==2
        raw_mz=np.asarray(sp['m/z array'],np.float32); raw_i=np.asarray(sp['intensity array'],np.float32)
        observed=rows(base/'observed_peaks.tsv.gz')
        assert np.array_equal(raw_mz,np.array([float(r['mz']) for r in observed],np.float32))
        assert np.array_equal(raw_i,np.array([float(r['intensity']) for r in observed],np.float32))
        keep=np.flatnonzero((raw_mz>=200)&(raw_mz<=2000)); keep=sorted(keep,key=lambda i:(float(raw_i[i]),int(i)))[-150:]; keep=sorted(keep,key=lambda i:float(raw_mz[i]))
        mz=raw_mz[keep].astype(float); intensity=raw_i[keep].astype(float)
        processed=rows(base/'processed_peaks.tsv')
        assert np.array_equal(raw_mz[keep],np.array([float(r['mz']) for r in processed],np.float32))
        assert np.array_equal(raw_i[keep],np.array([float(r['intensity']) for r in processed],np.float32))
        h=hits(row['peptide'],float(row['calcmass']),int(row['charge']),mz,intensity)
        ann=rows(base/'matched_fragments.tsv')
        signatures=lambda rr:{(r['series'],int(r['ordinal']),int(r['charge'])) for r in rr}
        assert signatures(h)==signatures(ann) and len(h)==int(row['matched_peaks']), (cid,len(h),row['matched_peaks'])
        for a in ann:
            x=next(v for v in h if (v['series'],v['ordinal'],v['charge'])==(a['series'],int(a['ordinal']),int(a['charge'])))
            assert abs(x['observed_mz']-float(a['observed_mz']))<.0001
        alternatives=rows(base/'bounded_sequence_alternatives.tsv'); nonoriginal=[]
        for a in alternatives:
            ah=hits(a['peptide'],float(a['calculated_mass']),int(row['charge']),mz,intensity)
            assert len(ah)==int(a['matched_ions']), (cid,a['peptide'])
            if a['original']=='False':nonoriginal.append(a)
        same=[c for c in calls if c['case_id']==cid and c['route']=='semi' and c['label']=='1' and float(c['peptide_q'])<=.01 and float(c['spectrum_q'])<=.01 and canonical(c['peptide'])==pep]
        coverage=len({x['cleavage'] for x in h})/(len(pep)-1)
        result.append(dict(case_id=cid,sample=sample,compartment=row['compartment'],origin=row['origin'],side=row['side'],peptide=row['peptide'],matched_ions=len(h),backbone_coverage=coverage,peptide_q=float(row['peptide_q']),spectrum_q=float(row['spectrum_q']),accepted_in_semi=len(same)>0,nonoriginal_alternatives=len(nonoriginal),nonoriginal_matching_at_least_original=sum(int(a['matched_ions'])>=len(h) for a in nonoriginal),N_flanking_H_inferred=row['side']=='N',parents=';'.join(sorted(set(parents)))))
        plots.append((row,mz,intensity,h)); print('SPECTRUM_CASE_PASS',cid,len(h),len(nonoriginal),bool(same),flush=True)
    write_table(out/'independent_case_review.tsv',result)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(math.ceil(len(plots)/2),2,figsize=(15,3*math.ceil(len(plots)/2)),squeeze=False)
    for ax,(row,mz,intensity,h) in zip(axes.flat,plots):
        y=intensity/max(intensity)*100; ax.vlines(mz,0,y,color='0.72',lw=.6)
        for j,a in enumerate(h):
            i=a['peak_index']; color='#2475A8' if a['series']=='b' else '#B64040'; ax.vlines(mz[i],0,y[i],color=color,lw=1)
            ax.text(mz[i],y[i]+3+(j%2)*5,a['series']+str(a['ordinal'])+('++' if a['charge']==2 else '+'),rotation=90,ha='center',va='bottom',fontsize=7,color=color)
        ax.set_title(row['case_id']+' | '+row['compartment']+' '+row['origin']+' '+row['side']+'\n'+row['peptide'],fontsize=9)
        ax.set_ylim(0,138); ax.set_xlabel('m/z'); ax.set_ylabel('Relative intensity (%)')
    fig.tight_layout();fig.savefig(out/'measured_case_review.png',dpi=170);fig.savefig(out/'measured_case_review.pdf');plt.close(fig)
    report=dict(status='PASS',phase=args.phase,cases=len(result),source_acquisitions=306,all_source_receipt_hashes_verified=True,independent_checks='Deterministic selections, exact indexed RAW peak arrays, preprocessing, float64 fragment/alternative arithmetic, local canonical source loci and same-scan route interpretation',accepted_already_in_semi=sum(r['accepted_in_semi'] for r in result),nonoriginal_alternatives=sum(r['nonoriginal_alternatives'] for r in result),scientific_review='PENDING_VISUAL_AND_WRITTEN_SYNTHESIS',source_receipt_sha256=sha(SOURCE/'COMPLETE.json'),worker_sha256=sha(__file__))
    (out/'COMPLETE.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':main()
