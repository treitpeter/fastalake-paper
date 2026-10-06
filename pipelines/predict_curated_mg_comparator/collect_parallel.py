"""Recount the four distinct supplement contrasts and preserve their input lineages."""
from pathlib import Path
import collections, csv, gzip, hashlib, importlib.util, json, os, re, sys
import numpy as np
import pandas as pd
B=Path(__file__).resolve().parent;R=B.parents[2];E=B.parent
csv.field_size_limit(64*1024*1024)
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8<<20),b''):h.update(b)
    return h.hexdigest()
def stamp(p):p=Path(p).resolve();return dict(path=str(p),bytes=p.stat().st_size,sha256=sha(p))
def load(p):return json.loads(Path(p).read_text())
def dump(p,j):Path(p).write_text(json.dumps(j,indent=2,allow_nan=False)+'\n')
def seq(p):return re.sub(r'\[[^\]]*\]','',p)
def accepted(p,sample):
    result=set();domains=collections.defaultdict(set)
    with Path(p).open() as f:
        for r in csv.DictReader(f,delimiter='\t'):
            assert Path(r['filename']).name==sample+'.mzML'
            if r['rank']=='1' and r['label']=='1' and float(r['peptide_q'])<=.01:
                p=seq(r['peptide']).replace('I','L');assert 5<=len(p)<=50
                result.add(p)
                for a in r['proteins'].split(';'):
                    if a.strip():domains[p].add('Human' if a.strip().startswith('sp|') else 'Microbial')
    assert all(v for v in domains.values()) and set(domains)==result
    classes=collections.Counter('Mixed' if len(v)>1 else next(iter(v)) for v in domains.values())
    assert sum(classes.values())==len(result)
    return result,classes
def best_targets(p,sample):
    out={}
    with Path(p).open() as f:
        for r in csv.DictReader(f,delimiter='\t'):
            assert Path(r['filename']).name==sample+'.mzML'
            if r['label']!='1':continue
            key=r['scannr'];v=(seq(r['peptide']),float(r['hyperscore']),float(r['peptide_q']))
            if key not in out or v[1]>out[key][1]:out[key]=v
    return out
def loss(sample,tasks,sources):
    a,b=[B/'results'/str(tasks[sample,x]['index']) for x in ['loss5_FL','loss5_MG']]
    fa=Path(tasks[sample,'loss5_MG']['fasta']['path'])
    dgfile=R/'canonical/final_revision/thread_selection_math/R3_digest.py'
    spec=importlib.util.spec_from_file_location('archived_digestion',dgfile);dg=importlib.util.module_from_spec(spec);spec.loader.exec_module(dg)
    rawseq=[s for _,s in dg.read_fasta(fa)];digest={p for s in rawseq for p in dg.digest(s)};concat='\n'.join(rawseq)
    aa,bb=[best_targets(p/'results.sage.tsv',sample) for p in [a,b]]
    counts=dict.fromkeys(['A_absent','B_not_cleavable','C_top_fails_q','D_displaced','E_no_target'],0)
    details=[]
    for scan,v in aa.items():
        if v[2]>.01:continue
        w=bb.get(scan)
        if w is not None and w[2]<=.01:continue
        if v[0] not in digest:key='A_absent' if v[0] not in concat else 'B_not_cleavable'
        elif w is None:key='E_no_target'
        elif w[0]==v[0]:key='C_top_fails_q'
        else:key='D_displaced'
        counts[key]+=1;details.append(dict(sample=sample,scan=scan,peptide=v[0],category=key,FL_peptide_q=v[2],comparator_peptide=w[0] if w else '',comparator_peptide_q=w[2] if w else None))
    sources.extend(stamp(p) for p in [a/'results.sage.tsv',b/'results.sage.tsv',fa,dgfile])
    return dict(sample=sample,total=sum(counts.values()),**counts),details

CONTEXT=None
def process_sample(s):
    a16,matched,tasks,bad,source,full=CONTEXT
    sources=[];ladder=[];domains=[];restricted=[];lossrows=[];lossdetail=[];deltas=[];ownsets={}
    t=a16[s]
    if s in matched:
        for arm in ['FL_cohortMG','MG','FL_public']:
            rec=t['sets'][arm];assert sha(rec['path'])==rec['sha256'];values=set(Path(rec['path']).read_text().splitlines())
            psmp=Path(t['configs'][arm]['path']).parent/'results.sage.tsv';actual,cl=accepted(psmp,s)
            assert actual==values,(s,arm,len(actual),len(values))
            ladder.append(dict(sample=s,arm=arm,peptides=len(actual)));sources.extend([stamp(psmp),stamp(rec['path'])])
            if arm in ['FL_cohortMG','MG']:domains.append(dict(sample=s,arm=arm,Human=cl['Human'],Microbial=cl['Microbial'],Mixed=cl['Mixed']))
        if full and s not in bad:
            p=source/'R4K_ownmetag'/(s+'.peptides.txt.gz')
            with gzip.open(p,'rt') as f:own={seq(x.strip()).replace('I','L') for x in f if x.strip()}
            sources.append(stamp(p))
        else:
            p=B/'results'/str(tasks[s,'own_curated']['index'])/'accepted_canonical.txt';own=set(p.read_text().splitlines());sources.append(stamp(p))
        ownsets[s]=own;ladder.append(dict(sample=s,arm='MG_own_curated',peptides=len(own)))
    if full and s not in bad:
        old=source/'R8A_loss_psm5'/(s+'.json');v=load(old);lossrows.append(dict(sample=s,total=v['n_fl_only'],**v['counts']));sources.append(stamp(old))
    else:
        v,dd=loss(s,tasks,sources);lossrows.append(v);lossdetail.extend(dd)
        old=load(source/'R8A_loss_psm5'/(s+'.json'));deltas.extend(dict(sample=s,category=k,old=old['counts'][k],new=v[k],difference=v[k]-old['counts'][k]) for k in old['counts'])
    if (s in bad and full) or not full:
        rp=B/'results'/str(tasks[s,'restricted_lake']['index'])/'results.sage.tsv'
    else:rp=R/'canonical/publication_data/arms/metalab_lake_run/cohorts/PREDICT/5_sage/razor_hash-acc'/s/'results.sage.tsv'
    if s in bad and full:
        mp=E/'mlab_both_stages/results'/s/'step2/results.sage.tsv';assert load(mp.parent.parent/'COMPLETE.json')['status']=='PASS'
    else:mp=R/'canonical/publication_data/arms/metalab_mag/PREDICT'/s/'step2_refined/results.sage.tsv'
    aa,_=accepted(rp,s);bb,_=accepted(mp,s);assert bb
    restricted.append(dict(sample=s,restricted_lake=len(aa),MLAB=len(bb),shared=len(aa&bb),gained=len(aa-bb),lost=len(bb-aa),count_ratio=len(aa)/len(bb),reference_recovery=len(aa&bb)/len(bb)))
    sources.extend([stamp(rp),stamp(mp)]);print('S6 RECOUNT',s,flush=True)
    return sources,ladder,domains,restricted,lossrows,lossdetail,deltas,ownsets

def main(phase):
    code=load(B/'COLLECT_PARALLEL_CODE.json');assert sha(__file__)==code['sha256']
    m=load(B/'manifest.json');assert sha(B/'manifest.json')==(B/'manifest.sha256').read_text().strip()
    full=phase=='full';wanted=m['tasks'] if full else m['tasks'][:4]
    for t in wanted:
        r=load(B/'results'/str(t['index'])/'COMPLETE.json');assert r['status']=='PASS'
        for rec in r['outputs']:assert sha(rec['path'])==rec['sha256']
        actual=load(B/'results'/str(t['index'])/'results.json')
        def subset(expected,observed):
            if isinstance(expected,dict):
                for k,v in expected.items():
                    if k=='output_paths':assert [Path(x).name for x in v]==[Path(x).name for x in observed[k]]
                    elif k not in ['version','output_directory','mzml_paths','fasta']:subset(v,observed[k])
            else:assert expected==observed,(expected,observed)
        subset(t['config'],actual)
    if full:assert load(B/'collection/pilot/COMPLETE.json')['status']=='PASS'
    out=B/'collection_parallel'/phase;out.mkdir(parents=True,exist_ok=False)
    a16p=R/'MM_review_analysis/recovery_2026-09-15/public_comparison_recovery/A16/inputs/manifest.json'
    a16={t['sample']:t for t in load(a16p)['tasks'] if t['cohort']=='PREDICT'}
    inventory=R/'reviews/supplementary_completion_2026-09-13/assignment_exports/PREDICT_296_MLAB/SEARCH_INPUTS.json'
    matched=sorted(t['sample'] for t in load(inventory));assert len(matched)==296
    pilot=m['tasks'][0]['sample'];names=sorted(a16) if full else [pilot]
    if not full:matched=[pilot]
    tasks={(t['sample'],t['arm']):t for t in m['tasks'] if full or t['pilot']}
    bad={t['sample'] for t in m['tasks'] if not t['pilot']};source=R/'canonical/final_revision/thread_selection_math/out'
    sources=[stamp(a16p),stamp(inventory),stamp(B/'manifest.json')];ladder=[];domains=[];restricted=[];lossrows=[];lossdetail=[];deltas=[];ownsets={}
    global CONTEXT
    CONTEXT=(a16,matched,tasks,bad,source,full)
    import concurrent.futures,multiprocessing
    with concurrent.futures.ProcessPoolExecutor(max_workers=4 if full else 2,mp_context=multiprocessing.get_context('fork')) as executor:
        for ss,ll,dd,rr,lr,ld,de,oo in executor.map(process_sample,names):
            sources.extend(ss);ladder.extend(ll);domains.extend(dd);restricted.extend(rr);lossrows.extend(lr);lossdetail.extend(ld);deltas.extend(de);ownsets.update(oo)
    tables={'curation_counts':ladder,'domain_counts':domains,'scoring_counts':lossrows,'scoring_repaired_details':lossdetail,'scoring_historical_differences':deltas,'restricted_resource':restricted}
    for name,rows in tables.items():pd.DataFrame(rows).to_csv(out/(name+'.tsv'),sep='\t',index=False)
    assert len(ladder)==4*len(matched) and len(domains)==2*len(matched) and len(restricted)==len(names)
    assert all(sum(v for k,v in r.items() if k not in ['sample','total'])==r['total'] for r in lossrows)
    assert len(names)==(306 if full else 1)
    totals={k:sum(r[k] for r in lossrows) for k in lossrows[0] if k!='sample'}
    result=dict(status='PASS',phase=phase,matched_acquisitions=len(matched),full_acquisitions=len(names),scoring_totals=totals,
        restricted_median_paired_ratio=float(np.median([r['count_ratio'] for r in restricted])),
        restricted_median_identity_recovery=float(np.median([r['reference_recovery'] for r in restricted])),
        sources=sources,code=stamp(__file__),scientific_review='PENDING',
        limitations=['Scoring comparator is the curated cohort-MG pool, not specimen MG.','Scoring screen retains historical exact-I/L and unmodified-mass digest convention; not full modification-aware searchability.','Restricted source pool remains archived305-FASTA input; no rebuilding from repaired genome selection.','Native search arms, thresholds and selection lineages are distinct; no FDR calibration.'])
    dump(out/'COMPLETE.json',result);print('SUPPLEMENT6 COLLECTION PASS',phase,totals,flush=True)
if __name__=='__main__':main(sys.argv[1])
