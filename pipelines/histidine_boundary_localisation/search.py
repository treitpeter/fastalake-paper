#!/usr/bin/env python3
"""Seven-point histidine follow-up: frozen per-acquisition coordinate evidence."""
import csv,gzip,hashlib,json,math,os,re,sys
from collections import Counter,defaultdict
from pathlib import Path
import ahocorasick
import numpy as np
import pandas as pd
RUN=Path(sys.argv[1]).resolve();sys.path.insert(0,str(RUN/'frozen'))
from lib import COL,clean,low,load_queries,meta,verify,dump
from peptide import canonical,bare
from context_helpers import fasta,context
from origin import classify_members,is_classifiable_accession
from histidine_helpers import AA,protein_key,cuts,compatible,background,windows,brute_validate

m=json.loads((RUN/'inputs/manifest.json').read_text());idx=int(sys.argv[2]);nt=len(m['tasks'])
phase='available' if idx<nt else 'unspecific';t=m['tasks'][idx%nt]
out=RUN/'results'/phase/t['cohort']/t['sample'];out.mkdir(parents=True,exist_ok=True)
assert not (out/'COMPLETE.json').exists(),'Refuse completed-output overwrite'
for x in m['code']+m['inputs']+[t['fasta'],t['prediction']]:verify(x)
sources=[];routes={};records=[];route_enzyme={};qc=Counter();reported=defaultdict(set)
fields=['route','peptide','bare','scan','score','charge','precursor_mz','peptide_q','spectrum_q',
        'matched_peaks','matched_intensity_pct','hyperscore','rt','ms2_intensity']
def saveframe(name,rows,columns=None):pd.DataFrame(rows,columns=columns).to_csv(out/name,sep='\t',index=False)
def add(route,p,row):
    routes.setdefault(route,defaultdict(list))[p].append(row)

prediction=[];eligible=[];masses={};beam=Counter()
with open(t['prediction']['path']) as f:
    for number,r in enumerate(csv.DictReader(f),1):
        beam[r.get('beam_rank','MISSING')]+=1;s=clean(r.get(COL,''))
        try:score=float(r.get('score','nan'))
        except ValueError:qc['invalid_score']+=1;continue
        if not math.isfinite(score) or score<=0:continue
        if not (9<=len(s)<=50):continue
        assert set(s)<=set(AA),s
        spec=int(float(r['spec_idx']));charge=float(r['prec_charge']);mz=float(r['prec_mz'])
        row=dict(route='prediction',peptide=s.replace('I','L'),bare=s,scan='prediction:'+str(spec),
                 score=score,charge=charge,precursor_mz=mz,row_number=number)
        prediction.append(row);masses[spec]=charge*(mz-1.007276466621)
        if not low(s):eligible.append((number,s,score))
assert beam=={'0':sum(beam.values())},beam
threshold=sorted((x[2] for x in eligible),reverse=True)[max(1,math.ceil(.3*len(eligible)))-1] if eligible else None
if eligible:
    expected,_,kept=load_queries(t['prediction']['path'])
    assert expected==Counter(s.replace('I','L') for _,s,v in eligible if v>=threshold)
if phase=='available':
    for r in prediction:
        s=r['bare'];p=r['peptide'];v=r['score']
        label=next(label for upper,label in [(0.5,'score_0_0.5'),(.7,'score_0.5_0.7'),(.9,'score_0.7_0.9'),(.95,'score_0.9_0.95'),(.99,'score_0.95_0.99'),(float('inf'),'score_0.99_plus')] if v<upper)
        selected=[label]
        if v>=.9:selected.append('denovo_score09')
        if not low(s):
            if v>=.9:selected.append('denovo_score09_eligible')
            if v>=threshold:selected.append('denovo_top30')
        for route in selected:add(route,p,r);route_enzyme[route]='unspecific'
        records.append(r)
searches=dict(t['routes']) if phase=='available' else {}
if phase=='unspecific':
    receipt=Path(t['unspecific_receipt']);c=json.loads(receipt.read_text());assert c['status']=='PASS'
    assert c['input_fasta']['sha256']==t['fasta']['sha256'];sources.append(meta(receipt))
    oo={Path(x['path']).name:x for x in c['outputs']}
    searches['unspecific']=dict(search=oo['results.sage.tsv'],config=oo['results.json'])
scans={};search_records=[]
for enzyme,source in searches.items():
    verify(source['search']);verify(source['config']);sources.extend(source.values())
    config=json.loads(Path(source['config']['path']).read_text())
    assert Path(config['database']['fasta']).resolve()==Path(t['fasta']['path']).resolve()
    assert [Path(p).resolve() for p in config['mzml_paths']]==[Path(t['raw']['path']).resolve()]
    assert config['database']['static_mods']=={'C':57.0215}
    assert config['database']['variable_mods']=={'M':[15.9949]} and config['database']['max_variable_mods']==2
    assert config['database']['peptide_min_mass']==500 and config['database']['peptide_max_mass']==5000
    assert config['database']['enzyme']['min_len']==5 and config['database']['enzyme']['max_len']==50
    assert config['database']['enzyme']['missed_cleavages']==2
    assert (config['database']['enzyme']['cleave_at']=='')==(enzyme=='unspecific')
    if enzyme!='unspecific':
        assert config['database']['enzyme']['cleave_at']=='KR' and config['database']['enzyme']['restrict']=='P'
        assert bool(config['database']['enzyme']['semi_enzymatic'])==(enzyme=='semi')
    peps=set()
    with open(source['search']['path']) as f:
        for r in csv.DictReader(f,delimiter='\t'):
            if r['rank']!='1':continue
            scan=int(re.search(r'(?:^| )scan=(\d+)',r['scannr'])[1])
            mass=float(r['expmass']);charge=int(r['charge'])
            if scan in scans:assert abs(scans[scan]['mass']-mass)<.01 and scans[scan]['charge']==charge
            scans[scan]=dict(mass=mass,charge=charge)
            if r['label']!='1' or float(r['peptide_q'])>.01:continue
            p=canonical(r['peptide']);s=bare(r['peptide']);peps.add(p)
            if not 9<=len(p)<=50:continue
            row=dict(route='sage_'+enzyme,peptide=p,bare=s,scan='scan:'+str(scan),
                score=np.nan,charge=charge,precursor_mz=mass/charge+1.007276466621,
                **{k:float(r[k]) for k in ['peptide_q','spectrum_q','matched_peaks','matched_intensity_pct','hyperscore','rt','ms2_intensity']})
            search_records.append(row)
            for route in ['sage_'+enzyme]+(['sage_'+enzyme+'_bothq'] if float(r['spectrum_q'])<=.01 else []):
                add(route,p,row);route_enzyme[route]=enzyme;reported[(route,p)].update(r['proteins'].split(';'))
    if enzyme=='tryptic':
        verify(t['baseline_peptides']);assert peps==set(Path(t['baseline_peptides']['path']).read_text().splitlines())
    if enzyme=='unspecific':
        verify(oo['peptide_q01.txt']);assert peps==set(Path(oo['peptide_q01.txt']['path']).read_text().splitlines())
mapping=[]
for offset in [0,-1,1,7]:
    pairs=[(mass,scans[sp+offset]['mass']) for sp,mass in masses.items() if sp+offset in scans]
    good=sum(abs(a-b)/b*1e6<=20 for a,b in pairs)
    mapping.append(dict(offset=offset,pairs=len(pairs),within20ppm=good,fraction=good/len(pairs) if pairs else 0))
best=max(mapping[:3],key=lambda x:x['fraction'])
joined=best['pairs']>=30 and best['fraction']>=.95 and best['fraction']-max(x['fraction'] for x in mapping if x is not best)>=.8
dump(out/'SCAN_MAPPING.json',dict(status='PASS' if joined else 'UNVERIFIED',candidates=mapping,offset=best['offset'] if joined else None))
if joined:
    for row in prediction:row['scan']='scan:'+str(int(row['scan'].split(':')[1])+best['offset'])
records.extend(search_records)
saveframe('evidence_rows.tsv.gz',records,fields+['row_number'])
if joined:
    dn=defaultdict(list)
    for r in prediction:dn[r['scan']].append(r)
    agreement=[]
    for r in search_records:
        for d in dn.get(r['scan'],[]):
            agreement.append(dict(scan=r['scan'],sage_route=r['route'],prediction=d['peptide'],sage=r['peptide'],
                score=d['score'],top30=bool(not low(d['bare']) and d['score']>=threshold),score09=d['score']>=.9,
                exact_IL=d['peptide']==r['peptide'],prediction_C=d['bare'][-1],sage_C=r['bare'][-1],
                spectrum_q=r['spectrum_q']))
    saveframe('same_spectrum.tsv.gz',agreement,['scan','sage_route','prediction','sage','score','top30','score09','exact_IL','prediction_C','sage_C','spectrum_q'])
universe=set().union(*(set(p) for p in routes.values())) if routes else set()
machine=ahocorasick.Automaton()
for p in sorted(universe):machine.add_word(p,p)
if universe:machine.make_automaton()
occ=defaultdict(list);proteins={};bykey=defaultdict(set);audit={};protein_cuts={}
for h,seq in fasta(t['fasta']['path'],audit):
    pid=h.split()[0];assert pid not in proteins;proteins[pid]=seq;key=protein_key(seq);bykey[key].add(pid)
    protein_cuts[pid]=cuts(seq)
    if universe:
        for end,p in machine.iter(seq.replace('I','L')):occ[p].append((pid,end+1-len(p),end+1))
assert audit['sha256']==t['fasta']['sha256']
if t.get('pilot'):
    # Independent brute-force enumeration on real source sequences, including C/M.
    chosen=sorted(proteins)[:2]+[p for p in sorted(proteins) if 'C' in proteins[p] and 'M' in proteins[p]][:1]
    for pid in chosen:
        for enzyme in ['tryptic','semi','unspecific']:
            for side in ['N','C']:brute_validate(proteins[pid][:100],enzyme,side)
bgcache={};site_rows=[];summary=[];peptide_rows=[];enrichment=[];motif=[]
for route,peps in routes.items():
    enzyme=route_enzyme[route];sites={};mapped=0;unambiguous=0
    for p,rr in peps.items():
        compatible_occ=[]
        for pid,start,end in occ[p]:
            c=protein_cuts[pid]
            if enzyme!='unspecific' and not (int(c[start+1:end].sum())<=2 and ((c[start] and c[end]) if enzyme=='tryptic' else (c[start] or c[end]))):continue
            compatible_occ.append((pid,start,end))
        members={x[0] for x in compatible_occ};locs={(protein_key(proteins[pid]),s,e) for pid,s,e in compatible_occ}
        origin=('Unknown' if any(not is_classifiable_accession(pid) for pid in members) else classify_members(members)) if members else 'NotInSelectedFASTA'
        if route.startswith('sage_'):assert members and reported[(route,p)]<=members,(route,p)
        mapped+=bool(members);unique=len(locs)==1;unambiguous+=unique
        peptide_rows.append(dict(route=route,peptide=p,length=len(p),origin=origin,parents=';'.join(sorted(members)),
            compatible_loci=len(locs),unambiguous_locus=unique,rows=len(rr),C_H=p.endswith('H'),C_nonKR=p[-1] not in 'KR'))
        for pid,start,end in compatible_occ:
            seq=proteins[pid];key=protein_key(seq)
            for side,pos in [('N',start),('C',end)]:
                if not 0<pos<len(seq):continue
                sk=(key,pos)
                if sk not in sites:sites[sk]=dict(N=set(),C=set(),N_scans=set(),C_scans=set(),N_unique=set(),C_unique=set(),
                    N_unique_scans=set(),C_unique_scans=set(),parents=set(),windows=set(),origin=set())
                v=sites[sk];v[side].add(p);v[side+'_scans'].update(r['scan'] for r in rr)
                if unique:v[side+'_unique'].add(p);v[side+'_unique_scans'].update(r['scan'] for r in rr)
                v['parents'].add(pid);v['windows'].add(windows(seq,[pos])[0]);v['origin'].add(origin)
    summary.append(dict(route=route,peptides=len(peps),rows=sum(len(x) for x in peps.values()),mapped=mapped,unambiguous=unambiguous))
    aggregates={}
    for (key,pos),v in sites.items():
        origin=next(iter(v['origin'])) if len(v['origin'])==1 else ('Unknown' if 'Unknown' in v['origin'] else 'Mixed')
        for policy in ['all_compatible','unambiguous']:
            nn=v['N_unique'] if policy=='unambiguous' else v['N'];cc=v['C_unique'] if policy=='unambiguous' else v['C']
            ns=v['N_unique_scans'] if policy=='unambiguous' else v['N_scans'];cs=v['C_unique_scans'] if policy=='unambiguous' else v['C_scans']
            independent=bool(ns and cs and len(ns|cs)>1)
            site_rows.append(dict(route=route,policy=policy,protein_key=key,position=pos,origin=origin,
                parents=';'.join(sorted(v['parents'])),motifs=';'.join(sorted(v['windows'])),
                N_peptides=';'.join(sorted(nn)),C_peptides=';'.join(sorted(cc)),N_scans=';'.join(sorted(ns)),C_scans=';'.join(sorted(cs)),
                N_count=len(nn),C_count=len(cc),both_independent_spectra=independent))
            if policy!='unambiguous':continue
            for side,pp in [('N',nn),('C',cc)]:
                if not pp:continue
                # Equal weight over native I/L sequence variants of one canonical protein.
                native_seqs=sorted({proteins[pid] for pid in v['parents']})
                for seq in native_seqs:
                    residue=seq[pos-1]
                    if residue not in AA or residue in 'KR':continue
                    weight=1/len(native_seqs);w=windows(seq,[pos])[0]
                    cachekey=(seq,enzyme,side)
                    if cachekey not in bgcache:bgcache[cachekey]=background(seq,enzyme,side)
                    bg=bgcache[cachekey]
                    for origin_class in ['All',origin]:
                        for population,bgname in [('all_nonKR','eligible'),('H_only','H'),('all_nonKR_internal_background','internal')]:
                            if population=='H_only' and residue!='H':continue
                            k=(side,origin_class,population)
                            if k not in aggregates:aggregates[k]=[np.zeros((8,20)),np.zeros((8,20)),0,0]
                            obs,exp,n,undefined=aggregates[k];bm=bg[bgname];den=bm.sum(axis=1)
                            aggregates[k][2]+=weight
                            if den[3]==0:aggregates[k][3]+=weight;continue
                            for j,a in enumerate(w):
                                if a in AA and den[j]>0:
                                    obs[j,AA.index(a)]+=weight;exp[j]+=weight*bm[j]/den[j]
    for (side,origin,pop), (obs,exp,n,undefined) in aggregates.items():
        for j in range(8):
            for ai,a in enumerate(AA):
                motif.append(dict(route=route,side=side,origin=origin,population=pop,position=j-4 if j<4 else j-3,
                    residue=a,observed=obs[j,ai],expected=exp[j,ai],observed_total=obs[j].sum(),expected_total=exp[j].sum(),
                    enrichment=obs[j,ai]/exp[j,ai] if exp[j,ai] else np.nan,sites=n,undefined_background=undefined))
        for ai,a in enumerate(AA):
            enrichment.append(dict(route=route,side=side,origin=origin,population=pop,residue=a,
                observed=obs[3,ai],expected=exp[3,ai],denominator=obs[3].sum(),sites=n,undefined_background=undefined,
                enrichment=obs[3,ai]/exp[3,ai] if exp[3,ai] else np.nan))
    print('HISTIDINE ROUTE',phase,t['sample'],route,len(peps),len(sites),flush=True)
saveframe('peptide_mapping.tsv.gz',peptide_rows)
saveframe('site_evidence.tsv.gz',site_rows)
saveframe('route_summary.tsv',summary,['route','peptides','rows','mapped','unambiguous'])
saveframe('source_enrichment.tsv',enrichment)
saveframe('motif_counts.tsv.gz',motif)
dump(out/'COMPLETE.json',dict(status='PASS',phase=phase,cohort=t['cohort'],sample=t['sample'],design=t['design'],
    pilot=bool(t.get('pilot')),production_threshold=threshold,beam_ranks=dict(beam),scan_mapping=joined,
    prediction_residue_confidence='NOT_EXPORTED',scientific_review='PENDING',sources=sources,
    code=m['code'],source_fasta=t['fasta'],prediction=t['prediction'],job_id=os.environ.get('SLURM_JOB_ID'),
    outputs=[meta(p) for p in out.iterdir() if p.is_file() and p.name!='COMPLETE.json']))
print('HISTIDINE ACQUISITION PASS',phase,t['sample'],flush=True)
