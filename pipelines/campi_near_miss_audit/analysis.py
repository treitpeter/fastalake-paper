#!/usr/bin/env python3
import csv,gzip,hashlib,json,math,os,pathlib,re,shutil,sys,time
from collections import Counter,defaultdict
import ahocorasick
from lib import clean,low,load_queries,sha,meta,verify,dump,table,command,COL
from peptide import canonical
PROTON=1.007276466
AA=set('ACDEFGHIKLMNPQRSTVWY')
CATS=['exact_IL','adjacent_swap','rearrangement','substitution','length_difference','other_disagreement','invalid_prediction','unpaired']
def now():return time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())
def statcheck(m):
    p=pathlib.Path(m['path']);s=p.stat()
    assert (s.st_size,s.st_mtime_ns)==(m['bytes'],m['mtime_ns']),'Changed source '+str(p)
def records(path,hash_result=None):
    h=hashlib.sha256();header=None;parts=[]
    with open(path,'rb') as f:
        for raw in f:
            h.update(raw);line=raw.strip()
            if line.startswith(b'>'):
                if header is not None:yield header,b''.join(parts).decode('ascii').upper().replace('I','L')
                header=line[1:].split()[0].decode('ascii');parts=[]
            elif line:
                assert header is not None,'Sequence before header'
                parts.append(line)
        if header is not None:yield header,b''.join(parts).decode('ascii').upper().replace('I','L')
    if hash_result is not None:hash_result['sha256']=h.hexdigest()
def automaton(seqs):
    a=ahocorasick.Automaton()
    for s in sorted(set(seqs)):
        if s:a.add_word(s,s)
    a.make_automaton();return a
def distance(a,b):
    if len(a)<len(b):a,b=b,a
    row=list(range(len(b)+1))
    for i,x in enumerate(a,1):
        nxt=[i]
        for j,y in enumerate(b,1):nxt.append(min(nxt[-1]+1,row[j]+1,row[j-1]+(x!=y)))
        row=nxt
    return row[-1]
def classify(a,b):
    if not a or any(c not in AA for c in a):return 'invalid_prediction','','',None
    if a==b:return 'exact_IL','','',0
    i=0
    while i<min(len(a),len(b)) and a[i]==b[i]:i+=1
    j=0
    while j<min(len(a),len(b))-i and a[-1-j]==b[-1-j]:j+=1
    x=a[i:len(a)-j];y=b[i:len(b)-j]
    if len(x)==len(y)==2 and x==y[::-1]:cat='adjacent_swap'
    elif len(a)==len(b) and Counter(a)==Counter(b):cat='rearrangement'
    elif len(a)!=len(b):cat='length_difference'
    elif len(x)==len(y)==1:cat='substitution'
    else:cat='other_disagreement'
    return cat,x,y,distance(a,b)
def read_sage(path):
    byscan={};allpeps=set();n=0
    with open(path) as f:
        for r in csv.DictReader(f,delimiter='\t'):
            if r['label']!='1' or int(r['rank'])!=1 or not float(r['peptide_q'])<=.01:continue
            p=canonical(r['peptide']);assert p and all(c in AA for c in p),r['peptide']
            match=re.search(r'(?:^| )scan=(\d+)(?:$| )',r['scannr']);assert match,r['scannr']
            scan=int(match.group(1));assert scan not in byscan,'Ambiguous accepted reference scan'
            r=dict(r,canonical=p);byscan[scan]=r;allpeps.add(p);n+=1
    return byscan,allpeps,n
def families(files):
    values=defaultdict(set)
    for path in files:
        with open(path) as f:
            for r in csv.DictReader(f,delimiter='\t'):
                v=frozenset(x.split('@')[0] for x in r['OG'].split(',') if '@2|Bacteria' in x)
                if v:values[r['accession']].add(v)
    conflicting={k for k,v in values.items() if len(v)>1}
    out={k:set(next(iter(v))) for k,v in values.items() if len(v)==1}
    return out,sorted(conflicting)
def audit_mapping(rows,refs):
    # One mass measurement per spectrum; never multiply mapping evidence by beams.
    masses={}
    for r in rows:
        i=int(r['spec_idx']);v=float(r['prec_mz'])*float(r['prec_charge'])-float(r['prec_charge'])*PROTON
        if i in masses:assert abs(masses[i]-v)<.001,'Inconsistent precursor within spectrum'
        masses[i]=v
    out=[]
    for off in [0,-1,1,7]:
        pairs=[(m,float(refs[i+off]['expmass'])) for i,m in masses.items() if i+off in refs]
        n=len(pairs);ok=sum(abs(a-b)/b*1e6<=20 for a,b in pairs)
        out.append(dict(offset=off,n_paired=n,n_mass_agree=ok,fraction=ok/n if n else 0,role='planted_wrong' if off==7 else 'candidate'))
    best=max(out[:3],key=lambda r:r['fraction'])
    assert best['n_paired']>=100 and best['fraction']>=.95,out
    assert best['fraction']-out[-1]['fraction']>=.8,out
    assert best['fraction']-max(r['fraction'] for r in out[:3] if r is not best)>=.8,out
    return best['offset'],out
def family_status(qstats,rf):
    nf=qstats['candidate_proteins'];known=qstats['annotated_candidates'];fs=set(qstats['candidate_families'])
    any_ok=bool(fs & rf) if rf else None
    exclusive=bool(nf and known==nf and fs and fs<=rf) if rf else None
    return any_ok,exclusive
def fixtures():
    assert classify('ACDEFGHLK','ACDEFGHLK')[0]=='exact_IL'
    assert classify('ACDFEGHLK','ACDEFGHLK')[0]=='adjacent_swap'
    assert classify('ACDEYGHLK','ACDEFGHLK')[0]=='substitution'
    assert classify('ACDEFGHL','ACDEFGHLK')[0]=='length_difference'
    assert distance('ACDFEGHLK','ACDEFGHLK')==2
    a=automaton(['ACDFEGHLK','ACDEFGHLK','RRRR'])
    assert {q for _,q in a.iter('MACDFEGHLKTTACDEFGHLK')}=={'ACDFEGHLK','ACDEFGHLK'}
    assert family_status(dict(candidate_proteins=2,annotated_candidates=1,candidate_families=['F1']),{'F1'})==(True,False)
    assert family_status(dict(candidate_proteins=1,annotated_candidates=1,candidate_families=['F2']),{'F1'})==(False,False)
    assert family_status(dict(candidate_proteins=1,annotated_candidates=1,candidate_families=['F1']),set())==(None,None)
    fake_rows=[dict(spec_idx=str(i),prec_charge='2',prec_mz=str((1000+i*13+2*PROTON)/2)) for i in range(200)]
    fake_ref={i:dict(expmass=str(1000+i*13)) for i in range(200)}
    assert audit_mapping(fake_rows,fake_ref)[0]==0
    try:audit_mapping(fake_rows,{i:dict(expmass='100') for i in range(200)})
    except AssertionError:pass
    else:raise AssertionError('Bad mapping fixture was accepted')
    return dict(status='PASS',tests=['sequence categories','edit distance','co-occurrence','unknown families','correct scan join','planted wrong joins','mass-mismatched join rejection'])
def prepare(run,m):
    for c in m['code']:verify(c)
    prepared=dict(samples=[],common={})
    for key,src in m['common'].items():
        statcheck(src);dst=run/'inputs'/key/pathlib.Path(src['path']).name;dst.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(src['path'],dst);record=meta(dst);assert record['sha256']==sha(src['path'])
        prepared['common'][key]=record
    for sm in m['samples']:
        out=dict(sample=sm['sample'])
        for key,src in sm.items():
            if key=='sample':continue
            statcheck(src)
            if key in ['candidates','hits']:out[key]=src;continue
            dst=run/'inputs'/sm['sample']/key/pathlib.Path(src['path']).name;dst.parent.mkdir(parents=True,exist_ok=True)
            shutil.copy2(src['path'],dst);record=meta(dst);assert record['sha256']==sha(src['path'])
            out[key]=record
        prepared['samples'].append(out)
    prepared['fixtures']=fixtures();prepared['utc']=now();prepared['status']='PASS'
    dump(run/'derived/PREPARED.json',prepared)
def sample(run,m,index):
    prepared=json.loads((run/'derived/PREPARED.json').read_text());assert prepared['status']=='PASS'
    sm=prepared['samples'][index];name=sm['sample'];out=run/'results'/name;out.mkdir(parents=True,exist_ok=True)
    for key,r in sm.items():
        if key!='sample':statcheck(r)
        if key not in ['sample','candidates','hits']:verify(r)
    for r in prepared['common'].values():verify(r)
    predpath=sm['predictions']['path']
    with open(predpath) as f:rows=list(csv.DictReader(f))
    refs,refpeps,nref=read_sage(sm['ref_search']['path'])
    fl,flpeps,nfl=read_sage(sm['fl_search']['path'])
    offset,mapping=audit_mapping(rows,refs);table(out/'scan_mapping.tsv',mapping)
    qs,filter_summary,kept=load_queries(predpath);retained_ids={i for i,_,_ in kept};qset=set(qs)
    # Independent production query-set oracle contains all clean sequences, including rejected ones.
    fixture=out/'all_clean_queries.fasta'
    with fixture.open('w') as f:
        for i,s in enumerate(sorted({clean(r.get(COL) or '') for r in rows})):
            if s:f.write('>q'+str(i)+'\n'+s+'\n')
    oracle=out/'production_query_check';oracle.mkdir(exist_ok=True)
    command([prepared['common']['extractor']['path'],'--peptides-dir',pathlib.Path(predpath).parent,'--database',fixture,'--output',oracle/'matched.fasta','--output-hits',oracle/'hits.tsv','--top-percent','0.30','--min-length','9','--max-length','50','--normalize-il','true','-t','2'],oracle/'command')
    oracle_set=set()
    with (oracle/'hits.tsv').open() as f:
        assert f.readline().strip()=='peptide\tprotein_id\tposition'
        for line in f:oracle_set.add(line.split('\t',1)[0])
    assert oracle_set==qset,'Production query reconstruction mismatch '+name
    ann,conflicts=families([prepared['common'][k]['path'] for k in ['annotation_blind','annotation_ref']])
    ra=automaton(refpeps);rparents=defaultdict(set);rfamilies=defaultdict(set);rannot=Counter()
    reference_hash={}
    for protein,seq in records(prepared['common']['reference_fasta']['path'],reference_hash):
        found={q for _,q in ra.iter(seq)}
        for q in found:
            rparents[q].add(protein)
            if protein in ann:rfamilies[q].update(ann[protein]);rannot[q]+=1
    assert reference_hash['sha256']==prepared['common']['reference_fasta']['sha256']
    assert set(rparents)==refpeps,'Accepted REF peptide absent from REF FASTA'
    pairrefs=defaultdict(set)
    for r in rows:
        ref=refs.get(int(r['spec_idx'])+offset)
        q=clean(r.get(COL) or '').replace('I','L')
        if ref and q:pairrefs[q].add(ref['canonical'])
    pairs={(q,r):dict(candidate_cooccurrence=0,selected_cooccurrence=0,compatible_candidate_examples=[],compatible_selected_examples=[]) for q,rr in pairrefs.items() for r in rr}
    qstats={q:dict(candidate_proteins=0,annotated_candidates=0,candidate_families=set(),selected_proteins=0) for q in qset}
    ref_final=Counter();a=automaton(qset|refpeps)
    def scan(path,is_final=False):
        digest={};n=0
        for protein,seq in records(path,digest):
            n+=1;found={q for _,q in a.iter(seq)};qfound=found&qset
            if is_final:
                for r in found&refpeps:ref_final[r]+=1
            fam=ann.get(protein,set())
            for q in qfound:
                st=qstats[q]
                if is_final:st['selected_proteins']+=1
                else:
                    st['candidate_proteins']+=1
                    if fam:st['annotated_candidates']+=1;st['candidate_families'].update(fam)
                for r in pairrefs.get(q,()):
                    if r in found:
                        ps=pairs[q,r];key='selected_cooccurrence' if is_final else 'candidate_cooccurrence';ps[key]+=1
                        ex='compatible_selected_examples' if is_final else 'compatible_candidate_examples'
                        if len(ps[ex])<3:ps[ex].append(protein)
        digest['entries']=n;return digest
    candidate_scan=scan(sm['candidates']['path'])
    final_scan=scan(sm['final_fasta']['path'],True)
    assert final_scan['sha256']==sm['final_fasta']['sha256']
    # Fresh query-to-protein counts must equal the archived production edge table.
    counts=Counter();last=None;seen=set();h=hashlib.sha256();occ=0
    with open(sm['hits']['path'],'rb') as f:
        line=f.readline();h.update(line);assert line.strip()==b'peptide\tprotein_id\tposition'
        for line in f:
            h.update(line);q,protein,pos=line.decode().rstrip('\n').split('\t');occ+=1
            assert q in qset
            if protein!=last:last=protein;seen=set()
            if q not in seen:counts[q]+=1;seen.add(q)
    assert counts==Counter({q:s['candidate_proteins'] for q,s in qstats.items() if s['candidate_proteins']}),'Archived hit counts differ from FASTA containment '+name
    for q,st in qstats.items():assert st['selected_proteins']<=st['candidate_proteins']
    for p in pairs.values():assert p['selected_cooccurrence']<=p['candidate_cooccurrence']
    statcheck(sm['candidates']);statcheck(sm['hits'])
    qrows=[dict(query=q,length=len(q),retained_prediction_rows=qs[q],**{k:(';'.join(sorted(v)) if isinstance(v,set) else v) for k,v in st.items()}) for q,st in sorted(qstats.items())]
    table(out/'query_lookup.tsv',qrows)
    rr=[dict(reference_peptide=r,reference_parents=len(rparents[r]),annotated_reference_parents=rannot[r],reference_families=';'.join(sorted(rfamilies[r])),final_proteins=ref_final[r],identified_fl_anywhere=int(r in flpeps)) for r in sorted(refpeps)]
    table(out/'reference_retention.tsv',rr)
    pr=[]
    for (q,r),v in sorted(pairs.items()):
        st=qstats.get(q);anyfam,exfam=family_status(st,rfamilies[r]) if st else (None,None)
        pr.append(dict(query=q,reference_peptide=r,query_in_retained_set=int(q in qset),candidate_proteins=st['candidate_proteins'] if st else '',annotated_candidates=st['annotated_candidates'] if st else '',family_any_compatible=anyfam,family_exclusively_compatible=exfam,reference_family_known=bool(rfamilies[r]),**{k:';'.join(vv) if isinstance(vv,list) else vv for k,vv in v.items()}))
    table(out/'query_reference_pairs.tsv',pr)
    fields=['sample','prediction_row','spec_idx','scan','beam_rank','score','prediction_raw','prediction_IL','sequence_with_bags','production_row_retained','query_in_retained_set','reference_accepted','reference_peptide','reference_modified','reference_reported_proteins','reference_peptide_q','reference_spectrum_q','also_spectrum_q_01','category','edit_distance','differing_prediction_block','differing_reference_block','candidate_proteins','annotated_candidates','reference_family_known','family_any_compatible','family_exclusively_compatible','same_candidate_contains_reference','same_selected_contains_reference','reference_present_final','reference_identified_fl_anywhere','reference_identified_fl_same_spectrum']
    sums=defaultdict(Counter);ranks=Counter();bagcount=0
    with gzip.open(out/'spectrum_prediction_fate.tsv.gz','wt',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader()
        for i,row in enumerate(rows,1):
            scanid=int(row['spec_idx'])+offset;ref=refs.get(scanid);raw=row.get(COL) or '';q=clean(raw).replace('I','L');ret=i in retained_ids;ranks[row.get('beam_rank','')]+=1;bagcount+=int('[' in row.get('sequence_with_bags',''))
            d=dict.fromkeys(fields,'');d.update(sample=name,prediction_row=i,spec_idx=row['spec_idx'],scan=scanid,beam_rank=row.get('beam_rank',''),score=row.get('score',''),prediction_raw=raw,prediction_IL=q,sequence_with_bags=row.get('sequence_with_bags',''),production_row_retained=int(ret),query_in_retained_set=int(q in qset),reference_accepted=int(ref is not None),category='unpaired')
            if ref:
                r=ref['canonical'];cat,x,y,ed=classify(q,r);st=qstats.get(q);anyfam,exfam=family_status(st,rfamilies[r]) if st else (None,None);ps=pairs.get((q,r),{})
                d.update(reference_peptide=r,reference_modified=ref['peptide'],reference_reported_proteins=ref['proteins'],reference_peptide_q=ref['peptide_q'],reference_spectrum_q=ref['spectrum_q'],also_spectrum_q_01=int(float(ref['spectrum_q'])<=.01),category=cat,edit_distance=ed if ed is not None else '',differing_prediction_block=x,differing_reference_block=y,reference_family_known=int(bool(rfamilies[r])),reference_present_final=int(ref_final[r]>0),reference_identified_fl_anywhere=int(r in flpeps),reference_identified_fl_same_spectrum=int(scanid in fl and fl[scanid]['canonical']==r))
                if st:d.update(candidate_proteins=st['candidate_proteins'],annotated_candidates=st['annotated_candidates'],family_any_compatible='' if anyfam is None else int(anyfam),family_exclusively_compatible='' if exfam is None else int(exfam),same_candidate_contains_reference=int(ps.get('candidate_cooccurrence',0)>0),same_selected_contains_reference=int(ps.get('selected_cooccurrence',0)>0))
                for pop in ['all_reference_paired']+(['production_retained_reference_paired'] if ret else [])+(['production_retained_reference_paired_both_q01'] if ret and float(ref['spectrum_q'])<=.01 else []):
                    z=sums[pop,cat];z['n']+=1
                    for k in ['reference_present_final','reference_identified_fl_anywhere','reference_identified_fl_same_spectrum','same_candidate_contains_reference','same_selected_contains_reference','family_any_compatible','family_exclusively_compatible']:
                        if d[k]!='':z[k+'_evaluable']+=1;z[k]+=int(d[k])
                    if st:z['any_lookup']+=int(st['candidate_proteins']>0);z['all_candidates_annotated']+=int(st['candidate_proteins']>0 and st['candidate_proteins']==st['annotated_candidates'])
            w.writerow(d)
    summary_rows=[dict(sample=name,population=pop,category=cat,**dict(z)) for (pop,cat),z in sorted(sums.items())]
    keys=sorted(set().union(*(set(x) for x in summary_rows)));table(out/'summary.tsv',summary_rows,keys)
    checks=dict(status='PASS',sample=name,utc=now(),job_id=os.environ.get('SLURM_JOB_ID'),prediction_rows=len(rows),retained_rows=len(kept),reference_psms=nref,fl_psms=nfl,paired_prediction_rows=sum(z['n'] for (p,c),z in sums.items() if p=='all_reference_paired'),recorded_beam_ranks=dict(ranks),rows_with_bags=bagcount,filter=filter_summary,scan_offset=offset,production_query_set='PASS',candidate_counts_equal_archived_hits=True,candidate_fasta=candidate_scan,archived_hits_sha256=h.hexdigest(),archived_hit_occurrences=occ,final_fasta=final_scan,annotation_conflicts=len(conflicts),reference_unique_peptides=len(refpeps),reference_peptides_with_family=sum(bool(rfamilies[r]) for r in refpeps),source_manifest_sha256=sha(run/'inputs/manifest.json'),code=[meta(p) for p in pathlib.Path(__file__).parent.glob('*.py')])
    checks['outputs']=[meta(p) for p in out.iterdir() if p.is_file() and p.name not in ['COMPLETE.json','all_clean_queries.fasta']]
    dump(out/'COMPLETE.json',checks)
def main():
    run=pathlib.Path(sys.argv[1]).resolve();m=json.loads((run/'inputs/manifest.json').read_text())
    for c in m['code']:verify(c)
    stage=sys.argv[2]
    if stage=='prepare':prepare(run,m)
    elif stage=='sample':sample(run,m,int(sys.argv[3]))
    elif stage=='report':
        from report import build
        build(run,m)
    else:raise ValueError(stage)
    marker=stage+('_'+str(sys.argv[3]) if stage=='sample' else '')
    dump(run/'status'/(marker+'.json'),dict(status='PASS',stage=stage,job_id=os.environ.get('SLURM_JOB_ID'),utc=now()))
    print('VALIDATED',marker,flush=True)
if __name__=='__main__':main()
