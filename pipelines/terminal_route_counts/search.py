#!/usr/bin/env python3
"""Audit actual local peptide contexts and all terminal residues across search routes."""
import csv,gzip,hashlib,json,math,os,sys
from collections import Counter,defaultdict
from pathlib import Path
RUN=Path(sys.argv[1]).resolve();sys.path.insert(0,str(RUN/'frozen'))
from lib import COL,clean,load_queries,meta,sha,verify,dump,table
from peptide import canonical,bare
from digest import auto,fixture
from origin import classify_members,is_classifiable_accession
from parameters import parameters
AA='ACDEFGHIKLMNPQRSTVWY'

def fasta(path,audit):
    digest=hashlib.sha256();header=None;parts=[];n=0
    with open(path,'rb') as f:
        for raw in f:
            digest.update(raw);line=raw.strip()
            if line.startswith(b'>'):
                if header is not None:
                    n+=1;yield header,b''.join(parts).decode('ascii').upper()
                header=line[1:].decode('ascii');parts=[]
            elif line:parts.append(line)
    if header is not None:
        n+=1;yield header,b''.join(parts).decode('ascii').upper()
    audit.update(sha256=digest.hexdigest(),entries=n)

def context(s,start,end):
    n=start==0 or (s[start-1] in 'KR' and s[start]!='P')
    c=end==len(s) or (s[end-1] in 'KR' and s[end]!='P')
    missed=sum(s[i] in 'KR' and s[i+1]!='P' for i in range(start,end-1))
    return (s[start-1] if start else 'PROTEIN_END',s[end] if end<len(s) else 'PROTEIN_END',n,c,missed)

def eligible(ctx,route):
    return route in ('unspecific','denovo') or (ctx[4]<=2 and ((ctx[2] and ctx[3]) if route=='tryptic' else (ctx[2] or ctx[3])))

def psmread(path):
    seqs=Counter();parents=defaultdict(set);qc=Counter()
    with open(path) as f:
        rd=csv.DictReader(f,delimiter='\t');assert {'rank','label','peptide_q','spectrum_q','peptide','proteins'}<=set(rd.fieldnames)
        for r in rd:
            qc['reported_rows']+=1
            if r['rank']!='1':continue
            q=float(r['peptide_q']);sq=float(r['spectrum_q']);assert math.isfinite(q) and 0<=q<=1 and math.isfinite(sq) and 0<=sq<=1
            assert r['label'] in ('1','-1')
            if q>.01:continue
            qc['decoy_psms_q01' if r['label']=='-1' else 'target_psms_q01']+=1
            if r['label']=='-1':continue
            s=bare(r['peptide']);p=canonical(r['peptide']);assert 5<=len(p)<=50 and s and p==s.replace('I','L')
            seqs[s]+=1;parents[p].update(x for x in r['proteins'].split(';') if x and not x.startswith('rev_'))
    return seqs,parents,dict(qc)

def checked_config(x,t):
    verify(x);c=json.loads(Path(x['path']).read_text())
    assert len(c['mzml_paths'])==1 and Path(c['mzml_paths'][0]).resolve()==Path(t['raw']['path']).resolve()
    assert Path(c['database']['fasta']).resolve()==Path(t['fasta']['path']).resolve()
    assert c.get('version','0.14.6')=='0.14.6'
    return c

def run(index):
    fixture();assert not context('MACDEK',1,6)[2];assert context('KACDEK',1,6)[2]
    m=json.loads((RUN/'inputs/manifest.json').read_text())
    for x in m['frozen_code']+[m['protocol']]:verify(x)
    t=m['tasks'][index];out=RUN/'results'/t['cohort']/t['sample'];out.mkdir(parents=True,exist_ok=True)
    assert not (out/'COMPLETE.json').exists(),'Refuse completed-output overwrite'
    for x in t['inputs']:verify(x)
    receipt=json.loads(Path(t['unspecific_receipt']).read_text());assert receipt['status']=='PASS' and receipt['sample']==t['sample']
    assert receipt['input_fasta']['sha256']==sha(t['fasta']['path'])
    source_receipt=meta(t['unspecific_receipt']);oo={Path(x['path']).name:x for x in receipt['outputs']}
    unspecific_config=oo['results.json'];checked_config(unspecific_config,t)
    route_inputs=dict(t['routes'],unspecific=dict(search=oo['results.sage.tsv'],config=unspecific_config))
    configs={r:checked_config(v['config'],t) for r,v in route_inputs.items()}
    search_parameters=[]
    for route,c in configs.items():
        p=parameters(c);e=p['database'].pop('enzyme');search_parameters.append(p)
        assert e['missed_cleavages']==2 and e['min_len']==5 and e['max_len']==50
        assert (e['cleave_at']=='')==(route=='unspecific')
        if route!='unspecific':assert e['cleave_at']=='KR' and e['restrict']=='P' and bool(e['semi_enzymatic'])==(route=='semi')
    assert all(p==search_parameters[0] for p in search_parameters),'Non-digestion search settings differ'
    forms={};reported={};qcs={};source_files=[source_receipt]
    for route,v in route_inputs.items():
        verify(v['search']);source_files.extend([meta(v['search']['path']),v['config']])
        forms[route],reported[route],qcs[route]=psmread(v['search']['path'])
    assert set(reported['unspecific'])==set(Path(oo['peptide_q01.txt']['path']).read_text().splitlines())
    verify(oo['peptide_q01.txt'])
    # The two independent recounters must agree before any biological summary.
    verify(t['baseline_peptides']);assert set(reported['tryptic'])==set(Path(t['baseline_peptides']['path']).read_text().splitlines())
    queries,filtering,kept=load_queries(t['prediction']['path']);kept_ids={i for i,_,_ in kept};dn=Counter();ranks=Counter()
    with open(t['prediction']['path']) as f:
        for i,r in enumerate(csv.DictReader(f),1):
            ranks[r.get('beam_idx',r.get('beam_index',r.get('beam_rank','FIELD_NOT_FOUND')))]+=1
            if i in kept_ids:dn[clean(r.get(COL) or '')]+=1
    collapsed=Counter()
    for seq,n in dn.items():collapsed[seq.replace('I','L')]+=n
    assert collapsed==queries
    forms['denovo']=dn;reported['denovo']={};qcs['denovo']=dict(filtering,beam_ranks=dict(ranks))
    sets={r:{s.replace('I','L') for s in ss} for r,ss in forms.items()};universe=set().union(*sets.values())
    machine=auto(universe);occ=defaultdict(lambda:defaultdict(set));fa={};audit={}
    for h,s in fasta(t['fasta']['path'],audit):
        pid=h.split()[0];assert pid not in fa;fa[pid]=s
        for e,p in machine.iter(s.replace('I','L')):occ[p][pid].add(context(s,e+1-len(p),e+1))
    assert audit['sha256']==receipt['input_fasta']['sha256']
    peptide_rows=[];properties={}
    for route,peps in sets.items():
        props={}
        for p in sorted(peps):
            available={pid:{c for c in cc if eligible(c,route)} for pid,cc in occ[p].items()}
            available={pid:cc for pid,cc in available.items() if cc};members=set(available)
            contexts=set().union(*available.values()) if available else set()
            if route!='denovo':
                assert members,(route,p,'No compatible actual local parent')
                assert reported[route][p] and reported[route][p]<=members,(route,p,'Reported target parent not digest compatible',reported[route][p]-members)
            readable=classify_members(members);unknown=sum(not is_classifiable_accession(pid) for pid in members)
            conservative='Unknown' if unknown else readable
            p1={c[0] for c in contexts};c_end={c[1]=='PROTEIN_END' for c in contexts}
            if not members:conservative=readable='NotInSelectedFASTA'
            row=dict(route=route,peptide=p,length=len(p),origin=conservative,readable_origin=readable,
                n_local_parents=len(members),n_unreadable_parents=unknown,parents=';'.join(sorted(members)),
                context_count=len(contexts),n_p1_options=';'.join(sorted(p1)),
                n_protein_end_any=int('PROTEIN_END' in p1),n_protein_end_all=int(p1=={'PROTEIN_END'}),
                c_protein_end_any=int(True in c_end),c_protein_end_all=int(c_end=={True}),
                fully_tryptic_any=int(any(c[2] and c[3] for c in contexts)),
                fully_tryptic_all=int(bool(contexts) and all(c[2] and c[3] for c in contexts)),
                fully_nontryptic_any=int(any(not c[2] and not c[3] for c in contexts)),
                fully_nontryptic_all=int(bool(contexts) and all(not c[2] and not c[3] for c in contexts)),
                n_histidine_internal_any=int('H' in p1),n_histidine_internal_all=int(p1=={'H'}),
                c_histidine_internal_any=int(p[-1]=='H' and False in c_end),
                c_histidine_internal_all=int(p[-1]=='H' and c_end=={False}))
            peptide_rows.append(row);props[p]=row
        properties[route]=props
    with gzip.open(out/'peptide_contexts.tsv.gz','wt') as f:
        w=csv.DictWriter(f,fieldnames=list(peptide_rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(peptide_rows)
    counts=Counter();denoms=Counter();summaries=[];route_pairs=[]
    for route,ss in forms.items():
        for lo in (5,9):
            if route=='denovo' and lo==5:continue
            for policy in ['conservative','readable_parent_sensitivity']:
                for unit in ['canonical_peptide','assigned_bare_sequence','accepted_psm_or_prediction_row']:
                    items=Counter({p:1 for p in sets[route]}) if unit=='canonical_peptide' else (Counter({s:1 for s in ss}) if unit=='assigned_bare_sequence' else ss)
                    for s,weight in items.items():
                        if len(s)<lo:continue
                        p=s.replace('I','L');v=properties[route][p];origin=v['origin' if policy=='conservative' else 'readable_origin']
                        for cls in ['All',origin]:
                            for term,residues in [('N_first',[s[0]]),('C_last',[s[-1]]),('N_cleavage_P1',v['n_p1_options'].split(';') if v['n_p1_options'] else ['UNMAPPED'])]:
                                key=(route,lo,policy,unit,cls,term);denoms[key]+=weight
                                for aa in residues:
                                    if unit=='canonical_peptide' and aa in ('I','L'):aa='I/L'
                                    counts[key+(aa,)]+=weight/len(residues)
            for cls in ['All','Human','Microbial','Mixed','Unknown','NotInSelectedFASTA']:
                rr=[v for p,v in properties[route].items() if len(p)>=lo and (cls=='All' or v['origin']==cls)]
                summaries.append(dict(route=route,min_length=lo,origin=cls,canonical_peptides=len(rr),
                    **{k:sum(v[k] for v in rr) for k in ['fully_tryptic_any','fully_tryptic_all','fully_nontryptic_any','fully_nontryptic_all','n_histidine_internal_any','n_histidine_internal_all','c_histidine_internal_any','c_histidine_internal_all','n_protein_end_any','c_protein_end_any']}))
            for other in sets:
                if other==route or (other=='denovo' and lo==5):continue
                a={p for p in sets[route] if len(p)>=lo};b={p for p in sets[other] if len(p)>=lo}
                route_pairs.append(dict(test=route,reference=other,min_length=lo,test_count=len(a),reference_count=len(b),overlap=len(a&b),gained=len(a-b),lost=len(b-a)))
    rows=[]
    for key,denom in sorted(denoms.items()):
        assert math.isclose(sum(v for k,v in counts.items() if k[:-1]==key),denom,rel_tol=1e-9)
        alphabet=set(AA)|{'PROTEIN_END','UNMAPPED'}|{k[-1] for k in counts if k[:-1]==key}
        if key[3]=='canonical_peptide':alphabet=(alphabet-{'I','L'})|{'I/L'}
        for aa in sorted(alphabet):
            n=counts[key+(aa,)];rows.append(dict(zip(['route','min_length','origin_policy','unit','origin','terminus'],key),residue=aa,count=n,denominator=denom,fraction=n/denom if denom else None))
    table(out/'terminal_counts.tsv',rows);table(out/'context_summary.tsv',summaries);table(out/'route_identity.tsv',route_pairs)
    dump(out/'COMPLETE.json',dict(status='PASS',cohort=t['cohort'],sample=t['sample'],design=t['design'],
        source_receipt=source_receipt,sources=source_files+[meta(t['prediction']['path']),meta(t['fasta']['path'])],
        search_counts={r:len(s) for r,s in sets.items()},qc=qcs,matched_search_chemistry=True,
        origin_rule='All compatible local parents; unreadable parents force Unknown. Readable-only historical rule is an explicitly labelled sensitivity.',
        n_p1_rule='Equal fractional weight across distinct eligible P1 residues; protein ends and absent contexts separate. Any/all context bounds also retained.',
        prediction_scope='Production top30% retained exported predictions, 9–50, not database-q-filtered',
        job_id=os.environ.get('SLURM_JOB_ID'),outputs=[meta(p) for p in out.iterdir() if p.is_file() and p.name!='COMPLETE.json']))
    print('TERMINAL PASS',t['cohort'],t['sample'],{r:len(s) for r,s in sets.items()},flush=True)

if __name__=='__main__':run(int(sys.argv[2]))
