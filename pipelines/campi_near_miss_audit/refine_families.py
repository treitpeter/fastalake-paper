#!/usr/bin/env python3
import csv,gzip,json,pathlib,shutil,sys,os,time
from collections import Counter,defaultdict
RUN=pathlib.Path(sys.argv[1]).resolve()
sys.path.insert(0,str(RUN/'frozen'))
from lib import dump,table,meta,sha,verify
from analysis import records,automaton,family_status
M=json.loads((RUN/'derived/PREPARED.json').read_text())
V=json.loads((RUN/'inputs/family_refinement_manifest.json').read_text())
def readfa_raw(p):
    out={};header=None;buf=[]
    for line in open(p):
        if line.startswith('>'):
            if header is not None:
                assert header not in out;out[header]=''.join(buf).upper()
            header=line[1:].split()[0];buf=[]
        elif line.strip():buf.append(line.strip())
    if header is not None:assert header not in out;out[header]=''.join(buf).upper()
    return out
def main():
    verify(V['source']);verify(V['frozen']);verify(V['script'])
    for r in M['common'].values():verify(r)
    refseq=readfa_raw(M['common']['reference_fasta']['path'])
    appendseq=readfa_raw(V['frozen']['path'])
    values=defaultdict(set);aliases=[]
    for key in ['annotation_blind','annotation_ref']:
        for r in csv.DictReader(open(M['common'][key]['path']),delimiter='\t'):
            fs=frozenset(x.split('@')[0] for x in r['OG'].split(',') if '@2|Bacteria' in x)
            acc=r['accession']
            if acc.startswith('CAMPI_REF__'):
                raw=acc[len('CAMPI_REF__'):]
                assert raw in refseq and acc in appendseq
                assert refseq[raw]==appendseq[acc],('Reference alias sequence mismatch',acc)
                aliases.append(dict(annotation_accession=acc,reference_accession=raw,exact_sequence_sha256=__import__('hashlib').sha256(refseq[raw].encode()).hexdigest(),family_annotated=bool(fs)))
                acc=raw
            if fs:values[acc].add(fs)
    assert len(aliases)>25000,'Reference alias mapping unexpectedly incomplete'
    conflicts={k for k,v in values.items() if len(v)>1}
    ann={k:set(next(iter(v))) for k,v in values.items() if len(v)==1}
    table(RUN/'results/reference_annotation_aliases.tsv',aliases)
    basecomplete=RUN/'results/COMPLETE.json'
    original=RUN/'status/FAMILY_REFINEMENT_ORIGINAL_COMPLETE.json'
    assert not original.exists(),'Refinement already attempted; inspect state'
    shutil.copy2(basecomplete,original)
    for sm in M['samples']:
        out=RUN/'results'/sm['sample'];old=json.loads((out/'COMPLETE.json').read_text())
        for f in old['outputs']:assert sha(f['path'])==f['sha256'],f['path']
        archive=out/'before_family_alias_audit';archive.mkdir()
        filenames=['COMPLETE.json','query_lookup.tsv','reference_retention.tsv','query_reference_pairs.tsv','spectrum_prediction_fate.tsv.gz','summary.tsv']
        for name in filenames:shutil.copy2(out/name,archive/name)
        qr=list(csv.DictReader((archive/'query_lookup.tsv').open(),delimiter='\t'))
        qs={r['query']:dict(candidate_proteins=int(r['candidate_proteins']),annotated_candidates=0,candidate_families=set(),selected_proteins=int(r['selected_proteins']),annotated_selected=0,selected_families=set()) for r in qr}
        h=__import__('hashlib').sha256();counts=Counter();last=None;seen=set()
        with open(sm['hits']['path'],'rb') as f:
            raw=f.readline();h.update(raw)
            for raw in f:
                h.update(raw);q,p,_=raw.decode().rstrip('\n').split('\t')
                if p!=last:last=p;seen=set()
                if q in seen:continue
                seen.add(q);counts[q]+=1
                if p in ann:qs[q]['annotated_candidates']+=1;qs[q]['candidate_families'].update(ann[p])
        assert h.hexdigest()==old['archived_hits_sha256'],'Source hit table changed'
        assert counts==Counter({q:r['candidate_proteins'] for q,r in qs.items() if r['candidate_proteins']})
        rr=list(csv.DictReader((archive/'reference_retention.tsv').open(),delimiter='\t'))
        refs={r['reference_peptide'] for r in rr};rf=defaultdict(set);rp=Counter();rna=Counter();a=automaton(refs)
        for p,s in refseq.items():
            for q in {q for _,q in a.iter(s.replace('I','L'))}:
                rp[q]+=1
                if p in ann:rf[q].update(ann[p]);rna[q]+=1
        a=automaton(qs)
        observed_selected=Counter()
        for p,s in records(sm['final_fasta']['path']):
            for q in {q for _,q in a.iter(s)}:
                observed_selected[q]+=1
                if p in ann:qs[q]['annotated_selected']+=1;qs[q]['selected_families'].update(ann[p])
        assert observed_selected==Counter({q:r['selected_proteins'] for q,r in qs.items() if r['selected_proteins']})
        for r in rr:
            q=r['reference_peptide'];assert rp[q]==int(r['reference_parents'])
            r['annotated_reference_parents']=rna[q];r['reference_families']=';'.join(sorted(rf[q]))
        table(out/'reference_retention.tsv',rr)
        for r in qr:
            q=r['query']
            for key in ['annotated_candidates','candidate_families','annotated_selected','selected_families']:
                value=qs[q][key];r[key]=';'.join(sorted(value)) if isinstance(value,set) else value
        table(out/'query_lookup.tsv',qr)
        def statuses(q,r):
            if q not in qs:return None,None,None,None
            st=qs[q];ac,ec=family_status(st,rf[r])
            asc,esc=family_status(dict(candidate_proteins=st['selected_proteins'],annotated_candidates=st['annotated_selected'],candidate_families=st['selected_families']),rf[r])
            return ac,ec,asc,esc
        pairs=list(csv.DictReader((archive/'query_reference_pairs.tsv').open(),delimiter='\t'))
        for r in pairs:
            q,p=r['query'],r['reference_peptide'];r['reference_family_known']=bool(rf[p])
            r['annotated_candidates']=qs[q]['annotated_candidates'] if q in qs else ''
            ac,ec,asc,esc=statuses(q,p)
            r.update(family_any_compatible=ac,family_exclusively_compatible=ec,selected_family_any_compatible=asc,selected_family_exclusively_compatible=esc)
        table(out/'query_reference_pairs.tsv',pairs)
        sums=defaultdict(Counter);metrics=['reference_present_final','reference_identified_fl_anywhere','reference_identified_fl_same_spectrum','same_candidate_contains_reference','same_selected_contains_reference','family_any_compatible','family_exclusively_compatible','selected_family_any_compatible','selected_family_exclusively_compatible']
        n=0
        with gzip.open(archive/'spectrum_prediction_fate.tsv.gz','rt') as fi,gzip.open(out/'spectrum_prediction_fate.tsv.gz','wt',newline='') as fo:
            rd=csv.DictReader(fi,delimiter='\t');fields=rd.fieldnames+['selected_family_any_compatible','selected_family_exclusively_compatible','annotated_selected_proteins']
            w=csv.DictWriter(fo,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader()
            for row in rd:
                n+=1;row.update(selected_family_any_compatible='',selected_family_exclusively_compatible='',annotated_selected_proteins='')
                if row['reference_accepted']=='1':
                    q,r=row['prediction_IL'],row['reference_peptide'];row['reference_family_known']=int(bool(rf[r]))
                    if q in qs:
                        row['annotated_candidates']=qs[q]['annotated_candidates'];row['annotated_selected_proteins']=qs[q]['annotated_selected']
                        vals=statuses(q,r)
                        for k,v in zip(metrics[5:],vals):row[k]='' if v is None else int(v)
                    pops=['all_reference_paired']
                    if row['production_row_retained']=='1':
                        pops.append('production_retained_reference_paired')
                        if row['also_spectrum_q_01']=='1':pops.append('production_retained_reference_paired_both_q01')
                    for pop in pops:
                        z=sums[pop,row['category']];z['n']+=1
                        for k in metrics:
                            if row[k]!='':z[k+'_evaluable']+=1;z[k]+=int(row[k])
                        if row['candidate_proteins']!='':
                            z['any_lookup']+=int(int(row['candidate_proteins'])>0)
                            z['all_candidates_annotated']+=int(int(row['candidate_proteins'])>0 and int(row['candidate_proteins'])==int(row['annotated_candidates']))
                w.writerow(row)
        assert n==old['prediction_rows']
        sr=[dict(sample=sm['sample'],population=p,category=c,**dict(z)) for (p,c),z in sorted(sums.items())]
        table(out/'summary.tsv',sr,sorted(set().union(*(set(r) for r in sr))))
        old['reference_peptides_with_family']=sum(bool(rf[r]) for r in refs);old['annotation_conflicts']=len(conflicts);old['family_alias_sequence_validation']='PASS';old['selected_family_audit']='PASS';old['refinement_job_id']=os.environ.get('SLURM_JOB_ID')
        old['outputs']=[meta(out/f) for f in filenames if f!='COMPLETE.json']+[meta(out/'scan_mapping.tsv')]
        dump(out/'COMPLETE.json',old)
        print('Family audit validated',sm['sample'],old['reference_peptides_with_family'],'/',len(refs),flush=True)
    import importlib.util
    spec=importlib.util.spec_from_file_location('report_v2',RUN/'frozen/report_v2.py');module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    module.build(RUN,json.loads((RUN/'inputs/manifest.json').read_text()))
    comp=json.loads((RUN/'results/COMPLETE.json').read_text());z=comp['nonexact_retained']
    extra='\n\nFamily alias reconciliation verified every prefixed reference accession against its exact, unmodified protein sequence. The initial direct-ID join and its original output tables remain archived. Selected-family compatibility was also measured explicitly after parsimony: '+str(z.get('selected_family_any_compatible',0))+'/'+str(z['n'])+' non-exact retained reference-paired rows have a selected reference-compatible family. Exclusive selected-family compatibility: '+str(z.get('selected_family_exclusively_compatible',0))+'/'+str(z['n'])+'.\n'
    with (RUN/'results/REPORT.md').open('a') as f:f.write(extra)
    comp['family_refinement']=dict(status='PASS',reference_aliases_verified=len(aliases),conflicting_families=len(conflicts),manifest=meta(RUN/'inputs/family_refinement_manifest.json'),job_id=os.environ.get('SLURM_JOB_ID'))
    comp['outputs']=[meta(p) for p in (RUN/'results').iterdir() if p.is_file() and p.name!='COMPLETE.json']
    dump(RUN/'results/COMPLETE.json',comp)
    dump(RUN/'status/family_refinement.json',dict(status='PASS',utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),job_id=os.environ.get('SLURM_JOB_ID'),aliases=len(aliases)))
if __name__=='__main__':main()
