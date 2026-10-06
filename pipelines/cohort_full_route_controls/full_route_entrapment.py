"""Execute a reservoir-entry control through the unclustered production recipe.

This runner calls the pinned production extractor, razor and Sage binaries. It
does not use the extractor's --entrapment option, which removes control hits from
the output. Evaluation classes live in a sidecar, not selection parameters.
Outputs are descriptive control-labelled discoveries, never a paired FDP.

Usage: python -m fasta_lake.full_route_entrapment --config experiment.json
       {evidence,sample} --cohort NAME [--index N]
"""
from __future__ import annotations
import argparse, collections, copy, csv, hashlib, json, os, re, shutil, subprocess, time
from pathlib import Path

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(4*1024*1024),b''):h.update(block)
    return h.hexdigest()

def fasta(path):
    acc=None;seq=[]
    with Path(path).open() as f:
        for line in f:
            if line.startswith('>'):
                if acc is not None:yield acc,''.join(seq)
                acc=line[1:].split()[0];seq=[]
            else:seq.append(line.strip())
    if acc is not None:yield acc,''.join(seq)

def write_json(path,obj):
    path=Path(path);tmp=path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(obj,indent=2));tmp.replace(path)

def command(args,folder):
    folder=Path(folder);folder.mkdir(parents=True,exist_ok=True)
    record={'argv':list(map(str,args)),'started_unix':time.time()}
    write_json(folder/'command.json',record)
    with (folder/'stdout.log').open('w') as out,(folder/'stderr.log').open('w') as err:
        result=subprocess.run(record['argv'],stdout=out,stderr=err)
    record.update(returncode=result.returncode,elapsed_seconds=time.time()-record['started_unix'])
    write_json(folder/'command.json',record)
    if result.returncode:raise RuntimeError('Command failed; see '+str(folder))

def control_map(path):
    with Path(path).open() as f:return {r['accession']:r['class'] for r in csv.DictReader(f,delimiter='\t')}

def audit_fasta(path,classes,output):
    counts=collections.Counter();retained=[];total=0
    for acc,seq in fasta(path):
        total+=1
        if acc in classes:
            counts[classes[acc]]+=1;retained.append((acc,classes[acc],len(seq)))
        elif acc.startswith('FLX'):raise ValueError('Unmapped control-like accession '+acc)
        else:counts['target']+=1
    with Path(output).open('w',newline='') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(['accession','class','length']);w.writerows(retained)
    return {'total_sequences':total,'counts':dict(counts),'fasta_sha256':sha(path),'control_accessions':str(output)}

def assert_binary(config,key):
    r=config['binaries'][key]
    if sha(r['path'])!=r['sha256']:raise RuntimeError('Pinned binary changed: '+key)
    return r['path']

def extractor(config,peptides,database,out,threads):
    return [assert_binary(config,'extractor'),'--peptides-dir',str(peptides),'--database',str(database),'--output',str(out/'selected.fasta'),'--output-hits',str(out/'hits.tsv'),'--top-percent','0.30','--min-length','9','--max-length','50','--normalize-il','true','-t',str(threads)]

def evidence(config,cohort):
    root=Path(config['output_root'])/cohort;out=root/'1_evidence';out.mkdir(parents=True,exist_ok=True)
    if (out/'COMPLETE.json').exists():raise FileExistsError('Evidence already complete; use a new experiment output')
    c=config['cohorts'][cohort];rows=list(csv.DictReader(open(c['sample_list']),delimiter='\t'))
    assert len(rows)==c['n'] and len({r['sample_id'] for r in rows})==c['n']
    assert sha(c['sample_list'])==c['sample_list_sha256'] and sha(config['controls_fasta'])==config['controls_sha256']
    pool=out/'predictions';pool.mkdir(exist_ok=True)
    for r in rows:
        p=Path(r['denovo_csv']);assert p.is_file()
        st=p.stat();assert (st.st_size,st.st_mtime_ns)==(int(r['denovo_bytes']),int(r['denovo_mtime_ns']))
        q=pool/(r['sample_id']+'_predictions.csv')
        if not q.exists():q.symlink_to(p)
    # The input is the ORIGINAL fixed stage-0 reservoir, never an evidence lake.
    original=Path(c['reservoir']);before=original.stat()
    assert '/0_lake/' in str(original),'Experiment must declare an unselected stage-0 input'
    augmented=root/'0_augmented.fasta'
    with augmented.open('wb') as dest:
        for source in [original,Path(config['controls_fasta'])]:
            with source.open('rb') as src:shutil.copyfileobj(src,dest,length=8*1024*1024)
    after=original.stat();assert (before.st_size,before.st_mtime_ns)==(after.st_size,after.st_mtime_ns)
    command(extractor(config,pool,augmented,out,16),out)
    classes=control_map(config['controls_map']);audit=audit_fasta(out/'selected.fasta',classes,out/'controls.tsv')
    stats=json.loads((out/'selected.stats.json').read_text())
    assert stats['peptide_length_range']==[9,50] and stats['score_filter']=='top_30%_per_sample'
    assert stats['proteins_searched']==c['reservoir_sequences']+config['control_sequences']
    audit.update(cohort=cohort,n_pooled=len(rows),original_reservoir=str(original),original_bytes=before.st_size,original_mtime_ns=before.st_mtime_ns,controls_sha256=config['controls_sha256'],extractor_stats=stats)
    write_json(out/'COMPLETE.json',audit)
    # The augmented file is exactly reproducible from the immutable source + bank.
    # Keep it until the entire experiment is complete; never discard candidate FASTAs.

def canonical(peptide):
    p=re.sub(r'\[[^\]]*\]','',peptide).replace('I','L')
    if not re.fullmatch('[A-Z]+',p):raise ValueError('Unsupported peptide notation: '+peptide)
    return p

def summarize_search(result,classes,out):
    # Store canonical accepted PSM memberships once; derive several q-threshold views.
    records=[]
    with Path(result).open() as f:
        for r in csv.DictReader(f,delimiter='\t'):
            if r['label']!='1':continue
            q=float(r['peptide_q'])
            if not q<=.05:continue
            p=canonical(r['peptide']);assert 5<=len(p)<=50
            labels=set()
            for acc in re.split('[;,]',r['proteins']):
                acc=acc.strip()
                if acc.startswith('rev_'):continue
                if acc.startswith('FLX') and acc not in classes:raise ValueError('Lost evaluation mapping: '+acc)
                labels.add(classes.get(acc,'target'))
            if not labels:raise ValueError('Target-labelled row without target membership')
            records.append((p,q,labels))
    summary=[];discovery=[]
    for threshold in [.001,.005,.01,.02,.05]:
        labels=collections.defaultdict(set)
        for p,q,ll in records:
            if q<=threshold:labels[p].update(ll)
        for floor in [5,9]:
            data={p:v for p,v in labels.items() if len(p)>=floor};count=collections.Counter()
            for p,ll in data.items():
                group=next(iter(ll)) if len(ll)==1 else 'shared'
                count[group]+=1
                if floor==5:discovery.append(dict(peptide=p,q_threshold=threshold,assignment=group,classes=';'.join(sorted(ll))))
            n=len(data)
            row={'q_threshold':threshold,'min_residues':floor,'n_canonical':n,'n_target_only':count['target'],'n_shared':count['shared']}
            for cls in ['shuffle','plant','thermophilic']:
                row['n_'+cls+'_only']=count[cls];row[cls+'_only_fraction']=count[cls]/n if n else None
            summary.append(row)
    for name,data in [('discoveries',discovery),('summary',summary)]:
        with (Path(out)/(name+'.tsv')).open('w',newline='') as f:
            fields=list(data[0]) if data else ['peptide','q_threshold','assignment','classes']
            w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(data)
    return summary

def sample(config,cohort,index):
    c=config['cohorts'][cohort];rows=list(csv.DictReader(open(c['sample_list']),delimiter='\t'));r=rows[index]
    st=Path(r['denovo_csv']).stat();assert (st.st_size,st.st_mtime_ns)==(int(r['denovo_bytes']),int(r['denovo_mtime_ns']))
    root=Path(config['output_root'])/cohort;ev=root/'1_evidence';assert (ev/'COMPLETE.json').is_file()
    out=root/'samples'/r['sample_id'];out.mkdir(parents=True,exist_ok=True)
    if (out/'COMPLETE.json').exists():raise FileExistsError('Sample already complete; refuse overwrite')
    classes=control_map(ev/'controls.tsv');selection=out/'3_selection';selection.mkdir(exist_ok=True)
    pool=selection/'predictions';pool.mkdir(exist_ok=True);link=pool/(r['sample_id']+'_predictions.csv')
    if not link.exists():link.symlink_to(r['denovo_csv'])
    command(extractor(config,pool,ev/'selected.fasta',selection,8),selection)
    sel=audit_fasta(selection/'selected.fasta',classes,selection/'controls.tsv')
    par=out/'4_razor';par.mkdir(exist_ok=True)
    # Match production: parsimony has no top-percent argument and can use all
    # qualifying predictions to assign evidence among candidates admitted above.
    args=[assert_binary(config,'parsimony'),'--strategy','razor','--tiebreak','hash-acc','--fasta',str(selection/'selected.fasta'),'--peptides',r['denovo_csv'],'--output-dir',str(par),'--sample',r['sample_id'],'--min-length','9','--max-length','50','--normalize-il','-t','8']
    command(args,par);db=par/(r['sample_id']+'_razor.fasta');assert db.is_file()
    razor=audit_fasta(db,classes,par/'controls.tsv')
    assert {r['accession'] for r in csv.DictReader((par/'controls.tsv').open(),delimiter='\t')} <= {r['accession'] for r in csv.DictReader((selection/'controls.tsv').open(),delimiter='\t')}
    search=out/'5_sage';search.mkdir(exist_ok=True)
    cfg=json.loads(Path(r['sage_config']).read_text());cfg['database']['fasta']=str(db);cfg['output_directory']=str(search);cfg['mzml_paths']=[r['mzml']]
    write_json(search/'config.json',cfg)
    command([assert_binary(config,'sage'),str(search/'config.json'),'--write-pin'],search)
    result=search/'results.sage.tsv';assert result.is_file()
    summary=summarize_search(result,classes,search)
    write_json(out/'COMPLETE.json',dict(cohort=cohort,sample=r['sample_id'],selection=sel,razor=razor,search_summary=summary,paired_estimate=None,full_procedure_calibration_claimed=False))

def main():
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('stage',choices=['evidence','sample']);p.add_argument('--cohort',required=True);p.add_argument('--index',type=int)
    a=p.parse_args();config=json.loads(Path(a.config).read_text())
    if a.stage=='evidence':evidence(config,a.cohort)
    else:sample(config,a.cohort,a.index)
if __name__=='__main__':main()
