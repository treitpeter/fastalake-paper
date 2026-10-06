from pathlib import Path
import argparse, csv, json, hashlib, subprocess, os, re, math, datetime, shutil, statistics, collections, tempfile
from lxml import etree as ET

W=Path(__file__).resolve().parent
ROOT=next(p for p in W.parents if (p/'00_START_HERE_V7.md').exists())
CAN=ROOT/'canonical'; ARMS=CAN/'publication_data/arms'
PY='/mpcdf/soft/SLE_15/packages/x86_64/anaconda/3/2023.03/bin/python'
EX=CAN/'rust/fasta_extractor_v2/target/release/fasta_extractor'
PAR=CAN/'rust/parsimony_engine/target/release/parsimony_engine'
SAGE=ROOT/'sample_specific_merge_all_raw/2026-09-11/bin/sage'
NS='http://psi.hupo.org/ms/mzml'; M='{'+NS+'}'
COHORTS=['HSTOOL','CAMPI_MOCK','PREDICT','CAPSCAN']
csv.field_size_limit(10_000_000)
def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def dump(p,j):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);tmp=p.with_suffix(p.suffix+'.tmp')
    tmp.write_text(json.dumps(j,indent=2)+'\n');tmp.replace(p)
def rows(p):return list(csv.DictReader(open(p),delimiter='\t'))
def write_rows(p,rr):
    p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rr[0]),delimiter='\t');w.writeheader();w.writerows(rr)
def command(cmd,prefix):
    prefix=Path(prefix);prefix.parent.mkdir(parents=True,exist_ok=True)
    dump(str(prefix)+'.command.json',dict(command=list(map(str,cmd)),utc=now()))
    with open(str(prefix)+'.out','w') as out,open(str(prefix)+'.err','w') as err:
        subprocess.run(list(map(str,cmd)),stdout=out,stderr=err,check=True)
def check_frozen():
    j=json.loads((W/'MANIFEST.json').read_text())
    for r in j['code']:assert sha(r['path'])==r['sha256'],r['path']
    return j
def setup():
    assert not (W/'MANIFEST.json').exists(),'already prepared'
    rdir=W.parents[1].parent/'main_methods_pass_2026-09-24/review_decisions'
    keepH={r['sample_id'] if 'sample_id' in r else r.get('acquisition_id',r.get('sample')) for r in rows(rdir/'M01_sample_rosters/healthy_campi/HSTOOL_acquisitions.tsv')}
    keepC={r['sample_id'] for r in rows(CAN/'publication_data/cohorts/CAPSCAN/conservative_331.tsv')}
    capmanifest=json.loads((ROOT/'MM_review_analysis/A08_4_5_capscan_public_only/inputs/manifest.json').read_text())
    cm={r['sample']:r for r in capmanifest['samples']}
    records=[];cohorts={}
    for co,listname,num in [('HSTOOL','hstool',58),('CAMPI_MOCK','campi_mock',15),('PREDICT','predict',306),('CAPSCAN','capscan',331)]:
        ro=rows(CAN/f'publication/cohorts/{listname}_sample_list.tsv')
        rr=[]
        for row in ro:
            s=row.get('sample',row.get('sample_id'))
            if co=='HSTOOL' and s not in keepH:continue
            if co=='CAPSCAN' and s not in keepC:continue
            csvpath=row.get('denovo_csv',row.get('alphanovo_csv'));mzml=row['mzml']
            arm={'HSTOOL':'noclust_2026-08-07','CAMPI_MOCK':'uniprotkb_2026_03_noclust_2026-09-09','PREDICT':'abl_nometag_noclust'}.get(co)
            if co=='CAPSCAN':
                base=ROOT/f'MM_review_analysis/A08_4_5_capscan_public_only/results/{s}/search'
                cache=ROOT/'MM_review_analysis/A08_4_5_capscan_public_only/evidence/selected.fasta'
                reservoir=Path(capmanifest['reservoir']['path'])
                csvpath=cm[s]['prediction']['path']
            else:
                base=ARMS/f'{arm}/cohorts/{co}/5_sage/razor_hash-acc/{s}'
                cache=ARMS/f'{arm}/cohorts/{co}/1_evidence/evidence_lake.fasta'
                stats=json.loads(cache.with_name('evidence_lake.stats.json').read_text());reservoir=Path(stats['database'])
            config=base/'results.json'
            cfg=json.loads(config.read_text());fasta=cfg['database']['fasta']
            for p in [csvpath,mzml,config,fasta,cache,reservoir]:assert Path(p).is_file() and Path(p).stat().st_size>0,str(p)
            with open(csvpath) as f:
                hd=next(csv.reader(f));assert {'spec_idx','score','prec_mz','prec_charge','peptide_prediction_detokenized_unmodified'}<=set(hd)
            r=dict(index=len(records)+len(rr),cohort=co,sample=s,csv=csvpath,mzml=mzml,config=str(config),config_sha256=sha(config),fasta=fasta,fasta_sha256=sha(fasta),csv_bytes=Path(csvpath).stat().st_size,csv_mtime_ns=Path(csvpath).stat().st_mtime_ns,fasta_bytes=Path(fasta).stat().st_size,production_cache=str(cache),reservoir=str(reservoir))
            rr.append(r)
        assert len(rr)==num and len({r['sample'] for r in rr})==num,(co,len(rr))
        order=sorted(rr,key=lambda r:r['fasta_bytes'])
        pilot={order[0]['sample'],order[len(order)//2]['sample'],order[-1]['sample']}
        for r in rr:r['pilot']=r['sample'] in pilot
        records.extend(rr);cohorts[co]=dict(n=num,reservoir=str(reservoir),production_cache=str(cache))
    assert len(records)==710
    for p in [EX,PAR,SAGE]:assert p.is_file() and os.access(p,os.X_OK),str(p)
    version=subprocess.check_output([str(SAGE),'--version'],text=True).strip();assert '0.14.6' in version,version
    codepaths=[Path(__file__),W/'PROTOCOL.md',EX,PAR,SAGE]
    for folder in ['logs','jobs','inputs','cache','outputs','results']:(W/folder).mkdir(exist_ok=True)
    write_rows(W/'tasks.tsv',records)
    dump(W/'MANIFEST.json',dict(utc=now(),cohorts=cohorts,tasks=records,sage_version=version,code=[dict(path=str(p),sha256=sha(p),bytes=p.stat().st_size) for p in codepaths]))
    print('PREPARED',len(records),'pilot',sum(r['pilot'] for r in records),flush=True)
def split(co):
    m=check_frozen();rr=[r for r in m['tasks'] if r['cohort']==co];summary=[]
    for r in rr:
        src=Path(r['csv']);assert src.stat().st_size==r['csv_bytes'] and src.stat().st_mtime_ns==r['csv_mtime_ns']
        paths={p:W/f"inputs/{co}/{p}/{r['sample']}_predictions.csv" for p in ['even','odd']}
        for p in paths.values():p.parent.mkdir(parents=True,exist_ok=True)
        fs={k:p.open('w',newline='') for k,p in paths.items()};counts=collections.Counter();ids=set()
        with src.open() as f:
            rd=csv.DictReader(f);ws={k:csv.DictWriter(v,fieldnames=rd.fieldnames) for k,v in fs.items()}
            for w in ws.values():w.writeheader()
            for row in rd:
                n=float(row['spec_idx']);assert n.is_integer() and n>=0
                scan=int(n);assert scan not in ids;ids.add(scan)
                assert float(row.get('beam_rank','0'))==0,'not top-ranked export'
                p='even' if scan%2==0 else 'odd';ws[p].writerow(row);counts[p]+=1
        for f in fs.values():f.close()
        assert min(counts['even'],counts['odd'])>0
        summary.append(dict(sample=r['sample'],csv_sha256=sha(src),counts=dict(counts),half_sha256={k:sha(p) for k,p in paths.items()}))
        print(co,r['sample'],dict(counts),flush=True)
    dump(W/f'inputs/{co}/COMPLETE.json',dict(cohort=co,n=len(rr),samples=summary,utc=now()))
def extract(csvdir,db,out,threads=8):
    out.mkdir(parents=True,exist_ok=True)
    command([EX,'--peptides-dir',csvdir,'--database',db,'--output',out/'selected.fasta','--output-hits',out/'hits.tsv','--top-percent','0.30','--min-length','9','--max-length','50','--normalize-il','true','-t',str(threads)],out/'extract')
    assert (out/'selected.fasta').stat().st_size>0
def cache(index):
    m=check_frozen();co=COHORTS[index//2];parity=['even','odd'][index%2]
    receipt=json.loads((W/f'inputs/{co}/COMPLETE.json').read_text())
    assert receipt['n']==m['cohorts'][co]['n']
    for r in receipt['samples']:assert sha(W/f"inputs/{co}/{parity}/{r['sample']}_predictions.csv")==r['half_sha256'][parity]
    out=W/f'cache/{co}/{parity}';extract(W/f'inputs/{co}/{parity}',m['cohorts'][co]['reservoir'],out,16)
    stats=json.loads((out/'selected.stats.json').read_text())
    dump(out/'COMPLETE.json',dict(cohort=co,parity=parity,stats=stats,fasta_sha256=sha(out/'selected.fasta'),input_receipt_sha256=sha(W/f'inputs/{co}/COMPLETE.json'),utc=now()))
def seqset(path):
    seqs=set();parts=[]
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                if parts:seqs.add(''.join(parts).upper())
                parts=[]
            else:parts.append(line.strip())
    if parts:seqs.add(''.join(parts).upper())
    return seqs
def select(r,parity,db,csvpath):
    out=W/f"outputs/{r['cohort']}/{r['sample']}/selection/{parity}";pool=out/'queries';pool.mkdir(parents=True,exist_ok=True)
    dest=pool/(r['sample']+'_predictions.csv')
    if not dest.exists():dest.symlink_to(csvpath)
    extract(pool,db,out/'retrieval')
    par=out/'parsimony';par.mkdir(exist_ok=True)
    command([PAR,'--strategy','razor','--tiebreak','hash-acc','--fasta',out/'retrieval/selected.fasta','--peptides',csvpath,'--output-dir',par,'--sample',r['sample'],'--min-length','9','--max-length','50','--normalize-il','-t','8'],par/'select')
    choices=list(par.glob('*.fasta'));assert len(choices)==1,[str(p) for p in choices]
    return choices[0]
def odd_mzml(src,csvpath,dest,receipt_path):
    pred={}
    with open(csvpath) as f:
        for row in csv.DictReader(f):pred[int(float(row['spec_idx']))]=(float(row['prec_mz']),int(float(row['prec_charge'])))
    count=0;seen=set();matched=0;ms1=0;odd=0;run_attrs=None;root_attrs=None;root_ns=None;headers=[];slist_attrs=None
    # First pass checks every prediction's native ID and precursor and determines the XML count.
    for event,e in ET.iterparse(src,events=['start','end'],huge_tree=True):
        if event=='start':
            if e.tag==M+'mzML':root_attrs=dict(e.attrib);root_ns=e.nsmap
            if e.tag==M+'run':run_attrs=dict(e.attrib)
            if e.tag==M+'spectrumList':slist_attrs=dict(e.attrib)
            continue
        if e.tag==M+'spectrum':
            cv={x.get('accession'):x.get('value') for x in e.iter(M+'cvParam')}
            level=int(cv['MS:1000511']);scan=int(re.search(r'scan=(\d+)',e.get('id')).group(1))
            if level==2:
                seen.add(scan)
                if scan in pred:
                    mz,ch=pred[scan];assert abs(float(cv['MS:1000744'])-mz)<.001,(scan,'m/z')
                    assert int(cv.get('MS:1000041','0'))==ch,(scan,'charge');matched+=1
            keep=level!=2 or scan%2==1
            count+=keep;ms1+=level==1;odd+=level==2 and scan%2==1
            e.clear()
            while e.getprevious() is not None:del e.getparent()[0]
        elif e.getparent() is not None and e.getparent().tag==M+'mzML' and e.tag!=M+'run':
            headers.append(ET.tostring(e))
    assert matched==len(pred) and set(pred)<=seen
    assert root_attrs is not None and run_attrs is not None and slist_attrs is not None
    slist_attrs['count']=str(count);written=0
    with ET.xmlfile(str(dest),encoding='utf-8') as xf:
        xf.write_declaration()
        with xf.element(M+'mzML',root_attrs,nsmap=root_ns):
            for h in headers:xf.write(ET.fromstring(h))
            with xf.element(M+'run',run_attrs):
                with xf.element(M+'spectrumList',slist_attrs):
                    for _,e in ET.iterparse(src,events=['end'],tag=M+'spectrum',huge_tree=True):
                        level=int(next(x.get('value') for x in e.findall(M+'cvParam') if x.get('accession')=='MS:1000511'))
                        scan=int(re.search(r'scan=(\d+)',e.get('id')).group(1))
                        if level!=2 or scan%2:
                            e.set('index',str(written));xf.write(e);written+=1
                        e.clear()
                        while e.getprevious() is not None:del e.getparent()[0]
    assert written==count and odd>0
    dump(receipt_path,dict(source=str(src),source_bytes=Path(src).stat().st_size,predictions=len(pred),predictions_matched=matched,original_ms2=len(seen),retained_ms1=ms1,retained_odd_ms2=odd,output_spectra=written,output_sha256=sha(dest),native_ids_preserved=True,utc=now()))
def accepted(path,out):
    acc=set();nrows=0;ndecoy=0;scans=set()
    with open(path) as f:
        for r in csv.DictReader(f,delimiter='\t'):
            scan=int(re.search(r'scan=(\d+)',r['scannr']).group(1));assert scan%2==1,'even scan entered search'
            nrows+=1;scans.add(scan);ndecoy+=int(r['label'])==-1
            if int(r['rank'])==1 and int(r['label'])==1 and float(r['peptide_q'])<=.01:
                pep=re.sub(r'\[[^\]]*\]','',r['peptide']).replace('I','L')
                assert pep.isalpha() and pep.isupper(),pep
                if 5<=len(pep)<=50:acc.add(pep)
    assert nrows>0 and ndecoy>0,'missing scores or decoys'
    out.write_text('\n'.join(sorted(acc))+'\n')
    return dict(psms=nrows,decoy_psms=ndecoy,scans=len(scans),peptides_5_50=len(acc),peptides_9_50=sum(len(p)>=9 for p in acc),accepted_sha256=sha(out)),acc
def worker(index):
    m=check_frozen();r=m['tasks'][index];out=W/f"outputs/{r['cohort']}/{r['sample']}";out.mkdir(parents=True,exist_ok=True)
    assert sha(r['config'])==r['config_sha256'] and sha(r['fasta'])==r['fasta_sha256']
    fastas={'FULL':Path(r['fasta'])};full_reproduced=None
    for p,arm in [('even','HALF_DISJOINT'),('odd','HALF_SELF')]:
        ca=W/f"cache/{r['cohort']}/{p}"
        cj=json.loads((ca/'COMPLETE.json').read_text());assert cj['cohort']==r['cohort'] and cj['parity']==p
        cp=W/f"inputs/{r['cohort']}/{p}/{r['sample']}_predictions.csv"
        fastas[arm]=select(r,p,ca/'selected.fasta',cp)
    if r['pilot']:
        full=select(r,'full_reproduction',r['production_cache'],r['csv'])
        a=seqset(full);b=seqset(r['fasta']);full_reproduced=a==b
        dump(out/'FULL_REPRODUCTION.json',dict(equal=full_reproduced,reproduced=len(a),original=len(b),only_reproduced=len(a-b),only_original=len(b-a)))
        assert full_reproduced,'FULL selection does not reproduce original sequence set'
    scratch=out/'odd_only.mzML'
    odd_mzml(r['mzml'],r['csv'],scratch,out/'PARITY.json')
    metrics={};sets={}
    for arm,fasta in fastas.items():
        search=out/arm;search.mkdir(exist_ok=True)
        cfg=json.loads(Path(r['config']).read_text());cfg.pop('version',None)
        cfg['database']['fasta']=str(fasta);cfg['mzml_paths']=[str(scratch)];cfg['output_directory']=str(search)
        dump(search/'config.json',cfg)
        command([SAGE,search/'config.json','--write-pin'],search/'sage')
        metrics[arm],sets[arm]=accepted(search/'results.sage.tsv',search/'ACCEPTED.txt')
        assert all(int(re.search(r'scan=(\d+)',s).group(1))%2 for s in []) # per-row check is in accepted()
    F=sets['FULL'];S=sets['HALF_SELF'];D=sets['HALF_DISJOINT']
    result=dict(index=index,cohort=r['cohort'],sample=r['sample'],pilot=r['pilot'],metrics=metrics,full_reproduced=full_reproduced,disjoint_over_self=len(D)/len(S) if S else None,disjoint_over_full=len(D)/len(F) if F else None,self_over_full=len(S)/len(F) if F else None,disjoint_in_self=len(D&S)/len(D) if D else None,disjoint_in_full=len(D&F)/len(D) if D else None,whole_procedure_fdr_control_established=False,utc=now())
    dump(out/'COMPLETE.json',result)
    scratch.unlink() # reproducible derivative only, after receipt + all successful searches
    print(json.dumps(result),flush=True)
def summarize(require_all=True):
    m=check_frozen();reports={}
    for co in COHORTS:
        tasks=[r for r in m['tasks'] if r['cohort']==co];valid=[];missing=[]
        for r in tasks:
            p=W/f"outputs/{co}/{r['sample']}/COMPLETE.json"
            if p.exists():valid.append(json.loads(p.read_text()))
            else:missing.append(r['sample'])
        if require_all:assert not missing,(co,len(missing))
        rep=dict(cohort=co,expected=len(tasks),completed=len(valid),missing=missing,complete=not missing,whole_procedure_fdr_control_established=False,utc=now())
        for k in ['disjoint_over_self','disjoint_over_full','self_over_full','disjoint_in_self','disjoint_in_full']:
            vals=[r[k] for r in valid if r[k] is not None];rep['median_'+k]=statistics.median(vals) if vals else None
        dump(W/f'results/{co}_SUMMARY.json',rep);reports[co]=rep
        if valid:write_rows(W/f'results/{co}_per_acquisition.tsv',[{k:v for k,v in r.items() if k!='metrics'}|{arm+'_'+k:v for arm,x in r['metrics'].items() for k,v in x.items() if k.startswith('peptides')} for r in valid])
    print(json.dumps(reports),flush=True)
if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('mode',choices=['setup','split','cache','worker','summary']);a.add_argument('index',nargs='?');x=a.parse_args()
    if x.mode=='setup':setup()
    elif x.mode=='split':split(COHORTS[int(x.index)])
    elif x.mode=='cache':cache(int(x.index))
    elif x.mode=='worker':worker(int(x.index))
    else:summarize()
