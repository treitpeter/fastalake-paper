"""Native abundance cross-omics; all outputs isolated beside this script."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
os.environ.setdefault('MPLBACKEND','Agg')
import argparse, csv, hashlib, json, math, re, time
from collections import Counter
from itertools import combinations
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.spatial.distance import pdist, squareform
from scipy.stats import rankdata, spearmanr

HERE=Path(__file__).resolve().parent
ROOT=Path(__import__('os').environ['FASTALAKE_DATA_DIR'])
OLD=ROOT/'analyses/enterotypes_2026-09-07'
MG=ROOT/'analyses/capscan_metag_analysis_2026-09-10'
A33=ROOT/'MM_review_analysis/A33_4_6_peptide_centric_downstream_2026-09-19'
CURRENT=A33/'results/crossomics_repaired__simplified/source_data'
DNA_RESULT=ROOT/'fable_5.1_review_2026-09-16/scoping_review_2026-09-19/figure_comments_2026-09-21/capscan_dna/results'
REG=['duodenum','jejunum','ileum','colon','stool']
SEED=20260907

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(8<<20),b''):h.update(b)
    return h.hexdigest()
def rec(path):
    p=Path(path).resolve();s=p.stat()
    return dict(path=str(p),bytes=s.st_size,mtime_ns=s.st_mtime_ns)
def check(r):
    p=Path(r['path']);s=p.stat()
    assert (s.st_size,s.st_mtime_ns)==(r['bytes'],r['mtime_ns']),('Changed input',str(p))
    return p
def dump(path,obj):
    Path(path).parent.mkdir(parents=True,exist_ok=True)
    Path(path).write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')
def read(path,**kw):return pd.read_csv(path,sep='\t',**kw)
def relative(frame):return frame.div(frame.sum(1).replace(0,np.nan),axis=0)
def manifest():return json.loads((HERE/'INPUT_MANIFEST.json').read_text())
def group_path(cohort):
    base=ROOT/('simplified_grouping_results/manuscript_cohorts_2026-09-10_2200' if cohort=='CAPSCAN' else 'MM_review_analysis/recovery_2026-09-15/biological_grouping_recovery/full')
    return base/'cohorts'/(cohort+'_noclust')/'reaggregated'

def prepare():
    assert not (HERE/'INPUT_MANIFEST.json').exists(),'Inspect existing manifest before replacing it'
    pred=read(OLD/'PREDICT_fullstack_primary_samples.tsv',index_col=0)
    assert len(pred)==77 and pred.idpatient.nunique()==77
    tpm=read(OLD/'TPM_input_manifest.tsv').set_index('MicrobPred_ID')
    items={'PREDICT':[],'CAPSCAN':[]}
    for mid in pred.index:
        p=Path(tpm.loc[mid,'path']);asm=p.with_name('predictedGenesFromAssembly.'+p.name.replace('_TPM_annotated.tsv','.fasta'))
        items['PREDICT'].append(dict(id=mid,tpm=rec(p),fasta=rec(asm),specimen=mid,participant=int(pred.loc[mid,'idpatient']),region='stool'))
    cap=read(DNA_RESULT/'DNA_donor_region_roster.tsv');donors=sorted(cap.subject_id.unique().tolist())
    assert donors==[1,2,3,4,6,8,9,10,11,12,13,14] and len(cap)==60
    design=read(OLD/'CAPSCAN_specimen_design.tsv',index_col=0);design.specimen=design.specimen.astype(str)
    met=read(OLD/'CAPSCAN_metabolite_intensity.tsv',index_col=0);met.index=met.index.astype(str)
    og=read(CURRENT/'CAPSCAN_noclust/microbial_OG.tsv',index_col=0).reindex(design.index)
    eligible_specimens=met.index.intersection(design.loc[og.sum(1)>0,'specimen'].unique())
    dspec=design.groupby('specimen').first().loc[eligible_specimens]
    dspec=dspec[dspec.subject_id.isin(donors)&dspec.group.isin(REG)]
    runs=read(MG/'manifests/runs.tsv');runs.specimen=runs.specimen.astype(str)
    dna_saved=read(MG/'aggregate/metaG_OG.tsv',index_col=0)
    good=set(dna_saved.index[dna_saved.sum(1)>0])
    for row in runs[runs.specimen.isin(dspec.index)&runs.run_accession.isin(good)].itertuples():
        source=MG/'results'/row.run_accession
        complete=source/('genes.complete.json' if (source/'genes.complete.json').exists() else 'all.complete.json')
        assert json.loads(complete.read_text())
        d=dspec.loc[row.specimen]
        items['CAPSCAN'].append(dict(id=row.run_accession,specimen=row.specimen,participant=int(d.subject_id),region=d.group,tpm=rec(source/'salmon_quant/quant.sf'),fasta=rec(source/'genes.faa'),completion=rec(complete)))
    assert set((r['participant'],r['region']) for r in items['CAPSCAN'])==set((d,r) for d in donors for r in REG)
    shared=[OLD/'PREDICT_fullstack_primary_samples.tsv',OLD/'CAPSCAN_specimen_design.tsv',OLD/'PREDICT_metabolite_intensity.tsv',OLD/'CAPSCAN_metabolite_intensity.tsv',DNA_RESULT/'DNA_donor_region_roster.tsv',DNA_RESULT/'capscan_profile_correlations_with_DNA.tsv',CURRENT/'crossomics/fullstack_profile_correlations.tsv']
    for co in items:
        g=group_path(co);receipt=json.loads((g/'SUMMARY.json').read_text())
        assert receipt['status']=='PASS' and all(receipt['checks'].values())
        shared.extend([g/'simplified_study_group_matrix.tsv.gz',g/'SUMMARY.json',CURRENT/(co+'_noclust')/'host_representative_intensity.tsv'])
    shared.append(CURRENT/'CAPSCAN_noclust/microbial_OG.tsv')
    obj=dict(created_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),code_sha256=sha(__file__),items=items,shared=[rec(p) for p in shared],capscan_donors=donors,sequence_coordinate='SHA256 of complete encoded amino-acid sequence, terminal stop removed, no similarity clustering',protocol_sha256=sha(HERE/'PROTOCOL.md'))
    dump(HERE/'INPUT_MANIFEST.json',obj)
    for co,rows in items.items():
        read_rows=[{k:v for k,v in r.items() if k not in ['tpm','fasta']} for r in rows]
        pd.DataFrame(read_rows).to_csv(HERE/(co+'_source_roster.tsv'),sep='\t',index=False)
    print(json.dumps({c:dict(sources=len(v),source_GB=sum(x[k]['bytes'] for x in v for k in ['tpm','fasta'])/1e9) for c,v in items.items()}),flush=True)

def fasta(path):
    key=None;parts=[]
    with Path(path).open() as f:
        for line in f:
            if line.startswith('>'):
                if key is not None:yield key,''.join(parts).upper().removesuffix('*')
                key=line[1:].split()[0];parts=[]
            else:parts.append(line.strip())
        if key is not None:yield key,''.join(parts).upper().removesuffix('*')

def source_vector(co,row,out):
    folder=out/'sources';folder.mkdir(parents=True,exist_ok=True)
    cache=folder/(row['id']+'.npz');receipt=folder/(row['id']+'.json')
    if receipt.exists():
        r=json.loads(receipt.read_text());assert r['source']==row and r['status']=='PASS'
        assert sha(cache)==r['output_sha256'];return cache,r
    cols=['prodigal_protein','metaG_tpm','metaT_tpm'] if co=='PREDICT' else ['Name','TPM']
    table=read(check(row['tpm']),usecols=cols).set_index(cols[0])
    assert table.index.is_unique and table.index.notna().all()
    table=table.apply(pd.to_numeric,errors='raise');missing=int(table.isna().sum().sum());table=table.fillna(0)
    assert np.isfinite(table.values).all() and (table.values>=0).all()
    values=table.values;ids={k:i for i,k in enumerate(table.index)};found=set();summed={};gene_hash={};seen=set()
    for gene,seq in fasta(check(row['fasta'])):
        assert gene not in seen,('Duplicate FASTA ID',gene);seen.add(gene)
        if gene not in ids:continue
        assert seq and re.fullmatch('[ACDEFGHIKLMNPQRSTVWYBXZJUO]+',seq),(gene,seq[:80])
        found.add(gene);h=hashlib.sha256(seq.encode()).hexdigest();gene_hash[gene]=h
        v=values[ids[gene]]
        if np.any(v>0):
            if h not in summed:summed[h]=v.copy()
            else:summed[h]+=v
    missing_genes=set(ids)-found
    if co=='CAPSCAN':assert not missing_genes,('TPM genes absent from verified Salmon reference',row['id'],len(missing_genes))
    excluded=table.loc[sorted(missing_genes)]
    excluded.to_csv(folder/(row['id']+'_unresolved_TPM_genes.tsv.gz'),sep='\t')
    excluded_sum=excluded.sum(0).values
    keys=sorted(summed);v=np.array([summed[h] for h in keys],float)
    assert np.allclose(v.sum(0)+excluded_sum,values.sum(0),rtol=1e-10,atol=1e-6)
    assert (v.sum(0)>0).all(),('Empty layer',row['id'])
    np.savez_compressed(cache,hashes=np.array(keys,dtype='S64'),values=v)
    # Bounded identity witness: same assembly ID across samples need not be the same sequence.
    witness={g:gene_hash[g] for g in table.index[:1000] if g in gene_hash}
    r=dict(status='PASS',source=row,genes=len(table),positive_genes=int(np.any(values>0,axis=1).sum()),distinct_positive_sequences=len(keys),abundance_sum=values.sum(0).tolist(),unresolved_gene_count=len(missing_genes),unresolved_abundance_sum=excluded_sum.tolist(),sequence_resolved_TPM_fraction=(v.sum(0)/values.sum(0)).tolist(),missing_value_count=missing,tpm_sha256=sha(check(row['tpm'])),fasta_sha256=sha(check(row['fasta'])),output_sha256=sha(cache),identity_witness=witness)
    if co=='CAPSCAN':
        check(row['completion'])
        root=Path(row['fasta']['path']).parent
        for filename,field in [('04_prodigal.done.json','fasta'),('06_salmon_quant.done.json','tpm')]:
            checkpoint=json.loads((root/filename).read_text())['output_sha256']
            matched=[digest for path,digest in checkpoint.items() if Path(path).resolve()==Path(row[field]['path']).resolve()]
            assert matched==[r[field+'_sha256']],('Checkpoint hash mismatch',row['id'],field)
    dump(receipt,r);print(co,row['id'],len(keys),'sequences',flush=True)
    return cache,r

def dense_control():
    # Independent dense Hellinger and balanced-regression checks for sparse geometry.
    x=np.random.default_rng(93).integers(0,8,size=(15,50)).astype(float);x[:,0]=0
    distances,quality=profile_distance(sparse.csr_matrix(x),'x',.1,capscan=True)
    keep=(x>0).mean(0)>=.1;y=x[:,keep];y=np.sqrt(y/y.sum(1,keepdims=True))
    design=np.eye(5)[np.tile(np.arange(5),3)];res=y-design @ np.linalg.lstsq(design,y,rcond=None)[0]
    assert np.allclose(distances,squareform(pdist(res)),atol=1e-12)

def pilot():
    m=manifest();assert m['code_sha256']==sha(__file__)
    results={}
    for co in ['CAPSCAN','PREDICT']:
        out=HERE/'pilot'/co;out.mkdir(parents=True,exist_ok=True)
        rows=m['items'][co];chosen=[rows[0],next(r for r in rows if r['participant']!=rows[0]['participant'])]
        rr=[source_vector(co,row,out)[1] for row in chosen]
        a,b=[r['identity_witness'] for r in rr];shared=set(a)&set(b)
        differing=sum(a[k]!=b[k] for k in shared)
        results[co]=dict(sources=[r['source']['id'] for r in rr],reused_local_names=len(shared),reused_names_different_sequences=differing,sequence_mapping_and_abundance_accounting='PASS',sequence_resolved_TPM_fraction=[r.get('sequence_resolved_TPM_fraction',[1.]) for r in rr])
        assert differing>0,('Pilot should demonstrate local gene naming',co)
    dense_control()
    dump(HERE/'pilot/COMPLETE.json',dict(status='PASS',code_sha256=sha(__file__),results=results,sparse_geometry_vs_dense_regression='PASS',job_id=os.environ.get('SLURM_JOB_ID')))
    print(json.dumps(results,indent=2),flush=True)

def gene_layers(co,rows,out,cell_index):
    # Sparse global sequence coordinates keep the complete native sequence repertoire feasible.
    lookup={};col_chunks=[];data_chunks=[];row_chunks=[];receipts=[]
    for i,row in enumerate(rows):
        cache,r=source_vector(co,row,out);receipts.append(r)
        z=np.load(cache);keys=z['hashes'];v=z['values'];cols=np.empty(len(keys),np.int64)
        for j,k in enumerate(keys):
            key=bytes(k)
            if key not in lookup:lookup[key]=len(lookup)
            cols[j]=lookup[key]
        col_chunks.append(cols);row_chunks.append(np.full(len(keys),i,np.int32));data_chunks.append(v)
    cols=np.concatenate(col_chunks);rr=np.concatenate(row_chunks);values=np.concatenate(data_chunks)
    shape=(len(rows),len(lookup));nfeature=len(lookup)
    with (out/'gene_sequence_features.tsv.gz').open('wb') as _unused:pass
    import gzip
    with gzip.open(out/'gene_sequence_features.tsv.gz','wt') as f:
        f.write('feature_index\tencoded_protein_sha256\n')
        for key,j in lookup.items():f.write(str(j)+'\t'+key.decode()+'\n')
    del lookup,col_chunks,data_chunks,row_chunks
    # Average runs within specimens, normalize each specimen, then average specimens in a cell.
    specimens=list(dict.fromkeys(r['specimen'] for r in rows));sp_index={s:i for i,s in enumerate(specimens)}
    nrun=Counter(r['specimen'] for r in rows)
    run_to_sp=sparse.csr_matrix(([1/nrun[r['specimen']] for r in rows],([sp_index[r['specimen']] for r in rows],range(len(rows)))),shape=(len(specimens),len(rows)))
    sp_meta={r['specimen']:(r['participant'],r['region']) for r in rows}
    if co=='CAPSCAN':
        cell_map={k:i for i,k in enumerate(cell_index)};ns=Counter(sp_meta[s] for s in specimens)
        sp_to_cell=sparse.csr_matrix(([1/ns[sp_meta[s]] for s in specimens],([cell_map[sp_meta[s]] for s in specimens],range(len(specimens)))),shape=(len(cell_index),len(specimens)))
    else:
        assert specimens==list(cell_index)
        sp_to_cell=sparse.eye(len(specimens),format='csr')
    layers={}
    for j,key in enumerate(['DNA','RNA'] if co=='PREDICT' else ['DNA']):
        raw=sparse.csr_matrix((values[:,j],(rr,cols)),shape=shape);raw.eliminate_zeros()
        spec=run_to_sp@raw;total=np.asarray(spec.sum(1)).ravel();assert (total>0).all()
        spec=sparse.diags(1/total)@spec
        layers[key]=(sp_to_cell@spec).tocsr();layers[key].eliminate_zeros()
        assert np.allclose(np.asarray(layers[key].sum(1)).ravel(),1)
        sparse.save_npz(out/(key+'_native_cells.npz'),layers[key])
    return layers,dict(n_sources=len(rows),n_specimens=len(specimens),n_sequence_features=nfeature,input_receipts=[dict(id=r['source']['id'],tpm_sha256=r['tpm_sha256'],fasta_sha256=r['fasta_sha256'],genes=r['genes'],positive_sequences=r['distinct_positive_sequences'],sequence_resolved_TPM_fraction=r.get('sequence_resolved_TPM_fraction',[1.]),unresolved_gene_count=r.get('unresolved_gene_count',0)) for r in receipts])

def other_layers(co,cell_index,out):
    base=group_path(co)
    matrix=read(base/'simplified_study_group_matrix.tsv.gz',index_col=0).fillna(0).T
    assert matrix.index.is_unique and matrix.columns.is_unique and np.isfinite(matrix.values).all() and (matrix.values>=0).all()
    host=matrix.loc[:,matrix.columns.str.contains('_HUMAN',regex=False)]
    expected=read(CURRENT/(co+'_noclust')/'host_representative_intensity.tsv',index_col=0)
    assert list(host.columns)==list(expected.columns) and list(host.index)==list(expected.index)
    assert np.allclose(host.values,expected.values,rtol=1e-10,atol=1e-10)
    nonhost=matrix.loc[:,~matrix.columns.str.contains('_HUMAN',regex=False)]
    assert np.allclose(host.sum(1)+nonhost.sum(1),matrix.sum(1))
    met=read(OLD/(co+'_metabolite_intensity.tsv'),index_col=0)
    met=np.log10(met.where(met>0));frames={}
    if co=='PREDICT':
        ids=matrix.index.to_series().str.extract(r'(MicrobPred_\d+)',expand=False);assert ids.notna().all()
        frames['Non-host protein']=nonhost.groupby(ids).mean().loc[cell_index]
        frames['Host protein']=host.groupby(ids).mean().loc[cell_index]
        frames['Metabolites']=met.loc[cell_index]
        pd.DataFrame({'acquisition':matrix.index,'specimen':ids}).query('specimen in @cell_index').to_csv(out/'protein_acquisition_roster.tsv',sep='\t',index=False)
    else:
        design=read(OLD/'CAPSCAN_specimen_design.tsv',index_col=0);design.specimen=design.specimen.astype(str);met.index=met.index.astype(str)
        og=read(CURRENT/'CAPSCAN_noclust/microbial_OG.tsv',index_col=0).reindex(design.index)
        valid=set(design.loc[og.sum(1)>0,'specimen']);met=met.loc[met.index.isin(valid)]
        dspec=design.groupby('specimen').first().loc[met.index]
        for key,frame in [('Non-host protein',nonhost),('Host protein',host)]:
            byspec=relative(frame.loc[design.index]).groupby(design.specimen).median()
            ix=byspec.index.intersection(dspec.index)
            frames[key]=byspec.loc[ix].groupby([dspec.loc[ix,'subject_id'],dspec.loc[ix,'group']]).mean().loc[cell_index]
        frames['Metabolites']=met.groupby([dspec.subject_id,dspec.group]).mean().loc[cell_index]
        keep=dspec.apply(lambda r:(r.subject_id,r.group) in set(cell_index),axis=1)
        dspec.loc[keep,['subject_id','group']].to_csv(out/'protein_metabolite_specimen_roster.tsv',sep='\t')
    result={}
    for key,frame in frames.items():
        assert list(frame.index)==list(cell_index)
        if key!='Metabolites':
            assert np.isfinite(frame.values).all() and (frame.sum(1)>0).all()
            result[key]=sparse.csr_matrix(frame.values)
            pd.Series(frame.columns,name='representative').to_csv(out/(key.replace(' ','_')+'_features.tsv.gz'),sep='\t',index=False)
        else:result[key]=frame.values
    return result

def profile_distance(x,key,prevalence,capscan=False):
    if key=='Metabolites':
        full=x.shape[1];keep=np.isfinite(x).mean(0)>=.7;y=x[:,keep]
        y=np.where(np.isfinite(y),y,np.nanmedian(y,axis=0));sd=y.std(0);y=y[:,sd>1e-8]
        y=(y-y.mean(0))/sd[sd>1e-8]
        assert np.isfinite(y).all() and y.shape[1]>1
        g=y@y.T;qc=dict(input_features=full,retained_features=y.shape[1],minimum_finite_fraction=.7)
    else:
        full=x.shape[1];counts=np.asarray((x>0).sum(0)).ravel()
        keep=(counts>=math.ceil(x.shape[0]*prevalence-1e-10))&(counts>0)
        y=x[:,keep].tocsr();tot=np.asarray(y.sum(1)).ravel();alltot=np.asarray(x.sum(1)).ravel()
        assert y.shape[1]>1 and (tot>0).all(),(key,prevalence,y.shape,float(tot.min()))
        fractions=tot/alltot;y=sparse.diags(1/tot)@y;y.data=np.sqrt(y.data)
        g=(y@y.T).toarray()
        qc=dict(input_features=full,retained_features=y.shape[1],minimum_prevalence=prevalence,retained_signal_fraction=fractions.tolist(),median_retained_signal_fraction=float(np.median(fractions)))
    if capscan:
        n=len(g);assert n%5==0
        region=np.eye(5)[np.tile(np.arange(5),n//5)];r=np.eye(n)-region @ np.linalg.pinv(region)
        g=r@g@r.T
    d2=np.diag(g)[:,None]+np.diag(g)[None,:]-2*g
    assert d2.min()>-1e-8
    d=np.sqrt(np.maximum(d2,0));np.fill_diagonal(d,0)
    assert np.allclose(d,d.T,atol=1e-12)
    return d,qc

def bh(values):
    v=np.asarray(values);o=np.argsort(v);q=np.empty(len(v));q[o]=np.minimum(1,np.minimum.accumulate((v[o]*len(v)/np.arange(1,len(v)+1))[::-1])[::-1]);return q

def correlations(distances,co):
    n=len(next(iter(distances.values())));ix=np.triu_indices(n,1)
    rng=np.random.default_rng(SEED)
    if co=='CAPSCAN':
        blocks=np.arange(n).reshape(-1,5);perms=np.array([blocks[rng.permutation(len(blocks))].ravel() for _ in range(1999)])
    else:perms=np.array([rng.permutation(n) for _ in range(1999)])
    ranks={k:rankdata(v[ix]) for k,v in distances.items()};norm={}
    for k,r in ranks.items():
        centered=r-r.mean();scale=np.linalg.norm(centered);assert scale>1e-10,(co,k,'No distance variation')
        norm[k]=centered/scale
    rows=[]
    for a,b in combinations(distances,2):
        rho=float(norm[a]@norm[b]);checkrho=float(spearmanr(distances[a][ix],distances[b][ix]).statistic)
        assert abs(rho-checkrho)<1e-12
        bmat=squareform(norm[b]);extreme=0
        for chunk in np.array_split(perms,20):
            null=bmat[chunk[:,ix[0]],chunk[:,ix[1]]]@norm[a]
            extreme+=int((abs(null)>=abs(rho)-1e-12).sum())
        rows.append(dict(cohort=co,layer_a=a,layer_b=b,rho=rho,p=(1+extreme)/2000,permutations=1999,n_profiles=n,n_participants=n//5 if co=='CAPSCAN' else n))
    table=pd.DataFrame(rows);table['q']=bh(table.p);return table

def run(co):
    m=manifest();assert m['code_sha256']==sha(__file__)
    p=json.loads((HERE/'pilot/COMPLETE.json').read_text());assert p['status']=='PASS' and p['code_sha256']==sha(__file__)
    for r in m['shared']:check(r)
    out=HERE/'results'/co;out.mkdir(parents=True,exist_ok=True)
    assert not (out/'COMPLETE.json').exists(),'Completed run exists; do not overwrite'
    if co=='PREDICT':index=read(OLD/'PREDICT_fullstack_primary_samples.tsv',index_col=0).index
    else:index=pd.MultiIndex.from_product([m['capscan_donors'],REG],names=['participant','region'])
    # Reuse the checked pilot source-vector outputs by copying their small deterministic cache.
    import shutil
    (out/'sources').mkdir(exist_ok=True)
    for f in (HERE/'pilot'/co/'sources').iterdir():shutil.copy2(f,out/'sources'/f.name)
    layers,geneqc=gene_layers(co,m['items'][co],out,index)
    layers.update(other_layers(co,index,out));allrows=[];quality={};dist_archive={}
    for label,prevalence in [('primary_prevalence10',.1),('sensitivity_all_features',0.)]:
        ds={};quality[label]={}
        for k,x in layers.items():
            ds[k],qc=profile_distance(x,k,prevalence,capscan=co=='CAPSCAN');quality[label][k]=qc
            dist_archive[label+'__'+k]=ds[k]
        table=correlations(ds,co);table.insert(1,'analysis',label);allrows.append(table)
        print(co,label,table[['layer_a','layer_b','rho','q']].to_string(index=False),flush=True)
    result=pd.concat(allrows,ignore_index=True);result.to_csv(out/'correlations.tsv',sep='\t',index=False)
    np.savez_compressed(out/'distance_matrices.npz',**dist_archive)
    dump(out/'feature_quality.json',quality)
    hostrow=result[(result.analysis=='primary_prevalence10')&(result.layer_a=='Host protein')&(result.layer_b=='Metabolites')].iloc[0]
    if co=='PREDICT':
        ref=read(CURRENT/'crossomics/fullstack_profile_correlations.tsv');ref=ref[(ref.cohort=='PREDICT')&(ref.layer_a=='host_metaP')&(ref.layer_b=='metabolites')].iloc[0]
    else:
        ref=read(DNA_RESULT/'capscan_profile_correlations_with_DNA.tsv');ref=ref[(ref.n_donors==12)&(ref.layer_a=='host_metaP')&(ref.layer_b=='metabolites')].iloc[0]
    assert abs(hostrow.rho-ref.rho)<1e-12,(co,hostrow.rho,ref.rho)
    assert abs(hostrow.p-ref.p)<1e-12,(co,hostrow.p,ref.p)
    dense_control()
    dump(out/'COMPLETE.json',dict(status='PASS',cohort=co,code_sha256=sha(__file__),manifest_sha256=sha(HERE/'INPUT_MANIFEST.json'),job_id=os.environ.get('SLURM_JOB_ID'),n_profiles=len(index),layers=list(layers),gene_profile=geneqc,host_metabolite_control=dict(rho=float(hostrow.rho),saved_rho=float(ref.rho),p=float(hostrow.p),saved_p=float(ref.p)),checks=['TPM genes mapped within source assembly or explicitly exported as unresolved','resolved plus unresolved TPM conserves input abundance','native host intensities equal current published source','host/non-host partition conserves all protein intensity','host-metabolite rho/p reproduce frozen result','sparse Hellinger/regional geometry agrees with dense control','all Spearman rhos independently checked'],outputs={f.name:sha(f) for f in out.iterdir() if f.is_file()}))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('command',choices=['prepare','pilot','run']);p.add_argument('--cohort',choices=['CAPSCAN','PREDICT']);a=p.parse_args()
    if a.command=='prepare':prepare()
    elif a.command=='pilot':pilot()
    else:run(a.cohort)
