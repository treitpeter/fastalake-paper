"""Frozen native multi-omics method comparison and exploratory MOFA."""
import os
for name in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[name]='1'
import sys,json,hashlib,gzip,argparse,math,time,subprocess
from pathlib import Path
from itertools import combinations
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import rankdata,spearmanr,pearsonr
from scipy.spatial.distance import squareform,pdist
HERE=Path(__file__).resolve().parent;PARENT=HERE.parent
sys.path.insert(0,str(PARENT));import native_pyramids as native
SEEDS=[20260923,20260924,20260925]

def dump(p,x):native.dump(p,x)
def sha(p):return native.sha(p)
def read(p,**kw):return pd.read_csv(p,sep='\t',**kw)
def freeze():
    assert not (HERE/'MANIFEST.json').exists()
    upstream={co:sha(PARENT/'results'/co/'COMPLETE.json') for co in ['PREDICT','CAPSCAN']}
    import importlib.metadata as md
    dump(HERE/'MANIFEST.json',dict(code_sha256=sha(__file__),protocol_sha256=sha(HERE/'PROTOCOL.md'),upstream=upstream,
        original_code_sha256=sha(PARENT/'native_pyramids.py'),versions={k:md.version(k) for k in ['mofapy2','numpy','scipy','pandas','h5py']},
        wuyts_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=HERE/'external/multiomics_Wuyts_2022',text=True).strip()))
def gate(pilot=False):
    m=json.loads((HERE/'MANIFEST.json').read_text());assert m['code_sha256']==sha(__file__)
    assert m['protocol_sha256']==sha(HERE/'PROTOCOL.md') and m['original_code_sha256']==sha(PARENT/'native_pyramids.py')
    for co,digest in m['upstream'].items():assert sha(PARENT/'results'/co/'COMPLETE.json')==digest
    if pilot:
        p=json.loads((HERE/'PILOT_COMPLETE.json').read_text());assert p['status']=='PASS' and p['code_sha256']==sha(__file__)
    return m

def cohort_index(co):
    if co=='PREDICT':
        d=read(native.OLD/'PREDICT_fullstack_primary_samples.tsv',index_col=0)
        return d.index,pd.DataFrame({'profile':d.index,'participant':d.idpatient.values,'region':'stool'})
    idx=pd.MultiIndex.from_product([native.manifest()['capscan_donors'],native.REG],names=['participant','region'])
    return idx,pd.DataFrame({'profile':[str(p)+':'+r for p,r in idx],'participant':[p for p,r in idx],'region':[r for p,r in idx]})

def projection(n):
    a=np.eye(5)[np.tile(np.arange(5),n//5)];return np.eye(n)-a @ np.linalg.pinv(a)
def relative(x):
    tot=np.asarray(x.sum(1)).ravel();assert (tot>0).all();return sparse.diags(1/tot)@x

def rank_gram(x):
    """Exact centered unit-rank Gram matrix, including tied zeros, without densifying."""
    x=x.tocsr(copy=True);x.eliminate_zeros();n,m=x.shape
    counts=np.diff(x.indptr).astype(float);base=-counts/2
    for i in range(n):
        s=slice(x.indptr[i],x.indptr[i+1]);k=int(counts[i])
        x.data[s]=rankdata(x.data[s])+(m-k-1)/2
    g=(x@x.T).toarray()-m*np.outer(base,base)
    norms=np.sqrt(np.diag(g));assert (norms>1e-10).all()
    g=g/norms[:,None]/norms[None,:]
    assert np.allclose(np.diag(g),1) and np.max(abs(g))<1+1e-10
    return g

def gram_distance(g,adjusted=False,squared_half=False):
    if adjusted:r=projection(len(g));g=r@g@r.T
    d2=np.diag(g)[:,None]+np.diag(g)[None,:]-2*g
    assert d2.min()>-1e-8
    d=np.maximum(d2,0)/2 if squared_half else np.sqrt(np.maximum(d2,0))
    np.fill_diagonal(d,0);return d

def met_matrix(x,recipe=False):
    keep=np.isfinite(x).mean(0)>=.7;y=x[:,keep].copy()
    if recipe:
        y=10**y;y[~np.isfinite(y)]=0;y/=y.sum(1,keepdims=True)
        for j in range(y.shape[1]):y[y[:,j]==0,j]=y[y[:,j]>0,j].min()
    else:
        y=np.where(np.isfinite(y),y,np.nanmedian(y,axis=0));sd=y.std(0);y=y[:,sd>1e-8];y=(y-y.mean(0))/sd[sd>1e-8]
    return y

def gram_for(x,key,variant,metric):
    if key=='Metabolites':
        y=met_matrix(x,variant=='wuyts_recipe')
        if metric=='rank':
            y=np.apply_along_axis(rankdata,1,y);y-=y.mean(1,keepdims=True);y/=np.linalg.norm(y,axis=1,keepdims=True)
        return y@y.T,dict(retained_features=y.shape[1])
    x=relative(x).tocsr();n,m=x.shape
    if variant=='prevalence10':keep=np.asarray((x>0).sum(0)).ravel()>=math.ceil(n*.1)
    elif variant=='wuyts_recipe':keep=np.asarray(x.mean(0)).ravel()>(1e-7 if key in ['DNA','RNA'] else 1e-5)
    else:keep=np.asarray((x>0).sum(0)).ravel()>0
    y=x[:,keep].tocsr();signal=np.asarray(y.sum(1)).ravel();assert (signal>0).all()
    qc=dict(input_features=m,retained_features=int(keep.sum()),retained_signal_fraction=signal.tolist())
    if metric=='rank':g=rank_gram(y)
    else:
        y=relative(y).tocsr();y.data=np.sqrt(y.data);g=(y@y.T).toarray()
    return g,qc

def rank_checks():
    x=np.random.default_rng(18).integers(0,9,(15,80)).astype(float);x[:,0:45]=0
    x[:,51:56]=3
    g=rank_gram(sparse.csr_matrix(x));assert np.allclose(g,spearmanr(x,axis=1).statistic,atol=1e-12)
    expected=1-spearmanr(x,axis=1).statistic
    assert np.allclose(gram_distance(g,squared_half=True),expected,atol=1e-12)
    y=np.apply_along_axis(rankdata,1,x);y-=y.mean(1,keepdims=True);y/=np.linalg.norm(y,axis=1,keepdims=True)
    y=projection(len(y))@y
    assert np.allclose(gram_distance(g,True,True),squareform(pdist(y,'sqeuclidean'))/2,atol=1e-12)

def load_layers(co,out):
    idx,meta=cohort_index(co);out.mkdir(parents=True,exist_ok=True)
    layers={k:sparse.load_npz(PARENT/'results'/co/(k+'_native_cells.npz')) for k in (['DNA','RNA'] if co=='PREDICT' else ['DNA'])}
    layers.update(native.other_layers(co,idx,out));return layers,meta

def select_mofa(x,key,adjusted):
    n=x.shape[0]
    if key in ['DNA','RNA']:
        counts=np.asarray((x>0).sum(0)).ravel();eligible=np.flatnonzero(counts>=math.ceil(.2*n))
        z=x[:,eligible].tocsr(copy=True);z.data=np.log2(1+1e6*z.data)
        mean=np.asarray(z.mean(0)).ravel();var=np.asarray(z.power(2).mean(0)).ravel()-mean**2
        if adjusted:
            means=np.vstack([np.asarray(z[r::5].mean(0)).ravel() for r in range(5)])
            var=np.asarray(z.power(2).mean(0)).ravel()-np.mean(means**2,axis=0)
        order=np.argsort(-var,kind='stable');chosen=order[var[order]>1e-10][:1000]
        y=z[:,chosen].toarray();indices=eligible[chosen]
        coverage=np.asarray(x[:,indices].sum(1)).ravel().tolist()
    else:
        a=x if key=='Metabolites' else relative(x).toarray()
        if key!='Metabolites':a=np.where(a>0,np.log2(np.maximum(a,1e-300)*1e6),np.nan)
        frac=.7 if key=='Metabolites' else .5
        eligible=np.flatnonzero(np.isfinite(a).mean(0)>=frac);z=a[:,eligible].copy()
        if adjusted:
            for r in range(5):z[r::5]-=np.nanmean(z[r::5],axis=0)
        var=np.nanvar(z,axis=0);order=np.argsort(-var,kind='stable');chosen=order[var[order]>1e-10][:1000]
        y=z[:,chosen];indices=eligible[chosen]
        coverage=None if key=='Metabolites' else np.asarray(relative(x)[:,indices].sum(1)).ravel().tolist()
    if adjusted and key in ['DNA','RNA']:y=projection(n)@y
    y=(y-np.nanmean(y,axis=0))/np.nanstd(y,axis=0)
    assert y.shape[1]>10 and (np.isfinite(y).sum(1)>0).all()
    return y,indices,dict(input_features=x.shape[1],eligible_features=len(eligible),selected_features=y.shape[1],
        selected_signal_fraction=coverage,missing_fraction=float(np.isnan(y).mean()))

def feature_names(co,wanted):
    found={}
    with gzip.open(PARENT/'results'/co/'gene_sequence_features.tsv.gz','rt') as f:
        next(f)
        for line in f:
            i,name=line.rstrip().split('\t');i=int(i)
            if i in wanted:found[i]='sha256:'+name
    assert len(found)==len(wanted);return found

def fit_model(co,adjustment,budget,seed,pilot=False):
    from mofapy2.run.entry_point import entry_point
    folder=HERE/'prepared'/co;data=np.load(folder/('mofa_'+adjustment+'.npz'))
    info=json.loads((folder/('mofa_'+adjustment+'.json')).read_text());views=list(info['views']);meta=read(folder/'profiles.tsv')
    mats=[];names=[]
    for view in views:
        y=data[view][:,:budget].copy();fn=info['views'][view]['features'][:budget]
        if pilot:
            y=y[:20,:100];fn=fn[:100];ok=(np.isfinite(y).mean(0)>.1)&(np.nanvar(y,axis=0)>1e-8);y=y[:,ok];fn=np.array(fn)[ok].tolist()
        mats.append([y]);names.append([view+'::'+str(v) for v in fn])
    ids=meta.profile.astype(str).tolist()[:len(mats[0][0])]
    out=HERE/('pilot_models' if pilot else 'mofa')/co/adjustment/('features'+str(budget))/('seed'+str(seed));out.mkdir(parents=True,exist_ok=True)
    assert not (out/'COMPLETE.json').exists(),'Completed fit must not be overwritten'
    k=3 if pilot else (10 if co=='PREDICT' else 5);iterations=50 if pilot else 2000
    ent=entry_point();ent.set_data_options(scale_views=True,center_groups=True)
    ent.set_data_matrix(mats,likelihoods=['gaussian']*len(views),views_names=views,groups_names=[co],samples_names=[ids],features_names=names)
    ent.set_model_options(factors=k,spikeslab_weights=True,ard_weights=True)
    ent.set_train_options(iter=iterations,convergence_mode='slow',seed=seed,quiet=True,freqELBO=5,startELBO=1)
    ent.build();ent.run();assert ent.model.trained
    ent.save(str(out/'model.hdf5'))
    exp=ent.model.getExpectations(only_first_moments=True);z=exp['Z'];w=exp['W'];assert np.isfinite(z).all()
    factors=['Factor'+str(i+1) for i in range(z.shape[1])]
    pd.DataFrame(z,index=ids,columns=factors).to_csv(out/'factor_scores.tsv',sep='\t',index_label='profile')
    r2=ent.model.calculate_variance_explained()[0];total=ent.model.calculate_variance_explained(total=True)[0]
    pd.DataFrame(r2,index=views,columns=factors).to_csv(out/'variance_explained.tsv',sep='\t',index_label='view')
    for view,features,weights in zip(views,names,w):
        pd.DataFrame(weights,index=features,columns=factors).to_csv(out/(view.replace(' ','_')+'_loadings.tsv.gz'),sep='\t',index_label='feature')
    elbo=np.array(ent.model.getTrainingStats()['elbo']);pd.Series(elbo,name='ELBO').to_csv(out/'training_ELBO.tsv',sep='\t',index_label='iteration')
    finite=elbo[np.isfinite(elbo)];assert len(finite)>1 and np.isfinite(total).all()
    dump(out/'COMPLETE.json',dict(status='PASS',cohort=co,adjustment=adjustment,budget=budget,seed=seed,pilot=pilot,
        factors=z.shape[1],n_profiles=len(ids),elbo_final=float(finite[-1]),converged_before_limit=len(elbo)<iterations,
        total_variance_explained=dict(zip(views,map(float,total))),code_sha256=sha(__file__),job_id=os.environ.get('SLURM_JOB_ID'),
        outputs={p.name:sha(p) for p in out.iterdir() if p.is_file()}))
    return out

def prepare_pilot():
    gate();rank_checks()
    for co in ['PREDICT','CAPSCAN']:
        out=HERE/'prepared'/co;layers,meta=load_layers(co,out);meta.to_csv(out/'profiles.tsv',sep='\t',index=False)
        for key,x in layers.items():
            if key not in ['DNA','RNA']:
                if key=='Metabolites':np.save(out/'Metabolites_log_profiles.npy',x)
                else:sparse.save_npz(out/(key.replace(' ','_')+'.npz'),x)
        indexsets={};arrays={};reports={};wanted=set()
        for adj in ['unadjusted']+(['region_adjusted'] if co=='CAPSCAN' else []):
            arrays[adj]={};indexsets[adj]={};reports[adj]={}
            for key,x in layers.items():
                y,indices,q=select_mofa(x,key,adj=='region_adjusted');arrays[adj][key]=y;indexsets[adj][key]=indices;reports[adj][key]=q
                if key in ['DNA','RNA']:wanted.update(map(int,indices))
        genes=feature_names(co,wanted)
        for adj in arrays:
            np.savez_compressed(out/('mofa_'+adj+'.npz'),**arrays[adj]);views={}
            for key,indices in indexsets[adj].items():
                if key in ['DNA','RNA']:names=[genes[int(i)] for i in indices]
                elif key=='Metabolites':names=read(native.OLD/(co+'_metabolite_intensity.tsv'),nrows=0).columns[1:][indices].tolist()
                else:names=read(out/(key.replace(' ','_')+'_features.tsv.gz')).representative.iloc[indices].tolist()
                views[key]=dict(features=names,source_column_indices=indices.tolist(),quality=reports[adj][key])
            dump(out/('mofa_'+adj+'.json'),dict(views=views,cohort=co,adjustment=adj))
        old=np.load(PARENT/'results'/co/'distance_matrices.npz')
        for key,x in layers.items():
            for variant,oldlabel in [('all_observed','sensitivity_all_features'),('prevalence10','primary_prevalence10')]:
                g,_=gram_for(x,key,variant,'hellinger');d=gram_distance(g,co=='CAPSCAN')
                assert np.allclose(d,old[oldlabel+'__'+key],atol=1e-10),(co,key,variant)
            if key in ['DNA','RNA']:
                small=x[:15,:2000].toarray();ok=np.std(small,axis=1)>0;small=small[ok]
                assert len(small)>2 and np.allclose(rank_gram(sparse.csr_matrix(small)),spearmanr(small,axis=1).statistic,atol=1e-12)
        del layers
        fit_model(co,'unadjusted',100,20260923,pilot=True)
        dump(out/'COMPLETE.json',dict(status='PASS',code_sha256=sha(__file__),checks=['native inputs frozen','original Hellinger distances reproduced','actual sparse rank geometry versus scipy Spearman'],outputs={p.name:sha(p) for p in out.iterdir() if p.is_file()}))
        print(co,'preparation and real-data pilot PASS',flush=True)
    dump(HERE/'PILOT_COMPLETE.json',dict(status='PASS',code_sha256=sha(__file__),manifest_sha256=sha(HERE/'MANIFEST.json'),job_id=os.environ.get('SLURM_JOB_ID')))

def load_prepared(co):
    folder=HERE/'prepared'/co;p=json.loads((folder/'COMPLETE.json').read_text());assert p['status']=='PASS' and p['code_sha256']==sha(__file__)
    for name,digest in p['outputs'].items():assert sha(folder/name)==digest
    layers={k:sparse.load_npz(PARENT/'results'/co/(k+'_native_cells.npz')) for k in (['DNA','RNA'] if co=='PREDICT' else ['DNA'])}
    for k in ['Non-host protein','Host protein']:layers[k]=sparse.load_npz(folder/(k.replace(' ','_')+'.npz'))
    layers['Metabolites']=np.load(folder/'Metabolites_log_profiles.npy');return layers

def correlate(ds,co,outer):
    n=len(next(iter(ds.values())));ix=np.triu_indices(n,1);rng=np.random.default_rng(native.SEED)
    blocks=np.arange(n).reshape(-1,5) if co=='CAPSCAN' else None
    perms=np.array([blocks[rng.permutation(len(blocks))].ravel() if blocks is not None else rng.permutation(n) for _ in range(1999)])
    norm={}
    for key,d in ds.items():
        v=rankdata(d[ix]) if outer=='spearman' else d[ix].copy();v-=v.mean();norm[key]=v/np.linalg.norm(v)
    rows=[]
    for a,b in combinations(ds,2):
        r=float(norm[a]@norm[b]);ind=(spearmanr if outer=='spearman' else pearsonr)(ds[a][ix],ds[b][ix]).statistic;assert abs(ind-r)<1e-12
        mat=squareform(norm[b]);null=mat[perms[:,ix[0]],perms[:,ix[1]]]@norm[a]
        rows.append(dict(cohort=co,layer_a=a,layer_b=b,coefficient=r,p=(1+int((abs(null)>=abs(r)-1e-12).sum()))/2000,
            p_positive_tail=(1+int((null>=r-1e-12).sum()))/2000,n_profiles=n,n_participants=n//5 if co=='CAPSCAN' else n))
    t=pd.DataFrame(rows);t['q']=native.bh(t.p);return t

def distances(co):
    gate(True);layers=load_prepared(co);out=HERE/'distances'/co;out.mkdir(parents=True,exist_ok=True)
    assert not (out/'COMPLETE.json').exists();allrows=[];archive={};quality={}
    for variant in ['all_observed','prevalence10','wuyts_recipe']:
        for metric in ['hellinger','rank']:
            if variant=='wuyts_recipe' and metric=='hellinger':continue
            grams={};quality[variant+'__'+metric]={}
            for key,x in layers.items():grams[key],quality[variant+'__'+metric][key]=gram_for(x,key,variant,metric)
            for adj in ['unadjusted']+(['region_adjusted'] if co=='CAPSCAN' else []):
                ds={key:gram_distance(g,adj=='region_adjusted',metric=='rank') for key,g in grams.items()}
                for key,d in ds.items():archive['__'.join([variant,metric,adj,key])]=d
                for outer in ['pearson','spearman']:
                    t=correlate(ds,co,outer);t.insert(1,'features',variant);t.insert(2,'distance',metric);t.insert(3,'adjustment',adj);t.insert(4,'outer_correlation',outer);allrows.append(t)
    pd.concat(allrows,ignore_index=True).to_csv(out/'comparisons.tsv',sep='\t',index=False)
    np.savez_compressed(out/'distance_matrices.npz',**archive);dump(out/'feature_quality.json',quality)
    dump(out/'COMPLETE.json',dict(status='PASS',code_sha256=sha(__file__),job_id=os.environ.get('SLURM_JOB_ID'),outputs={p.name:sha(p) for p in out.iterdir() if p.is_file()}))

def full_mofa(co,adj):
    gate(True)
    for budget in [1000,500]:
        for seed in SEEDS:fit_model(co,adj,budget,seed)
    dump(HERE/'mofa'/co/adj/'COMPLETE.json',dict(status='PASS',cohort=co,adjustment=adj,code_sha256=sha(__file__),job_id=os.environ.get('SLURM_JOB_ID'),fits=6))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('command',choices=['freeze','pilot','distance','mofa']);p.add_argument('--cohort',choices=['PREDICT','CAPSCAN']);p.add_argument('--adjustment',default='unadjusted',choices=['unadjusted','region_adjusted']);a=p.parse_args()
    if a.command=='freeze':freeze()
    elif a.command=='pilot':prepare_pilot()
    elif a.command=='distance':distances(a.cohort)
    else:full_mofa(a.cohort,a.adjustment)
