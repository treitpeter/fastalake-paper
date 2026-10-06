#!/usr/bin/env python3
"""Reconstruct identified, quantified and transferred HSTOOL evidence separately."""
import csv,gzip,json,math,os,sys
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np
import pandas as pd
csv.field_size_limit(1<<30)
RUN=Path(sys.argv[1]).resolve();sys.path.insert(0,str(RUN/'frozen'))
from lib import dump,meta,sha,verify,table
from peptide import canonical
from parameters import parameters
from origin import classify_members,is_classifiable_accession
from digest import fasta

def key(p,c):return p.replace('I','L')+'|charge='+str(c)
def protkey(p):return ';'.join(sorted(set(p.split(';'))))

def evidence(path,names):
    ids={n:set() for n in names};features={n:set() for n in names};pq={};peps_by_group=defaultdict(set);qc=Counter()
    # Capture every native protein-set q value, not only peptide-q passing PSMs.
    with open(path) as f:
        r=csv.DictReader(f,delimiter='\t');assert {'rank','label','peptide_q','protein_q','peptide','proteins','filename','charge'}<=set(r.fieldnames)
        for row in r:
            qc['psm_rows']+=1
            if row['rank']!='1' or row['label']!='1':continue
            name=Path(row['filename']).name;assert name in names
            q=float(row['peptide_q']);p_q=float(row['protein_q']);assert math.isfinite(q) and math.isfinite(p_q) and 0<=q<=1 and 0<=p_q<=1
            group=protkey(row['proteins'])
            if group in pq:assert pq[group]==p_q,('Native protein q varies within a search',group,pq[group],p_q)
            pq[group]=p_q
            if q>.01:continue
            p=canonical(row['peptide']);assert 5<=len(p)<=50
            ids[name].add(p);features[name].add(row['peptide']);peps_by_group[group].add(p)
    return ids,features,pq,peps_by_group,dict(qc)

def readsets(spec):
    if 'baseline_reconstruction' in spec:
        b=spec['baseline_reconstruction']
        for k in ['receipt','new','gained','lost']:verify(b[k])
        a=set(Path(b['new']['path']).read_text().splitlines());gain=set(Path(b['gained']['path']).read_text().splitlines());lost=set(Path(b['lost']['path']).read_text().splitlines())
        assert gain<=a and not a&lost;old=(a-gain)|lost
        assert len(old)==b['expected_old'] and len(a)==b['expected_new'];return old
    verify(spec);return set(Path(spec['path']).read_text().splitlines())

def run(index):
    m=json.loads((RUN/'inputs/manifest_revision2.json').read_text());pilot=index<5;t=m['arms'][index%5]
    for x in m['frozen_code']+[m['protocol']]:verify(x)
    out=RUN/'results'/('pilot_'+t['name'] if pilot else t['name']);out.mkdir(parents=True,exist_ok=True)
    assert not (out/'COMPLETE.json').exists()
    samples=m['samples'];selected=samples[:1] if pilot else samples;selected_names=[s+'.mzML' for s in selected]
    mapping=pd.read_csv(m['mapping']['path'],sep='\t',dtype=str).set_index('canonical_peptide').study_group.to_dict();verify(m['mapping'])
    aliases=defaultdict(set)
    with gzip.open(m['union_provenance']['path'],'rt') as fh:
        for r in csv.DictReader(fh,delimiter='\t'):aliases[r['union_id']].add(r['original_accession'])
    verify(m['union_provenance'])
    class_table=pd.read_csv(m['local_classification']['path'],sep='\t',usecols=['protein_group','classification'])
    classes={protkey(p):c for p,c in zip(class_table.protein_group,class_table.classification)}
    verify(m['local_classification'])
    origin_mats={c:{} for c in ['Human','Microbial','Mixed','Unknown']}
    records=t['records'][:1] if pilot and len(t['records'])>1 else t['records']
    allmat=defaultdict(dict);idmat=defaultdict(dict);strictmat=defaultdict(dict);q01mat=defaultdict(dict);q05mat=defaultdict(dict)
    ids_by_sample={};native_frames=[];native_labels={};coverage=[];identity=[];inputs=[];quality=[];transfer=[]
    native_meta={};reference_parameters=None;source_mappings={}
    for rec in records:
        for x in rec['inputs']:verify(x)
        config=json.loads(Path(rec['config']['path']).read_text());all_names=[Path(p).name for p in config['mzml_paths']]
        assert len(all_names)==len(set(all_names));assert set(all_names)<=set(s+'.mzML' for s in samples)
        assert all(Path(p).resolve()==Path(m['raws'][Path(p).stem]['path']).resolve() for p in config['mzml_paths'])
        pp=parameters(config)
        assert pp==m['search_parameters'],'Search identification settings differ from baseline'
        assert config['quant']['lfq'] and config['quant']['lfq_settings']['combine_charge_states']
        actual_names=[n for n in all_names if n in selected_names];assert actual_names
        ids,modified,pq,native_peps,qc=evidence(rec['psm']['path'],set(all_names))
        cols=['peptide','charge','proteins','q_value','score','spectral_angle']+all_names
        f=pd.read_csv(rec['lfq']['path'],sep='\t',usecols=cols)
        assert f.charge.eq(-1).all(),'Unexpected noncombined-charge LFQ feature'
        f['key']=[key(p,c) for p,c in zip(f.peptide,f.charge)];f['canonical']=f.peptide.map(canonical);f['native_group']=f.proteins.map(protkey)
        fa_audit={};fa={h.split()[0]:seq for h,seq in fasta(rec['fasta']['path'],fa_audit)}
        assert fa_audit['sha256']==rec['fasta']['sha256']
        rowclasses=[];group_classes={}
        for peptide,group in set(zip(f['canonical'],f['native_group'])):
            members=group.split(';')
            assert all(pid in fa and peptide in fa[pid] for pid in members),('LFQ peptide/parent containment',peptide,group)
            expanded=set().union(*(aliases[pid] if pid.startswith('FLU_') else {pid} for pid in members))
            assert expanded
            group_classes[group]=classes.get(protkey(';'.join(expanded)),classify_members(expanded))
        rowclasses=f.native_group.map(group_classes).to_numpy()
        raw=f[actual_names].to_numpy(float);positive=np.isfinite(raw)&(raw>0);v=np.where(positive,raw,0.)
        q=f.q_value.to_numpy(float);assert np.isfinite(q).all() and (q>=0).all() and (q<=1).all()
        gate=np.stack([f.canonical.isin(ids[n]).to_numpy() for n in actual_names],axis=1)
        # Modification-specific same-acquisition evidence is retained as a stricter sensitivity.
        strict=np.stack([f.peptide.isin(modified[n]).to_numpy() for n in actual_names],axis=1)
        assert not (strict & ~gate).any()
        for mat,mask in [(allmat,np.ones_like(gate)),(idmat,gate),(strictmat,strict),(q01mat,np.repeat((q<=.01)[:,None],len(actual_names),axis=1)),(q05mat,np.repeat((q<=.05)[:,None],len(actual_names),axis=1))]:
            frame=pd.DataFrame(np.where(mask,v,0.),index=f['key'],columns=[Path(n).stem for n in actual_names]).groupby(level=0).sum()
            for sample in frame:mat[sample]=frame[sample]
        for cls,store in origin_mats.items():
            vals=np.where(gate & (rowclasses==cls)[:,None],v,0.)
            origin_frame=pd.DataFrame(vals,index=f['key'],columns=[Path(n).stem for n in actual_names]).groupby(level=0).sum()
            for sample in origin_frame:store[sample]=origin_frame[sample]
        for j,n in enumerate(actual_names):
            sample=Path(n).stem;ids_by_sample[sample]=ids[n]
            if sample in rec['expected_sets']:assert ids[n]==readsets(rec['expected_sets'][sample]),('Independent PSM recount differs',t['name'],sample)
            coverage.append(dict(sample=sample,arm=t['name'],identified_canonical_peptides=len(ids[n]),lfq_source_rows=len(f),
                positive_source_cells=int(positive[:,j].sum()),same_file_id_source_cells=int((positive[:,j]&gate[:,j]).sum()),
                same_file_modified_id_source_cells=int((positive[:,j]&strict[:,j]).sum()),
                positive_without_same_file_id=int((positive[:,j]&~gate[:,j]).sum()),
                pooled_q01_positive=int((positive[:,j]&(q<=.01)).sum()),pooled_q05_positive=int((positive[:,j]&(q<=.05)).sum()),
                negative_cells=int((raw[:,j]<0).sum()),nonfinite_cells=int((~np.isfinite(raw[:,j])).sum())))
            for qlabel,cut in [('no_lfq_q_filter',1.),('pooled_q01',.01),('pooled_q05',.05)]:
                for state,mask in [('same_acquisition_id',gate[:,j]),('no_same_acquisition_id',~gate[:,j])]:
                    ok=positive[:,j] & (q<=cut) & mask
                    transfer.append(dict(arm=t['name'],sample=sample,lfq_gate=qlabel,state=state,source_cells=int(ok.sum()),intensity=float(v[ok,j].sum()),confidence_scope='pooled_feature_not_cell' if t['corrected'] else 'old_exported_lfq_q_unreliable'))
        # Native Sage output retains parent-protein sets. No FastaLake representative
        # cover is used here. Quantification is done separately by directLFQ.
        if t['corrected']:
            native_q=f.native_group.map(pq);assert native_q.notna().all(),'LFQ native group has no Sage picked-protein score'
            valid=native_q.le(.01).to_numpy() & f.canonical.isin(set().union(*ids.values())).to_numpy()
            for mode,keep in [('same_file_id',valid[:,None]&gate),('pooled_q01',valid[:,None]&(q<=.01)[:,None]),('pooled_q05',valid[:,None]&(q<=.05)[:,None])]:
                nv=np.where(keep,v,0.);native=f[['native_group','peptide','charge']].copy()
                native['ion']=native.peptide+'|charge='+native.charge.astype(str)
                for j,n in enumerate(actual_names):native[Path(n).stem]=nv[:,j]
                native=native.groupby(['native_group','ion'],sort=True)[[Path(n).stem for n in actual_names]].sum()
                native=native[native.gt(0).any(axis=1)]
                assert native.index.get_level_values('ion').is_unique,'A native ion maps to more than one native protein set'
                np.testing.assert_allclose(native.sum().to_numpy(),nv.sum(axis=0),rtol=1e-12)
                native.index.names=['protein','ion'];dest=out/(mode+'.native_sage_sets.aq_reformat.tsv');native.reset_index().to_csv(dest,sep='\t',index=False)
                native_meta[mode]=dict(groups=native.index.get_level_values('protein').nunique(),ions=len(native),positive_cells=int(native.gt(0).sum().sum()),input=meta(dest))
            table(out/'native_protein_evidence.tsv',[dict(protein_set=g,sage_protein_q=pq[g],accepted_canonical_peptides=len(native_peps[g])) for g in sorted(set(f.native_group))])
        inputs.extend(rec['inputs']);quality.append(dict(record=rec['id'],**qc,lfq_rows=len(f)))
    assert set(idmat)==set(selected)
    outputs=[]
    for name,mat in [('all_positive',allmat),('same_file_id',idmat),('same_file_modified_id',strictmat),('pooled_q01',q01mat),('pooled_q05',q05mat)]:
        frame=pd.DataFrame(mat).reindex(columns=selected).fillna(0.);frame=frame[frame.gt(0).any(axis=1)]
        assert frame.index.is_unique and np.isfinite(frame.to_numpy()).all() and (frame.to_numpy()>=0).all()
        np.savez_compressed(out/(name+'.features.npz'),keys=frame.index.to_numpy(str),samples=np.array(selected),values=frame.to_numpy())
        # Projection onto the PRE-EXISTING dictionary isolates changes in measured
        # features. It is explicitly restricted; unmatched peptides are counted.
        groups=pd.Series([mapping.get(canonical(k.rsplit('|charge=',1)[0])) for k in frame.index],index=frame.index)
        known=groups.notna();g=frame.loc[known].groupby(groups[known]).sum()
        np.testing.assert_allclose(g.sum().to_numpy()+frame.loc[~known].sum().to_numpy(),frame.sum().to_numpy(),rtol=1e-12)
        np.savez_compressed(out/(name+'.fixed_groups.npz'),keys=g.index.to_numpy(str),samples=np.array(selected),values=g.to_numpy())
        source_mappings[name]=dict(features=len(frame),mapped_features=int(known.sum()),unmapped_features=int((~known).sum()),
            mapped_positive_cells=int(frame.loc[known].gt(0).sum().sum()),unmapped_positive_cells=int(frame.loc[~known].gt(0).sum().sum()),groups=len(g))
    origin_group_arrays={}
    for cls,mat in origin_mats.items():
        frame=pd.DataFrame(mat).reindex(columns=selected).fillna(0.);groups=pd.Series([mapping.get(canonical(k.rsplit('|charge=',1)[0])) for k in frame.index],index=frame.index)
        known=groups.notna();g=frame.loc[known].groupby(groups[known]).sum();g=g[g.gt(0).any(axis=1)]
        origin_group_arrays[cls]=g
        np.savez_compressed(out/('same_file_id.fixed_groups.'+cls+'.npz'),keys=g.index.to_numpy(str),samples=np.array(selected),values=g.to_numpy())
    if t['name']=='original' and not pilot:
        z=np.load(m['baseline_local_groups']['path']);verify(m['baseline_local_groups']);assert set(z['samples'])==set(selected)
        for cls in ['Human','Microbial']:
            actual=origin_group_arrays[cls].reindex(index=z['groups'],columns=z['samples']).fillna(0.).to_numpy()
            np.testing.assert_allclose(actual,z[cls],rtol=1e-10,atol=1e-6,err_msg='Original local-classified group quantities do not reproduce the paper baseline')
    with gzip.open(out/'accepted_identities.tsv.gz','wt') as f:
        w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['sample','peptide'])
        for sample,ps in ids_by_sample.items():w.writerows((sample,p) for p in sorted(ps))
    table(out/'coverage.tsv',coverage);table(out/'transfer_coverage.tsv',transfer)
    dump(out/'COMPLETE.json',dict(status='PASS',arm=t['name'],pilot=pilot,acquisitions=len(selected),corrected_lfq_confidence=t['corrected'],
        quality=quality,dictionary_projection=source_mappings,native=native_meta,inputs=inputs,
        exclusions_explicit=True,job_id=os.environ.get('SLURM_JOB_ID'),
        outputs=[meta(p) for p in out.iterdir() if p.is_file() and p.name!='COMPLETE.json']))
    print('HSTOOL QUANT EXTRACT PASS',t['name'],len(selected),source_mappings['same_file_id'],flush=True)

if __name__=='__main__':run(int(sys.argv[2]))
