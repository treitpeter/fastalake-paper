#!/usr/bin/env python3
"""Predeclared participant-level terminal and selection-threshold comparisons."""
import hashlib,itertools,json,math,os,sys
from pathlib import Path
import numpy as np
import pandas as pd
RUN=Path(sys.argv[1]).resolve();sys.path.insert(0,str(RUN/'frozen'))
from lib import dump,meta,verify
m=json.loads((RUN/'inputs/manifest.json').read_text());index=int(sys.argv[2]);pilot=index==0
for x in m['code']+m['inputs']:verify(x)
mode=['pilot','denovo','CAPSCAN','PREDICT'][index];out=RUN/'results'/mode;out.mkdir(exist_ok=True);assert not (out/'COMPLETE.json').exists()
REGIONS=['duodenum','jejunum','ileum','colon','stool'];AA=set('ACDEFGHIKLMNPQRSTVWY');NONKR=AA-{'K','R'}
receipts=[];rows=[];raw_design=[];production_checks=[]
tasks=m['tasks'][:3] if pilot else [t for t in m['tasks'] if mode=='denovo' or t['cohort']==mode]
keys=['cohort','route','min_length','unit','origin','terminus','metric']
for t in tasks:
    d=dict(cohort=t['cohort'],sample=t['sample'],**t['design']);raw_design.append(d)
    a=Path(t['denovo_receipt']);c=json.loads(a.read_text());assert c['status']=='PASS';receipts.append(meta(a))
    output={Path(x['path']).name:x for x in c['outputs']};verify(output['terminal_counts.tsv']);f=pd.read_csv(output['terminal_counts.tsv']['path'],sep='\t')
    if mode in ['CAPSCAN','PREDICT'] or (pilot and Path(t['routes_receipt']).exists()):
        p=Path(t['routes_receipt']);q=json.loads(p.read_text());assert q['status']=='PASS';receipts.append(meta(p))
        outputs={Path(x['path']).name:x for x in q['outputs']};verify(outputs['terminal_counts.tsv']);g=pd.read_csv(outputs['terminal_counts.tsv']['path'],sep='\t')
        kk=['route','min_length','origin_policy','unit','origin','terminus','residue']
        left=f[f.route.eq('denovo')].sort_values(kk).reset_index(drop=True)
        right=g[g.route.eq('denovo')].sort_values(kk).reset_index(drop=True)
        pd.testing.assert_frame_equal(left,right,check_dtype=False,check_exact=False,rtol=1e-9,atol=1e-10)
        production_checks.append(t['sample']);f=pd.concat([f,g[~g.route.eq('denovo')]],ignore_index=True)
    f=f[f.origin_policy.eq('conservative')]
    for labels,gg in f.groupby(['route','min_length','unit','origin','terminus'],sort=False):
        route,lo,unit,origin,terminus=labels;assert gg.denominator.nunique()==1
        counts=gg.set_index('residue')['count'];den=float(gg.denominator.iloc[0])
        residue_rows=~counts.index.isin(['K','R','PROTEIN_END','UNMAPPED']);nk=float(counts[residue_rows].sum());h=float(counts.get('H',0))
        values={'nonKR_fraction_all':nk/den if den else np.nan,'H_fraction_all':h/den if den else np.nan,'H_fraction_nonKR':h/nk if nk else np.nan}
        if unit=='assigned_bare_sequence':
            for aa in sorted(NONKR):values['residue_'+aa+'_fraction_nonKR']=float(counts.get(aa,0))/nk if nk else np.nan
            values['other_residue_fraction_nonKR']=float(counts[residue_rows & ~counts.index.isin(NONKR)].sum())/nk if nk else np.nan
        for metric,value in values.items():rows.append(dict(d,route=route,min_length=lo,unit=unit,origin=origin,terminus=terminus,metric=metric,value=value,denominator=den,nonKR_denominator=nk))
    print('TERMINAL STATS INPUT',t['cohort'],t['sample'],flush=True)
f=pd.DataFrame(rows);f.to_csv(out/'per_acquisition_endpoints.tsv.gz',sep='\t',index=False)
design=pd.DataFrame(raw_design);design.to_csv(out/'acquisition_design.tsv',sep='\t',index=False)
f['participant']=f.participant.astype(str);valid=f[f.participant.ne('')]
sp=valid.groupby(keys+['participant','region','specimen_id'],dropna=False).agg(value=('value','mean'),acquisitions=('sample','nunique'),defined_acquisitions=('value','count')).reset_index()
pa=sp.groupby(keys+['participant','region'],dropna=False).agg(value=('value','mean'),specimens=('specimen_id','nunique'),defined_specimens=('value','count')).reset_index()
sp.to_csv(out/'per_specimen_endpoints.tsv.gz',sep='\t',index=False);pa.to_csv(out/'per_participant_region_endpoints.tsv.gz',sep='\t',index=False)
pa.groupby(keys+['region']).agg(mean=('value','mean'),median=('value','median'),participants=('participant','nunique'),defined_participants=('value','count')).reset_index().to_csv(out/'region_summary.tsv',sep='\t',index=False)
def seed(key):return int.from_bytes(hashlib.sha256(('|'.join(map(str,key))+'|20260914').encode()).digest()[:4],'little')
null_cache={};boot_cache={}
def assess(delta,key,bootstrap=2000):
    delta=np.asarray(delta,float);delta=delta[np.isfinite(delta)];n=len(delta)
    result=dict(n_participants=n,mean_difference=float(delta.mean()) if n else np.nan,median_difference=float(np.median(delta)) if n else np.nan,n_positive=int((delta>0).sum()),n_negative=int((delta<0).sum()),p=1.,estimable=n>=8,CI95_low=np.nan,CI95_high=np.nan)
    if n<8:return result
    if n not in null_cache:
        if n<=16:
            z=np.arange(2**n,dtype=np.uint32)[:,None];null_cache[n]=((z>>np.arange(n))&1).astype(np.int8)*2-1
        else:null_cache[n]=np.random.default_rng(seed(['null',n])).choice(np.array([-1,1],np.int8),size=(200000,n))
    means=null_cache[n]@delta/n;hits=int((abs(means)>=abs(delta.mean())-1e-12).sum())
    result['p']=hits/len(means) if n<=16 else (hits+1)/(len(means)+1)
    result['test']='exact_sign_flip' if n<=16 else '200000_sign_flips_plus_one'
    if (n,bootstrap) not in boot_cache:boot_cache[n,bootstrap]=np.random.default_rng(seed(['bootstrap',n,bootstrap])).multinomial(n,np.full(n,1/n),size=bootstrap)
    ci=np.quantile(boot_cache[n,bootstrap]@delta/n,[.025,.975]);result.update(CI95_low=float(ci[0]),CI95_high=float(ci[1]))
    return result
tests=[];differences=[]
primary=pa[pa.origin.isin(['All','Human','Microbial'])]
if not pilot:
    cap=primary[primary.cohort.eq('CAPSCAN')] if mode!='denovo' else primary.iloc[:0]  # Regional families wait for all Sage routes.
    for labels,gg in cap.groupby(keys,sort=False):
        info=dict(zip(keys,labels));is_residue=info['metric'].startswith('residue_')
        if info['metric']=='other_residue_fraction_nonKR':continue
        if is_residue and (info['unit']!='assigned_bare_sequence' or info['min_length']!=9):continue
        if not is_residue and info['unit']=='assigned_bare_sequence':continue
        pv=gg.pivot(index='participant',columns='region',values='value').reindex(columns=REGIONS)
        pv['proximal']=pv[['duodenum','jejunum']].mean(axis=1,skipna=False);pv['distal']=pv[['ileum','colon']].mean(axis=1,skipna=False)
        pv['capsule']=pv[['proximal','distal']].mean(axis=1,skipna=False)
        contrasts=list(itertools.combinations(REGIONS,2))+[('proximal','stool'),('distal','stool'),('capsule','stool'),('proximal','distal')]
        for a,b in contrasts:
            pair=pv[[a,b]].dropna();delta=pair[a]-pair[b]
            family=('all18_nonKR_composition' if is_residue else 'terminal_fraction')+'__length_'+str(info['min_length'])
            row=dict(info,contrast=a+'__minus__'+b,family=family,**assess(delta,[*labels,a,b],2000 if is_residue else 50000));tests.append(row)
            differences.extend(dict(info,contrast=row['contrast'],participant=str(pid),difference=float(v)) for pid,v in delta.items())
    # Threshold contrasts pair the same participant and region. PREDICT ALL uses
    # specimen means, not equal weights for clinical groups or raw acquisitions.
    all_patient=sp.groupby(keys+['participant'],dropna=False).value.mean().reset_index();all_patient['region']='ALL_observed_specimens'
    threshold_data=pd.concat([primary,all_patient[all_patient.origin.isin(['All','Human','Microbial'])]],ignore_index=True)
    threshold_data=threshold_data[threshold_data.metric.isin(['nonKR_fraction_all','H_fraction_all','H_fraction_nonKR']) & threshold_data.unit.isin(['canonical_peptide','accepted_psm_or_prediction_row'])]
    comparisons=[('historical_all_lengths','denovo_highconfidence',0),('same_9_50_lengths','denovo_highconfidence',9),('same_production_eligibility','denovo_highconfidence_production_eligible',9)]
    tk=['cohort','unit','origin','terminus','metric','participant','region']
    prod=threshold_data[threshold_data.route.eq('denovo')&threshold_data.min_length.eq(9)]
    for label,reference,lo in comparisons:
        ref=threshold_data[threshold_data.route.eq(reference)&threshold_data.min_length.eq(lo)]
        matched=prod.merge(ref,on=tk,how='outer',suffixes=('_top30','_reference'),indicator=True,validate='1:1')
        matched['difference']=matched.value_top30-matched.value_reference
        matched.to_csv(out/('threshold_pairs_'+label+'.tsv.gz'),sep='\t',index=False)
        for labels,gg in matched.groupby([x for x in tk if x!='participant'],sort=False):
            info=dict(zip([x for x in tk if x!='participant'],labels));family='selection_'+label+'__'+info['cohort']
            tests.append(dict(info,route='top30_vs_'+reference,min_length=9,contrast=label,family=family,**assess(gg.difference,[label,*labels],50000)))
    tf=pd.DataFrame(tests)
    def bh(p):
        p=np.asarray(p,float);order=np.argsort(p);v=np.empty(len(p));v[order]=np.minimum(1,np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]);return v
    tf['q']=tf.groupby('family').p.transform(bh);tf.to_csv(out/'complete_statistical_families.tsv',sep='\t',index=False)
    pd.DataFrame(differences).to_csv(out/'participant_paired_differences.tsv.gz',sep='\t',index=False)
    tf.groupby('family').agg(tests=('p','size'),estimable=('estimable','sum'),q_lt_05=('q',lambda v:int((v<.05).sum()))).reset_index().to_csv(out/'statistical_family_summary.tsv',sep='\t',index=False)
(out/'REPORT.md').write_text('''# Terminal residues and de novo selection sensitivity

Fractions are calculated within acquisition, averaged within canonical specimen,
then within participant/region. Paired differences retain the same participants.
Proximal and distal summaries require both constituent regions; the capsule
summary averages those two compartments. PREDICT ALL averages observed specimens
within patient and does not assert a clinical-group effect. Missing denominators
remain undefined. High-confidence all-length, common 9-50 and common production-
eligibility comparisons separate combined filtering from the score-cutoff effect.

All declared regional and compartment contrasts and all eighteen non-K/R residues
are retained. Primary inference uses conservative All/Human/Microbial origin and
canonical-peptide/prediction-or-PSM-row units; eighteen-residue composition uses
assigned bare sequences to retain separate I/L labels. Unknown/mixed/absent-context
profiles remain in the descriptive source tables. Boundary/unmapped contexts are
explicit in the source counts. N-first is the peptide residue, not the preceding
cleavage P1. Neither terminal profile establishes a protease mechanism.

Tests use paired sign flips: exact up to sixteen participants, otherwise 200,000
draws with plus-one correction. At least eight paired participants are required;
non-estimable planned tests remain in the family with p=1. BH families include all
routes/origins/termini/contrasts in the stated endpoint/length family. Participant
bootstrap intervals use 50,000 resamples for the compact terminal endpoints and
2,000 for full residue-composition sensitivities. Tests and pointwise intervals
are exploratory and assume sign symmetry under the null. No pooled-peptide
pseudoreplication or biological independence of the de novo/Sage methods is claimed.
Scientific review, route corroboration, figures and prose remain pending.
''')
dump(out/'COMPLETE.json',dict(status='PASS',pilot=pilot,mode=mode,acquisitions=len(tasks),
    production_route_independently_reproduced=len(production_checks),scientific_review='PENDING',source_receipts=receipts,
    job_id=os.environ.get('SLURM_JOB_ID'),outputs=[meta(p) for p in out.iterdir() if p.is_file() and p.name!='COMPLETE.json']))
print('TERMINAL STATISTICS PASS',mode,flush=True)
