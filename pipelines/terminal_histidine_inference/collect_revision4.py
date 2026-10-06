#!/usr/bin/env python3
"""Participant-level source enrichment, motifs and complementary-site controls."""
import itertools,json,os,sqlite3,sys
from collections import defaultdict
from pathlib import Path
import numpy as np
import pandas as pd
RUN=Path(sys.argv[1]).resolve();sys.path.insert(0,str(RUN/'frozen'))
from lib import dump,meta,verify
from histidine_helpers import AA,REGIONS,bh,signflip,protein_key,background,opportunities
from context_helpers import fasta
m=json.loads((RUN/'inputs/manifest_revision4.json').read_text());stage=int(sys.argv[2])
pilot=stage==0;cohort=['PILOT','CAPSCAN','PREDICT','CAPSCAN','PREDICT','CAPSCAN'][stage]
phases=['available','unspecific'] if stage in [0,3,4] else ['available']
out=RUN/'results_revision4'/['pilot','available_CAPSCAN','available_PREDICT','full_CAPSCAN','full_PREDICT','statistics_smoke_CAPSCAN'][stage]
out.mkdir(exist_ok=True);assert not (out/'COMPLETE.json').exists()
for x in m['code']+m['inputs']:verify(x)
primary=['denovo_score09','denovo_top30','sage_semi','sage_semi_bothq','sage_unspecific','sage_unspecific_bothq']
tasks=m['tasks'][:4] if pilot else m['statistics_smoke_tasks'] if stage==5 else [t for t in m['tasks'] if t['cohort']==cohort]
enrich=[];motifs=[];reciprocal=[];concentration=[];overlap=[];mapping=[];receipts=[];unresolved_P1=[];P1_coverage=[]
db=sqlite3.connect(out/'boundary_recurrence.sqlite');db.execute('PRAGMA journal_mode=OFF')
site_columns=['route','policy','protein_key','position','origin','motifs','N_count','C_count','both_independent_spectra','N_scans','C_scans']
def load_output(c,name):
    x=next(x for x in c['outputs'] if Path(x['path']).name==name);verify(x)
    try:return pd.read_csv(x['path'],sep='\t',keep_default_na=False)
    except pd.errors.EmptyDataError:return pd.DataFrame()
def stamp(frame,t):
    for k,v in dict(cohort=t['cohort'],sample=t['sample'],**t['design']).items():frame[k]=v
    return frame
for t in tasks:
    siteparts=[];peptideparts=[]
    for phase in phases:
        p=Path(t['receipts'][phase]);c=json.loads(p.read_text());assert c['status']=='PASS';receipts.append(meta(p))
        assert c['sample']==t['sample'] and c['source_fasta']['sha256']==t['fasta']['sha256']
        f=load_output(c,'source_enrichment.tsv')
        if len(f):
            f=f[f.origin.isin(['All','Human','Microbial']) & f.population.eq('all_nonKR')].copy()
            enrich.append(stamp(f,t))
        f=load_output(c,'motif_counts.tsv.gz')
        if len(f) and t['cohort']=='CAPSCAN':
            f=f[f.route.isin(primary)&f.origin.isin(['All','Human','Microbial'])&f.population.isin(['all_nonKR','H_only'])].copy()
            motifs.append(stamp(f,t))
        f=load_output(c,'route_summary.tsv');mapping.append(stamp(f,t))
        f=load_output(c,'peptide_mapping.tsv.gz')
        if len(f):peptideparts.append(f[f.route.isin(primary)].copy())
        x=next(x for x in c['outputs'] if Path(x['path']).name=='site_evidence.tsv.gz');verify(x)
        try:
            ss=pd.read_csv(x['path'],sep='\t',usecols=site_columns,keep_default_na=False)
        except (pd.errors.EmptyDataError,ValueError):
            # Empty acquisition files were emitted without columns by the frozen first worker.
            if Path(x['path']).stat().st_size>100:raise
            ss=pd.DataFrame(columns=site_columns)
        ss=ss[ss.route.isin(primary)&ss.policy.eq('unambiguous')&((ss.N_count>0)|(ss.C_count>0))].copy()
        if len(ss):siteparts.append(ss)
    if peptideparts:
        pp=pd.concat(peptideparts,ignore_index=True);sets={route:set(g.peptide) for route,g in pp.groupby('route')}
        for dn in ['denovo_score09','denovo_top30']:
            for sage in [r for r in primary if r.startswith('sage_') and r in sets]:
                if dn not in sets:continue
                for label,peps in [('denovo_only',sets[dn]-sets[sage]),('sage_only',sets[sage]-sets[dn]),('agreeing',sets[dn]&sets[sage])]:
                    for a in AA:
                        overlap.append(dict(cohort=t['cohort'],sample=t['sample'],**t['design'],comparison=dn+'__'+sage,
                            population=label,terminal_residue=a,count=sum(p[-1]==a for p in peps),denominator=len(peps),
                            IL_note='I/L collapsed as L; identity-overlap composition is not 20 distinguishable residues'))
    if not siteparts:continue
    ss=pd.concat(siteparts,ignore_index=True)
    # P1 is unresolved for X/nonstandard source symbols. Native I/L alternatives
    # alone are compatible with the canonical full-protein key and H classification.
    ss['P1']=ss.motifs.map(lambda text: ''.join(sorted({w[3] for w in text.split(';')})))
    resolved=ss.P1.map(lambda a:bool(a) and set(a)<=set(AA) and (len(a)==1 or set(a)==set('IL')))
    ss['P1_class']=np.where(~resolved,'UNRESOLVED_P1',np.where(ss.P1.eq('H'),'H',np.where(ss.P1.isin(['K','R']),'KR','OTHER_STANDARD_NONKR')))
    for labels,g in ss.groupby(['route','P1_class'],sort=False):
        P1_coverage.append(dict(cohort=t['cohort'],sample=t['sample'],**t['design'],route=labels[0],P1_class=labels[1],sites=len(g),N_supported_sites=int((g.N_count>0).sum()),C_supported_sites=int((g.C_count>0).sum())))
    if (~resolved).any():
        unresolved_P1.append(stamp(ss.loc[~resolved].copy(),t))
    ss=ss[resolved & ~ss.P1.isin(['K','R'])].copy();ss['is_H']=ss.P1.eq('H')
    recurrence=ss[['route','protein_key','position','origin','P1','N_count','C_count','both_independent_spectra']].copy()
    stamp(recurrence,t).to_sql('sites',db,if_exists='append',index=False)
    # Exact boundary agreement between methods, in addition to peptide agreement.
    for dn in ['denovo_score09','denovo_top30']:
        for sage in [x for x in primary if x.startswith('sage_')]:
            a=ss[ss.route.eq(dn)];b=ss[ss.route.eq(sage)]
            if a.empty or b.empty:continue
            for side in ['N','C']:
                aset=set(zip(a.loc[a[side+'_count']>0,'protein_key'],a.loc[a[side+'_count']>0,'position']))
                bset=set(zip(b.loc[b[side+'_count']>0,'protein_key'],b.loc[b[side+'_count']>0,'position']))
                overlap.append(dict(cohort=t['cohort'],sample=t['sample'],**t['design'],comparison=dn+'__'+sage,
                    population='exact_unambiguous_'+side+'_boundary',count=len(aset&bset),denominator=len(aset),
                    sage_denominator=len(bset),terminal_residue='ALL_NONKR'))
    verify(t['fasta']);seqs={};audit={}
    for header,seq in fasta(t['fasta']['path'],audit):seqs.setdefault(protein_key(seq),set()).add(seq)
    assert audit['sha256']==t['fasta']['sha256'];bgcache={};opcache={}
    for route,gg in ss.groupby('route',sort=False):
        enzyme='semi' if route.startswith('sage_semi') else 'unspecific'
        for side in ['N','C']:
            g=gg[gg[side+'_count']>0]
            for origin in ['All','Human','Microbial']:
                z=g if origin=='All' else g[g.origin.eq(origin)]
                counts=z.groupby('protein_key').agg(sites=('position','size'),H=('is_H','sum')).reset_index()
                expected=[]
                for row in counts.itertuples():
                    prob=[]
                    for seq in seqs[row.protein_key]:
                        key=(seq,enzyme,side)
                        if key not in bgcache:bgcache[key]=background(seq,enzyme,side)
                        v=bgcache[key]['eligible'][3];prob.append(v[AA.index('H')]/v.sum() if v.sum() else np.nan)
                    expected.append(row.sites*np.mean(prob))
                counts['expected_H']=expected
                # Rank sources by all non-K/R site count, never by the H enrichment outcome.
                counts=counts.sort_values(['sites','protein_key'],ascending=[False,True])
                for remove in [0,1,5]:
                    remain=counts.iloc[remove:];den=int(remain.sites.sum());obs=int(remain.H.sum());exp=float(remain.expected_H.sum())
                    concentration.append(dict(cohort=t['cohort'],sample=t['sample'],**t['design'],route=route,side=side,origin=origin,
                        removed_top_proteins=remove,excluded_protein_keys=';'.join(counts.iloc[:remove].protein_key),proteins=len(remain),sites=den,H=obs,expected_H=exp,
                        H_fraction=obs/den if den else np.nan,expected_fraction=exp/den if den else np.nan,
                        enrichment=obs/exp if exp else np.nan,top1_site_share=float(counts.sites.iloc[0]/counts.sites.sum()) if len(counts) else np.nan))
        # Matched observed non-H boundaries in the same protein, opportunity class
        # and observed independent-spectrum-depth bin. This is a descriptive control.
        z=gg.copy();availability=[]
        for row in z.itertuples():
            key=(row.protein_key,enzyme)
            if key not in opcache:
                ns=set();cs=set()
                for seq in seqs[row.protein_key]:ns.update(opportunities(seq,enzyme,'N'));cs.update(opportunities(seq,enzyme,'C'))
                opcache[key]=(ns,cs)
            availability.append(str(int(row.position in opcache[key][0]))+str(int(row.position in opcache[key][1])))
        z['opportunity']=availability
        def depth(row):
            n=len({x for x in (str(row.N_scans)+';'+str(row.C_scans)).split(';') if x})
            return '1' if n<=1 else '2_4' if n<=4 else '5_plus'
        z['depth']=z.apply(depth,axis=1)
        for origin in ['All','Human','Microbial']:
            g=z if origin=='All' else z[z.origin.eq(origin)];matched=[];htotal=int(g.is_H.sum())
            for _,cell in g.groupby(['protein_key','opportunity','depth'],sort=False):
                hh=cell[cell.is_H];other=cell[~cell.is_H]
                if len(hh) and len(other):
                    other_rate=float(other.both_independent_spectra.mean())
                    matched.extend((int(v),other_rate) for v in hh.both_independent_spectra)
            reciprocal.append(dict(cohort=t['cohort'],sample=t['sample'],**t['design'],route=route,origin=origin,
                H_sites=htotal,matched_H_sites=len(matched),unmatched_H_sites=htotal-len(matched),
                H_both_rate=np.mean([x[0] for x in matched]) if matched else np.nan,
                matched_nonH_both_rate=np.mean([x[1] for x in matched]) if matched else np.nan))
    print('HISTIDINE COLLECT',t['cohort'],t['sample'],flush=True)

def write(name,rows):
    f=pd.concat(rows,ignore_index=True) if rows and isinstance(rows[0],pd.DataFrame) else pd.DataFrame(rows)
    f.to_csv(out/name,sep='\t',index=False);return f
ef=write('acquisition_enrichment.tsv.gz',enrich);mf=write('acquisition_motif_counts.tsv.gz',motifs)
cf=write('protein_concentration_sensitivities.tsv.gz',concentration);rf=write('matched_reciprocal_controls.tsv',reciprocal)
write('method_overlap.tsv.gz',overlap);write('mapping_coverage.tsv',mapping)
write('unresolved_P1_sites.tsv.gz',unresolved_P1);write('P1_classification_coverage.tsv',P1_coverage)
if db.execute("SELECT count(*) FROM sqlite_master WHERE type='table' AND name='sites'").fetchone()[0]:
    query='''SELECT cohort,route,protein_key,position,origin,P1,count(distinct sample) acquisitions,
      count(distinct CASE WHEN participant!='' THEN participant END) participants,
      count(distinct CASE WHEN N_count>0 THEN sample END) N_acquisitions,
      count(distinct CASE WHEN C_count>0 THEN sample END) C_acquisitions,
      count(distinct CASE WHEN N_count>0 AND participant!='' THEN participant END) N_participants,
      count(distinct CASE WHEN C_count>0 AND participant!='' THEN participant END) C_participants,
      sum(both_independent_spectra) acquisitions_with_distinct_spectra_both_sides
      FROM sites GROUP BY cohort,route,protein_key,position,origin,P1'''
    for i,chunk in enumerate(pd.read_sql_query(query,db,chunksize=100000)):
        chunk['cross_acquisition_both_sides']=(chunk.N_acquisitions>0)&(chunk.C_acquisitions>0)&(chunk.acquisitions>1)
        chunk['cross_donor_both_sides']=(chunk.N_participants>0)&(chunk.C_participants>0)&(chunk.participants>1)
        # Both sides may occur only in the same donor while another has opposite?
        # Distinct-donor existence follows union>1 with nonempty side sets.
        chunk.to_csv(out/'boundary_recurrence.tsv',sep='\t',index=False,mode='w' if i==0 else 'a',header=i==0)
db.close()

tests=[]
def ratio_interval(g):
    if not {'O','E'}<=set(g.columns):return {}
    a=g[['O','E']].replace([np.inf,-np.inf],np.nan).dropna().values;n=len(a)
    if not n:return {}
    means=a.mean(axis=0);r=dict(observed_fraction_mean=means[0],expected_fraction_mean=means[1],enrichment_ratio=means[0]/means[1] if means[1] else np.nan)
    if n>=8:
        weights=np.random.default_rng(20261114+n).multinomial(n,np.full(n,1/n),size=20000)
        v=weights@a/n;ratios=np.divide(v[:,0],v[:,1],out=np.full(len(v),np.nan),where=v[:,1]>0)
        valid=ratios[np.isfinite(ratios)];r['ratio_bootstrap_defined']=len(valid)
        if len(valid):r.update(enrichment_CI95_low=np.quantile(valid,.025),enrichment_CI95_high=np.quantile(valid,.975))
    return r

def participant(frame,key,value,all_specimens=False):
    f=frame.copy();f[value]=pd.to_numeric(f[value],errors='coerce');f=f[f.participant.astype(str).ne('')]
    sp=f.groupby(key+['cohort','participant','region','specimen_id'],dropna=False)[value].mean().reset_index()
    if all_specimens:return sp.groupby(key+['cohort','participant'],dropna=False)[value].mean().reset_index().assign(region='ALL_observed_specimens')
    return sp.groupby(key+['cohort','participant','region'],dropna=False)[value].mean().reset_index()
if len(ef):
    for k in ['observed','expected','denominator']:ef[k]=pd.to_numeric(ef[k],errors='coerce')
    ef['O']=ef.observed/ef.denominator.replace(0,np.nan);ef['E']=ef.expected/ef.denominator.replace(0,np.nan)
    key=['route','side','origin','residue'];a=participant(ef,key,'O');b=participant(ef,key,'E')
    pe=a.merge(b,on=key+['cohort','participant','region'],validate='1:1');pe['difference']=pe.O-pe.E
    pe['enrichment']=pe.O/pe.E.replace(0,np.nan);pe.to_csv(out/'participant_region_enrichment.tsv',sep='\t',index=False)
    if not pilot:
        # All represented route/origin/residue cells retained; unavailable cells
        # are added below as planned p=1 hypotheses for the primary grid.
        routes=['score_0_0.5','score_0.5_0.7','score_0.7_0.9','score_0.9_0.95','score_0.95_0.99','score_0.99_plus','denovo_score09','denovo_score09_eligible','denovo_top30','sage_tryptic','sage_tryptic_bothq']+(['sage_semi','sage_semi_bothq'] if cohort=='CAPSCAN' else [])+(['sage_unspecific','sage_unspecific_bothq'] if 'unspecific' in phases else []);origins=['All','Human','Microbial'];regions=REGIONS if cohort=='CAPSCAN' else ['ALL_observed_specimens']
        if cohort=='CAPSCAN':pp=pe
        else:
            ap=participant(ef,key,'O',True);bp=participant(ef,key,'E',True)
            pp=ap.merge(bp,on=key+['cohort','participant','region'],validate='1:1');pp['difference']=pp.O-pp.E
        lookup={labels:g for labels,g in pp.groupby(key+['region'],sort=False)}
        for route,side,origin,residue,region in itertools.product(routes,['N','C'],origins,AA,regions):
            g=lookup.get((route,side,origin,residue,region),pd.DataFrame({'difference':[]}))
            tests.append(dict(endpoint='source_enrichment',family='source_P1_all_residues',route=route,side=side,origin=origin,
                residue=residue,contrast=region,**signflip(g.difference),**ratio_interval(g)))
        h=pe[pe.residue.eq('H')]
        if cohort=='CAPSCAN':
            spatial_cells={labels:g for labels,g in h.groupby(['route','side','origin'],sort=False)}
            for labels in itertools.product(routes,['N','C'],origins):
                g=spatial_cells.get(labels)
                for metric in ['difference','enrichment']:
                    pv=g.pivot(index='participant',columns='region',values=metric).reindex(columns=REGIONS) if g is not None else pd.DataFrame(columns=REGIONS,dtype=float)
                    pv['proximal']=pv[['duodenum','jejunum']].mean(axis=1,skipna=False);pv['distal']=pv[['ileum','colon']].mean(axis=1,skipna=False)
                    for x,y in list(itertools.combinations(REGIONS,2))+[('proximal','stool'),('distal','stool'),('proximal','distal')]:
                        d=pv[x]-pv[y];tests.append(dict(endpoint='spatial_H_'+metric,family='H_spatial_all_routes',route=labels[0],side=labels[1],origin=labels[2],contrast=x+'__'+y,**signflip(d)))
                    slopes=[]
                    for _,row in pv[REGIONS[:4]].iterrows():
                        if row.notna().all():slopes.append(np.polyfit([1,2,3,4],row.values,1)[0])
                    tests.append(dict(endpoint='ordinal_fluid_H_'+metric,family='H_spatial_all_routes',route=labels[0],side=labels[1],origin=labels[2],contrast='ordinal1_to4_complete_fluid_regions',**signflip(slopes)))
if len(mf):
    for k in ['observed','expected','observed_total','expected_total']:mf[k]=pd.to_numeric(mf[k],errors='coerce')
    mf['difference']=mf.observed/mf.observed_total.replace(0,np.nan)-mf.expected/mf.expected_total.replace(0,np.nan)
    key=['route','side','origin','population','position','residue'];pm=participant(mf,key,'difference')
    pm.to_csv(out/'participant_motif_differences.tsv.gz',sep='\t',index=False)
    # Full matrices are retained; motif inference excludes the fixed H-conditioning P1.
    if not pilot and cohort=='CAPSCAN':
        cells={labels:g for labels,g in pm.groupby(key,sort=False)}
        motif_routes=['denovo_score09','denovo_top30','sage_semi','sage_semi_bothq']+(['sage_unspecific','sage_unspecific_bothq'] if 'unspecific' in phases else [])
        for labels in itertools.product(motif_routes,['N','C'],['All','Human','Microbial'],['all_nonKR','H_only'],[-4,-3,-2,-1,1,2,3,4],AA):
            if labels[3]=='H_only' and labels[4]==-1:continue
            g=cells.get(labels)
            pv=g.pivot(index='participant',columns='region',values='difference').reindex(columns=REGIONS) if g is not None else pd.DataFrame(columns=REGIONS,dtype=float)
            for compartment,columns in [('proximal',REGIONS[:2]),('distal',REGIONS[2:4]),('stool',['stool'])]:
                d=pv[columns].mean(axis=1,skipna=False)
                tests.append(dict(endpoint='motif',family='all_motif_positions_residues',**dict(zip(key,labels)),contrast=compartment,**signflip(d)))
if len(cf):
    cf['difference']=cf.H_fraction-cf.expected_fraction
    pc=participant(cf,['route','side','origin','removed_top_proteins'],'difference')
    pc.to_csv(out/'participant_concentration_sensitivity.tsv',sep='\t',index=False)
if len(rf):
    rf['difference']=rf.H_both_rate-rf.matched_nonH_both_rate
    pr=participant(rf,['route','origin'],'difference');pr.to_csv(out/'participant_reciprocal_controls.tsv',sep='\t',index=False)
    if not pilot:
        reciprocal_cells={labels:g for labels,g in pr.groupby(['route','origin','region'],sort=False)}
        reciprocal_routes=[route for route in primary if route in routes]
        for labels in itertools.product(reciprocal_routes,origins,REGIONS if cohort=='CAPSCAN' else ['ALL_observed_specimens']):
            if cohort=='PREDICT':
                # Match the all-specimen patient weighting used for source tests.
                pg=participant(rf,['route','origin'],'difference',True)
                g=pg[(pg.route==labels[0])&(pg.origin==labels[1])]
            else:g=reciprocal_cells.get(labels,pd.DataFrame({'difference':[]}))
            tests.append(dict(endpoint='matched_reciprocal',family='all_reciprocal_controls',route=labels[0],origin=labels[1],contrast=labels[2],**signflip(g.difference)))
if tests:
    tf=pd.DataFrame(tests);tf['q']=tf.groupby('family').p.transform(bh);tf.to_csv(out/'complete_statistical_families.tsv',sep='\t',index=False)
    tf.groupby('family').agg(tests=('p','size'),estimable=('estimable','sum'),q_lt05=('q',lambda x:int((x<.05).sum()))).reset_index().to_csv(out/'statistical_family_summary.tsv',sep='\t',index=False)
if stage==5:
    assert len({t['design']['participant'] for t in tasks})==8
    estimable=tf[tf.estimable.eq(True)]
    for family in ['source_P1_all_residues','H_spatial_all_routes','all_motif_positions_residues']:
        assert len(estimable[estimable.family.eq(family)])>0,'Inferential branch not exercised: '+family
    assert len(estimable[estimable.endpoint.str.startswith('ordinal_fluid')])>0,'Ordinal branch not exercised'
    assert tf.get('enrichment_CI95_low',pd.Series(dtype=float)).notna().any(),'Ratio-bootstrap branch not exercised'
(out/'REPORT.md').write_text('''# Histidine defence collection

Source-normalized site and motif evidence is averaged within specimen and then
participant/region. H-source tests use observed minus source-expected fractions;
ratios are retained separately. Intervals from signflip() are donor-bootstrap
intervals for mean differences, not ratio confidence intervals. CapScan anatomy
and PREDICT observed clinical-region mixtures remain separate. All tests are
exploratory sign-symmetry tests; n<8 is non-estimable/p=1. Full families and null
outcomes remain in source tables. Primary enrichment uses uniquely localized
sites; ambiguous evidence is retained in the underlying A27 acquisition outputs.

Reciprocal controls match the same source protein, N/C opportunity and observed
spectrum-depth bin. Their denominator is observed non-K/R boundaries, not all
proteolytic events. Concentration controls remove sources ranked by all non-K/R
site count. Same-protein source backgrounds condition on a de novo-selected FASTA;
they do not remove selection, detectability, digestion, handling or ex-vivo effects.
Unresolved P1 source symbols are exported with complete classification counts and
excluded from resolved H/non-H concentration, reciprocity and boundary overlap.
Taxonomic specificity, raw-MS covariates, measured spectrum review and final
scientific interpretation remain distinct checks. No figure or prose is finalized.
''')
dump(out/'COMPLETE.json',dict(status='PASS',pilot=pilot,cohort=cohort,phases=phases,acquisitions=len(tasks),
    scientific_review='PENDING',validation_only=stage in [0,5],source_receipts=receipts,job_id=os.environ.get('SLURM_JOB_ID'),
    outputs=[meta(p) for p in out.iterdir() if p.is_file() and p.name!='COMPLETE.json']))
print('HISTIDINE COLLECTION PASS',cohort,phases,flush=True)
