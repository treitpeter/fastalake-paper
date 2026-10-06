#!/usr/bin/env python3
"""Compare fixed biology outputs and matched CapScan attribution opportunities."""
import json,os,sys
from pathlib import Path
import numpy as np
import pandas as pd
RUN=Path(sys.argv[1]).resolve();sys.path.insert(0,str(RUN/'frozen'))
from lib import dump,meta,verify

def read(p):return pd.read_csv(p,sep='\t')
def write(df,p):df.to_csv(p,sep='\t',index=False)
def aggregate(v,d,keys):
    groups=d.groupby(keys,sort=True).indices
    a=np.stack([v[:,ix].mean(axis=1) for ix in groups.values()],axis=1)
    return a,pd.DataFrame(list(groups),columns=keys)

def region_state(folder,eligible):
    d=folder/'source_data/Figure5';z=np.load(d/'function_species.npz')
    ac=read(d/'acquisition_composition_local_candidate.tsv').set_index('sample_id').loc[z['samples']].reset_index()
    v=z['values']/ac.total_LFQ.to_numpy();ko=z['ko'];species=z['species']
    vs,sm=aggregate(v,ac,['participant','region','specimen_id']);vr,cm=aggregate(vs,sm,['participant','region'])
    nonhost=~np.isin(species,['Host','Mixed host/microbial'])
    resolved=nonhost & ~np.isin(species,['No species annotation','Conflicting species','Unassigned'])
    old=read(d/'contributor_region_cells.tsv').set_index(['KO','participant','region'])
    rows=[]
    for k in sorted(eligible):
        ix=np.flatnonzero((ko==k)&resolved);allix=np.flatnonzero((ko==k)&nonhost)
        totals=vr[allix].sum(axis=0);rr=vr[ix].sum(axis=0)
        for j,c in cm.iterrows():
            key=(k,c.participant,c.region);leader='';gap=np.nan;fraction=rr[j]/totals[j] if totals[j]>0 else np.nan
            if rr[j]>0:
                order=sorted(ix,key=lambda i:(-vr[i,j],species[i]));leader=str(species[order[0]])
                gap=(vr[order[0],j]-(vr[order[1],j] if len(order)>1 else 0))/rr[j]
                prior=old.loc[key];assert prior.leading_species==leader
                assert np.isclose(prior.gap,gap,rtol=1e-10,atol=1e-12)
                assert np.isclose(prior.species_resolved_fraction,fraction,rtol=1e-10,atol=1e-12)
            else:assert key not in old.index
            rows.append(dict(KO=k,participant=int(c.participant),region=c.region,
                nonhost_function_fraction=totals[j],resolved_fraction_total=rr[j],species_resolved_fraction=fraction,
                leading_species=leader,leader_gap=gap,strict=bool(rr[j]>0 and fraction>=.2 and gap>=.05),
                state='resolved' if rr[j]>0 else ('unresolved_positive_signal' if totals[j]>0 else 'not_detected')))
    return pd.DataFrame(rows)

m=json.loads((RUN/'inputs/collector_manifest.json').read_text())
for x in m['receipts']+m['code']:verify(x)
for x in m['receipts']:
    c=json.loads(Path(x['path']).read_text());assert c['status']=='PASS'
    for y in c['outputs']:verify(y)
out=RUN/'results';folders={method:out/('CAPSCAN__'+method) for method in ['original','simplified']}
elig={method:set(read(p/'source_data/Figure5/contributor_cells.tsv').KO) for method,p in folders.items()}
common=elig['original']&elig['simplified'];union=elig['original']|elig['simplified']
assert common
write(pd.DataFrame([dict(KO=k,original_eligible=k in elig['original'],simplified_eligible=k in elig['simplified'],common_eligible=k in common) for k in sorted(union)]),out/'contributor_eligibility.tsv')
states={method:region_state(p,union) for method,p in folders.items()}
paired=states['original'].merge(states['simplified'],on=['KO','participant','region'],suffixes=('_original','_simplified'),validate='1:1',how='outer',indicator=True)
assert paired['_merge'].eq('both').all();paired=paired.drop(columns='_merge')
paired['common_eligible_function']=paired.KO.isin(common)
paired['both_resolved']=paired.state_original.eq('resolved')&paired.state_simplified.eq('resolved')
paired['leader_changed']=np.where(paired.both_resolved,paired.leading_species_original.ne(paired.leading_species_simplified),np.nan)
paired['both_strict']=paired.strict_original & paired.strict_simplified
write(paired,out/'contributor_method_cells.tsv')
transitions=paired.groupby(['common_eligible_function','state_original','state_simplified']).size().rename('cells').reset_index()
write(transitions,out/'contributor_resolution_transitions.tsv')
comparisons=[];donors=[];summaries=[]
for scope,cut in [('common_functions',paired.common_eligible_function),('eligible_union_sensitivity',pd.Series(True,index=paired.index))]:
    pp=paired[cut].copy()
    for threshold,ok in [('resolved',pp.both_resolved),('strict_both_methods',pp.both_strict)]:
        rr=pp[ok]
        ds=rr.groupby(['participant','region']).leader_changed.agg(['mean','size']).reset_index().rename(columns={'mean':'leader_switch_fraction','size':'cells'})
        ds['scope']=scope;ds['threshold']=threshold;donors.append(ds)
        dv=ds.groupby('participant').leader_switch_fraction.mean()
        summaries.append(dict(endpoint='between_method_leader_switch',scope=scope,threshold=threshold,cells=len(rr),participants=len(dv),participant_mean=float(dv.mean()),pooled_cell_mean=float(rr.leader_changed.mean())))
    stool=pp[pp.region.eq('stool')].drop(columns='region')
    cp=pp[pp.region.ne('stool')].merge(stool,on=['KO','participant'],suffixes=('_capsule','_stool'),validate='m:1')
    for threshold,ok in [('resolved',cp.both_resolved_capsule & cp.both_resolved_stool),('strict_both_methods',cp.both_strict_capsule & cp.both_strict_stool)]:
        cc=cp[ok].copy()
        for method in ['original','simplified']:
            cc['agreement_'+method]=cc['leading_species_'+method+'_capsule'].eq(cc['leading_species_'+method+'_stool'])
        dd=cc.groupby(['participant','region'])[['agreement_original','agreement_simplified']].mean().reset_index()
        dd['cells']=cc.groupby(['participant','region']).size().to_numpy();dd['scope']=scope;dd['threshold']=threshold
        comparisons.append(dd)
        for method in ['original','simplified']:
            dv=dd.groupby('participant')['agreement_'+method].mean()
            summaries.append(dict(endpoint='capsule_stool_agreement_'+method,scope=scope,threshold=threshold,cells=len(cc),participants=len(dv),participant_mean=float(dv.mean()),pooled_cell_mean=float(cc['agreement_'+method].mean())))
        cc['scope']=scope;cc['threshold']=threshold
        write(cc,out/('capsule_stool_cells_'+scope+'_'+threshold+'.tsv'))
write(pd.concat(donors,ignore_index=True),out/'between_method_leader_switch_donors.tsv')
write(pd.concat(comparisons,ignore_index=True),out/'matched_capsule_stool_donors.tsv')
write(pd.DataFrame(summaries),out/'matched_leader_summary.tsv')
# Retain every original inferential result and align only explicit matching endpoints.
definitions=[('PREDICT','Figure4','clinical_paired_statistics.tsv',['metric']),
 ('PREDICT','Figure4','holobiont_EDS_context_tests.tsv',['endpoint','model']),
 ('CAPSCAN','Figure5','target_function_tests.tsv',['KO','metric','capsule_region']),
 ('CAPSCAN','Figure5','target_depth_host_sensitivity_local_candidate.tsv',['KO','region']),
 ('CAPSCAN','Figure5','BSH_individual_metabolite_tests.tsv',['KO','metabolite']),
 ('CAPSCAN','Figure5','own_stool_donors.tsv',['participant']),
 ('CAPSCAN','Figure5','variance_partition.tsv',['population'])]
stats=[]
for co,fig,name,keys in definitions:
    a=read(out/(co+'__original')/'source_data'/fig/name);b=read(out/(co+'__simplified')/'source_data'/fig/name)
    joined=a.merge(b,on=keys,suffixes=('_original','_simplified'),how='outer',validate='1:1',indicator=True)
    dest=out/(co+'_'+name);write(joined,dest)
    stats.append(dict(cohort=co,table=name,original_rows=len(a),simplified_rows=len(b),matched_rows=int(joined['_merge'].eq('both').sum()),original_only=int(joined['_merge'].eq('left_only').sum()),simplified_only=int(joined['_merge'].eq('right_only').sum())))
write(pd.DataFrame(stats),out/'statistical_endpoint_alignment.tsv')
lines=['# Grouping robustness: fixed biology and matched CapScan contributors','',
 'Both simplified baseline replays reproduced the predeclared manuscript tables. The original rule was replayed on the same acquisitions with the same annotation definitions, local candidate composition, participant aggregation and model code.', '',
 f'CapScan contributor eligibility: {len(elig["original"]):,} original, {len(elig["simplified"]):,} simplified, {len(common):,} common functions. Method-specific eligibility is retained as a sensitivity result.', '',
 '| Matched endpoint | Threshold | Cells | Participants | Equal-participant mean |','|---|---|---:|---:|---:|']
for r in summaries:
    if r['scope']=='common_functions':lines.append(f"| {r['endpoint']} | {r['threshold']} | {r['cells']:,} | {r['participants']} | {r['participant_mean']:.2%} |")
lines += ['', 'Capsule/stool comparisons require resolution in both compartments under both grouping methods. Strict comparisons additionally require ≥20% species-resolved signal and ≥5% leader margin in all four method/compartment cells. The estimator averages functions within participant/region, then regions within participant, then participants. Pooled cell means are saved separately.', '',
 'Not-detected and positive-but-unresolved transitions are explicit. Representative taxon labels are not proof of unique organism origin. Overlapping function/region observations are not treated as independent biological participants.', '',
 'The Hartout H02 workstream supplies the separate HSTOOL technical/workflow leader benchmark; it is not subtracted from CapScan biological differences. Statistical interpretation, plots and manuscript integration remain pending.']
(out/'REPORT.md').write_text('\n'.join(lines)+'\n')
outputs=[meta(p) for p in out.iterdir() if p.is_file() and p.name!='AGGREGATE_VALIDATION.json']
dump(out/'AGGREGATE_VALIDATION.json',dict(status='PASS',replays=4,baseline_reproductions=2,common_eligible_functions=len(common),all_region_cells_paired=True,job_id=os.environ.get('SLURM_JOB_ID'),outputs=outputs,scientific_review='PENDING'))
print('BIOLOGY MATCHED COLLECTOR PASS',flush=True)
