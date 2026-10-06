#!/usr/bin/env python3
"""Declared acquisition -> specimen -> participant terminal-residue summaries."""
import json,os,sys
from pathlib import Path
import numpy as np
import pandas as pd
RUN=Path(sys.argv[1]).resolve();sys.path.insert(0,str(RUN/'frozen'))
from lib import dump,meta,verify
index=int(sys.argv[2]);pilot=index<2;co=['CAPSCAN','PREDICT'][index%2]
m=json.loads((RUN/'inputs/manifest.json').read_text());cm=json.loads((RUN/'inputs/collector_manifest.json').read_text())
for x in cm['code']+[cm['source_manifest']]:verify(x)
out=RUN/'results'/(('pilot_summary_' if pilot else 'summary_')+co);out.mkdir(parents=True,exist_ok=True);assert not (out/'COMPLETE.json').exists()
receipts=[];profiles=[];context=[];overlap=[];design=[]
tasks=[t for t in m['tasks'] if t['cohort']==co]
if pilot:tasks=tasks[:1]
for t in tasks:
    folder=RUN/'results'/co/t['sample'];c=json.loads((folder/'COMPLETE.json').read_text());assert c['status']=='PASS'
    for x in c['outputs']:verify(x)
    receipts.append(meta(folder/'COMPLETE.json'));d=dict(sample=t['sample'],**t['design']);design.append(d)
    for name,collection in [('terminal_counts.tsv',profiles),('context_summary.tsv',context),('route_identity.tsv',overlap)]:
        f=pd.read_csv(folder/name,sep='\t',keep_default_na=False)
        for k,v in d.items():f[k]=str(v)
        collection.append(f)
design=pd.DataFrame(design);design.to_csv(out/'acquisition_design.tsv',sep='\t',index=False)
f=pd.concat(profiles,ignore_index=True);ctx=pd.concat(context,ignore_index=True);ov=pd.concat(overlap,ignore_index=True)
keys=['sample','route','min_length','origin_policy','unit','origin','terminus']
assert not f.duplicated(keys+['residue']).any()
f['count']=f['count'].astype(float);f['denominator']=f['denominator'].astype(float)
known_residues=set('ACDEFGHIKLMNPQRSTVWY')|{'I/L'};nonkr=known_residues-{'K','R'}
nd=f[f.residue.isin(nonkr)].groupby(keys)['count'].sum().rename('nonkr_denominator').reset_index()
f=f.merge(nd,on=keys,validate='m:1');f['fraction_all']=f['count']/f.denominator.replace(0,np.nan)
f['fraction_nonkr']=np.where(f.residue.isin(nonkr),f['count']/f.nonkr_denominator.replace(0,np.nan),np.nan)
# Source task tables preserve zeros for every residue within every observed origin.
# Context summaries additionally preserve absent origin/route cells with zero peptides.
f.to_csv(out/'per_acquisition_residues.tsv.gz',sep='\t',index=False)
ctx.to_csv(out/'per_acquisition_context_bounds.tsv.gz',sep='\t',index=False)
ov.to_csv(out/'per_acquisition_route_identity.tsv.gz',sep='\t',index=False)
missing=design[design.participant.eq('')];missing.to_csv(out/'participant_exclusions.tsv',sep='\t',index=False)
population=f[f.participant.ne('')].copy()
endpoint=['route','min_length','origin_policy','unit','origin','terminus','residue']
sk=endpoint+['participant','region','specimen_id']
sp=population.groupby(sk,dropna=False).agg(fraction_all=('fraction_all','mean'),fraction_nonkr=('fraction_nonkr','mean'),n_acquisitions=('sample','nunique'),n_defined_all=('fraction_all','count'),n_defined_nonkr=('fraction_nonkr','count')).reset_index()
sp.to_csv(out/'per_specimen_residues.tsv.gz',sep='\t',index=False)
pk=endpoint+['participant','region']
pa=sp.groupby(pk,dropna=False).agg(fraction_all=('fraction_all','mean'),fraction_nonkr=('fraction_nonkr','mean'),n_specimens=('specimen_id','nunique'),n_defined_all=('fraction_all','count'),n_defined_nonkr=('fraction_nonkr','count')).reset_index()
pa.to_csv(out/'per_participant_region_residues.tsv.gz',sep='\t',index=False)
summary=pa.groupby(endpoint+['region']).agg(mean_fraction_all=('fraction_all','mean'),median_fraction_all=('fraction_all','median'),mean_fraction_nonkr=('fraction_nonkr','mean'),median_fraction_nonkr=('fraction_nonkr','median'),participants=('participant','nunique'),defined_participants_all=('fraction_all','count'),defined_participants_nonkr=('fraction_nonkr','count')).reset_index()
summary.to_csv(out/'residue_summary.tsv',sep='\t',index=False)
primary=pa[pa.origin_policy.eq('conservative') & pa.unit.eq('canonical_peptide') & pa.residue.eq('H')]
primary.to_csv(out/'histidine_participant_regions.tsv',sep='\t',index=False)
paired_n=0
if co=='CAPSCAN':
    stool=primary[primary.region.eq('stool')].drop(columns='region')
    caps=primary[primary.region.ne('stool')]
    pairs=caps.merge(stool,on=endpoint+['participant'],suffixes=('_capsule','_stool'),how='outer',validate='m:1',indicator=True)
    for metric in ['fraction_all','fraction_nonkr']:pairs[metric+'_difference']=pairs[metric+'_capsule']-pairs[metric+'_stool']
    pairs.to_csv(out/'histidine_capsule_stool_all_pair_states.tsv',sep='\t',index=False)
    matched=pairs[pairs['_merge'].eq('both')]
    # Average observed capsule regions within participant before population summaries.
    dk=endpoint+['participant'];donors=matched.groupby(dk).agg(delta_all=('fraction_all_difference','mean'),delta_nonkr=('fraction_nonkr_difference','mean'),regions=('region','nunique')).reset_index()
    donors.to_csv(out/'histidine_capsule_stool_per_participant.tsv',sep='\t',index=False)
    ds=donors.groupby(endpoint).agg(mean_delta_all=('delta_all','mean'),median_delta_all=('delta_all','median'),mean_delta_nonkr=('delta_nonkr','mean'),median_delta_nonkr=('delta_nonkr','median'),participants=('participant','nunique'),defined_participants_all=('delta_all','count'),defined_participants_nonkr=('delta_nonkr','count')).reset_index()
    ds.to_csv(out/'histidine_capsule_stool_summary.tsv',sep='\t',index=False);paired_n=donors.participant.nunique()
report=f'''# {co} terminal-residue analysis

All {len(design)} acquisition summaries passed their independent peptide recount,
actual-parent and digestion checks. There are {design.participant.replace('',np.nan).nunique()}
recorded participants; {len(missing)} acquisitions without participant IDs remain
in acquisition summaries. CapScan paired summaries include {paired_n} participants.

Use assigned-bare-sequence profiles to display all 18 assigned non-K/R residues;
canonical profiles explicitly combine I/L. N_first is the peptide's first residue.
N_cleavage_P1 is the preceding parent residue and weights ambiguous eligible
contexts equally. Protein ends and unmapped predictions remain separate. Histidine
counts are shown relative to all observations and to the non-K/R denominator.
Any/all parent-context bounds distinguish compatible from unambiguous internal cuts.

Primary source attribution requires every local compatible parent to be readable;
otherwise it is Unknown. Mixed source and the historical readable-only sensitivity
are retained. These labels do not establish unique organism or protease origin.
Every route, residue and declared denominator is retained. De novo comparisons
use 9–50 and their production selection, without claiming database-q equivalence.

Fractions are averaged within specimen, then participant/region, then participants.
Capsule/stool differences use matched participant cells and average capsule regions
within participant. Missing denominators are undefined, not zero. No pooled PSM
inference or post-selected residue p values. Scientific review, plots and manuscript
integration remain pending; do not call histidine enrichment before that review.
'''
(out/'REPORT.md').write_text(report)
dump(out/'COMPLETE.json',dict(status='PASS',pilot=pilot,cohort=co,acquisitions=len(design),participants=int(design.participant.replace('',np.nan).nunique()),missing_participant_acquisitions=len(missing),paired_capscan_participants=int(paired_n),source_receipts=receipts,job_id=os.environ.get('SLURM_JOB_ID'),scientific_review='PENDING',outputs=[meta(p) for p in out.iterdir() if p.is_file() and p.name!='COMPLETE.json']))
print('TERMINAL COHORT SUMMARY PASS',co,len(design),flush=True)
