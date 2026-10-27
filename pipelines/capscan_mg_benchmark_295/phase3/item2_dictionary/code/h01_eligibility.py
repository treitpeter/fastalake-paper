#!/usr/bin/env python3
"""Reproduce saved leader summaries and explicitly bounded eligibility sensitivities."""
import hashlib,json,sys
from pathlib import Path
import numpy as np
import pandas as pd

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8<<20),b''): h.update(b)
    return h.hexdigest()

def read(p): return pd.read_csv(p,sep='\t',dtype={'participant':str},keep_default_na=False)
def mean(pairs):
    if pairs.empty:return None
    return float(pairs.groupby(['participant','region']).same_leader.mean().groupby('participant').mean().mean())

def run(run,stage):
    m=json.loads((run/'manifest.json').read_text())
    for x in m['sources']:assert sha(x['path'])==x['sha256'],x['path']
    out=run/stage;out.mkdir(exist_ok=False)
    base=Path(m['data'])
    c=read(base/'contributor_cells.tsv.gz');s=read(base/'contributor_region_cells.tsv.gz')
    archived=read(base/'contributor_agreement_cells.tsv.gz')
    assert not c.duplicated(['KO','participant']).any()
    assert not s.duplicated(['KO','participant','region']).any()
    eligible=set(c.groupby('KO').participant.nunique().loc[lambda x:x>=12].index)
    assert len(eligible)==1365 and len(c)==19564
    stool=s[s.region.eq('stool')].drop(columns='region')
    p=s[s.region.ne('stool') & s.KO.isin(eligible)].merge(stool,on=['KO','participant'],suffixes=('_capsule','_stool'),validate='m:1')
    p['same_leader']=p.leading_species_capsule.eq(p.leading_species_stool)
    p['well_resolved']=p[['species_resolved_fraction_capsule','species_resolved_fraction_stool']].min(axis=1).ge(.2)&p[['gap_capsule','gap_stool']].min(axis=1).ge(.05)
    keys=['KO','participant','region']
    pp=p.sort_values(keys).reset_index(drop=True);aa=archived.sort_values(keys).reset_index(drop=True)
    assert pp[keys].equals(aa[keys])
    for col in ['same_leader','well_resolved']:
        av=aa[col] if aa[col].dtype==bool else aa[col].astype(str).str.lower().eq('true')
        assert np.array_equal(pp[col],av),col
    strict=p[p.well_resolved].copy()
    original={'primary':mean(p),'strict_pair_filter':mean(strict)}
    assert abs(original['primary']-.4728904409733759)<1e-12
    assert round(100*original['strict_pair_filter'],1)==53.2
    summary={'status':'PASS','stage':stage,'source_rule_reproduced':True,'original':original,'primary_kos':len(eligible),'primary_cells':len(c),'primary_pairs':len(p),'strict_pairs':len(strict)}
    if stage=='pilot':
        # Independent arithmetic fixture distinguishes cell pooling from donor balance.
        toy=pd.DataFrame([('a','colon',True),('a','colon',True),('b','colon',False)],columns=['participant','region','same_leader'])
        assert mean(toy)==.5 and toy.same_leader.mean()==2/3
    else:
        # Two explicit interpretations, neither claimed to be the author's intended rule.
        # A: donor-pooled primary cells must now satisfy the same quality thresholds.
        strict_c=c[c.species_resolved_fraction.ge(.2)&c.gap.ge(.05)]
        keepA=set(strict_c.groupby('KO').participant.nunique().loc[lambda x:x>=12].index)
        # B: require twelve donors with at least one quality-passing capsule/stool pair.
        keepB=set(strict.groupby('KO').participant.nunique().loc[lambda x:x>=12].index)
        cohorts={'primary':p,'strict_pair_filter':strict,'strict_donor_cell_eligibility':strict[strict.KO.isin(keepA)],'strict_paired_donor_eligibility':strict[strict.KO.isin(keepB)]}
        rows=[];donors=[]
        for rule,f in cohorts.items():
            d=f.groupby(['participant','region']).agg(agreement=('same_leader','mean'),n_pairs=('same_leader','size')).reset_index()
            d['rule']=rule;donors.append(d)
            rows.append(dict(rule=rule,donor_balanced_agreement=mean(f),pooled_cell_agreement=float(f.same_leader.mean()) if len(f) else None,n_pairs=len(f),n_kos=f.KO.nunique(),n_donors=f.participant.nunique(),n_donor_regions=len(d)))
        p['strict_donor_cell_eligible_KO']=p.KO.isin(keepA)
        p['strict_paired_donor_eligible_KO']=p.KO.isin(keepB)
        p.to_csv(out/'pair_accounting.tsv.gz',sep='\t',index=False)
        c.assign(strict_donor_cell=c.species_resolved_fraction.ge(.2)&c.gap.ge(.05)).to_csv(out/'donor_eligibility.tsv.gz',sep='\t',index=False)
        pd.DataFrame(rows).to_csv(out/'rule_comparison.tsv',sep='\t',index=False)
        pd.concat(donors).to_csv(out/'donor_region_comparison.tsv',sep='\t',index=False)
        summary['sensitivities']=rows
        (out/'REPORT.md').write_text('# Contributor eligibility audit\n\nThe historical 47.3% and 53.2% reproduce from archived region cells. The strict result filters primary-eligible pairs and does not reapply donor eligibility. Two explicit alternative eligibility interpretations are supplied as sensitivities, not replacements or claims about author intent. Original frozen LFQ and simplified grouping lineage; not corrected joint LFQ. Participant-mean agreement is primary; pooled cells are secondary.\n\n'+pd.DataFrame(rows).to_csv(sep='\t',index=False)+'\n')
    summary['outputs']=[{'path':str(x),'sha256':sha(x)} for x in sorted(out.iterdir())]
    (out/'COMPLETE.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary),flush=True)
if __name__=='__main__':run(Path(sys.argv[1]),sys.argv[2])
