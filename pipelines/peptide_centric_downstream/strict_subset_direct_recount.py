"""Recount CapScan leading-species agreement (Fig. 5f) under direct peptide assignment only.

The v30 manuscript removed the original protein grouping. The published Fig. 5f values
(45.6% all eligible, 60.0% strict) had been restricted to functions eligible under both
grouping rules and comparisons resolved under both. This script recomputes them with the
direct-assignment criteria alone (46.4% and 57.6%) and reproduces the published values as
controls.

Input: the run directory of ``fig5ef_run.py`` (argument 1), which holds
``results/contributor_method_cells.tsv`` and ``results/contributor_eligibility.tsv``.
Output (written to the current directory): ``recount.json`` and four
``*_per_participant.tsv`` tables.

Estimator: comparisons (function x participant x intestinal region against the same
participant's stool) are averaged within participant x region, then over regions within
participant, then over the 15 participants. Statistical unit: participant, n = 15.
Strict = species-resolved fraction >= 0.20 and leader margin >= 0.05 in both the
intestinal-region and the stool cell (columns ``strict_simplified`` / ``both_strict``).

Code as it ran (LAB_ROOT/.../final/v30/grouping_recount/recount.py, Slurm job 6074394);
only this docstring and the function docstrings were added.
"""
import sys, json
import pandas as pd, numpy as np
R=sys.argv[1]
p=pd.read_csv(R+'/results/contributor_method_cells.tsv',sep='\t')
el=pd.read_csv(R+'/results/contributor_eligibility.tsv',sep='\t')
simp=set(el.KO[el.simplified_eligible]); common=set(el.KO[el.common_eligible])
p['res_o']=p.state_original.eq('resolved'); p['res_s']=p.state_simplified.eq('resolved')
def run(pp, okfun, label):
    """Mean participant agreement of capsule and stool leading species.

    Parameters
    ----------
    pp : DataFrame
        Method cells restricted to the eligible functions.
    okfun : callable
        Row filter on the capsule-stool merged table (which comparisons count).
    label : str
        Printed label.

    Returns
    -------
    dict
        mean (participant unit, n = 15), participants, cells, functions and the
        per-participant table.
    """
    st=pp[pp.region.eq('stool')].drop(columns='region')
    cp=pp[pp.region.ne('stool')].merge(st,on=['KO','participant'],suffixes=('_capsule','_stool'),validate='m:1')
    cc=cp[okfun(cp)].copy()
    cc['agree']=cc.leading_species_simplified_capsule.eq(cc.leading_species_simplified_stool)
    dd=cc.groupby(['participant','region']).agree.mean().reset_index()
    dv=dd.groupby('participant').agree.mean()
    nf=cc.groupby('participant').KO.nunique()
    ncells=cc.groupby('participant').size()
    print(f"{label}: participants={len(dv)} mean={dv.mean():.6f} cells={len(cc)} functions_total={cc.KO.nunique()}")
    return dict(label=label,mean=float(dv.mean()),participants=int(len(dv)),cells=int(len(cc)),functions=int(cc.KO.nunique()),
                per_participant=pd.DataFrame(dict(agreement=dv,functions=nf,comparisons=ncells)).reset_index())
B=lambda s: (lambda c: c[s+'_capsule'].astype(bool)&c[s+'_stool'].astype(bool))
out={}
# controls: published definitions
out['ctrl_resolved_both']=run(p[p.KO.isin(common)], lambda c:(c.res_o_capsule&c.res_s_capsule&c.res_o_stool&c.res_s_stool),'CONTROL resolved, common functions, both rules')
out['ctrl_strict_both']=run(p[p.KO.isin(common)], B('both_strict'),'CONTROL strict, common functions, both rules')
# direct only
out['direct_resolved']=run(p[p.KO.isin(simp)], lambda c:(c.res_s_capsule&c.res_s_stool),'DIRECT resolved, direct-eligible functions')
out['direct_strict']=run(p[p.KO.isin(simp)], B('strict_simplified'),'DIRECT strict, direct-eligible functions')
for k,v in out.items():
    v['per_participant'].to_csv(f'{k}_per_participant.tsv',sep='\t',index=False)
    print(k); print(v['per_participant'].to_string(index=False))
json.dump({k:{a:b for a,b in v.items() if a!='per_participant'} for k,v in out.items()},open('recount.json','w'),indent=1)
