import os
"""Figure 2 panel h source: exact member-set sharing of pooled blind-FL peptides, derived from the reviewed A14 per-acquisition peptide_membership tables (arm FL, 15 acquisitions)."""
import hashlib,json
from pathlib import Path
import pandas as pd
import numpy as np
ROOT=Path(os.environ['FASTALAKE_DATA_DIR']);A14=next((ROOT/'MM_review_analysis').glob('A14_*'))/'results';OUT=Path(__file__).resolve().parent/'source_data'
rs=pd.read_csv(A14/'recovery_summary.tsv',sep='\t');expect=int(rs[rs.arm.eq('FL')&rs.min_length.eq(5)].pooled_peptides.iloc[0])
mem={};srcs={};acq=0
for d in sorted(p for p in A14.iterdir() if p.is_dir() and (p/'peptide_membership.tsv').exists()):
    f=d/'peptide_membership.tsv';srcs[str(f.relative_to(ROOT))]=hashlib.sha256(f.read_bytes()).hexdigest();acq+=1
    t=pd.read_csv(f,sep='\t',keep_default_na=False);t=t[t.arms.str.split(';').apply(lambda a:'FL' in a)]
    for pep,m in zip(t.peptide,t.members):
        assert 5<=len(pep)<=50
        if pep in mem:assert mem[pep]==m,(pep,mem[pep],m)
        else:mem[pep]=m
assert acq==15,acq
assert len(mem)==expect,(len(mem),expect)
lens=pd.DataFrame({'peptide':list(mem.keys()),'members':list(mem.values())});lens['length']=lens.peptide.str.len();lens['class']=np.where(lens.members.eq(''),'No member',np.where(lens.members.str.contains(';'),'Shared by members','One member'))
lt=lens.groupby(['class','length']).size().rename('n_peptides').reset_index();lt['pct_within_class']=100*lt.n_peptides/lt.groupby('class').n_peptides.transform('sum');lt.to_csv(OUT/'length_by_member_class_a14.tsv',sep='\t',index=False)
lensum=lens.groupby('class').length.agg(['size','median']).reset_index()
c=pd.Series(list(mem.values())).value_counts()
rows=[dict(members=k if k else '(none of the eight)',degree=0 if not k else k.count(';')+1,n_peptides=int(v)) for k,v in c.items()]
u=pd.DataFrame(rows).sort_values('n_peptides',ascending=False).reset_index(drop=True);u.insert(0,'rank',range(1,len(u)+1))
u['members']=u.members.str.replace(';','|')
u.to_csv(OUT/'upset_intersections_a14.tsv',sep='\t',index=False)
nz=u[u.degree.gt(0)];u['drawn']=u.degree.eq(1)|u.index.isin(nz[nz.degree.gt(1)].head(6).index);assert u.drawn.sum()==14 and u[u.degree.eq(1)].shape[0]==8;u.to_csv(OUT/'upset_intersections_a14.tsv',sep='\t',index=False);top=u[u.drawn]
info=dict(acquisitions=acq,pooled_fl_peptides=len(mem),nonempty_combinations=int(len(nz)),peptides_in_member_sets=int(nz.n_peptides.sum()),no_member=int(u[u.degree.eq(0)].n_peptides.sum()),drawn_peptides=int(top.n_peptides.sum()),length_classes=lensum.to_dict('records'),drawn_rule='all eight single-member sets and the six largest shared sets',sources=srcs)
json.dump(info,open(OUT/'upset_intersections_a14.json','w'),indent=2);print(json.dumps({k:v for k,v in info.items() if k!='sources'},indent=1));print(top.to_string())
