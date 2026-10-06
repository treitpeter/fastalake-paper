import pandas as pd,numpy as np,glob,json,sys
from pathlib import Path
W=Path(__file__).resolve().parents[1];V=W/'verify'
def load(mode):return pd.concat([pd.read_csv(p,sep='\t',dtype={'pair_bitscore':float}) for p in sorted(glob.glob(str(V/f'prep_{mode}_out/shard_*.tsv.gz')))],ignore_index=True)
d=load('real');d['ratio']=d.pair_bitscore.fillna(0)/d.native_score
res=dict(rows=len(d),genomes=int(d.genome.nunique()),status=d.status.value_counts().to_dict(),
  ratio_quantiles={str(q):round(float(d.ratio.quantile(q)),4) for q in [.001,.005,.01,.02,.03,.05,.1,.5,.9]})
for t in [0.2,0.5,0.8]:res[f'fail_at_{t}']=dict(n=int((d.ratio<t).sum()),fraction=float((d.ratio<t).mean()))
# impossibility check: reported bitscore cannot exceed what the attached sequence could reach (about 2.2 bits per residue upper bound for BLOSUM62 self-score is ~5 bits/res; use self-agnostic: ratio)
d['num']=d.accession.str.split('_').str[1].astype(int)
seed=dict(zip(d.accession,d.seed_ortholog))
nb=[(seed.get(f"{g}_{n-1:05d}")==s) or (seed.get(f"{g}_{n+1:05d}")==s) for g,n,s in zip(d.genome,d.num,d.seed_ortholog)]
d['seed_equals_adjacent_gene_seed']=nb
f=d.ratio<0.5
res['fail_seed_equals_adjacent']=float(d.loc[f,'seed_equals_adjacent_gene_seed'].mean());res['pass_seed_equals_adjacent']=float(d.loc[~f,'seed_equals_adjacent_gene_seed'].mean())
res['fail_query_shorter_than_half_seed']=float((d.loc[f,'query_length']<0.5*d.loc[f,'seed_length']).mean());res['pass_query_shorter_than_half_seed']=float((d.loc[~f,'query_length']<0.5*d.loc[~f,'seed_length']).mean())
pg=d.groupby('genome').apply(lambda g:(g.ratio<0.5).mean());res['per_genome_fail_fraction_quantiles']={str(q):round(float(pg.quantile(q)),4) for q in [.1,.5,.9,.99,1]}
res['genomes_fail_over_20pct']=int((pg>0.2).sum())
ctl={a:d.set_index('accession').loc[a,['seed_ortholog','native_score','pair_bitscore','ratio']].to_dict() for a in ['MGYG000000644_00789','MGYG000000644_00788','MGYG000000543_00549']}
res['positive_controls']=ctl
c=load('control');c['ratio']=c.pair_bitscore.fillna(0)/c.native_score
res['planted_negative_control']=dict(rows=len(c),shards='0-3',pass_at_0_5=int((c.ratio>=0.5).sum()),pass_fraction_at_0_5=float((c.ratio>=0.5).mean()))
# verification table for the projection code
d[['genome','accession','seed_ortholog','native_score','pair_bitscore','ratio','status','seed_equals_adjacent_gene_seed','query_length','seed_length','best_seed_in_genome_set','best_seed_bitscore']].to_csv(V/'native_row_verification.tsv.gz',sep='\t',index=False)
json.dump(res,open(V/'VERIFICATION_SUMMARY.json','w'),indent=1,default=str);print(json.dumps(res,indent=1,default=str))
