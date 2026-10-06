from pathlib import Path
import json,hashlib,csv,math,collections
import pandas as pd,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
H=Path(__file__).resolve().parent;O=H/'evidence';O.mkdir(exist_ok=True)
# Source-count presentation requested by ML11: actual inventories, no new searches.
src=Path('MM_review_analysis/A01_4_1_database_sizes/results/acquisition_sizes.tsv')
t=pd.read_csv(src,sep='\t');assert len(t)==710 and (t.entries<20431).all()
cols=['cohort','sample','entries','reservoir_entries'];t[cols].to_csv(O/'database_sizes_710.tsv',sep='\t',index=False)
order=['CAMPI','HSTOOL','PREDICT','CapScan']
fig,ax=plt.subplots(figsize=(7.2,4.2))
rng=np.random.default_rng(20260928)
for i,c in enumerate(order):
 z=t[t.cohort==c];ax.scatter(i+rng.uniform(-.14,.14,len(z)),z.entries,s=8,alpha=.35,color='#1675a9',linewidths=0)
 ax.plot([i-.2,i+.2],[z.entries.median()]*2,color='#093e5a',lw=2)
 pool=int(z.reservoir_entries.iloc[0]);assert z.reservoir_entries.nunique()==1
 ax.scatter([i],[pool],marker='D',s=40,color='#c67621',zorder=4)
 ax.text(i,pool*1.25,f'{pool/1e6:.1f} M',ha='center',fontsize=8)
ax.axhline(20431,ls='--',c='#555',lw=1,label='Human reference: 20,431 entries')
ax.set_yscale('log');ax.set_ylim(150,1.2e9);ax.set_xticks(range(4),['CAMPI\nn=15','Healthy stool\nn=58','PREDICT\nn=306','CapScan\nn=331'])
ax.set_ylabel('Protein entries (log scale)');ax.spines[['top','right']].set_visible(False)
ax.scatter([],[],s=20,color='#1675a9',label='Selected acquisition FASTAs (median line)')
ax.scatter([],[],s=35,color='#c67621',marker='D',label='Full cohort reservoir')
ax.legend(frameon=False,fontsize=8,loc='center left')
fig.tight_layout();fig.savefig(O/'Database_size_context.pdf');fig.savefig(O/'Database_size_context.png',dpi=180);fig.savefig(O/'Database_size_context.svg');plt.close(fig)
(O/'DATABASE_SIZE_NOTE.md').write_text('All 710 selected FASTAs and actual full reservoirs are shown on one log scale. PREDICT/CapScan reservoir diamonds include cohort assemblies; the public-only benchmarking reservoir is 233,497,702 entries. CAMPI uses the broad-UniProt reservoir. This is sequence-entry count, not a measured peptide-index/search-cost axis. An external UHGP90 release is not substituted for the actual reservoirs; adding that optional comparator requires a pinned version/count.')
# Complete terminal-control report from every receipt, not array status alone.
T=H/'analyses/terminal_scramble'
m=json.loads((T/'search_manifest.json').read_text());rows=[]
for i,task in enumerate(m['tasks']):
 p=T/'end_to_end'/str(i)/'COMPLETE.json';d=json.loads(p.read_text());assert d['status']=='PASS' and d['arm']==task['arm']
 assert d['manifest_sha256']==hashlib.sha256((T/'search_manifest.json').read_bytes()).hexdigest()
 for q in d['accepted_counts']:rows.append(dict(index=i,arm=d['arm'],seed=d['seed'],selected_proteins=d['selected_proteins'],parsimony_proteins=d['parsimony_proteins'],**q))
e=pd.DataFrame(rows);e.to_csv(O/'scramble_end_to_end.tsv',sep='\t',index=False)
main=e[(e.gate=='q_only')&(e.q==.01)&(e.min_length==5)&(e.max_length==50)]
assert len(main)==11 and main.query("arm != 'authentic'").n_peptides.eq(0).all()
assert main.query("arm == 'authentic'").n_peptides.iloc[0]==29663
new=pd.read_csv(T/'terminal_preserved_by_length.tsv',sep='\t')
old=pd.read_csv(T/'archived_all_position_by_length.tsv',sep='\t')
old=old[old.arm.str.startswith('SIM')].copy();old['seed']=old.seed.astype(int)
joined=new.merge(old[['seed','length','distinct_queries','matching_queries','fraction_matched']],on=['seed','length'],suffixes=('_terminal','_all'),validate='one_to_one')
joined['difference_pp']=100*(joined.match_fraction-joined.fraction_matched)
joined.to_csv(O/'scramble_lookup_comparison.tsv',sep='\t',index=False)
j=joined.groupby('length').agg(terminal_min=('match_fraction','min'),terminal_median=('match_fraction','median'),terminal_max=('match_fraction','max'),all_min=('fraction_matched','min'),all_median=('fraction_matched','median'),all_max=('fraction_matched','max')).reset_index()
j.to_csv(O/'scramble_lookup_summary.tsv',sep='\t',index=False)
# Explicit adverse confidence cases and independent table arithmetic review.
C=H/'analyses/confidence';st=json.loads((C/'summary/STATUS.json').read_text());assert st['status']=='PASS'
a=pd.read_csv(C/'summary/ordinary_all_pairs.tsv',sep='\t')
assert len(a)==578*36 and a['index'].nunique()==578
assert ((a.n_fl==a.shared+a.fl_only)&(a.n_mg==a.shared+a.mg_only)).all()
f=pd.read_csv(C/'summary/fusion_all.tsv',sep='\t')
assert len(f)==30*36
pooled=f.groupby(['seed','gate','q','min_length','max_length'])[['n_canonical','paired_numerator','n_unpaired_accepted']].sum()
published=pd.read_csv(C/'summary/fusion_summary.tsv',sep='\t').set_index(['seed','gate','q','min_length','max_length'])
assert np.allclose(((pooled.paired_numerator+pooled.n_unpaired_accepted)/pooled.n_canonical).sort_index(),published.diagnostic.sort_index(),equal_nan=True)
sent=[]
for p in (C/'recovery').glob('*/COMPLETE.json'):
 for task in json.loads(p.read_text())['inputs']:
  for arm,d in task['inputs'].items():
   if d['counts'].get('pep_unset_rows'):sent.append(dict(index=task['index'],arm=arm,**d['counts']))
pd.DataFrame(sent).to_csv(O/'unset_pep_cases.tsv',sep='\t',index=False)
zero=a[(a.gate=='pep01')&(a.q==.01)&(a.min_length==5)&(a.max_length==50)&((a.n_fl==0)|(a.n_mg==0))]
zero.to_csv(O/'pep01_zero_yield_cases.tsv',sep='\t',index=False)
# Data-only scientific report preserves distinctions and does not promote thresholds.
report={'database_size_source_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'selected_databases':710,'scramble_complete_searches':11,'scramble_control_searches':10,'all_control_q01_counts_zero':True,'fresh_authentic_count':29663,'ordinary_pairs':578,'fusion_searches':30,'unset_pep_acquisitions':len(sent),'pep01_zero_yield_acquisitions':len(zero),'pooled_fusion_arithmetic_verified':True,'full_workflow_fdr_established':False}
(O/'VALIDATION.json').write_text(json.dumps(report,indent=2))
print(json.dumps(report,indent=2));print(j.head(5).to_string(index=False))

