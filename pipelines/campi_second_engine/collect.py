"""Combine validated acquisition outputs and independently verify displayed counts."""
import collections,csv,hashlib,itertools,json,math,os,sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(sys.argv[1]).resolve();OUT=R/'results'
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8<<20),b''):h.update(b)
 return h.hexdigest()
def dump(p,v):p.write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
def completion(p,names=None):
 c=json.loads((p/'COMPLETE.json').read_text());assert c['status']=='PASS'
 for x in c['outputs']:
  f=Path(x['path'])
  if names is None or f.name in names:assert sha(f)==x['sha256'],f
 return c
def readset(p):
 lines=p.read_text().splitlines();assert len(lines)==len(set(lines));assert all(x and x.isalpha() and x==x.upper() and 'I' not in x for x in lines);return set(lines)
def savefig(fig,name):
 fig.tight_layout()
 for ext in ['png','pdf','svg']:fig.savefig(OUT/(name+'.'+ext),dpi=180)
 plt.close(fig)
def finish(scope,files,counts):
 dump(OUT/'AGGREGATE_VALIDATION.json',dict(status='PASS',scope=scope,counts=counts,job_id=os.getenv('SLURM_JOB_ID'),source_manifest_sha256=sha(R/'inputs/manifest.json') if (R/'inputs/manifest.json').exists() else sha(R/'inputs/run_manifest.json'),outputs={p.name:sha(p) for p in files},scientific_review='pending'))
def support():
 m=json.loads((R/'inputs/manifest.json').read_text());assert len(m['tasks'])==710;rows=[];hist=[]
 for t in m['tasks']:
  p=OUT/t['cohort']/t['sample'];c=completion(p,{'protein_support.tsv'});x=pd.read_csv(p/'protein_support.tsv',sep='\t');assert len(x)==c['populations']['all_selected']['n'];assert not x.protein.duplicated().any();assert (x.production_prediction_sequences>0).all() and (x.parsimony_prediction_sequences>0).all()
  for pop,xx in [('all_selected',x),('identified_accessions',x[x.identified_accession==1]),('native_protein_q01',x[x.native_protein_q01==1])]:
   record=dict(cohort=t['cohort'],sample=t['sample'],population=pop,n=len(xx),one_top_prediction=int((xx.production_prediction_sequences==1).sum()),one_parsimony_prediction=int((xx.parsimony_prediction_sequences==1).sum()),with_FL_corroborated_prediction=int((xx.predictions_corroborated_by_FL_peptides>0).sum()))
   if t['cohort']=='CAMPI':record.update(with_REF_correct_top=int((xx.same_spectrum_REF_correct_top>0).sum()),with_REF_correct_all=int((xx.same_spectrum_REF_correct_all>0).sum()))
   assert record['n']==c['populations'][pop]['n']
   for key in c['populations'][pop]:assert record[key]==c['populations'][pop][key],(t['sample'],pop,key)
   for key in ['production_prediction_sequences','parsimony_prediction_sequences','production_prediction_spectra','parsimony_prediction_spectra','unique_among_selected_proteins','shared_among_selected_proteins']:
    for n,v in xx[key].value_counts().items():hist.append(dict(cohort=t['cohort'],sample=t['sample'],population=pop,support_definition=key,support_count=int(n),protein_acquisition_entries=int(v)))
    record[key+'_median']=float(xx[key].median()) if len(xx) else None
   for key in list(record):
    if key.startswith(('one_','with_')):record[key+'_fraction']=record[key]/len(xx) if len(xx) else None
   rows.append(record)
 d=pd.DataFrame(rows);h=pd.DataFrame(hist);d.to_csv(OUT/'acquisition_support_summary.tsv',sep='\t',index=False);h.to_csv(OUT/'support_histograms.tsv',sep='\t',index=False)
 cols=[c for c in d if c.endswith('_fraction')];g=d.groupby(['cohort','population'])[cols].agg(['median','min','max']);g.to_csv(OUT/'cohort_support_summary.tsv',sep='\t')
 fig,ax=plt.subplots(1,2,figsize=(10,4));a=d.query('population=="identified_accessions"');cohorts=['CAMPI','HSTOOL','PREDICT','CAPSCAN']
 ax[0].boxplot([100*a[a.cohort==c].one_top_prediction_fraction.dropna() for c in cohorts],labels=cohorts,showfliers=False);ax[0].set(ylabel='Identified accessions with one lookup prediction (%)',title='Support in the production top 30%')
 z=a[a.cohort=='CAMPI'];ax[1].plot(range(len(z)),100*z.with_REF_correct_top_fraction,'o',label='At least one REF-agreeing lookup prediction');ax[1].plot(range(len(z)),100*z.with_REF_correct_all_fraction,'s',label='At least one REF-agreeing parsimony prediction');ax[1].set(xlabel='CAMPI acquisition',ylabel='Identified accessions (%)');ax[1].legend(frameon=False,fontsize=8);savefig(fig,'prediction_support')
 lines=['# Prediction support — all 710 primary acquisitions','', 'Every selected accession has reconstructed lookup and parsimony evidence. The histograms distinguish distinct predicted sequences from supporting spectra, and selected entries from accessions assigned accepted peptides and the native protein-q subset. Counts across acquisitions are repeated observations of accessions, not study-wide distinct proteins.','', 'The principal question is how many identified accessions have at least one independently corroborated prediction. CAMPI supplies same-spectrum agreement with the accepted REF search; this is conditional reference agreement, not whole-model benchmark recall. The other cohorts supply agreement with the downstream FastaLake search, which is supporting evidence with shared data and selection, not independent truth.','', '| Cohort | Acquisitions | Median fraction of identified accessions with one lookup prediction | Median with one parsimony prediction |','|---|---:|---:|---:|']
 for cohort,x in a.groupby('cohort'):lines.append(f'| {cohort} | {len(x)} | {100*x.one_top_prediction_fraction.median():.2f}% | {100*x.one_parsimony_prediction_fraction.median():.2f}% |')
 lines+=['',f'Across the 15 CAMPI acquisitions, the median fraction of identified accessions with at least one same-spectrum REF-agreeing production lookup prediction is {100*z.with_REF_correct_top_fraction.median():.2f}%; using all qualifying parsimony predictions it is {100*z.with_REF_correct_all_fraction.median():.2f}%.','', 'These counts describe available supporting predictions. They do not establish the minimal number causally required by a rerun of protein selection. One shared prediction can support several homologous proteins, and an accession with one supporting prediction is not necessarily species-specific. The af3ac89b benchmark provenance audit remains a separate requirement.']
 (OUT/'REPORT.md').write_text('\n'.join(lines)+'\n');finish('Acquisition support counts independently reconciled with protein tables',list(OUT.glob('*.tsv'))+[OUT/'REPORT.md'],{'acquisitions':710})
def overlaps():
 m=json.loads((R/'inputs/manifest.json').read_text());assert len(m['tasks'])==652;rows=[];curves=[];params=[]
 for t in m['tasks']:
  p=OUT/t['cohort']/t['sample'];names={a+'_'+str(lo)+'_50.txt' for a in t['arms'] for lo in [5,9]}|{'q_length_counts.tsv','pairwise_overlaps.tsv','parameter_comparison.json'};c=completion(p,names);sets={(a,lo):readset(p/(a+'_'+str(lo)+'_50.txt')) for a in t['arms'] for lo in [5,9]}
  for a in t['arms']:assert sets[a,9]<=(sets[a,5])
  for ref in c['comparisons']:
   a,b,lo=ref['arm_a'],ref['arm_b'],ref['min_length'];x,y=sets[a,lo],sets[b,lo]
   computed=dict(n_a=len(x),n_b=len(y),intersection=len(x&y),a_only=len(x-y),b_only=len(y-x))
   for k,v in computed.items():assert v==ref[k]
   rows.append(dict(cohort=t['cohort'],sample=t['sample'],**ref))
  q=pd.read_csv(p/'q_length_counts.tsv',sep='\t')
  for _,row in q[q.q_threshold==.01].iterrows():assert len(sets[row.arm,int(row.min_length)])==row.target_canonical
  q['cohort']=t['cohort'];q['sample']=t['sample'];curves.append(q);pc=json.loads((p/'parameter_comparison.json').read_text());params.append(dict(cohort=t['cohort'],sample=t['sample'],identical_search_parameters=pc['identical_search_parameters'],versions=json.dumps(pc['versions'],sort_keys=True)))
 d=pd.DataFrame(rows);d.to_csv(OUT/'all_pairwise_overlaps.tsv',sep='\t',index=False);pd.concat(curves).to_csv(OUT/'all_q_length_counts.tsv',sep='\t',index=False);pd.DataFrame(params).to_csv(OUT/'search_parameter_agreement.tsv',sep='\t',index=False)
 summary=[]
 for k,x in d.groupby(['cohort','arm_a','arm_b','min_length']):
  z=dict(zip(['cohort','arm_a','arm_b','min_length'],k));z['acquisitions']=len(x);z['more_peptides_arm_a']=int((x.n_a>x.n_b).sum());z['same_peptides_count']=int((x.n_a==x.n_b).sum());z['fewer_peptides_arm_a']=int((x.n_a<x.n_b).sum())
  for field in ['n_a','n_b','intersection','a_only','b_only','a_recovery','b_recovery','ratio_a_to_b']:z.update({field+'_median':x[field].median(),field+'_min':x[field].min(),field+'_max':x[field].max()})
  summary.append(z)
 s=pd.DataFrame(summary);s.to_csv(OUT/'cohort_comparison_summary.tsv',sep='\t',index=False)
 for cohort,expected in [('PREDICT',296),('CAPSCAN',282)]:assert len(d.query('cohort==@cohort and arm_a=="FL_cohortMG" and arm_b=="MG" and min_length==5'))==expected
 fig,ax=plt.subplots(1,2,figsize=(10,4))
 for j,cohort in enumerate(['PREDICT','CAPSCAN']):
  x=d.query('cohort==@cohort and arm_a=="FL_cohortMG" and arm_b=="MG" and min_length==5');ax[j].scatter(x.n_b,x.n_a,s=9,alpha=.5);lim=max(x.n_a.max(),x.n_b.max());ax[j].plot([0,lim],[0,lim],':',color='grey');ax[j].set(xlabel='Specimen-metagenome peptides',ylabel='FastaLake peptides',title=cohort+' — '+str(len(x))+' matched acquisitions')
 savefig(fig,'matched_metagenome_yield')
 lines=['# Identity recovery — 652 acquisition recounts','', 'The full five-arm CAMPI comparison, PREDICT public-only and cohort-metagenome comparisons, and CapScan cohort-metagenome comparison have been recounted. Identity intersections and both directions of disagreement are retained, alongside peptide-count ratios. A larger list does not imply that every comparator peptide was recovered.','', '| Cohort / comparison | Matched acquisitions | More / same / fewer peptides | Median comparator identities recovered | Median count ratio |','|---|---:|---:|---:|---:|']
 for _,x in s[s.min_length==5].iterrows():lines.append(f'| {x.cohort}: {x.arm_a} / {x.arm_b} | {int(x.acquisitions)} | {int(x.more_peptides_arm_a)} / {int(x.same_peptides_count)} / {int(x.fewer_peptides_arm_a)} | {100*x.b_recovery_median:.2f}% | {x.ratio_a_to_b_median:.3f} |')
 lines+=['', 'Acceptance: rank-one target PSM rows, peptide q ≤ 0.01, modifications and I/L collapsed to canonical sequences. Both 5–50 and nested 9–50 populations are provided; the five q thresholds are not independently recalibrated confidence claims. All arms are tied to the same physical acquisition. Per-acquisition parameter differences are retained and must be considered in interpretation.','', 'The 652 units are 15 CAMPI, 306 PREDICT and 331 CapScan acquisitions. Missing specimen metagenomes remain missing, yielding 296 and 282 matched comparisons. The newly selected public-only CapScan searches are a separate dependent report. CAMPI member-specific peptide assignment, CapScan donor-level interpretation, and public-only PREDICT source provenance remain required supplements to these count and identity results.']
 (OUT/'REPORT.md').write_text('\n'.join(lines)+'\n');finish('Exact identity sets and q=.01 counts reconciled',list(OUT.glob('*.tsv'))+[OUT/'REPORT.md'],{'acquisitions':652,'PREDICT_matched_MG':296,'CAPSCAN_matched_MG':282})
def engine():
 m=json.loads((R/'inputs/manifest.json').read_text());assert len(m['tasks'])==45;rows=[];sets={};a10=R.parent/'A10_4_3_4_5_identity_recovery/results/CAMPI'
 for t in m['tasks']:
  p=OUT/t['id'];c=completion(p,{'accepted_5_50.txt','accepted_9_50.txt','accepted_canonical_peptides.tsv','pin_peptide_identity_audit.json'});assert c['sample']==t['sample'] and c['arm']==t['arm'];sage_arm='FL' if t['arm']=='FL_UNIPROT' else t['arm'];completion(a10/t['sample'],{sage_arm+'_5_50.txt',sage_arm+'_9_50.txt'})
  for lo,key in [(5,'accepted_canonical_peptides'),(9,'accepted_length9_50')]:
   x=readset(p/('accepted_'+str(lo)+'_50.txt'));assert len(x)==c[key];y=readset(a10/t['sample']/(sage_arm+'_'+str(lo)+'_50.txt'));sets[t['sample'],t['arm'],lo]=x
   rows.append(dict(sample=t['sample'],arm=t['arm'],min_length=lo,MSFragger=len(x),Sage=len(y),intersection=len(x&y),MSFragger_only=len(x-y),Sage_only=len(y-x),Sage_identity_recovered=len(x&y)/len(y) if y else None))
 d=pd.DataFrame(rows);d.to_csv(OUT/'between_engine_identities.tsv',sep='\t',index=False);comp=[]
 for sample in sorted({t['sample'] for t in m['tasks']}):
  for lo in [5,9]:
   x=sets[sample,'FL_UNIPROT',lo]
   for arm in ['REF','MG']:
    y=sets[sample,arm,lo];comp.append(dict(sample=sample,comparator=arm,min_length=lo,FL=len(x),reference=len(y),intersection=len(x&y),FL_only=len(x-y),reference_only=len(y-x),reference_recovery=len(x&y)/len(y) if y else None,count_ratio=len(x)/len(y) if y else None))
 c=pd.DataFrame(comp);c.to_csv(OUT/'within_MSFragger_comparisons.tsv',sep='\t',index=False);fig,ax=plt.subplots(1,2,figsize=(10,4))
 for j,arm in enumerate(['REF','MG']):
  z=c.query('comparator==@arm and min_length==5');ax[j].plot(z.reference,z.FL,'o');lim=max(z.reference.max(),z.FL.max());ax[j].plot([0,lim],[0,lim],':',color='grey');ax[j].set(xlabel=arm+' peptides',ylabel='FastaLake peptides',title='MSFragger + Percolator')
 savefig(fig,'second_engine')
 lines=['# Independent search engine — 45 completed searches','', 'All 15 CAMPI acquisitions were searched with the actual selected FastaLake FASTA, deposited REF, and metagenome FASTA using MSFragger 4.1 and Percolator 3.06.5. The comparison within the second engine is the primary test of portability; agreement with Sage is an additional identity analysis.','', '| Comparator | Acquisitions where FL gives more peptides | Median count ratio | Median comparator identity recovery |','|---|---:|---:|---:|']
 for arm,z in c.query('min_length==5').groupby('comparator'):lines.append(f'| {arm} | {int((z.FL>z.reference).sum())}/15 | {z.count_ratio.median():.3f} | {100*z.reference_recovery.median():.2f}% |')
 lines+=['', 'Peptide confidence is computed independently per acquisition and engine at q ≤ 0.01. Native MSFragger PIN precursor-charge suffixes were checked against charge features and removed from the peptide identity before Percolator peptide competition. Decoy labels, scores and other features are retained. Exact target FASTAs and input spectra are recorded by hash. Reverse-protein decoy generation is explicit.','', 'Search chemistry and main mass/digestion settings follow the source analyses, but engine preprocessing, scoring and confidence estimation remain different. A portable FASTA does not establish identical results or whole-procedure statistical calibration. Both 5–50 and nested 9–50 canonical peptide populations are supplied. Community member-specific interpretation remains in the CAMPI identity analysis.']
 (OUT/'REPORT.md').write_text('\n'.join(lines)+'\n');finish('Second-engine sets reconciled with accepted counts; exact cross-engine intersections',list(OUT.glob('*.tsv'))+[OUT/'REPORT.md'],{'searches':45})
if __name__=='__main__':{'support':support,'overlaps':overlaps,'engine':engine}[sys.argv[2]]()
