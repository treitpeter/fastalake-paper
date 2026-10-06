"""Extend the frozen CAMPI clustering comparison; never overwrite prior results."""
import csv,gzip,hashlib,json,os,re,subprocess,sys,time
from pathlib import Path
import select as unused_select
# Import local scientific selector without shadowing Python's select module.
import importlib.util
BASE=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('cluster_selection',BASE/'frozen/selection.py');selection=importlib.util.module_from_spec(spec);spec.loader.exec_module(selection)
sha=selection.sha

def peptides(path):
 out=set()
 with Path(path).open() as f:
  for r in csv.DictReader(f,delimiter='\t'):
   if r['label']=='1' and 0<=float(r['peptide_q'])<=.01:
    p=re.sub(r'\[[^\]]*\]','',r['peptide']).replace('I','L')
    if 5<=len(p)<=50:
     assert re.fullmatch('[A-Z]+',p)
     assert int(r.get('rank','1'))==1,'Source comparison eligibility has changed: rank>1 accepted'
     out.add(p)
 return out

def run(index):
 m=json.loads((BASE/'manifest.json').read_text());t=m['tasks'][index];out=BASE/'results'/t['sample']
 assert not out.exists(),str(out)
 for rec in m['frozen_code']+m['lineage']+t['search_sources']:
  p=Path(rec['path']);st=p.stat();assert st.st_size==rec['bytes'] and st.st_mtime_ns==rec['mtime_ns'],str(p)
  if 'sha256' in rec:assert sha(p)==rec['sha256'],str(p)
 assert subprocess.check_output([m['sage'],'--version'],text=True).strip()=='sage 0.14.6'
 # Reproduce selection, including pilot99% gate, before doing any new search.
 selection.run(BASE,index)
 producer=json.loads((out/'COMPLETE.json').read_text());assert producer['status']=='PASS' and producer['baseline_sequence_set_reproduced']
 prior=Path(t['prior_search']);receipt=json.loads((prior/'COMPLETE.json').read_text());assert receipt['status']=='PASS'
 for rec in receipt['outputs']:
  if rec['path'].endswith(('results.sage.tsv','accepted_canonical_5_50.txt','config.json')):assert sha(rec['path'])==rec['sha256']
 baseline=peptides(prior/'unclustered/results.sage.tsv');saved=set((prior/'unclustered/accepted_canonical_5_50.txt').read_text().splitlines());assert baseline==saved
 archive=peptides(t['archived_search'])
 if t['sample']=='S02':
  assert baseline-archive=={'LEKDLLDGMALLSDDSYSVPK'} and archive-baseline=={'KVTAPYKYPR'}
 else:assert baseline==archive
 reference=set(Path(t['reference_peptides']).read_text().splitlines())
 oldsets={a:peptides(prior/a/'results.sage.tsv') for a in ['cluster_097','cluster_099']}
 for arm,s in oldsets.items():assert s==set((prior/arm/'accepted_canonical_5_50.txt').read_text().splitlines())
 cfg0=json.loads(Path(t['config']).read_text());cfg0.pop('output_paths',None)
 basecfg=json.loads((prior/'unclustered/config.json').read_text());assert cfg0==basecfg,'Original source settings differ from fixed denominator'
 sets={'unclustered':baseline,**oldsets};rows=[]
 for threshold in m['new_thresholds']:
  arm='cluster_'+threshold.replace('.','');fasta=out/(arm+'_parsimony')/(t['sample']+'_razor.fasta');dest=out/'searches'/arm;dest.mkdir(parents=True)
  cfg=json.loads(json.dumps(cfg0));cfg['database']['fasta']=str(fasta);(dest/'config.json').write_text(json.dumps(cfg,indent=2)+'\n')
  argv=[m['sage'],str(dest/'config.json'),'-o',str(dest),'--write-pin','--disable-telemetry-i-dont-want-to-improve-sage'];started=time.time()
  with (dest/'search.log').open('w') as f:r=subprocess.run(argv,env=dict(os.environ,RAYON_NUM_THREADS='2'),stdout=f,stderr=subprocess.STDOUT)
  with (out/'SEARCH_COMMANDS.jsonl').open('a') as f:f.write(json.dumps({'argv':argv,'returncode':r.returncode,'seconds':time.time()-started})+'\n')
  assert r.returncode==0
  sets[arm]=peptides(dest/'results.sage.tsv');(dest/'accepted_canonical_5_50.txt').write_text(''.join(p+'\n' for p in sorted(sets[arm])))
 baseline_sequences=list(selection.fasta(t['final_fasta']));nb=len(baseline_sequences)
 import ahocorasick
 union=set.union(reference,*sets.values());auto=ahocorasick.Automaton()
 for pep in union:auto.add_word(pep,pep)
 auto.make_automaton();availability={}
 for arm,peps in sets.items():
  fasta=Path(t['final_fasta']) if arm=='unclustered' else (Path(t['prior_h04']) if arm in oldsets else out)/(arm+'_parsimony')/(t['sample']+'_razor.fasta')
  seqs=list(selection.fasta(fasta));found={pep for _,seq in seqs for _,pep in auto.iter(seq.replace('I','L'))};assert peps<=found,arm;availability[arm]=found
  loss=baseline-peps;gain=peps-baseline
  rows.append({'sample':t['sample'],'arm':arm,'identity':1.0 if arm=='unclustered' else float('0.'+arm.split('_')[1][1:]),'protein_entries':len(seqs),'protein_entry_ratio':len(seqs)/nb,'peptides':len(peps),'baseline_peptides':len(baseline),'peptide_count_ratio':len(peps)/len(baseline),'peptide_reduction_pct':100*(1-len(peps)/len(baseline)),'baseline_peptide_identities_retained':len(peps&baseline),'baseline_identity_recovery':len(peps&baseline)/len(baseline),'REF_peptides':len(reference),'REF_recovered':len(peps&reference),'REF_recovery':len(peps&reference)/len(reference),'lost':len(loss),'gained':len(gain),'lost_sequence_absent':len(loss-found),'lost_sequence_available':len(loss&found)})
 with (out/'identification_summary.tsv').open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
 with gzip.open(out/'peptide_fates.tsv.gz','wt') as f:
  w=csv.writer(f,delimiter='\t');w.writerow(['peptide','REF',*sets,*[a+'_available' for a in sets]])
  for p in sorted(union):w.writerow([p,int(p in reference),*[int(p in s) for s in sets.values()],*[int(p in availability[a]) for a in sets]])
 if t['sample']=='S02':
  baselines={'archive':archive,'fixed_contemporaneous':baseline}
  for label in ['repeat_1','repeat_2']:baselines[label]=peptides(Path(m['s02_diagnostic'])/label/'results.sage.tsv')
  sensitivity=[{'baseline':name,'arm':arm,'count_ratio':len(peps)/len(bs),'gained':len(peps-bs),'lost':len(bs-peps)} for name,bs in baselines.items() for arm,peps in sets.items()]
  (out/'S02_BASELINE_SENSITIVITY.json').write_text(json.dumps(sensitivity,indent=2)+'\n')
 outputs=[out/'identification_summary.tsv',out/'peptide_fates.tsv.gz',out/'SEARCH_COMMANDS.jsonl',*sorted((out/'searches').rglob('results.sage.tsv')),*sorted((out/'searches').rglob('accepted_canonical_5_50.txt'))]
 result={'status':'PASS','sample':t['sample'],'job_id':os.environ.get('SLURM_JOB_ID'),'new_searches':len(m['new_thresholds']),'manifest_sha256':sha(BASE/'manifest.json'),'baseline_fixed_before_extension':True,'baseline_archive_exact':t['sample']!='S02','s02_exception_preserved':t['sample']=='S02','pilot_99_selected_and_representative_sequences_exact':True if t['sample'] in m['pilot_samples'] else None,'outputs':[{'path':str(p),'sha256':sha(p)} for p in outputs]}
 (out/'SEARCH_COMPLETE.json').write_text(json.dumps(result,indent=2)+'\n');print(t['sample'],'PASS',flush=True)
if __name__=='__main__':
 m=json.loads((BASE/'manifest.json').read_text())
 if sys.argv[1]=='pilot':
  for index in m['pilot_indices']:run(index)
  (BASE/'PILOT_COMPLETE.json').write_text(json.dumps({'status':'PASS','samples':m['pilot_samples'],'job_id':os.environ.get('SLURM_JOB_ID'),'manifest_sha256':sha(BASE/'manifest.json')},indent=2)+'\n')
 else:
  assert json.loads((BASE/'PILOT_COMPLETE.json').read_text())['status']=='PASS';run(int(sys.argv[1]))
