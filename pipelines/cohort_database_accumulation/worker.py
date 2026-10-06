from pathlib import Path
import sys,json,hashlib,csv,os,gc,re,collections
import numpy as np
H=Path(__file__).resolve().parent
sys.path.insert(0,str(H.parent/'confidence'))
from peptide import canonical
D=np.dtype([('seq','S32'),('header','S32'),('sid','u2')]);DR=np.dtype(D.descr+[('local','S32')])
def dump(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(8<<20),b''):h.update(b)
 return h.hexdigest()
def fasta(path,sid,raw=False):
 stat=Path(path).stat();h=hashlib.sha256();rows=[];acc=None;chunks=[]
 def add():
  if acc is None:return
  seq=b''.join(chunks);assert seq and acc
  sh=hashlib.sha256(seq).digest();hh=hashlib.sha256(acc).digest()
  row=(sh,hh,sid)
  if raw:
   i=acc.find(b'k141_');local=acc[i:] if i>=0 else acc;row+=(hashlib.sha256(local).digest(),)
  rows.append(row)
 with open(path,'rb') as f:
  for line in f:
   h.update(line)
   if line.startswith(b'>'):add();acc=line[1:].split()[0];chunks=[]
   else:chunks.append(line.strip())
 add();st=Path(path).stat();assert (stat.st_size,stat.st_mtime_ns)==(st.st_size,st.st_mtime_ns)
 arr=np.array(rows,dtype=DR if raw else D);assert len(arr)>0
 return arr,{'path':path,'sha256':h.hexdigest(),'bytes':st.st_size,'mtime_ns':st.st_mtime_ns,'records':len(arr)}
def peptides(info):
 path=info['results'];st=Path(path).stat();assert (st.st_size,st.st_mtime_ns)==(info['results_bytes'],info['results_mtime_ns'])
 out=set();h=hashlib.sha256();rows=0
 with open(path,'rb') as f:
  line=f.readline();h.update(line);names=line.decode().rstrip('\r\n').split('\t');ix={s:i for i,s in enumerate(names)}
  for line in f:
   h.update(line);v=line.decode().rstrip('\r\n').split('\t');assert len(v)==len(names);rows+=1
   if v[ix['rank']]!='1' or v[ix['label']]!='1' or float(v[ix['peptide_q']])>.01:continue
   p=canonical(v[ix['peptide']]);assert re.fullmatch('[A-Z]+',p)
   if 5<=len(p)<=50:out.add(p)
 for length in [5,9]:
  e=info['expected_sets'][str(length)];assert sha(e['path'])==e['sha256'];expected=set(Path(e['path']).read_text().splitlines());assert {p for p in out if len(p)>=length}==expected,(path,length)
 assert (st.st_size,st.st_mtime_ns)==(Path(path).stat().st_size,Path(path).stat().st_mtime_ns)
 return out,{'path':path,'sha256':h.hexdigest(),'bytes':st.st_size,'mtime_ns':st.st_mtime_ns,'reported_psm_rows':rows,'baseline_sets_reproduced':[5,9]}
def table(path,rows):
 with path.open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def incidence(arr,key):
 # One occurrence per distinct identity per specimen. Full fixed-width SHA256 bytes.
 pairs=np.empty(len(arr),dtype=[('key','S32'),('sid','u2')]);pairs['key']=arr[key];pairs['sid']=arr['sid'];pairs=np.unique(pairs)
 ids,starts,freq=np.unique(pairs['key'],return_index=True,return_counts=True)
 first=pairs['sid'][starts].copy();del pairs;return ids,freq,first
def audit(arr,key,label):
 ids,sfreq,_=incidence(arr,key);_,hfreq,_=incidence(arr,'header')
 pairs=np.empty(len(arr),dtype=[('key','S32'),('header','S32')]);pairs['key']=arr[key];pairs['header']=arr['header'];pairs=np.unique(pairs)
 seqkeys,alias=np.unique(pairs['key'],return_counts=True);assert np.array_equal(ids,seqkeys)
 _,collision=np.unique(pairs['header'],return_counts=True);assert len(collision)==len(hfreq)
 sharedseq=sfreq>=2;sharedhdr=hfreq>=2;na=int(np.count_nonzero(sharedseq&(alias>1)));nc=int(np.count_nonzero(sharedhdr&(collision>1)))
 return {'namespace':label,'distinct_sequences':len(ids),'sequences_shared_across_specimens':int(sharedseq.sum()),'shared_sequences_with_multiple_headers':na,'alias_rate_among_shared_sequences':na/int(sharedseq.sum()) if sharedseq.sum() else None,'alias_rate_among_all_sequences':na/len(ids),'distinct_headers':len(hfreq),'headers_shared_across_specimens':int(sharedhdr.sum()),'shared_headers_with_multiple_sequences':nc,'collision_rate_among_shared_headers':nc/int(sharedhdr.sum()) if sharedhdr.sum() else None,'collision_rate_among_all_headers':nc/len(hfreq)}
def curves(arr,key,n):
 ids,freq,first=incidence(arr,key);hist=np.bincount(freq,minlength=n+1);cum=np.cumsum(np.bincount(first,minlength=n));p=np.ones(n+1);expected=[]
 for k in range(1,n+1):
  p*=np.maximum(0,(n-np.arange(n+1)-k+1)/(n-k+1));expected.append(float(np.dot(hist,1-p)))
 assert cum[-1]==len(ids) and abs(expected[-1]-len(ids))<.001
 assert np.all(np.diff(cum)>=0) and np.all(np.diff(expected)>=-1e-6)
 return ids,freq,cum,expected

def main():
 cohort,mode=sys.argv[1:];manifest=json.loads((H/'manifest.json').read_text());rows=[r for r in manifest['rows'] if r['cohort']==cohort];specimens=sorted(set(r['design']['specimen_id'] for r in rows))
 if mode=='pilot':specimens=specimens[:3];rows=[r for r in rows if r['design']['specimen_id'] in specimens]
 assert mode in ['pilot','full'];out=H/mode/cohort;out.mkdir(parents=True,exist_ok=False);sids={s:i for i,s in enumerate(specimens)}
 arrays={a:[] for a in ['FL_public','MG']};rawarrays=[];rawseen=set();recs=[];acq=[];pools={a:[set() for _ in specimens] for a in arrays};checks=[]
 human_path=H.parents[3]/'MM_review_analysis/A01_4_1_database_sizes/inputs/human_swissprot_canonical_archived_2026_03.fasta'
 human,hr=fasta(str(human_path),0);humanseq=np.unique(human['seq']);del human
 for i,r in enumerate(rows):
  sid=sids[r['design']['specimen_id']]
  for arm,info in r['arms'].items():
   arr,rec=fasta(info['fasta'],sid);arrays[arm].append(arr);recs.append(rec)
   ids,pr=peptides(info);recs.append(pr);pools[arm][sid].update(ids)
   acq.append({'cohort':cohort,'sample':r['sample'],'specimen_id':r['design']['specimen_id'],'subject_id':r['design']['subject_id'],'group':r['design']['group'],'arm':arm,'fasta_records':len(arr),'distinct_sequences':len(np.unique(arr['seq'])),'peptide_ids':len(ids)})
   if arm=='MG':
    raw,rr=fasta(info['raw_fasta'],sid,True);recs.append(rr);actual=np.unique(arr['seq']);wanted=np.union1d(raw['seq'],humanseq)
    assert np.array_equal(actual,wanted),(r['sample'],'MG raw+human set mismatch');checks.append({'sample':r['sample'],'MG_original_plus_human_exact_set_equal':True})
    key=(sid,info['raw_fasta'])
    if key not in rawseen:rawarrays.append(raw);rawseen.add(key)
  print(cohort,mode,i+1,'/',len(rows),r['sample'],flush=True)
 table(out/'per_acquisition.tsv',acq);dump(out/'INPUT_RECEIPTS.json',recs+[hr]);dump(out/'MG_SEQUENCE_PRESERVATION.json',checks)
 curve=[];audits=[];specrows=[]
 for arm in ['FL_public','MG']:
  arr=np.concatenate(arrays.pop(arm));gc.collect();print('AGGREGATE',cohort,arm,len(arr),flush=True)
  keys,freq,cum,expected=curves(arr,'seq',len(specimens));np.savez_compressed(out/(arm+'_sequence_incidence.npz'),keys=keys,frequency=freq);del keys,freq
  for sid,s in enumerate(specimens):
   rr=[r for r in rows if r['design']['specimen_id']==s];curve.append({'cohort':cohort,'arm':arm,'k':sid+1,'specimen_id':s,'cumulative_sequences_fixed_order':int(cum[sid]),'mean_cumulative_sequences_random_order':expected[sid],'new_sequences_fixed_order':int(cum[sid]-(cum[sid-1] if sid else 0)),'peptide_ids_in_this_specimen':len(pools[arm][sid])})
   specrows.append({'cohort':cohort,'arm':arm,'specimen_id':s,'subject_id':rr[0]['design']['subject_id'],'acquisitions':len(rr),'peptide_ids':len(pools[arm][sid])})
  audits.append(audit(arr,'seq',arm+'_as_searched'));del arr;gc.collect()
 raw=np.concatenate(rawarrays);rawarrays.clear();gc.collect();print('RAW HEADER AUDIT',cohort,len(raw),flush=True)
 audits.append(audit(raw,'seq','MG_original_prefixed_headers'));raw['header']=raw['local'];audits.append(audit(raw,'seq','MG_assembly_local_headers'));del raw;gc.collect()
 table(out/'accumulation.tsv',curve);table(out/'per_specimen.tsv',specrows);table(out/'header_audit.tsv',audits)
 dump(out/'COMPLETE.json',{'status':'PASS','cohort':cohort,'mode':mode,'acquisitions':len(rows),'specimens':len(specimens),'participants':len(set(r['design']['subject_id'] for r in rows)),'code_sha256':sha(__file__),'manifest_sha256':sha(H/'manifest.json'),'job_id':os.getenv('SLURM_JOB_ID'),'all_original_MG_plus_human_sets_equal':True,'all_current_peptide_sets_reproduced':True,'header_rate_denominator':'shared across specimens; also all-distinct denominator provided','spectral_candidates':'not measured; definition awaits clarification'})
 print('COMPLETE',cohort,mode,flush=True)
if __name__=='__main__':main()
