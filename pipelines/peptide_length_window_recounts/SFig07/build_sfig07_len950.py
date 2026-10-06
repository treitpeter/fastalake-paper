"""Supplementary Figure 7 at 9-50 residues: panels a and b rebuilt, c and d unchanged.

Inputs
- work/hstool/identity_per_acquisition_len.tsv: hstool_recount.py (Slurm srun 6072769) from A21/A22
  accepted identities, 9-50 = accepted 5-50 identities with >= 9 residues. Replaces
  identity_per_acquisition.tsv for panel a (only the joint_fileq column is read from it).
- repo identity_comparison.tsv rows with min_length 9 (A21, upstream): the other panel a columns.
- work/clin/part_*.tsv: clin_recount.py (Slurm array 6072755) accepted peptides per clinical
  acquisition, own vs common FASTA, 5-50 and 9-50. Aggregated here exactly as summarize_v2.py
  (acquisitions without participant excluded; mean within specimen, then within participant;
  100*(common/original-1)). Pseudonyms from the repo table by matching the 5-50 participant means.
  Panel b y (directLFQ groups) is unchanged.
The figure is drawn by panels_len950.render (repo panels.py with MIN_LENGTH_PANEL_A = 9).
Panels c, d must equal the base pixel for pixel; otherwise only the a and b boxes are pasted into the base.
Usage: build_sfig07_len950.py OUTDIR"""
import os  # release: paths come from environment variables
import glob,json,sys,hashlib,shutil
from pathlib import Path
import numpy as np,pandas as pd
from PIL import Image
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE))
import panels_len950 as P
LAB=Path(os.environ['FASTALAKE_DATA_DIR'])
REPO=LAB/'release_prep/fastalake-paper_v28';SD=REPO/'source_data/SFig07'
BASE_PNG=REPO/'figures/SFig07/expected/SFig07_v28f.png'
E=LAB/'MM_review_analysis/end_of_day_2026-09-15'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=True);sd=out/'source_data_len950';sd.mkdir(exist_ok=True)
for f in SD.glob('*.tsv'):shutil.copy(f,sd/f.name)
# panel a: per-acquisition table at 9-50 (same columns as identity_per_acquisition.tsv)
L=pd.read_csv(HERE/'work/hstool/identity_per_acquisition_len.tsv',sep='\t')
L[L.min_length.eq(9)].drop(columns='min_length').to_csv(sd/'identity_per_acquisition.tsv',sep='\t',index=False)
# control: 5-50 recount equals the repo table
r5=L[L.min_length.eq(5)].drop(columns='min_length').reset_index(drop=True);repo=pd.read_csv(SD/'identity_per_acquisition.tsv',sep='\t')
m=repo.merge(r5,on=['arm','sample'],suffixes=('','_r'))
checks=[dict(check='hstool 5-50 recount equals repo identity_per_acquisition (common arms)',passed=bool(len(m)==4*58 and all((m[c]==m[c+'_r']).all() for c in ['original','test','overlap'])))]
# panel b: participant aggregation
c=pd.concat([pd.read_csv(f,sep='\t') for f in sorted(glob.glob(str(HERE/'work/clin/part_*.tsv')))]);assert len(c)==637 and c['sample'].is_unique
fz=json.loads((E/'clinical_expansion/inputs/FROZEN.json').read_text());man=json.loads((Path(fz['run'])/'manifest.json').read_text())
rows=[]
for co,cohort in [('PREDICT','PREDICT'),('CAPSCAN','CAPSCAN')]:
    tasks=[t for t in man['tasks'] if t['cohort']==co]
    md=pd.DataFrame([dict(sample=t['sample'],participant=str(t['design']['participant']),specimen=str(t['design']['specimen_id']),region=str(t['design']['region'])) for t in tasks])
    d=c[c.cohort.eq(co)].merge(md,on='sample',validate='one_to_one');assert len(d)=={'PREDICT':306,'CAPSCAN':331}[co]
    paired=pd.read_csv(E/'common_cache/comparison/results/full'/co/'paired_acquisition_comparison.tsv',sep='\t').set_index('sample')
    checks.append(dict(check=f'{co} 5-50 recount equals paired_acquisition_comparison',passed=bool((d.set_index('sample').original_5==paired.original_peptides.reindex(d['sample'])).all() and (d.set_index('sample').common_5==paired.common_peptides.reindex(d['sample'])).all())))
    dp=d[d.participant.ne('')];cols=['original_5','common_5','original_9','common_9']
    spec=dp.groupby(['participant','region','specimen'],sort=True)[cols].mean().reset_index();don=spec.groupby('participant',sort=True)[cols].mean().reset_index()
    don['cohort']=cohort;rows.append(don)
don=pd.concat(rows)
pat=pd.read_csv(SD/'common_search_participants.tsv',sep='\t')
key=lambda a,b:(round(a,6),round(b,6))
look={(r.cohort,)+key(r.original_5,r.common_5):r for r in don.itertuples()}
assert len(look)==len(don)==len(pat)
newx=[];o9=[];c9=[]
for r in pat.itertuples():
    q=look[(r.cohort,)+key(r.original_peptides,r.common_peptides)]
    o9.append(q.original_9);c9.append(q.common_9);newx.append(100*(q.common_9/q.original_9-1))
checks.append(dict(check='participant 5-50 means reproduce repo peptides_percent_difference',passed=bool(np.allclose([100*(look[(r.cohort,)+key(r.original_peptides,r.common_peptides)].common_5/look[(r.cohort,)+key(r.original_peptides,r.common_peptides)].original_5-1) for r in pat.itertuples()],pat.peptides_percent_difference,rtol=0,atol=1e-9))))
pat9=pat.copy();pat9['original_peptides']=o9;pat9['common_peptides']=c9;pat9['peptides_difference']=pat9.common_peptides-pat9.original_peptides;pat9['peptides_percent_difference']=newx
for col in ['gained_peptides','lost_peptides']:pat9[col]=np.nan   # not recounted at 9-50; not plotted
pat9.to_csv(sd/'common_search_participants.tsv',sep='\t',index=False)
print(json.dumps(checks,indent=1));assert all(x['passed'] for x in checks)
res=P.render(sd,out)
print(res['checks'],res['font'],res['axes_bbox_px'])
new=np.array(Image.open(out/'SFig07.png').convert('RGB')).astype(int);base=np.array(Image.open(BASE_PNG).convert('RGB')).astype(int)
assert new.shape==base.shape,(new.shape,base.shape)
diff=(new!=base).any(2);boxes=res['axes_bbox_px']
mask=np.zeros_like(diff)
for k in 'ab':l,t,r_,b=boxes[k];mask[t:b,l:r_]=True
rec=dict(base_png=str(BASE_PNG),base_png_sha256=sha(BASE_PNG),panel_boxes_px=boxes,differing_outside_a_b=int((diff&~mask).sum()),differing_inside_a_b=int((diff&mask).sum()),checks=checks+res['checks'])
if rec['differing_outside_a_b']:
    comp=base.copy();comp[mask]=new[mask];Image.fromarray(comp.astype(np.uint8)).save(out/'SFig07_pasted.png',dpi=(300,300));rec['pasted']=True
Image.fromarray(np.where(diff[...,None],[255,0,0],new).astype(np.uint8)).save(out/'SFig07_diff.png')
(out/'SFig07_receipt.json').write_text(json.dumps(rec,indent=1,default=str)+'\n');print(json.dumps({k:v for k,v in rec.items() if k!='checks'},indent=1,default=str))
