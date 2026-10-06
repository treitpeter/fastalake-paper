"""fixed_database_evidence_ablation: plot corrected ablation. Moved from the original analysis as it ran, apart from path scrubbing; see README.md."""
from pathlib import Path
import json,hashlib,shutil
import pandas as pd,numpy as np
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
B=Path(__file__).resolve().parents[1];R=B.parents[1];O=B/'figures/H05_H09_evidence_contribution';O.mkdir(parents=True,exist_ok=True);(O/'source_data').mkdir(exist_ok=True)
assert json.loads((B/'evidence_ablation/full_PREDICT/INDEPENDENT_AUDIT.json').read_text())['status']=='PASS'
for name in ['Regular.ttf','Bold.ttf']:font_manager.fontManager.addfont(str(R/'MM_review_analysis/v7_fusion_2026-09-15/render_tools/fonts'/name))
plt.rcParams.update({'font.family':'FastaLakeSoehne','font.size':7,'axes.labelsize':7,'axes.titlesize':8,'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
ds=[]
for co,p,n in [('PREDICT',B/'evidence_ablation/full_PREDICT/acquisition_summary.tsv',296),('CapScan',R/'Hartout_edits/H09_evidence_ablation/full_CAPSCAN/acquisition_summary.tsv',282)]:
 d=pd.read_csv(p,sep='\t');assert len(d)==n and d['sample'].is_unique;d['cohort']=co;d['intensity_removed_percent']=100*d.fraction_intensity_removed;d['KO_support_lost_percent']=100*(1-d.shared_KOs/d.all_KOs);d['BrayCurtis']=d.normalized_KO_profile_BrayCurtis_all_vs_shared;ds.append(d);shutil.copy2(p,O/'source_data'/(co+'_acquisition_summary.tsv'))
d=pd.concat(ds);d.to_csv(O/'source_data/plot_values.tsv',sep='\t',index=False)
fig,axes=plt.subplots(1,3,figsize=(178/25.4,70/25.4));fig.subplots_adjust(left=.075,right=.99,bottom=.18,top=.84,wspace=.42)
for ai,(ax,col,title) in enumerate(zip(axes,['intensity_removed_percent','KO_support_lost_percent','BrayCurtis'],['Intensity removed (%)','KEGG orthologues losing support (%)','KO-profile Bray–Curtis distance'])):
 for i,(co,color) in enumerate([('PREDICT','#28688f'),('CapScan','#4d8362')]):
  v=d[d.cohort.eq(co)][col].to_numpy();x=i+np.random.default_rng(20260915+i).uniform(-.20,.20,len(v));ax.scatter(x,v,s=4,color=color,alpha=.40,edgecolors='none');q=np.quantile(v,[.25,.5,.75]);ax.plot([i,i],q[[0,2]],color='black',lw=1.4);ax.plot([i-.12,i+.12],[q[1],q[1]],color='black',lw=1.4)
 ax.set_xticks([0,1],['PREDICT\n296 acquisitions','CapScan\n282 acquisitions']);ax.set_xlim(-.5,1.5);ax.set_ylim(bottom=0);ax.set_title(title,fontsize=7);ax.text(-.18,1.14,'abc'[ai],transform=ax.transAxes,fontweight='bold',fontsize=10)
for ext in ['pdf','svg','png']:fig.savefig(O/('H05_H09_evidence_contribution.'+ext),dpi=300)
fig.savefig(O/'REVIEW.png',dpi=180);plt.close(fig)
(O/'CAPTION.md').write_text('Contribution of additional peptide evidence to measured functional profiles. a–c, Every matched acquisition is shown for positive label-free quantification intensity removed, detected Kyoto Encyclopedia of Genes and Genomes orthologues (KOs) losing all support, and normalised KO-profile Bray–Curtis distance after removing features whose peptide was not identified by the same-acquisition metagenome search. Black marks show medians and interquartile ranges. PREDICT uses the corrected 306-acquisition study dictionary restricted to 296 metagenome-matched acquisitions; CapScan retains its reviewed 282-acquisition dictionary. The dictionary and annotations are fixed within each all-evidence/shared-only contrast. KO intensity is split equally among assigned terms, with unannotated intensity retained in the total denominator. Median intensity removed is 23.89% and 30.40%; median KO support lost is 10.81% and 15.52%, for PREDICT and CapScan, respectively. These are acquisition summaries, not independent patient estimates. No new hypothesis tests are shown. Absence from the metagenome identification set does not imply absence from its database, correct sequence assignment or unique species origin. The shared-only profile is not the metagenome workflow profile; this ablation does not refit clinical models.\n')
(O/'VALIDATION.json').write_text(json.dumps(dict(status='NUMERICAL_RENDER_PASS_VISUAL_REVIEW_REQUIRED',code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),outputs=[dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in O.rglob('*') if p.is_file() and p.name!='VALIDATION.json']),indent=2)+'\n')
print('CORRECTED ABLATION FIGURE RENDERED',flush=True)
