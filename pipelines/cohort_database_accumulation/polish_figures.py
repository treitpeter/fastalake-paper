"""Restyle the accumulation figures written by plot.py for review display. No numbers are computed."""
from pathlib import Path
import csv,json,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from docx import Document
from docx.shared import Inches
H=Path(__file__).resolve().parent;O=H/'figures'
fig,axs=plt.subplots(1,2,figsize=(13,5));labels=['FastaLake\nas searched','Metagenome\nas searched','Metagenome\noriginal full\nheaders','Metagenome\nassembly-local\nnames'];keys=['FL_public_as_searched','MG_as_searched','MG_original_prefixed_headers','MG_assembly_local_headers']
for ax,c in zip(axs,['PREDICT','CAPSCAN']):
 lookup={r['namespace']:r for r in csv.DictReader((H/'full'/c/'header_audit.tsv').open(),delimiter='\t')};x=np.arange(4)
 for j,(metric,label,color) in enumerate([('alias_rate_among_shared_sequences','Shared sequence has multiple headers','#009E73'),('collision_rate_among_shared_headers','Shared header has multiple sequences','#CC79A7')]):
  vals=[float(lookup[k][metric])*100 if lookup[k][metric] else np.nan for k in keys];bars=ax.bar(x+(j-.5)*.36,vals,width=.36,color=color,label=label)
  for b,v in zip(bars,vals):
   label='N/A' if not np.isfinite(v) else ('<0.01%' if 0<v<.01 else (f'{v:.1f}%' if v else '0%'))
   ax.text(b.get_x()+b.get_width()/2,(v if np.isfinite(v) else 0)+1.7,label,ha='center',fontsize=8)
 ax.set_xticks(x,labels,fontsize=9);ax.set_ylim(0,111);ax.set_ylabel('Inconsistency among shared identities (%)');ax.set_title('CapScan' if c=='CAPSCAN' else c);ax.grid(axis='y',alpha=.15)
fig.legend(*axs[0].get_legend_handles_labels(),loc='lower center',bbox_to_anchor=(.5,-.07),ncol=2,frameon=False,fontsize=9);fig.tight_layout()
for ext in ['pdf','svg','png']:fig.savefig(O/('header_inconsistency_audit.'+ext),dpi=220,bbox_inches='tight')
plt.close(fig)
caption=(O/'CAPTIONS.md').read_text().strip()+' N/A in the header panel means no header was shared across specimens, so the shared-header denominator is zero. The as-searched metagenome alias rates are small but nonzero: 53 shared sequences have multiple headers in each cohort; these are explicitly labelled <0.01%, not rounded to zero. Candidate measurements are now provided in the separate candidate-competition figure proposal.'
(O/'CAPTIONS.md').write_text(caption+'\n')
doc=Document();doc.add_heading('FastaLake — empirical database accumulation and header audit',0);doc.add_paragraph('Figure proposals for Peter’s review. Full recounts completed for both cohorts; candidate measurements are provided separately. Figures have not been inserted into the manuscript.');doc.add_paragraph(caption)
for name in ['sequence_accumulation_dual_axis_linear','sequence_accumulation_dual_axis_log','sequence_accumulation_order_averaged','header_inconsistency_audit']:
 doc.add_heading(name.replace('_',' '),1);doc.add_picture(str(O/(name+'.png')),width=Inches(6.4))
doc.save(O/'Sequence_accumulation_figure_proposals.docx')
# Give candidate panel styles an explicit legend; retain all extreme points.
p=H/'candidates/plot.py';s=p.read_text();s=s.replace("fig.tight_layout();files=[]", "fig.text(.5,.005,'Solid: scored candidates (left axis). Dotted: peptide IDs in that acquisition (right axis).',ha='center',fontsize=10)\nfig.tight_layout(rect=(0,.03,1,1));files=[]")
p.write_text(s)
report=json.loads((O/'COMPLETE.json').read_text());report['candidate_per_spectrum']='Completed separately: candidates/full and candidates/figures';report['files']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in O.iterdir() if p.suffix in ['.png','.pdf','.svg','.docx']};(O/'COMPLETE.json').write_text(json.dumps(report,indent=2)+'\n')
print('Header labels clarified; nonzero rates and undefined denominators preserved.')
