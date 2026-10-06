"""capscan_metagenome_yield_join: build figures v15. Moved from the original analysis as it ran, apart from path scrubbing; see README.md."""
from pathlib import Path
import json,hashlib,shutil,collections,re
import fitz,pandas as pd,numpy as np,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
W=Path(__file__).resolve().parent;O=W.parent;B=O/'release_candidate';P=W/'coauthor_package';F=P/'figures';S=P/'supplementary_figures';T=P/'source_tables'
for d in [P,F,S,T]:d.mkdir(parents=True,exist_ok=True)
RD=O.parent.parent/'main_methods_pass_2026-09-24/review_decisions'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
font=Path('STYLEGUIDE_DIR/soehne-buch.ttf');font_manager.fontManager.addfont(str(font));family=font_manager.FontProperties(fname=str(font)).get_name()
for n in [1,2,3,5,6]:
    for ext in ['pdf','png']:shutil.copy2(B/f'figures/Figure{n}.{ext}',F/f'Figure{n}.{ext}')
def export(doc,path,width=None):
    if width:
        src=doc;doc=fitz.open();w=width/25.4*72;p=doc.new_page(width=w,height=src[0].rect.height*w/src[0].rect.width);p.show_pdf_page(p.rect,src,0)
    if path.exists():path.unlink()
    doc.save(path,garbage=4,deflate=True);p=fitz.open(path)[0];p.get_pixmap(dpi=300,alpha=False).save(path.with_suffix('.png'));p.get_pixmap(dpi=110,alpha=False).save(path.with_name(path.stem+'_preview.png'))
    return {'sha256':sha(path),'png_sha256':sha(path.with_suffix('.png')),'width_mm':p.rect.width/72*25.4,'height_mm':p.rect.height/72*25.4}
# Preserve all a-g artwork through the last complete antibiotic label; replace i by h.
old=fitz.open(B/'figures/Figure4_native.pdf');upper=315.3
builder=RD/'M20_crossomics_populations/prevalence_comparison/integration/build_figures.py';ns={'__file__':str(builder)};exec(compile(builder.read_text().split('# Replace only the current cross-omics panels')[0],str(builder),'exec'),ns);ns['OUT']=F
plt.rcParams.update({'font.family':family,'pdf.fonttype':42,'svg.fonttype':'none'})
trip=fitz.open(ns['triptych']('PREDICT',old[0].rect.width,157,'h'))
doc=fitz.open();p=doc.new_page(width=old[0].rect.width,height=upper+8+157);p.show_pdf_page(fitz.Rect(0,0,p.rect.width,upper),old,0,clip=fitz.Rect(0,0,p.rect.width,upper));p.show_pdf_page(fitz.Rect(0,upper+8,p.rect.width,p.rect.height),trip,0)
oldvals=json.loads((B/'figures/CROSSOMICS_VALUES_UNCHANGED.json').read_text());assert all(x in oldvals for x in ns['values'])
report={'Figure4':export(doc,F/'Figure4.pdf',183),'Figure4_current_crossomics_values_preserved':len(ns['values']),'unchanged_main_figures':[1,2,3,5,6]}
# Heatmap moved intact; the former panel letter h is unnecessary in a standalone supplement.
eds=fitz.open();p=eds.new_page(width=old[0].rect.width,height=500.9-upper);p.show_pdf_page(p.rect,old,0,clip=fitz.Rect(0,upper,old[0].rect.width,500.9));p.add_redact_annot(fitz.Rect(0,0,10,12),fill=(1,1,1));p.apply_redactions(images=0,graphics=0)
for box in p.search_for('Enterotype Dysbiosis Score across 85 patients'):p.add_redact_annot(box,fill=(1,1,1))
p.apply_redactions(images=0,graphics=0)
report['EDS_moved']=export(eds,S/'Supplementary_Figure_08.pdf',180)
# All retained supplementary figures get their new first-citation numbers.
mapping={1:1,2:2,3:3,4:4,5:5,6:6,7:9,8:10,9:11,10:12,11:13}
for oldn,newn in mapping.items():
    for ext in ['pdf','png']:shutil.copy2(B/f'supplementary_figures/Supplementary_Figure_{oldn:02}.{ext}',S/f'Supplementary_Figure_{newn:02}.{ext}')
# Remove Wuyts from the displayed e/f grid and its underlying table; retain other methods.
original=RD/'M20_crossomics_populations/supplement_figure_rebuild/build_crossomics_supplement.py';code=original.read_text();code=code.replace("OUT = Path(__file__).resolve().parent/'output'; OUT.mkdir(exist_ok=True)","OUT = Path(__file__).resolve().parent/'s7_rebuild'; OUT.mkdir(exist_ok=True)")
needle="meth = pd.read_csv(SUP/'method_comparison/SOURCE_DATA_method_comparisons.tsv',sep='\\t')"
assert needle in code;code=code.replace(needle,needle+"\nmeth = meth[meth.features.ne('wuyts_recipe')].copy()\nmeth.to_csv(OUT/'SOURCE_DATA_method_comparisons_retained.tsv',sep='\\t',index=False)")
code=code.replace("fig.suptitle('Cross-cohort functional and native-profile robustness',fontsize=11,fontweight='bold')",'')
(W/'rebuild_s7.py').write_text(code)
exec(compile(code,str(original),'exec'),{'__file__':str(W/'rebuild_s7.py')})
# Reuse original a-d pixels and vector objects; splice only the lower method grids.
before=fitz.open(B/'supplementary_figures/Supplementary_Figure_07.pdf');after=fitz.open(W/'s7_rebuild/CROSSOMICS_ROBUSTNESS.pdf');doc=fitz.open();p=doc.new_page(width=before[0].rect.width,height=before[0].rect.height);split=472
p.show_pdf_page(fitz.Rect(0,0,p.rect.width,split),before,0,clip=fitz.Rect(0,0,p.rect.width,split));p.show_pdf_page(fitz.Rect(0,split,p.rect.width,p.rect.height),after,0,clip=fitz.Rect(0,split,p.rect.width,p.rect.height))
report['Supplement9_Wuyts_removed']=export(doc,S/'Supplementary_Figure_09.pdf')
for file in (W/'s7_rebuild').glob('*.tsv'):shutil.copy2(file,T/file.name)
methods=pd.read_csv(T/'SOURCE_DATA_method_comparisons_retained.tsv',sep='\t');assert not methods.features.str.contains('wuyts',case=False).any();report['retained_distance_estimates']=len(methods)
# Assembly/yield association: specimen-balanced scatter plus the reported acquisition-level tertiles.
d=pd.read_csv(RD/'MM6_capsule_metagenome_2026-09-25/capscan_assembly_dna_yield_gain_by_acquisition.tsv',sep='\t');d['gain_exact']=100*(d.FL_peptides/d.MG_peptides-1)
spec=d.groupby(['biosample','subject','region'],as_index=False).agg(assembly_proteins=('assembly_proteins','median'),gain_percent=('gain_exact','median'),acquisitions=('sample','size'))
caps=d[d.region.ne('stool')].copy();caps['assembly_tertile']=pd.qcut(caps.assembly_proteins,3,labels=['Smallest','Middle','Largest']);summary=caps.groupby('assembly_tertile',observed=True).agg(median_gain=('gain_exact','median'),n=('sample','size'),min_assembly=('assembly_proteins','min'),max_assembly=('assembly_proteins','max')).reset_index()
assert np.allclose(summary.median_gain,[80.1655,68.594,50.386],atol=.002)
spec.to_csv(T/'capsule_assembly_gain_per_specimen.tsv',sep='\t',index=False);caps.to_csv(T/'capsule_assembly_gain_acquisitions.tsv',sep='\t',index=False);summary.to_csv(T/'capsule_assembly_gain_tertiles.tsv',sep='\t',index=False)
plt.rcParams.update({'font.family':family,'font.size':8,'axes.titlesize':9,'pdf.fonttype':42})
fig,axes=plt.subplots(1,2,figsize=(7.2,3.5),gridspec_kw={'width_ratios':[1.5,1]});colors={'duodenum':'#2878a3','jejunum':'#60a5bf','ileum':'#8ba777','colon':'#df8f44','stool':'#777777'}
for reg,c in colors.items():
    q=spec[spec.region==reg];axes[0].scatter(q.assembly_proteins,q.gain_percent,s=11,c=c,label=reg.capitalize(),alpha=.8,linewidths=0)
axes[0].set(xscale='log',xlabel='Matched metagenome proteins per specimen',ylabel='FastaLake peptide gain over MG (%)');axes[0].set_title('a  Specimen-level association',loc='left');axes[0].legend(frameon=False,fontsize=7,ncol=2)
axes[0].set_ylim(0,spec.gain_percent.max()*1.2);axes[0].set_yscale('symlog',linthresh=50);axes[0].set_yticks([0,50,100,200,500,1000,2000],[0,50,100,200,500,1000,2000])
axes[1].bar(range(3),summary.median_gain,color='#2878a3',width=.65);axes[1].set_xticks(range(3),summary.assembly_tertile);axes[1].set(xlabel='Capsule assembly-size tertile',ylabel='Median acquisition-level gain (%)',ylim=(0,100));axes[1].set_title('b  Capsule acquisitions',loc='left')
for i,r in summary.iterrows():axes[1].text(i,r.median_gain+2,f'{r.median_gain:.1f}%',ha='center');axes[1].text(i,4,f'n = {r.n}',ha='center',color='white',fontsize=7)
for ax in axes:ax.spines[['top','right']].set_visible(False)
fig.tight_layout();fig.savefig(S/'Supplementary_Figure_07.pdf');fig.savefig(S/'Supplementary_Figure_07.png',dpi=300);fig.savefig(S/'Supplementary_Figure_07_preview.png',dpi=130);plt.close(fig)
report['capsule_panel']={'specimens':len(spec),'capsule_acquisitions':len(caps),'tertiles':summary.to_dict(orient='records'),'association_not_causal':True}
# Integrate the independently derived fusion decomposition as a continuation of S4.
for ext in ['pdf','png']:shutil.copy2(B/f'supplementary_additions/Supplementary_Figure_4_decomposition_proposed.{ext}',S/f'Supplementary_Figure_04_continued.{ext}')
for file in ['fusion_decomposition_current375.tsv','null_control_by_acquisition_seed.tsv','null_control_summary.tsv']:shutil.copy2(B/'supplementary_additions'/file,T/file)
report['supplement_renumbering']=mapping;(W/'FIGURE_CHANGES.json').write_text(json.dumps(report,indent=2));print('FIGURES_V15_PASS',report['capsule_panel'],flush=True)
