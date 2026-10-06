"""Supplementary 9, 10 and 19 drawn from the currently selected A33 tables.

Conceptions approved 20 September 2026 (ITEM8_RENDERING_DECISION_2026-09-20.md). Supplementary 9
and 19 had no drawing code (typography derivatives of older images); Supplementary 10 keeps the
layout of MM_Hartout_Figure_review/code/build_figures.py technical() with current sources, and its
panel c takes the four H01 eligibility definitions from render/Supp10/s10c_check. Supplementary 19c
is not drawn: the agreement analysis is shown once, in Supplementary 10c.
"""
import os,sys,json,hashlib,re
from pathlib import Path
os.environ['MPLBACKEND']='Agg';os.environ.setdefault('FL_STYLE_DIR','STYLEGUIDE_DIR')
ROOT=Path('LAB_ROOT');sys.path.insert(0,str(ROOT/'figure_style'))
import fl_style as fl,matplotlib.pyplot as plt,numpy as np,pandas as pd
fl.init();plt.rcParams['savefig.bbox']=None;plt.rcParams.update({'axes.labelcolor':fl.INK,'text.color':fl.INK,'xtick.color':fl.INK,'ytick.color':fl.INK,'font.size':6.5,'axes.labelsize':6.5,'xtick.labelsize':6,'ytick.labelsize':6,'legend.fontsize':6})
OUT=Path(__file__).resolve().parent;A33=OUT.parents[1];RES=A33/'results'
assert A33.name=='A33_4_6_peptide_centric_downstream_2026-09-19' and OUT.parent.name=='render',OUT
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest();MMW=1/25.4
def passed(sel):
    c=json.loads((RES/sel/'COMPLETE.json').read_text());assert c.get('status')=='PASS',sel
    return {'COMPLETE.json':sha(RES/sel/'COMPLETE.json'),'RUN_SELECTION.json':sha(RES/sel/'RUN_SELECTION.json') if (RES/sel/'RUN_SELECTION.json').exists() else None}
class Fig:
    def __init__(s,name,w,h):s.name,s.w,s.h=name,w,h;s.fig=plt.figure(figsize=(w*MMW,h*MMW));s.src={};s.num={};s.sel={}
    def ax(s,x,y,w,h):return s.fig.add_axes([x/s.w,1-(y+h)/s.h,w/s.w,h/s.h])
    def read(s,p,**kw):p=Path(p);s.src[str(p.relative_to(ROOT))]=sha(p);return pd.read_csv(p,sep='\t',**kw)
    def letter(s,x,y,l):s.fig.text(x/s.w,1-y/s.h,l,fontsize=8,color=fl.INK,va='top',ha='left')
    def save(s):
        d=OUT/s.name;d.mkdir(exist_ok=True);s.fig.savefig(d/f'{s.name}.pdf');s.fig.savefig(d/f'{s.name}.png',dpi=300);s.fig.savefig(d/f'{s.name}.svg')
        outs={f'{s.name}.{e}':sha(d/f'{s.name}.{e}') for e in ('pdf','png','svg')}
        json.dump(dict(sources=s.src,selections=s.sel,numbers=s.num,outputs_sha256=outs,size_mm=[s.w,s.h],code_sha256=sha(__file__)),open(d/f'{s.name}_PANELS.json','w'),indent=2,default=str);plt.close(s.fig)

TIERS=[('TechnicalReplicates','Technical\nreplicates'),('WorkflowReplicates','Workflow\nreplicates'),('OnePiece','Spatial\nsampling'),('Intra','Intra-\nindividual'),('Interindividual','Inter-\nindividual')]
# display names of the published twelve-function mosaic (figure_restoration Figure3 e_mosaic)
NAMES={'K00134':'Glyceraldehyde-3-P DH','K01006':'Pyruvate phosphate dikinase','K03737':'Pyruvate oxidoreductase','K02358':'Elongation factor Tu','K02355':'Elongation factor G','K01176':'Alpha-amylase','K00262':'Glutamate dehydrogenase','K04043':'Chaperone DnaK','K00656':'Formate C-acetyltransferase','K10117':'Raffinose-family transporter','K01803':'Triosephosphate isomerase','K00927':'Phosphoglycerate kinase'}
PAL=[fl.FL,fl.METAG,fl.GREEN,fl.VERM,fl.SKY,fl.PURPLE,fl.YELLOW,'#1d2433','#8c6d46','#6fa8a1','#b56576','#9b9b3a'];OTHER='#c9c9c9';NA='#ffffff'
def spname(x):
    g,_,sp=x.partition(' ');g=re.sub(r'_[A-Z]+$','',g);sp=re.sub(r'_[A-Z]+$','',sp);return (g[0]+'. '+sp) if sp and not sp.startswith('sp') else g+' '+sp

def supp9():
    F=Fig('S09',180,112);S=RES/'HSTOOL__simplified/source_data/Figure3';F.sel['HSTOOL__simplified']=passed('HSTOOL__simplified')
    t=F.read(S/'function_tier_summary.tsv');summ=json.loads((S/'SUMMARY.json').read_text());F.src[str((S/'SUMMARY.json').relative_to(ROOT))]=sha(S/'SUMMARY.json')
    e=t[t.eligible_in_all_five_designs.eq(1)];assert e.ko.nunique()==summ['n_common_functions']==582,e.ko.nunique()
    # ---- a: rarefied distinct leaders per design ----
    a=F.ax(16,7,84,36);fl.style(a);F.letter(2,2,'a');col='expected_distinct_leaders_in_8_resolved_acquisitions'
    for i,(tier,_) in enumerate(TIERS):
        v=e[e.tier.eq(tier)][col].astype(float);assert len(v)==582,(tier,len(v));q1,m,q3=np.percentile(v,[25,50,75])
        assert abs(m-summ['rarefied_medians'][tier])<1e-9,(tier,m)
        a.plot([i,i],[q1,q3],color=fl.FL,lw=1.2,zorder=2);a.scatter([i],[m],s=16,color=fl.FL,zorder=3);a.text(i+.16,m,f'{m:.1f}',va='center',ha='left',fontsize=5.8)
        F.num[f'a:{tier}']=dict(median=float(m),q1=float(q1),q3=float(q3),n_functions=int(len(v)))
    a.set_xticks(range(5),[l for _,l in TIERS]);a.set_xlim(-.5,4.7);a.set_ylim(0,6.5);a.set_ylabel('Distinct leading species per function\n(expected in 8 resolved acquisitions)')
    # ---- b: leader mosaic, the published twelve functions ----
    top=t[t.main_top12.eq(1)];kos=list(top[top.tier.eq('TechnicalReplicates')].sort_values('mean_intensity_pct_equal_three_designs',ascending=False).ko)
    assert set(kos)==set(NAMES) and len(kos)==12,kos
    pa=F.read(S/'per_function_acquisition.tsv',usecols=['tier','sample','ko','resolved_fraction','leading_species','n_resolved_contributors'],keep_default_na=False)
    pa=pa[pa.ko.isin(kos)];des=F.read(S/'design.tsv');assert len(des)==58
    # colours by species identity: Figure 5f's species keep their Figure 5f colour; the most frequent
    # remaining leaders here take the Figure 5f colours that no species of this panel uses.
    mr=F.read(RES/'CAPSCAN__simplified/source_data/Figure5/mosaic_donor_regions.tsv');F.sel['CAPSCAN__simplified']=passed('CAPSCAN__simplified')
    f5=list(mr[mr.status.eq('Resolved')].species.value_counts().index[:12]);f5col=dict(zip(f5,PAL))
    lead=pa[pa.leading_species.ne('')].leading_species.value_counts();colmap={sp:f5col[sp] for sp in lead.index if sp in f5col}
    free=[c for c in PAL if c not in colmap.values()]
    for sp in lead.index:
        if sp not in colmap and free:colmap[sp]=free.pop(0)
    b=F.ax(44,56,112,40);b.axis('off');F.letter(2,50,'b');x=0;gap=1.4;look=pa.set_index(['ko','sample']);xS={}
    for tier,lab in TIERS:
        ss=sorted(des[des.group.eq(tier)].sample_id);x0=x
        for s_ in ss:
            for i,ko in enumerate(kos):
                r=look.loc[(ko,s_)];c=NA if r.leading_species=='' else colmap.get(r.leading_species,OTHER)
                b.add_patch(plt.Rectangle((x,i-.5),1,1,facecolor=c,edgecolor='#e6e6e6' if c==NA else 'white',lw=.2))
            x+=1
        b.text((x0+x)/2,-1.0,lab.replace('\n',' ').replace('- ','-')+f'\nn = {len(ss)}',ha='center',va='bottom',fontsize=5.5,linespacing=1.1)
        xS[tier]=x+.9;x+=2.6+gap
    for i,ko in enumerate(kos):
        b.text(-.8,i,NAMES[ko],ha='right',va='center',fontsize=5.5)
        for tier,_ in TIERS:
            d=pa[pa.ko.eq(ko)&pa.tier.eq(tier)&pa.leading_species.ne('')];n=int(d.leading_species.nunique());b.text(xS[tier],i,str(n),ha='center',va='center',fontsize=5.3,color=fl.INK)
            ref=int(top[(top.ko==ko)&(top.tier==tier)].n_distinct_leaders.iloc[0]);assert n==ref,(ko,tier,n,ref);F.num[f'b:S:{ko}:{tier}']=n
    for tier,_ in TIERS:b.text(xS[tier],-1.0,'S',ha='center',va='bottom',fontsize=5.5)
    b.set_xlim(-.5,x);b.set_ylim(11.6,-2.2)
    F.num['b:technical_median_S']=float(np.median([F.num[f'b:S:{k}:TechnicalReplicates'] for k in kos]));F.num['b:interindividual_median_S']=float(np.median([F.num[f'b:S:{k}:Interindividual'] for k in kos]))
    F.num['b:functions_with_multiple_technical_leaders']=sorted(k for k in kos if F.num[f'b:S:{k}:TechnicalReplicates']>1)
    # abundance bars
    ab=F.ax(160,56,14,40);fl.style(ab);v=[float(top[(top.ko==k)&(top.tier=='TechnicalReplicates')].mean_intensity_pct_equal_three_designs.iloc[0]) for k in kos]
    F.num['b:mean_LFQ_pct']={k:float(x) for k,x in zip(kos,v)};ab.barh(range(12),v,color=fl.GREY,height=.62);ab.set_ylim(11.6,-2.2);ab.set_yticks([]);ab.spines['left'].set_visible(False);ab.set_xlabel('Mean LFQ (%)');ab.set_xlim(0,2.4);ab.set_xticks([0,1,2])
    # legend
    shown=[sp for sp in lead.index if sp in colmap][:12];lg=F.ax(44,99,112,13);lg.axis('off');lg.set_xlim(0,112);lg.set_ylim(13,0)
    short={sp:spname(sp) for sp in shown};dup={n for n in short.values() if list(short.values()).count(n)>1}
    # two GTDB taxa can shorten to one name; those keep their full GTDB label so the key stays one-to-one
    items=[(colmap[sp],sp.replace('_',' ') if short[sp] in dup else short[sp],'white') for sp in shown]+[(OTHER,'Other species','white'),(NA,'Unresolved','#bdbdbd')]
    for j,(c,name,ec) in enumerate(items):
        cx=(j%4)*28;cy=1.5+(j//4)*3.1;lg.add_patch(plt.Rectangle((cx,cy-.9),2.2,1.8,facecolor=c,edgecolor=ec,lw=.4));lg.text(cx+3.2,cy,name,va='center',ha='left',fontsize=5.5,style='italic' if name not in('Other species','Unresolved') else 'normal')
    F.num['b:legend']=[n for _,n,_ in items];F.num['b:functions']=kos;F.save();return F

def supp10():
    F=Fig('S10',183,66);d=F.read(RES/'S10_technical_leaders/full/SUMMARY.tsv');F.sel['S10_technical_leaders/full']={'COMPLETE.json':sha(RES/'S10_technical_leaders/full/COMPLETE.json')}
    assert json.loads((RES/'S10_technical_leaders/full/COMPLETE.json').read_text()).get('status')=='PASS'
    for k,(tier,xl) in enumerate([('TechnicalReplicates','Technical injections'),('WorkflowReplicates','Workflow preparations')]):
        a=F.ax(24+k*52,6,34,42);fl.style(a);F.letter((2,70)[k],2,'ab'[k]);labs=[];i=0
        for subset,short in [('all_eligible','All eligible'),('author_top12','Historical 12')]:
            for gate,glab in [('primary','Primary'),('well_resolved','Strict')]:
                q=d[d.tier.eq(tier)&d.subset.eq(subset)&d['filter'].eq(gate)].set_index('method');assert set(q.index)=={'original','simplified'},q
                for off,method,colr in [(-.12,'original',fl.GREY),(.12,'simplified',fl.FL)]:
                    val=100*float(q.loc[method,'equal_KO_agreement']);a.scatter(val,i+off,color=colr,s=14,label=method.capitalize() if i==0 else None,zorder=3);F.num[f'{"ab"[k]}:{subset}:{gate}:{method}']=round(val,2)
                labs.append(short+'\n'+glab);i+=1
        a.set_yticks(range(4),labs if k==0 else ['']*4);a.set_ylim(3.6,-.6);a.set_xlim(0,103);a.set_xlabel(xl+'\nEqual-KO same-leader pairs (%)')
        if k==0:a.legend(frameon=False,loc='lower left',handletextpad=.2,borderaxespad=0)
    c=F.read(OUT.parent/'Supp10/s10c_check/s10c_rule_comparison.tsv');F.src[str((OUT.parent/'Supp10/s10c_check/S10C_CHECK.json').relative_to(ROOT))]=sha(OUT.parent/'Supp10/s10c_check/S10C_CHECK.json')
    assert json.loads((OUT.parent/'Supp10/s10c_check/S10C_CHECK.json').read_text())['status']=='PASS';F.sel['CAPSCAN__simplified']=passed('CAPSCAN__simplified')
    s=c[c.source.eq('A33_CAPSCAN__simplified')].set_index('rule');order=['primary','strict_pair_filter','strict_donor_cell_eligibility','strict_paired_donor_eligibility'];assert list(s.index)==order and (s.n_donors==14).all() and (s.n_donor_regions==56).all()
    EXP={'primary':(45.03,1401,53617),'strict_pair_filter':(56.18,1219,22512),'strict_donor_cell_eligibility':(61.46,358,11330),'strict_paired_donor_eligibility':(62.70,154,6466)}
    a=F.ax(144,6,34,42);fl.style(a);F.letter(116,2,'c')
    for i,r in enumerate(order):
        val=100*float(s.loc[r,'donor_balanced_agreement']);nk,npair=int(s.loc[r,'n_kos']),int(s.loc[r,'n_pairs']);assert (round(val,2),nk,npair)==EXP[r],(r,val,nk,npair)
        a.scatter(val,i,color=fl.FL,s=14,zorder=3);a.text(2,i+.34,f'{nk:,} KOs · {npair:,} pairs',fontsize=5.2,va='center',ha='left',color=fl.INK);F.num[f'c:{r}']=dict(agreement_pct=round(val,2),n_kos=nk,n_pairs=npair,n_participants=14,n_participant_regions=56)
    a.set_yticks(range(4),['Primary','Strict pair filter','Strict participant cells','Strict paired participants']);a.set_ylim(3.6,-.6);a.set_xlim(0,103);a.set_xlabel('Capsule–stool same-leader\nagreement (%), simplified rule')
    F.save();return F

def supp19():
    F=Fig('S19',120,52);S=RES/'CAPSCAN__simplified/source_data/Figure5';F.sel['CAPSCAN__simplified']=passed('CAPSCAN__simplified')
    c=F.read(S/'contributor_cells.tsv',dtype={'participant':str});summ=json.loads((S/'SUMMARY.json').read_text())['contributors'];F.src[str((S/'SUMMARY.json').relative_to(ROOT))]=sha(S/'SUMMARY.json')
    assert len(c)==summ['eligible_cells']==20352 and c.KO.nunique()==summ['eligible_functions']==1416 and not c.duplicated(['KO','participant']).any()
    assert sorted(c.participant.unique(),key=int)==[str(i) for i in range(1,16)] and c.groupby('KO').participant.nunique().min()>=12
    F.num['n_cells']=int(len(c));F.num['n_functions']=int(c.KO.nunique());F.num['n_participants']=int(c.participant.nunique());F.num['participant_15_cells']=int(c.participant.eq('15').sum())
    for k,(colname,xl,key) in enumerate([('leader_share','Leading-species share of\nspecies-resolved intensity (%)','leader_share_median'),('species_resolved_fraction','Species-resolved share of\nfunction intensity (%)','median_resolved_fraction')]):
        a=F.ax(16+k*58,7,40,30);fl.style(a);F.letter(2+k*58,2,'ab'[k]);v=100*c[colname].astype(float);assert v.between(0,100.0000001).all();m=float(v.median());assert abs(m/100-summ[key])<1e-12,(colname,m)
        edges=np.arange(0,102.5,2.5);cnt,_=np.histogram(v,bins=edges);assert int(cnt.sum())==len(v);F.num['ab'[k]+':hist_edges_pct']=[float(x) for x in edges];F.num['ab'[k]+':hist_counts']=[int(x) for x in cnt]
        a.hist(v,bins=edges,color=fl.FL,lw=0);a.axvline(m,color=fl.INK,lw=.9);a.text(m+2 if m<50 else m-2,.97,f'Median {m:.1f}%',transform=a.get_xaxis_transform(),ha='left' if m<50 else 'right',va='top',fontsize=5.8)
        a.set_xlim(0,100);a.set_xlabel(xl);a.set_ylabel('Function–participant cells' if k==0 else '');F.num['ab'[k]+':median_pct']=round(m,4)
    F.save();return F

if __name__=='__main__':
    which=sys.argv[1:] or ['S09','S10','S19']
    for name,fn in (('S09',supp9),('S10',supp10),('S19',supp19)):
        if name in which:f=fn();print(name,json.dumps({k:v for k,v in f.num.items() if not k.startswith('b:S:')},default=str)[:1800])
