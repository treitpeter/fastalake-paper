"""Build the approved post hoc EDS examples from frozen patient-matched inputs.

Writes only EDS_approved in this passage-review branch. No discovery screen is
rerun or restricted to the displayed examples. Saved full-family q values remain
authoritative; the small additional adjustment decomposition is descriptive.
"""
from pathlib import Path
import hashlib
import itertools
import json
import shutil
from datetime import datetime, timezone
import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
import fitz
from PIL import Image

D = Path(__file__).resolve().parent
B = D.parent/'B'   # CR003 follow-up shadow: analysis tree of this variant
R = Path(__import__('os').environ['FASTALAKE_DATA_DIR'])
O = D/'EDS_approved'
O.mkdir(exist_ok=True)
S = O/'source_data'
S.mkdir(exist_ok=True)
F = O/'figures'
F.mkdir(exist_ok=True)
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p, v):
    p.write_text(json.dumps(v, indent=2, ensure_ascii=False)+'\n')
def wr(df, name):
    df.to_csv(S/name, sep='\t', index=False)

manifest = json.loads((B/'eds_refit/INPUTS.json').read_text())
for k in ['order', 'heatmap', 'features', 'metabolites', 'KO', 'annotation']:
    assert sha(manifest[k]['path']) == manifest[k]['sha256'], k
read = lambda k, **kw: pd.read_csv(manifest[k]['path'], sep='\t', **kw)
order = read('order')
legacy = np.load(manifest['heatmap']['path'])
order = order.set_index('sample_id').loc[legacy['samples']].reset_index()
assert len(order) == order.patient_id.nunique() == 85
meta = read('metabolites', index_col=0).loc[order.microbiome_id]
meta.index = order.sample_id
ko = read('KO', index_col=0).drop(index='FunctionallyUnassigned', errors='ignore').loc[:, order.sample_id].T
assert ko.index.equals(meta.index)
pairs = [
    ('HILIC_FT23003', 'K01812', 'Glucuronic acid', 'Glucuronate isomerase', 'R01482'),
    ('HILIC_FT24070', 'K01667', 'Tryptophan', 'Tryptophanase', 'R00673'),
    ('C18_F2423', 'K01442', 'Glycochenodeoxycholic acid', 'Bile salt hydrolase', 'R03975'),
]
prior = pd.read_csv(B/'eds_refit/augmented_pair_screen.tsv', sep='\t')
eds = pd.read_csv(B/'eds_refit/all_EDS_associations.tsv', sep='\t')
names = ['host', 'depth', 'antibiotics', 'clinical']
subsets = [c for n in range(5) for c in itertools.combinations(names, n)]
cats = {k: pd.get_dummies(order[k], drop_first=True).to_numpy(float) for k in ['antibiotic_state', 'clinical_group']}
def blocks(idx):
    return {'host': rankdata(order.Human.to_numpy()[idx])[:, None],
            'depth': rankdata(np.log2(order.pg_parsimony.to_numpy()[idx]))[:, None],
            'antibiotics': cats['antibiotic_state'][idx], 'clinical': cats['clinical_group'][idx]}
def estimate(x, y, block, subset):
    Z = np.column_stack([np.ones(len(x))]+[block[k] for k in subset])
    ranks = np.column_stack([rankdata(x), rankdata(y)])
    residuals = ranks-Z @ np.linalg.lstsq(Z, ranks, rcond=None)[0]
    rho = np.corrcoef(residuals.T)[0, 1]
    U, _, _ = np.linalg.svd(Z, full_matrices=False)
    U = U[:, :np.linalg.matrix_rank(Z)]
    check = ranks-U@(U.T@ranks)
    assert np.isclose(rho, np.corrcoef(check.T)[0, 1], atol=1e-10)
    return float(rho)

masks = {m: (np.isfinite(meta[m]) & (meta[m]>0) & np.isfinite(ko[k]) & (ko[k]>0)).to_numpy() for m,k,*_ in pairs}
common = np.logical_and.reduce(list(masks.values()))
rows, summary, rosters, features, raw = [], [], [], [], []
for m, k, ml, kl, reaction in pairs:
    saved = prior[(prior.layer=='KO') & (prior.metabolite==m) & (prior.protein_feature==k)].copy()
    assert len(saved) == 6
    wr(saved, f'{m}_{k}_full_family_results.tsv')
    original = saved[(saved.normalization=='native') & (saved.model=='unadjusted')].iloc[0]
    adjusted = saved[(saved.normalization=='native') & (saved.model=='host_depth_antibiotic_clinical')].iloc[0]
    eds_rows = eds[(eds.normalization=='native') & (eds.model=='unadjusted') &
                   (((eds.layer=='Metabolite') & (eds.feature==m)) | ((eds.layer=='KO') & (eds.feature==k)))]
    assert len(eds_rows)==2 and (eds_rows.q<.05).all()
    summary.append(dict(metabolite=m, metabolite_label=ml, KO=k, enzyme_label=kl,
                        reaction=reaction, reaction_url=f'https://www.kegg.jp/entry/{reaction}',
                        n=int(original.n), rho=float(original.rho), full_family_q=float(original.q_augmented),
                        joint_adjusted_rho=float(adjusted.rho), joint_adjusted_full_family_q=float(adjusted.q_augmented)))
    for scope, mask in [('pair_complete', masks[m]), ('common_84', common)]:
        idx = np.flatnonzero(mask)
        for i in idx:
            rosters.append(dict(metabolite=m, KO=k, scope=scope, sample_id=order.sample_id.iloc[i], patient_id=order.patient_id.iloc[i]))
        for norm in ['native', 'annotated_relative']:
            x = ko[k].to_numpy(float)
            if norm=='annotated_relative':
                x = x/ko.sum(axis=1).to_numpy(float)
            x, y = x[idx], meta[m].to_numpy(float)[idx]
            estimates = [estimate(x, y, blocks(idx), subset) for subset in subsets]
            assert np.isclose(estimates[0], spearmanr(x, y).statistic, atol=1e-12)
            if scope=='pair_complete':
                saved_norm = saved[saved.normalization==norm]
                assert np.isclose(estimates[0], saved_norm[saved_norm.model=='unadjusted'].rho.iloc[0], atol=1e-12)
                global_blocks = {key: value[idx] for key,value in blocks(np.arange(85)).items()}
                replay = estimate(x, y, global_blocks, tuple(names))
                assert np.isclose(replay, saved_norm[saved_norm.model=='host_depth_antibiotic_clinical'].rho.iloc[0], atol=1e-12)
            for subset, rho in zip(subsets, estimates):
                rows.append(dict(metabolite=m, KO=k, scope=scope, normalization=norm,
                                 model='+'.join(subset) if subset else 'none', n=len(idx), rho=rho,
                                 delta_from_unadjusted=rho-estimates[0]))
    for feature, label, layer, vals in [(m, ml, 'Metabolite', meta[m]), (k, f'{kl} ({k})', 'KO', ko[k])]:
        vals = vals.to_numpy(float)
        positive = np.isfinite(vals) & (vals>0)
        row = eds_rows[eds_rows.feature==feature].iloc[0]
        assert np.isclose(spearmanr(vals[positive], order.dysbiosis_score.to_numpy()[positive]).statistic, row.rho, atol=1e-12)
        features.append(dict(feature=feature, label=label, layer=layer, n_positive=int(positive.sum()),
                             EDS_rho=float(row.rho), EDS_full_family_q=float(row.q)))
        raw.append(vals)

assert common.sum()==84
models = pd.DataFrame(rows)
assert len(models)==192
wr(models, 'separate_adjustments_all_models.tsv')
wr(pd.DataFrame(rosters), 'pair_patient_rosters.tsv')
wr(pd.DataFrame(summary), 'selected_pair_summary.tsv')
wr(pd.DataFrame(features), 'focused_features.tsv')
wr(order, 'patient_order_and_covariates.tsv')
raw = np.array(raw)
pos = np.where(np.isfinite(raw) & (raw>0), raw, np.nan)
zz = (pos-np.nanmean(pos, axis=1)[:, None])/np.nanstd(pos, axis=1)[:, None]
np.savez_compressed(S/'focused_values.npz', raw=raw, z=zz,
                    features=np.array([f['feature'] for f in features]), samples=order.sample_id.to_numpy(dtype=str))
long = pd.DataFrame(raw.T, columns=[f['feature'] for f in features])
wr(pd.concat([order.reset_index(drop=True), long], axis=1), 'focused_raw_measurements.tsv')
for name in ['augmented_pair_screen.tsv', 'all_EDS_associations.tsv']:
    shutil.copy2(B/'eds_refit'/name, S/name)
legacy_features = read('features')
keep = np.flatnonzero(legacy_features.n_positive_85>0)
assert len(keep)==19
shutil.copy2(manifest['features']['path'], S/'legacy_features.tsv')
shutil.copy2(manifest['heatmap']['path'], S/'legacy_85_values.npz')
for name in ['heatmap_values.npz', 'heatmap_columns.tsv']:
    shutil.copy2(B/'eds_refit/input_data'/name, S/('legacy_259_'+name))

font = R.parent/'ADAPT/SHARED_PLOTTING_FUNCTIONS_STYLEGUIDE/soehne-buch.ttf'
font_manager.fontManager.addfont(font)
plt.rcParams.update({'font.family':font_manager.FontProperties(fname=font).get_name(), 'font.size':6.4,
                     'axes.linewidth':.5, 'pdf.fonttype':42, 'savefig.bbox':None})
cm = LinearSegmentedColormap.from_list('fl', ['#0072B2', 'white', '#E69F00'])
cm.set_bad('#dedede')
labels = [f['label'] for f in features]
oldlabels = [s.replace('Dimethylsuccinic acid', 'Dimethylsuccinic / methylglutaric / adipic acid') for s in legacy_features.iloc[keep].label]
def tracks(fig, left, right, top, height, fontsize):
    tt = [(order.dysbiosis_score.to_numpy(), 'EDS (0–1)', 'Greys', 1),
          (100*order.Human.to_numpy(), 'Human signal (0–100%)', LinearSegmentedColormap.from_list('human', ['white','#CC79A7']), 100),
          (order.antibiotic_state.map({'Off':0,'On':1}).to_numpy(), 'Antibiotics (dark: on)', ListedColormap(['#ededed','#333333']), 1)]
    for i,(v,l,c,m) in enumerate(tt):
        ax=fig.add_axes([left,top-(i+1)*height,right-left,height*.72])
        ax.imshow(v[None,:],cmap=c,vmin=0,vmax=m,aspect='auto',interpolation='nearest')
        ax.set(xticks=[],yticks=[0],yticklabels=[l]);ax.tick_params(length=0,labelsize=fontsize)
        for sp in ax.spines.values():sp.set_visible(False)
def matrix(fig, values, labs, rect, seps, fontsize, xticks=True):
    ax=fig.add_axes(rect)
    im=ax.imshow(np.clip(values,-2,2),cmap=cm,vmin=-2,vmax=2,aspect='auto',interpolation='nearest')
    ax.set_yticks(range(len(labs)),labs);ax.tick_params(axis='y',length=0,pad=3,labelsize=fontsize)
    if xticks:
        ax.set_xticks([0,21,42,63,84],['1','22','43','64','85'])
        ax.set_xlabel('85 patients ordered by metagenomic EDS',labelpad=3,fontsize=fontsize)
        ax.tick_params(axis='x',labelsize=fontsize)
    else:ax.set_xticks([])
    for sep in seps:ax.axhline(sep,color='#666666',lw=.6)
    for sp in ax.spines.values():sp.set_visible(False)
    return im
def colorbar(fig, im, rect, fontsize):
    cb=fig.colorbar(im,cax=fig.add_axes(rect),ticks=[-2,0,2])
    cb.set_label('Within-feature z score',labelpad=2,fontsize=fontsize)
    cb.ax.tick_params(length=2,pad=1,labelsize=fontsize)
def save(fig,name):
    fig.savefig(F/(name+'.pdf'));fig.savefig(F/(name+'.png'),dpi=180);plt.close(fig)

fig=plt.figure(figsize=(518.7401733398438/72,215/72))
fig.text(.006,.965,'h',fontsize=9,weight='bold',va='top')
tracks(fig,.385,.929,.97,.052,6.4)
im=matrix(fig,zz,labels,[.385,.135,.544,.665],[1.5,3.5],6.6)
colorbar(fig,im,[.945,.28,.011,.40],6.4)
save(fig,'PREDICT_EDS_focused_panel')

fig=plt.figure(figsize=(720/72,470/72))
fig.text(.006,.985,'a  Selected biochemical examples',fontsize=9,weight='bold',va='top')
tracks(fig,.385,.929,.945,.023,7)
im=matrix(fig,zz,labels,[.385,.65,.544,.218],[1.5,3.5],7,False)
fig.text(.006,.604,'b  Earlier EDS-selected overview',fontsize=9,weight='bold',va='top')
matrix(fig,legacy['z'][keep],oldlabels,[.385,.085,.544,.49],[4.5,10.5],7)
colorbar(fig,im,[.945,.32,.011,.40],7)
save(fig,'PREDICT_EDS_supplement')

# Recompose from the original vector panels. Wrapping the already composited
# PDF in further clipped XObjects loses its effective clipping on this renderer.
base=R/'review/v15_integration/integration_2026-09-29/four_author_integration_2026-10-02/actual_FAR_261002/FAR23_Figure4_applied/Figure4.pdf'
src=fitz.open(base)
original_page=src[0]
cut=fitz.Rect(0,331,original_page.rect.width,546)
out=fitz.open(); page=out.new_page(width=original_page.rect.width,height=original_page.rect.height)
rects=[fitz.Rect(33,203,169,323),fitz.Rect(171,203,321,331),cut,fitz.Rect(0,549,page.rect.width,page.rect.height)]
clips=[original_page.rect]
for rect in rects:
    pieces=[]
    for clip in clips:
        inter=clip & rect
        if inter.is_empty:
            pieces.append(clip);continue
        for piece in [fitz.Rect(clip.x0,clip.y0,clip.x1,inter.y0),fitz.Rect(clip.x0,inter.y1,clip.x1,clip.y1),
                      fitz.Rect(clip.x0,inter.y0,inter.x0,inter.y1),fitz.Rect(inter.x1,inter.y0,clip.x1,inter.y1)]:
            if not piece.is_empty:pieces.append(piece)
    clips=pieces
for clip in clips:page.show_pdf_page(clip,src,0,clip=clip)
panel_paths=[B/'figures/panel_e.pdf',B/'figures/panel_f.pdf',F/'PREDICT_EDS_focused_panel.pdf',B/'figures/panel_i.pdf']
for path,rect in zip(panel_paths,rects):
    with fitz.open(path) as panel:page.show_pdf_page(rect,panel,0)
out.save(F/'Figure4_EDS_updated.pdf',garbage=4,deflate=True)
# Reopen the saved PDF before rendering so its final resource references are used.
out.close();out=fitz.open(F/'Figure4_EDS_updated.pdf');page=out[0]
im=Image.open(B/'figures/Figure4_revised.png').convert('RGB')   # CR003 follow-up: this variant's build_figures raster (panels e,f,i) instead of the frozen review_package copy
before=np.asarray(im).copy()
pix=page.get_pixmap(matrix=fitz.Matrix(im.width/page.rect.width,im.height/page.rect.height),alpha=False)
rendered=Image.frombytes('RGB',[pix.width,pix.height],pix.samples)
box=(0,int(np.floor(331*im.height/page.rect.height)),im.width,int(np.ceil(546*im.height/page.rect.height)))
im.paste(rendered.crop(box),box[:2])
after=np.asarray(im)
assert np.array_equal(before[:box[1]],after[:box[1]]) and np.array_equal(before[box[3]:],after[box[3]:])
im.save(F/'Figure4_EDS_updated.png')
page.get_pixmap(dpi=110,alpha=False).save(F/'Figure4_EDS_preview.png')
for name in ['panel_h.pdf','panel_h.png','Supplementary_EDS_heatmap.pdf','Supplementary_EDS_heatmap.png']:
    shutil.copy2(B/'figures'/name,F/('previous_'+name))

dump(O/'BUILD_VALIDATION.json',dict(
    utc=datetime.now(timezone.utc).isoformat(),status='NUMERIC_PASS_VISUAL_REVIEW_PENDING',
    user_decision='ok i like', selection='Post hoc approval of three interpretable annotated metabolite–KO pairs; not a new or independently validated signature.',
    patients=85,common_pair_patients=int(common.sum()),focused_rows=6,legacy_rows=19,
    original_full_screen_preserved=True,prior_native_and_relative_correlations_reproduced=True,
    descriptive_adjustment_models=len(models),independent_projection_checked=True,
    new_p_values_or_FDR_claims=False,all_other_Figure4_pixels_unchanged=True,
    source_manifest=manifest,script_sha256=sha(__file__),
    full_family_tables={n:sha(S/n) for n in ['augmented_pair_screen.tsv','all_EDS_associations.tsv']},
    source_files={p.name:sha(p) for p in S.iterdir() if p.is_file()},
    figure_files={p.name:sha(p) for p in F.iterdir() if p.is_file()},
    related_MM_requests=['R20','R24','R27']))
print(models[(models.scope=='pair_complete') & (models.normalization=='native') & models.model.isin(['none',*names])][['metabolite','KO','model','n','rho']].to_string(index=False))
print('Build complete:',O)
