from common import *
import importlib.util
mode=sys.argv[1];assert mode in ['baseline','refit'];verify()
for p,h in load(H/'NESTED_INPUTS.json').items():assert sha(p)==h,p
assert sha(__file__)==load(H/'NESTED_CODE.json')['nested_eds.py']
if mode=='refit':gate(H/'nested_baseline')
out=H/('nested_'+mode);out.mkdir(exist_ok=False)
spec=importlib.util.spec_from_file_location('nested_source',N/'code/analyze_score.py');a=importlib.util.module_from_spec(spec);spec.loader.exec_module(a)
a.ROOT=R;a.A19=N/'inputs';a.H=out;a.D=out
oldinputs=a.inputs
if mode=='refit':
 def replacement():
  meta,frames,hist=oldinputs();new=frame(C/'profiles/PREDICT_noclust','OG_other_unresolved').loc[meta.sample_id]
  cov=rd(C/'clinical_refit/acquisition_accounting.tsv').set_index('sample_id').loc[meta.sample_id]
  assert np.array_equal(cov.patient_id.astype(int),meta.idpatient.astype(int))
  for grouping in frames:
   # Both models use the same accepted peptide evidence and revised function allocation.
   # Grouping sensitivity survives only through its frozen original/simplified depth.
   meta[grouping+'_Human']=cov.Human.to_numpy()
   meta[grouping+'_Microbial']=(cov.Microbial_species+cov.Microbial_above_species).to_numpy()
   meta[grouping+'_Mixed']=cov.Other_annotated_origins.to_numpy();meta[grouping+'_Unknown']=cov.Unresolved_origin.to_numpy()
   frames[grouping]=(new.to_numpy(float),np.column_stack([np.log10(meta[grouping+'_pg_parsimony']),cov.Human]))
  return meta,frames,hist
 a.inputs=replacement
# Original source verifies cached nested predictions against direct sklearn fits.
a.main(pilot=False)
checks={}
if mode=='baseline':
 for name in ['model_metrics','heldout_scores','permutation_results','selected_alphas','outer_folds','increment_over_baseline','phenotype_associations']:
  got=rd(out/(name+'.tsv'));old=rd(N/'full/source_data'/(name+'.tsv'));assert list(got)==list(old) and got.shape==old.shape,name
  errors={}
  for col in got:
   if pd.api.types.is_numeric_dtype(old[col]):
    assert np.array_equal(got[col].isna(),old[col].isna()),(name,col)
    err=np.nanmax(abs(got[col]-old[col]),initial=0);assert np.allclose(got[col],old[col],rtol=1e-7,atol=1e-9,equal_nan=True),(name,col,err);errors[col]=float(err)
   else:assert got[col].fillna('').equals(old[col].fillna('')),(name,col)
  checks[name]=errors
else:
 old=rd(H/'nested_baseline/model_metrics.tsv');new=rd(out/'model_metrics.tsv');keys=['grouping','seed','model'];wr(old.merge(new,on=keys,suffixes=('_previous','_revised')),out/'method_comparison.tsv')
 for f in ['outer_folds.tsv']:
  assert rd(out/f).equals(rd(H/'nested_baseline'/f)),f
receipt(out,mode=mode,patients=85,full_eight_test_family=True,three_split_seeds=True,baseline_replay=checks,source_code_sha256=sha(N/'code/analyze_score.py'),driver_sha256=sha(__file__),scope='Internally evaluated; no external calibration. Peptide functional evidence is common across grouping controls; original and simplified group depth retained separately.')
