"""Audit the completed two-seed collection; aggregate existing diagnostics without new search rules."""
from pathlib import Path
import csv, datetime, hashlib, json, math, re, subprocess, sys, zipfile
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[2]
OUT=Path(__file__).resolve().parent
EXP=ROOT/"analyses/fusion_entrapment_2026-09-07"
FIG=ROOT/"FIGURATOR/ap_style_revision_2026-09-07/final"
RELEASE=ROOT/"release_prep/github_biorxiv_2026-09-07"
def sha(p):
 h=hashlib.sha256()
 with Path(p).open("rb") as f:
  for b in iter(lambda:f.read(8*1024*1024),b""):h.update(b)
 return h.hexdigest()
def save(name,j):(OUT/name).write_text(json.dumps(j,indent=2,default=str)+"\n")
def tsv(name,d):d.to_csv(OUT/name,sep="\t",index=False)
d=pd.read_csv(OUT/"snapshot/per_acquisition.tsv",sep="\t")
keys=["cohort","sample","seed","q_threshold","min_residues"]
assert not d.duplicated(keys).any()
assert len(d)==14200 and (d.n_canonical>=0).all()
assert d.loc[d.n_canonical.eq(0), 'paired_upper_diagnostic_with_unpaired_worst_case'].isna().all()
assert (d.n_target_only+d.n_entrapment_only+d.n_shared==d.n_canonical).all()
assert np.allclose((d.paired_numerator+d.n_unpaired_accepted)/d.n_canonical,d.paired_upper_diagnostic_with_unpaired_worst_case,equal_nan=True)
assert np.allclose((2*d.n_paired_entrapment+d.n_unpaired_accepted)/d.n_canonical,d.combined_upper_diagnostic_with_unpaired_worst_case,equal_nan=True)
assert not d.whole_procedure_fdr_control_established.any()
cfg=json.loads((EXP/"config.json").read_text())
expected={c:v["n"] for c,v in cfg["cohorts"].items()}
counts=d.groupby(["cohort","seed"]).agg(completed=("sample","nunique")).reset_index()
counts["expected"]=counts.cohort.map(expected)
tsv("completeness.tsv",counts)
# Verify every collected row against its completion manifest and separate summary.
manifest_checks=0;config_checks=0
original_rows={}
for c,v in cfg["cohorts"].items():
 original_rows[c]={r["sample_id"]:r for r in csv.DictReader(Path(v["sample_list"]).open(),delimiter="\t")}
for (c,sample,seed),g in d.groupby(["cohort","sample","seed"],sort=False):
 base=EXP/"runs_v2"/c/str(seed)/"samples"/sample
 done=json.loads((base/"COMPLETE.json").read_text())
 assert (done["cohort"],done["sample"],done["seed"])==(c,sample,seed)
 a=g.drop(columns=["cohort","sample","seed"]).reset_index(drop=True)
 b=pd.DataFrame(done["summary"]);z=pd.read_csv(base/"summary.tsv",sep="\t")
 pd.testing.assert_frame_equal(a,b,check_dtype=False,rtol=1e-12,atol=1e-14)
 pd.testing.assert_frame_equal(a,z,check_dtype=False,rtol=1e-12,atol=1e-14)
 original=json.loads(Path(original_rows[c][sample]["sage_config"]).read_text())
 actual=json.loads((base/"search_input/sage.json").read_text())
 for x in (original,actual):
  x.pop("output_directory",None);x.pop("mzml_paths",None);x["database"].pop("fasta",None)
 assert actual==original,(c,sample,seed)
 manifest_checks+=1;config_checks+=1
print("Verified",manifest_checks,"completion/summary records and frozen search configurations",flush=True)
def summarize(frame,group):
 rows=[]
 for k,g in frame.groupby(group):
  if not isinstance(k,tuple):k=(k,)
  r=dict(zip(group,k))
  for col in ["n_canonical","n_target_only","n_entrapment_only","n_shared","n_unpaired_accepted","n_paired_entrapment","paired_numerator"]:
   r[col]=int(g[col].sum())
  n=r["n_canonical"];u=r["n_unpaired_accepted"];e=r["n_paired_entrapment"]
  r.update(acquisitions=len(g),paired_pct=100*(r["paired_numerator"]+u)/n,combined_pct=100*(2*e+u)/n,unknown_pct=100*u/n,entrapment_only_pct=100*r["n_entrapment_only"]/n,paired_median_pct=100*g.paired_upper_diagnostic_with_unpaired_worst_case.median(),paired_q25_pct=100*g.paired_upper_diagnostic_with_unpaired_worst_case.quantile(.25),paired_q75_pct=100*g.paired_upper_diagnostic_with_unpaired_worst_case.quantile(.75),acquisitions_paired_above_1pct=int(g.paired_upper_diagnostic_with_unpaired_worst_case.gt(.01).sum()))
  rows.append(r)
 return pd.DataFrame(rows)
summary=summarize(d,["cohort","seed","q_threshold","min_residues"])
tsv("cohort_threshold_summary.tsv",summary)
q=summary[summary.q_threshold.eq(.01)];tsv("cohort_q001.tsv",q)
common=d.groupby(["cohort","sample"]).seed.nunique(); common=set(common[common.eq(2)].index)
matched=d[[ (c,s) in common for c,s in zip(d.cohort,d["sample"]) ]]
ms=summarize(matched,["cohort","seed","q_threshold","min_residues"]);tsv("matched_seed_summary.tsv",ms)
basecols=["n_canonical","n_entrapment_only","n_shared","n_unpaired_accepted","n_paired_entrapment","paired_numerator"]
short=d[d.min_residues.eq(5)].set_index(["cohort","sample","seed","q_threshold"])[basecols]-d[d.min_residues.eq(9)].set_index(["cohort","sample","seed","q_threshold"])[basecols]
assert (short>=0).all().all()
short["paired_upper_diagnostic_with_unpaired_worst_case"]=(short.paired_numerator+short.n_unpaired_accepted)/short.n_canonical
short["n_target_only"]=short.n_canonical-short.n_entrapment_only-short.n_shared
sh=summarize(short.reset_index(),["cohort","seed","q_threshold"]);tsv("short_5_to_8_summary.tsv",sh)
tsv("acquisition_outliers_q001.tsv",d[d.q_threshold.eq(.01)].sort_values("paired_upper_diagnostic_with_unpaired_worst_case",ascending=False))
# Independent raw arithmetic on the highest 9-50 diagnostic acquisition in each cohort/seed.
# Outcome-selected checks audit difficult cases; they do not estimate population error.
raw_audit=[]
sel=d[d.q_threshold.eq(.01)&d.min_residues.eq(9)].sort_values("paired_upper_diagnostic_with_unpaired_worst_case").groupby(["cohort","seed"]).tail(1)
for _,row in sel.iterrows():
 base=EXP/"runs_v2"/row.cohort/str(row.seed)/"samples"/row["sample"]
 done=json.loads((base/"COMPLETE.json").read_text())
 result=base/"5_sage/results.sage.tsv";digest=sha(result);assert digest==done["results_sha256"]
 p=pd.read_csv(result,sep="\t",usecols=["label","peptide","peptide_q","proteins"]);p=p[p.label.eq(1)].copy()
 p["canonical"]=p.peptide.str.replace(r"\[[^\]]*\]","",regex=True).str.replace("I","L",regex=False)
 assert p.canonical.str.fullmatch("[A-Z]+").all()
 score=p.groupby("canonical").peptide_q.min()
 pairs=pd.read_csv(base/"search_input/peptide_pairs.tsv",sep="\t")
 assert pairs.target.is_unique and pairs.entrapment.is_unique and not(set(pairs.target)&set(pairs.entrapment))
 pairs=pairs[pairs.target.str.len().ge(9)]
 accepted=set(score[score.le(.01)&(score.index.str.len()>=9)&(score.index.str.len()<=50)].index)
 e_acc=pairs.entrapment.isin(accepted);t_acc=pairs.target.isin(accepted)
 ne=int(e_acc.sum());eo=int((e_acc&~t_acc).sum())
 eb=int((e_acc&t_acc&(pairs.entrapment.map(score)<=pairs.target.map(score))).sum())
 unknown=len(accepted-set(pairs.target)-set(pairs.entrapment))
 num=ne+eo+2*eb
 assert (len(accepted),unknown,ne,num)==(row.n_canonical,row.n_unpaired_accepted,row.n_paired_entrapment,row.paired_numerator)
 raw_audit.append(dict(cohort=row.cohort,seed=int(row.seed),sample=row["sample"],results_sha256=digest,n=len(accepted),unknown=unknown,paired_entrapment=ne,entrapment_only_pairs=eo,entrapment_better_or_tied=eb,paired_numerator=num))
 print("Raw audit passed",row.cohort,row.seed,row["sample"],flush=True)
save("RAW_AUDIT.json",raw_audit)

save("VALIDATION.json",dict(completed_searches=manifest_checks,expected_searches=1420,unchanged_search_configurations=config_checks,independent_raw_recounts=len(raw_audit),whole_procedure_fdr_control_established=False))
assert manifest_checks==1420
print(q[["cohort","seed","min_residues","acquisitions","paired_pct"]].to_string(index=False))
