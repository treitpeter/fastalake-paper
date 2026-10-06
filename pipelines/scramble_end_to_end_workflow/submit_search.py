from pathlib import Path
import json,subprocess,hashlib,py_compile
H=Path('review/v15_integration/analyses/terminal_scramble').resolve();A=Path('MM_review_analysis/A02_4_2_lookup_specificity').resolve()
m=json.loads((H/'manifest.json').read_text());c=json.loads(Path('00_PACKAGE_2026-09-17/cleanup_2026-09-22/campi_fusion_rerun/config.json').read_text())
tasks=[]
for arm,seed,pred,union in [('authentic',0,A/'inputs/authentic.csv',A/'derived/lookup/matches.fasta')]+[(arm,seed,(H/('Cterm_'+str(seed))/(m['sample']+'_predictions.csv') if arm=='terminal' else A/'derived/predictions'/('SIM_'+str(seed))/(m['sample']+'_predictions.csv')),(H/'full_reservoir/matches.fasta' if arm=='terminal' else A/'derived/lookup/matches.fasta')) for seed in m['seeds'] for arm in ['all_position','terminal']]:
 assert pred.exists()
 tasks.append(dict(arm=arm,seed=seed,predictions=str(pred),predictions_sha256=hashlib.sha256(pred.read_bytes()).hexdigest(),candidate_union=str(union)))
cfg=Path('canonical/publication_data/arms/fig1b_shuffle_control_2026-09-01/REAL/sage/config.json').resolve()
out=dict(tasks=tasks,binaries={k:c['binaries'][k] for k in ['extractor','parsimony','sage']},sage_config=str(cfg),sage_config_sha256=hashlib.sha256(cfg.read_bytes()).hexdigest())
(H/'search_manifest.json').write_text(json.dumps(out,indent=2));py_compile.compile(str(H/'search_worker.py'),doraise=True)
(H/'END_TO_END_PROTOCOL.md').write_text("""# Full five-seed paired scramble experiment
User requested full execution on 28 September, not pilot-only delivery. Run one fresh authentic arm and all five all-position/fixed-C-terminal seed pairs: eleven complete selection/parsimony/Sage searches on the same worked acquisition. Both control families are searched anew using the same pinned binaries/configuration; preserve historical comparisons independently.

Each seed re-extracts its exact candidates from a complete union returned by that family's full-reservoir lookup. This is lossless: every protein matching any seed query was retained in that union, and per-seed extraction again applies the production filters. Do not use a top-hit subset. Fixed-C-terminal union depends on completed full scan. Stage counts, selected FASTAs, canonical accepted identities, q/PEP sensitivities and zeros are saved. In addition to ordinary q<=.01, expose 5–8 and 9–50 residues and PEP<=.01.

One serial Slurm array, eight CPUs/24GB/90min per task, nice10000. Original raw input/configuration/predictions and old outputs remain untouched. No new null control should be called FDR calibration: random scrambled queries are not the distribution of actual de novo errors, and an incidental genuine target can be recovered. The endpoint is whether terminal preservation changes chance lookup, protein selection or final target-peptide acceptance.
""")
s=json.loads((H/'SUBMISSIONS.json').read_text());assert 'end_to_end' not in s
cmd=['sbatch','--parsable','--job-name=FL_ML_Cterm_search','--partition=p.hpcl8,p.hpcl92,p.hpcl93','--nice=10000','--cpus-per-task=8','--mem=24G','--time=01:30:00','--array=0-10%1','--dependency=afterok:'+s['lookup']['job_id'],'--output='+str(H/'logs/end_to_end_%A_%a.out'),'--error='+str(H/'logs/end_to_end_%A_%a.err'),'--wrap','/mpcdf/soft/SLE_15/packages/x86_64/anaconda/3/2023.03/bin/python '+str(H/'search_worker.py')+' $SLURM_ARRAY_TASK_ID']
r=subprocess.run(cmd,text=True,capture_output=True);r.check_returncode();j=r.stdout.strip().split(';')[0];s['end_to_end']=dict(job_id=j,command=cmd);(H/'SUBMISSIONS.json').write_text(json.dumps(s,indent=2));print('FULL END TO END',j,'11 searches')

