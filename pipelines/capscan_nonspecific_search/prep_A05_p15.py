"""Build an A05 (non-specific Sage 0.14.6) run directory for the 16 participant-15 capsules.
Frozen A05 code is copied unchanged; only the task list is new. Configs are derived from each
acquisition's own production tryptic results.json with exactly the A05 transformation
(verified against the frozen A05 config of 1249: cleave_at '', restrict null, c_terminal true,
semi_enzymatic false, lfq false, Dari raw mzML, output_directory)."""
import os
import json, shutil, sys, pathlib, hashlib
import pandas as pd
F6 = pathlib.Path(sys.argv[1]); RUN = F6/'A05'
SRC = pathlib.Path(os.environ['FASTALAKE_DATA_DIR'], 'MM_review_analysis/A05_4_7_unspecific')
sys.path.insert(0, str(SRC/'frozen')); from lib import meta
CAP = pathlib.Path(os.environ['FASTALAKE_DATA_DIR'], 'canonical/publication_data/arms/noclust_2026-08-07/cohorts/CAPSCAN')
RAW = pathlib.Path(os.environ['FASTALAKE_CAPSCAN_MZML_DIR'])
p15 = pd.read_csv(F6.parent/'roster/participant15_added.tsv', sep='\t')
old = json.loads((SRC/'inputs/manifest.json').read_text())
(RUN/'frozen').mkdir(parents=True, exist_ok=False); (RUN/'inputs').mkdir(); (RUN/'results').mkdir(); (RUN/'logs').mkdir()
for f in ['search.py', 'lib.py', 'peptide.py']: shutil.copy2(SRC/'frozen'/f, RUN/'frozen'/f)
for a, b in zip(sorted(old['frozen_code'], key=lambda r: r['path']), sorted((RUN/'frozen').glob('*.py'))):
    assert hashlib.sha256(b.read_bytes()).hexdigest() == a['sha256'] and pathlib.Path(a['path']).name == b.name

def transform(cfg, fasta, mzml, outdir):
    c = json.loads(json.dumps(cfg))
    c['database']['enzyme'].update(cleave_at='', restrict=None, c_terminal=True, semi_enzymatic=False)
    c['database']['fasta'] = fasta; c['quant']['lfq'] = False; c['mzml_paths'] = [mzml]
    c['output_directory'] = outdir; c.pop('output_paths', None); c.pop('version', None)
    return c

# control: the transformation reproduces the frozen A05 config of the first task exactly
t0 = old['tasks'][0]; s0 = t0['sample']
ctrl = transform(json.loads((CAP/'5_sage/razor_hash-acc'/s0/'results.json').read_text()), t0['fasta']['path'], t0['raw']['path'],
                 str(SRC/'results'/t0['id']/'search'))
assert ctrl == json.loads(pathlib.Path(t0['config']['path']).read_text()), 'config transformation does not reproduce frozen A05 config'

tasks = []
for _, r in p15.iterrows():
    s = r.sample_id; tid = s + '__unspecific'
    fasta = str(CAP/'4_parsimony/razor_hash-acc'/s/'parsimony.fasta'); mzml = str(RAW/(s + '.mzML'))
    src = CAP/'5_sage/razor_hash-acc'/s/'results.json'
    cfg = transform(json.loads(src.read_text()), fasta, mzml, str(RUN/'results'/tid/'search'))
    cp = RUN/'inputs'/(tid + '.json'); cp.write_text(json.dumps(cfg, indent=2) + '\n')
    tasks.append(dict(id=tid, sample=s, arm='unspecific', raw=meta(mzml, False), fasta=meta(fasta, False), source_config=meta(src),
                      config=meta(cp), baseline_search=meta(CAP/'5_sage/razor_hash-acc'/s/'results.sage.tsv', False),
                      semi_search=meta(CAP/'7_semitryptic/razor_hash-acc'/s/'results.sage.tsv', False),
                      design=dict(sample_id=s, cohort='CAPSCAN', subject_id=str(int(r.subject_id)), group=r.region, group_ordinal=str(int(r.region_ordinal)),
                                  aliquot='pellet', aliquot_main_id=str(int(r.workbook_Main)), aliquot_pellet_id=s.split('_T760_')[1], aliquot_super_id=str(int(r.workbook_Super)),
                                  subject_source='C2 fixed_roster_347.tsv (join fix)', group_source='C2 fixed_roster_347.tsv (join fix)'),
                      require_nontryptic=True))
m = dict(analysis='A05_4_7_unspecific__C2_participant15', register=old['register'], source_manifest=old['source_manifest'], tasks=tasks,
         sage=old['sage'], actual_version=old['actual_version'], frozen_code=[meta(p) for p in sorted((RUN/'frozen').glob('*.py'))],
         protocol=old['protocol'], copied_from=str(SRC), note='Frozen A05 code unchanged; 16 participant-15 tasks added by the C2 join fix')
(RUN/'inputs/manifest.json').write_text(json.dumps(m, indent=2) + '\n')
print('PREPARED', len(tasks))
