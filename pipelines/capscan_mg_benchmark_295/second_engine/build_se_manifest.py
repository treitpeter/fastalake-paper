#!/usr/bin/env python3
"""(e) second-engine (MSFragger 4.1 + Percolator 3.06.5) MG tasks for mg295.
inputs/manifest.json = the production second_engine_all manifest with every top-level key unchanged and
tasks replaced by: [0] the production CAPSCAN MG task of the control acquisition (copied verbatim),
[1..13] one MG task per participant-15 acquisition, built by the production rule (build_manifest.py /
append_p15.py): FASTA and source_config from that acquisition's MG Sage search (mg295/mg_full/out/<s>/results.json),
raw = the same mzML as that acquisition's production FL_public second-engine task (1610-1625).
FL_public second-engine searches for these 13 already exist (second_engine_all tasks 1610-1625); not rerun.
usage: build_se_manifest.py [sample ...]   (default: all 13; a subset for the pilot)"""
import csv, json, os, pathlib, sys
S = pathlib.Path('LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/second_engine_all'); RUN = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(S / 'frozen'))
from lib import meta, dump
m = json.loads((S / 'inputs/manifest.json').read_text()); assert len(m['tasks']) == 1626
CONTROL = '20210513_EXPL0_PeTr_SA_T760_1311'
ctl = next(t for t in m['tasks'] if t['cohort'] == 'CAPSCAN' and t['arm'] == 'MG' and t['sample'] == CONTROL)
fl = {t['sample']: t for t in m['tasks'] if t['cohort'] == 'CAPSCAN' and t['arm'] == 'FL_public'}
acq = list(csv.DictReader(open(RUN.parent / 'inputs/ACQ13.tsv'), delimiter='\t'))
only = sys.argv[1:]
new = []
for a in acq:
    s = a['sample']
    if only and s not in only: continue
    cfgp = pathlib.Path(a['mg_search_dir']) / 'results.json'
    cnt = json.loads((cfgp.parent / 'COUNTS.json').read_text()); assert cnt['sample'] == s
    cfg = json.loads(cfgp.read_text())
    fasta = os.path.realpath(cfg['database']['fasta']); assert fasta == os.path.realpath(a['mg_db']), (s, fasta)
    assert len(cfg['mzml_paths']) == 1 and os.path.realpath(cfg['mzml_paths'][0]) == os.path.realpath(fl[s]['raw']['path'])
    raw = meta(fl[s]['raw']['path']); assert raw == fl[s]['raw'], (s, 'mzML differs from the FL_public second-engine task')
    new.append(dict(id='CAPSCAN__' + s + '__MG', cohort='CAPSCAN', sample=s, arm='MG', fasta=meta(fasta), raw=raw, source_config=meta(cfgp)))
m2 = dict(m, tasks=[ctl] + new)
m2['mg295'] = dict(note='task 0 = production CAPSCAN MG task of the control acquisition (verbatim); tasks 1.. = participant-15 MG tasks', production_manifest=meta(S / 'inputs/manifest.json'))
dump(RUN / 'inputs/manifest.json', m2)
print('TASKS', len(m2['tasks']), [t['id'] for t in m2['tasks']])
