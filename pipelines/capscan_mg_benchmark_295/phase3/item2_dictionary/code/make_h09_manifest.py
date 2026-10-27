#!/usr/bin/env python3
"""Write runs/<R>/h09/manifest.json in the format of Hartout_edits/H09_evidence_ablation/manifest.json (CAPSCAN task only).
Inputs are this run's grouping output (features) and its export (study_group_long, candidate_annotations, metaG peptide_compatibility)."""
import json, sys
from pathlib import Path
run = Path(sys.argv[1]).resolve(); n = int(sys.argv[2])
names = (run / 'inputs/acquisitions.txt').read_text().split(); assert len(names) == n
ex = run / 'export'
paths = dict(features=run / 'cohorts/FL/direct/feature_assignments.tsv.gz',
             original_long=ex / f'CAPSCAN_{n}_FL/study_group_long.tsv.gz',
             annotations=ex / f'CAPSCAN_{n}_FL/candidate_annotations.tsv.gz',
             mg_peptides=ex / f'CAPSCAN_{n}_metaG/peptide_compatibility.tsv.gz')
srcs = []
for p in paths.values():
    st = p.stat(); srcs.append(dict(path=str(p), bytes=st.st_size, mtime_ns=st.st_mtime_ns))
task = dict(cohort='CAPSCAN', samples=sorted(names), **{k: str(v) for k, v in paths.items()}, sources=srcs)
h = run / 'h09'; h.mkdir(exist_ok=True)
(h / 'manifest.json').write_text(json.dumps(dict(tasks=[task]), indent=1) + '\n')
print('h09 manifest', h / 'manifest.json')
