#!/usr/bin/env python3
"""Pick benchmark tasks: per cohort x arm, the task with the largest FASTA and the task with the largest mzML."""
import json, sys, pathlib
RUN = pathlib.Path(sys.argv[1]); m = json.loads((RUN / 'inputs/manifest.json').read_text())
groups = {}
for i, t in enumerate(m['tasks']): groups.setdefault((t['cohort'], t['arm']), []).append(i)
pick = []
for k, ix in sorted(groups.items()):
    for key in ('fasta', 'raw'):
        j = max(ix, key=lambda i: m['tasks'][i][key]['bytes'])
        if j not in pick: pick.append(j)
        print(k[0], k[1], 'max', key, j, m['tasks'][j]['sample'], m['tasks'][j]['fasta']['bytes'], m['tasks'][j]['raw']['bytes'], file=sys.stderr)
print(','.join(map(str, pick)))
