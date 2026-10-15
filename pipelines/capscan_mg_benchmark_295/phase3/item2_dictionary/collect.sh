#!/bin/bash
# Collection (small reads only): run after all jobs in JOBS.txt finished.
I=$(cd "$(dirname "$0")" && pwd)
sacct -j $(grep -oE '[0-9]{7}' $I/JOBS.txt | paste -sd,) --format=JobID,JobName%18,State,MaxRSS,Elapsed -P | grep -v extern
python3 -c "import json;d=json.load(open('$I/controls/CONTROL_282.json'));print('CONTROL_282',d['status'],d['n_checks'],[c['item'] for c in d['checks'] if not c['identical']])"
cat $I/controls/NEGATIVE_CONTROL.json
cat $I/runs/ctl282/values/medians_282.tsv $I/runs/r295/values/medians_295.tsv
python3 -c "import json;[print(r, json.dumps(json.load(open('$I/runs/'+r+'/comparison/SUMMARY.json'))['median_quantified_groups']), json.dumps(json.load(open('$I/runs/'+r+'/comparison/SUMMARY.json'))['comparisons'])) for r in ['ctl282','r295']]"
python3 $I/code/make_edits.py
