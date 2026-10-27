#!/usr/bin/env python3
"""(d) MG_ODD for the spectrum-split comparison (Supp. Fig. 6d,e), via the production iteration3/disjoint
worker.py copied byte-identically. Only the roster differs (FL_ROSTER=roster_mg295.tsv: production 578 rows
+ 13 rows i=578..590). worker.py's own __main__ always runs ['odd','full']; 'full' needs archived frozen MG
sets, which the 13 do not have (their same-spectrum MG search is job (b), mg_full/). So this driver calls
worker.main(i, ['odd']) for the 13 and worker.main(i, ['odd','full']) for the control row (335).
usage: run_disjoint.py <roster i>   or   run_disjoint.py <list file> (SLURM_ARRAY_TASK_ID indexes it)"""
import os, sys
from pathlib import Path
os.environ['FL_ROSTER'] = 'roster_mg295.tsv'
sys.path.insert(0, str(Path(__file__).resolve().parent))
import worker
a = sys.argv[1]
i = int(a) if a.isdigit() else int(open(Path(__file__).resolve().parent / a).read().split()[int(os.environ['SLURM_ARRAY_TASK_ID'])])
worker.main(i, ['odd', 'full'] if i < 578 else ['odd'])
