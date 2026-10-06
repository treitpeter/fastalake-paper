from pathlib import Path
import os,sys,json,subprocess
O=Path(__file__).resolve().parent.parent
mode=sys.argv[1]
if mode=='pilot':
 idx=json.loads((O/'pilots.json').read_text())[int(os.environ['SLURM_ARRAY_TASK_ID'])]
 cmd=[sys.executable,str(O/'code/run_search.py'),str(idx),'--pilot']
elif mode=='standard':cmd=[sys.executable,str(O/'code/run_search.py'),os.environ['SLURM_ARRAY_TASK_ID']]
elif mode=='joint':cmd=[sys.executable,str(O/'code/run_search.py'),'88']
else:raise ValueError(mode)
subprocess.run(cmd,check=True)
