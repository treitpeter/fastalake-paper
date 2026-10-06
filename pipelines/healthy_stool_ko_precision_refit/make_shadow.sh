#!/bin/bash
# Build a read-only shadow of biology_refit (B) and the hstool driver (C) whose profile inputs point at
# W/profiles/<variant>. Analysis scripts are byte-identical copies; only common.py path roots and the
# hstool H/NEW/out roots are redirected. Frozen baselines and groups are symlinked (read-only use).
set -euo pipefail
V=$1
W=${FASTALAKE_DATA_DIR}/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native
CR=${FASTALAKE_DATA_DIR}/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01
BR=$CR/biology_refit
S=$W/downstream/$V; mkdir -p $S/C $S/B
ln -sfn $W/profiles/$V $S/C/profiles
ln -sfn $CR/clinical_refit $S/C/clinical_refit
for f in crossomics.py capscan_functions.py eds_biochemistry.py nested_eds.py individual_functions.py build_groups.py run.sh PROTOCOL.txt CODE.json INPUTS.json CROSS_INPUTS.json CROSS_CODE.json CAPSCAN_INPUTS.json CAPSCAN_CODE.json EDS_INPUTS.json EDS_CODE.json NESTED_INPUTS.json NESTED_CODE.json INDIVIDUAL_INPUTS.json INDIVIDUAL_CODE.json; do cp -p $BR/$f $S/B/$f; done
for d in groups cross_baseline capscan_baseline eds_baseline nested_baseline; do ln -sfn $BR/$d $S/B/$d; done
python3 - "$S" "$BR" "$CR" <<'PY'
import sys,re
S,BR,CR=sys.argv[1:]
src=open(BR+'/common.py').read()
old_head="H=Path(__file__).resolve().parent; C=H.parent; J=C.parents[1]"
assert old_head in src
new_head=("H=Path(__file__).resolve().parent; C=H.parent/'C'; J=Path('%s').parents[1]\n"
 "# SHADOW: C/profiles -> corrected profiles; every other root is the real frozen tree.\n"
 "_CREAL=Path('%s')") % (CR,CR)
src=src.replace(old_head,new_head)
for a,b in [("O=C.parent/","O=_CREAL.parent/"),("E=C.parent/","E=_CREAL.parent/")]:
    assert a in src;src=src.replace(a,b)
old_verify=" for p,h in load(H/'CODE.json').items():assert sha(H/p)==h,p\n for p,h in load(H/'INPUTS.json').items():assert sha(p)==h,p"
assert old_verify in src
src=src.replace(old_verify," for p,h in load(H/'CODE.json').items():\n  if p!='common.py':assert sha(H/p)==h,p\n for p,h in load(H/'INPUTS.json').items():\n  if '/manuscript_completion_2026-10-01/profiles/' not in p:assert sha(p)==h,p")
open(S+'/B/common.py','w').write(src)
h=open(CR+'/hstool_precision.py').read()
a="H=Path(__file__).resolve().parent\n";assert a in h
h=h.replace(a,"H=Path('%s')\nSHADOW=Path(__file__).resolve().parent\n"%CR)
a="NEW=H/'profiles/HSTOOL_noclust'";assert a in h;h=h.replace(a,"NEW=SHADOW/'C/profiles/HSTOOL_noclust'")
a="out=H/('hstool_precision_'+mode)";assert a in h;h=h.replace(a,"out=SHADOW/('hstool_precision_'+mode)")
open(S+'/hstool_precision.py','w').write(h)
PY
cd $S && diff <(sed 's/^/ /' $BR/common.py) <(sed 's/^/ /' B/common.py) > common.diff || true
diff $CR/hstool_precision.py hstool_precision.py > hstool.diff || true
echo shadow $S
