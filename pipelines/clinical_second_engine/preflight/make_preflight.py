#!/usr/bin/env python3
"""Build preflight RUN dirs: CAMPI reproduction (S02 FL_UNIPROT, S03 MG) from the
A07 manifest, and two negative controls (missing mzML; wrong FASTA sha256)."""
import copy, json, os, pathlib, sys
OUT = pathlib.Path(sys.argv[1]).resolve(); P = OUT / 'preflight'
sys.path.insert(0, str(OUT / 'frozen'))
from lib import meta, dump  # noqa
A07 = pathlib.Path(os.environ['FASTALAKE_DATA_DIR']) / 'MM_review_analysis/A07_4_2_second_engine'
src = json.loads((A07 / 'inputs/manifest.json').read_text())
frozen = [meta(OUT / 'frozen' / n) for n in ['search.py', 'lib.py', 'closed_fragger.params', 'peptide.py']]
assert [f['sha256'] for f in frozen] == [f['sha256'] for f in src['frozen_code']]
pick = [t for t in src['tasks'] if t['id'] in ('S02__FL_UNIPROT', 'S03__MG')]
assert len(pick) == 2
for t in pick:
    t['raw'] = meta(t['raw']['path'])  # adds sha256
    t['cohort'] = 'CAMPI'
def mk(name, tasks):
    """Write one preflight run directory with its manifest and a symlink to the frozen code."""
    run = P / name; (run / 'inputs').mkdir(parents=True, exist_ok=True)
    if not (run / 'frozen').exists(): os.symlink(OUT / 'frozen', run / 'frozen')
    m = {k: v for k, v in src.items() if k != 'tasks'}; m['frozen_code'] = frozen; m['tasks'] = tasks
    dump(run / 'inputs/manifest.json', m)
mk('campi_repro', pick)
neg1 = copy.deepcopy(pick[0]); neg1['raw']['path'] = neg1['raw']['path'].replace('S02.mzML', 'S02_DOES_NOT_EXIST.mzML'); neg1['id'] = 'NEG_missing_mzml'
mk('neg_missing_mzml', [neg1])
neg2 = copy.deepcopy(pick[0]); neg2['fasta']['sha256'] = '0' * 64; neg2['id'] = 'NEG_wrong_fasta_hash'
mk('neg_wrong_fasta_hash', [neg2])
print('PREFLIGHT MANIFESTS OK')
