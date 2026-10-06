#!/usr/bin/env python3
"""Healthy-stool replay: Fig 3f, Supplementary 9 and 10, and the leader analyses.

Three modes.

  baseline    runs the UNPATCHED published producer through this driver and requires it to
              reproduce the published tables. This gates the other two: if the plumbing
              cannot reproduce what shipped, a peptide-centric number from the same
              plumbing means nothing.
  simplified  the patched sibling, species labels from the simplified grouping rule
  original    the patched sibling, species labels from the original grouping rule

The published producer is read, never modified. Supplementary 10's own statistics are not
touched at all: its computation consumes a per-acquisition leader table, so it is re-run
afterwards from the frozen h02_technical_leaders.py with these tables as input.
"""
import json
import os
import shutil
import sys
from pathlib import Path
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
PACKAGE = HERE.parent
ROOT = PACKAGE.parents[1]
MSCODE = ROOT / 'release_prep/integration_2026-09-14/repo/manuscript/current/code'
PUBLISHED = ROOT / 'release_prep/integration_2026-09-14/repo/manuscript/current/source_data/Figure3'
GROUPING = ROOT / 'simplified_grouping_results/manuscript_cohorts_2026-09-10_2200'
LOCALPKG = ROOT / 'PeTr3_simplified_grouping_2026-09-11/source_data/HSTOOL_noclust'
ARM = ROOT / 'canonical/publication_data/arms/noclust_2026-08-07'
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(MSCODE))
import selection as sl
import compare as cmp

# Every table the producer must emit. A table that fails to appear is a failure, which the
# earlier loop could not express because it iterated whatever happened to exist.
TABLES = ['pair_summary_local_candidate.tsv', 'pyramid_pairwise_local_candidate.tsv',
          'per_function_acquisition.tsv', 'function_tier_summary.tsv', 'proposed_top12.tsv']
CODE = ['revised/hstool_pc.py', 'revised/selection.py', 'revised/compare.py',
        'revised/common_pepcentric.py', 'revised/derive_hstool_pc.py']


def _superseded_compare(new, old, keys):
    a = pd.read_csv(new, sep='\t')
    b = pd.read_csv(old, sep='\t')
    keys = [k for k in keys if k in a.columns and k in b.columns]
    entry = dict(table=Path(new).name, revised_rows=len(a), published_rows=len(b), keys=keys)
    if not keys:
        return entry
    j = a.merge(b, on=keys, how='outer', suffixes=('_revised', '_published'), indicator=True)
    both = j['_merge'].eq('both')
    entry.update(matched=int(both.sum()), revised_only=int(j['_merge'].eq('left_only').sum()),
                 published_only=int(j['_merge'].eq('right_only').sum()), columns={})
    for c in a.columns:
        if c in keys or not pd.api.types.is_numeric_dtype(a[c]):
            continue
        if f'{c}_revised' not in j or f'{c}_published' not in j:
            continue
        x = j.loc[both, f'{c}_revised'].to_numpy(float)
        y = j.loc[both, f'{c}_published'].to_numpy(float)
        ok = np.isfinite(x) & np.isfinite(y)
        if not ok.any():
            continue
        d = np.abs(x[ok] - y[ok])
        entry['columns'][c] = dict(
            n=int(ok.sum()), max_abs=float(d.max()),
            median_abs=float(np.median(d)),
            identical=bool(np.allclose(x[ok], y[ok], rtol=1e-8, atol=1e-10)),
            revised_median=float(np.median(x[ok])), published_median=float(np.median(y[ok])))
    return entry


def run(name):
    sel = sl.resolve(name)
    mode = sel['mode']
    assert sel['driver'] == 'hstool_pc', f'{name} is not an hstool_pc analysis'
    assert mode in ('baseline', 'simplified', 'original'), mode
    out = Path(sel['out_dir'])
    seal = sl.seal(sel, CODE)
    assert not (out / 'COMPLETE.json').exists(), 'this replay already completed'
    data = out / 'source_data'
    core = data / 'HSTOOL_noclust'
    core.mkdir(parents=True, exist_ok=True)
    shutil.copy2(GROUPING / 'cohorts/HSTOOL_noclust/reaggregated/per_acquisition.tsv',
                 core / 'per_acquisition.tsv')
    shutil.copy2(LOCALPKG / 'local_candidate_group_intensities.npz',
                 core / 'local_candidate_group_intensities.npz')
    import common
    import common_pepcentric as pc
    common.ROOT = ROOT
    common.RUN = GROUPING
    common.DATA = data
    common.ARM = ARM
    pc.INPUTS = Path(sel['matrices_dir'])
    pc.COHORT = 'HSTOOL'
    pc.RULE = 'simplified' if mode == 'baseline' else mode
    if mode == 'baseline':
        import derive_hstool as producer                      # published, unmodified
    else:
        import derive_hstool_pc as producer                   # patched sibling
    producer.main()
    figure = data / 'Figure3'
    missing = [t for t in TABLES if not (figure / t).exists()]
    assert not missing, f'the producer did not emit {missing}'
    diffs = []
    for tname in TABLES:
        # The published package gzips the leader table and leaves the others plain. The
        # earlier gate looked only for the plain name, so per_function_acquisition.tsv was
        # silently skipped and the reproduction claim never covered the table S09 and S10
        # depend on.
        candidates = [PUBLISHED / tname, PUBLISHED / (tname + '.gz')]
        base = next((c for c in candidates if c.exists()), None)
        if base is None:
            diffs.append(dict(table=tname, status='NO_MATCHED_BASELINE',
                              looked_for=[str(c) for c in candidates]))
            continue
        entry = cmp.compare(figure / tname, base, name=tname)
        entry['baseline'] = str(base)
        diffs.append(entry)
    (out / 'BEFORE_AFTER.json').write_text(json.dumps(
        dict(selection=name, tables=diffs), indent=2, default=str) + chr(10))
    gate = None
    if mode == 'baseline':
        # Every published table must come back identical, or nothing downstream is trusted.
        bad = []
        for e in diffs:
            if e.get('status') == 'NO_MATCHED_BASELINE':
                bad.append(f"{e['table']}: no published baseline to reproduce")
                continue
            if not e.get('comparison_is_substantive'):
                bad.append(f"{e['table']}: nothing was actually compared")
            if e.get('added_rows') or e.get('removed_rows'):
                bad.append(f"{e['table']}: {len(e['added_rows'])} added, "
                           f"{len(e['removed_rows'])} removed rows")
            for c, v in e.get('changed_labels', {}).items():
                bad.append(f"{e['table']}:{c} {v['n']} changed labels")
            for c, n in e.get('changed_missingness', {}).items():
                bad.append(f"{e['table']}:{c} missingness changed in {n} rows")
            for c, v in e.get('columns', {}).items():
                if 'identical' in v and not v['identical']:
                    bad.append(f"{e['table']}:{c} max_abs={v['max_abs']:.3g}")
                elif v.get('status'):
                    bad.append(f"{e['table']}:{c} {v['status']}")
        gate = dict(reproduced=not bad, failures=bad, tables_required=TABLES)
        assert not bad, 'BASELINE REPRODUCTION FAILED: ' + '; '.join(bad[:8])
    else:
        # A dependent run needs the gate to have actually passed, not merely to have run.
        greceipt = PACKAGE / 'results' / ('HSTOOL__baseline' if sel['role'] == 'primary'
                   else 'sensitivity_ungated/HSTOOL__baseline') / 'COMPLETE.json'
        assert greceipt.exists(), f'the reproduction gate has not run: {greceipt}'
        g = json.loads(greceipt.read_text())
        assert g.get('status') == 'PASS' and g.get('baseline_gate', {}).get('reproduced'), \
            f'the reproduction gate did not pass: {greceipt}'
        gate = dict(depends_on=str(greceipt), gate_reproduced=True)
    (out / 'COMPLETE.json').write_text(json.dumps(dict(
        status='PASS', selection=name, mode=mode, cohort='HSTOOL',
        species_rule=sel['species_rule'], class_rule=sel['class_rule'],
        input_lineage=sel['input_lineage'], run_selection=str(seal),
        projection='group-centric full-intensity (published)' if mode == 'baseline'
                   else 'peptide-centric equal split',
        species_label_rule=pc.RULE,
        baseline_gate=gate, tables=TABLES,
        job_id=os.environ.get('SLURM_JOB_ID')), indent=2, default=float) + '\n')
    print('HSTOOL REPLAY PASS', name, flush=True)


if __name__ == '__main__':
    run(sys.argv[1])
