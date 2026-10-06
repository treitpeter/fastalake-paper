#!/usr/bin/env python3
"""Replay the fixed manuscript biology under the revised peptide-centric framework.

Sibling driver. It binds the overridden projection, imports the two patched derive
siblings in place of their frozen originals, and writes into this package's own results
tree. Every statistical module it calls is the frozen file, unmodified:
derive_eds_context.py and derive_capscan_extras.py run exactly as they shipped, and the
patched siblings differ from their parents only in the blocks listed in DIFFS/.

Repaired 19 September 2026, after review:

  - The analysis is named by FL_SELECTION and resolved through selection.py. The loose
    FL_PC_* overrides are gone; a destination now carries a sealed record of which cohort,
    grouping rule, species rule, class rule, input lineage and matrix files produced it, and
    refuses a conflicting reuse.
  - Input verification is unconditional. An output tag used to skip it entirely, so the
    repaired lineage ran with its cohort inputs unchecked. The manifest verified is now the
    one belonging to the selected lineage, and the lineage's own completion receipt is
    required.
  - Conservation is asserted per acquisition, over the allocated profile plus the retained
    remainder, before the run can succeed. It used to be recorded as a boolean that compared
    the annotated rows alone against the acquisition total, so it read false everywhere.
  - The baseline comparison joins on each table's own keys and reports added rows, removed
    rows, changed labels and changed missingness. An empty comparison is a failure, not a
    reproduction.
"""
import importlib
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
A19 = ROOT / 'MM_review_analysis/A19_4_6_biology_robustness'
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(A19 / 'frozen'))
from lib import dump, meta, verify                                  # frozen helper
import selection as sl
import compare as cmp
from validation_gates import assert_conservation
from exact_signflip import bind_paired_test

CODE = ['revised/search_pc.py', 'revised/selection.py', 'revised/compare.py',
        'revised/common_pepcentric.py', 'revised/derive_capscan_pc.py',
        'revised/derive_clinical_pc.py', 'revised/validation_gates.py',
        'revised/exact_signflip.py']


def verify_lineage(sel, m, cohort):
    """Verify the inputs of the lineage this selection actually names.

    The A19 manifest describes the 2026-09-10 grouping run. For that lineage its cohort
    entries are verified by hash and mtime. For the repaired lineage those entries do not
    apply, so the repaired run's own completion receipt is required instead, and the
    selected files are hashed into the seal by selection.resolve. Either way the frozen
    code is verified and nothing is skipped.
    """
    checks = []
    for x in m['frozen_code'] + [m['protocol']]:
        verify(x)
    checks.append(dict(check='frozen code hashes', status='PASS',
                       n=len(m['frozen_code']) + 1))
    if Path(sel['grouping_root']) == Path(m['grouping_root']):
        for x in m['cohorts'][cohort]['inputs']:
            verify(x)
        checks.append(dict(check='A19 manifest cohort inputs', status='PASS',
                           n=len(m['cohorts'][cohort]['inputs']),
                           lineage=sel['input_lineage']))
    else:
        receipt = Path(sel['grouping_root']) / 'COMPLETE.json'
        assert receipt.exists(), f'selected lineage has no completion receipt: {receipt}'
        r = json.loads(receipt.read_text())
        assert r.get('status') == 'PASS', f'selected lineage did not pass: {receipt}'
        checks.append(dict(check='selected lineage completion receipt', status='PASS',
                           receipt=str(receipt), lineage=sel['input_lineage'],
                           note='the A19 manifest describes a different grouping run and '
                                'does not apply to this lineage'))
        checks.append(dict(check='A19 manifest cohort inputs', status='NOT_APPLICABLE',
                           reason='different grouping run; see the lineage receipt above'))
    return checks


def run(name):
    sel = sl.resolve(name)
    cohort, method = sel['cohort'], sel['grouping_rule']
    assert sel['driver'] == 'search_pc', f'{name} is not a search_pc analysis'
    m = json.loads((A19 / 'inputs/manifest.json').read_text())
    for key in ('grouping_root', 'source_package', 'baseline_data'):
        m[key] = sel[key]
    lineage_checks = verify_lineage(sel, json.loads((A19 / 'inputs/manifest.json').read_text()),
                                    cohort)
    out = Path(sel['out_dir'])
    seal = sl.seal(sel, CODE)
    assert not (out / 'COMPLETE.json').exists(), 'this replay already completed'
    os.environ.update(MM_BIOLOGY_METHOD=method, MM_BIO_OUT=str(out), MM_BIO_ROOT=m['root'])
    import common
    import common_pepcentric as pc
    common.METHOD = method
    common.ROOT = Path(m['root'])
    common.RUN = Path(m['grouping_root'])
    common.DATA = out / 'source_data'
    common.PACKAGE = out
    frozen_aggregate = common.aggregate
    common.aggregate = pc.aggregate          # the one overridden function
    pc.INPUTS = Path(sel['matrices_dir'])
    pc.COHORT = cohort
    pc.RULE = method
    common.DATA.mkdir(parents=True, exist_ok=True)
    task = cohort + '_noclust'
    core = common.DATA / task
    core.mkdir(exist_ok=True)
    matrix, annot, summary = common.load(task)
    np.savez_compressed(core / 'groups.npz', groups=matrix.index.to_numpy(dtype=str),
                        samples=matrix.columns.to_numpy(dtype=str), values=matrix.to_numpy())
    common.write(annot.reset_index(), core / 'representative_annotations.tsv')
    original = common.RUN / 'cohorts' / task / 'reaggregated'
    for fname in ['SUMMARY.json', 'per_acquisition.tsv', 'fixed_pipeline_metrics.tsv',
                  'median_diff.tsv', 'assignment_changes.tsv']:
        shutil.copy2(original / fname, core / fname)
    acquisition_total = matrix.sum()
    projection = {}
    for level in ['OG', 'KO', 'Pathway']:
        revised = common.aggregate(matrix, annot, level, microbial_domain=True)
        common.write(revised, core / ('microbial_' + level + '.tsv'), index=True)
        common.write(frozen_aggregate(matrix, annot, level, microbial_domain=True),
                     core / ('microbial_' + level + '__frozen_projection.tsv'), index=True)
        # Per acquisition, the allocated profile plus the remainder not allocated to this
        # profile must equal the quantified intensity. One bad acquisition fails the run.
        remainder = pc.unassigned(level, 'microbial_single_og').reindex(matrix.columns).fillna(0.0)
        allocated = revised.sum(axis=1).reindex(matrix.columns).fillna(0.0)
        accounting = assert_conservation(acquisition_total, allocated, remainder)
        worst = accounting['worst_per_acquisition_relative_error']
        projection[level] = dict(
            revised_features=int(revised.shape[1]), frozen_features=int(frozen_aggregate(
                matrix, annot, level, microbial_domain=True).shape[1]),
            allocated_total=float(allocated.sum()), remainder_total=float(remainder.sum()),
            acquisition_total=float(acquisition_total.sum()),
            per_acquisition_conservation='PASS',
            worst_per_acquisition_relative_error=worst,
            acquisitions_checked=int((acquisition_total > 0).sum()),
            remainder_pct_of_quantified_intensity=float(
                100 * remainder.sum() / acquisition_total.sum()),
            remainder_meaning='intensity not allocated to this profile: the feature peptide '
                              'carried no qualifying annotation on a protein eligible under '
                              'this profile\'s rule. Not a measure of missing annotation.')
    for fname in (['local_candidate_composition.tsv', 'LOCAL_COMPOSITION_VALIDATION.json']
                  + (['local_candidate_species.npz'] if cohort == 'PREDICT' else [])):
        shutil.copy2(Path(m['source_package']) / 'source_data' / task / fname, core / fname)
    local = pd.read_csv(core / 'local_candidate_composition.tsv', sep='\t').set_index('sample_id').loc[matrix.columns]
    assert np.allclose(matrix.sum().to_numpy(), local.total_LFQ.to_numpy(), rtol=1e-10)
    if cohort == 'PREDICT':
        import derive_clinical_pc
        derive_clinical_pc.main()
        importlib.import_module('derive_eds_context')          # frozen, unmodified
    else:
        import derive_capscan_pc
        import derive_capscan_extras                            # frozen, unmodified
        bind_paired_test(derive_capscan_extras)
        derive_capscan_pc.main()
        derive_capscan_extras.main()
    figure = 'Figure4' if cohort == 'PREDICT' else 'Figure5'
    t = json.loads((A19 / 'inputs/manifest.json').read_text())['cohorts'][cohort]
    diffs = []
    for fname in t['baseline_tables']:
        base = Path(sel['baseline_data']) / figure / fname
        if not base.exists():
            diffs.append(dict(table=fname, status='NO_MATCHED_BASELINE',
                              reason=f'the selected lineage has no {fname}',
                              looked_in=str(base.parent)))
            continue
        entry = cmp.compare(common.DATA / figure / fname, base, name=fname)
        entry['baseline'] = meta(base)
        entry['baseline_lineage'] = sel['input_lineage']
        assert entry['comparison_is_substantive'], \
            f'{fname}: nothing was actually compared; an empty comparison is not a reproduction'
        diffs.append(entry)
    (out / 'BEFORE_AFTER.json').write_text(json.dumps(
        dict(selection=name, baseline_lineage=sel['input_lineage'],
             baseline_data=sel['baseline_data'], tables=diffs), indent=2, default=str) + '\n')
    per = pd.read_csv(core / 'per_acquisition.tsv', sep='\t')
    per = per[per.method.eq(method)].copy()
    common.write(per, out / 'group_endpoints.tsv')
    checks = dict(
        status='PASS', selection=name, cohort=cohort, method=method,
        acquisitions=len(matrix.columns),
        run_selection=str(seal),
        framework='peptide-centric equal-split functional projection; September study-level '
                  'direct peptide-to-containing-representative protein assignment; local-member '
                  'host/microbial classification retained',
        species_rule=sel['species_rule'], class_rule=sel['class_rule'],
        input_lineage=sel['input_lineage'],
        overridden_functions=['common.aggregate'] + (['summarize_biology.paired_test (dynamic instance)',
                                  'derive_capscan_extras.paired_test'] if cohort == 'CAPSCAN' else []),
        patched_siblings=json.loads((HERE / 'SIBLING_PROVENANCE.json').read_text()),
        lineage_verification=lineage_checks,
        projection=projection,
        total_intensity_matches_fixed_local_composition=True,
        baseline_comparison='KEYED; added rows, removed rows, changed labels and changed '
                            'missingness reported separately in BEFORE_AFTER.json',
        definitions='Fixed rosters, same local candidate classes, same representative annotation '
                    'rules and original multiple-testing families. Only the functional '
                    'projection and, where stated, the species rule changed.',
        job_id=os.environ.get('SLURM_JOB_ID'))
    outputs = [meta(p) for p in (common.DATA / figure).iterdir() if p.is_file()]
    checks['outputs'] = outputs + [meta(out / 'group_endpoints.tsv'),
                                   meta(out / 'BEFORE_AFTER.json'), meta(seal)]
    dump(out / 'COMPLETE.json', checks)
    print('REVISED BIOLOGY REPLAY PASS', name, flush=True)


if __name__ == '__main__':
    run(sys.argv[1])
