"""Extended task list for the job-6048085 Sage input-coverage audit.

The audit code is sage_input_integrity.py, copied byte-identical from
ROOT/work/sage_error/ (COPY_SHA256.txt). Only the task list (searches()) is
replaced; per-search reading, mzML checks, flag rule and outputs are unchanged.

Scopes (argv[1]):
  all      original task list, unchanged (reproduction of job 6048085)
  current  current search roots:
             primary_public_current  PREDICT 306 (PRIMARY_INPUTS, as in 6048085)
                                     + CAPSCAN 347 (final/public347/A08_347, fixed 347 roster)
             healthy_stool_58        canonical HSTOOL 5_sage, exact HSTOOL58 roster
             public_route_semi / public_route_unspecific   (as in 6048085; 917)
             public_route_unspecific_participant15         C2_participant15/fig6/A05 (16)
             canonical_assembly_inclusive CAPSCAN          (as in 6048085)
  plant    negative control: 2 real current searches + planted truncated inputs
           listed in argv[3] (JSON list of search dicts)
Output: argv[2] (TSVs + SUMMARY.json, as the original).
"""
import csv, glob, json, sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import sage_input_integrity as S  # noqa: E402

FL = S.FL
ROOT = FL / ('review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/'
             'biology_refit/claude_review_2026-10-03')
A08_347 = ROOT / 'final/public347/A08_347/results'
ROSTER347 = ROOT / 'coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/roster/fixed_roster_347.tsv'
A05 = ROOT / 'coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A05/results'
HSTOOL58 = FL / 'simplified_grouping/publication/analyses/simplified_grouping_2026-09-10/hstool_samples.txt'
ORIG = S.searches


def current():
    base = ORIG('all')
    out = [s for s in base if s['lineage'] == 'primary_public_current' and s['cohort'] == 'PREDICT']
    with open(ROSTER347) as f:
        cap = [r['sample_id'] for r in csv.DictReader(f, delimiter='\t')]
    assert len(cap) == 347 == len(set(cap))
    out += [dict(lineage='primary_public_current', cohort='CAPSCAN', sample=s, repaired=False,
                 dir=str(A08_347 / s / 'search')) for s in cap]
    hs = [l.strip() for l in open(HSTOOL58) if l.strip()]
    assert len(hs) == 58 == len(set(hs))
    canon_hs = {s['sample']: s for s in base if s['lineage'] == 'canonical_assembly_inclusive' and s['cohort'] == 'HSTOOL'}
    out += [dict(canon_hs[s], lineage='healthy_stool_58') for s in hs]
    out += [s for s in base if s['lineage'].startswith('public_route_')]
    for rj in sorted(glob.glob(str(A05 / '*__unspecific/search/results.json'))):
        d = Path(rj).parent
        out.append(dict(lineage='public_route_unspecific_participant15', cohort='CAPSCAN',
                        sample=d.parent.name.split('__')[0], repaired=None, dir=str(d)))
    out += [s for s in base if s['lineage'] == 'canonical_assembly_inclusive' and s['cohort'] == 'CAPSCAN']
    return out


def searches(scope):
    if scope == 'current':
        return current()
    if scope == 'plant':
        real = [s for s in current() if s['lineage'] == 'healthy_stool_58'][:2]
        return real + json.loads(Path(sys.argv[3]).read_text())
    return ORIG(scope)


S.searches = searches

if __name__ == '__main__':
    S.main()
