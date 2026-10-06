#!/usr/bin/env python3
"""Append the 16 CapScan participant-15 acquisitions to the second-engine manifest as
FastaLake-only tasks (arm FL_public, no MG pair) on their public-sequence FASTAs from
final/public347/A08_347. Same rule as build_manifest.py: FASTA and mzML are those in the
Sage results.json of that acquisition's public347 search.

  append_p15.py RUN                 apply: back up inputs/manifest.json once
                                    (inputs/manifest_1610_before_p15.json), append tasks
                                    1610-1625 atomically, append manifest.tsv rows
  append_p15.py RUN --dry OUT [S..] write the would-be manifest to OUT only (preflight),
                                    optionally for the listed samples only

Guards: the first 1610 tasks and every other manifest key stay identical; each
public347 sample must have COMPLETE.json status PASS; the mzML must be the one the
existing CAPSCAN FL_cohortMG task for that acquisition uses.
"""
import csv, json, os, pathlib, sys
RUN = pathlib.Path(sys.argv[1]).resolve()
sys.path.insert(0, str(RUN / 'frozen'))
from lib import meta, dump  # noqa: E402

ROOT = RUN.parent.parent
A08 = ROOT / 'final/public347/A08_347/results'
ROSTER347 = ROOT / 'coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/roster/fixed_roster_347.tsv'
N0 = 1610


def main():
    """Append the 16 participant-15 FL_public tasks to inputs/manifest.json (``--dry`` writes to a scratch inputs dir)."""
    dry = '--dry' in sys.argv
    mp = RUN / 'inputs/manifest.json'
    m = json.loads(mp.read_text())
    if len(m['tasks']) != N0:
        raise SystemExit(f'manifest has {len(m["tasks"])} tasks, expected {N0}; refusing to append twice')
    cap331 = {t['sample'] for t in m['tasks'] if t['cohort'] == 'CAPSCAN' and t['arm'] == 'FL_public'}
    with open(ROSTER347) as f:
        roster = [r['sample_id'] for r in csv.DictReader(f, delimiter='\t')]
    p15 = [s for s in roster if s not in cap331]
    assert len(cap331) == 331 and len(roster) == 347 and len(p15) == 16, (len(cap331), len(roster), len(p15))
    if dry and len(sys.argv) > 4:
        only = sys.argv[4:]
        assert set(only) <= set(p15), only
        p15 = [s for s in p15 if s in only]
    cohortmg_raw = {t['sample']: t['raw'] for t in m['tasks'] if t['cohort'] == 'CAPSCAN' and t['arm'] == 'FL_cohortMG'}
    new = []
    for s in p15:
        c = json.loads((A08 / s / 'COMPLETE.json').read_text())
        assert c['status'] == 'PASS' and c['sample'] == s, s
        cfgp = A08 / s / 'search/results.json'
        cfg = json.loads(cfgp.read_text())
        fasta = os.path.realpath(cfg['database']['fasta'])
        assert fasta == os.path.realpath(A08 / s / 'parsimony' / f'{s}_razor.fasta'), (s, fasta)
        assert len(cfg['mzml_paths']) == 1
        raw = os.path.realpath(cfg['mzml_paths'][0])
        assert pathlib.Path(raw).stem == s or pathlib.Path(raw).stem.startswith(s), (s, raw)
        rm = meta(raw)
        assert rm['path'] == cohortmg_raw[s]['path'] and rm['sha256'] == cohortmg_raw[s]['sha256'], (s, 'mzML differs from FL_cohortMG task')
        new.append(dict(id='CAPSCAN__' + s + '__FL_public', cohort='CAPSCAN', sample=s, arm='FL_public',
                        fasta=meta(fasta), raw=rm, source_config=meta(cfgp)))
    ids = {t['id'] for t in m['tasks']}
    assert not ids & {t['id'] for t in new}
    m2 = dict(m, tasks=m['tasks'] + new)
    m2['appended'] = dict(note='Tasks 1610-1625: CapScan participant-15 acquisitions, FastaLake-only (FL_public, no MG pair), '
                               'public-sequence FASTAs from final/public347/A08_347; appended by code/append_p15.py',
                          first_index=N0, n=len(new), roster=meta(ROSTER347))
    assert m2['tasks'][:N0] == m['tasks'] and {k: v for k, v in m2.items() if k not in ('tasks', 'appended')} == \
        {k: v for k, v in m.items() if k != 'tasks'}
    if dry:
        out = pathlib.Path(sys.argv[sys.argv.index('--dry') + 1])
        dump(out / 'manifest.json', m2)
        print('DRY', len(new), [t['id'] for t in new])
        return
    assert len(new) == 16
    bak = RUN / 'inputs/manifest_1610_before_p15.json'
    if not bak.exists():
        bak.write_bytes(mp.read_bytes())
    dump(mp, m2)
    with open(RUN / 'manifest.tsv', 'a') as f:
        w = csv.writer(f, delimiter='\t', lineterminator='\n')
        for i, t in enumerate(new, N0):
            for role in ['fasta', 'raw', 'source_config']:
                w.writerow([i, t['id'], t['cohort'], t['arm'], t['sample'], role, t[role]['path'], t[role]['bytes'], t[role]['sha256']])
    with open(RUN / 'MANIFEST_SUMMARY.md', 'a') as f:
        f.write(f'\nAppended {len(new)} tasks ({N0}-{N0 + len(new) - 1}): CAPSCAN FL_public for the participant-15 acquisitions, '
                'FastaLake-only, FASTA = `final/public347/A08_347/results/{s}/parsimony/{s}_razor.fasta` (code/append_p15.py). '
                'CAPSCAN FL_public is now 347; the paired MG comparison stays 282.\n')
    print('APPENDED', len(new), N0, N0 + len(new) - 1)


if __name__ == '__main__':
    main()
