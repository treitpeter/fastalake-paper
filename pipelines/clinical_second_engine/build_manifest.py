#!/usr/bin/env python3
"""Stage the second-engine inputs for PREDICT, CapScan and healthy stool.

Every FASTA and mzML is taken from the Sage results.json of the search that the
main benchmark used for that acquisition and arm (database.fasta, mzml_paths).
Writes inputs/manifest.json (consumed unchanged by frozen/search.py from the
CAMPI run A07_4_2_second_engine), manifest.tsv and MANIFEST_SUMMARY.md.
Run on Slurm (hashes ~1 TB).
"""
import csv, json, os, sys, pathlib, collections
from multiprocessing import Pool
OUT = pathlib.Path(sys.argv[1]).resolve()
sys.path.insert(0, str(OUT / 'frozen'))
from lib import meta, sha, dump  # noqa: E402

FL = os.environ['FASTALAKE_DATA_DIR']  # LAB_ROOT layout
ARMS = FL + '/canonical/publication_data/arms'
SD = (FL + '/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/'
      'biology_refit/claude_review_2026-10-03/final/figure_sd/source_data/package/FastaLake_Source_Data')
A07 = FL + '/MM_review_analysis/A07_4_2_second_engine'
A08 = FL + '/MM_review_analysis/A08_4_5_capscan_public_only/results'
ROSTER347 = (FL + '/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/'
             'biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/roster/fixed_roster_347.tsv')
HSTOOL58 = FL + '/simplified_grouping/publication/analyses/simplified_grouping_2026-09-10/hstool_samples.txt'

CONFIG = {
    ('PREDICT', 'FL_public'): ARMS + '/abl_nometag_noclust/cohorts/PREDICT/5_sage/razor_hash-acc/{s}/results.json',
    ('PREDICT', 'MG'): ARMS + '/metag_hashed_human/PREDICT/{s}/results.json',
    ('CAPSCAN', 'FL_public'): A08 + '/{s}/search/results.json',
    ('CAPSCAN', 'FL_cohortMG'): ARMS + '/noclust_2026-08-07/cohorts/CAPSCAN/5_sage/razor_hash-acc/{s}/results.json',
    ('CAPSCAN', 'MG'): ARMS + '/metag_hashed_human/CAPSCAN/{s}/results.json',
    ('HSTOOL', 'FL'): ARMS + '/noclust_2026-08-07/cohorts/HSTOOL/5_sage/razor_hash-acc/{s}/results.json',
}
EXPECT_FASTA = {
    ('PREDICT', 'FL_public'): ARMS + '/abl_nometag_noclust/cohorts/PREDICT/4_parsimony/razor_hash-acc/{s}/parsimony.fasta',
    ('PREDICT', 'MG'): ARMS + '/metag_hashed_human/PREDICT/_db/{s}.metag_hashed_plus_human.fasta',
    ('CAPSCAN', 'FL_public'): A08 + '/{s}/parsimony/{s}_razor.fasta',
    ('CAPSCAN', 'FL_cohortMG'): ARMS + '/noclust_2026-08-07/cohorts/CAPSCAN/4_parsimony/razor_hash-acc/{s}/parsimony.fasta',
    ('CAPSCAN', 'MG'): ARMS + '/metag_hashed_human/CAPSCAN/_db/{s}.metag_hashed_plus_human.fasta',
    ('HSTOOL', 'FL'): ARMS + '/noclust_2026-08-07/cohorts/HSTOOL/4_parsimony/razor_hash-acc/{s}/parsimony.fasta',
}


def rows(path):
    """Read a tab-separated file as a list of dicts."""
    with open(path) as f:
        return list(csv.DictReader(f, delimiter='\t'))


def bench_roster(fname, cohort, test, ref):
    """Sorted acquisition roster of one Fig. 4 benchmark comparison (cohort, test arm, reference arm, min length 5)."""
    return sorted({r['sample'] for r in rows(SD + '/Figure4/' + fname)
                   if r['cohort'] == cohort and r['test_arm'] == test and r['reference_arm'] == ref and r['min_length'] == '5'})


def main():
    """Build the task list (acquisition x arm), hash every FASTA and mzML, write the manifests."""
    predict296 = bench_roster('benchmark_0_per_acquisition.tsv', 'PREDICT', 'FL_public', 'MG')
    cap_fl331 = bench_roster('benchmark_1_per_acquisition.tsv', 'CAPSCAN', 'FL_public', 'FL_cohortMG')
    cap_mg282 = bench_roster('benchmark_1_per_acquisition.tsv', 'CAPSCAN', 'FL_public', 'MG')
    cap347 = [r['sample_id'] for r in rows(ROSTER347)]
    hs58 = [l.strip() for l in open(HSTOOL58) if l.strip()]
    assert (len(predict296), len(cap_fl331), len(cap_mg282), len(cap347), len(hs58)) == (296, 331, 282, 347, 58)
    assert len(set(cap347)) == 347 and set(cap_fl331) <= set(cap347) and set(cap_mg282) <= set(cap_fl331)
    # Every CapScan roster acquisition with a matched MG database (benchmark directory set)
    mg_dirs = {p.name for p in pathlib.Path(ARMS + '/metag_hashed_human/CAPSCAN').iterdir() if p.is_dir() and p.name != '_db'}
    cap_mg = sorted(s for s in cap347 if s in mg_dirs and os.path.exists(CONFIG['CAPSCAN', 'MG'].format(s=s)))
    assert cap_mg == cap_mg282, (len(cap_mg), sorted(set(cap_mg) ^ set(cap_mg282)))
    plan = ([('PREDICT', 'FL_public', s) for s in predict296] + [('PREDICT', 'MG', s) for s in predict296]
            + [('CAPSCAN', 'FL_public', s) for s in cap_fl331] + [('CAPSCAN', 'FL_cohortMG', s) for s in cap347]
            + [('CAPSCAN', 'MG', s) for s in cap_mg282] + [('HSTOOL', 'FL', s) for s in hs58])
    staged = []
    for cohort, arm, s in plan:
        cfgp = CONFIG[cohort, arm].format(s=s)
        c = json.load(open(cfgp))
        fasta = os.path.realpath(c['database']['fasta'])
        assert fasta == os.path.realpath(EXPECT_FASTA[cohort, arm].format(s=s)), (cohort, arm, s, fasta)
        assert len(c['mzml_paths']) == 1
        raw = os.path.realpath(c['mzml_paths'][0])
        assert pathlib.Path(raw).stem == s or pathlib.Path(raw).stem.startswith(s), (s, raw)
        staged.append((cohort, arm, s, os.path.realpath(cfgp), fasta, raw))
    # one raw per acquisition across arms
    byacq = collections.defaultdict(set)
    for cohort, arm, s, cfgp, fasta, raw in staged:
        byacq[cohort, s].add(raw)
    assert all(len(v) == 1 for v in byacq.values()), [k for k, v in byacq.items() if len(v) > 1][:5]
    files = sorted({x for t in staged for x in t[3:]})
    with Pool(int(os.environ.get('SLURM_CPUS_PER_TASK', '8'))) as pool:
        metas = dict(zip(files, pool.map(meta, files, chunksize=1)))
    # Benchmark-comparator cross-check: PREDICT MG FASTA sha256 recorded by the Fig. 4d analysis
    summ = json.load(open(SD + '/Figure4/SUMMARY.json'))
    mg_checked = 0
    for cohort, arm, s, cfgp, fasta, raw in staged:
        if cohort == 'PREDICT' and arm == 'MG':
            assert summ['input_hashes'][s]['MG_fasta_sha256'] == metas[fasta]['sha256'], s
            mg_checked += 1
    assert mg_checked == 296
    src = json.load(open(A07 + '/inputs/manifest.json'))
    tasks = [dict(id=cohort + '__' + s + '__' + arm, cohort=cohort, sample=s, arm=arm,
                  fasta=metas[fasta], raw=metas[raw], source_config=metas[cfgp])
             for cohort, arm, s, cfgp, fasta, raw in staged]
    assert len({t['id'] for t in tasks}) == len(tasks)
    frozen = [meta(OUT / 'frozen' / n) for n in ['search.py', 'lib.py', 'closed_fragger.params', 'peptide.py']]
    for a, b in zip(frozen, src['frozen_code']):
        assert a['sha256'] == b['sha256'], (a['path'], 'differs from CAMPI A07 frozen code')
    m = dict(analysis='second_engine_all', source_pipeline=A07, source_manifest=meta(A07 + '/inputs/manifest.json'),
             frozen_code=frozen, jar=src['jar'], percolator=src['percolator'], java=src['java'], seed=src['seed'],
             revision=src['revision'], created_utc=__import__('datetime').datetime.utcnow().isoformat() + 'Z',
             rosters=dict(PREDICT_benchmark=SD + '/Figure4/benchmark_0_per_acquisition.tsv',
                          CAPSCAN_benchmark=SD + '/Figure4/benchmark_1_per_acquisition.tsv',
                          CAPSCAN_347=meta(ROSTER347), HSTOOL_58=meta(HSTOOL58)),
             tasks=tasks)
    for k in ['jar', 'percolator', 'java']:
        assert sha(m[k]['path']) == m[k]['sha256'], k
    dump(OUT / 'inputs/manifest.json', m)
    with open(OUT / 'manifest.tsv', 'w') as f:
        w = csv.writer(f, delimiter='\t', lineterminator='\n')
        w.writerow(['index', 'id', 'cohort', 'arm', 'sample', 'role', 'path', 'bytes', 'sha256'])
        for i, t in enumerate(tasks):
            for role in ['fasta', 'raw', 'source_config']:
                w.writerow([i, t['id'], t['cohort'], t['arm'], t['sample'], role, t[role]['path'], t[role]['bytes'], t[role]['sha256']])
    cnt = collections.Counter((t['cohort'], t['arm']) for t in tasks)
    L = ['# Second-engine manifest summary', '',
         'Inputs per task are the FASTA and mzML recorded in the Sage results.json of the main-benchmark search for that acquisition and arm; sha256 of every file in manifest.tsv.', '',
         '| cohort | arm | acquisitions | FASTA source |', '|---|---|---:|---|']
    for (c, a), n in sorted(cnt.items()):
        L.append(f'| {c} | {a} | {n} | `{EXPECT_FASTA[c, a].replace(FL, "FastaLake")}` |')
    L += ['', f'Total tasks: {len(tasks)}. Unique mzML files: {len({t["raw"]["path"] for t in tasks})}.',
          f'Paired FL/MG acquisitions: PREDICT {len(predict296)}, CapScan {len(cap_mg282)} (of the 347-acquisition roster; the other 65 have no specimen-matched MG database, including all 16 participant-15 acquisitions).',
          'PREDICT MG FASTA sha256 equal to the Fig. 4d benchmark record for all 296.']
    (OUT / 'MANIFEST_SUMMARY.md').write_text('\n'.join(L) + '\n')
    print('MANIFEST', len(tasks), dict(cnt))


if __name__ == '__main__':
    main()
