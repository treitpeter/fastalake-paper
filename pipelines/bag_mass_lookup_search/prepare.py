"""Freeze the native AlphaNovo/FastaLake experiment; no submission in this step."""
from pathlib import Path
from collections import Counter
from datetime import datetime, timezone
import csv
import hashlib
import json
import subprocess
import xml.etree.ElementTree as ET

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
REPO = ROOT / 'worktrees/bags'
OLD = HERE.parent / 'full_reservoir_2026-09-23'
PYTHON = ROOT.parent / 'AlphaForge/alphanovo/an_window_af3ac89b/.venv/bin/python'

def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

def meta(path):
    path = Path(path).resolve()
    return dict(path=str(path), bytes=path.stat().st_size, sha256=sha(path))

def dump(name, value):
    (HERE / name).write_text(json.dumps(value, indent=2) + '\n')

assert not (HERE / 'manifest.json').exists()
tests = ET.parse(HERE / 'FASTALAKE_TESTS.xml').getroot().find('testsuite')
assert tests.get('errors') == tests.get('failures') == '0'
assert int(tests.get('tests')) - int(tests.get('skipped')) >= 14
assert subprocess.check_output(['git', '-C', str(REPO), 'status', '--porcelain'], text=True).strip() == ''
revision = subprocess.check_output(['git', '-C', str(REPO), 'rev-parse', 'HEAD'], text=True).strip()
m = json.loads((OLD / 'manifest.json').read_text())
m['name'] = 'alphanovo_native_paper_rules_all_cohorts_2026-09-23'
m['output'] = str(HERE / 'outputs')
m['alphanovo_source'] = str(HERE / 'alphanovo_source')
m['alphanovo_revision'] = '93d647cf1d996864c3367d5f04e880b5df08d977'
m['comparison_branch_reviewed_not_submitted'] = '04ab8a4e6722969bae66ea470d703eda327b4b9c'
m['fasta_lake_revision'] = revision
m['python'] = str(PYTHON)
m['pilot_candidate_limit_seconds'] = 18 * 3600
m['paired_pilot_indices'] = [max((i for i,t in enumerate(m['tasks']) if t['cohort'] == c),
                               key=lambda i: m['tasks'][i]['predictions']['bytes'])
                              for c in m['expected_cohort_counts']]
m['pilot_selection'] = 'Largest archived prediction file in each cohort, without consulting new outcomes'
m['code'] = [meta(p) for p in sorted((REPO / 'fasta_lake').rglob('*.py'))]
m['code'] += [meta(REPO / p) for p in ['tools/run_bags.py', 'BAGS.md', 'rust/fasta_extractor_v2/src/bin/bag_extract.rs',
                                      'rust/fasta_extractor_v2/Cargo.lock', 'tests/test_bags.py']]
m['code'] += [meta(p) for p in sorted((HERE / 'alphanovo_source/alphanovo').rglob('*.py'))]
m['code'].append(meta(HERE / 'alphanovo_source/alphanovo/constants/vocabulary.yaml'))
m['bag_engine'] = meta(REPO / 'rust/fasta_extractor_v2/target/release/bag_extract')
m['sage'] = json.loads((OLD / 'search_manifest.json').read_text())['sage']
m['policy'] = dict(method='Exact, order bags, native two-substitution mass alternatives, union, two scrambled controls',
                   input='Archived rank-zero AlphaNovo; original FastaLake filters and top30% with ties',
                   below_cutoff_parsimony='Original exact support only',
                   mass_modification_policy='Disable D->N and E->Q per supplied paper metaproteome screen; native CAM alternatives separate',
                   candidate_space='Complete reservoir; all hit proteins; no baseline-context exclusion',
                   selection='One original IL-normalized prediction sequence per evidence unit; full graph fixed counts; FNV accession tie',
                   clinical_reservoirs='Assembly-inclusive; public-only benchmark remains separate',
                   error_controls='Scrambled-query retrieval diagnostics, not novel-subset FDP calibration',
                   speed='Joint scan time of all six comparisons; no per-arm speedup claim')
for obsolete in ['order_engine','frozen_files','validation','protocol','source_manifests']:
    m.pop(obsolete, None)
headers = Counter()
for t in m['tasks']:
    for key in ['predictions','baseline','search_config','archived_search']:
        rec = t[key]; stat = Path(rec['path']).stat()
        assert stat.st_size == rec['bytes'] and stat.st_mtime_ns == rec['mtime_ns'], (t['sample'], key)
    for rec in t['mzml']:
        stat = Path(rec['path']).stat()
        assert stat.st_size == rec['bytes'] and stat.st_mtime_ns == rec['mtime_ns']
    with open(t['predictions']['path']) as f:
        fields = next(csv.reader(f))
    assert {'sequence_with_bags','score','spec_idx','beam_rank','peptide_prediction_detokenized',
            'peptide_prediction_detokenized_unmodified','peptide_prediction'}.issubset(fields)
    headers[t['cohort']] += 1
assert dict(headers) == m['expected_cohort_counts']
for folder in ['logs', 'jobs', 'outputs']:
    (HERE / folder).mkdir(exist_ok=True)
m['created_utc'] = datetime.now(timezone.utc).isoformat()
dump('manifest.json', m)
dump('SOURCE_REVIEW.json', dict(source_paper=meta(ROOT / '_Nature_Biotech__AlphaNovo (1).docx'),
     original_roster=meta(OLD / 'manifest.json'), fasta_lake_revision=revision,
     alphanovo_revision=m['alphanovo_revision'], native_tests=dict(entrapment_mass_order=101, newer_mass_bag=59),
     fasta_lake_tests=dict(passed=14, skipped=10, reason='Optional directLFQ absent; this workflow uses sum'),
     inspected_paper_paragraphs=['94–106','130–141','157–160','181–185'],
     tests_receipt=meta(HERE / 'FASTALAKE_TESTS.xml')))
dump('PREFLIGHT.json', dict(status='PASS_INPUT_STAT_SCHEMA_CODE_TESTS', cohorts=dict(headers), tasks=len(m['tasks']),
                           pilot_indices=m['paired_pilot_indices'], code_files=len(m['code']),
                           reservoir_hash_check='Complete streamed hash and expected count in every candidate job'))
print(json.dumps(dict(tasks=len(m['tasks']),pilots=m['paired_pilot_indices'],revision=revision),indent=2))
