"""Read the frozen cohort with only explicitly recorded, checksum-pinned exclusions."""
import csv
import hashlib
import json
from pathlib import Path

BASE = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def table(path):
    with Path(path).open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def selected_runs(base=BASE):
    folder = Path(base) / 'manifests'
    original = table(folder / 'runs.tsv')
    if not (folder / 'selection.json').exists():
        return original
    selection = json.loads((folder / 'selection.json').read_text())
    for name in ['runs.tsv', 'active_runs.tsv', 'excluded_runs.tsv']:
        if sha(folder / name) != selection['file_sha256'][name]:
            raise RuntimeError('Changed selection input: ' + name)
    active = table(folder / 'active_runs.tsv')
    excluded = table(folder / 'excluded_runs.tsv')
    names = [r['run_accession'] for r in original]
    excluded_names = set(selection['excluded_runs'])
    if len(set(names)) != len(names) or not excluded_names.issubset(names):
        raise RuntimeError('Invalid original or excluded run identities')
    expected = [r for r in original if r['run_accession'] not in excluded_names]
    if active != expected or len(active) != selection['active_run_count'] or len(original) != selection['original_run_count']:
        raise RuntimeError('Active cohort differs from the authorized selection')
    if {r['run_accession'] for r in excluded} != excluded_names or len(excluded) != len(excluded_names):
        raise RuntimeError('Exclusion table differs from the authorized selection')
    if not all(r.get('reason') for r in excluded) or not selection.get('authorization'):
        raise RuntimeError('Missing exclusion reason or authorization')
    return active


def provenance(base=BASE):
    rows = selected_runs(base)
    folder = Path(base) / 'manifests'
    if not (folder / 'selection.json').exists():
        return dict(active_run_count=len(rows), original_run_count=len(rows), excluded_runs=[],
                    active_manifest_sha256=sha(folder / 'runs.tsv'), selection_sha256=None)
    value = json.loads((folder / 'selection.json').read_text())
    return dict(active_run_count=len(rows), original_run_count=value['original_run_count'],
                excluded_runs=value['excluded_runs'], matched_acquisitions_removed=value['matched_acquisitions_removed'],
                matched_specimens_removed=value['matched_specimens_removed'],
                active_manifest_sha256=sha(folder / 'active_runs.tsv'), selection_sha256=sha(folder / 'selection.json'))


def require_input_gate(base=BASE):
    state = json.loads((Path(base) / 'ANALYSIS_INPUTS_COMPLETE.json').read_text())
    if state.get('state') != 'complete' or state.get('selection') != provenance(base):
        raise RuntimeError('Missing or stale verification of selected analysis inputs')
    return state
