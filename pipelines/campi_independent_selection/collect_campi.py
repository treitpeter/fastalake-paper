"""Recount the common-reference CAMPI comparison and retain all diagnostics."""
import csv
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent
PROJECT = ROOT.parents[1]
sys.path.insert(0, str(PROJECT / 'canonical/publication/lib'))
from peptide import peptide_set


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_table(path, rows):
    with path.open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def main():
    base = ROOT / 'campi_ref_v1'
    output = base / 'review'
    output.mkdir(exist_ok=True)
    ordinary = base / 'ordinary'
    assert json.loads((ordinary / 'COMPLETE.json').read_text())['searches'] == 6
    matrix = json.loads((ordinary / 'SEARCH_MATRIX.json').read_text())
    sets, rows = {}, []
    for item in matrix:
        assert sha(item['result']) == item['result_sha256']
        sets[item['evaluated'], item['arm']] = peptide_set(item['result'])
    for item in matrix:
        sample, arm = item['evaluated'], item['arm']
        peptides = sets[sample, arm]
        reference = sets[sample, 'fixed']
        for low, high in ((5, 50), (9, 50), (5, 8)):
            retained = {p for p in peptides if low <= len(p) <= high}
            fixed = {p for p in reference if low <= len(p) <= high}
            rows.append(dict(evaluated=sample, arm=arm, min_residues=low, max_residues=high,
                             canonical_target_peptides=len(retained), recovered_fixed=len(retained & fixed),
                             lost_fixed=len(fixed - retained), additional_vs_fixed=len(retained - fixed),
                             q_threshold=.01, unit='canonical_peptides_at_peptide_q',
                             interpretation='yield_and_identity_not_truth'))
    write_table(output / 'ordinary_counts.tsv', rows)
    diagnostics = []
    for seed in (20260908, 20260909):
        directory = base / str(seed)
        assert json.loads((directory / 'COMPLETE.json').read_text())['searches'] == 6
        with (directory / 'diagnostics.tsv').open() as stream:
            diagnostics.extend(csv.DictReader(stream, delimiter='\t'))
    assert len(diagnostics) == 120
    write_table(output / 'fusion_diagnostics.tsv', diagnostics)
    report = dict(ordinary_searches=6, augmented_searches=12, independent_biological_communities=1,
                  independent_spectral_acquisitions=2, control_seeds=2,
                  independence='evaluated acquisition excluded from ALL selection stages in independent arm',
                  reservoir='common SIHUMIx reference, not the full production reservoir',
                  uncertainty='two acquisitions and nested seeds; no biological significance or FDR certification',
                  calibration_established=False,
                  outputs={n: sha(output / n) for n in ('ordinary_counts.tsv', 'fusion_diagnostics.tsv')},
                  collector_sha256=sha(__file__))
    (output / 'VALIDATION.json').write_text(json.dumps(report, indent=2))
    lines = ['# CAMPI common-reference developmental comparison', '',
             'All 18 searches completed: six ordinary and six per fusion-control seed.', '',
             'The fixed reference, self-selected and independently selected databases are searched against the same evaluated acquisition. Each independent database was rebuilt from the other acquisition at cohort selection, acquisition selection and razor. S01 and S02 are separate technical acquisitions of one community; they are not independent biological communities.', '',
             '| Acquisition | Arm | Target peptides (9–50) | Recovered fixed identities | Additional identities |',
             '|---|---|---:|---:|---:|']
    for row in rows:
        if row['min_residues'] == 9:
            lines.append(f"| {row['evaluated']} | {row['arm']} | {row['canonical_target_peptides']} | {row['recovered_fixed']} | {row['additional_vs_fixed']} |")
    lines += ['', 'These are yields and identity overlaps at peptide q ≤ 0.01. Additional identities are not automatically true discoveries; missing reference-search identities are not automatically errors.', '',
              '| Acquisition | Arm | Seed | Paired diagnostic, 9–50 (%) |', '|---|---|---|---:|']
    for row in diagnostics:
        if float(row['q_threshold']) == .01 and row['min_residues'] == '9':
            value = 100 * float(row['paired_upper_diagnostic_with_unpaired_worst_case'])
            lines.append(f"| {row['evaluated']} | {row['arm']} | {row['seed']} | {value:.3f} |")
    lines += ['', 'Fusion diagnostics include the augmented target/control discovery population. They are assumption-dependent and cannot be substituted for the ordinary-target true error rate. Both full 5–50 and 9–50 populations and all five thresholds are retained in the source table. A value above 1% is not by itself evidence that true error exceeds 1%.', '',
              'This common-reference experiment isolates selection using a manageable, known candidate source. It does not establish calibration for the 233-million-sequence production reservoir, larger cohort evidence pools, protein groups or study-wide unions. The separate generating-spectrum experiment is needed to count actual errors.']
    (output / 'REPORT.md').write_text('\n'.join(lines) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
