"""Independent accounting audit of frozen fusion searches; no FDR certification.

Uses only the Python standard library. Does not import the original estimator.
Q-based output is reproduced exactly. Score-based ordering is a sensitivity
analysis of the same accepted set, not a newly validated acceptance procedure.
"""
import argparse
import collections
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import subprocess
import time

ROOT = Path(__file__).resolve().parent
THRESHOLDS = (0.001, 0.005, 0.01, 0.02, 0.05)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    tmp.replace(path)


def write_tsv(path, rows):
    rows = list(rows)
    if not rows:
        raise ValueError('Refusing a schema-free empty table')
    with Path(path).open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
        w.writeheader()
        w.writerows(rows)


def source_commit():
    dirty = subprocess.check_output(
        ['git', '-C', str(ROOT), 'status', '--porcelain', '--untracked-files=no'], text=True)
    if dirty.strip():
        raise RuntimeError('Commit source before running: ' + dirty)
    return subprocess.check_output(['git', '-C', str(ROOT), 'rev-parse', 'HEAD'], text=True).strip()


def canonical(peptide):
    value = re.sub(r'\[[^\]]*\]', '', peptide).replace('I', 'L')
    if not re.fullmatch('[A-Z]+', value):
        raise ValueError('Unrecognised peptide notation: ' + peptide)
    return value


def fasta(path):
    header, chunks = None, []
    with Path(path).open() as f:
        for line in f:
            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(chunks)
                header, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.strip())
    if header is not None:
        yield header, ''.join(chunks)


def intervals(sequence):
    # Independently enumerate full tryptic segments and their joins, including
    # protein-terminal fragments and initiator methionine clipping.
    cuts = [0]
    for index in range(len(sequence)):
        if sequence[index] in 'KR' and (
                index == len(sequence) - 1 or sequence[index + 1] != 'P'):
            cuts.append(index + 1)
    if cuts[-1] != len(sequence):
        cuts.append(len(sequence))
    for left in range(len(cuts) - 1):
        for n_segments in (1, 2, 3):
            right = left + n_segments
            if right >= len(cuts):
                continue
            a, b = cuts[left], cuts[right]
            if 5 <= b - a <= 50:
                yield a, b
            if a == 0 and sequence.startswith('M') and 5 <= b - 1 <= 50:
                yield 1, b


def reconstruct_pairs(path):
    forward, reverse = collections.defaultdict(set), collections.defaultdict(set)
    records = iter(fasta(path))
    count = 0
    for t_acc, t_seq in records:
        try:
            e_acc, e_seq = next(records)
        except StopIteration as error:
            raise ValueError('Missing protein partner') from error
        if not t_acc.startswith('T_') or e_acc != 'E_' + t_acc[2:]:
            raise ValueError('Search FASTA is not paired T_/E_ records')
        t_seq, e_seq = t_seq.replace('I', 'L'), e_seq.replace('I', 'L')
        if len(t_seq) != len(e_seq) or list(intervals(t_seq)) != list(intervals(e_seq)):
            raise ValueError('Control digestion differs from target')
        count += 1
        for a, b in intervals(t_seq):
            t, e = t_seq[a:b], e_seq[a:b]
            if not re.fullmatch('[ACDEFGHKLMNPQRSTVWY]+', t + e):
                continue
            if sorted(t) != sorted(e):
                raise ValueError('Control peptide composition differs')
            forward[t].add(e)
            reverse[e].add(t)
    if count == 0:
        raise ValueError('Empty search FASTA')
    shared = forward.keys() & reverse.keys()
    pairs = {}
    for t, ee in forward.items():
        if t in shared or len(ee) != 1:
            continue
        e = next(iter(ee))
        if e not in shared and len(reverse[e]) == 1:
            pairs[t] = e
    reasons = {}
    paired = pairs.keys() | set(pairs.values())
    for p in (forward.keys() | reverse.keys()) - paired:
        if p in shared:
            reasons[p] = 'shared_target_control_sequence'
        elif p in forward:
            reasons[p] = ('partner_shared' if any(e in shared for e in forward[p])
                          else 'nonbijective_mapping')
        else:
            reasons[p] = ('partner_shared' if any(t in shared for t in reverse[p])
                          else 'nonbijective_mapping')
    return pairs, reasons, {'protein_pairs': count, 'valid_pairs': len(pairs),
                            'shared_digest_sequences': len(shared)}


def read_search(path):
    peptides = {}
    n_rows = 0
    with Path(path).open() as f:
        for row in csv.DictReader(f, delimiter='\t'):
            n_rows += 1
            if row['label'] != '1':
                continue
            p = canonical(row['peptide'])
            q, score = float(row['peptide_q']), float(row['sage_discriminant_score'])
            if not math.isfinite(q) or not 0 <= q <= 1 or not math.isfinite(score):
                raise ValueError('Invalid q-value or score')
            record = peptides.setdefault(p, {'q': math.inf, 'score': -math.inf,
                                             'classes': set(), 'forms': set()})
            record['q'] = min(record['q'], q)
            record['score'] = max(record['score'], score)
            record['forms'].add(row['peptide'])
            for accession in re.split('[;,]', row['proteins']):
                if accession.startswith('T_'):
                    record['classes'].add('target')
                elif accession.startswith('E_'):
                    record['classes'].add('entrapment')
                elif not accession.startswith('rev_'):
                    raise ValueError('Unknown protein class: ' + accession)
            if not record['classes']:
                raise ValueError('Target row has no evaluated protein class')
    if n_rows == 0:
        raise ValueError('Empty search output')
    return peptides, n_rows


def decompose(peptides, pairs, reasons, q, minimum, maximum):
    accepted = {p for p, r in peptides.items()
                if r['q'] <= q and minimum <= len(p) <= maximum}
    paired = pairs.keys() | set(pairs.values())
    unknown = accepted - paired
    counters = collections.Counter()
    for p in unknown:
        counters[reasons.get(p, 'outside_audited_digest')] += 1
    ne = alone = q_better = q_tied = score_better = score_tied = discordant = 0
    for t, e in pairs.items():
        if e not in accepted:
            continue
        ne += 1
        if t not in accepted:
            alone += 1
            continue
        eq, tq = peptides[e]['q'], peptides[t]['q']
        es, ts = peptides[e]['score'], peptides[t]['score']
        q_better += eq < tq
        q_tied += eq == tq
        score_better += es > ts
        score_tied += es == ts
        discordant += (eq < tq and es < ts) or (eq > tq and es > ts)
    n = len(accepted)
    numerator = ne + alone + 2 * (q_better + q_tied)
    score_numerator = ne + alone + 2 * (score_better + score_tied)
    row = dict(q_threshold=q, min_residues=minimum, max_residues=maximum,
               n_canonical=n,
               n_target_only=sum(peptides[p]['classes'] == {'target'} for p in accepted),
               n_entrapment_only=sum(peptides[p]['classes'] == {'entrapment'} for p in accepted),
               n_shared=sum(len(peptides[p]['classes']) > 1 for p in accepted),
               n_unpaired_accepted=len(unknown), n_paired_entrapment=ne,
               n_entrapment_alone=alone, n_entrapment_q_strict_better=q_better,
               n_entrapment_q_tied=q_tied, n_entrapment_score_strict_better=score_better,
               n_entrapment_score_tied=score_tied, n_q_score_order_discordant=discordant,
               paired_numerator=numerator, score_order_numerator_sensitivity=score_numerator,
               paired_diagnostic=(numerator + len(unknown)) / n if n else None,
               score_order_diagnostic_sensitivity=(score_numerator + len(unknown)) / n if n else None,
               paired_component_fraction=numerator / n if n else None,
               unpaired_component_fraction=len(unknown) / n if n else None,
               q_tie_component_fraction=2 * q_tied / n if n else None,
               n_multiple_peptidoforms=sum(len(peptides[p]['forms']) > 1 for p in accepted))
    for reason in ('shared_target_control_sequence', 'partner_shared',
                   'nonbijective_mapping', 'outside_audited_digest'):
        row['n_unknown_' + reason] = counters[reason]
    if sum(row['n_' + key] for key in ('target_only', 'entrapment_only', 'shared')) != n:
        raise ValueError('Membership accounting does not partition discoveries')
    return row


def run_task(task, output):
    start = time.time()
    commit = source_commit()
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    base = Path(task['directory'])
    result = base / '5_sage/results.sage.tsv'
    digest = sha(result)
    if digest != task['expected_results_sha256']:
        raise ValueError('Raw search hash differs from frozen completion manifest')
    pairs, reasons, pairing = reconstruct_pairs(base / 'search_input/search.fasta')
    with (base / 'search_input/peptide_pairs.tsv').open() as f:
        saved_rows = list(csv.DictReader(f, delimiter='\t'))
    saved_pairs = {r['target']: r['entrapment'] for r in saved_rows}
    if len(saved_pairs) != len(saved_rows) or pairs != saved_pairs:
        raise ValueError('Independent pairing does not match original pairing')
    peptides, n_rows = read_search(result)
    rows = [decompose(peptides, pairs, reasons, q, low, high)
            for q in THRESHOLDS for low, high in ((5, 50), (9, 50), (5, 8))]
    completion = json.loads((base / 'COMPLETE.json').read_text())
    checks = ('n_canonical', 'n_target_only', 'n_entrapment_only', 'n_shared',
              'n_unpaired_accepted', 'n_paired_entrapment', 'paired_numerator')
    for original in completion['summary']:
        row = next(r for r in rows if r['q_threshold'] == original['q_threshold']
                   and r['min_residues'] == original['min_residues'] and r['max_residues'] == 50)
        for key in checks:
            if row[key] != original[key]:
                raise ValueError(f'Independent recount differs: {key}: {row[key]} != {original[key]}')
    for row in rows:
        row.update(cohort=task['cohort'], sample=task['sample'], seed=task['seed'])
    write_tsv(output / 'decomposition.tsv', rows)
    paired_members = pairs.keys() | set(pairs.values())
    accepted_rows = [dict(peptide=p, peptide_q=r['q'], score=r['score'],
                          assignment=';'.join(sorted(r['classes'])),
                          pairing_status='paired' if p in paired_members
                          else reasons.get(p, 'outside_audited_digest'))
                     for p, r in peptides.items() if r['q'] <= 0.05]
    if accepted_rows:
        write_tsv(output / 'accepted_sequences.tsv', accepted_rows)
    else:
        (output / 'accepted_sequences.tsv').write_text(
            'peptide\tpeptide_q\tscore\tassignment\tpairing_status\n')
    write_json(output / 'COMPLETE.json', dict(
        task=task, source_commit=commit, seconds=time.time() - start, input_rows=n_rows,
        independent_accounting_matches=True, independent_final_pairing_matches=True,
        scientific_calibration='not_established_by_accounting', pairing=pairing,
        inputs={str(p): sha(p) for p in (result, base / 'COMPLETE.json',
                                       base / 'search_input/search.fasta',
                                       base / 'search_input/peptide_pairs.tsv')},
        outputs={p.name: sha(p) for p in output.glob('*.tsv')}))
    print(json.dumps({'task': task['index'], 'seconds': time.time() - start,
                      'accounting': 'PASS'}), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--manifest', type=Path, default=ROOT / 'inputs/fusion_tasks.json')
    parser.add_argument('--index', required=True, type=int)
    parser.add_argument('--output-root', type=Path, default=ROOT / 'results/audit')
    args = parser.parse_args()
    tasks = json.loads(args.manifest.read_text())
    run_task(tasks[args.index], args.output_root / str(args.index))
