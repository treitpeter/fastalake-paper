#!/usr/bin/env python3
"""Frozen CAMPI stage and reference-member audit; never perform a search."""
import csv
import hashlib
import json
import math
import os
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import ahocorasick

RUN = Path(sys.argv[1]).resolve()
sys.path.insert(0, str(RUN / 'frozen'))
from lib import clean, low, load_queries, meta, sha, verify, dump, table, COL


def read_table(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter='\t'))


def fasta(path, audit):
    h, header, parts, n = hashlib.sha256(), None, [], 0
    with open(path, 'rb') as f:
        for raw in f:
            h.update(raw)
            line = raw.strip()
            if line.startswith(b'>'):
                if header is not None:
                    n += 1
                    yield header, b''.join(parts).decode('ascii').upper().replace('I', 'L')
                header, parts = line[1:].decode('ascii'), []
            elif line:
                assert header is not None
                parts.append(line)
    if header is not None:
        n += 1
        yield header, b''.join(parts).decode('ascii').upper().replace('I', 'L')
    audit.update(sha256=h.hexdigest(), entries=n)


def auto(peptides):
    a = ahocorasick.Automaton()
    for p in sorted(peptides):
        a.add_word(p, p)
    a.make_automaton()
    return a


def compatible(seq, start, end):
    # Exact KR/P production digestion: source enzyme.rs has no initiator-M clipping.
    n = start == 0 or (seq[start - 1] in 'KR' and seq[start] != 'P')
    c = end == len(seq) or (seq[end - 1] in 'KR' and seq[end] != 'P')
    missed = sum(seq[i] in 'KR' and seq[i + 1] != 'P' for i in range(start, end - 1))
    return n and c and missed <= 2


def scan(path, peptides):
    machine = auto(peptides)
    present, digest, audit = Counter(), Counter(), {}
    for _, seq in fasta(path, audit):
        found, eligible = set(), set()
        for end, p in machine.iter(seq):
            found.add(p)
            if compatible(seq, end + 1 - len(p), end + 1):
                eligible.add(p)
        present.update(found)
        digest.update(eligible)
    return present, digest, audit


def fixture():
    assert compatible('KACDEKAA', 1, 6)
    assert not compatible('MACDEKAA', 1, 6)
    assert not compatible('KPACDEKAA', 1, 7)
    assert compatible('PACDEK', 0, 6)
    assert not compatible('ACDEKP', 0, 5)
    assert compatible('ACDEK', 0, 5)
    assert not compatible('AKAKAKAK', 0, 8)


def run(index):
    fixture()
    m = json.loads((RUN / 'inputs/manifest.json').read_text())
    for x in m['frozen_code'] + [m['protocol'], m['members'], m['historical_flags']]:
        verify(x)
    t = m['tasks'][index]
    out = RUN / 'results' / t['sample']
    out.mkdir(parents=True, exist_ok=True)
    assert not (out / 'COMPLETE.json').exists(), 'Completed result already exists'
    for x in t['inputs_to_verify']:
        verify(x)
    sm = t['selection_inputs']
    for k in ('predictions', 'candidates', 'final_fasta'):
        p = Path(sm[k]['path']); st = p.stat()
        assert (st.st_size, st.st_mtime_ns) == (sm[k]['bytes'], sm[k]['mtime_ns'])
    sets = {arm: set(Path(x['path']).read_text().splitlines()) for arm, x in t['peptide_sets'].items()}
    union = set().union(*sets.values())
    for p in union:
        assert 5 <= len(p) <= 50 and 'I' not in p
    member_map = read_table(m['members']['path'])
    members = [r['member'] for r in member_map]
    members_for_peptide = defaultdict(set)
    protein_counts, ref_audit = Counter(), {}
    machine = auto(union)
    for header, seq in fasta(m['reference_fasta']['path'], ref_audit):
        labels = [r['member'] for r in member_map if ' OS=' + r['campi_os_string'] + ' ' in header + ' ']
        assert len(labels) <= 1, header
        if not labels:
            protein_counts['not_one_of_eight'] += 1
            continue
        member = labels[0]; protein_counts[member] += 1
        for _, p in machine.iter(seq):
            members_for_peptide[p].add(member)
    assert ref_audit['sha256'] == m['reference_fasta']['sha256']
    for r in member_map:
        assert protein_counts[r['member']] == int(r['reference_proteins']), (r, protein_counts)
    old = {r['peptide']: r for r in read_table(m['historical_flags']['path'])}
    old_compared, old_absent = 0, []
    for p in union:
        if p not in old:
            old_absent.append(p); continue
        actual = members_for_peptide[p]
        assert len(actual) == int(old[p]['n_members']), (p, actual, old[p])
        assert actual == set(filter(None, old[p]['members'].split(';'))), (p, actual, old[p])
        old_compared += 1
    table(out / 'member_counts.tsv', [dict(sample=t['sample'], arm=arm, min_length=lo,
          member=member, specific_peptides=sum(len(members_for_peptide[p]) == 1 and member in members_for_peptide[p]
          for p in peps if len(p) >= lo)) for arm, peps in sets.items() for lo in (5, 9) for member in members])
    table(out / 'peptide_membership.tsv', [dict(peptide=p, n_members=len(members_for_peptide[p]),
          members=';'.join(sorted(members_for_peptide[p])), arms=';'.join(a for a in sets if p in sets[a])) for p in sorted(union)])
    kept, filtering, _ = load_queries(sm['predictions']['path'])
    assert filtering == t['selection_check']['filter']
    all_predictions, qualifying, beam_ranks = set(), set(), Counter()
    with open(sm['predictions']['path']) as f:
        for row in csv.DictReader(f):
            p = clean(row.get(COL) or '').replace('I', 'L')
            if p: all_predictions.add(p)
            beam_ranks[row.get('beam_idx', row.get('beam_index', row.get('beam_rank', 'FIELD_NOT_FOUND')))] += 1
            try: score = float(row.get('score') or '')
            except ValueError: continue
            if 9 <= len(p) <= 50 and score > 0 and not low(clean(row.get(COL) or '')):
                qualifying.add(p)
    refs = sets['REF']
    candidate_present, candidate_digest, candidate_audit = scan(sm['candidates']['path'], refs)
    final_present, final_digest, final_audit = scan(sm['final_fasta']['path'], refs)
    assert candidate_audit == t['selection_check']['candidate_fasta']
    assert final_audit == t['selection_check']['final_fasta']
    previous = {r['reference_peptide']: r for r in read_table(t['reference_retention']['path'])}
    assert refs <= set(previous)
    rows, stages = [], Counter()
    for p in sorted(refs):
        assert final_present[p] == int(previous[p]['final_proteins']), p
        assert (p in sets['FL']) == bool(int(previous[p]['identified_fl_anywhere'])), p
        assert final_present[p] <= candidate_present[p], p
        assert final_digest[p] <= candidate_digest[p], p
        if p in sets['FL']:
            stage = 'recovered_FL'
            assert final_digest[p] > 0, ('Accepted FL peptide incompatible with audited digest', p)
        elif final_digest[p]: stage = 'final_digest_compatible_not_accepted'
        elif final_present[p]: stage = 'final_sequence_present_digest_incompatible'
        elif candidate_digest[p]: stage = 'candidate_digest_compatible_removed_by_selection'
        elif candidate_present[p]: stage = 'candidate_sequence_present_digest_incompatible_absent_final'
        else: stage = 'absent_from_acquisition_candidates'
        stages[stage] += 1
        rows.append(dict(sample=t['sample'], peptide=p, length=len(p), stage=stage,
            direct_prediction_any=int(p in all_predictions), direct_prediction_qualifying=int(p in qualifying),
            direct_prediction_lookup=int(p in kept), candidate_parents=candidate_present[p],
            candidate_digest_parents=candidate_digest[p], final_parents=final_present[p],
            final_digest_parents=final_digest[p], accepted_FL=int(p in sets['FL']),
            accepted_FL_R=int(p in sets['FL_R']), accepted_FL_R_M=int(p in sets['FL_R_M']),
            n_members=len(members_for_peptide[p]), members=';'.join(sorted(members_for_peptide[p]))))
    assert sum(stages.values()) == len(refs)
    table(out / 'reference_gap.tsv', rows)
    table(out / 'stage_counts.tsv', [dict(sample=t['sample'], min_length=lo, stage=stage,
        peptides=sum(r['stage'] == stage and r['length'] >= lo for r in rows),
        no_direct_prediction=sum(r['stage'] == stage and r['length'] >= lo and not r['direct_prediction_any'] for r in rows))
        for lo in (5, 9) for stage in sorted(stages)])
    dump(out / 'COMPLETE.json', dict(status='PASS', sample=t['sample'], job_id=os.environ.get('SLURM_JOB_ID'),
        reference_peptides=len(refs), counts={a: len(s) for a, s in sets.items()}, stages=dict(stages),
        source_manifest_sha256=sha(RUN / 'inputs/manifest.json'),
        candidate_fasta=candidate_audit, final_fasta=final_audit, reference_fasta=ref_audit,
        historical_member_labels_agree=old_compared, peptides_absent_historical_map=len(old_absent),
        protein_counts=dict(protein_counts), prediction_filter=filtering, observed_beam_field_counts=dict(beam_ranks),
        outputs=[meta(p) for p in out.iterdir() if p.is_file() and p.name != 'COMPLETE.json']))
    print('PASS', t['sample'], dict(stages), flush=True)


if __name__ == '__main__':
    run(int(sys.argv[2]))
