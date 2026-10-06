"""Meaningful checks for sequence boundaries, context and classification."""
from pathlib import Path
import tempfile

from analyze import canonical, category, locate, spectrum_states

assert canonical("[+42.0106]-AC[+57.0215]ILK") == "ACLLK"
records = [("p1", "IACD"), ("p2", "EFGL"), ("p3", "ACDACD"), ("p4", "ACDEFG")]
queries = {"ACD", "DEFG", "ACDEFG", "LACD", "CDE", "DEFGL"}
with tempfile.TemporaryDirectory(prefix="bag_gainloss_test_") as folder:
    path = Path(folder) / "input.fasta"
    path.write_text("".join(f">{header}\n{sequence[:2]}\n{sequence[2:]}\n" for header, sequence in records))
    proteins, contexts, memberships = locate(path, queries)
    assert len(proteins) == len(records)
    for peptide in queries:
        expected_contexts, expected_members = set(), set()
        for header, original in records:
            sequence = original.replace("I", "L")
            for offset in range(len(sequence) - len(peptide) + 1):
                if sequence[offset:offset + len(peptide)] == peptide:
                    end = offset + len(peptide)
                    expected_contexts.add((sequence[offset-1] if offset else "^", sequence[end] if end < len(sequence) else "$"))
                    expected_members.add(header)
        assert contexts[peptide] == expected_contexts
        assert memberships[peptide] == expected_members
assert category("gained", False, None) == "sequence_added"
assert category("lost", False, None) == "sequence_removed"
assert category("lost", True, 0.02, True) == "sequence_retained_above_q_cutoff"
assert category("gained", True, None) == "sequence_retained_no_reported_target"
assert category("gained", True, None, True) == "sequence_retained_nonfinite_q"
for args in (("gained", False, 0.02, True), ("lost", True, 0.005, True)):
    try:
        category(*args)
    except AssertionError:
        pass
    else:
        raise AssertionError("Inconsistent category input accepted")
states = spectrum_states(range(5), "ACD", {0: ("ACD", 1, 0.02), 1: ("EFG", 1, 0.001),
                                          2: ("EFG", 1, 0.04), 3: ("EFG", -1, 0.01)})
assert set(states.values()) == {1} and sum(states.values()) == 5
print("PASS: canonicalization, exhaustive sequence/context matching, record boundaries, categories, inconsistent-input rejection and same-spectrum states")
