"""C2 summary for one arm: run A (originals) vs run B (originals + swap-only variants).
usage: summarize.py <set_dir_all> <runA_dir> <runB_dir> <parsimony_fasta_B or -> <parsimony_fasta_A or -> <out_json> [a02_reference_proteins_matched]
(a) distinct proteins by originals, (b) additional proteins by swaps only, (c) swap-only proteins in the razor cover."""
import os
import csv, json, sys, collections, os
set_dir, runA, runB, parsB, parsA, out = sys.argv[1:7]
ref = int(sys.argv[7]) if len(sys.argv) > 7 else None
originals = set(open(os.path.join(set_dir, "originals.txt")).read().split())
swap_only = set(open(os.path.join(set_dir, "swap_only.txt")).read().split())
vmap = collections.defaultdict(set)
for r in csv.DictReader(open(os.path.join(set_dir, "variant_map.tsv")), delimiter="\t"):
    if r["equals_an_original"] == "0":
        vmap[r["variant"]].add(r["origin_query"])

def hits(d):
    """Read an extractor hits.tsv: protein ids matched per query."""
    pep2prot = collections.defaultdict(set)
    with open(os.path.join(d, "hits.tsv")) as f:
        rd = csv.reader(f, delimiter="\t"); next(rd)
        for pep, prot, pos in rd:
            pep2prot[pep].add(prot)
    return pep2prot

A = hits(runA); B = hits(runB)
protA = set().union(*A.values()) if A else set()
protB = set().union(*B.values()) if B else set()
unknown = [p for p in B if p not in originals and p not in swap_only]
assert not unknown, f"{len(unknown)} hit peptides in run B belong to neither set, e.g. {unknown[:3]}"
protB_orig = set().union(*(B[p] for p in B if p in originals)) if any(p in originals for p in B) else set()
protB_swap = set().union(*(B[p] for p in B if p in swap_only)) if any(p in swap_only for p in B) else set()
swap_only_prot = protB - protB_orig
consistency = {"runA_proteins_equal_runB_original_hit_proteins": protA == protB_orig,
               "B_minus_A_equals_swap_only": (protB - protA) == swap_only_prot}
matched_orig_A = {p for p in A}
unmatched_orig = originals - matched_orig_A
rescued = {q for v in B if v in swap_only for q in vmap[v] if q in unmatched_orig}
by_len = {}
for L in sorted({len(v) for v in swap_only} | {len(q) for q in originals}):
    vs = [v for v in swap_only if len(v) == L]
    vh = [v for v in vs if v in B]
    pro = set().union(*(B[v] for v in vh)) if vh else set()
    oq = [q for q in originals if len(q) == L]
    oh = [q for q in oq if q in A]
    by_len[L] = {"originals": len(oq), "originals_with_hit": len(oh), "swap_variants": len(vs), "swap_variants_with_hit": len(vh),
                 "swap_only_proteins_at_this_length": len(pro & swap_only_prot)}

def pars_ids(p):
    """Accessions of a parsimony razor FASTA."""
    if p == "-" or not os.path.exists(p): return None
    ids = set()
    for line in open(p):
        if line.startswith(">"): ids.add(line[1:].split()[0])
    return ids
PB = pars_ids(parsB); PA = pars_ids(parsA)
S = {"set_dir": set_dir, "runA": runA, "runB": runB,
     "n_originals": len(originals), "n_swap_only_variants": len(swap_only),
     "a_proteins_by_originals": len(protA), "b_additional_proteins_by_swaps_only": len(swap_only_prot),
     "b_over_a_fraction": round(len(swap_only_prot) / len(protA), 4) if protA else None,
     "proteins_run_B_total": len(protB),
     "originals_with_exact_hit": len(matched_orig_A), "originals_without_exact_hit": len(unmatched_orig),
     "originals_without_exact_hit_rescued_by_a_swap": len(rescued),
     "swap_variants_with_hit": sum(1 for v in B if v in swap_only),
     "swap_only_proteins_by_variant_length": by_len, "consistency": consistency}
if PB is not None:
    S["c_parsimony_runB"] = {"fasta": parsB, "proteins_in_razor_cover": len(PB), "swap_only_candidates_in_cover": len(PB & swap_only_prot),
                             "original_candidates_in_cover": len(PB & protA), "cover_ids_not_in_candidates": len(PB - protB)}
if PA is not None:
    S["c_parsimony_runA"] = {"fasta": parsA, "proteins_in_razor_cover": len(PA), "cover_ids_not_in_candidates": len(PA - protA)}
if ref is not None:
    S["positive_control_vs_A02_archived_proteins_matched"] = {"a02": ref, "here": len(protA), "equal": ref == len(protA)}
json.dump(S, open(out, "w"), indent=2)
print(json.dumps({k: v for k, v in S.items() if k != "swap_only_proteins_by_variant_length"}, indent=1))
print("length\toriginals\torig_with_hit\tswap_variants\tswap_with_hit\tswap_only_proteins")
for L, d in by_len.items():
    print(f"{L}\t{d['originals']}\t{d['originals_with_hit']}\t{d['swap_variants']}\t{d['swap_variants_with_hit']}\t{d['swap_only_proteins_at_this_length']}")
