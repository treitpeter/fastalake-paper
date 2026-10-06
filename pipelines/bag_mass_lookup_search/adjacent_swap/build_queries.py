"""C2 query sets. Adjacent-swap definition is A04's (frozen/benchmark.py line 42):
   {q[:i]+q[i+1]+q[i]+q[i+2:] for i in range(len(q)-1) if q[i]!=q[i+1]}
Inputs: the A02 retained, I/L-collapsed distinct query sets (authentic REAL; scrambled SIM seed 20260913).
Writes one extractor input directory per set (queries.csv with generic `sequence,score` columns, score 1.0),
a variant map for attribution, and QUERIES.json with every count."""
import os
import csv, json, os, random, sys, hashlib
R = os.environ["FASTALAKE_OUT_DIR"]
OUT = R + "/C2/queries"
PILOT_N, PILOT_SEED = 1000, 20260917

def swaps(q):
    """Every single adjacent transposition of q, skipping equal neighbours (A04 definition)."""
    return {q[:i] + q[i + 1] + q[i] + q[i + 2:] for i in range(len(q) - 1) if q[i] != q[i + 1]}

def load(p):
    """Read a frozen query table and return its I/L-collapsed distinct queries."""
    qs = [r["query"].replace("I", "L") for r in csv.DictReader(open(p), delimiter="\t")]
    assert len(qs) == len(set(qs)), "queries not distinct after I/L collapse"
    assert all(9 <= len(q) <= 50 for q in qs)
    return sorted(qs)

def write_set(name, originals, swap_only, vmap_rows):
    """Write one query set as an extractor input CSV."""
    d = os.path.join(OUT, name); os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, "queries.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["sequence", "score"])
        for q in originals: w.writerow([q, "1.0"])
        for v in swap_only: w.writerow([v, "1.0"])
    with open(os.path.join(d, "originals.txt"), "w") as f:
        f.write("\n".join(originals) + "\n")
    with open(os.path.join(d, "swap_only.txt"), "w") as f:
        f.write("\n".join(swap_only) + ("\n" if swap_only else ""))
    with open(os.path.join(d, "variant_map.tsv"), "w", newline="") as f:
        w = csv.writer(f, delimiter="\t"); w.writerow(["variant", "origin_query", "length", "equals_an_original"])
        w.writerows(vmap_rows)
    return {"dir": d, "n_originals": len(originals), "n_swap_only": len(swap_only), "n_rows": len(originals) + len(swap_only)}

def build(arm, qs, name_orig, name_all):
    """Build the originals and originals-plus-swaps query sets of one arm."""
    orig = set(qs)
    var = {}; vmap = []
    for q in qs:
        for v in sorted(swaps(q)):
            var.setdefault(v, []).append(q)
    swap_only = sorted(v for v in var if v not in orig)
    n_eq_orig = sum(1 for v in var if v in orig)
    for v in sorted(var):
        for q in var[v]:
            vmap.append([v, q, len(v), int(v in orig)])
    info = {"arm": arm, "n_queries": len(qs), "sum_L_minus_1": sum(len(q) - 1 for q in qs),
            "distinct_swap_variants": len(var), "swap_variants_equal_to_an_original": n_eq_orig,
            "swap_only_variants": len(swap_only),
            "mean_swap_only_per_query": round(len(swap_only) / len(qs), 3)}
    info["set_orig"] = write_set(name_orig, qs, [], [])
    info["set_all"] = write_set(name_all, qs, swap_only, vmap)
    return info

REAL = load(R + "/inputs/REAL_queries.tsv")
SIM = load(R + "/inputs/SIM_20260913_queries.tsv")
rng = random.Random(PILOT_SEED); P = sorted(rng.sample(REAL, PILOT_N))
rng = random.Random(PILOT_SEED); PS = sorted(rng.sample(SIM, PILOT_N))
J = {"pilot_seed": PILOT_SEED, "pilot_n": PILOT_N, "swap_definition": "A04 frozen/benchmark.py line 42: single adjacent-residue transposition, identical variants skipped",
     "pilot": {"REAL": build("REAL_pilot", P, "P_REAL_orig", "P_REAL_all"), "SIM": build("SIM_pilot", PS, "P_SIM_orig", "P_SIM_all")},
     "full": {"REAL": build("REAL", REAL, "REAL_orig", "REAL_all"), "SIM": build("SIM", SIM, "SIM_orig", "SIM_all")}}
for k, d in [("pilot_REAL", P), ("pilot_SIM", PS)]:
    J[k + "_sha256_of_list"] = hashlib.sha256("\n".join(d).encode()).hexdigest()
json.dump(J, open(OUT + "/QUERIES.json", "w"), indent=2)
print(json.dumps({k: {a: {x: y for x, y in v.items() if not x.startswith("set_")} for a, v in J[k].items()} for k in ("pilot", "full")}, indent=1))
