#!/usr/bin/env python
"""Figure 4d spectrum audit, verbatim logic of R8A_loss_psm5.sbatch (canonical/final_revision/
thread_selection_math), applied to two rank-5 Sage searches of the same acquisition.

usage: audit.py FL_results.sage.tsv[.gz] MG_results.sage.tsv[.gz] MG.fasta SAMPLE OUT.json

Definitions (unchanged from R8A):
  best target per spectrum   = highest-hyperscore TARGET row (label == 1) among the reported ranks
  FL-only spectrum           = FL best target has peptide_q <= 0.01 and the MG best target is absent
                               or has peptide_q > 0.01
  peptide                    = residue letters only (modifications stripped), exact I/L
  A_absent        FL peptide not in the tryptic digest of the MG FASTA and not a substring of any MG protein
  B_not_cleavable FL peptide is a substring of an MG protein but not a tryptic peptide of it
  C_top_fails_q   MG best target is the same peptide but fails q
  D_displaced     MG best target is a different peptide (and the FL peptide is in the MG digest)
  E_no_target     MG reported no target row for the spectrum although the peptide is in its digest
"""
import csv, gzip, hashlib, json, os, sys, time
sys.path.insert(0, os.environ.get("FASTALAKE_R3_DIGEST_DIR", os.path.dirname(os.path.abspath(__file__))))
import R3_digest as DG
csv.field_size_limit(1 << 30)

def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(8 << 20), b""):
            h.update(b)
    return h.hexdigest()

def opn(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)

def best_target(path):
    """Per spectrum: the best-scoring TARGET row and whether it passes q <= 0.01."""
    out = {}; rows = 0; maxrank = 0
    with opn(path) as fh:
        for x in csv.DictReader(fh, delimiter="\t"):
            rows += 1
            try:
                maxrank = max(maxrank, int(x["rank"]))
                if int(x["label"]) != 1: continue
                k = (x["filename"], x["scannr"]); sc = float(x["hyperscore"])
                q = float(x["peptide_q"])
            except (ValueError, KeyError): continue
            p = "".join(c for c in x["peptide"] if c.isalpha())
            if k not in out or sc > out[k][1]: out[k] = (p, sc, q)
    return out, rows, maxrank

def main(flp, mgp, fa, sample, outp):
    t0 = time.time()
    dg = set(); concat = []; nprot = 0
    for a, s in DG.read_fasta(fa):
        nprot += 1; dg.update(DG.digest(s)); concat.append(s)
    concat = "\n".join(concat)
    fl, fl_rows, fl_maxrank = best_target(flp)
    mg, mg_rows, mg_maxrank = best_target(mgp)
    c = {"A_absent": 0, "B_not_cleavable": 0, "C_top_fails_q": 0, "D_displaced": 0, "E_no_target": 0}
    n_fl_q01 = n_shared = 0
    for k, v in fl.items():
        if v[2] > 0.01: continue
        n_fl_q01 += 1
        o = mg.get(k)
        if o is not None and o[2] <= 0.01:
            n_shared += 1; continue
        pep = v[0]
        if pep not in dg:
            c["A_absent" if pep not in concat else "B_not_cleavable"] += 1
        elif o is None:
            c["E_no_target"] += 1
        elif o[0] == pep:
            c["C_top_fails_q"] += 1
        else:
            c["D_displaced"] += 1
    n = sum(c.values())
    assert n == n_fl_q01 - n_shared, (n, n_fl_q01, n_shared)
    n_mg_q01 = sum(1 for v in mg.values() if v[2] <= 0.01)
    res = dict(sample=sample, n_fl_only=n, counts=c,
               fractions={k: (v / float(n) if n else None) for k, v in c.items()},
               n_fl_spectra_q01=n_fl_q01, n_mg_spectra_q01=n_mg_q01, n_shared_spectra_q01=n_shared,
               n_mg_only=n_mg_q01 - n_shared,
               fl_rows=fl_rows, fl_max_rank=fl_maxrank, mg_rows=mg_rows, mg_max_rank=mg_maxrank,
               comparator_fasta=dict(path=fa, proteins=nprot, digest_peptides=len(dg), bytes=os.path.getsize(fa), sha256=sha(fa)),
               inputs=dict(FL_results=dict(path=flp, sha256=sha(flp)), MG_results=dict(path=mgp, sha256=sha(mgp))),
               seconds=round(time.time() - t0, 1))
    json.dump(res, open(outp, "w"), indent=1, sort_keys=True)
    print("  %s n=%d  A %.1f%% B %.1f%% C %.1f%% D %.1f%% E %.1f%%  (FL q01 %d, MG q01 %d, shared %d)" % (
        sample, n, *[100 * c[k] / max(1, n) for k in ("A_absent", "B_not_cleavable", "C_top_fails_q", "D_displaced", "E_no_target")],
        n_fl_q01, n_mg_q01, n_shared))

if __name__ == "__main__":
    main(*sys.argv[1:6])
