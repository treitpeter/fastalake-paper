"""Historical de novo reader: for one acquisition, counts exported predictions with score >= 0.9 (HICONF) and the number whose C-terminal residue is not K or R, for all lengths. Reproduces the original manuscript counts that A26 must match before any new route is trusted."""
import csv
SEQ_COL='peptide_prediction_detokenized_unmodified'
SCORE_COL='score'
HICONF=0.9
TRYPTIC={'K','R'}
def denovo_one(args):
    sample, path = args
    n_hi = non_hi = n_all = non_all = 0
    with open(path, newline="") as fh:
        r = csv.DictReader(fh)
        if SEQ_COL not in (r.fieldnames or []):
            raise KeyError(f"{path}: no {SEQ_COL!r}")
        for row in r:
            seq = "".join(ch for ch in (row.get(SEQ_COL) or "")
                          if ch.isascii() and ch.isupper())
            if not seq:
                continue
            try:
                sc = float(row.get(SCORE_COL) or "nan")
            except ValueError:
                continue
            if not (sc > 0.0):
                continue
            nk = seq[-1] not in TRYPTIC
            n_all += 1
            non_all += nk
            if sc >= HICONF:
                n_hi += 1
                non_hi += nk
    return dict(sample=sample, dn_n=n_hi, dn_non=non_hi,
                dn_rate=(non_hi / n_hi) if n_hi else None,
                dn_n_ungated=n_all,
                dn_rate_ungated=(non_all / n_all) if n_all else None)
