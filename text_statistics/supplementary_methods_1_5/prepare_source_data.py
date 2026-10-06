#!/usr/bin/env python3
"""Build ``source_data/supplementary_methods_1_5/`` from the original producer outputs.

This script is the record of how every table in ``source_data/supplementary_methods_1_5/``
was derived. It reads only the original producer outputs under the lab project directory
named by the environment variable ``FASTALAKE_LAB_ROOT`` (written ``LAB_ROOT`` in
``PROVENANCE.md``) and writes only into ``source_data/supplementary_methods_1_5/`` and
``prepare_receipt.json`` beside this script. It is not needed to run the notebook.

Most tables are byte copies or column subsets of producer outputs. Four tables are
intermediate results recomputed here with the original producer functions so that the
statistics can be re-run from small tables:

* ``allocation_capscan_own_stool_distances_347.tsv``: participant-by-participant mean
  Hellinger distances behind the CapScan own-stool permutation test (code of
  ``alloc_capscan347.py``, run 6064925);
* ``allocation_predict_distances.tsv``: the 77 x 77 Euclidean distance matrices behind the
  PREDICT profile-distance correlations (``check_prevalence.py:transform`` as called by
  ``alloc_sensitivity.py:predict_crossomics``, run 6048208);
* ``allocation_hstool_ko_pairs.tsv``: per-pair Pearson correlations behind the healthy-stool
  KO precision summary (``hstool_precision.py:pairs`` as called by
  ``alloc_sensitivity.py:hstool_precision``, run 6048208);
* ``mgnify_seed_ortholog_check_by_genome.tsv``: per-genome counts of the 9,368,708-row
  seed-ortholog verification table.

Each recomputed table is checked against the producer's own summary before it is written.

Run (on Slurm, about 10 GB and 10 minutes)::

    FASTALAKE_LAB_ROOT=<lab project directory> python prepare_source_data.py
"""
import ast
import csv
import gzip
import hashlib
import json
import os
import re
import shutil
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.spatial.distance import cdist, pdist, squareform
from scipy.stats import rankdata

UNIT = Path(__file__).resolve().parent
REPO = UNIT.parents[1]
SD = REPO / "source_data" / UNIT.name
LAB = Path(os.environ["FASTALAKE_LAB_ROOT"]).resolve()

AF = LAB / "review/v15_integration/jess_integration_2026-09-29/analysis_followup"
C = AF / "manuscript_completion_2026-10-01"
B = C / "biology_refit"
F = B / "claude_review_2026-10-03"
J = AF.parent
GW = AF / "github_workflow_reconciliation_2026-09-30/resumed_audit/annotation_repair"
C2 = F / "coauthor_comment_check_v25_2026-10-05/compute/C2_participant15"
ALLOC = F / "iteration1/proposals_methods/alloc_sensitivity"
A3 = F / "final/reruns_347/A3_allocation_sensitivity"

RECEIPT = {"origins": {}, "outputs": {}, "checks": {}}
FLOAT = "%.17g"


def sha(p):
    """Return the sha256 hex digest of a file."""
    h = hashlib.sha256()
    with Path(p).open("rb") as f:
        for b in iter(lambda: f.read(8 << 20), b""):
            h.update(b)
    return h.hexdigest()


def rel(p):
    """Return a path relative to LAB_ROOT, written LAB_ROOT/..."""
    return "LAB_ROOT/" + Path(p).resolve().relative_to(LAB).as_posix()


def origin(p):
    """Record and return an origin file."""
    RECEIPT["origins"][rel(p)] = sha(p)
    return p


def write(df, name):
    """Write a tab-separated table into source_data and record its sha256."""
    out = SD / name
    df.to_csv(out, sep="\t", index=False, float_format=FLOAT, lineterminator="\n")
    RECEIPT["outputs"][name] = sha(out)
    print("wrote", name, len(df), flush=True)


def copy(src, name):
    """Copy a producer file byte for byte."""
    shutil.copyfile(origin(src), SD / name)
    RECEIPT["outputs"][name] = sha(SD / name)
    print("copied", name, flush=True)


def import_functions(path, names, env):
    """Execute only the named top-level functions of a producer script (as the producer did)."""
    tree = ast.parse(Path(path).read_text())
    defs = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name in names]
    assert {n.name for n in defs} == set(names), (path, names)
    exec(compile(ast.Module(body=defs, type_ignores=[]), str(path), "exec"), env)
    return env


# ---------------------------------------------------------------- simple copies and subsets

def parser_audit():
    """Parser audit, 726 prediction files (A7 rerun of M07, jobs 6064887/6064888)."""
    copy(F / "final/reruns_347/A7_other/parser_audit/fixed726/ACQUISITION_SUMMARY.tsv",
         "parser_audit_acquisitions_726.tsv")


def completeness():
    """Sage input completeness (job 6065843); path columns dropped."""
    src = origin(F / "final/completeness_653/out/search_coverage.tsv")
    d = pd.read_csv(src, sep="\t", dtype=str, keep_default_na=False)
    d = d[d.lineage.isin(["primary_public_current", "healthy_stool_58"])]
    d = d.drop(columns=["mzml_cfg", "dir"])
    assert len(d) == 711
    write(d, "search_input_completeness.tsv")


def hstool97():
    """Original 97-acquisition healthy-stool study-level lookup (production chain stage 1)."""
    ev = LAB / "canonical/publication_data/cohorts/HSTOOL/1_evidence"
    st = json.loads(origin(ev / "evidence_lake.stats.json").read_text())
    st["database"] = Path(st["database"]).name
    rows = [dict(key=k, value=json.dumps(v) if isinstance(v, list) else v) for k, v in st.items()]
    write(pd.DataFrame(rows), "hstool97_stage1_stats.tsv")
    files = sorted((ev / "_pooled_csvs").iterdir())
    RECEIPT["origins"][rel(ev / "_pooled_csvs")] = "directory listing"
    pf = pd.DataFrame(dict(prediction_file=[p.name for p in files], bytes=[p.stat().st_size for p in files]))
    assert len(pf) == 97
    write(pf, "hstool97_pooled_prediction_files.tsv")


def origin_rows():
    """Annotated-origin rows and pooled LFQ intensity (inputs of partial_human.py)."""
    sets = {
        "HSTOOL": GW / "annotated_origin_consensus_2026-10-01/full/HSTOOL_noclust",
        "PREDICT": GW / "annotated_origin_consensus_2026-10-01/full/PREDICT_noclust",
        "CAPSCAN331": GW / "annotated_origin_consensus_2026-10-01/full/CAPSCAN_noclust",
        "CAPSCAN347": C2 / "fig5/ra/annotation_repair/annotated_origin_consensus_2026-10-01/full/CAPSCAN_noclust",
    }
    expect = {"HSTOOL": (58, "0.120"), "PREDICT": (306, "1.964"), "CAPSCAN331": (331, "0.604"), "CAPSCAN347": (347, "0.613")}
    for name, d in sets.items():
        tx = json.loads(origin(d / "taxonomy_rows.json").read_text())
        m = sparse.load_npz(origin(d / "taxonomy_evidence_signal.npz")).tocsr()
        r = np.asarray(m.sum(1)).ravel()
        h = np.array([t["supported_lineage"][-1:] == ["s__Homo sapiens"] for t in tx])
        miss = np.array([t["state"] == "annotated_parent_lca" for t in tx])
        n_acq, pct = expect[name]
        assert m.shape[1] == n_acq
        assert "%.3f" % (100 * r[h & miss].sum() / m.sum()) == pct, name
        df = pd.DataFrame(dict(row=range(len(tx)), rank=[t["rank"] for t in tx], state=[t["state"] for t in tx],
                               lineage=[";".join(t["supported_lineage"]) for t in tx], human_species=h,
                               signal_pooled=r, acquisitions=m.shape[1]))
        write(df, f"origin_rows_{name}.tsv")


def mgnify():
    """Per-genome counts of the MGnify seed-ortholog verification (summarize_verification.py)."""
    src = origin(F / "iteration1/fix_eggnog_native/verify/native_row_verification.tsv.gz")
    d = pd.read_csv(src, sep="\t", usecols=["genome", "native_score", "pair_bitscore", "ratio", "status",
                                             "seed_equals_adjacent_gene_seed", "query_length", "seed_length"])
    ratio = d.pair_bitscore.fillna(0) / d.native_score
    assert np.allclose(ratio, d.ratio, rtol=0, atol=1e-12)
    d["fail"] = ratio < 0.5
    d["no_alignment"] = d.status.eq("pair_no_alignment")
    d["short"] = d.query_length < 0.5 * d.seed_length
    d["adj"] = d.seed_equals_adjacent_gene_seed.astype(bool)
    g = d.groupby("genome")
    out = pd.DataFrame(dict(
        rows=g.size(), no_alignment=g.no_alignment.sum(), fail_ratio_below_half=g.fail.sum(),
        fail_query_shorter_than_half_seed=g.apply(lambda x: int((x.fail & x.short).sum())),
        fail_seed_equals_adjacent=g.apply(lambda x: int((x.fail & x.adj).sum())),
        pass_query_shorter_than_half_seed=g.apply(lambda x: int((~x.fail & x.short).sum())),
        pass_seed_equals_adjacent=g.apply(lambda x: int((~x.fail & x.adj).sum())),
    )).reset_index()
    summ = json.loads(origin(F / "iteration1/fix_eggnog_native/verify/VERIFICATION_SUMMARY.json").read_text())
    assert out.rows.sum() == summ["rows"] and len(out) == summ["genomes"]
    assert out.fail_ratio_below_half.sum() == summ["fail_at_0.5"]["n"]
    assert out.no_alignment.sum() == summ["status"]["pair_no_alignment"]
    write(out, "mgnify_seed_ortholog_check_by_genome.tsv")


def annotation_coverage():
    """MGnify gene-sequence-verified share of compatible proteins per study arm."""
    rows = []
    for d in sorted((GW / "evidence_repaired/full").iterdir()):
        c = json.loads(origin(d / "COMPLETE.json").read_text())
        rows.append(dict(arm=d.name, roster="v28d", acquisitions=c["acquisitions"], compatible_parents=c["compatible_parents"],
                         **{"native_" + k: v for k, v in c["native_eggnog_checks"].items()}))
    c = json.loads(origin(C2 / "fig5/ra/annotation_repair/evidence_repaired/full/CAPSCAN_noclust/COMPLETE.json").read_text())
    rows.append(dict(arm="CAPSCAN_noclust", roster="347", acquisitions=c["acquisitions"], compatible_parents=c["compatible_parents"],
                     **{"native_" + k: v for k, v in c["native_eggnog_checks"].items()}))
    write(pd.DataFrame(rows), "annotation_mgnify_coverage_by_arm.tsv")


def campi_supplement():
    """CAMPI reference proteins annotated in the supplementary eggNOG-mapper run (job 6034651)."""
    d = GW / "supplemental_native_extract"
    s = pd.read_csv(origin(d / "native_sequences.tsv.gz"), sep="\t", dtype=str)
    f = pd.read_csv(origin(d / "native_functions.tsv.gz"), sep="\t", dtype=str)
    c = json.loads(origin(d / "COMPLETE.json").read_text())
    s = s[["accession", "sequence_sha256_IL", "sequence_verified", "terminal_stop_removed"]].copy()
    s["has_eggnog_record"] = s.accession.isin(set(f.accession))
    assert len(s) == c["queries"] and s.has_eggnog_record.sum() == c["annotation_records"] == len(f)
    write(s, "campi_supplemental_annotation_queries.tsv")


UNIPROT = re.compile(r"^(sp|tr)\|([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:-\d+)?\|")
CONTIG = re.compile(r"^k\d+_\d+")


def campi_targets():
    """Records of the two supplied CAMPI target FASTAs searched by Sage."""
    rows = []
    for db, fn in [("reference", "campi_ref_targets.faa"), ("MG", "campi_mo_targets.faa")]:
        p = origin(LAB / "canonical/publication_data/arms/campi_db" / fn)
        head, seq, i = None, [], 0

        def flush():
            s = "".join(seq)
            tok = head[1:].split()[0] if head[1:].split() else ""
            kind = "uniprot" if UNIPROT.match(tok) else "assembly_contig" if CONTIG.match(tok) else "other"
            rows.append(dict(database=db, record=i, header_token=tok, header_kind=kind,
                             header_has_trailing_whitespace=head != head.rstrip(), length=len(s),
                             sequence_sha256=hashlib.sha256(s.encode()).hexdigest()))
        with open(p) as fh:
            for line in fh:
                line = line.rstrip("\n")
                if line.startswith(">"):
                    if head is not None:
                        flush(); i += 1
                    head, seq = line, []
                else:
                    seq.append(line.strip())
            flush()
    write(pd.DataFrame(rows), "campi_target_fasta_records.tsv")


# ---------------------------------------------------------------- allocation sensitivity

def alloc_coverage():
    """Per-acquisition allocated intensity under each allocation rule (coverage.tsv of each run)."""
    for name, p in [("HSTOOL", ALLOC / "out/profiles/HSTOOL_noclust/coverage.tsv"),
                    ("PREDICT", ALLOC / "out/profiles/PREDICT_noclust/coverage.tsv"),
                    ("CAPSCAN347", A3 / "run347/profiles/CAPSCAN_noclust/coverage.tsv")]:
        copy(p, f"allocation_coverage_{name}.tsv")
    copy(A3 / "allocated_fraction_medians_347.tsv", "allocation_fraction_medians_producer.tsv")
    copy(A3 / "capscan_own_stool_331_vs_347.tsv", "allocation_capscan_own_stool_producer.tsv")
    copy(ALLOC / "out/predict_crossomics_protein_rows.tsv", "allocation_predict_crossomics_producer.tsv")
    copy(ALLOC / "out/hstool_KO_precision.tsv", "allocation_hstool_ko_precision_producer.tsv")


def frame(root, arm, key):
    """Profile matrix (acquisitions x features) as written by alloc_sensitivity.py:build_arm."""
    a = json.loads(origin(root / "profiles" / arm / (key + "_axes.json")).read_text())
    v = sparse.load_npz(origin(root / "profiles" / arm / (key + ".npz")))
    return pd.DataFrame(v.toarray(), index=a["samples"], columns=a["features"])


def own_stool_distances():
    """Participant x participant mean Hellinger distances (alloc_capscan347.py, run 6064925)."""
    cfg = json.loads(origin(A3 / "fixed347.json").read_text())
    design = pd.read_csv(origin(C2 / "roster/design_347.tsv"), sep="\t").set_index("sample_id")
    root = A3 / "run347"
    ARM = "CAPSCAN_noclust"
    REG = ["duodenum", "jejunum", "ileum", "colon", "stool"]
    pos = frame(root, ARM, "OG_other_unresolved__equal").reindex(design.index)
    roster = design.index[pos.sum(1) > 0]
    obsr = design.loc[roster].groupby("participant").region.nunique()
    parts = sorted(obsr.index[obsr.eq(5)]); n = len(parts)
    assert n == 15 and cfg["mode"] == "fixed347"
    ix = pd.MultiIndex.from_product([parts, REG], names=["participant", "region"])

    def cells(raw):
        d = design.loc[roster]; v = raw.loc[roster]; ok = v.sum(1) > 0
        relv = v.div(v.sum(1).where(ok), axis=0)
        spec = relv.groupby([d.participant, d.region, d.specimen_id]).mean()
        return spec.groupby(level=[0, 1]).mean().loc[ix]

    def dd_of(x, prev):
        keep = (x > 0).mean(0) >= prev; v = x[:, keep]; h = np.sqrt(v / v.sum(1, keepdims=True)); cube = h.reshape(n, 5, -1)
        return np.stack([cdist(cube[i, :4], cube[:, 4]).mean(0) for i in range(n)]), int(keep.sum())

    rows = []
    for v in ["equal", "unique", "complete", "unique_complete"]:
        x = cells(frame(root, ARM, f"OG_other_unresolved__{v}").reindex(design.index)).to_numpy()
        for prev in [0.0, 0.1, 0.2]:
            dd, nf = dd_of(x, prev)
            rows += [dict(variant=v, prevalence=prev, features=nf, intestinal_participant=parts[i],
                          stool_participant=parts[j], mean_hellinger_distance=dd[i, j]) for i in range(n) for j in range(n)]
    cube = cells(frame(root, ARM, "OG_other_unresolved__equal").reindex(design.index)).to_numpy().reshape(n, 5, -1).copy()
    cube[:, 4] = cube[np.roll(np.arange(n), 1), 4]
    dd, nf = dd_of(cube.reshape(5 * n, -1), 0.0)
    rows += [dict(variant="NEGATIVE_CONTROL_stool_rotated", prevalence=0.0, features=nf, intestinal_participant=parts[i],
                  stool_participant=parts[j], mean_hellinger_distance=dd[i, j]) for i in range(n) for j in range(n)]
    write(pd.DataFrame(rows), "allocation_capscan_own_stool_distances_347.tsv")


def predict_distances():
    """77 x 77 distance matrices of the PREDICT functional profiles (alloc_sensitivity.py:predict_crossomics)."""
    X = J / "source_data/current/Crossomics"
    OLDE = LAB / "analyses/enterotypes_2026-09-07"
    env = dict(np=np, pd=pd, pdist=pdist, squareform=squareform, rankdata=rankdata, combinations=combinations,
               SEED=20260907, NPERM=1999, results=[], coverage=[], rosters=[], matrices={})
    import_functions(origin(X / "provenance/check_prevalence.py"), ["relative", "logmet", "transform", "evaluate"], env)

    def rd(p, **kw):
        return pd.read_csv(origin(p), sep="\t", **kw)

    def bframe(folder, key):
        a = json.loads(origin(folder / (key + "_axes.json")).read_text()); v = sparse.load_npz(origin(folder / (key + ".npz")))
        return pd.DataFrame(v.toarray(), index=a["samples"], columns=a["features"])

    def collapse(d):
        ids = d.index.to_series().str.extract(r"(MicrobPred_\d+)", expand=False); assert ids.notna().all()
        return d.groupby(ids).mean()

    first = rd(OLDE / "PREDICT_fullstack_primary_samples.tsv", index_col=0); assert len(first) == 77
    unit = {k: f"U{i + 1:02d}" for i, k in enumerate(first.index)}       # pseudonyms in analysis order
    host = collapse(bframe(B / "groups/PREDICT", "supported_human_groups"))
    met = env["logmet"](rd(OLDE / "PREDICT_metabolite_intensity.tsv", index_col=0))
    rows = []
    iu = np.triu_indices(77, 1)

    def emit(level, layer, variant, mat):
        dist = pdist(mat)
        rows.extend(dict(level=level, layer=layer, variant=variant, unit_i=f"U{a + 1:02d}", unit_j=f"U{b + 1:02d}",
                         distance=x) for a, b, x in zip(iu[0], iu[1], dist))

    for lev in ["KO", "OG", "Pathway"]:
        fs = dict(metaG=rd(OLDE / ("metaG_OG_TPM.tsv" if lev == "OG" else "PREDICT_metaG_" + lev + ".tsv"), index_col=0),
                  metaT=rd(OLDE / ("metaT_OG_TPM.tsv" if lev == "OG" else "PREDICT_metaT_" + lev + ".tsv"), index_col=0),
                  host_metaP=host, metabolites=met)
        for k, v in fs.items():
            emit(lev, k, "all", env["transform"](v.loc[first.index], k, 0.0, "PREDICT", lev, "displayed77"))
        for v in ["equal", "unique", "complete", "unique_complete", "NEGATIVE_CONTROL_protein_patients_permuted"]:
            vv = "equal" if v.startswith("NEGATIVE") else v
            prot = collapse(frame(ALLOC / "out", "PREDICT_noclust", f"{lev}_other_unresolved__{vv}"))
            if v.startswith("NEGATIVE"):
                rng = np.random.default_rng(20261003); sub = prot.loc[first.index]
                prot = pd.DataFrame(sub.to_numpy()[rng.permutation(len(sub))], index=sub.index, columns=sub.columns)
            emit(lev, "microbial_metaP", v, env["transform"](prot.loc[first.index], "microbial_metaP", 0.0, "PREDICT", lev, "displayed77"))
        print("PREDICT distances", lev, flush=True)
    assert all(u.startswith("U") for u in unit.values())
    write(pd.DataFrame(rows), "allocation_predict_distances.tsv")


def hstool_pairs():
    """Per-pair Pearson r of log2 KO intensity (hstool_precision.py:pairs via alloc_sensitivity.py)."""
    hp = import_functions(origin(C / "hstool_precision.py"), ["measure", "pairs", "summarize"],
                          dict(np=np, pd=pd, METRICS=["r_per_pair", "cv_per_pair_pct", "n_features_per_pair"]))
    design = pd.read_csv(origin(C / "hstool_precision_refit/design.tsv"), sep="\t", dtype=str)
    tiers = ["TechnicalReplicates", "WorkflowReplicates", "OnePiece", "Intra", "Interindividual"]
    groups = {t: design.loc[design.group.eq(t), "sample_id"].tolist() for t in tiers}
    pairlist = [(t, a, b) for t in tiers for a, b in combinations(groups[t], 2)]; assert len(pairlist) == 316
    samples = design.sample_id.tolist(); out = []
    for pname in ["all", "other_unresolved"]:
        for v in ["equal", "unique", "complete", "unique_complete"]:
            f = frame(ALLOC / "out", "HSTOOL_noclust", f"KO_{pname}__{v}").T.reindex(columns=samples)
            out.append(hp["pairs"](f, pairlist).assign(profile=pname, variant=v))
    rng = np.random.default_rng(20261003)
    f = frame(ALLOC / "out", "HSTOOL_noclust", "KO_other_unresolved__equal").T.reindex(columns=samples)
    sh = f.apply(lambda c: pd.Series(rng.permutation(c.to_numpy()), index=c.index))
    out.append(hp["pairs"](sh, pairlist).assign(profile="other_unresolved", variant="NEGATIVE_CONTROL_feature_shuffle"))
    out = pd.concat(out, ignore_index=True)
    ref = pd.read_csv(ALLOC / "out/hstool_KO_precision.tsv", sep="\t")
    s = out.groupby(["profile", "variant", "design"], sort=False)[["r_per_pair", "cv_per_pair_pct", "n_features_per_pair"]].median().reset_index()
    m = s.merge(ref, on=["profile", "variant", "design"], suffixes=("", "_ref"))
    assert len(m) == len(ref) == 45
    for c in ["r_per_pair", "cv_per_pair_pct", "n_features_per_pair"]:
        assert np.allclose(m[c], m[c + "_ref"], rtol=1e-12, atol=1e-12), c
    RECEIPT["checks"]["hstool_pairs_reproduce_producer_medians"] = True
    write(out[["profile", "variant", "design", "i", "j", "r_per_pair", "cv_per_pair_pct", "n_features_per_pair"]],
          "allocation_hstool_ko_pairs.tsv")


def main():
    SD.mkdir(parents=True, exist_ok=True)
    steps = [parser_audit, completeness, hstool97, origin_rows, annotation_coverage, campi_supplement, campi_targets,
             alloc_coverage, own_stool_distances, predict_distances, hstool_pairs, mgnify]
    only = sys.argv[1:]
    for s in steps:
        if not only or s.__name__ in only:
            s()
    for name in ["Fig2/alphanovo_accuracy_per_acquisition.tsv", "Fig2/reservoir_protein_presence_by_proteome.tsv",
                 "SFig02/bags_bounded48_per_acquisition.tsv", "SFig02/bags_bounded48_orderings.tsv", "Fig1/stage1_stats.json"]:
        RECEIPT["checks"]["shared:" + name] = sha(REPO / "source_data" / name)
    name = "prepare_receipt.json" if not only else "prepare_receipt_" + "_".join(only) + ".json"
    (UNIT / name).write_text(json.dumps(RECEIPT, indent=1, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
