"""Freeze every input of task C with sha256 (bytes + mtime_ns for the 61 GB reservoir, whose
full sha256 is computed inside the C2 pilot job on Slurm and written to inputs/reservoir.sha256.json)."""
import os
import hashlib, json, os, sys, datetime
R = os.environ["FASTALAKE_OUT_DIR"]
FL = os.environ["FASTALAKE_DATA_DIR"]

def sha(p):
    """sha256 hex digest of a file."""
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 24), b""):
            h.update(b)
    return h.hexdigest()

def meta(p, digest=True):
    """Size, mtime_ns and (optionally) sha256 of an input file."""
    st = os.stat(p)
    d = {"path": p, "bytes": st.st_size, "mtime_ns": st.st_mtime_ns}
    if digest:
        d["sha256"] = sha(p)
    return d

M = {
    "task": "R8 task C: C1 AlphaNovo reproduction comparison; C2 adjacent-swap lookup pilot on the full HSTOOL reservoir",
    "written_utc": datetime.datetime.utcnow().isoformat() + "Z",
    "canonical_commit": "66f1e6ae233f6cc2bd60199214a205566ea11f78",
    "C1": {
        "july_predictions": meta(os.environ["FASTALAKE_PHD_DIR"] + "/Peter_DIA/Peter_Sophia_FlashLFQ/mzML_files/alphanovo_inference/20191024_QX7_PeTr_LC02_SA_MicrobPred_10_191025211701_predictions.csv"),
        "repro_predictions_job5950191": meta(FL + "/fable_5.1_review_2026-09-16/runs/an_repro/predictions/20191024_QX7_PeTr_LC02_SA_MicrobPred_10_191025211701_predictions.csv"),
        "compare_py_original": meta(FL + "/fable_5.1_review_2026-09-16/runs/an_repro/compare.py"),
        "compare_py_copy": meta(R + "/scripts/compare.py"),
        "infer_log_5950191": meta(FL + "/fable_5.1_review_2026-09-16/runs/an_repro/infer_5950191.out"),
        "python": os.environ.get("FASTALAKE_PYTHON", "python") + " (pandas 2.3.3)",
    },
    "C2": {
        "sample": "20181231_QX7_PhGe_LC95_SA_WorkflowReplicates_2",
        "real_queries_source": meta(FL + "/MM_review_analysis/A02_4_2_lookup_specificity/derived/REAL_queries.tsv"),
        "real_queries_frozen": meta(R + "/inputs/REAL_queries.tsv"),
        "sim_queries_source": meta(FL + "/MM_review_analysis/A02_4_2_lookup_specificity/derived/SIM_20260913_queries.tsv"),
        "sim_queries_frozen": meta(R + "/inputs/SIM_20260913_queries.tsv"),
        "reservoir": meta(FL + "/canonical/publication_data/cohorts/HSTOOL/0_lake/HSTOOL_lake.fasta", digest=False),
        "extractor_source": meta(FL + "/canonical/rust/fasta_extractor_v2/target/release/fasta_extractor"),
        "extractor_frozen": meta(R + "/bin/fasta_extractor"),
        "parsimony_source": meta(FL + "/canonical/rust/parsimony_engine/target/release/parsimony_engine"),
        "parsimony_frozen": meta(R + "/bin/parsimony_engine"),
        "a02_lookup_reference": meta(FL + "/MM_review_analysis/A02_4_2_lookup_specificity/inputs/REAL_sample_lake.stats.json"),
        "a02_accounting": meta(FL + "/MM_review_analysis/A02_4_2_lookup_specificity/status/slurm_accounting.tsv"),
        "pilot_seed": 20260917,
        "pilot_n": 1000,
        "scramble_seed_used": 20260913,
    },
}
# pairs that must be byte-identical
for a, b in [("C1.compare_py_original", "C1.compare_py_copy"), ("C2.real_queries_source", "C2.real_queries_frozen"),
             ("C2.sim_queries_source", "C2.sim_queries_frozen"), ("C2.extractor_source", "C2.extractor_frozen"),
             ("C2.parsimony_source", "C2.parsimony_frozen")]:
    ga = M[a.split(".")[0]][a.split(".")[1]]; gb = M[b.split(".")[0]][b.split(".")[1]]
    assert ga["sha256"] == gb["sha256"], (a, b)
M["identical_pairs_checked"] = 5
out = R + "/inputs/manifest.json"
with open(out, "w") as f:
    json.dump(M, f, indent=2)
print("wrote", out)
