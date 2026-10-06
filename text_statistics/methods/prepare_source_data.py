#!/usr/bin/env python3
"""Extract the source data of the Methods statistics unit from the lab analysis tree.

Run once, by the authors, on Slurm. It is not part of the notebook. It documents exactly
how every table in ``source_data/methods/`` was cut from the original outputs. It needs
read access to the lab project directory, given as the environment variable
``FASTALAKE_LAB_ROOT`` (written ``LAB_ROOT`` in PROVENANCE.md); every input path below is
relative to it.

Tables written (see ``source_data/methods/README.md`` for the data dictionary)
------------------------------------------------------------------------------
reservoir_sources.tsv
    The SOURCES BY PRIORITY block of the four lake-builder statistics files (public,
    PREDICT cohort, CapScan cohort and CAMPI reservoirs), one row per source.
campi_comparator_records.tsv
    One row per record of the two supplied CAMPI comparator FASTAs (reference and MG), as
    frozen for the comparator arms: accession, length and SHA-256 of the sequence.
union_contributions.tsv
    Independent recount of the exact-sequence unions of the selected FASTAs (healthy
    stool 58, PREDICT 306, CapScan 331 and 347). One row per FASTA: records, FASTA
    sha256 and the number of sequences that FASTA adds to the union when FASTAs are read
    in acquisition-name order. The union size is the sum of the additions. The recount
    re-implements the identity rule of the original ``build_union`` (exact full sequence,
    whitespace removed, upper case, I/L kept distinct) without importing it, and is
    checked against the original ``union_summary.json`` files.
crossfit_lda_fallback_347.tsv
    Roster index, cohort and ``lda_fallback`` of the 578 matched acquisitions of the
    10-fold cross-fit rerun on the 347-acquisition CapScan roster.
split_spectra_mg_per_acquisition.tsv
    Source Data package sheet ``Supplementary_Figure_06/metagenome_same_spectra/
    per_acquisition.tsv`` without the acquisition file name and without the PREDICT
    specimen-visit code (column ``group`` is kept for CapScan sampling regions only).
split_spectra_mg_summary_original.tsv
    Byte-identical copy of the package sheet ``.../metagenome_same_spectra/summary.tsv``
    (the original bootstrap summary, used as the expected values).
predict_reconversion.tsv
    The 306 PREDICT public-sequence searches of the Sage input-completeness check: whether
    the searched mzML was the reconverted file, and the end-of-file check of the original
    mzML in the deposit folder.
compute_nodes.tsv
    ``scontrol show node`` records (CPUs, memory, features, generic resources) of the
    nodes that ran the CAMPI reservoir build and the CAMPI AlphaNovo inference.
crossfit_interpreter.tsv
    Python and NumPy versions of the interpreter named in the cross-fit worker launcher.
    The interpreter is passed as ``FASTALAKE_CROSSFIT_PYTHON``.

The controlled clinical input ``paired_antibiotic_profiles.tsv`` is copied unchanged to
``--controlled-out`` (never into the repository).

Usage (on Slurm)::

    FASTALAKE_LAB_ROOT=... FASTALAKE_CROSSFIT_PYTHON=... \\
        python text_statistics/methods/prepare_source_data.py \\
        --out source_data/methods --controlled-out <dir> --receipt <file.json>
"""

import argparse
import csv
import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
from pathlib import Path

import pandas as pd

#: Analysis roots, relative to LAB_ROOT.
MC = Path("review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01")
CR = MC / "biology_refit/claude_review_2026-10-03"
PKG = CR / "final/figure_sd/source_data/package/FastaLake_Source_Data"
A01 = Path("MM_review_analysis/A01_4_1_database_sizes/inputs")
ARMS = Path("canonical/publication_data/arms")
MERGE = Path("sample_specific_merge_all_raw/2026-09-11")

#: Lake-builder statistics files, one per reservoir.
RESERVOIR_STATS = {
    "public": A01 / "HSTOOL_reservoir_stats.txt",
    "PREDICT_cohort": A01 / "PREDICT_reservoir_stats.txt",
    "CapScan_cohort": A01 / "CAPSCAN_reservoir_stats.txt",
    "CAMPI": ARMS / "uniprotkb_2026_03_noclust_2026-09-09/cohorts/CAMPI_MOCK/0_lake/CAMPI_MOCK_lake_stats.txt",
}

#: Supplied CAMPI comparator FASTAs as frozen for the comparator arms, with their sha256.
COMPARATOR_FASTAS = {
    "reference": (ARMS / "uniprotkb_2026_03_noclust_campi_ref_append_2026-09-09/inputs/campi_ref_targets.faa",
                  "9d4def2d1de49d61e88f163e8c456da891b407d77e8959f58da9dbead6d9fce9"),
    "MG": (ARMS / "uniprotkb_2026_03_noclust_campi_ref_metag_append_2026-09-09/inputs/campi_mo_targets_fixed.faa",
           "364c64150e23b8eeabbeaf9c0158bac915f87645dfbbc18710a34af6aea6c8d0"),
}

#: Union rosters: manifest task id, original union summary, extra FASTAs.
UNION_TASKS = {
    "HSTOOL58": ("HSTOOL_noclust", MERGE / "cohorts/HSTOOL_noclust/union_summary.json"),
    "PREDICT306": ("PREDICT_noclust", MERGE / "cohorts/PREDICT_noclust/union_summary.json"),
    "CAPSCAN331": ("CAPSCAN_noclust", MERGE / "cohorts/CAPSCAN_noclust/union_summary.json"),
    "CAPSCAN347": ("CAPSCAN_noclust", CR / "final/reruns_347/A7_other/union/fixed347/union_summary.json"),
}
P15_ROSTER = CR / "coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/roster/participant15_added.tsv"
P15_FASTA_DIR = ARMS / "noclust_2026-08-07/cohorts/CAPSCAN/4_parsimony/razor_hash-acc"
UNION_DECOY_TAG = "rev_"   # database.decoy_tag of the union task's sage.json

CROSSFIT347 = CR / "final/public347/crossfit_capscan347/run/PER_ACQUISITION.tsv"
SPLIT_MG = PKG / "Supplementary_Figure_06/metagenome_same_spectra/per_acquisition.tsv"
SPLIT_MG_SUMMARY = PKG / "Supplementary_Figure_06/metagenome_same_spectra/summary.tsv"
COVERAGE = CR / "final/completeness_653/out/search_coverage.tsv"
ENDINGS = CR / "final/second_engine_all/truncation_check/predict_mzml_endings.tsv"
CLIN_PAIRED = (MC / "clinical_refit/paired_antibiotic_profiles.tsv",
               "b99dd07a5cf75150f7cb6520eedcc816cefcd6e562dea9b457c309dfa95f27ec")
CAMPI_LAKE_JOB_NODE = "hpcl9202"        # NodeList of job 5895248 (source_data/SFig15/compute_task_records.tsv)
CAMPI_INFERENCE_ARRAY = "5450700"       # AlphaNovo CAMPI inference array (source_data/SFig15)


def sha256(path):
    """Return the hex sha256 of a file."""
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 23), b""):
            h.update(block)
    return h.hexdigest()


def lab_relative(path_text):
    """Return a manifest path as a path relative to LAB_ROOT (the part after 'FastaLake/')."""
    i = path_text.find("/FastaLake/")
    if i < 0:
        raise ValueError("path outside the lab project directory")
    return Path(path_text[i + len("/FastaLake/"):])


def fasta_records(path):
    """Yield ``(header, sequence)``; sequence has whitespace removed and is upper case."""
    header, parts = None, []
    with open(path, "rt") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(parts).upper()
                header, parts = line[1:], []
            else:
                parts.append("".join(line.split()))
    if header is not None:
        yield header, "".join(parts).upper()


def write(df, path):
    """Write a tab-separated table without index."""
    df.to_csv(path, sep="\t", index=False)


def reservoir_sources(lab):
    """Parse the SOURCES BY PRIORITY block of each lake-builder statistics file."""
    rows = []
    pat = re.compile(r"^(\S+)\s+(\d+)\s+(\d+)\s+([\d,]+)\s+([\d,]+)\s+([\d,]+)\s+([\d.]+)s$")
    for reservoir, rel in RESERVOIR_STATS.items():
        text = (lab / rel).read_text().splitlines()
        threads = int(next(l for l in text if l.startswith("Threads:")).split()[1])
        block = text[text.index("SOURCES BY PRIORITY") + 3:]
        for line in block:
            if line.startswith("---"):
                break
            m = pat.match(line.strip())
            assert m, line
            rows.append(dict(reservoir=reservoir, source=m[1], priority=int(m[2]), files=int(m[3]),
                             input_records=int(m[4].replace(",", "")), unique_sequences=int(m[5].replace(",", "")),
                             duplicate_records=int(m[6].replace(",", "")), source_seconds=float(m[7]),
                             builder_threads=threads))
    return pd.DataFrame(rows)


def comparator_records(lab):
    """One row per record of the frozen CAMPI comparator FASTAs."""
    rows = []
    for db, (rel, expected) in COMPARATOR_FASTAS.items():
        assert sha256(lab / rel) == expected, rel
        for k, (header, seq) in enumerate(fasta_records(lab / rel)):
            rows.append(dict(database=db, record_index=k, accession=header.split()[0], length=len(seq),
                             sequence_sha256=hashlib.sha256(seq.encode("ascii")).hexdigest()))
    return pd.DataFrame(rows)


def union_rosters(lab):
    """Return ``{roster: [(acquisition, lab-relative FASTA path), ...]}``."""
    manifest = json.loads((lab / MERGE / "manifest.json").read_text())
    tasks = {t["id"]: t for t in manifest["tasks"]}
    out = {}
    for roster, (task_id, _) in UNION_TASKS.items():
        items = [(lab_relative(f["path"]).parent.name, lab_relative(f["path"])) for f in tasks[task_id]["fastas"]]
        if roster == "CAPSCAN347":
            with open(lab / P15_ROSTER) as fh:
                added = [r["sample_id"] for r in csv.DictReader(fh, delimiter="\t")]
            assert len(added) == 16
            items += [(s, P15_FASTA_DIR / s / "parsimony.fasta") for s in added]
        out[roster] = sorted(items)
    return out


def union_contributions(lab):
    """Recount each union and return the per-FASTA table and a check dictionary."""
    rows, checks = [], {}
    for roster, items in union_rosters(lab).items():
        seen = set()
        n_records = 0
        for order, (acq, rel) in enumerate(items):
            before = len(seen)
            records = 0
            for header, seq in fasta_records(lab / rel):
                assert seq and re.fullmatch("[A-Z*]+", seq), (rel, header)
                assert UNION_DECOY_TAG not in header.split()[0], (rel, header)
                seen.add(hashlib.sha256(seq.encode("ascii")).digest())
                records += 1
            rows.append(dict(roster=roster, order=order, acquisition=acq, records=records,
                             new_unique_sequences=len(seen) - before, fasta_sha256=sha256(lab / rel)))
            n_records += records
        summary = json.loads((lab / UNION_TASKS[roster][1]).read_text())
        checks[roster] = dict(fastas=len(items), records=n_records, unique=len(seen),
                              original_records=summary["input_records"], original_unique=summary["unique_sequences"],
                              equal=(n_records == summary["input_records"] and len(seen) == summary["unique_sequences"]))
    return pd.DataFrame(rows), checks


def crossfit_fallback(lab):
    """Roster index, cohort and LDA fallback of the 347-roster cross-fit rerun."""
    x = pd.read_csv(lab / CROSSFIT347, sep="\t", dtype=str, keep_default_na=False)
    return x[["i", "cohort", "lda_fallback"]]


def split_mg(lab):
    """Package sheet without acquisition names and without PREDICT specimen-visit codes."""
    x = pd.read_csv(lab / SPLIT_MG, sep="\t", dtype=str, keep_default_na=False)
    x.loc[x.cohort.eq("PREDICT"), "group"] = ""
    return x.drop(columns=["sample"])


def predict_reconversion(lab):
    """PREDICT public-sequence searches: reconverted flag and original-file end check."""
    c = pd.read_csv(lab / COVERAGE, sep="\t", dtype=str, keep_default_na=False)
    c = c[c.lineage.eq("primary_public_current") & c.cohort.eq("PREDICT")]
    keep = ["sample", "repaired", "closing_tag", "nul_in_tail", "declared_spectra", "max_scan",
            "max_scan_fraction", "software"]
    c = c[keep].rename(columns={"sample": "acquisition", "repaired": "searched_reconverted_mzml",
                                "closing_tag": "searched_mzml_closing_tag", "software": "searched_mzml_software"})
    e = pd.read_csv(lab / ENDINGS, sep="\t", header=None, names=["acquisition", "original_mzml_bytes",
                                                                "original_mzml_ending"], dtype=str)
    out = c.merge(e, on="acquisition", how="left", validate="one_to_one")
    assert len(out) == 306 and out.original_mzml_ending.notna().all()
    return out.sort_values("acquisition")


def compute_nodes(lab_sfig15):
    """scontrol records of the CAMPI reservoir node and the CAMPI inference nodes."""
    sacct = pd.read_csv(lab_sfig15, sep="\t", dtype=str)
    jobs = sacct[sacct.JobID.str.match(rf"^{CAMPI_INFERENCE_ARRAY}_\d+$")]
    nodes = [(CAMPI_LAKE_JOB_NODE, "CAMPI reservoir build, job 5895248")]
    nodes += [(n, f"CAMPI AlphaNovo inference, array {CAMPI_INFERENCE_ARRAY}") for n in sorted(set(jobs.NodeList))]
    rows = []
    for node, role in nodes:
        text = subprocess.run(["scontrol", "show", "node", node], capture_output=True, text=True, check=True).stdout
        kv = dict(re.findall(r"(\w+)=(\S*)", text))
        rows.append(dict(node=node, role=role, CPUTot=kv.get("CPUTot"), RealMemory_MB=kv.get("RealMemory"),
                         AvailableFeatures=kv.get("AvailableFeatures"), Gres=kv.get("Gres"),
                         queried_utc=pd.Timestamp.now(tz='UTC').isoformat()))
    return pd.DataFrame(rows)


def crossfit_interpreter():
    """Python and NumPy versions of the cross-fit worker interpreter."""
    py = os.environ["FASTALAKE_CROSSFIT_PYTHON"]
    out = subprocess.run([py, "-c", "import sys, numpy; print(sys.version.split()[0]); print(numpy.__version__)"],
                         capture_output=True, text=True, check=True).stdout.split()
    return pd.DataFrame([dict(launcher="iteration4/crossfit/design/run.sbatch", python=out[0], numpy=out[1])])


def main():
    """Write all tables and a receipt with input and output sha256."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--controlled-out", required=True)
    ap.add_argument("--receipt", required=True)
    ap.add_argument("--sfig15-sacct", required=True, help="source_data/SFig15/campi_alphanovo_inference_sacct.tsv")
    a = ap.parse_args()
    lab = Path(os.environ["FASTALAKE_LAB_ROOT"])
    out, ctl = Path(a.out), Path(a.controlled_out)
    out.mkdir(parents=True, exist_ok=True)
    ctl.mkdir(parents=True, exist_ok=True)
    receipt = dict(inputs={}, outputs={})
    for rel in list(RESERVOIR_STATS.values()) + [CROSSFIT347, SPLIT_MG, SPLIT_MG_SUMMARY, COVERAGE, ENDINGS] + \
            [r for _, r in UNION_TASKS.values()]:
        receipt["inputs"][str(rel)] = sha256(lab / rel)
    write(reservoir_sources(lab), out / "reservoir_sources.tsv")
    write(comparator_records(lab), out / "campi_comparator_records.tsv")
    uc, checks = union_contributions(lab)
    write(uc, out / "union_contributions.tsv")
    receipt["union_checks"] = checks
    write(crossfit_fallback(lab), out / "crossfit_lda_fallback_347.tsv")
    write(split_mg(lab), out / "split_spectra_mg_per_acquisition.tsv")
    shutil.copyfile(lab / SPLIT_MG_SUMMARY, out / "split_spectra_mg_summary_original.tsv")
    write(predict_reconversion(lab), out / "predict_reconversion.tsv")
    write(compute_nodes(a.sfig15_sacct), out / "compute_nodes.tsv")
    write(crossfit_interpreter(), out / "crossfit_interpreter.tsv")
    src, expected = CLIN_PAIRED
    assert sha256(lab / src) == expected
    shutil.copyfile(lab / src, ctl / "paired_antibiotic_profiles.tsv")
    for p in sorted(out.glob("*.tsv")):
        receipt["outputs"][p.name] = sha256(p)
    Path(a.receipt).write_text(json.dumps(receipt, indent=1) + "\n")
    print(json.dumps(checks, indent=1))


if __name__ == "__main__":
    main()
