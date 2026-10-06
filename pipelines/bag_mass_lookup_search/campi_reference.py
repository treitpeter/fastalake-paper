"""Corroborate S06 changes against its archived reference-database Sage search.

This is the same acquisition/engine, not independent molecular ground truth.
CAMPI reference surrogate strains prevent interpreting absence as false discovery.
"""
from pathlib import Path
from collections import Counter
import csv
import json
import os

from analyze import HERE, SOURCES, dump, identity, locate, read, search, tsv


def main():
    out = HERE / "reference_crosscheck"
    out.mkdir(exist_ok=False)
    assert read(HERE / "results/VALIDATION.json")["status"] == "PASS"
    inputs = read(HERE / "INPUTS.json")
    entry = next(item for item in inputs["cases"] if item["cohort"] == "CAMPI" and item["sample"] == "S06")
    case = Path(entry["directory"])
    root = HERE.parents[3]
    previous = read(root / "MM_review_analysis/A03_4_2_near_misses/inputs/manifest.json")
    task = next(item for item in previous["samples"] if item["sample"] == "S06")
    for key in ("ref_search", "ref_config"):
        source = task[key]
        actual = identity(source["path"])
        assert actual["bytes"] == source["bytes"] and actual["mtime_ns"] == source["mtime_ns"]
    reference_dir = Path(task["ref_search"]["path"]).parent
    config = read(reference_dir / "results.json")
    assert config["version"] == "0.14.6"
    reference_fasta = Path(config["database"]["fasta"])
    native_config = read(case / "searches/exact/config.json")
    assert [str(Path(p).resolve()) for p in config["mzml_paths"]] == [str(Path(p).resolve()) for p in native_config["mzml_paths"]]
    for settings in (config, native_config):
        for key in ("version", "output_paths", "output_directory"):
            settings.pop(key, None)
        settings["database"].pop("fasta")
    assert config == native_config, "Search settings differ beyond database/output"
    reference = search(task["ref_search"]["path"])
    source_rows = HERE / "results/PEPTIDE_MECHANISMS.tsv"
    identity(source_rows)
    with source_rows.open() as handle:
        rows = [r for r in csv.DictReader(handle, delimiter="\t") if r["cohort"] == "CAMPI"]
    _, contexts, _ = locate(reference_fasta, {r["peptide_IL"] for r in rows})
    native = {arm: search(case / "searches" / arm / "results.sage.tsv") for arm in ("exact", "bags", "mass", "combined")}
    # Stable spectrum IDs must refer to overlapping observations in this same mzML.
    shared = set(reference["spectra"]) & set(native["exact"]["spectra"])
    assert len(shared) > 10000, "Reference spectrum join insufficient"
    details, summary = [], Counter()
    for row in rows:
        upper, lower = row["contrast"].split("_vs_")
        own = upper if row["change"] == "gained" else lower
        peptide = row["peptide_IL"]
        ref_accepted = peptide in reference["accepted"]
        same = different = no_accepted = 0
        for spectrum in native[own]["accepted"][peptide]:
            counterpart = reference["spectra"].get(spectrum)
            if counterpart and counterpart[1] == 1 and counterpart[2] <= 0.01:
                if counterpart[0] == peptide:
                    same += 1
                else:
                    different += 1
            else:
                no_accepted += 1
        detail = dict(contrast=row["contrast"], change=row["change"], category=row["category"], peptide_IL=peptide,
                      reference_fasta_contains=bool(contexts[peptide]), reference_accepted_any_spectrum=ref_accepted,
                      same_spectrum_reference_agreements=same, different_accepted_reference_assignments=different,
                      no_accepted_reference_assignment=no_accepted)
        details.append(detail)
        key = (row["contrast"], row["change"], row["category"])
        summary[key + ("total",)] += 1
        summary[key + ("reference_fasta_contains",)] += bool(contexts[peptide])
        summary[key + ("reference_accepted_any_spectrum",)] += ref_accepted
        summary[key + ("any_same_spectrum_agreement",)] += same > 0
        summary[key + ("any_different_accepted_assignment",)] += different > 0
    table = []
    for key in sorted({k[:3] for k in summary}):
        table.append(dict(contrast=key[0], change=key[1], category=key[2], **{name: summary[key + (name,)] for name in
                         ("total", "reference_fasta_contains", "reference_accepted_any_spectrum", "any_same_spectrum_agreement", "any_different_accepted_assignment")}))
    tsv(out / "PEPTIDES.tsv", details)
    tsv(out / "SUMMARY.tsv", table)
    for record in SOURCES.values():
        stat = Path(record["path"]).stat()
        assert (stat.st_size, stat.st_mtime_ns) == (record["bytes"], record["mtime_ns"])
    dump(out / "SOURCE_IDENTITIES.json", list(SOURCES.values()))
    dump(out / "VALIDATION.json", dict(status="PASS", job_id=os.environ.get("SLURM_JOB_ID"),
         sample="S06", shared_reported_spectra=len(shared), peptide_comparisons=len(details),
         interpretation="Same-acquisition same-engine reference-database agreement; not independent truth or FDR calibration; reference contains surrogate strains."))
    print(json.dumps(dict(status="PASS", shared_spectra=len(shared), output=str(out))))


if __name__ == "__main__":
    main()
