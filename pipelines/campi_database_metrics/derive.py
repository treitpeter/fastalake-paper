#!/usr/bin/env python3
"""Run the existing Figure 2 derivations with the five explicitly named database arms."""

import argparse
import csv
import hashlib
import importlib
import json
import subprocess
import sys
from pathlib import Path

CANON = Path(__file__).resolve().parents[3]
SRC = CANON / "publication/analyses/fig2_revision_2026-08-11/src"
sys.path.insert(0, str(SRC))
import derive_fig2_revision as D  # noqa: E402


def sha(path):
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            value.update(chunk)
    return value.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    out = args.out.resolve()
    if out.exists():
        raise FileExistsError(f"Refusing to overwrite existing output: {out}")
    arms = CANON / "publication_data/arms"
    names = {
        "fastalake": ("uniprotkb_2026_03_noclust_2026-09-09", "Blind UniProt noclust"),
        "ref_append": ("uniprotkb_2026_03_noclust_campi_ref_append_2026-09-09", "+CAMPI REF"),
        "ref_metag_append": (
            "uniprotkb_2026_03_noclust_campi_ref_metag_append_2026-09-09",
            "+CAMPI REF+MG",
        ),
    }
    D.ARMS = {
        key: (arms / name / "cohorts/CAMPI_MOCK/5_sage/razor_hash-acc", label, name, key)
        for key, (name, label) in names.items()
    }
    D.ARMS.update(
        {
            "campi_ref": (arms / "campi_ref_solo/CAMPI_MOCK", "CAMPI REF", "REF", "off"),
            "campi_mo": (arms / "campi_mo_fixed/CAMPI_MOCK", "CAMPI MG", "MG", "off"),
        }
    )
    for key in ["ref_append", "ref_metag_append"]:
        validation = json.loads(
            (arms / names[key][0] / "SEARCH_COMPARISON_VALIDATED.json").read_text()
        )
        if validation["status"] != "PASS" or validation["n_acquisitions"] != 15:
            raise ValueError(f"Incomplete paired validation for {key}")
    inputs = {
        D.MEMBER_FASTA,
        D.MEMBERS_TSV,
        CANON / "publication/cohorts/campi_mock_sample_list.tsv",
    }
    database_rows = []
    for arm, (path, label, _, _) in D.ARMS.items():
        for sample in D.SAMPLES:
            files = [path / sample / f for f in ["config.json", "results.sage.tsv", "lfq.tsv"]]
            for file in files:
                if not file.is_file() or not file.stat().st_size:
                    raise ValueError(f"Missing input {file}")
            inputs.update(files)
            config = json.loads(files[0].read_text())
            fasta = Path(config["database"]["fasta"]).resolve()
            inputs.add(fasta)
            with fasta.open() as stream:
                count = sum(line.startswith(">") for line in stream)
            database_rows.append([arm, label, sample, str(fasta), count])
    out.mkdir(parents=True)
    D.TAB = out
    provenance = {
        "commit": subprocess.check_output(
            ["git", "-C", str(CANON), "rev-parse", "HEAD"], text=True
        ).strip(),
        "arms": {k: {"label": v[1], "search_root": str(v[0])} for k, v in D.ARMS.items()},
        "input_sha256": {str(p): sha(p) for p in sorted(inputs)},
        "script_sha256": {
            str(p): sha(p)
            for p in [
                Path(__file__),
                SRC / "derive_fig2_revision.py",
                SRC / "derive_upset.py",
                SRC / "derive_length_specificity.py",
            ]
        },
        "scope": "Five arms; 15 acquisitions; existing Figure 2 peptide/LFQ definitions.",
    }
    (out / "PROVENANCE.json").write_text(json.dumps(provenance, indent=2) + "\n")
    with (out / "database_sizes.tsv").open("w") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(["arm", "arm_label", "sample", "searched_fasta", "database_records"])
        writer.writerows(database_rows)
    D.main()
    for module in ["derive_upset", "derive_length_specificity"]:
        importlib.import_module(module).main()
    # Historical comments contain fixed labels; actual source paths/arms are above.
    for name in ["upset_intersections.tsv", "upset_member_totals.tsv", "length_specificity.tsv"]:
        path = out / name
        if path.exists():
            value = path.read_text().replace("four arms", "five arms")
            value = value.replace(
                "noclust_2026-08-07", "see PROVENANCE.json for arm-specific builds"
            )
            path.write_text(value)
    (out / "DERIVATION_COMPLETE.json").write_text(
        json.dumps(
            {
                "status": "PASS",
                "n_arms": 5,
                "n_acquisitions_per_arm": 15,
                "note": "Source tables only; final figure and manuscript pending.",
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
