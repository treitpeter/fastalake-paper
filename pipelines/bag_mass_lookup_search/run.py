"""Annotation replays and explicit raw-spectrum coverage for frozen comparisons."""
from pathlib import Path
from collections import defaultdict
import csv
import json
import math
import os
import subprocess
import sys
import time
import traceback

from pyteomics import mzml
from ions import peptide_ions, legacy_ics, backbone, relation, site_metrics, self_test

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from analyze import canonical, identity, dump, tsv

ARMS = ("exact", "bags", "mass", "combined")


def rows(path):
    with Path(path).open() as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def row_index(path):
    values = {}
    for row in rows(path):
        assert int(row["rank"]) == 1
        key = (row["filename"], row["scannr"])
        assert key not in values
        values[key] = row
    return values


def main(index):
    manifest = json.loads((HERE / "MANIFEST.json").read_text())
    for path, expected in manifest["code_sha256"].items():
        assert identity(path)["sha256"] == expected, path
    assert identity(manifest["sage"]["path"])["sha256"] == manifest["sage"]["sha256"]
    self_test()
    case = manifest["cases"][index]
    source_stats = {}
    for expected in case["frozen"]:
        path = Path(expected["path"])
        assert identity(path)["sha256"] == expected["sha256"], str(path)
        source_stats[str(path)] = (path.stat().st_size, path.stat().st_mtime_ns)
    source = Path(case["directory"])
    output = HERE / "runs" / str(index)
    output.mkdir(parents=True, exist_ok=False)
    try:
        task = json.loads((source / "task.json").read_text())
        assert len(task["mzml"]) == 1
        raw = Path(task["mzml"][0]["path"])
        stat = raw.stat()
        assert (stat.st_size, stat.st_mtime_ns) == (task["mzml"][0]["bytes"], task["mzml"][0]["mtime_ns"])
        frozen_rows = [r for r in rows(HERE.parent / "results/PEPTIDE_MECHANISMS.tsv")
                       if (r["cohort"], r["sample"]) == (case["cohort"], case["sample"])]
        changed = {r["peptide_IL"] for r in frozen_rows}
        original = {arm: row_index(source / "searches" / arm / "results.sage.tsv") for arm in ARMS}
        wanted = {key for arm in ARMS for key, row in original[arm].items()
                  if int(row["label"]) == 1 and float(row["peptide_q"]) <= 0.01 and canonical(row["peptide"]) in changed}
        assert wanted and {name for name, _ in wanted} == {raw.name}
        seen, replays, results, ions, annotations = {}, {}, {}, {}, {}
        env = dict(os.environ, RAYON_NUM_THREADS="4", OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
        selection = json.loads((source / "selection.json").read_text())
        for arm in ARMS:
            fa = source / (arm + ".fasta")
            fingerprint = identity(fa)["sha256"]
            assert fingerprint == selection["arms"][arm]["fasta"]["sha256"]
            run = output / arm
            if fingerprint in seen:
                run.symlink_to(seen[fingerprint], target_is_directory=True)
            else:
                run.mkdir()
                settings = json.loads((source / "searches" / arm / "config.json").read_text())
                settings["database"]["fasta"] = str(fa)
                settings["output_directory"] = str(run)
                dump(run / "config.json", settings)
                args = [manifest["sage"]["path"], str(run / "config.json"), "--write-pin", "--annotate-matches", "--disable-telemetry-i-dont-want-to-improve-sage"]
                dump(run / "COMMAND.json", args)
                with (run / "stdout.log").open("w") as log, (run / "stderr.log").open("w") as error:
                    subprocess.run(["/usr/bin/time", "-v", *args], check=True, stdout=log, stderr=error, env=env)
                seen[fingerprint] = arm
            replay = row_index(run / "results.sage.tsv")
            assert set(replay) == set(original[arm]), "Replay spectrum roster differs"
            raw_differences = 0
            q_differences = 0
            max_q_delta = 0.0
            for key, before in original[arm].items():
                after = replay[key]
                if (before["peptide"], before["label"]) != (after["peptide"], after["label"]) or not math.isclose(float(before["hyperscore"]), float(after["hyperscore"]), rel_tol=1e-7, abs_tol=1e-5):
                    raw_differences += 1
                delta = abs(float(before["peptide_q"]) - float(after["peptide_q"]))
                if delta > 1e-8:
                    q_differences += 1
                    max_q_delta = max(max_q_delta, delta)
            replays[arm] = dict(psms=len(replay), raw_assignment_differences=raw_differences,
                                peptide_q_differences=q_differences, max_q_delta=max_q_delta)
            dump(output / "REPLAY_COMPARISON.json", replays)
            assert raw_differences == 0, "Annotation replay changed top-ranked sequence/label/raw score"
            ids = {}
            for key in wanted & replay.keys():
                row = original[arm][key]
                sequence, neutral, theoretical = peptide_ions(row["peptide"])
                assert abs(neutral - float(row["calcmass"])) < 0.02, "Modification/mass parser disagrees with Sage"
                ions[(arm, key)] = (sequence, theoretical)
                results[(arm, key)] = dict(arm=arm, filename=key[0], scannr=key[1], peptide=row["peptide"],
                     label=int(row["label"]), peptide_q=float(row["peptide_q"]), spectrum_q=float(row["spectrum_q"]),
                     hyperscore=float(row["hyperscore"]), matched_peaks=int(row["matched_peaks"]), length=len(sequence))
                ids[replay[key]["psm_id"]] = key
            for fragment in rows(run / "matched_fragments.sage.tsv"):
                if fragment["psm_id"] in ids:
                    key = ids[fragment["psm_id"]]
                    ion = (fragment["fragment_type"], int(fragment["fragment_ordinals"]), int(fragment["fragment_charge"]))
                    if ion[0] in ("b", "y"):
                        seq, theoretical = ions[(arm, key)]
                        assert ion in theoretical, "Unexpected fragment charge or ordinal"
                        assert abs(theoretical[ion] - float(fragment["fragment_mz_calculated"])) < 0.02
                        annotations.setdefault((arm, key), set()).add(ion)
            for key in wanted & replay.keys():
                seq, _ = ions[(arm, key)]
                matched = annotations.get((arm, key), set())
                results[(arm, key)]["sage_backbone_coverage"] = len(backbone(matched, len(seq))) / (len(seq)-1)
                results[(arm, key)]["sage_unique_by_coverage"] = len({(s, n) for s, n, z in matched}) / (2 * (len(seq)-1))
            print(json.dumps(dict(stage="REPLAY_COMPLETE", arm=arm, **replays[arm])), flush=True)
        remaining = {scan for _, scan in wanted}
        with mzml.MzML(str(raw), use_index=False, iterative=True, read_schema=False) as spectra:
            for spectrum in spectra:
                scan = spectrum["id"]
                if scan not in remaining:
                    continue
                assert int(spectrum["ms level"]) == 2
                key = (raw.name, scan)
                observed = spectrum["m/z array"]
                for arm in ARMS:
                    if (arm, key) in results:
                        seq, theoretical = ions[(arm, key)]
                        coverage, hits = legacy_ics(theoretical, observed, len(seq))
                        results[(arm, key)]["legacy_ics_raw_z1_20ppm"] = coverage
                        results[(arm, key)]["raw_backbone_coverage_z1"] = len(backbone(hits, len(seq))) / (len(seq)-1)
                remaining.remove(scan)
                if not remaining:
                    break
        assert not remaining, f"Missing raw spectra: {len(remaining)}"
        paired = []
        accepted_keys = defaultdict(list)
        for arm in ARMS:
            for key in wanted & original[arm].keys():
                row = original[arm][key]
                if int(row["label"]) == 1 and float(row["peptide_q"]) <= 0.01:
                    accepted_keys[(arm, canonical(row["peptide"]))].append(key)
        for row in frozen_rows:
            upper, lower = row["contrast"].split("_vs_")
            own, other = (upper, lower) if row["change"] == "gained" else (lower, upper)
            peptide = row["peptide_IL"]
            for key in accepted_keys[(own, peptide)]:
                a = results[(own, key)]
                b = results.get((other, key))
                pair = dict(cohort=case["cohort"], sample=case["sample"], contrast=row["contrast"], change=row["change"],
                            category=row["category"], filename=key[0], scannr=key[1], own_arm=own, other_arm=other,
                            own_peptide=a["peptide"], other_peptide=b["peptide"] if b else "",
                            other_label=b["label"] if b else "", counterpart_status="unreported" if b is None else "reported_top1_only",
                            relation="unreported", pairwise_site_ions="", pairwise_diagnostic_site_ions="")
                for metric in ("hyperscore", "legacy_ics_raw_z1_20ppm", "sage_backbone_coverage", "sage_unique_by_coverage", "peptide_q", "spectrum_q"):
                    pair["own_"+metric] = a[metric]
                    pair["other_"+metric] = b[metric] if b else ""
                    pair["delta_"+metric] = a[metric]-b[metric] if b else ""
                if b:
                    sa, ia = ions[(own, key)]
                    sb, ib = ions[(other, key)]
                    relationship, positions = relation(sa, sb, b["label"])
                    pair["relation"] = relationship
                    if positions:
                        spans, diagnostic = site_metrics(annotations.get((own, key), set()), ia, ib, len(sa), positions)
                        pair["pairwise_site_ions"] = spans
                        pair["pairwise_diagnostic_site_ions"] = diagnostic
                paired.append(pair)
        tsv(output / "PSM_COVERAGE.tsv", [dict(cohort=case["cohort"], sample=case["sample"], **r) for r in results.values()])
        tsv(output / "PAIRED_EVIDENCE.tsv", paired)
        after = raw.stat()
        assert (stat.st_size, stat.st_mtime_ns) == (after.st_size, after.st_mtime_ns)
        for path, expected in source_stats.items():
            stat_after = Path(path).stat()
            assert expected == (stat_after.st_size, stat_after.st_mtime_ns), path
        dump(output / "VALIDATION.json", dict(status="PASS", job_id=os.environ.get("SLURM_JOB_ID"),
             cohort=case["cohort"], sample=case["sample"], raw_spectra=len(wanted), scored_psms=len(results),
             paired_observations=len(paired), replays=replays, raw_input=task["mzml"][0],
             limits="Top-one comparator is not guaranteed to be the corresponding canonical peptide; local sequence differences are candidates only; no variant-subset FDR calibration."))
        print(json.dumps(dict(status="PASS", index=index, spectra=len(wanted))), flush=True)
    except BaseException as error:
        dump(output / "FAILURE.json", dict(error=str(error), traceback=traceback.format_exc()))
        raise


if __name__ == "__main__":
    main(int(sys.argv[1]))
