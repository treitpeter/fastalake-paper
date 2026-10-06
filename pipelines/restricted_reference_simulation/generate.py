"""Construct matched mixtures and independent spectral perturbations; no search."""

import copy
import importlib.metadata
import json
import os
import random
import re
import sys
import time

from common import (
    ACQUISITIONS,
    N_SPECTRA,
    NOISE,
    PREVIOUS,
    REFERENCE,
    digest_intervals,
    execute_main,
    fasta,
    location,
    module,
    run,
    save,
    sha,
    table,
)


def choose_mixture(unit):
    from pyteomics.mass import calculate_mass

    rng = random.Random(unit["seed"])
    proteins = dict(fasta(REFERENCE))
    valid = sorted(
        a
        for a, seq in proteins.items()
        if len(seq) >= 150 and re.fullmatch("[ACDEFGHIKLMNPQRSTVWY]+", seq)
    )
    omitted_n = 64 if unit["regime"] == "baseline" else 128
    chosen = rng.sample(valid, 256 + omitted_n)
    retained, omitted = chosen[:256], set(chosen[256:])
    search = {a: seq for a, seq in proteins.items() if a not in omitted}
    digest = {
        seq[b:e].replace("I", "L") for seq in search.values() for b, e in digest_intervals(seq)
    }
    n_present = 1536 if unit["regime"] == "baseline" else 1024
    rows, used = [], set()
    for kind, accessions, count in [
        ("present", retained, n_present),
        ("omitted", sorted(omitted), N_SPECTRA - n_present),
    ]:
        candidates = {}
        for accession in accessions:
            seq = proteins[accession]
            for begin, end in digest_intervals(seq):
                peptide = seq[begin:end]
                canonical = peptide.replace("I", "L")
                if (
                    not 7 <= len(peptide) <= 30
                    or canonical in used
                    or (kind == "omitted" and canonical in digest)
                ):
                    continue
                charges = [
                    z
                    for z in (2, 3)
                    if 400
                    <= calculate_mass(sequence=peptide, charge=z)
                    + peptide.count("C") * 57.021464 / z
                    <= 1000
                ]
                if charges:
                    candidates.setdefault(canonical, (peptide, accession, charges))
        if len(candidates) < count:
            raise ValueError(f"Insufficient {kind} candidates: {len(candidates)} < {count}")
        for canonical in rng.sample(sorted(candidates), count):
            peptide, accession, charges = candidates[canonical]
            used.add(canonical)
            sites = [str(i + 1) for i, aa in enumerate(peptide) if aa == "C"]
            rows.append(
                dict(
                    sequence=peptide,
                    generating_accession=accession,
                    source_kind=kind,
                    charge=rng.choice(charges),
                    mods=";".join(["Carbamidomethyl@C"] * len(sites)),
                    mod_sites=";".join(sites),
                    sid=len(rows),
                )
            )
    assert len(rows) == len(used) == N_SPECTRA
    assert all(
        (r["sequence"].replace("I", "L") not in digest)
        for r in rows
        if r["source_kind"] == "omitted"
    )
    return rows, search, omitted


def main(phase, index):
    import numpy as np
    import pandas as pd
    import torch
    from peptdeep.settings import global_settings
    from pyteomics import mzml as parser
    from pyteomics.mass import calculate_mass

    unit, output = location(phase, index)
    if phase == "validation":
        from common import ROOT

        gate = json.loads((ROOT / "runs/development/DEVELOPMENT_GATE.json").read_text())
        assert gate["ready_for_reserved_validation"], "Development gate is required"
    output.mkdir(parents=True, exist_ok=False)
    started = time.time()
    for sub in ("truth", "generator", "spectra"):
        (output / sub).mkdir()
    selected, reservoir, omitted = choose_mixture(unit)
    database = output / "generator/search_reservoir.fasta"
    with database.open("w") as handle:
        for accession, seq in reservoir.items():
            handle.write(f">{accession}\n{seq}\n")
    table(output / "truth/mixture.tsv", selected)
    global_settings["torch_device"]["device_type"] = "cpu"
    global_settings["thread_num"] = 4
    torch.set_num_threads(4)
    from peptdeep.pretrained_models import ModelManager

    manager = ModelManager(mask_modloss=True, device="cpu")
    frame = pd.DataFrame(selected)
    frame["nce"], frame["instrument"] = 30.0, "Lumos"
    predicted = manager.predict_all(
        frame,
        predict_items=["rt", "ms2"],
        frag_types=["b_z1", "b_z2", "y_z1", "y_z2"],
        multiprocessing=False,
    )
    precursors = predicted["precursor_df"].sort_values(["rt_pred", "sid"])
    mz, intensity = predicted["fragment_mz_df"], predicted["fragment_intensity_df"]
    noise = module("published_noise", NOISE)
    # Reuse the parser-tested mzML writer, without running the old generator.
    sys.path.insert(0, str(PREVIOUS))
    sys.path.insert(0, str(PREVIOUS.parent / "fusion_entrapment_2026-09-07"))
    writer = module("previous_mzml_writer", PREVIOUS / "pilot/generate.py")
    clean, truth = [], []
    for i, row in enumerate(precursors.to_dict("records")):
        begin, end = int(row["frag_start_idx"]), int(row["frag_stop_idx"])
        masses = mz.iloc[begin:end].values.ravel()
        values = intensity.iloc[begin:end].values.ravel()
        keep = (masses >= 200) & (masses <= 2000) & (values > 0.001)
        peaks = sorted(zip(masses[keep], values[keep] * 1e6))
        assert len(peaks) >= 4 and np.isfinite(values).all()
        expected = (
            calculate_mass(sequence=row["sequence"], charge=row["charge"])
            + row["sequence"].count("C") * 57.021464 / row["charge"]
        )
        assert abs(expected - row["precursor_mz"]) < 0.002, "Independent precursor mass check"
        clean.append(
            dict(
                params=dict(
                    TITLE=f"opaque_{i}",
                    SCANS=str(2 * i + 2),
                    PEPMASS=f"{row['precursor_mz']:.8f}",
                    CHARGE=f"{row['charge']}+",
                    RTINSECONDS=str((1 + i * 0.1) * 60),
                ),
                peaks=peaks,
            )
        )
        truth.append(
            dict(
                scannr=f"controllerType=0 controllerNumber=1 scan={2 * i + 2}",
                generating_sequence=row["sequence"],
                generating_accession=row["generating_accession"],
                source_kind=row["source_kind"],
                sid=row["sid"],
            )
        )
    clean_path = output / "generator/clean.mgf"
    noise.write_mgf(copy.deepcopy(clean), clean_path)
    stress = unit["regime"] == "stress"
    settings = dict(
        loss_rate=0.35 if stress else 0.25,
        jitter_ppm=12 if stress else 8,
        amp_jitter=0.35 if stress else 0.25,
        noise_low=80 if stress else 50,
        noise_high=15 if stress else 8,
        keep_top_k=12 if stress else 10,
        precursor_sigma_ppm=6 if stress else 5,
        precursor_tail_prob=0.4 if stress else 0.3,
        precursor_laplace_b=20 if stress else 15,
        isotope_shift_prob=0.12 if stress else 0.08,
        chimera_prob=0,
        mz_tol=0.02,
    )
    for acq_index, acquisition in enumerate(ACQUISITIONS):
        noisy = output / "generator" / f"{acquisition}_noisy.mgf"
        args = [
            sys.executable,
            NOISE,
            "-i",
            clean_path,
            "-o",
            noisy,
            "--seed",
            unit["seed"] * 100 + acq_index + 1,
        ]
        for key, value in settings.items():
            args.extend(["--" + key.replace("_", "-"), value])
        run(args, output / "generator" / f"{acquisition}_noise")
        perturbed = noise.parse_mgf(noisy)
        assert len(perturbed) == N_SPECTRA
        spectra = []
        for i, spectrum in enumerate(perturbed):
            assert spectrum["params"]["TITLE"] == f"opaque_{i}"
            spectra.append(
                dict(
                    pm=float(spectrum["params"]["PEPMASS"].split()[0]),
                    charge=int(spectrum["params"]["CHARGE"].rstrip("+")),
                    rt=float(spectrum["params"]["RTINSECONDS"]) / 60,
                    peaks=spectrum["peaks"],
                )
            )
        mzml = output / "spectra" / f"{acquisition}.mzML"
        writer.mzml(mzml, spectra)
        table(
            output / "truth" / f"{acquisition}.tsv",
            [dict(filename=str(mzml), **row) for row in truth],
        )
        with parser.MzML(str(mzml)) as reader:
            seen = {r["id"] for r in reader if r["ms level"] == 2 and len(r["m/z array"]) >= 15}
        assert seen == {r["scannr"] for r in truth}
    assert sha(output / "generator/donor_noisy.mgf") != sha(output / "generator/test_noisy.mgf")
    save(
        output / "GENERATED.json",
        dict(
            unit=unit,
            spectra_per_acquisition=N_SPECTRA,
            reservoir_proteins=len(reservoir),
            omitted_proteins=len(omitted),
            noise_settings=settings,
            matched_generating_peptides=True,
            conditionally_independent_noise=True,
            generator_versions={
                p: importlib.metadata.version(p) for p in ["peptdeep", "alphabase", "torch"]
            },
            elapsed_seconds=time.time() - started,
            slurm_job=os.getenv("SLURM_JOB_ID"),
            files={str(p.relative_to(output)): sha(p) for p in output.rglob("*") if p.is_file()},
        ),
    )
    print(unit, "generated with complete independent truth", flush=True)


if __name__ == "__main__":
    execute_main(main)
