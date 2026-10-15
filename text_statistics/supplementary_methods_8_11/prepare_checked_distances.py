"""Rebuild own-stool distances from the corrected Figure 5 OG profiles.

Run after placing the Figure 5 profile inputs in build/zenodo/Fig5. The
participant/region averaging and distances are identical to Fig5.stat_cd_profiles.
The original own_stool_distance_matrices.tsv remains a historical input.
"""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

from figures.Fig5 import panels

REPO = Path(__file__).resolve().parents[2]


def prepare():
    """Write the 675 distances after verifying the corrected participant summaries."""
    sd = REPO / "source_data/Fig5"
    profiles = REPO / "build/zenodo/Fig5"
    for name in ("OG_other_unresolved.npz", "OG_other_unresolved_axes.json"):
        assert hashlib.sha256((profiles / name).read_bytes()).hexdigest() == panels.ZENODO_FILES[name], name
    design = panels.load_design(sd)
    og = panels.frame(profiles, "OG_other_unresolved")
    assert set(design.index) == set(og.index)
    roster = design.index[og.loc[design.index].sum(1) > 0]
    counts = design.loc[roster].groupby("participant").region.nunique()
    participants = sorted(counts.index[counts.eq(5)])
    assert len(participants) == 15
    index = pd.MultiIndex.from_product([participants, panels.REG], names=["participant", "region"])
    d = design.loc[roster]
    values = og.loc[roster]
    relative = values.div(values.sum(1), axis=0)
    cells = relative.groupby([d.participant, d.region, d.specimen_id]).mean().groupby(level=[0, 1]).mean().loc[index]
    x = cells.to_numpy()
    reference = pd.read_csv(sd / "own_stool_participants.tsv", sep="\t", float_precision="round_trip")
    rows = []
    for prevalence in panels.PREVALENCES:
        kept = x[:, (x > 0).mean(0) >= prevalence]
        cube = np.sqrt(kept / kept.sum(1, keepdims=True)).reshape(15, 5, -1)
        distances = np.stack([cdist(cube[i, :4], cube[:, 4]).mean(0) for i in range(15)])
        own = np.diag(distances)
        other = (distances.sum(1) - own) / 14
        ref = reference[reference.prevalence.eq(prevalence)].set_index("participant").loc[participants]
        np.testing.assert_allclose(own, ref.own_stool_distance, rtol=1e-12, atol=0)
        np.testing.assert_allclose(other, ref.other_stool_mean_distance, rtol=1e-12, atol=0)
        for i, participant in enumerate(participants):
            for j, stool in enumerate(participants):
                rows.append(dict(prevalence=prevalence, capsule_participant=participant,
                                 stool_participant=stool, mean_distance=repr(float(distances[i, j]))))
    output = REPO / "source_data/supplementary_methods_8_11/own_stool_distance_matrices_checked.tsv"
    pd.DataFrame(rows).to_csv(output, sep="\t", index=False)
    receipt = {"rows": len(rows), "sha256": hashlib.sha256(output.read_bytes()).hexdigest(),
               "profile_inputs": {name: panels.ZENODO_FILES[name] for name in
                                  ("OG_other_unresolved.npz", "OG_other_unresolved_axes.json")},
               "reference": "source_data/Fig5/own_stool_participants.tsv",
               "reference_sha256": hashlib.sha256((sd / "own_stool_participants.tsv").read_bytes()).hexdigest()}
    output.with_suffix(".json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    prepare()
