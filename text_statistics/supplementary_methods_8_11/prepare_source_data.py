"""Write source_data/supplementary_methods_8_11/ from the original analysis outputs.

Run once by the authors on Slurm. The notebook never calls this script; it reads only
``source_data/``. The script documents exactly how each file in
``source_data/supplementary_methods_8_11/`` was made. Original locations are under
``LAB_ROOT``, the lab project directory, given by the environment variable
``FASTALAKE_LAB_ROOT``.

Files written (all tab-separated)
---------------------------------
own_stool_distance_matrices.tsv
    For each OG prevalence filter (0, 0.1, 0.2), the 15 x 15 matrix of mean Euclidean
    distances between the four capsule-region Hellinger profiles of one participant and
    the stool profile of every participant. Computed here with the unchanged code of the
    original producer (``capscan_fig5bcd.py``, functions ``frame`` and ``cells`` and the
    ``dd`` line) from its own inputs; the script asserts that the diagonal and row means
    equal the producer's ``own_stool_participants.tsv`` exactly.
capscan_dna_runs.tsv
    One row per CapScan metagenome run: participant, sampling group, whether the run is
    in the functional DNA aggregate and whether it is in the native gene-sequence profile.
predict_native_feature_summary.tsv
    PREDICT DNA and RNA rows of the native cross-omics feature summary (copied; CapScan
    rows are from the superseded 14-participant population and protein rows predate the
    revised functional allocation, so both are dropped).
eds_prediction_*.tsv
    Outputs of the exploratory held-out EDS prediction (Source Data package copy, which
    is pseudonymised), restricted to the columns the statistics need. The patient roster,
    which holds clinical covariates, is not copied.
nonspecific_search_summary.tsv
    One row per cohort: number of completed non-specific Sage searches and the search
    settings, each checked to be identical in every acquisition's configuration.
hstool_grouping_receipt.tsv
    The counts that the healthy-stool grouping receipts record (key, value).

Usage (on Slurm, about 2 min, 4 GB)::

    FASTALAKE_LAB_ROOT=... python text_statistics/supplementary_methods_8_11/prepare_source_data.py
"""

import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.spatial.distance import cdist

UNIT = Path(__file__).resolve().parent
REPO = UNIT.parents[1]
SD = REPO / "source_data" / UNIT.name
LAB = Path(os.environ["FASTALAKE_LAB_ROOT"])

MC = LAB / "review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01"
BR = MC / "biology_refit"
C2 = BR / "claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15"
PKG = BR / "claude_review_2026-10-03/final/figure_sd/source_data/package/FastaLake_Source_Data"
RECON = LAB / "review/v15_integration/jess_integration_2026-09-29/cohort_reconciliation"
METAG = LAB / "analyses/capscan_metag_analysis_2026-09-10"
A05 = LAB / "MM_review_analysis/A05_4_7_unspecific"
A13 = LAB / "MM_review_analysis/A13_4_7_predict_unspecific"
HSREC = (LAB / "00_PACKAGE_2026-09-17/cleanup_2026-09-22/v10_mm_philip_comment_audit_2026-09-22/nine_point_followup/"
         "manuscript_final_pass_2026-09-23/mm_v11_revision_2026-09-23/reappraisal_2026-09-24/supplements/"
         "revision_2026-09-24/delivery/source_data/compute_records/HSTOOL_grouping_COMPLETE.json")
HSKO = LAB / "review/v15_integration/jess_integration_2026-09-29/source_data/current/HSTOOL_KO_refit/upstream_receipts/COMPLETE.json"

REG = ['duodenum', 'jejunum', 'ileum', 'colon', 'stool']
INPUTS = {}


def sha(p):
    """Return the sha256 of file ``p`` and record it in INPUTS."""
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    INPUTS[str(Path(p).relative_to(LAB))] = h.hexdigest()
    return h.hexdigest()


def rd(p, **kw):
    sha(p)
    return pd.read_csv(p, sep="\t", **kw)


def wr(d, name):
    d.to_csv(SD / name, sep="\t", index=False, lineterminator="\n")


# ---------------------------------------------------------------- own-stool distances
def frame(folder, key):
    """Producer function (capscan_fig5bcd.py), unchanged."""
    sha(folder / (key + '_axes.json')); sha(folder / (key + '.npz'))
    a = json.loads((folder / (key + '_axes.json')).read_text()); v = sparse.load_npz(folder / (key + '.npz'))
    return pd.DataFrame(v.toarray(), index=a['samples'], columns=a['features'])


def own_stool():
    """Recompute the producer's participant x participant distance matrix ``dd``."""
    prof = C2 / "fig5/mc/profiles/CAPSCAN_noclust"
    design = rd(C2 / "roster/design_347.tsv").set_index('sample_id')
    og = frame(prof, 'OG_other_unresolved'); assert set(design.index) == set(og.index)
    roster = design.index[og.loc[design.index].sum(1) > 0]
    obs = design.loc[roster].groupby('participant').region.nunique()
    parts = sorted(obs.index[obs.eq(5)]); n = len(parts)
    ix = pd.MultiIndex.from_product([parts, REG], names=['participant', 'region'])
    d = design.loc[roster]

    def cells(v):
        rel = v.div(v.sum(1), axis=0)
        return rel.groupby([d.participant, d.region, d.specimen_id]).mean().groupby(level=[0, 1]).mean().loc[ix]

    x = cells(og.loc[roster]).to_numpy()
    ref = rd(C2 / "fig5/fig5bcd_347/own_stool_participants.tsv", float_precision="round_trip")
    rows = []
    for prev in [0., .1, .2]:
        keep = (x > 0).mean(0) >= prev; v = x[:, keep]; h = np.sqrt(v / v.sum(1, keepdims=True)); cube = h.reshape(n, 5, -1)
        dd = np.stack([cdist(cube[i, :4], cube[:, 4]).mean(0) for i in range(n)])
        own = np.diag(dd); other = (dd.sum(1) - own) / (n - 1)
        r = ref[ref.prevalence == prev].set_index('participant').loc[parts]
        assert np.array_equal(own, r.own_stool_distance.to_numpy()), prev
        assert np.array_equal(other, r.other_stool_mean_distance.to_numpy()), prev
        for i in range(n):
            for j in range(n):
                rows.append(dict(prevalence=prev, capsule_participant=parts[i], stool_participant=parts[j],
                                 mean_distance=repr(float(dd[i, j]))))
    wr(pd.DataFrame(rows), "own_stool_distance_matrices.tsv")


# ---------------------------------------------------------------- DNA runs
def dna_runs():
    runs = rd(METAG / "manifests/runs.tsv")
    agg = {}
    for lev in ['OG', 'KO', 'Pathway']:
        a = rd(METAG / f"aggregate/metaG_{lev}.tsv", index_col=0)
        agg[lev] = set(a.index[a.sum(1) > 0])
    assert agg['OG'] == agg['KO'] == agg['Pathway']
    native = RECON / "crossomics/native/current331_plus16_all/COMPLETE.json"; sha(native)
    nat = json.loads(native.read_text())
    native_ids = {r['id'] for r in nat['gene_profile']['input_receipts']}
    out = runs[['run_accession', 'subject_id', 'group']].rename(columns={'subject_id': 'participant', 'group': 'region'})
    out['functional_aggregate'] = out.run_accession.isin(agg['OG'])
    out['native_profile'] = out.run_accession.isin(native_ids)
    assert native_ids <= set(out.run_accession)
    assert sorted(out.loc[out.native_profile, 'participant'].unique().tolist()) == sorted(nat['participants'])
    wr(out.sort_values('run_accession'), "capscan_dna_runs.tsv")


# ---------------------------------------------------------------- PREDICT native features
def native_features():
    f = LAB / "review/v15_integration/jess_integration_2026-09-29/source_data/current/Crossomics/full/native_SOURCE_DATA_feature_summary.tsv"
    pk = PKG / "current/Crossomics/full/native_SOURCE_DATA_feature_summary.tsv"
    assert sha(f) == sha(pk)
    d = rd(f, float_precision='round_trip')
    # DNA and RNA rows only: the protein rows of this summary predate the revised functional
    # allocation (the protein layers are reported by Supplementary Fig. 11 instead).
    wr(d[(d.cohort == 'PREDICT') & d.layer.isin(['DNA', 'RNA'])], "predict_native_feature_summary.tsv")


# ---------------------------------------------------------------- nested EDS prediction
def eds_prediction():
    N = PKG / "Biology_refit/nested_refit"
    m = rd(N / "model_metrics.tsv", float_precision='round_trip'); wr(m, "eds_prediction_model_metrics.tsv")
    wr(rd(N / "permutation_results.tsv", float_precision='round_trip'), "eds_prediction_permutations.tsv")
    wr(rd(N / "increment_over_baseline.tsv", float_precision='round_trip'), "eds_prediction_increment.tsv")
    wr(rd(N / "selected_alphas.tsv", float_precision='round_trip'), "eds_prediction_selected_alphas.tsv")
    h = rd(N / "heldout_scores.tsv", float_precision='round_trip')
    blocks = [g.reset_index(drop=True) for _, g in h.groupby(['grouping', 'seed', 'model'], sort=False)]
    for g in blocks:   # every block lists the 85 patients in the same order
        assert g.idpatient.tolist() == blocks[0].idpatient.tolist()
        assert np.array_equal(g.ET_dysbiosis.to_numpy(), blocks[0].ET_dysbiosis.to_numpy())
    h['row'] = np.tile(np.arange(85), len(blocks))
    wr(h[['row', 'idpatient', 'ET_dysbiosis', 'grouping', 'seed', 'model', 'score']].rename(
        columns={'idpatient': 'patient_pseudonym', 'ET_dysbiosis': 'dysbiosis_score'}), "eds_prediction_heldout.tsv")
    fo = rd(N / "outer_folds.tsv")
    wr(fo[['seed', 'fold', 'row']], "eds_prediction_outer_folds.tsv")
    fa = rd(N / "fit_boundary_audit.tsv")
    fa['n_train'] = fa.train.str.count(';') + 1; fa['n_test'] = fa.test.str.count(';') + 1
    wr(fa[['seed', 'model', 'outer_fold', 'stage', 'n_train', 'n_test', 'features']], "eds_prediction_fit_audit.tsv")
    dr = rd(N / "distance_reference_folds.tsv", float_precision='round_trip')
    wr(dr[['grouping', 'seed', 'fold', 'training_EDS_cutoff', 'n_reference']], "eds_prediction_distance_reference.tsv")


# ---------------------------------------------------------------- non-specific searches
FIELDS = {
    'cleave_at': lambda c: c['database']['enzyme']['cleave_at'],
    'semi_enzymatic': lambda c: c['database']['enzyme']['semi_enzymatic'],
    'min_len': lambda c: c['database']['enzyme']['min_len'],
    'max_len': lambda c: c['database']['enzyme']['max_len'],
    'peptide_min_mass': lambda c: c['database']['peptide_min_mass'],
    'peptide_max_mass': lambda c: c['database']['peptide_max_mass'],
    'precursor_tol_ppm': lambda c: c['precursor_tol']['ppm'][1],
    'fragment_tol_ppm': lambda c: c['fragment_tol']['ppm'][1],
    'static_mod_C': lambda c: c['database']['static_mods'].get('C'),
    'variable_mod_M': lambda c: c['database']['variable_mods'].get('M', [None])[0],
    'max_variable_mods': lambda c: c['database']['max_variable_mods'],
    'lfq': lambda c: c['quant']['lfq'],
}


def nonspecific():
    rows = []
    for cohort, roots in [('CAPSCAN', [A05, C2 / "fig6/A05"]), ('PREDICT', [A13])]:
        values = {k: set() for k in FIELDS}; done = 0; tasks = 0
        for root in roots:
            for cfg in sorted((root / "inputs").glob("*__unspecific.json")):
                tasks += 1
                c = json.loads(cfg.read_text())
                for k, f in FIELDS.items():
                    values[k].add(json.dumps(f(c)))
                res = root / "results" / cfg.stem / "COMPLETE.json"
                if res.is_file() and json.loads(res.read_text()).get('status') == 'PASS':
                    done += 1
        row = dict(cohort=cohort, configured_searches=tasks, completed_searches=done)
        for k, v in values.items():
            assert len(v) == 1, (cohort, k, v)   # identical in every acquisition's configuration
            row[k] = json.loads(next(iter(v)))
        rows.append(row)
    wr(pd.DataFrame(rows), "nonspecific_search_summary.tsv")


# ---------------------------------------------------------------- healthy-stool grouping receipt
def hstool_receipt():
    g = json.loads(HSREC.read_text()); sha(HSREC)
    k = json.loads(HSKO.read_text()); sha(HSKO)
    det = {d['method']: d for d in g['detection']}
    rows = [('acquisitions', g['acquisitions']),
            ('quantified_protein_groups_direct_peptide_assignment', det['simplified']['detected_anywhere']),
            ('quantified_protein_groups_original_grouping', det['original']['detected_anywhere']),
            ('eligible_lfq_features', g['eligible_lfq_features']),
            ('ko_refit_quantified_features', k['quantified_features']),
            ('ko_refit_positive_acquisitions', k['positive_acquisitions'])]
    wr(pd.DataFrame(rows, columns=['key', 'value']), "hstool_grouping_receipt.tsv")


if __name__ == "__main__":
    SD.mkdir(parents=True, exist_ok=True)
    own_stool(); dna_runs(); native_features(); eds_prediction(); nonspecific(); hstool_receipt()
    out = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(SD.glob("*.tsv"))}
    print(json.dumps(dict(inputs=INPUTS, outputs=out), indent=1))
