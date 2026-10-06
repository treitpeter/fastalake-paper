"""Peptide-centric equal-split projection, as a sibling of the frozen A19 implementation.

Nothing under A19_4_6_biology_robustness/frozen/ is read for anything but import, and
nothing there is modified. This module supplies one replacement function, aggregate(),
with the frozen signature and return shape, and the loader the patched sibling derive
modules use for the projections the frozen code performs inline.

The frozen rule gave each study group's full intensity to every distinct annotation its
representative carried, so totals across annotations were not additive. The revised rule
allocates at the level of the accepted quantified LFQ feature: the feature's canonical
peptide fixes the set of compatible proteins in that cohort's study search space, the
distinct annotations of those proteins at the level form the target set, and the feature
intensity is divided equally among them. Intensity is therefore conserved, and features
whose peptide carries no qualifying annotation are retained explicitly as functionally
unassigned rather than dropped in silence.

Annotation vocabularies are unchanged: they are read from the frozen module's own
annotation_terms, which the matrix builder copied verbatim.
"""
from pathlib import Path
import json
import os
import numpy as np
import pandas as pd

PACKAGE = Path(__file__).resolve().parents[1]
INPUTS = Path(os.environ.get('FL_PC_MATRICES', str(PACKAGE / 'inputs' / 'matrices')))
# The frozen helpers the patched siblings load by path. They are read, never written.
FROZEN_A19 = PACKAGE.parents[1] / 'MM_review_analysis/A19_4_6_biology_robustness/frozen'
UNASSIGNED = 'FunctionallyUnassigned'

# Set by the driver before any frozen code runs.
COHORT = None
# Which grouping rule's study-group representative taxa label the species dimension.
# The function x sample allocation itself is grouping-invariant by construction: it depends
# only on the accepted feature, its canonical peptide and that peptide's compatible
# proteins. Only the species dimension can differ between the two rules.
RULE = 'simplified'

# Peptide-level analogues of the two group-level microbial rules in the frozen code.
#   microbial_single_og  <- common.aggregate(microbial_domain=True), which kept groups
#                           with exactly one bacterial/archaeal OG token.
#   microbial_broad      <- the inline rule in derive_capscan.py and derive_clinical.py,
#                           group.startswith(('MGYG','GMGC','GMSC')) or any OG token,
#                           excluding '_HUMAN'.
#   all                  <- no domain restriction.
VARIANTS = ('all', 'microbial_single_og', 'microbial_broad')

_cache = {}


def cohort():
    assert COHORT in ('HSTOOL', 'PREDICT', 'CAPSCAN'), 'set common_pepcentric.COHORT first'
    return COHORT


def matrix_path(level, variant, coh=None):
    assert level in ('KO', 'OG', 'Pathway'), level
    assert variant in VARIANTS, variant
    return INPUTS / f'{coh or cohort()}_{level}_{variant}.tsv.gz'


def profile(level, variant, coh=None):
    """function x sample allocation, as built. Rows include the remainder accounting row."""
    # The driver may repoint INPUTS at a different lineage; a cache keyed without the
    # directory would serve the previous lineage's matrix under the new selection.
    key = (str(INPUTS), coh or cohort(), level, variant)
    if key not in _cache:
        frame = pd.read_csv(matrix_path(level, variant, coh), sep='\t', index_col=0)
        assert frame.index.is_unique and frame.columns.is_unique
        assert np.isfinite(frame.to_numpy()).all() and (frame.to_numpy() >= 0).all()
        _cache[key] = frame
    return _cache[key]


def unassigned(level, variant, coh=None):
    """Per-sample intensity not allocated to this profile.

    The feature's canonical peptide carried no qualifying annotation on any protein eligible
    under this profile's rule. A feature annotated only on ineligible proteins lands here
    too, so this is an accounting remainder for this profile and not a measurement of
    missing biological annotation.
    """
    frame = profile(level, variant, coh)
    if UNASSIGNED not in frame.index:
        return pd.Series(0.0, index=frame.columns)
    return frame.loc[UNASSIGNED]


def functions(level, variant, coh=None):
    """function x sample allocation over annotations only, with the remainder row removed.

    The remainder is an accounting quantity, not a function: leaving it in would put a
    single very large pseudo-feature into every relative, Hellinger and distance transform
    downstream. It is reported separately by unassigned() and in the audit. Conservation is
    therefore allocated plus remainder, never allocated alone.
    """
    frame = profile(level, variant, coh)
    return frame.drop(index=[UNASSIGNED], errors='ignore')


def joint(coh=None, rule=None):
    """KO x species x sample allocation. Species is the representative taxon of the study
    group the given grouping rule gave the feature, so contributor labels stay
    representative-based and no LCA is introduced."""
    rule = rule or RULE
    assert rule in ('simplified', 'original'), rule
    stem = 'KO_species' if rule == 'simplified' else 'KO_species_original'
    path = INPUTS / f'{coh or cohort()}_{stem}.npz'
    with np.load(path, allow_pickle=False) as z:
        return dict(ko=z['ko'], species=z['species'], values=z['values'], samples=z['samples'])


def stats(coh=None):
    return json.loads((INPUTS / f'{(coh or cohort())}_variant_stats.json').read_text())


def aggregate(matrix, annotations, level, mode='split', microbial_domain=False):
    """Replacement for frozen common.aggregate, same signature and same return shape.

    Returns a sample x annotation frame, indexed by matrix.columns and sorted on columns,
    exactly as the frozen function did, so every frozen consumer reads it unchanged. The
    group matrix and the representative annotations are no longer the source of the
    projection; they are used only to fix the sample axis, which is asserted to agree.
    """
    assert mode == 'split', 'the revised framework is equal-split only; mode=%r' % mode
    variant = 'microbial_single_og' if microbial_domain else 'all'
    frame = functions(level, variant)
    extra = set(frame.columns) - set(matrix.columns)
    assert not extra, f'{level}/{variant}: samples absent from the group matrix: {sorted(extra)[:5]}'
    # The reverse direction matters more: a sample present in the group matrix but absent
    # here would be zero-filled and fed silently into the distance and Mantel transforms.
    # One CapScan acquisition is legitimately absent, the empty quantified one, so the
    # check names what it tolerates instead of allowing any number of them.
    missing = [s for s in matrix.columns if s not in set(frame.columns)]
    empty = [s for s in missing if float(matrix[s].sum()) == 0.0]
    assert missing == empty, (f'{level}/{variant}: acquisitions with signal are absent from '
                              f'the peptide-centric matrix: {sorted(set(missing) - set(empty))[:5]}')
    out = frame.T.reindex(index=matrix.columns).fillna(0.0)
    out.index.name = None
    return out.sort_index(axis=1)


def matrix_path_joint(coh=None, rule=None):
    """Path of the joint allocation actually used, for provenance records."""
    rule = rule or RULE
    stem = 'KO_species' if rule == 'simplified' else 'KO_species_original'
    return INPUTS / f'{coh or cohort()}_{stem}.npz'


def joint_long(coh=None, rule=None, drop_unassigned=True):
    """The joint allocation as the long sample/KO/representative_taxon/intensity frame the
    frozen healthy-stool leader code expects from KO_taxon_long.tsv.gz."""
    z = joint(coh, rule)
    ko, species, values, samples = z['ko'], z['species'], z['values'], z['samples']
    keep = np.flatnonzero(values.sum(axis=1) > 0)
    rows = []
    for i in keep:
        if drop_unassigned and ko[i] == UNASSIGNED:
            continue
        v = values[i]
        nz = np.flatnonzero(v > 0)
        if not len(nz):
            continue
        rows.append(pd.DataFrame(dict(sample=samples[nz], KO=ko[i],
                                      representative_taxon=species[i], intensity=v[nz])))
    out = pd.concat(rows, ignore_index=True)
    # 'Host' and 'No species annotation' are the frozen labels for unresolvable origin; the
    # leader code treats 'Unassigned' as unresolved, so restore that spelling for them.
    out['representative_taxon'] = out.representative_taxon.replace(
        {'No species annotation': 'Unassigned', 'Host': 'Unassigned'})
    return out


def hstool_local_candidate_ko():
    """Peptide-centric KO profiles for Figure 3f, one per local-candidate class.

    features x acquisitions, exactly the orientation the frozen code produced with
    aggregate(...).T. The class shares come from the published
    local_candidate_group_intensities.npz, so the classification itself is reused.
    """
    out = {}
    for cls in ['all', 'host', 'microbial']:
        frame = pd.read_csv(INPUTS / f'HSTOOL_KO_local_candidate_{cls}.tsv.gz',
                            sep='\t', index_col=0)
        out[cls] = frame.drop(index=[UNASSIGNED], errors='ignore')
    return out
