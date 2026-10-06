"""SIBLING of release_prep/integration_2026-09-14/repo/manuscript/current/code/derive_hstool.py. The frozen file is not modified.

Reapply the five-design analysis and contributor eligibility to direct assignments."""
from collections import Counter, defaultdict
from itertools import combinations
import json
import math
import numpy as np
import pandas as pd
from common import ROOT, RUN, DATA, ARM, TIERS, load, write, annotation_terms, aggregate, record
import common_pepcentric as pc


def main():
    out = DATA / 'Figure3'
    out.mkdir(exist_ok=True)
    matrix, annot, summary = load('HSTOOL_noclust')
    design_path = ARM / '_figroot/publication/figures/Fig4/tables/design.tsv'
    design = pd.read_csv(design_path, sep='\t', dtype=str).fillna('')
    design = design[design.sample_id.isin(matrix.columns)].copy()
    assert len(design) == 58
    samples = {t: design.loc[design.group.eq(t), 'sample_id'].tolist() for t in TIERS}
    assert [len(samples[t]) for t in TIERS] == [12, 15, 10, 11, 10]
    write(design, out / 'design.tsv')
    q = pd.read_csv(DATA / 'HSTOOL_noclust/per_acquisition.tsv', sep='\t')
    q = q.merge(design[['sample_id', 'group']], left_on='sample', right_on='sample_id', validate='m:1')
    write(q, out / 'group_evidence_per_acquisition.tsv')
    write(q.groupby(['method', 'group']).median(numeric_only=True).reset_index(), out / 'group_evidence_by_tier.tsv')
    # Keep classifications for both reviewed reporting choices explicit. The
    # representative-only comparison is a named sensitivity until attribution is settled.
    host = matrix.index.str.contains('_HUMAN', regex=False)
    microbial = ~host & (matrix.index.str.startswith(('MGYG', 'GMGC', 'GMSC')) |
                        annot.OG.str.contains(r'@(?:2\|Bacteria|2157\|Archaea)', regex=True).to_numpy())
    local = np.load(DATA / 'HSTOOL_noclust/local_candidate_group_intensities.npz')
    frames = {'all': matrix}
    for key, field in [('host', 'Human'), ('microbial', 'Microbial')]:
        mm = pd.DataFrame(local[field], index=local['groups'], columns=local['samples']).reindex(index=matrix.index, columns=matrix.columns)
        frames[key] = mm.loc[mm.sum(axis=1) > 0]
    # Revised framework. The protein-group space is unchanged, being a protein-level
    # quantity. The KO space is the peptide-centric equal-split projection, with each
    # allocated share carried into the Human and Microbial fractions its own study group
    # holds in that acquisition, so the published local-candidate classification is reused
    # and not recomputed.
    ko_frames = pc.hstool_local_candidate_ko()
    pairs = []
    for cls, m in frames.items():
        ko = ko_frames[cls].reindex(columns=m.columns).fillna(0.0)
        for space, table in [('groups', m), ('KO', ko)]:
            for tier in TIERS:
                for a, b in combinations(samples[tier], 2):
                    v = table[[a, b]].to_numpy()
                    v = v[(v > 0).all(axis=1)]
                    corr = np.corrcoef(np.log2(v).T)[0, 1] if len(v) > 2 else np.nan
                    cv = np.median(v.std(axis=1, ddof=0) / v.mean(axis=1)) * 100 if len(v) else np.nan
                    pairs.append(dict(space=space, design=tier, cls=cls, i=a, j=b,
                                      r_per_pair=corr, cv_per_pair_pct=cv, n_features_per_pair=len(v)))
    pairs = pd.DataFrame(pairs)
    write(pairs, out / 'pyramid_pairwise_local_candidate.tsv')
    write(pairs.groupby(['space', 'design', 'cls']).median(numeric_only=True).reset_index(), out / 'pair_summary_local_candidate.tsv')

    # Revised framework. The KO x taxon table is the peptide-centric equal-split
    # allocation instead of the group-centric one. Species labels remain the representative
    # taxon of the study group the grouping rule gave the feature, so no taxonomic rule
    # changes. Acquisition totals are identical under either grouping rule.
    source = pc.matrix_path_joint()
    taxa = pc.joint_long()
    totals = matrix.sum()
    taxa['normalised'] = taxa.intensity / taxa['sample'].map(totals)
    tier_of = design.set_index('sample_id')['group'].to_dict()
    old = ROOT / 'FIGURATOR/Figure3/review_2026-09-07/functional'
    labels = pd.read_csv(old / 'per_function_acquisition.tsv', sep='\t')[['ko', 'label']].drop_duplicates('ko').set_index('ko').label.to_dict()
    records = []
    for (sample, ko), frame in taxa.groupby(['sample', 'KO']):
        resolved = frame[frame.representative_taxon.ne('Unassigned')].sort_values(['normalised', 'representative_taxon'], ascending=[False, True])
        total = frame.normalised.sum()
        rs = resolved.normalised.sum()
        top = resolved.iloc[0] if len(resolved) else None
        second = resolved.iloc[1].normalised if len(resolved) > 1 else 0
        records.append(dict(tier=tier_of[sample], sample=sample, ko=ko, label=labels.get(ko, ko),
                            normalised_function_intensity=total, resolved_fraction=rs/total if total else 0,
                            leading_species=top.representative_taxon if top is not None else '',
                            leader_share_of_resolved=top.normalised/rs if rs else 0,
                            leader_margin_over_runnerup=(top.normalised-second)/rs if rs else 0,
                            n_resolved_contributors=len(resolved)))
    write(pd.DataFrame(records), out / 'per_function_acquisition.tsv')
    lookup = {(r['sample'], r['ko']): r for r in records}
    all_kos = sorted({r['ko'] for r in records})
    keep = {t: {k for k in all_kos if sum(bool(lookup.get((s, k), {}).get('leading_species')) for s in samples[t]) >= math.ceil(.8 * len(samples[t]))} for t in TIERS}
    common = set.intersection(*(keep[t] for t in TIERS))
    def abundance(k):
        return np.mean([np.mean([lookup.get((s, k), {}).get('normalised_function_intensity', 0) for s in samples[t]]) for t in TIERS[2:]])
    # Preserve all twelve author-selected functions; order by the newly measured
    # abundance using the same equal weighting of biological designs.
    selected = pd.read_csv(old / 'proposed_top12.tsv', sep='\t').ko.tolist()
    selected.sort(key=lambda k: (-abundance(k), k))
    rows = []
    for tier in TIERS:
        for ko in sorted(keep[tier]):
            rr = [lookup[s, ko] for s in samples[tier] if lookup.get((s, ko), {}).get('leading_species')]
            counts = Counter(r['leading_species'] for r in rr)
            n = len(rr)
            rare = sum(1-math.comb(n-v, 8)/math.comb(n, 8) if n-v >= 8 else 1 for v in counts.values())
            rows.append(dict(tier=tier, ko=ko, label=labels.get(ko, ko), eligible_in_all_five_designs=int(ko in common),
                             n_observed_resolved=n, n_distinct_leaders=len(counts),
                             expected_distinct_leaders_in_8_resolved_acquisitions=rare,
                             mean_intensity_pct_equal_three_designs=100*abundance(ko),
                             species_names='; '.join(sorted(counts)), main_top12=int(ko in selected)))
    stats = pd.DataFrame(rows)
    write(stats, out / 'function_tier_summary.tsv')
    write(pd.DataFrame([dict(rank=i+1, ko=k, label=labels.get(k, k), mean_intensity_pct=100*abundance(k)) for i, k in enumerate(selected)]), out / 'proposed_top12.tsv')
    result = dict(n_acquisitions=58, n_common_functions=len(common), selected_functions=selected,
                  rarefied_medians=stats[stats.eligible_in_all_five_designs.eq(1)].groupby('tier').expected_distinct_leaders_in_8_resolved_acquisitions.median().to_dict(),
                  grouping_rule_for_species_labels=pc.RULE,
                  source_scope='Peptide-centric equal-split functional projection over the frozen direct peptide assignment; production representative KO and taxon annotation; unchanged acquisition designs and 80% eligibility / eight-observation rarefaction.')
    (out / 'SUMMARY.json').write_text(json.dumps(result, indent=2))
    (out / 'PROVENANCE.json').write_text(json.dumps([record(source), record(design_path), record(old / 'proposed_top12.tsv')], indent=2))
    print(json.dumps(result), flush=True)


if __name__ == '__main__':
    main()
