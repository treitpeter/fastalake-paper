"""Freeze exact run/specimen/acquisition joins and analysis configuration."""
import hashlib
import json
from pathlib import Path
import pandas as pd

BASE = Path(__file__).resolve().parent
OLD = BASE.parent / 'enterotypes_2026-09-07'
DOWNLOAD = BASE.parent / 'capscan_metag_availability_2026-09-09'


def main():
    w = pd.read_csv(DOWNLOAD / 'ena_run_table_2026-09-10.tsv', sep='\t')
    w = w.loc[w.library_strategy.eq('WGS')].copy()
    alias = w.sample_alias.str.extract(r'^(?P<specimen>\d+)_(?P<donor>\d+)_S(?P<set>[^_]+)_T(?P<type>.+)_meta$')
    assert alias.notna().all().all(), 'Unparsed ENA aliases'
    w['specimen'] = alias.specimen.astype(int).astype(str)
    w['subject_id'] = alias.donor.astype(int)
    w['published_type'] = alias.type
    pub = pd.read_csv(OLD / 'CAPSCAN_published_shotgun_rows.tsv', sep='\t')
    retained = pub.loc[pub.drop_meta.eq(False)].copy()
    retained['specimen'] = retained.Main.astype(int).astype(str)
    assert w.specimen.is_unique and retained.specimen.is_unique
    assert set(w.specimen) == set(retained.specimen) and len(w) == 302
    w = w.merge(retained[['specimen', 'Subject', 'Type', 'meta_samplename', 'SeqRun',
                          'reads_meta_filtered', 'poor_seq_meta', 'drop_meta', 'duplicate_meta']],
                on='specimen', validate='one_to_one')
    assert w.subject_id.eq(w.Subject).all() and w.published_type.eq(w.Type).all()
    groups = {'Capsule 1': 'duodenum', 'Capsule 2': 'jejunum', 'Capsule 3': 'ileum',
              'Capsule 4': 'colon', 'Stool': 'stool', 'Saliva': 'saliva'}
    w['group'] = w.Type.map(groups)
    assert w.group.notna().all()
    files = pd.read_csv(DOWNLOAD / 'wgs_files.tsv', sep='\t')
    byrun = files.groupby('run_accession')
    w['task_index'] = w.run_accession.map(byrun.task_index.first())
    for mate in [1, 2]:
        f = files.loc[files.relative_path.str.endswith('_%d.fastq.gz' % mate)].set_index('run_accession')
        w['read%d' % mate] = w.run_accession.map(f.relative_path.map(lambda x: str(DOWNLOAD / x)))
        w['md5_%d' % mate] = w.run_accession.map(f.md5)
        w['bytes_%d' % mate] = w.run_accession.map(f.bytes)
    design = pd.read_csv(OLD / 'CAPSCAN_specimen_design.tsv', sep='\t', dtype={'specimen': str})
    assert design.sample_id.is_unique
    acquisitions = design.merge(w[['specimen', 'run_accession', 'sample_accession', 'subject_id', 'group']],
                                on='specimen', how='left', suffixes=('', '_ena'), validate='many_to_one')
    matched = acquisitions.run_accession.notna()
    assert acquisitions.loc[matched, 'subject_id'].eq(acquisitions.loc[matched, 'subject_id_ena']).all()
    assert acquisitions.loc[matched, 'group'].eq(acquisitions.loc[matched, 'group_ena']).all()
    acquisitions['status'] = matched.map({True: 'exact_specimen_match', False: 'no_public_retained_WGS'})
    w = w.sort_values('task_index')
    w.to_csv(BASE / 'manifests/runs.tsv', sep='\t', index=False)
    acquisitions.to_csv(BASE / 'manifests/acquisitions.tsv', sep='\t', index=False)
    acquisitions.loc[~matched].to_csv(BASE / 'manifests/excluded_acquisitions.tsv', sep='\t', index=False)
    summary = dict(runs=len(w), donors=int(w.subject_id.nunique()), compartments=w.group.value_counts().to_dict(),
                   all_runs_exactly_match_published_retained_samples=True,
                   proteomics_acquisitions=len(design), matched_acquisitions=int(matched.sum()),
                   matched_specimens=int(acquisitions.loc[matched, 'specimen'].nunique()),
                   excluded_acquisitions=int((~matched).sum()),
                   excluded_specimens=int(acquisitions.loc[~matched, 'specimen'].nunique()))
    inputs = [DOWNLOAD / 'ena_run_table_2026-09-10.tsv', DOWNLOAD / 'wgs_files.tsv',
              OLD / 'CAPSCAN_published_shotgun_rows.tsv', OLD / 'CAPSCAN_specimen_design.tsv',
              OLD / 'published_FKM_reference.npz', BASE / 'mOTUs_3.0.0_GTDB_tax.tsv']
    summary['input_sha256'] = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}
    (BASE / 'manifests/identity_audit.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
