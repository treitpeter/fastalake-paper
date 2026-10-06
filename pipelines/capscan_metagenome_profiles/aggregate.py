"""Collect completed DNA profiles; preserve unresolved taxonomic/functional mass."""
import argparse
import csv
import hashlib
import json
import re
import subprocess
import sys
import shutil
from run_selection import selected_runs, provenance, require_input_gate
from collections import defaultdict
from pathlib import Path
import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parent
OLD = BASE.parent / 'enterotypes_2026-09-07'


def function_terms(row):
    domains = [x for x in str(row.get('eggNOG_OGs', '')).split(',')
               if re.fullmatch(r'[^,]+@(?:2\|Bacteria|2157\|Archaea)', x)]
    if len(domains) != 1:
        return {'OG': [], 'KO': [], 'Pathway': []}
    return {'OG': domains,
            'KO': sorted(set(re.findall(r'K\d{5}', str(row.get('KEGG_ko', ''))))),
            'Pathway': sorted(x for x in set(re.findall(r'map\d{5}', str(row.get('KEGG_Pathway', ''))))
                              if not x.startswith('map011'))}


def read_annotations(path):
    with path.open() as f:
        for line in f:
            if line.startswith('#query\t'):
                reader = csv.DictReader(f, fieldnames=line[1:].rstrip('\n').split('\t'), delimiter='\t')
                out = {}
                for r in reader:
                    if r['query'].startswith('#'):
                        continue
                    if r['query'] in out:
                        raise ValueError('Duplicate annotation gene')
                    out[r['query']] = r
                return out
    raise ValueError('Missing eggNOG annotation header')


def taxonomy(path, mapping):
    genera = defaultdict(float)
    species = {}
    total = unresolved = bacterial = archaeal = 0.
    with path.open() as f:
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            fields = line.rstrip().split('\t')
            value = float(fields[-1])
            if not np.isfinite(value) or value < 0:
                raise ValueError('Invalid taxonomic abundance')
            match = re.search(r'(?:ref|meta|ext)_mOTU_v\d+_\d+', fields[0])
            key = match.group() if match else fields[0]
            if key in species:
                raise ValueError('Duplicate mOTU')
            species[key] = value
            total += value
            tax = mapping.get(key, [])
            if tax and tax[0] == 'd__Bacteria':
                bacterial += value
            if tax and tax[0] == 'd__Archaea':
                archaeal += value
            # Use only an explicit genus rank, never a genus embedded in an ambiguity label.
            if len(tax) >= 6 and tax[5].startswith('g__') and len(tax[5]) > 3:
                genera[tax[5][3:]] += value
            else:
                unresolved += value
    if not np.isclose(sum(genera.values()) + unresolved, total):
        raise ValueError('Taxonomic mass accounting failed')
    return dict(genera), species, dict(total_profile_mass=total,
            unresolved_genus_mass=unresolved, bacterial_mass=bacterial, archaeal_mass=archaeal)


def gene_profiles(root):
    q = pd.read_csv(root / 'salmon_quant/quant.sf', sep='\t').set_index('Name')
    if not q.index.is_unique or not np.isfinite(q.to_numpy()).all() or (q < 0).any().any():
        raise ValueError('Invalid gene quantification table')
    ann = read_annotations(root / 'genes.emapper.annotations')
    if set(ann) - set(q.index):
        raise ValueError('Annotation genes absent from Salmon reference')
    profiles = {level: defaultdict(float) for level in ['OG', 'KO', 'Pathway']}
    mapped = {level: 0. for level in profiles}
    for gene, value in q.TPM.items():
        terms = function_terms(ann.get(gene, {}))
        for level, labels in terms.items():
            if labels:
                mapped[level] += value
                for label in labels:
                    profiles[level][label] += value / len(labels)
    for level in profiles:
        if not np.isclose(sum(profiles[level].values()), mapped[level]):
            raise ValueError('Functional mass accounting failed')
    combined = q.join(pd.DataFrame.from_dict(ann, orient='index').drop(columns=['query'], errors='ignore'))
    combined.to_csv(root / 'gene_abundance.tsv.gz', sep='\t')
    meta = json.loads((root / 'salmon_quant/aux_info/meta_info.json').read_text())
    audit = dict(genes=len(q), annotated_genes=len(ann), total_gene_TPM=float(q.TPM.sum()),
                 assigned_fragments=float(q.NumReads.sum()),
                 mapped_TPM_by_level=mapped, percent_mapped=meta.get('percent_mapped'),
                 num_processed=meta.get('num_processed'), num_mapped=meta.get('num_mapped'))
    return {k: dict(v) for k, v in profiles.items()}, audit


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--allow-partial', action='store_true', help='Smoke testing only; list all missing run/stage outputs')
    p.add_argument('--output', type=Path, default=BASE / 'aggregate')
    a = p.parse_args()
    runs = pd.DataFrame(selected_runs(BASE)).set_index('run_accession')
    if (BASE / 'manifests/selection.json').exists() and not a.allow_partial:
        require_input_gate(BASE)
    taxmap = {}
    with (BASE / 'mOTUs_3.0.0_GTDB_tax.tsv').open() as f:
        for line in f:
            fields = line.rstrip('\n').split('\t'); taxmap[fields[0]] = fields[1:]
    matrices = {k: {} for k in ['genus', 'mOTU', 'OG', 'KO', 'Pathway']}
    missing, audit = [], []
    for run in runs.index:
        root = BASE / 'results' / run
        item = dict(run_accession=run, specimen=runs.loc[run, 'specimen'])
        if (root / '02_motus.done.json').exists():
            g, s, t = taxonomy(root / 'motus.tsv', taxmap)
            matrices['genus'][run] = g; matrices['mOTU'][run] = s
            item.update(t)
        else:
            missing.append(dict(run_accession=run, stage='taxonomy'))
        if (root / 'genes.complete.json').exists() or (root / 'all.complete.json').exists():
            f, coverage = gene_profiles(root)
            for level in f:
                matrices[level][run] = f[level]
            item.update(coverage)
        else:
            missing.append(dict(run_accession=run, stage='genes'))
        if len(item) > 2:
            audit.append(item)
    a.output.mkdir(parents=True, exist_ok=True)
    for name in ['selection.json', 'active_runs.tsv', 'excluded_runs.tsv']:
        source = BASE / 'manifests' / name
        if source.exists():
            shutil.copy2(source, a.output / name)
    pd.DataFrame(missing, columns=['run_accession', 'stage']).to_csv(a.output / 'missing_profiles.tsv', sep='\t', index=False)
    if missing and not a.allow_partial:
        raise RuntimeError('Incomplete profiles; see missing_profiles.tsv. No combined results produced.')
    for level, values in matrices.items():
        if not values:
            continue
        table = pd.DataFrame.from_dict(values, orient='index').fillna(0.)
        table.index.name = 'run_accession'
        table.to_csv(a.output / ('metaG_%s.tsv' % level), sep='\t')
    pd.DataFrame(audit).to_csv(a.output / 'coverage.tsv', sep='\t', index=False)
    if matrices['genus']:
        subprocess.run([sys.executable, str(OLD / 'project_genus_reference.py'),
                        '--input', str(a.output / 'metaG_genus.tsv'), '--output', str(a.output / 'metaG_reference_FKM.tsv'),
                        '--assay', 'metagenomic_genus'], check=True)
    status = dict(complete=not missing, missing_run_stages=len(missing),
                  selection=provenance(BASE),
                  rows_by_level={k: len(v) for k, v in matrices.items()},
                  semantics='Relative DNA gene abundance within assembled gene references; not RNA expression, absolute copies, or pathway flux.',
                  taxonomy='mOTUs 3.0.3; official GTDB R207 mapping; exact genus names retained, including suffixes.',
                  reference='Released FKM projection with mapped mass reported; not the online XGBoost model or a capsule diagnostic.',
                  manifest_sha256=hashlib.sha256((BASE / 'manifests/runs.tsv').read_bytes()).hexdigest())
    (a.output / 'aggregation_status.json').write_text(json.dumps(status, indent=2))
    print(json.dumps(status, indent=2))


if __name__ == '__main__':
    main()
