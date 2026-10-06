"""Prepare taxonomic and functional DNA profiles from one verified paired WGS run."""
import argparse
import csv
import datetime
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from run_selection import selected_runs, provenance, require_input_gate

BASE = Path(__file__).resolve().parent
DOWNLOAD = BASE.parent / 'capscan_metag_availability_2026-09-09'
MOTUS = BASE / 'vendor/mOTUs_3.0.3/motus/motus.py'
DB = MOTUS.parent / 'db_mOTU'
EGGDATA = Path('LAB_ROOT/.._Microbpred/Microbpred/eggnog-mapper/data')
EGGPY = os.path.expandvars('${FASTALAKE_EGGNOG_ENV}/bin/python')
EMAPPER = os.path.expandvars('${FASTALAKE_EGGNOG_ENV}/bin/emapper.py')


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def command(args, log):
    args = [str(x) for x in args]
    print(datetime.datetime.now(datetime.timezone.utc).isoformat(), json.dumps(args), flush=True)
    with log.open('a') as f:
        f.write(json.dumps(args) + '\n'); f.flush()
        subprocess.run(args, stdout=f, stderr=subprocess.STDOUT, check=True)


def stage(root, name, args, outputs):
    marker = root / (name + '.done.json')
    signature = hashlib.sha256(json.dumps([str(x) for x in args]).encode()).hexdigest()
    if marker.exists():
        old = json.loads(marker.read_text())
        if old['command_sha256'] != signature or any(not p.exists() for p in outputs):
            raise RuntimeError('Checkpoint mismatch for ' + name)
        for p in outputs:
            if sha(p) != old['output_sha256'][str(p)]:
                raise RuntimeError('Changed checkpoint output: ' + str(p))
        return
    if any(p.exists() for p in outputs):
        raise RuntimeError('Uncheckpointed stage output requires inspection: ' + name)
    command(args, root / (name + '.log'))
    if any(not p.exists() or not p.stat().st_size for p in outputs):
        raise RuntimeError('Missing or empty stage output: ' + name)
    marker.write_text(json.dumps(dict(command_sha256=signature,
                                     output_sha256={str(p): sha(p) for p in outputs}), indent=2))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--index', required=True, type=int)
    parser.add_argument('--stage', choices=['taxonomy', 'genes', 'all'], default='all')
    parser.add_argument('--threads', type=int, default=8)
    parser.add_argument('--assembly-memory-gb', type=int, default=24)
    parser.add_argument('--preflight', action='store_true')
    a = parser.parse_args()
    rows = selected_runs(BASE)
    chosen = [r for r in rows if int(r['task_index']) == a.index]
    if len(chosen) != 1:
        raise ValueError('Expected one exact run for index')
    row = chosen[0]
    root = BASE / 'results' / row['run_accession']
    programs = ['fastp', 'bwa', 'samtools']
    if a.stage in ['genes', 'all']:
        programs += ['megahit', 'prodigal', 'salmon', 'diamond']
    missing = [x for x in programs if not shutil.which(x)]
    if a.stage in ['taxonomy', 'all'] and not (DB / 'db_mOTU_versions').exists():
        missing.append(str(DB / 'db_mOTU_versions'))
    if a.stage in ['genes', 'all']:
        missing += [str(EGGDATA / x) for x in ['eggnog.db', 'eggnog_proteins.dmnd'] if not (EGGDATA / x).is_file()]
    marker = DOWNLOAD / 'status' / (row['run_accession'] + '.json')
    if not marker.exists() or json.loads(marker.read_text())['state'] != 'verified':
        missing.append('checksum-verified input pair')
    for mate in [1, 2]:
        path = Path(row['read%d' % mate])
        if not path.is_file() or path.stat().st_size != int(row['bytes_%d' % mate]):
            missing.append(str(path))
    report = dict(run=row['run_accession'], stage=a.stage, missing=missing,
                  tools={x: shutil.which(x) for x in programs})
    print(json.dumps(report, indent=2), flush=True)
    if missing:
        raise RuntimeError('Preflight incomplete')
    if a.preflight:
        return
    if (BASE / 'manifests/selection.json').exists():
        require_input_gate(BASE)
    if not os.environ.get('SLURM_JOB_ID'):
        raise RuntimeError('Run profiling within a bounded Slurm allocation')
    root.mkdir(parents=True, exist_ok=True)
    with (root / 'worker.lock').open('w') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        (root / 'input_provenance.json').write_text(json.dumps(dict(run=row,
                  manifest_sha256=sha(BASE / 'manifests/runs.tsv'), selection=provenance(BASE),
                  runner_sha256=sha(__file__), slurm_job_id=os.environ['SLURM_JOB_ID']), indent=2))
        r1, r2 = root / 'clean_1.fastq.gz', root / 'clean_2.fastq.gz'
        stage(root, '01_fastp', ['fastp', '-i', row['read1'], '-I', row['read2'],
              '-o', r1, '-O', r2, '--detect_adapter_for_pe', '--thread', a.threads,
              '--json', root / 'fastp.json', '--html', root / 'fastp.html'],
              [r1, r2, root / 'fastp.json', root / 'fastp.html'])
        if a.stage in ['taxonomy', 'all']:
            stage(root, '02_motus', [sys.executable, '-m', 'motus.motus', 'profile', '-f', r1, '-r', r2,
                  '-n', row['run_accession'], '-t', a.threads, '-db', DB,
                  '-M', root / 'motus_mgc.tsv', '-o', root / 'motus.tsv'],
                  [root / 'motus.tsv', root / 'motus_mgc.tsv'])
        if a.stage in ['genes', 'all']:
            assembly = root / 'assembly'
            stage(root, '03_megahit', ['megahit', '-1', r1, '-2', r2, '-t', a.threads,
                  '-m', a.assembly_memory_gb * 1024 ** 3, '-o', assembly], [assembly / 'final.contigs.fa'])
            stage(root, '04_prodigal', ['prodigal', '-i', assembly / 'final.contigs.fa',
                  '-p', 'meta', '-a', root / 'genes.faa', '-d', root / 'genes.fna',
                  '-o', root / 'genes.gff', '-f', 'gff'],
                  [root / 'genes.faa', root / 'genes.fna', root / 'genes.gff'])
            index = root / 'salmon_index'
            stage(root, '05_salmon_index', ['salmon', 'index', '-t', root / 'genes.fna',
                  '-i', index, '-p', a.threads, '--keepDuplicates'], [index / 'versionInfo.json'])
            stage(root, '06_salmon_quant', ['salmon', 'quant', '-i', index, '-l', 'IU',
                  '-1', r1, '-2', r2, '--meta', '--validateMappings', '-p', a.threads,
                  '-o', root / 'salmon_quant'], [root / 'salmon_quant/quant.sf', root / 'salmon_quant/aux_info/meta_info.json'])
            stage(root, '07_eggnog', [EGGPY, EMAPPER, '-i', root / 'genes.faa',
                  '--output', 'genes', '--output_dir', root, '--data_dir', EGGDATA,
                  '--cpu', a.threads, '-m', 'diamond', '--tax_scope', 'prokaryota_broad'],
                  [root / 'genes.emapper.annotations'])
        (root / (a.stage + '.complete.json')).write_text(json.dumps(dict(
            run=row['run_accession'], stage=a.stage, completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()), indent=2))


if __name__ == '__main__':
    main()
