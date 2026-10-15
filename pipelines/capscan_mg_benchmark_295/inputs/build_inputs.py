#!/usr/bin/env python3
"""mg295 inputs: the 13 participant-15 acquisitions whose specimen has an active (completed) WGS run,
the 3 that stay out, and the control acquisition (an existing 282 pair).

Rule (identical to the 282): an acquisition is MG-matched iff its specimen (workbook Main) has one of
the 301 active runs (analyses/capscan_metag_analysis_2026-09-10/aggregate/active_runs.tsv).
Assembly = canonical/publication_data/sources/assemblies_prefixed/CAPSCAN/<SAMN>.psa_megahit.prodigal.faa
(the production chain: prefixed .faa -> rehash_metag_headers.py -> + human Swiss-Prot).

Writes ACQ13.tsv, EXCLUDED3.tsv, CONTROL.tsv. File sha256s are added by db/build_mg_fasta.py (Slurm).
"""
import csv, json, os, sys
from pathlib import Path

LAB = Path('LAB_ROOT')
R = LAB / 'release_prep/fastalake-paper_v28'
FINAL = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
ADDED = R / 'pipelines/capscan_participant_roster/participant15_added.tsv'
ACTIVE = LAB / 'analyses/capscan_metag_analysis_2026-09-10/aggregate/active_runs.tsv'
EXCL = LAB / 'analyses/capscan_metag_analysis_2026-09-10/aggregate/excluded_runs.tsv'
ENA = LAB / 'analyses/capscan_metag_availability_2026-09-09/ena_run_table_2026-09-10.tsv'
PREF = LAB / 'canonical/publication_data/sources/assemblies_prefixed/CAPSCAN'
HASHED = LAB / 'canonical/publication_data/sources/assemblies_hashed/CAPSCAN'
MAN347 = FINAL / 'public347/crossfit_capscan347/inputs/MANIFEST_347.json'
A08 = FINAL / 'public347/A08_347/results'
SS = FINAL / 'public347/spectrum_split347/outputs/CAPSCAN'
SE = FINAL / 'second_engine_all/inputs/manifest.json'
ROSTER = FINAL.parent / 'iteration3/disjoint/roster.tsv'
MGSOLO_MAN = LAB / 'canonical/publication/arms/manifests/metag_solo_CAPSCAN.tsv'
CONTROL = '20210513_EXPL0_PeTr_SA_T760_1311'   # roster i=335; also the public347 cross-fit worker control


def tsv(p):
    return list(csv.DictReader(open(p), delimiter='\t'))


def write(p, rows):
    with open(p, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        w.writeheader(); w.writerows(rows)


added = tsv(ADDED); assert len(added) == 16 and {r['subject_id'] for r in added} == {'15'}
active = {r['specimen']: r for r in tsv(ACTIVE)}; assert len(active) == 301
excluded = tsv(EXCL)
ena = tsv(ENA)
man = json.loads(MAN347.read_text()); tasks = {t['sample']: t for t in man['tasks'] if t['cohort'] == 'CAPSCAN'}
se = json.loads(SE.read_text())
se_fl = {t['sample']: t for t in se['tasks'] if t['cohort'] == 'CAPSCAN' and t['arm'] == 'FL_public'}
roster = tsv(ROSTER); assert len(roster) == 578
assert not {r['sample_id'] for r in added} & {r['sample'] for r in roster}

rows, out3 = [], []
for a in added:
    s, main = a['sample_id'], a['workbook_Main']
    r = active.get(main)
    ena_rows = [e for e in ena if e['sample_alias'].split('_')[0] == main]
    if r is None:
        samn = sorted({e['sample_accession'] for e in ena_rows})
        wgs = [e['run_accession'] for e in ena_rows if e['library_strategy'] == 'WGS']
        ex = [x["run_accession"] + ": " + x["reason"] for x in excluded if x["run_accession"] in wgs]
        asm = [str(PREF / f'{n}.psa_megahit.prodigal.faa') for n in samn if (PREF / f'{n}.psa_megahit.prodigal.faa').exists()]
        out3.append(dict(sample=s, region=a['region'], participant=a['subject_id'], workbook_Main=main,
                         biosamples=';'.join(samn), wgs_runs=';'.join(wgs) or 'none',
                         reason=('WGS run excluded from the 301 active runs: ' + ';'.join(json.dumps(x) for x in ex)) if wgs else
                                'no WGS run in PRJNA822660 (amplicon only); no assembly',
                         assembly_present=';'.join(asm) or 'none'))
        continue
    samn = r['sample_accession']
    faa = PREF / f'{samn}.psa_megahit.prodigal.faa'; assert faa.is_file(), faa
    hashed = HASHED / f'{samn}.psa_megahit.prodigal.fasta'; assert hashed.is_file(), hashed
    t = tasks[s]
    parity = json.loads((SS / s / 'PARITY.json').read_text())
    a08 = json.loads((A08 / s / 'COMPLETE.json').read_text()); assert a08['status'] == 'PASS'
    assert r['group'] == a['region'], (s, r['group'], a['region'])
    rows.append(dict(
        sample=s, region=a['region'], participant=a['subject_id'], workbook_Main=main, run=r['run_accession'],
        biosample=samn, assembly_faa=str(faa), production_hashed_fasta=str(hashed),
        mg_db=str(FINAL / f'mg295/db/_db/{s}.metag_hashed_plus_human.fasta'),
        mzml=t['mzml'], mzml_real=os.path.realpath(t['mzml']), csv=t['csv'], split_index=t['index'],
        odd_mzml_sha256=parity['output_sha256'], fl_disjoint_config=str(SS / s / 'HALF_DISJOINT/config.json'),
        mg_search_dir=str(FINAL / f'mg295/mg_full/out/{s}'),
        se_raw=se_fl[s]['raw']['path'], se_raw_sha256=se_fl[s]['raw']['sha256'],
        fl_a08_peptides_5_50=a08['counts']['peptide_q01'], fl_a08_peptides_9_50=a08['counts']['length9_50']))
assert len(rows) == 13 and len(out3) == 3, (len(rows), len(out3))
assert len({r['biosample'] for r in rows}) == 13
write(OUT / 'ACQ13.tsv', rows)
write(OUT / 'EXCLUDED3.tsv', out3)

# control: an existing 282 pair (roster row 335)
ro = next(r for r in roster if r['sample'] == CONTROL)
solo = next(r for r in tsv(MGSOLO_MAN) if r['sample'] == CONTROL)
samn = Path(solo['fasta']).name.split('.')[0]
write(OUT / 'CONTROL.tsv', [dict(sample=CONTROL, roster_i=ro['i'], region=ro['group'], participant=ro['subject_id'], biosample=samn,
                                 assembly_faa=solo['fasta'], production_hashed_fasta=str(HASHED / f'{samn}.psa_megahit.prodigal.fasta'),
                                 production_mg_db=ro['mg_fasta'], production_mg_config=ro['mg_frozen_config'],
                                 metag_solo_config=str(LAB / f'canonical/publication_data/arms/metag_solo/CAPSCAN/{CONTROL}/config.json'),
                                 mzml=ro['mg_mzml_frozen'], frozen_n_mg_5=ro['frozen_n_mg_5'], frozen_n_mg_9=ro['frozen_n_mg_9'],
                                 mg_expected_5_sha256=ro['mg_expected_5_sha256'], mg_expected_9_sha256=ro['mg_expected_9_sha256'])])
print('ACQ13', len(rows), 'EXCLUDED3', [(r['sample'], r['wgs_runs']) for r in out3])
