"""Read only the first MiB of each RAW; retain embedded acquisition filename."""
from pathlib import Path
from urllib.parse import unquote
import collections, csv, hashlib, json, re

O = Path(__file__).resolve().parent
rows = list(csv.DictReader((O / 'checksum_name_crosswalk.tsv').open(), delimiter='\t'))
evidence = []
for row in rows:
    with Path(row['mzml']).open() as f:
        xml = f.read(100000)
    match = re.search(r'<sourceFile [^>]*location="file:/*([^"]+)"', xml)
    raw = Path('/' + unquote(match[1])) if match else None
    row['raw_location_recorded_in_mzml'] = str(raw or '')
    if raw and not raw.is_file():
        adjacent = Path(row['mzml']).parent.parent / 'raw' / raw.name
        if adjacent.is_file():
            raw = adjacent
    row['raw'] = str(raw or '')
    row['header_original_names'] = ''
    row['header_specimen_number_changed'] = ''
    row['header_status'] = 'RAW_UNAVAILABLE'
    if raw and raw.is_file():
        with raw.open('rb') as f:
            data = f.read(1 << 20)
        strings = []
        for offset in (0, 1):
            text = data[offset:len(data)-((len(data)-offset)%2)].decode('utf-16le', errors='ignore')
            strings.extend(s for s in re.findall(r'[\x20-\x7e]{6,}', text)
                           if 'originalRawFileNameWithoutExtension' in s)
        originals = sorted(set(name for s in strings for name in
            re.findall(r'"originalRawFileNameWithoutExtension":"([^"]+)"', s)))
        evidence.append(dict(sample=row['sample'], raw=str(raw),
                             header_sha256=hashlib.sha256(data).hexdigest(), strings=strings))
        row['header_original_names'] = ';'.join(originals)
        row['header_status'] = 'PARSED' if originals else 'NO_ORIGINAL_NAME'
        number = re.search(r'T760_(\d+)', row['sample'])
        other = [re.search(r'T760_(\d+)', x) for x in originals]
        if number and other and all(other):
            row['header_specimen_number_changed'] = str(any(x[1] != number[1] for x in other))
with (O / 'all_RAW_header_names.tsv').open('w') as f:
    writer = csv.DictWriter(f, fieldnames=rows[0].keys(), delimiter='\t')
    writer.writeheader()
    writer.writerows(rows)
(O / 'all_RAW_header_evidence.json').write_text(json.dumps(evidence, indent=2))
summary = dict(status='HEADER_AUDIT_COMPLETE', n=len(rows),
               header_status=dict(collections.Counter(r['header_status'] for r in rows)),
               specimen_number_differences=[r for r in rows if r['header_specimen_number_changed']=='True'],
               scope='Embedded original acquisition names; not independent proof of biological sample identity.')
(O / 'HEADER_AUDIT.json').write_text(json.dumps(summary, indent=2))
print(summary['header_status'])
for r in summary['specimen_number_differences']:
    print(r['sample'], r['header_original_names'], r['current_331'], sep='\t')
