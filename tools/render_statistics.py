#!/usr/bin/env python3
"""Build the statistics index from the unit records, without duplicating results."""
import argparse
import csv
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]


def render():
    lines = [
        '# Statistics', '',
        'Statistical records are kept alongside the analysis that produces each result. '
        'Use the links below for the test, statistical unit, sample size, repeated-measure '
        'structure, correction family, estimate, uncertainty, seed and implementation. '
        'Descriptive summaries and historical comparisons are labelled in those records.', '',
        'Acquisitions are nested within specimens and participants. Check the stated '
        'sampling unit for each result: acquisition counts, participant comparisons and '
        'participant-block bootstrap intervals describe different quantities.', '',
        '## Updated selection tests and error controls', '',
        'The [updated-analysis tables](../source_data/v42_updates/README.md) contain '
        'peptide/precursor-group held-out testing, paired-entrapment controls, the current '
        'native-scan comparison and expanded validation. Both gains and losses are retained. '
        'The [change notes](../CHANGELOG.md) summarize the changes awaiting the next bioRxiv revision.', '',
        '```bash', 'python tools/check_updated_analyses.py --bootstrap', '```', '',
        'This repeats the group-fold and entrapment summary calculations and their '
        '10,000 participant-cluster bootstrap draws (seed 20261004). It also checks the '
        'current native-scan medians. Raw searches and nested model fits are upstream steps.', '',
        '## Records by analysis', '',
        '| Analysis | Statistical records | Notebook and source references |',
        '|---|---:|---|',
    ]
    for folder in ('figures', 'tables', 'text_statistics'):
        for path in sorted((REPO / folder).glob('*/STATS.tsv')):
            count = sum(1 for _ in csv.DictReader(path.open(), delimiter='\t'))
            rel = path.parent.relative_to(REPO)
            lines.append(f'| {path.parent.name} | [{count} records](../{rel}/STATS.tsv) | [Open analysis](../{rel}/) |')
    lines += ['', 'The [reviewer guide](REVIEWER_GUIDE.md) links registered manuscript values '
              'to their calculations. Each unit’s `OPEN_ITEMS.md` describes remaining gaps. '
              'Controlled-input calculations require the [documented access route](DATA_AVAILABILITY.md).', '',
              'Regenerate this index with `python tools/render_statistics.py`.', '']
    return '\n'.join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    path = REPO / 'docs/STATISTICS.md'
    expected = render()
    if args.check:
        assert path.read_text() == expected, 'Run python tools/render_statistics.py'
        print('Statistics index matches unit records')
    else:
        path.write_text(expected)


if __name__ == '__main__':
    main()
