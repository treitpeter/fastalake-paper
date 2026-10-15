# Statistics

Statistical records are kept alongside the analysis that produces each result. Use the links below for the test, statistical unit, sample size, repeated-measure structure, correction family, estimate, uncertainty, seed and implementation. Descriptive summaries and historical comparisons are labelled in those records.

Acquisitions are nested within specimens and participants. Check the stated sampling unit for each result: acquisition counts, participant comparisons and participant-block bootstrap intervals describe different quantities.

## Updated selection tests and error controls

The [updated-analysis tables](../source_data/v42_updates/README.md) contain peptide/precursor-group held-out testing, paired-entrapment controls, the current native-scan comparison and expanded validation. Both gains and losses are retained. The [change notes](../CHANGELOG.md) summarize the changes awaiting the next bioRxiv revision.

```bash
python tools/check_updated_analyses.py --bootstrap
```

This repeats the group-fold and entrapment summary calculations and their 10,000 participant-cluster bootstrap draws (seed 20261004). It also checks the current native-scan medians. Raw searches and nested model fits are upstream steps.

## Records by analysis

| Analysis | Statistical records | Notebook and source references |
|---|---:|---|
| Fig1 | [9 records](../figures/Fig1/STATS.tsv) | [Open analysis](../figures/Fig1/) |
| Fig2 | [4 records](../figures/Fig2/STATS.tsv) | [Open analysis](../figures/Fig2/) |
| Fig3 | [19 records](../figures/Fig3/STATS.tsv) | [Open analysis](../figures/Fig3/) |
| Fig4 | [51 records](../figures/Fig4/STATS.tsv) | [Open analysis](../figures/Fig4/) |
| Fig5 | [33 records](../figures/Fig5/STATS.tsv) | [Open analysis](../figures/Fig5/) |
| Fig6 | [112 records](../figures/Fig6/STATS.tsv) | [Open analysis](../figures/Fig6/) |
| SFig01 | [0 records](../figures/SFig01/STATS.tsv) | [Open analysis](../figures/SFig01/) |
| SFig02 | [7 records](../figures/SFig02/STATS.tsv) | [Open analysis](../figures/SFig02/) |
| SFig03 | [13 records](../figures/SFig03/STATS.tsv) | [Open analysis](../figures/SFig03/) |
| SFig04 | [8 records](../figures/SFig04/STATS.tsv) | [Open analysis](../figures/SFig04/) |
| SFig05 | [21 records](../figures/SFig05/STATS.tsv) | [Open analysis](../figures/SFig05/) |
| SFig06 | [8 records](../figures/SFig06/STATS.tsv) | [Open analysis](../figures/SFig06/) |
| SFig07 | [9 records](../figures/SFig07/STATS.tsv) | [Open analysis](../figures/SFig07/) |
| SFig08 | [0 records](../figures/SFig08/STATS.tsv) | [Open analysis](../figures/SFig08/) |
| SFig09 | [5 records](../figures/SFig09/STATS.tsv) | [Open analysis](../figures/SFig09/) |
| SFig10 | [14 records](../figures/SFig10/STATS.tsv) | [Open analysis](../figures/SFig10/) |
| SFig11 | [252 records](../figures/SFig11/STATS.tsv) | [Open analysis](../figures/SFig11/) |
| SFig12 | [3 records](../figures/SFig12/STATS.tsv) | [Open analysis](../figures/SFig12/) |
| SFig13 | [5 records](../figures/SFig13/STATS.tsv) | [Open analysis](../figures/SFig13/) |
| SFig14 | [0 records](../figures/SFig14/STATS.tsv) | [Open analysis](../figures/SFig14/) |
| SFig15 | [0 records](../figures/SFig15/STATS.tsv) | [Open analysis](../figures/SFig15/) |
| SData1 | [0 records](../tables/SData1/STATS.tsv) | [Open analysis](../tables/SData1/) |
| STab1 | [7 records](../tables/STab1/STATS.tsv) | [Open analysis](../tables/STab1/) |
| STab2 | [2 records](../tables/STab2/STATS.tsv) | [Open analysis](../tables/STab2/) |
| STab3 | [8 records](../tables/STab3/STATS.tsv) | [Open analysis](../tables/STab3/) |
| STab4 | [17 records](../tables/STab4/STATS.tsv) | [Open analysis](../tables/STab4/) |
| STab5 | [10 records](../tables/STab5/STATS.tsv) | [Open analysis](../tables/STab5/) |
| STab6 | [1 records](../tables/STab6/STATS.tsv) | [Open analysis](../tables/STab6/) |
| abstract_discussion | [6 records](../text_statistics/abstract_discussion/STATS.tsv) | [Open analysis](../text_statistics/abstract_discussion/) |
| methods | [42 records](../text_statistics/methods/STATS.tsv) | [Open analysis](../text_statistics/methods/) |
| supplementary_methods_1_5 | [92 records](../text_statistics/supplementary_methods_1_5/STATS.tsv) | [Open analysis](../text_statistics/supplementary_methods_1_5/) |
| supplementary_methods_6_7 | [43 records](../text_statistics/supplementary_methods_6_7/STATS.tsv) | [Open analysis](../text_statistics/supplementary_methods_6_7/) |
| supplementary_methods_8_11 | [238 records](../text_statistics/supplementary_methods_8_11/STATS.tsv) | [Open analysis](../text_statistics/supplementary_methods_8_11/) |

The [reviewer guide](REVIEWER_GUIDE.md) links registered manuscript values to their calculations. Each unit’s `OPEN_ITEMS.md` describes remaining gaps. Controlled-input calculations require the [documented access route](DATA_AVAILABILITY.md).

Regenerate this index with `python tools/render_statistics.py`.
