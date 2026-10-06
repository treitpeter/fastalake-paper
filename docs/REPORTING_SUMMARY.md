# Reporting summary (Nature Portfolio), filled from the manuscript and this repository

Headings are those of the Nature Portfolio Reporting Summary (life sciences form). Each answer is taken from the accepted manuscript (v28d) with the line quoted, or from the code in this repository (file named). Where neither answers the item the entry reads AUTHOR INPUT NEEDED. This is a working draft for the authors; it is not signed. The form asks that no field be answered "not applicable"; where an item does not apply, the authors must write the form text asked for in its help text.

## Statistics

For all statistical analyses, confirm that the following items are present in the figure legend, table legend, main text, or Methods section.

### The exact sample size (n) for each experimental group/condition, given as a discrete number and unit of measurement

Every row of `docs/STATISTICS.md` (from the unit `STATS.tsv` files) gives n and the statistical unit (acquisition, specimen, participant, patient, peptide, protein group, sequence, species, region, synthetic mixture). Cohort sizes in the text: 15 CAMPI, 58 healthy-stool, 306 PREDICT and 347 CapScan acquisitions (726 in all).

AUTHOR INPUT NEEDED: per-test n for feature-dependent and histidine strata that the text gives only as families; check the `n` column against the text.

### A statement on whether measurements were taken from distinct samples or whether the same sample was measured repeatedly

Column `distinct_or_repeated` of every `STATS.tsv` row states what is repeated within what (for example acquisitions within specimens, specimens within participants, three entrapment seeds per acquisition). The healthy-stool series is designed with 12 technical replicates, 15 independent preparations, 10 spatial aliquots, 11 longitudinal samples and 10 inter-individual measurements.

> "technical replicates" (v28d line 99)

### The statistical test(s) used AND whether they are one- or two-sided

Column `test` and `sided` of every `STATS.tsv` row. Permutation tests are two-sided with plus-one P values.

> "Tests used 1,999 permutations (seed 20260907), two-sided plus-one P values" (v28d line 815)

### A description of all covariates tested

Column `covariates` of `STATS.tsv`. Mixed models include participant or patient as the random level. Patient-level clinical covariates of PREDICT are controlled inputs (`controlled/CONTROLLED_INPUTS.md`) and are not in the repository.

AUTHOR INPUT NEEDED: a list of covariates for every model, where the text does not name them.

### A description of any assumptions or corrections, such as tests of normality and adjustment for multiple comparisons

Columns `correction`, `family` and `family_size` of `STATS.tsv`. Resampling keeps the stated pairing and hierarchy.

AUTHOR INPUT NEEDED: mixed-model assumption checks; Holm family size; sizes of the EDS and nine function-correlation families; correction for the two Wilcoxon endpoints (priorities 3 and 4 of the statistics audit, see `docs/STATISTICS.md`).

> "Bootstrap resampling uses the stated participant/patient or synthetic-mixture unit and preserves pairing." (v28d line 887)

> "No correction family was limited to displayed or significant comparisons." (v28d line 893)

### A full description of the statistical parameters including central tendency (e.g. means) or other basic estimates (e.g. regression coefficient) AND variation (e.g. standard deviation) or associated estimates of uncertainty (e.g. confidence intervals)

Columns `estimate`, `uncertainty` and `uncertainty_construction` of `STATS.tsv`.

AUTHOR INPUT NEEDED: the own-stool confidence interval construction, draw count and seed if absent from the text.

> "Boxes, median and interquartile range (IQR); whiskers, 1.5 IQR." (v28d line 127)

> "Means and 95% bootstrap intervals were calculated over whole mixtures" (v28d line 733)

### For null hypothesis testing, the test statistic (e.g. F, t, r) with confidence intervals, effect sizes, degrees of freedom and P value noted

Columns `statistic`, `df`, `exact_P`, `effect_size` and `effect_size_how` of `STATS.tsv`; inequalities such as q<0.05 are kept as printed.

AUTHOR INPUT NEEDED: mixed-model likelihood-ratio and Wald statistics with df; Wilcoxon statistics.

> "Region 4 exceeded stool by a smaller 4.7 percentage points" (v28d line 149)

### For Bayesian analysis, information on the choice of priors and Markov chain Monte Carlo settings

No Bayesian analysis is described in the manuscript.

AUTHOR INPUT NEEDED: confirm that none was used and enter the form text for that case.

### For hierarchical and complex designs, identification of the appropriate level for tests and full reporting of outcomes

Participants and patients are the resampling and permutation unit; see the design hierarchy in `docs/STATISTICS.md`.

> "All prevalence, population and distance-method results are retained." (v28d line 839)

### Estimates of effect sizes (e.g. Cohen's d, Pearson's r), indicating how they were calculated

Column `effect_size_how` of `STATS.tsv`.

> "For each acquisition, peptide gain was the percentage change in accepted peptides with FastaLake relative to the matched metagenome search." (v28d line 229)

Seeds and draws (not a form heading, but required by the statistics policy): columns `draws` and `seed` of `STATS.tsv`.

> "50,000 resamples of the 42 patients (seed 20260907)" (v28d line 795)

> "10,000 participant-level bootstrap resamples (seed 20261004)" (v28d line 777)

## Software and code

Policy information about availability of computer code

### Data collection

Q Exactive HF-X, HF and Plus and Orbitrap Exploris 480 instruments; conversion with ProteoWizard 3.0.24025 and ThermoRawFileParser (versions in the Methods).

AUTHOR INPUT NEEDED: acquisition-control software and version; healthy-stool RAW-derived settings; PREDICT collection and preprocessing, which the text marks as pending.

> "Exploris 480" (v28d line 205)

> "ProteoWizard" (v28d line 759)

### Data analysis

FastaLake (MIT, MannLabs/fasta-lake release named at submission); this paper repository (`environment/README.md`: CPython 3.12.4, 66 pinned packages); AlphaNovo checkpoint named in Methods; Sage v0.14.7 (binary reports 0.14.6) and a corrected build for the joint healthy-stool searches; eggNOG-mapper 2.1.12; MSFragger 4.1; Percolator 3.06.5; FDRBench 0.0.4; directLFQ 0.3.3. Package versions of each statistics function are in the `software_version` column of `STATS.tsv`.

AUTHOR INPUT NEEDED: release, commit and Zenodo DOI of FastaLake and of this repository; MMseqs2 and AlphaPeptDeep versions if inputs depend on them. AlphaNovo code is not included in the repository (its checkpoint and exports are inputs, see `docs/DATA_AVAILABILITY.md`).

> "Sage v0.14.7" (v28d line 233)

> "directLFQ" (v28d line 437)

For manuscripts utilizing custom algorithms or software that are central to the research but not yet described in published literature, software must be made available to reviewers. We strongly encourage code deposition in a community repository (e.g. GitHub).

Code access: this repository (reviewers receive it with the submission, see `docs/CODE_AVAILABILITY.md`).

## Data

Policy information about availability of data

Data availability statement: proposed text in `docs/DATA_AVAILABILITY.md`; the manuscript statement is at v28d line 301.

> "PXD023217" (v28d line 303)

AUTHOR INPUT NEEDED: the PXD accessions for healthy-stool and PREDICT proteomics, the accession for AlphaNovo exports and selected FASTAs, the PREDICT metagenome, metatranscriptome and metabolite access route, and the access terms (reason, conditions, contact, time frame) for controlled PREDICT data.

## Research involving human participants, their data, or biological material

Note that full information on the approval of the study protocol must also be provided in the manuscript.

### Reporting on sex and gender

Not stated in the accepted manuscript. AUTHOR INPUT NEEDED: sex and gender variables, how collected, distributions, and whether analysed.

### Reporting on race, ethnicity, or other socially relevant groupings

Not stated in the accepted manuscript. AUTHOR INPUT NEEDED: whether collected or analysed, categories, or the reason they are not reported.

### Population characteristics

PREDICT: patients with decompensated cirrhosis or acute-on-chronic liver failure (306 acquisitions, 77 complete patients in the clinical analyses). CapScan: 15 healthy participants. Healthy-stool variability series: repeated specimens and 10 different-donor measurements. Patient-level covariates are controlled inputs (`controlled/CONTROLLED_INPUTS.md`).

AUTHOR INPUT NEEDED: age and other characteristics, eligibility, total number of healthy volunteers (donor overlap unknown).

### Recruitment

Secondary analysis of existing studies; PREDICT stool is attributed to the MUCOSA-PREDICT cohort at Debrecen.

AUTHOR INPUT NEEDED: recruitment, inclusion and exclusion procedures and any selection bias of the reused cohorts.

### Ethics oversight

PREDICT: local ethics approval and written informed consent; CapScan: WIRB-Copernicus Group study and informed consent. CAMPI is a defined bacterial community.

AUTHOR INPUT NEEDED: full committee name and approval number for each cohort, and for the healthy-stool volunteers.

> "35281" (v28d line 211)

> "1186513" (v28d line 211)

## Field-specific reporting

Life sciences (the form requested). Behavioural and social sciences and Ecological, evolutionary and environmental sciences: not selected.

## Life sciences study design

All studies must disclose on these points even when the disclosure is negative.

### Sample size

All available acquisitions of the reused datasets were analysed; sizes are in the Results and in `figures/*/README.md`. No sample size calculation is described.

AUTHOR INPUT NEEDED: whether any calculation was made, and the rationale if not.

### Data exclusions

Exclusions are described in the text (39 differential-centrifugation acquisitions excluded from the variability series, incomplete profiles, complete-case populations, positive-detection eligibility, missing-PEP handling); the roster is Supplementary Data 1 (`tables/SData1/`). Exclusions are never limited to displayed results.

AUTHOR INPUT NEEDED: which criteria were fixed before the data were examined.

> "39 differential-centrifugation acquisitions" (v28d line 501)

### Replication

Replication structure is in the text; all 72 ordinary simulation outcomes and non-significant controls are retained.

AUTHOR INPUT NEEDED: an explicit statement on attempts to reproduce the principal findings.

> "All 72 ordinary outcomes were retained" (v28d line 733)

### Randomization

Computational randomisation (cross-fit folds, label permutations, entrapment and scramble controls) is seeded; seeds are in `STATS.tsv`.

AUTHOR INPUT NEEDED: whether acquisition order or sample allocation was randomised, or why not.

### Blinding

No statement in the manuscript.

AUTHOR INPUT NEEDED: whether investigators were blinded at any stage, or why not.

## Reporting for specific materials, systems and methods

Antibodies, eukaryotic cell lines, palaeontology and archaeology, animals and other organisms, plants, ChIP-seq, flow cytometry, MRI-based neuroimaging: no such material or method is described in the manuscript. AUTHOR INPUT NEEDED: tick the form boxes; the draft does not infer a negative from silence.

### Clinical data

Policy information about clinical studies

### Clinical trial registration

NCT03056612 is stated for PREDICT (v28d line 211). The present work is a secondary analysis.

AUTHOR INPUT NEEDED: confirm.

> "NCT03056612" (v28d line 211)

### Study protocol

Source studies are cited; no protocol location is given.

AUTHOR INPUT NEEDED: protocol location and deviations.

### Data collection

PREDICT site and cohort, instruments and sample preparation are described; collection dates and recruitment periods are not.

AUTHOR INPUT NEEDED: where, when and by whom clinical measurements were collected.

### Outcomes

The outcomes analysed here are proteomic identification outcomes, paired antibiotic-state differences and exploratory molecular associations.

AUTHOR INPUT NEEDED: separate any prespecified source-study outcome.

### Dual use research of concern

Hazards: no hazard is described in the manuscript. AUTHOR INPUT NEEDED: tick the form boxes.

## Disagreements with the manuscript-only draft

The draft `REPORTING_SUMMARY_DRAFT.md` (written from the manuscript alone, before this repository existed) is the base of this file. Where the repository answers an item, this file differs from it as follows.

1. Exact n and units: the draft marks per-test n as missing; every `STATS.tsv` row now has n and unit (1,033 rows, none empty in those two columns).
2. Software for statistical implementations: the draft lists software as missing for 124 register rows; every `STATS.tsv` row now names the implementation and version. The manuscript text still does not state them.
3. Seeds and draws: the draft lists seeds as missing for 36 rows; `STATS.tsv` gives seed and draws wherever a stochastic procedure was used (other rows are exact or deterministic). The text still lacks several.
4. Data analysis versions: the draft asks for FastaLake and study-analysis versions; the repository pins CPython 3.12.4 and 66 packages (`environment/`), but release, commit and DOI remain placeholders.
5. Source Data: the draft cannot certify the contents of the Source Data package; the repository checks every number against its source table (`docs/REVIEWER_GUIDE.md`, 1,331 numbers: 1,232 MATCH, 38 v28f changes, 61 OPEN).
6. Cohort size: the draft and the text use 331 CapScan acquisitions in places; the repository units use 347 (see `docs/MANUSCRIPT_CHANGES_FOR_V28F.md`).
7. Items the repository cannot answer and the draft also leaves open (sex and gender, ethics approval numbers, blinding, sample size rationale, accessions) stay AUTHOR INPUT NEEDED.

