# Data availability

Public source datasets are linked below. File counts reflect the records
available on 12 October 2026.

| Data | Public record |
|---|---|
| CAMPI mass spectrometry | [PRIDE PXD023217](https://www.ebi.ac.uk/pride/archive/projects/PXD023217) |
| CapScan mass spectrometry | [PRIDE PXD038906](https://www.ebi.ac.uk/pride/archive/projects/PXD038906) |
| CAMPI nucleotide data | [ENA PRJEB42466](https://www.ebi.ac.uk/ena/browser/view/PRJEB42466); four public WGS runs were returned |
| CapScan nucleotide data | [BioProject PRJNA822660](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA822660); 302 WGS and 349 amplicon runs were returned by ENA |
| PREDICT nucleotide data | [ENA PRJEB89898](https://www.ebi.ac.uk/ena/browser/view/PRJEB89898); 390 public runs, with library strategies WGS (89) and OTHER (301) |
| CapScan metabolomics | Metabolomics Workbench [ST002073](https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&StudyID=ST002073), [ST002075](https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&StudyID=ST002075), [ST002407](https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&StudyID=ST002407), [ST002409](https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&StudyID=ST002409) and [ST002411](https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&StudyID=ST002411) |

The ENA results alone do not establish that all manuscript metagenomes and
metatranscriptomes have been deposited. In particular, PRJEB42466's returned run
list contains WGS entries only; the CAMPI transcriptomic file route needs reconciliation.
The PREDICT entries labelled OTHER must be matched to the study manifest before
being described as completed metatranscriptomic deposits.

## New study deposits

The healthy-stool mass spectrometry package is staged for PRIDE, including all
58 restored RAW files and their validated conversions. A completed PRIDE deposit
and accession are still required. AlphaNovo prediction exports for 726 acquisitions
and the selected search FASTA packages are staged; their archival record and DOI
are still required.

Source tables shipped in this repository are in `source_data/`. Larger inputs and
controlled files are described by the unit documentation and
[controlled input list](../controlled/CONTROLLED_INPUTS.md). The final Source Data
archive and its DOI are pending. The proposed data
licence is CC BY 4.0, subject to the authors' confirmation of redistribution rights.

## PREDICT access

PREDICT mass spectrometry, clinical covariates and metabolite data follow the
consortium's controlled data-hub route. Contact Jonel Trebicka (contact details in
the article) to initiate the required data transfer agreement and confirm access
conditions.

For five PREDICT acquisitions, the access package must distinguish both validated
conversion sets: ThermoRawFileParser 2.0.0.0 for searches and 1.3.4 for spectrum-split
and cross-fitting analyses. These versions must not be substituted for one another.
Public analysis tables use study pseudonyms; their correspondence to raw clinical
identifiers remains controlled.
