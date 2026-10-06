# Data availability (proposed statement)

Proposed text for the manuscript. Placeholders are marked; no accession or DOI has been invented. Public accessions below are those already in the v28d Data availability section (line 301) and are transcribed, not validated here.

> The healthy-stool mass spectrometry data have been deposited to the ProteomeXchange Consortium via the PRIDE partner repository with the dataset identifier [PXD accession to be added]. The PREDICT mass spectrometry data have been deposited under [PXD accession to be added] [access conditions to be added]. For five PREDICT acquisitions (MicrobPred_13, 44, 90, 115 and 181) the deposited mzML files are the validated reconversions used in the reported analyses. The CAMPI and CapScan mass spectrometry data reanalysed here are available under PXD023217 and PXD038906. CAMPI assemblies and metagenomic and metatranscriptomic reads are available from the European Nucleotide Archive under PRJEB42466. CapScan 16S rRNA and shotgun metagenomic reads are available from NCBI under BioProject PRJNA822660, and CapScan metabolomics data from the Metabolomics Workbench under ST002073, ST002075, ST002407, ST002409 and ST002411. [PREDICT metagenome, metatranscriptome and metabolite access route to be added.] AlphaNovo prediction exports for all 726 acquisitions and the selected search FASTA files are available at [accession to be added]. Source data are provided with this paper (Zenodo, DOI [to be added]; file list in source_data/ZENODO_MANIFEST.tsv); they underlie every figure panel, table and numerical statement and include complete test families, non-significant results and sensitivity analyses. Source data tables in the paper repository (source_data/) are released under CC BY 4.0. PREDICT participant and specimen codes are replaced by random pseudonyms used consistently across all files. Patient-level clinical covariates of PREDICT are not public: they are available to qualified researchers on request to the corresponding authors, subject to the ethics approval and data protection terms of the PREDICT study [route, committee, approval number and expected response time to be decided by the authors].

## What is public and what is controlled

| item | status |
|---|---|
| source tables behind every figure, table and number | in `source_data/` and the Zenodo package (CC BY 4.0, pending author confirmation) |
| PREDICT pseudonyms (`PREDICT_P###`, `PREDICT_A###`) | public; the mapping to raw identifiers is controlled |
| PREDICT clinical covariates (age, sex, MELD, CRP, antibiotic exposure, visit dates, tube identifiers) | controlled; listed in `controlled/CONTROLLED_INPUTS.md` with file, columns and sha256 |
| PREDICT study visit codes in Source Data sheets | listed as controlled in the unit `OPEN_ITEMS.md`; blanked in the repository copies |
| raw mass spectrometry files | PRIDE (accessions above, placeholders for healthy-stool and PREDICT) |
| files above 20 MB not in the repository | listed with sha256 in `source_data/ZENODO_MANIFEST.tsv` and in unit entries with prefix `ZENODO` |

## Open

- All bracketed placeholders (PXD accessions, AlphaNovo export and FASTA accession, Zenodo DOI).
- The access route for controlled PREDICT data is left to the authors; the repository only documents what is controlled and how code behaves without it.
- `source_data/ZENODO_MANIFEST.tsv` is empty and is to be filled from the unit `ZENODO` entries at release.
