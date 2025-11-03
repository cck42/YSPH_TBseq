## TB-seq: *M tuberculosis* Amplicon Sequencing

Created by Chaney Kalinich, Freddy Gonzalez, and Seth Redmond (Grubaugh Lab)

This repo contains data, figures, and code related to the development and validation of a whole-genome amplicon sequencing panel for *M. tuberculosis*.

### Data

Metadata for all specimens sequenced is available in the 'metadata' folder. This includes sample type, Ct, concentration (in GE/uL), dilution from template (if applicable), primer scheme concentration, and sequencing date.

### Results

Most results can be found in the 'results' folder: 
- 'alignment_calling' : Summary statistics of reference-guided genome assembly for each specimen sequenced
- 'metagenomics_czid' : Results from submitting raw reads to the CZID mNGS Illumina pipeline to characterize microbial composition in clinical specimens
- 'pangenome' : Pangenome assembly and *in-silico* amplicon prediction results

Results from PrimalScheme (e.g., the primer schemes used) can be found in 'schemes.' 

Figures, and code to generate said figures, can be found in 'figures' or 'supplemental.'

### Pipelines

There are four data analysis / processing pipelines:

- 'alignment_calling' for data processing, primer clipping, variants calling and QC metrics
- 'make_pangenome' pangenome construction for primer design and specificity analyses
- 'clean_human' removing human reads for data archival
- 'nextstrain' nextstrain build
    
