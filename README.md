# spruce-embryogenesis-miRNA

Analysis code for:

**miRNA dynamics in Norway spruce embryogenesis**
Street, N.R., Rizvi, A.Z., Kochakarn, T., Lind Björs, J., Stojkovič, K., Canovi, C. & Egertsdotter, U.
*Physiologia Plantarum* (DOI assigned on publication).

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24307186.svg)](https://doi.org/10.5281/zenodo.24307186)

This repository contains the R scripts that reproduce the analyses, tables and figures of the article.
The input data and the processed tables are archived separately (see **Data**); download them into
the repository folder before running the scripts.

## Overview
MicroRNAs (miRNAs) were profiled by small RNA sequencing across zygotic embryogenesis (ZE) in Norway
spruce (*Picea abies*), with zygotic embryos (Zem) and female gametophytes (FG) dissected from the same
seeds, and across somatic embryogenesis (SE). The scripts test differential expression between adjacent
stages for miRNAs and mRNAs, map the somatic stages onto the zygotic series, test the expression support
for predicted miRNA-target pairs, and run Gene Ontology (GO) enrichment of the supported targets.

## Data
- **Processed data and results**: SciLifeLab Data Repository, https://doi.org/10.17044/scilifelab.35037526.
  Download `data.tar.gz` and unpack it in the repository folder:

      tar -xzf data.tar.gz

  This creates `data/` (counts, library metadata, expression matrices, miRNA catalogue, predicted
  target pairs, GO annotation) and `results/` (the tables written by the scripts, and the processed
  tables that the figure scripts read). `FILE_DESCRIPTIONS.txt` in the archive describes every file.
- **Raw sequencing reads**: European Nucleotide Archive (ENA), accession PRJEB72619.

## Processing of the raw reads
The folder `upstream/` holds the scripts that processed the raw reads into the deposited inputs
(read trimming, miRNA prediction and annotation, the replicate check, counting, target prediction and
mRNA quantification). They were run on a computing cluster and are provided as run; see
`upstream/README.md`.

## Software
All scripts were run with R 4.4.3 and Bioconductor 3.20. `environment.yml` lists the exact package
versions. To create the environment with conda or mamba:

    conda env create -f environment.yml
    conda activate spruce-embryogenesis-mirna

## Running the scripts
Run each script with `Rscript` from the repository folder; all paths are relative to it.
`REPRODUCE.md` lists the script behind every figure and table, and the order in which to run them.
The analysis scripts overwrite the corresponding tables in `results/`, and the figure scripts write a
PDF and a PNG to `figures/`.

## Citation
Please cite the article, the data item (https://doi.org/10.17044/scilifelab.35037526) and, for the
code, version 1.0.0 archived at Zenodo (https://doi.org/10.5281/zenodo.24307186). `CITATION.cff` gives the
citation in machine-readable form.

## Licence
The code is released under the MIT licence (see `LICENSE`). The data item is licensed under CC BY 4.0.
