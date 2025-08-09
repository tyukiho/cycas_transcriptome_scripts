# RNA-seq analysis scripts for "Intermediate Evolutionary State of Motile Sperm and Pollen Tubes in the Extant Gymnosperm Cycas revoluta"

This repository contains R scripts used for the analysis in the study:

> Yukiho Toyama, Satohiro Okuda, Takamasa Suzuki, Tetsuya Higashiyama. "Intermediate Evolutionary State of Motile Sperm and Pollen Tubes in the Extant Gymnosperm Cycas revoluta." 
> in preparation

## Contents
- `script1.R`: Ortholog Group enrichment analysis using Fisher's exact test
- `script2.R`: Semantic similarity analysis and 2D visualization of GO terms
- `script3.R`: Upset plot visualization

## Requirements
- R version: 4.2.0
- Required R packages:  
  - GO.db (v3.16.0) 
  - AnnotationDbi (v1.60.2)
  - ggplot2 (v3.5.2)  
  - ggrepel (v0.9.6)
  - dplyr (v1.1.4)
  - tidyr (v1.3.1)
  - UpSetR (v1.4.0)
  - tibble (v3.3.0)

## How to run
Place input files in the `data/` directory.
 - Example: `data/counts_matrix.tsv`, `data/sample_metadata.tsv`

## Zenodo
https://doi.org/10.5281/zenodo.16784754
