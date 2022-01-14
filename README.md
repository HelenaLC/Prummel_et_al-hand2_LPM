This repository contains all code necessary to reproduce all analyses and figures from our manuscript:

> **Hand2 delineates mesothelium progenitors and is reactivated in mesothelioma**  
*Prummel KD, Crowell HL, Nieuwenhuize S, Brombacher EC, Daetwyler S, Soneson C, Kresoja-Rakic J, Kocere A, Ronner M, Ernst A, Labbaf Z, Clouthier DE, Firulli AB, Sánchez-Iranzo H, Naganathan SR, O’Rourke R, Raz E, Mercader N, Burger A, Felley-Bosco E, Huisken J, Robinson MD & Mosimann C*

A browesable `workflowr`<sup>[1](#f1)</sup> website of the analysis can be viewed [HERE](https://htmlpreview.github.io/?https://github.com/HelenaLC/Pummel_et_al-hand2_LPM/blob/master/docs/index.html).

## Prerequisites

The code in this repository was developed using **R v4.0.2** and **Bioconductor v3.11**. For installation of the required libraries, we'll fist install the `r BiocStyle::Biocpkg("BiocManager")` package:

```r
install.packages("BiocManager")
```

Versions of R and Bioconductor that are currently being run should be checked via:

```r
version
BiocManager::version()
```

Finally, the code chunk below will install all package dependencies:

```r
pkgs <- c("biomaRt", "circlize", "clustree", "ComplexHeatmap", 
    "cowplot", "dplyr", "ggplot2", "matrixStats", "RColorBrewer", 
    "readxl", "reshape2", "scater", "scran", "scales", 
    "SingleCellExperiment", "Seurat", "tidyr", "viridis")
BiocManager::install(pkgs, ask = FALSE)
```

## Contents

It follows a brief description of the directories and scripts contained herein. The contents of each `analysis/*.Rmd` script may be viewed on the landing page of the `workflowr` website (link above), but are still included here for completeness.

### Subdirectories

directory | description
----------|------------
data | raw and metadata required to run `analysis/*.Rmd` scripts
code | helper and utility scripts
analysis | .Rmd scripts to perform analyses and built the website
output | intermediate outputs (e.g. `.rds` objects) <br> generated by `analysis/*.Rmd` scripts
figures | stores all visual outputs of `analysis/*.Rmd` scripts (as `.pdf`)
docs | stores all `.html` (and other) outputs to generate the website
`_workflowr.yml` | mainly used to set a seed for consistent <br> random number generation across all analyses

### `/code`

`utils.R`

- reads in cannonical marker genes used by various scripts
- defines a cluster color palette used throughout all figures
- wrapper to compute trimmed 0-1 scales expression values
- custom function for dot plots (colored by expression and sized by expression frequency)

`merge_zUMIs.R`

- reads and merged UMI counts across 4 plates with Ensembl IDs as rownames,  
and plate IDs and barcodes as column names, separated by "_"

`Rehrauer18_heatmap.R`

- generates Figure 7B heatmap of bulk RNA-seq data of LPM-associated genes,  
up- (Msln, Wt1) or downregulated (Nf2, Bap1) mesothelioma genes, and neg. control (Tubb4a) 

### `/analysis`

`index.Rmd`  

- table of contents with links to corresponding analyses

`0-preprocessing.Rmd`

- calculation of QC metric
- outlier removal and gene and cell filtering
- normalization (log-transformed library-size-normalized counts)

`1-clustering.Rmd`

- `Seurat` clustering across a sequence of resolutions
- exploration of cluster numbers and stability 
- visualization of DE genes and reduced dimensions (t-SNE and UMAP)

`2-annotation.Rmd` 

- manual annotation of 15 `Seurat` clusters into 6 major subpopulations
- dot plots, heatmaps & UMAPs of canonical markers & hox genes
- DE analysis on hand2+ super-cluster merging hand2+ subpopulations

## References

<a name="f1">[1]</a>:
John Blischak, Peter Carbonetto and Matthew Stephens (2019).  
workflowr: A Framework for Reproducible and Collaborative Data Science.  
R package version 1.4.0. https://CRAN.R-project.org/package=workflowr
