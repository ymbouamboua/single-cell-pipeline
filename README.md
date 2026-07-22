Single-cell RNA-seq Automated Pipeline
================
Yvon Mbouamboua
2026-06-05

## Overview

An end-to-end automated pipeline for single-cell RNA-seq analysis from
raw 10x Genomics matrices to integrated Seurat objects, differential
expression analysis, pathway enrichment, and publication-ready reports.

## Features

- Automatic loading of 10x Genomics datasets (MTX or H5)
- Optional demultiplexing integration (Demuxafy)
- Automated QC filtering
- Optional ambient RNA removal (DecontX / SoupX)
- Batch correction and integration (Harmony, RPCA, CCA)
- Differential expression analysis
- Pathway enrichment analysis
- HTML reporting
- Interactive execution with R
- Fully reproducible execution with Snakemake
- HPC/SLURM compatible

## Project Structure

``` text
project/
├── master_summary.csv
├── Snakefile
├── run_pipeline.R
├── README.Rmd
│
├── config/
│   ├── pipeline_config.yaml
│   └── contrasts.csv
│
├── scripts/
│   ├── 01_load_seurat.R
│   ├── 02_add_demux.R
│   ├── 03_run_qc.R
│   ├── 04_decontaminate.R
│   ├── 05_integrate.R
│   └── 06_run_deg.R
│
├── logs/
│
└── results/
    ├── 01_loaded/
    ├── 02_demuxed/
    ├── 03_qc/
    ├── 04_clean/
    ├── 05_integrated/
    └── 06_deg/
```

## Workflows

### Basic Workflow (Default)

``` text
Load
  ↓
QC
  ↓
Integration
```

Scripts:

- 01_load_seurat.R
- 03_run_qc.R
- 05_integrate.R

``` bash
Rscript run_pipeline.R --workflow basic
```

### Demultiplex Workflow

``` text
Load
  ↓
Demultiplex
  ↓
QC
  ↓
Integration
```

Scripts:

- 01_load_seurat.R
- 02_add_demux.R
- 03_run_qc.R
- 05_integrate.R

``` bash
Rscript run_pipeline.R --workflow demux
```

### Full Workflow

``` text
Load
  ↓
Demultiplex (optional)
  ↓
QC
  ↓
Decontamination
  ↓
Integration
  ↓
DEG
```

Scripts:

- 01_load_seurat.R
- 02_add_demux.R
- 03_run_qc.R
- 04_decontaminate.R
- 05_integrate.R
- 06_run_deg.R

``` bash
Rscript run_pipeline.R --workflow full
```

## Quick Start

``` bash
Rscript run_pipeline.R --workflow basic
```

``` bash
Rscript run_pipeline.R --workflow full
```

``` bash
Rscript run_pipeline.R --workflow basic --sample GSM7476098
```

## Run Individual Steps

``` bash
Rscript run_pipeline.R --step load
Rscript run_pipeline.R --step qc
Rscript run_pipeline.R --step integrate
Rscript run_pipeline.R --step deg
```

## Snakemake

``` bash
snakemake -n
snakemake --cores 8
snakemake --cores 8 --use-conda
```

## master_summary.csv

Required columns:

| Column          | Description                      |
|-----------------|----------------------------------|
| SampleID        | Sample identifier                |
| Species         | human or mouse                   |
| CellRanger_Dir  | Path to 10x matrix directory     |
| Demuxafy_Dir    | Optional Demuxafy directory      |
| BulkMapping_TSV | Optional donor mapping           |
| Min_Feat        | Minimum detected genes           |
| Min_UMI         | Minimum UMIs                     |
| Max_Mito        | Maximum mitochondrial percentage |
| MAD_N           | MAD multiplier                   |
| Rm_Dbl          | Run scDblFinder                  |
| Dbl_Score       | Doublet threshold                |

## Integration Methods

- Harmony (recommended)
- RPCA
- CCA
- Merge

``` yaml
integration_method: Harmony
```

## Differential Expression

``` yaml
deg_method: pseudo_bulk
```

Supported:

- pseudo_bulk
- FindMarkers
- both

## Reproducibility

``` r
renv::snapshot()
```

Store:

- pipeline_config.yaml
- master_summary.csv
- contrasts.csv

for complete reproducibility.
