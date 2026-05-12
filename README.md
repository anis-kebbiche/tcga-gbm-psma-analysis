# TCGA GBM PSMA Analysis Pipeline

This repository contains an automated R pipeline to download, process, and analyze RNA-seq and clinical data for Glioblastoma (GBM) and Lower-Grade Glioma (LGG) from the TCGA database. 

The analysis specifically focuses on PSMA (FOLH1) expression, generating survival analyses, differential expression (DESeq2), and scoring for hypoxia and angiogenesis signatures.

## How to Use This Pipeline

To reproduce the analysis, you need to run the two R scripts in the following order:

### Step 1: Data Preparation
**File:** `Data_preparation_pipeline_TCGA_GBM_LGG.R`

Run this script first. It will automatically connect to the TCGA database via `TCGAbiolinks`, download the required clinical and RNA-seq (STAR - Counts) data, filter the patients, and generate the clean `.csv` matrices needed for the analysis.

**⚠️ IMPORTANT DATA PREREQUISITE:**
This script requires a master Excel tracking file to be placed in the same root directory before running. This file is a supplementary dataset from the following publication:
* **Article:** Impact of an evolving classification system on diffuse glioma repositories: experience from the Sydney brain tumour bank
* **Authors:** Laveniya Satgunaseelan et al.
* **Link/DOI:** https://doi.org/10.1007/s11060-026-05470-1
* **Action Required:** Download the supplementary material 3 (named `11060_2026_5470_MOESM3_ESM.xlsx`), ensure the filename matches exactly, and place it in the same folder as this repository's scripts.

### Step 2: Main Analysis
**File:** `TCGA_GBM_PSMA_analysis_final.R`

Once the data preparation is complete, run this script. It loads the cleaned datasets and performs the full downstream analysis, including:
* Propensity Score Matching (PSM) for Initial vs. Recurrent tumors.
* Differential Gene Expression (DESeq2) and Volcano plots.
* Survival analysis (Kaplan-Meier & Cox Proportional Hazards).
* Gene Ontology (GO) and g:Profiler enrichment analysis.
* Z-score heatmaps for Hallmark Hypoxia and Angiogenesis gene signatures.

All plots (`.png`/`.pdf`) and statistical summary tables (`.csv`/`.txt`) will be automatically exported to a dedicated output folder.

## Requirements & Packages
Both scripts include auto-installers that will check for and install missing packages. Key dependencies include:
* **CRAN:** `dplyr`, `ggplot2`, `survival`, `survminer`, `MatchIt`, `gprofiler2`, `patchwork`
* **Bioconductor:** `TCGAbiolinks`, `DESeq2`, `clusterProfiler`, `org.Hs.eg.db`, `SummarizedExperiment`

## ⚠️ Data Storage Note
To keep this repository lightweight and comply with GitHub's file size limits, the raw TCGA matrices and clinical datasets downloaded during Step 1 are not tracked in this repository (via `.gitignore`). Anyone cloning this repository simply needs to run Step 1 to generate the data locally.
