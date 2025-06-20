# rnaseq-nf

## Introduction

`rnaseq-nf` contains a modular, containerized (Docker) RNA-Seq analysis pipeline implemented using `Nextflow`. It is designed to process bulk RNA-Seq data from raw FASTQ files through quality control, quantification, and downstream analysis including differential expression (DE), gene set enrichment analysis (GSEA), and over-representation analysis (ORA).

## Pipeline Summary

Input:  Paired-end FASTQ files + project metadata/config
Output: Quality reports, expression quantification, DEG tables, GSEA/ORA results

Main Steps:
  1. Quality control (FastQC)
  2. Adapter and quality trimming (Fastp)
  3. Quantification (Kallisto)
  4. Post-trimming QC summary (MultiQC)
  5. Differential expression analysis (edgeR) # DESeq2 under contruction
  6. GSEA (Gene Set Enrichment Analysis)
  7. ORA (Over-Representation Analysis)

## Directory Structure

```bash
rna_seq_pipeline/
├── main.nf                               # Main Nextflow pipeline script
├── nextflow.config                       # System-wide Nextflow configuration
├── modules/                              # Modular Nextflow process scripts
├── containers/                           # Dockerfiles
├── data/
│   └── project_name/
│       ├── raw/                          # Raw FASTQ files
│       ├── config/
│       │   └── project_config.yaml       # Project-specific YAML config
│       ├── metadata/
│       │   └── sample_metadata.tsv       # Tab-separated metadata file
│       └── results/
│           ├── fastqc_pre/               # FASTQC results (before trimming)
│           ├── fastp_trim/               # Trimming FASTQ files
│           ├── fastqc_post/              # FASTQC results (after trimming)
│           ├── multiqc/                  # Aggregated QC reports
│           ├── kallisto_quant/           # Transcript quantification
│           ├── de_analysis/              # Differential expression analysis
│           ├── gsea_results/             # GSEA enrichment results
│           └── ora_results/              # Over-representation analysis results
└── README.md
```

## Install

1. Install Nextflow
2. Install Docker
3. Create Docker Image
4. Create Reference Files

## Customize project metadata/config

1. sample_metadata.tsv
2. project_config.yaml

## Usage


### gmt resouce
https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.5.1/