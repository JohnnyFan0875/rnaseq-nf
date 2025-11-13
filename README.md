# rnaseq-nf

## Introduction

`rnaseq-nf` contains a modular, containerized (Docker) RNA-Seq analysis pipeline implemented using `Nextflow`. It is designed to process bulk RNA-Seq data from raw FASTQ files through quality control, quantification, and downstream analysis including differential expression (DE), gene set enrichment analysis (GSEA), and over-representation analysis (ORA).

## Pipeline Summary

Input: Paired-end FASTQ files + project metadata/config files  
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
├── main.nf                                   # Main Nextflow pipeline script
├── nextflow.config                           # System-wide Nextflow configuration
├── modules/                                  # Modular Nextflow process scripts
├── images/                                   # Dockerfiles
├── scripts/                                  # Bash and R scripts
├── reference/                                # GTF, reference fasta, kallisto index files
├── data/
│   └── project_name/
│       ├── raw/   
│       │   ├── test_control_R1.fastq.gz      # Raw FASTQ files
│       │   ├── test_control_R2.fastq.gz      # Raw FASTQ files
│       │   ├── test_treatment_R2.fastq.gz    # Raw FASTQ files
│       │   └── test_treatment_R2.fastq.gz    # Raw FASTQ files
│       ├── config/
│       │   └── project_config.yaml           # Project-specific YAML config
│       ├── metadata/
│       │   └── sample_metadata.tsv           # Tab-separated metadata file
│       └── results/
│           ├── fastqc_pre/                   # FASTQC results (before trimming)
│           ├── fastp_trim/                   # Trimming FASTQ files
│           ├── fastqc_post/                  # FASTQC results (after trimming)
│           ├── multiqc/                      # Aggregated QC reports
│           ├── kallisto_quant/               # Transcript quantification
│           ├── de_analysis/                  # Differential expression analysis
│           ├── gsea_results/                 # GSEA enrichment results
│           └── ora_results/                  # Over-representation analysis results
└── README.md
```

## Install

```bash
bash scripts/install.sh
```

## Customize project metadata/config

1. Copy project config/metadata files

   ```bash
    project_name="project_name" # replace your desired name
   ```

   ```bash
    mkdir -p data/$project_name/{config,metadata,raw}
    cp template/project_config.yaml data/$project_name/project_config.yaml
    cp template/sample_metadata.tsv data/$project_name/sample_metadata.tsv
    unset project_name
   ```

2. Prepare raw data

   - Location: data/`project_name`/raw
   - Require **pair-end** fastq.gz files

3. Modify sample_metadata.csv

   - Location: data/`project_name`/sample_metadata.csv
   - Format: fastq_R1,fastq_R2,sample_name,group_name,replicate
   - Example: test_R1.fastq.gz,test_R2.fastq.gz,test,control,1

4. Modify project_config.yaml

   - Location: data/`project_name`/project_config.yaml

## Usage

```bash
nextflow run main.nf --project_name <project_name> -with-docker [-resume]
```

- `-resume` (optional)
  - Resume the workflow from where it left off in a previous run.
  - It skips completed tasks and avoids re-running steps.

## Others

### GMT Resouces

- [https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.5.1/](https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.5.1/)
