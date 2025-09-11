# rnaseq-nf

## Introduction

`rnaseq-nf` contains a modular, containerized (Docker) RNA-Seq analysis pipeline implemented using `Nextflow`. It is designed to process bulk RNA-Seq data from raw FASTQ files through quality control, quantification, and downstream analysis including differential expression (DE), gene set enrichment analysis (GSEA), and over-representation analysis (ORA).

## Pipeline Summary

Input: Paired-end FASTQ files + project metadata/config  
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

1. Install Nextflow (>=24.10.5) from [official website](https://www.nextflow.io/docs/latest/install.html#install-nextflow)
2. Install Docker (>=28.0.1) from [official website](https://docs.docker.com/desktop/)
3. Install required package
   ```bash
   apt install kallisto gffread
   ```
4. Create Docker Image
   ```bash
   docker build -t de_analysis -f containers/de_analysis.dockerfile containers
   docker build -t enrichment_analysis -f containers/enrichment.dockerfile containers
   ```
5. Create Reference Files

   1. kallisto index

   ```bash
   # Download GTF file
   wget -P reference/ https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz

   # Download FASTA file
   wget -P reference/ https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

   # Unzip gzip files
   gunzip reference/Homo_sapiens.GRCh38.113.gtf.gz
   gunzip reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

   # Extract transcript sequence
   gffread -w reference/transcripts.GRCh38.113.fa -g reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz reference/Homo_sapiens.GRCh38.113.gtf.gz

   # Build Kallisto index
   kallisto index -i reference/transcripts.GRCh38.113.fa.idx reference/transcripts.GRCh38.113.fa

   # Remove unrelated files
   rm reference/transcripts.GRCh38.113.fa reference/Homo_sapiens.GRCh38.113.gtf
   ```

   - alternatives: [kb ref](https://www.kallistobus.tools/kb_usage/kb_ref/)

   2. tx2gene.tsv (optional)  
      The reference file tx2gene.tsv is already in reference folder. If you need the latest tx2gene file, you can run

      ```bash
      Rscript scripts/generate_tx2gene.R
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

   - Put your raw fastq.gz files into `data/<project_name>/raw`
   - Require pair-end fastq.gz files
   - For example,
     - control_rep1_R1.fastq.gz, control_rep1_R2.fastq.gz
     - control_rep2_R1.fastq.gz, control_rep2_R2.fastq.gz
     - control_rep3_R1.fastq.gz, control_rep3_R2.fastq.gz
     - treatment_rep1_R1.fastq.gz, treatment_rep1_R2.fastq.gz
     - treatment_rep2_R1.fastq.gz, treatment_rep2_R2.fastq.gz
     - treatment_rep3_R1.fastq.gz, treatment_rep3_R2.fastq.gz

3. Modify sample_metadata.tsv

   - Location: data/<project_name>/sample_metadata.tsv
   - sample_id = name of raw fastq file without extension
   - Example:

     | sample_id      | group     | replicate |
     | -------------- | --------- | --------- |
     | control_rep1   | control   | 1         |
     | control_rep2   | control   | 2         |
     | control_rep3   | control   | 3         |
     | treatment_rep1 | treatment | 1         |
     | treatment_rep2 | treatment | 2         |
     | treatment_rep3 | treatment | 3         |

4. Modify project_config.yaml

   - Location: data/<project_name>/project_config.yaml
   - Example:

   ```text
    # project information
    metadata_file: ./data/<project_name>/metadata/sample_metadata.tsv
    comparisons:
    - name: treatment_vs_control
        control: control
        non_control: treatment

    # quantification
    kallisto_index: ./reference/transcriptome_index_v13.idx
    kallisto_threads: 4

    # differential expression analysis
    de_method: "edgeR"
   ```

   - comparisons: names is used for output file name
   - comparisons: control/non-control is referred to `sample_metadata.tsv` **group** column

## Usage

```bash
nextflow run main.nf --project_name <project_name> -with-docker -resume
```

- `-resume`: optional, resume the workflow from where it left off in a previous run. It skips completed tasks and avoids re-running steps

## Others

### GMT Resouces

- [https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.5.1/](https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.5.1/)
