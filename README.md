# rnaseq-nf

## Introduction

`rnaseq-nf` is a modular, containerized RNA-seq pipeline built with Nextflow DSL2.
It processes paired-end FASTQ files through QC, trimming, quantification,
differential expression, and optional enrichment analysis.

## Pipeline Summary

Input:

- Paired-end FASTQ files
- Sample metadata CSV
- Params YAML (`-params-file`)

Output:

- FastQC / fastp / MultiQC reports
- Kallisto quantification
- Differential expression results
- Optional GSEA / ORA results

Main steps:

1. FastQC (pre-trim)
2. fastp trimming
3. FastQC (post-trim)
4. Kallisto quantification
5. MultiQC aggregation
6. Differential expression (edgeR)
7. Optional GSEA / ORA

## Project Structure

```bash
rna_seq_pipeline/
├── main.nf
├── nextflow.config
├── conf/
│   ├── base.config
│   ├── modules.config
│   ├── slurm.config
│   └── test.config
├── modules/
├── scripts/
├── params/
│   ├── test.yaml
│   └── *.yaml
├── template/
│   ├── params_template.yaml
│   └── metadata_template.csv
├── data/
│   └── <project_name>/
│       ├── raw/
│       ├── metadata/
│       │   └── sample_metadata.csv
│       └── results/
└── reference/
```

## Requirements

- Nextflow `>= 24.10.0`
- One execution environment:
  - Docker (recommended), or
  - Singularity, or
  - Local binaries installed (`fastqc`, `fastp`, `kallisto`, `multiqc`, `Rscript`)

## Install

```bash
bash scripts/install.sh
```

Reference index helper:

- [scripts/kallisto_index.sh](./scripts/kallisto_index.sh)

## Prepare Input Data

1. Create project folders and copy templates

```bash
project_name="my_project"
mkdir -p data/${project_name}/{raw,metadata,results}
cp template/params_template.yaml params/${project_name}.yaml
cp template/metadata_template.csv data/${project_name}/metadata/sample_metadata.csv
unset project_name
```

2. Put FASTQ files in:

- `data/<project_name>/raw/`

3. Edit metadata CSV:

- `data/<project_name>/metadata/sample_metadata.csv`

Metadata columns (required):

- `sample_id`
- `fastq1`
- `fastq2`
- `group`

Example row:

```csv
sample_control_1,sample_control_1_R1.fastq.gz,sample_control_1_R2.fastq.gz,control
```

## Configure Parameters

Edit your copied params file:

- `params/<project_name>.yaml`

Minimum required params:

- `project_name`
- `outdir`
- `comparisons`

Important notes:

- If `metadata_file` is omitted, default is:
  - `data/<project_name>/metadata/sample_metadata.csv`
- `de_method` currently supports `edgeR` in this repo.

## Usage

### Standard run (Docker)

```bash
nextflow run main.nf -profile docker -params-file params/<project_name>.yaml
```

### Resume run

```bash
nextflow run main.nf -profile docker -params-file params/<project_name>.yaml -resume
```

### Local run (no container)

```bash
nextflow run main.nf -profile local -params-file params/<project_name>.yaml
```

## Test Dataset

This repo already includes test data under `data/test` and a test params file at
`params/test.yaml`.

Run test:

```bash
nextflow run main.nf -profile docker,test -params-file params/test.yaml
```

Notes:

- `conf/test.config` is intentionally minimal and mainly reserved for optional
  profile-level overrides (resource / executor / container behavior).
- Dataset paths and biological params should stay in `params/test.yaml`.

## Profiles

- `docker`: run with Docker containers
- `singularity`: run with Singularity
- `local`: run on host without container
- `slurm`: adds HPC SLURM executor settings (typically with `-profile singularity,slurm`)
- `test`: optional profile-specific overrides for testing

## Outputs

All outputs are written under `outdir` defined in your params file, typically:

- `fastqc_pre/`
- `fastp_trim/`
- `fastqc_post/`
- `kallisto_quant/`
- `multiqc/`
- `de_analysis/`
- `enrichment_analysis/` (when enrichment is enabled)
