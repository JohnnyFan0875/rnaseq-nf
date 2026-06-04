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
│       ├── logs/
│       ├── metadata/
│       │   └── sample_metadata.csv
│       └── results/
└── reference/
```

## Requirements

- Nextflow `>= 24.10.0`
- One execution environment:
  - Docker (recommended), or
  - Local binaries installed (`fastqc`, `fastp`, `kallisto`, `multiqc`, `Rscript`)

## Prepare Input Data

### Metadata

Template file:

- [template/sample_metadata.csv](/mnt/d/shared/rna_seq_pipeline/template/sample_metadata.csv:1)

Example row:

```csv
sample_control_1,/data/project/raw/sample_control_1_R1.fastq.gz,/data/project/raw/sample_control_1_R2.fastq.gz,control
```

## Configure Parameter file

Template file:

- [template/params_template.yaml](/mnt/d/shared/rna_seq_pipeline/template/params_template.yaml:1)

Minimum required params:

- `project_name`
- `outdir`
- `comparisons`

## Usage

### Standard run (Docker)

```bash
nextflow run main.nf \
  -profile docker \
  -params-file <param_file>.yaml
```

### Resume run

```bash
nextflow run main.nf \
  -profile docker \
  -params-file <param_file>.yaml \
  -resume
```

### Local run

```bash
nextflow run main.nf \
  -profile local \
  -params-file <param_file>.yaml
```

## Test Dataset

```bash
nextflow run main.nf \
  -profile docker,test \
  -params-file params/test.yaml
```

## Profiles

- `docker`: run with Docker containers
- `local`: run on host without container
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
