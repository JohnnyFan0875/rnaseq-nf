# RNA-seq Workflow Directory Structure

```bash
rna_seq_pipeline/
├── main.nf
├── nextflow.config
├── modules/
│   ├── fastqc.nf
│   └── fastp.nf 
├── envs/
├── data/
│   └── project_1/
│       ├── raw/
│       │   ├── control_rep1_R1.fastq.gz
│       │   ├── control_rep1_R2.fastq.gz
│       │   ├── control_rep2_R1.fastq.gz
│       │   ├── control_rep2_R2.fastq.gz
│       │   ├── control_rep3_R1.fastq.gz
│       │   ├── control_rep3_R2.fastq.gz
│       │   ├── treatment_rep1_R1.fastq.gz
│       │   ├── treatment_rep1_R2.fastq.gz
│       │   ├── treatment_rep2_R1.fastq.gz
│       │   ├── treatment_rep1_R2.fastq.gz
│       │   ├── treatment_rep3_R1.fastq.gz
│       │   └── treatment_rep3_R2.fastq.gz
│       ├── config/
│       │   └── config.yaml
│       ├── metadata/
│       │   └── sample_metadata.tsv
│       └── results/
│           ├── fastqc_pre/
│           ├── fastp_trim/
│           │   ├── control_rep1_R1.trimmed.fastq.gz
│           │   ├── control_rep1_R2.trimmed.fastq.gz
│           │   ├── control_rep2_R1.trimmed.fastq.gz
│           │   ├── control_rep2_R2.trimmed.fastq.gz
│           │   ├── control_rep3_R1.trimmed.fastq.gz
│           │   ├── control_rep3_R2.trimmed.fastq.gz
│           │   ├── treatment_rep1_R1.trimmed.fastq.gz
│           │   ├── treatment_rep1_R2.trimmed.fastq.gz
│           │   ├── treatment_rep2_R1.trimmed.fastq.gz
│           │   ├── treatment_rep1_R2.trimmed.fastq.gz
│           │   ├── treatment_rep3_R1.trimmed.fastq.gz
│           │   └── treatment_rep3_R2.trimmed.fastq.gz
│           ├── fastqc_post/
│           └── multiqc/

### gmt resouce
https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.5.1/