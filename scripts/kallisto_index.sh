#!/bin/bash

# Define file paths and names
GTF_GZ="reference/Homo_sapiens.GRCh38.113.gtf.gz"
FASTA_GZ="reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
GTF="reference/Homo_sapiens.GRCh38.113.gtf"
FASTA="reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
TRANSCRIPTS="reference/transcripts.GRCh38.113.fa"
KALLISTO_IDX="reference/transcripts.GRCh38.113.fa.idx"

# Create reference directory if it doesn't exist
mkdir -p reference

# Download the GTF and FASTA files from Ensembl FTP
wget -P reference/ https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz
wget -P reference/ https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# Decompress the downloaded files while keeping the original compressed files
gunzip -k "$GTF_GZ"
gunzip -k "$FASTA_GZ"

# Generate transcript fasta file from GTF and genome fasta using gffread in Docker
docker run -it --rm -v "$(pwd)/reference/:/mnt/" zavolab/gffread:0.11.7-slim \
    gffread -w /mnt/$(basename "$TRANSCRIPTS") -g /mnt/$(basename "$FASTA") /mnt/$(basename "$GTF")

# Build kallisto index from the generated transcripts fasta using kallisto Docker image
docker run -it --rm -v "$(pwd)/reference/:/mnt/" quay.io/biocontainers/kallisto:0.50.1--h6de1650_2 \
    kallisto index -i /mnt/$(basename "$KALLISTO_IDX") /mnt/$(basename "$TRANSCRIPTS")

# Clean up temporary files
rm "$FASTA"
rm "$GTF_GZ"
rm "$FASTA_GZ"
