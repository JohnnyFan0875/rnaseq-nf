#!/bin/bash

# Define file paths
GTF_GZ="reference/Homo_sapiens.GRCh38.113.gtf.gz"
FASTA_GZ="reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
GTF="reference/Homo_sapiens.GRCh38.113.gtf"
FASTA="reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
TRANSCRIPTS="reference/transcripts.GRCh38.113.fa"
KALLISTO_IDX="reference/transcripts.GRCh38.113.fa.idx"

# Create reference directory if not exists
mkdir -p reference

# Download GTF if not exists
if [ ! -f "$GTF_GZ" ] && [ ! -f "$GTF" ]; then
    wget -P reference/ https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz
else
    echo "GTF file already downloaded."
fi

# Download FASTA if not exists
if [ ! -f "$FASTA_GZ" ] && [ ! -f "$FASTA" ]; then
    wget -P reference/ https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
else
    echo "FASTA file already downloaded."
fi

# Unzip GTF if uncompressed GTF not exists but gz exists
if [ ! -f "$GTF" ] && [ -f "$GTF_GZ" ]; then
    gunzip -k "$GTF_GZ"
else
    echo "GTF file already uncompressed."
fi

# Unzip FASTA if uncompressed FASTA not exists but gz exists
if [ ! -f "$FASTA" ] && [ -f "$FASTA_GZ" ]; then
    gunzip -k "$FASTA_GZ"
else
    echo "FASTA file already uncompressed."
fi

# Extract transcript sequence if output not exists
if [ ! -f "$TRANSCRIPTS" ]; then
    gffread -w "$TRANSCRIPTS" -g "$FASTA" "$GTF"
else
    echo "Transcript FASTA already exists."
fi

# Build Kallisto index if not exists
if [ ! -f "$KALLISTO_IDX" ]; then
    kallisto index -i "$KALLISTO_IDX" "$TRANSCRIPTS"
else
    echo "Kallisto index already exists."
fi

# alternative: kb ref (https://www.kallistobus.tools/kb_usage/kb_ref/)