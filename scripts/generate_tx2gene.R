#!/usr/bin/env Rscript

# Load required library
if (!requireNamespace("biomaRt", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("biomaRt")
}
library(biomaRt)

# Connect to Ensembl BioMart for human genes
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Retrieve transcript IDs and corresponding gene IDs
tx2gene <- getBM(
  attributes = c("ensembl_transcript_id", "ensembl_gene_id"),
  mart = mart
)

# Optionally remove version numbers from transcript IDs if present
tx2gene$ensembl_transcript_id <- sub("\\..*$", "", tx2gene$ensembl_transcript_id)
tx2gene$ensembl_gene_id <- sub("\\..*$", "", tx2gene$ensembl_gene_id)

# Write to a tab-separated file
write.table(
  tx2gene,
  file = "reference/tx2gene.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

cat("Transcript-to-gene mapping file 'tx2gene.tsv' created successfully.\n")
