#!/usr/bin/env Rscript

library(optparse)

option_list <- list(
  make_option("--output", type = "character", default = "reference/tx2gene.tsv")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Load required library
if (!requireNamespace("biomaRt", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("biomaRt")
}
library(biomaRt)

output_dir <- dirname(opt$output)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

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
  file = opt$output,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

cat(sprintf("Transcript-to-gene mapping file created successfully: %s\n", opt$output))
