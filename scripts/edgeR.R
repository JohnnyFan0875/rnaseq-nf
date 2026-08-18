#!/usr/bin/env Rscript

library(optparse)
library(tximport)
library(edgeR)
library(readr)
library(jsonlite)
library(org.Hs.eg.db)
library(AnnotationDbi)

option_list <- list(
  make_option("--input_mode", type = "character"),
  make_option("--quant_type", type = "character"),
  make_option("--expression_inputs", type = "character"),
  make_option("--metadata", type = "character"),
  make_option("--comparisons", type = "character"),
  make_option("--out_dir", type = "character"),
  make_option("--gtf", type = "character"),
  make_option("--script_dir", type = "character")
)

opt <- parse_args(OptionParser(option_list = option_list))
opt$expression_inputs <- strsplit(opt$expression_inputs, ",")[[1]]
cat('Starting edgeR differential expression analysis\n')

source(file.path(opt$script_dir, "plot_utils_de.R"))
source(file.path(opt$script_dir, "gtf_utils.R"))

# Read metadata
metadata <- read_csv(opt$metadata)
rownames(metadata) <- metadata$sample_id

# Parse comparisons JSON string
comparisons <- fromJSON(opt$comparisons)
comparisons <- split(comparisons, seq(nrow(comparisons)))

# Build the transcript-to-gene map directly from the reference GTF
transcript_gene_map <- build_transcript_gene_map_from_gtf(opt$gtf)

if (opt$input_mode == "tximport") {
  sample_ids_from_files <- sapply(opt$expression_inputs, function(x) {
    parts <- strsplit(x, "/")[[1]]
    parts[length(parts) - 1]
  })

  if (!all(sample_ids_from_files %in% metadata$sample_id)) {
    stop("Some sample IDs from expression_inputs do not match metadata$sample_id")
  }

  files <- setNames(opt$expression_inputs, sample_ids_from_files)
  txi <- tximport(files, type = opt$quant_type, tx2gene = transcript_gene_map)
  metadata <- metadata[colnames(txi$counts), ]
  count_matrix <- txi$counts
} else if (opt$input_mode == "gene_counts") {
  if (length(opt$expression_inputs) != 1) {
    stop("gene_counts mode expects exactly one featureCounts matrix")
  }

  counts_tbl <- read.delim(opt$expression_inputs[[1]], comment.char = "#", check.names = FALSE)
  count_columns <- setdiff(colnames(counts_tbl), c("Geneid", "Chr", "Start", "End", "Strand", "Length"))
  sample_ids_from_files <- sub("\\.(Aligned\\.sortedByCoord\\.out|sorted)\\.bam$", "", basename(count_columns))

  if (!all(sample_ids_from_files %in% metadata$sample_id)) {
    stop("Some sample IDs from featureCounts columns do not match metadata$sample_id")
  }

  count_matrix <- as.matrix(counts_tbl[, count_columns, drop = FALSE])
  rownames(count_matrix) <- counts_tbl$Geneid
  colnames(count_matrix) <- sample_ids_from_files
  metadata <- metadata[colnames(count_matrix), ]
} else {
  stop(sprintf("Unsupported input_mode: %s", opt$input_mode))
}

# Extract gene IDs from rownames
gene_ids <- rownames(count_matrix)

# Map Ensembl IDs to gene symbols
gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = gene_ids,
                       column = "SYMBOL",
                       keytype = "ENSEMBL",
                       multiVals = "first")

# Create a data frame for annotation
annotation_df <- data.frame(gene_id = gene_ids, gene_symbol = gene_symbols, row.names = gene_ids)

# Create output directory
if (!dir.exists(opt$out_dir)) dir.create(opt$out_dir, recursive = TRUE)

# Prepare DGEList for all samples to generate exploratory plots
group_all <- factor(metadata$group)
y_all <- DGEList(counts = count_matrix, group = group_all)

# Filter lowly expressed genes, and normalize
design_all <- model.matrix(~group_all)
keep_all <- filterByExpr(y_all, design_all)
y_all <- y_all[keep_all, , keep.lib.sizes=FALSE]

y_all <- calcNormFactors(y_all)
logCPM_all <- cpm(y_all, log=TRUE)

# Generate exploratory plots
plot_mds(y_all, group_all, file.path(opt$out_dir, "MDS_plot.png"))
plot_pca(logCPM_all, group_all, file.path(opt$out_dir, "PCA_plot.png"))
plot_heatmap(logCPM_all, group_all, file.path(opt$out_dir, "Heatmap_top500_var_genes.png"), top_n = 500)

# Run edgeR
for (cmp in comparisons) {

  cat(sprintf("Running comparison: %s\n", cmp$name))
  
  # Create group factor
  group <- factor(metadata$group, levels = c(cmp$control, cmp$non_control))
  design <- model.matrix(~group)
  
  y <- DGEList(counts = count_matrix, group = group)
  
  # Filter lowly expressed genes
  keep <- filterByExpr(y, design)
  y <- y[keep, , keep.lib.sizes=FALSE]
  
  y <- calcNormFactors(y)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit)
  res <- topTags(qlf, n = Inf)$table

  # add annotation
  res$gene_symbol <- annotation_df[rownames(res), "gene_symbol"]
  
  # Write results to CSV file
  res <- cbind(gene_id = rownames(res), res)
  rownames(res) <- NULL

  write.csv(res, file = file.path(opt$out_dir, paste0(cmp$name, "_edgeR_results.csv")), row.names = FALSE)

  # Compute volcano plot
  plot_enhanced_volcano(res, out_file = file.path(opt$out_dir, paste0(cmp$name, "_EnhancedVolcano.png")), title = paste("Volcano Plot:", cmp$name))
  
  # MA plot
  baseMean <- rowMeans(cpm(y, log=FALSE))

  plot_ma(res, file.path(opt$out_dir, paste0(cmp$name, "_MAplot.png")), baseMean = baseMean)

}

cat('edgeR analysis completed\n')
