#!/usr/bin/env Rscript

library(optparse)
library(tximport)
library(edgeR)
library(readr)
library(jsonlite)
library(org.Hs.eg.db)
library(AnnotationDbi)

option_list <- list(
  make_option("--quant_dir", type = "character"),
  make_option("--metadata", type = "character"),
  make_option("--comparisons", type = "character"),
  make_option("--out_dir", type = "character"),
  make_option("--ref_dir", type = "character"),
  make_option("--script_dir", type = "character")
)

opt <- parse_args(OptionParser(option_list = option_list))
cat('Starting edgeR differential expression analysis\n')

source(file.path(opt$script_dir, "plot_utils_de.R"))

# Read metadata
metadata <- read_tsv(opt$metadata)
rownames(metadata) <- metadata$sample_id

# Construct file paths for abundance files
files <- setNames(
  file.path(opt$quant_dir, metadata$sample_id, "abundance.tsv"),
  metadata$sample_id
)

# Parse comparisons JSON string
comparisons <- fromJSON(opt$comparisons)
comparisons <- split(comparisons, seq(nrow(comparisons)))

# Load tx2gene mapping
tx2gene <- read_tsv(file.path(opt$ref_dir, "tx2gene.tsv"), col_names = c("transcript_id", "gene_id")) # reconfirm position

# Import transcript-level estimates
txi <- tximport(files, type = "kallisto", tx2gene = tx2gene)

# Extract gene IDs from rownames
gene_ids <- rownames(txi$counts)

# Map Ensembl IDs to gene symbols
gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = gene_ids,
                       column = "SYMBOL",
                       keytype = "ENSEMBL",
                       multiVals = "first")

# Create a data frame for annotation
annotation_df <- data.frame(gene_id = gene_ids, gene_symbol = gene_symbols, row.names = rownames(txi$counts))

# Create output directory
if (!dir.exists(opt$out_dir)) dir.create(opt$out_dir, recursive = TRUE)

# Prepare DGEList for all samples to generate exploratory plots
group_all <- factor(metadata$group)
y_all <- DGEList(counts = txi$counts, group = group_all)

# Filter lowly expressed genes, and normalize
design_all <- model.matrix(~group_all)
keep_all <- filterByExpr(y_all, design_all)
y_all <- y_all[keep_all, , keep.lib.sizes=FALSE]

y_all <- calcNormFactors(y_all)
logCPM_all <- cpm(y_all, log=TRUE)

# Generate exploratory plots
plot_mds(y_all, group_all, file.path(opt$out_dir, "MDS_plot.png"))
plot_pca(logCPM_all, group_all, file.path(opt$out_dir, "PCA_plot.png"))
plot_heatmap(logCPM_all, group_all, file.path(opt$out_dir, "Heatmap_top500_var_genes.png"))

# Run edgeR
for (cmp in comparisons) {

  cat(sprintf("Running comparison: %s\n", cmp$name))
  
  # Create group factor
  group <- factor(metadata$group, levels = c(cmp$control, cmp$non_control))
  design <- model.matrix(~group)
  
  y <- DGEList(counts = txi$counts, group = group)
  
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

  write.csv(res, file = file.path(opt$out_dir, paste0(cmp$name, "_edgeR_results.csv")), row.names = FALSE)}

  # Compute volcano plot
  plot_enhanced_volcano(res, out_file = file.path(opt$out_dir, paste0(cmp$name, "_EnhancedVolcano.png")), title = paste("Volcano Plot:", cmp$name))
  
  # MA plot
  baseMean <- rowMeans(cpm(y, log=FALSE))

  plot_ma(res, file.path(opt$out_dir, paste0(cmp$name, "_MAplot.png")), baseMean = baseMean)

cat('edgeR analysis completed\n')
