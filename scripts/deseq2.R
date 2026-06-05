#!/usr/bin/env Rscript

library(optparse)
library(tximport)
library(DESeq2)
library(readr)
library(jsonlite)
library(org.Hs.eg.db)
library(AnnotationDbi)

option_list <- list(
  make_option("--quant_files", type = "character"),
  make_option("--metadata", type = "character"),
  make_option("--comparisons", type = "character"),
  make_option("--out_dir", type = "character"),
  make_option("--gtf", type = "character"),
  make_option("--script_dir", type = "character")
)

opt <- parse_args(OptionParser(option_list = option_list))
opt$quant_files <- strsplit(opt$quant_files, ",")[[1]]
cat("Starting DESeq2 differential expression analysis\n")

source(file.path(opt$script_dir, "plot_utils_de.R"))
source(file.path(opt$script_dir, "gtf_utils.R"))

metadata <- read_csv(opt$metadata)
rownames(metadata) <- metadata$sample_id

sample_ids_from_files <- sapply(opt$quant_files, function(x) {
  parts <- strsplit(x, "/")[[1]]
  parts[length(parts) - 1]
})

if (!all(sample_ids_from_files %in% metadata$sample_id)) {
  stop("Some sample IDs from quant_files do not match metadata$sample_id")
}

files <- setNames(opt$quant_files, sample_ids_from_files)
comparisons <- fromJSON(opt$comparisons)
comparisons <- split(comparisons, seq(nrow(comparisons)))

transcript_gene_map <- build_transcript_gene_map_from_gtf(opt$gtf)
txi <- tximport(files, type = "kallisto", tx2gene = transcript_gene_map)
metadata <- metadata[colnames(txi$counts), ]

gene_ids <- rownames(txi$counts)
gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = gene_ids,
                       column = "SYMBOL",
                       keytype = "ENSEMBL",
                       multiVals = "first")
annotation_df <- data.frame(gene_id = gene_ids, gene_symbol = gene_symbols, row.names = gene_ids)

if (!dir.exists(opt$out_dir)) dir.create(opt$out_dir, recursive = TRUE)

dds_all <- DESeqDataSetFromTximport(txi = txi, colData = metadata, design = ~ group)
keep_all <- rowSums(counts(dds_all) >= 10) >= 2
dds_all <- dds_all[keep_all, ]
vsd <- vst(dds_all, blind = TRUE)
vsd_mat <- assay(vsd)

plot_pca(vsd_mat, metadata$group, file.path(opt$out_dir, "PCA_plot.png"))
plot_heatmap(vsd_mat, metadata$group, file.path(opt$out_dir, "Heatmap_top500_var_genes.png"), top_n = 500)

sample_dists <- dist(t(vsd_mat))
mds_coords <- cmdscale(sample_dists, k = 2)
mds_df <- data.frame(MDS1 = mds_coords[, 1], MDS2 = mds_coords[, 2], group = metadata$group)
png(file.path(opt$out_dir, "MDS_plot.png"), width = 800, height = 600)
plot(mds_df$MDS1, mds_df$MDS2,
     col = as.numeric(factor(mds_df$group)),
     pch = 20,
     xlab = "MDS1",
     ylab = "MDS2",
     main = "Multidimensional Scaling (MDS) Plot")
legend("topleft",
       legend = levels(factor(mds_df$group)),
       col = seq_along(levels(factor(mds_df$group))),
       pch = 20,
       border = NA)
dev.off()

for (cmp in comparisons) {
  cat(sprintf("Running comparison: %s\n", cmp$name))

  cmp_metadata <- metadata[metadata$group %in% c(cmp$control, cmp$non_control), , drop = FALSE]
  cmp_metadata$group <- factor(cmp_metadata$group, levels = c(cmp$control, cmp$non_control))

  cmp_counts <- txi$counts[, rownames(cmp_metadata), drop = FALSE]
  cmp_abundance <- txi$abundance[, rownames(cmp_metadata), drop = FALSE]
  cmp_length <- txi$length[, rownames(cmp_metadata), drop = FALSE]

  cmp_txi <- list(counts = cmp_counts, abundance = cmp_abundance, length = cmp_length, countsFromAbundance = txi$countsFromAbundance)

  dds <- DESeqDataSetFromTximport(txi = cmp_txi, colData = cmp_metadata, design = ~ group)
  keep <- rowSums(counts(dds) >= 10) >= 2
  dds <- dds[keep, ]
  dds <- DESeq(dds)

  res <- results(dds, contrast = c("group", cmp$non_control, cmp$control))
  res_df <- as.data.frame(res)
  res_df$logCPM <- log2(rowMeans(counts(dds, normalized = TRUE)) + 1)
  res_df$gene_symbol <- annotation_df[rownames(res_df), "gene_symbol"]
  res_df <- cbind(gene_id = rownames(res_df), res_df)
  rownames(res_df) <- NULL

  colnames(res_df)[colnames(res_df) == "log2FoldChange"] <- "logFC"
  colnames(res_df)[colnames(res_df) == "pvalue"] <- "PValue"
  colnames(res_df)[colnames(res_df) == "padj"] <- "FDR"

  write.csv(res_df, file = file.path(opt$out_dir, paste0(cmp$name, "_DESeq2_results.csv")), row.names = FALSE)

  plot_enhanced_volcano(
    res_df,
    out_file = file.path(opt$out_dir, paste0(cmp$name, "_EnhancedVolcano.png")),
    title = paste("Volcano Plot:", cmp$name)
  )

  plot_ma(
    res_df,
    file.path(opt$out_dir, paste0(cmp$name, "_MAplot.png")),
    baseMean = res_df$baseMean
  )
}

cat("DESeq2 analysis completed\n")
