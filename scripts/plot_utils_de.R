library(ggplot2)
library(matrixStats)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(EnhancedVolcano)

plot_mds <- function(dge_list, group, out_file) {
  png(out_file, width=800, height=600)
  plotMDS(dge_list, col=as.numeric(group), pch=20, main="Multidimensional Scaling (MDS) Plot")
  legend("topleft", legend=levels(group), col=seq_along(levels(group)), pch=20, border = NA)
  dev.off()
}

plot_pca <- function(log_cpm_matrix, group, out_file) {
  pca <- prcomp(t(log_cpm_matrix), center=TRUE, scale.=TRUE)
  percentVar <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 1)
  pca_data <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], group = group)
  p <- ggplot(pca_data, aes(PC1, PC2, color=group, label=group)) +
    geom_point(size=3) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    ggtitle("PCA Plot") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(out_file, p, width=8, height=6, bg = "white")
}

plot_heatmap <- function(log_cpm_matrix, group, out_file, top_n=100) {

  names(group) <- colnames(log_cpm_matrix)
  var_genes <- order(rowVars(log_cpm_matrix), decreasing=TRUE)[1:top_n]
  mat <- log_cpm_matrix[var_genes, , drop = FALSE]

  group <- factor(group[colnames(mat)])
  n_groups <- length(levels(group))
  group_colors <- setNames(brewer.pal(max(3, n_groups), "Pastel1")[1:n_groups], levels(group))
  
  annotation_col <- HeatmapAnnotation(
    Group = group,
    col = list(Group = group_colors),
    na_col = "grey"
  )

  col_fun <- colorRamp2(c(min(mat), mean(mat), max(mat)), c("blue", "white", "red"))
  
  png(out_file, width=1000, height=800)

  ht <- Heatmap(mat,
                name = "expression",
                top_annotation = annotation_col,
                show_row_names = FALSE,
                cluster_rows = TRUE,
                cluster_columns = TRUE,
                col = col_fun,
                column_title = paste("Top", top_n, "Variable Genes Heatmap"))
  
  draw(ht)
  
  dev.off()
}

plot_enhanced_volcano <- function(res_df, out_file,
                                  logFC_col = "logFC",
                                  pval_col = "FDR",
                                  gene_label_col = "gene_symbol",
                                  pCutoff = 0.05,
                                  FCcutoff = 1,
                                  top_n_labels = 20,
                                  title = "Volcano Plot") {

  # Ensure p-value column exists and is numeric
  if (!all(c(logFC_col, pval_col, gene_label_col) %in% colnames(res_df))) {
    stop("Required columns not found in results data frame for volcano plot")
  }
  
  # Optionally select top genes to label
  sig_genes <- res_df[res_df[[pval_col]] < pCutoff & abs(res_df[[logFC_col]]) >= FCcutoff, ]
  top_labels <- head(sig_genes[order(sig_genes[[pval_col]]), gene_label_col], top_n_labels)
  
  # Create volcano plot
  p <- EnhancedVolcano(res_df,
                       lab = res_df[[gene_label_col]],
                       x = logFC_col,
                       y = pval_col,
                       pCutoff = pCutoff,
                       FCcutoff = FCcutoff,
                       selectLab = top_labels,
                       title = title,
                       xlab = bquote(~Log[2]~ "fold change"),
                       ylab = bquote(~-Log[10]~ "p-value"),
                       gridlines.major = FALSE,
                       gridlines.minor = FALSE)
  
  ggsave(out_file, p, width = 10, height = 8, bg = "white")
}

plot_ma <- function(res_df, out_file, 
                    logFC_threshold = 1, 
                    FDR_threshold = 0.05, 
                    baseMean = NULL, 
                    gene_label_col = "gene_symbol", 
                    pval_col = "FDR",
                    logFC_col = "logFC") {
  
  # Ensure p-value column exists and is numeric
  if (!all(c(logFC_col, pval_col, gene_label_col) %in% colnames(res_df))) {
    stop("Required columns not found in results data frame for MA plot")
  }

  if (is.null(baseMean)) {
    stop("Please provide baseMean (mean expression) vector for MA plot.")
  }

  # Define regulation status with descriptive labels
  regulation <- rep("Not Significant", nrow(res_df))
  is_sig <- res_df[[pval_col]] < FDR_threshold & abs(res_df[[logFC_col]]) >= logFC_threshold
  regulation[is_sig & res_df[[logFC_col]] > 0] <- "Upregulated"
  regulation[is_sig & res_df[[logFC_col]] < 0] <- "Downregulated"

  # Create data frame for plotting
  ma_df <- data.frame(
    A = log2(baseMean + 1),  # average expression (A)
    M = res_df$logFC,        # log fold change (M)
    Regulation = factor(regulation, levels = c("Not Significant", "Upregulated", "Downregulated")),   
    gene_symbol = res_df$gene_symbol
  )

  colors <- c("Not Significant" = "grey", "Upregulated" = "red", "Downregulated" = "blue")
  
  p <- ggplot(ma_df, aes(x = A, y = M, color = Regulation)) +
    geom_point(alpha = 0.6, size = 1) +
    scale_color_manual(values = colors) +
    geom_hline(yintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed", color = "blue") +
    labs(title = "MA Plot",
         x = "Average Log2 Expression (A)",
         y = "Log2 Fold Change (M)",
         color = NULL) +
    theme_minimal() +
    theme(legend.position = "top")
  
  ggsave(out_file, p, width = 8, height = 6, bg = "white")
}