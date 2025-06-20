#!/usr/bin/env Rscript

library(optparse)
library(jsonlite)
library(tidyverse)
library(data.table)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(ReactomePA)
library(DOSE)
library(magrittr)

option_list <- list(
  make_option("--deg_dir", type="character"),
  make_option("--metadata", type="character"),
  make_option("--out_dir", type="character"),
  make_option("--script_dir", type="character"),
  make_option("--comparisons", type="character"),
  make_option("--de_method", type = "character"),
  make_option("--organism", type = "character", default = "human")
)

opt <- parse_args(OptionParser(option_list=option_list))
cat('Starting ora analysis\n')

source(file.path(opt$script_dir, "plot_utils_enrichment.R"))

# Read metadata
metadata <- read_tsv(opt$metadata)
rownames(metadata) <- metadata$sample_id

# Parse comparisons JSON string
comparisons <- fromJSON(opt$comparisons) %>% split(seq(nrow(.)))

# Function to convert gene symbols to Entrez IDs
convert_to_entrez <- function(genes) {
  gene_df <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  unique(gene_df$ENTREZID)
}

# # for testing in Rstudio
# yaml_data <- yaml::read_yaml("D:\\shared\\rna_seq_pipeline\\data\\project_1\\config\\project_config.yaml")
# json_data <- jsonlite::toJSON(yaml_data, pretty = TRUE, auto_unbox = TRUE)
# opt <- list(
#   deg_dir = "D:\\shared\\rna_seq_pipeline\\data\\project_1\\results\\de_analysis\\DEG_edgeR",
#   metadata = "D:\\shared\\rna_seq_pipeline\\data\\project_1\\metadata\\sample_metadata.tsv",
#   out_dir = "D:\\shared\\rna_seq_pipeline\\data\\project_1\\results\\ora_results",
#   script_dir = "D:\\shared\\rna_seq_pipeline\\scripts",
#   comparisons = json_data,
#   de_method = "edgeR",
#   organism = "human"
# )
# parsed_list <- fromJSON(opt$comparisons)
# comparisons_df <- parsed_list$comparisons
# comparisons <- split(comparisons_df, seq(nrow(comparisons_df)))
# 
# source(file.path(opt$script_dir, "plot_utils_enrichment.R"))
# 
# # Read metadata
# metadata <- read_tsv(opt$metadata)
# rownames(metadata) <- metadata$sample_id
# 
# # Function to convert gene symbols to Entrez IDs
# convert_to_entrez <- function(genes) {
#   gene_df <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
#   unique(gene_df$ENTREZID)
# }

set.seed(123)

# Run GSEA
for (cmp in comparisons) {
  
  deg_file <- file.path(opt$deg_dir, paste0(cmp$name, "_", opt$de_method, "_results.csv"))
  deg_data <- read.csv(deg_file)

  sig_genes <- deg_data %>%
    filter(!is.na(PValue), !is.na(gene_symbol)) %>%
    filter(PValue < 0.05) %>%
    pull(gene_symbol) %>%
    unique()

  # Convert gene symbols to Entrez IDs
  entrez_genes <- convert_to_entrez(sig_genes)

  # Define universe genes (all genes tested)
  universe_genes <- convert_to_entrez(deg_data$gene_symbol)

  # GO Biological Process ORA
  ora_go <- enrichGO(gene = entrez_genes,
                     universe = universe_genes,
                     OrgDb = org.Hs.eg.db,
                     ont = "BP",
                     keyType = "ENTREZID",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.2,
                     minGSSize = 10,
                     maxGSSize = 500,
                     readable = TRUE)
  
  fwrite(as.data.frame(ora_go), file.path(opt$out_dir, paste0(cmp$name, "_ora_GO_results.csv")))  
  
  # KEGG pathway ORA
  ora_kegg <- enrichKEGG(gene = entrez_genes,
                         universe = universe_genes,
                         organism = "hsa",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.2,
                         minGSSize = 10,
                         maxGSSize = 500)
  
  fwrite(as.data.frame(ora_kegg), file.path(opt$out_dir, paste0(cmp$name, "_ora_KEGG_results.csv")))  
  
  # Reactome pathway ORA
  ora_reactome <- enrichPathway(gene = entrez_genes,
                                universe = universe_genes,
                                organism = "human",
                                pAdjustMethod = "BH",
                                pvalueCutoff = 0.05,
                                qvalueCutoff = 0.2,
                                minGSSize = 10,
                                maxGSSize = 500,
                                readable = TRUE)
  
  fwrite(as.data.frame(ora_reactome), file.path(opt$out_dir, paste0(cmp$name, "_ora_REACTOME_results.csv")))
  
  # Disease Ontology ORA
  ora_do <- enrichDO(gene = entrez_genes,
                     universe = universe_genes,
                     ont = "HDO",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.2,
                     minGSSize = 10,
                     maxGSSize = 500,
                     readable = TRUE)
  
  fwrite(as.data.frame(ora_do), file.path(opt$out_dir, paste0(cmp$name, "_ora_DO_results.csv")))
  
  # ORA plots
  ora_list <- list(GO = ora_go,
                   KEGG = ora_kegg,
                   REACTOME = ora_reactome,
                   DO = ora_do)
  
  for (db_type in names(ora_list)) {

    ora_res <- ora_list[[db_type]]
    
    # ORA dotplot
    dot_output <- file.path(opt$out_dir, paste0(cmp$name, "_ora_", db_type, "_dotplot.png"))
    tryCatch({
      p_dotplot <- plot_dotplot(ora_res, analysis = 'ora')      
      ggsave(filename = dot_output, plot = p_dotplot, width = 10, height = 6)
    }, error = function(e) {
      fail_file <- sub("\\.png$", ".fail", dot_output)
      file.create(fail_file)
      message("Failed to create ora dotplot for ", cmp$name, " ", db_type, ": ", e$message)
    })
    
    # ORA barplot
    bar_output <- file.path(opt$out_dir, paste0(cmp$name, "_ora_", db_type, "_barplot.png"))
    tryCatch({
      p_barplot <- plot_barplot(ora_res)      
      ggsave(filename = bar_output, plot = p_barplot, width = 10, height = 6)
    }, error = function(e) {
      fail_file <- sub("\\.png$", ".fail", bar_output)
      file.create(fail_file)
      message("Failed to create ora barplot for ", cmp$name, " ", db_type, ": ", e$message)
    })
    
    # ORA cnetplot
    cnetplot_output <- file.path(opt$out_dir, paste0(cmp$name, "_ora_", db_type, "_cnetplot.png"))
    tryCatch({
      p_cnetplot <- plot_cnetplot(ora_res)
      ggsave(filename = cnetplot_output, plot = p_cnetplot, width = 20, height = 20)
    }, error = function(e) {
      fail_file <- sub("\\.png$", ".fail", cnetplot_output)
      file.create(fail_file)
      message("Failed to create ora cnetplot for ", cmp$name, " ", db_type, ": ", e$message)
    })

  }
}

cat("ORA analysis completed\n")
