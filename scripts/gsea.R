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
cat('Starting gsea analysis\n')

source(file.path(opt$script_dir, "plot_utils_enrichment.R"))

# Read metadata
metadata <- read_tsv(opt$metadata)
rownames(metadata) <- metadata$sample_id

# Parse comparisons JSON string
comparisons <- fromJSON(opt$comparisons) %>% split(seq(nrow(.)))

# for testing in Rstudio
# yaml_data <- yaml::read_yaml("D:\\shared\\rna_seq_pipeline\\data\\project_1\\config\\project_config.yaml")
# json_data <- jsonlite::toJSON(yaml_data, pretty = TRUE, auto_unbox = TRUE)
# opt <- list(
#   deg_dir = "D:\\shared\\rna_seq_pipeline\\data\\project_1\\results\\de_analysis\\DEG_edgeR",
#   metadata = "D:\\shared\\rna_seq_pipeline\\data\\project_1\\metadata\\sample_metadata.tsv",
#   out_dir = "D:\\shared\\rna_seq_pipeline\\data\\project_1\\results\\gsea_results",
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

set.seed(123)

# Run GSEA
for  (cmp in comparisons) {
  
  deg_file <- file.path(opt$deg_dir, paste0(cmp$name, "_", opt$de_method, "_results.csv"))
  deg_data <- read.csv(deg_file)

  ranked_gene_list <- deg_data %>%
    filter(!is.na(logFC), !is.na(PValue), !is.na(gene_symbol)) %>%
    filter(logCPM > 1) %>%
    mutate(rank_metric = -log10(PValue) * sign(logFC)) %>%
    { setNames(.$rank_metric, .$gene_symbol) } %>%
    sort(decreasing = TRUE) %>%
    { .[!duplicated(names(.))] }

  gene_df <- bitr(names(ranked_gene_list),
                  fromType = "SYMBOL", 
                  toType = "ENTREZID", 
                  OrgDb = org.Hs.eg.db)

  ranked_gene_list <- ranked_gene_list[gene_df$SYMBOL]
  names(ranked_gene_list) <- gene_df$ENTREZID
  ranked_gene_list <- ranked_gene_list[!is.na(names(ranked_gene_list))]
  ranked_gene_list <- ranked_gene_list[!duplicated(names(ranked_gene_list))]

  # Run GSEA for GO Biological Process
  gsea_go <- gseGO(geneList = ranked_gene_list,
                   OrgDb = org.Hs.eg.db,
                   ont = "BP",
                   keyType = "ENTREZID",
                   minGSSize = 10,
                   maxGSSize = 500,
                   pvalueCutoff = 0.05,
                   pAdjustMethod = "BH",
                   nPerm = 1000,
                   verbose = FALSE)

  fwrite(as.data.frame(gsea_go), file.path(opt$out_dir, paste0(cmp$name, "_gsea_GO_results.csv")))
  
  # Run GSEA for KEGG
  gsea_kegg <- gseKEGG(geneList = ranked_gene_list,
                       organism = "hsa",
                       minGSSize = 10,
                       maxGSSize = 500,
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "BH",
                       nPerm = 1000,
                       verbose = FALSE)

  fwrite(as.data.frame(gsea_kegg), file.path(opt$out_dir, paste0(cmp$name, "_gsea_KEGG_results.csv")))
  
  # Run GSEA for Reactome
  gsea_reactome <- gsePathway(geneList = ranked_gene_list,
                              organism = "human",
                              minGSSize = 10,
                              maxGSSize = 500,
                              pvalueCutoff = 0.05,
                              pAdjustMethod = "BH",
                              nPerm = 1000,
                              verbose = FALSE)

  fwrite(as.data.frame(gsea_reactome), file.path(opt$out_dir, paste0(cmp$name, "_gsea_REACTOME_results.csv")))  

  # Run GSEA for Disease Ontology (DO)
  gsea_do <- gseDO(geneList = ranked_gene_list,
                   minGSSize = 10,
                   maxGSSize = 500,
                   pvalueCutoff = 0.05,
                   pAdjustMethod = "BH",
                   nPerm = 1000,
                   verbose = FALSE)

  # Save results
  fwrite(as.data.frame(gsea_do), file.path(opt$out_dir, paste0(cmp$name, "_gsea_DO_results.csv")))

  # GSEA plots
  gsea_list <- list(GO = gsea_go,
                    KEGG = gsea_kegg,
                    REACTOME = gsea_reactome,
                    DO = gsea_do)
  
  for (db_type in names(gsea_list)) {

    gsea_res <- gsea_list[[db_type]]
    
    # GSEA dotplot
    dot_output <- file.path(opt$out_dir, paste0(cmp$name, "_gsea_", db_type, "_dotplot.png"))
    tryCatch({
      p_dotplot <- plot_dotplot(gsea_res, analysis = 'gsea')      
      ggsave(filename = dot_output, plot = p_dotplot, width = 10, height = 6)
    }, error = function(e) {
      fail_file <- sub("\\.png$", ".fail", dot_output)
      file.create(fail_file)
      message("Failed to create gsea dotplot for ", cmp$name, " ", db_type, ": ", e$message)
    })

    # GSEA gseaplot
    gseaplot_output <- file.path(opt$out_dir, paste0(cmp$name, "_gsea_", db_type, "_gseaplot.png"))
    tryCatch({
      p_gseaplot <- plot_gseaplot(gsea_res, geneSetID = 1)
      ggsave(filename = gseaplot_output, plot = p_gseaplot, width = 6, height = 6)
    }, error = function(e) {
      fail_file <- sub("\\.png$", ".fail", gseaplot_output)
      file.create(fail_file)
      message("Failed to create gseaplot for ", cmp$name, " ", db_type, ": ", e$message)
    })

    # GSEA cnetplot
    cnetplot_output <- file.path(opt$out_dir, paste0(cmp$name, "_gsea_", db_type, "_cnetplot.png"))
    tryCatch({
      p_cnetplot <- plot_cnetplot(gsea_res)
      ggsave(filename = cnetplot_output, plot = p_cnetplot, width = 12, height = 12)
    }, error = function(e) {
      fail_file <- sub("\\.png$", ".fail", cnetplot_output)
      file.create(fail_file)
      message("Failed to create gsea cnetplot for ", cmp$name, " ", db_type, ": ", e$message)
    })

    # GSEA ridgeplot
    ridgeplot_output <- file.path(opt$out_dir, paste0(cmp$name, "_gsea_", db_type, "_ridgeplot.png"))
    tryCatch({
      p_ridgeplot <- plot_ridgeplot(gsea_res)
      ggsave(filename = ridgeplot_output, plot = p_ridgeplot, width = 12, height = 10)
    }, error = function(e) {
      fail_file <- sub("\\.png$", ".fail", ridgeplot_output)
      file.create(fail_file)
      message("Failed to create gsea ridgeplot for ", cmp$name, " ", db_type, ": ", e$message)
    })

    # GSEA emap
    emap_output <- file.path(opt$out_dir, paste0(cmp$name, "_gsea_", db_type, "_emap.png"))
    tryCatch({
      p_emap <- plot_emap(gsea_res)  
      ggsave(filename = emap_output, plot = p_emap, width = 12, height = 10)
    }, error = function(e) {
      fail_file <- sub("\\.png$", ".fail", emap_output)
      file.create(fail_file)
      message("Failed to create gsea enrichment map for ", cmp$name, " ", db_type, ": ", e$message)
    })
    
    # GSEA treeplot
    treeplot_output <- file.path(opt$out_dir, paste0(cmp$name, "_gsea_", db_type, "_treeplot.png"))
    tryCatch({
      p_treeplot <- plot_treeplot(gsea_res)
      ggsave(filename = treeplot_output, plot = p_treeplot, width = 15, height = 10)
    }, error = function(e) {
      fail_file <- sub("\\.png$", ".fail", treeplot_output)
      file.create(fail_file)
      message("Failed to create gsea treeplot for ", cmp$name, " ", db_type, ": ", e$message)
    })
    
    # GSEA pathview (only for KEGG database)
    if (db_type == 'KEGG') {
      
      top_pathway <- gsea_res %>%
        filter(p.adjust < 0.05) %>%
        arrange(p.adjust) %>%
        head(5) 
      
      for (pathway_id in top_pathway$ID) {
        pathway_id <- sub("hsa", "", pathway_id)
        pathview_output <- file.path(opt$out_dir, paste0(cmp$name, "_gsea_", db_type, "_pathview_", pathway_id, ".png"))
        
        tryCatch({
          pathview_dir <- plot_pathview(ranked_gene_list, pathway_id, cmp_name = cmp$name, output_dir = opt$out_dir)
          file.copy(file.path(pathview_dir, paste0('hsa', pathway_id, '.colored.png')), pathview_output)
        }, error = function(e) {
          fail_file <- sub("\\.png$", ".fail", pathview_output)
          file.create(fail_file)
          message("Failed to create gsea pathview ", pathway_id, " for ", cmp$name, " ", db_type, ": ", e$message)
        })
      }
    }
  }
}

cat("GSEA analysis completed\n")
