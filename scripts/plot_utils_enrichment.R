library(enrichplot)
library(ggplot2)
library(DOSE)
library(GOSemSim)
library(pathview)

plot_dotplot <- function(enrichment_result, title = "Enriched Pathways", showCategory = 20, analysis = 'gsea') {
  
  if (is.null(enrichment_result) || nrow(as.data.frame(enrichment_result)) == 0) {
    stop("Enrichment result object is empty or invalid.")
  }
  
  if (analysis == 'gsea') {
    p <- dotplot(enrichment_result, showCategory = 10, split = ".sign", label_format = 50, title = title) +
      facet_grid(. ~ .sign) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5))
  } else if (analysis == 'ora') {
    p <- dotplot(enrichment_result, showCategory = 10, label_format = 50, title = title) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5))
  }
  
  
  return(p)
}

plot_gseaplot <- function(enrichment_result, geneSetID = 1) {
  
  if (is.null(enrichment_result) || nrow(as.data.frame(enrichment_result)) == 0) {
    stop("Enrichment result object is empty or invalid.")
  }
  
  p <- gseaplot2(enrichment_result, 
                 geneSetID = geneSetID, 
                 title = enrichment_result@result$Description[geneSetID], 
                 color = "red")
  
  return(p)
}

plot_cnetplot <- function(enrichment_result, 
                          key_type = "ENTREZID", 
                          OrgDb = org.Hs.eg.db, 
                          showCategory = 5, 
                          title = "Gene-Concept Network") {

  if (is.null(enrichment_result) || nrow(as.data.frame(enrichment_result)) == 0) {
    stop("Enrichment result object is empty or invalid.")
  }
  
  readable <- DOSE::setReadable(enrichment_result, OrgDb = OrgDb, keyType = key_type)
  
  p <- cnetplot(readable, 
                categorySize = "pvalue", 
                showCategory = showCategory, 
                node_label = "all",
                foldChange = NULL) +
       theme_bw() +
       ggtitle(title) +
       theme(panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             plot.title = element_text(hjust = 0.5))
  
  return(p)
}

plot_ridgeplot <- function(enrichment_result, showCategory = 20, title = "Ridgeplot") {
  
  if (is.null(enrichment_result) || nrow(as.data.frame(enrichment_result)) == 0) {
    stop("Enrichment result object is empty or invalid.")
  }
  
  p <- ridgeplot(enrichment_result, showCategory = showCategory, fill = "p.adjust") +
    ggtitle(title) +
    xlab("Enrichment Score") + 
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  
  return(p)
}

plot_emap <- function(enrichment_result) {
  
  if (is.null(enrichment_result) || nrow(as.data.frame(enrichment_result)) == 0) {
    stop("Enrichment result object is empty or invalid.")
  }
  
  emap <- pairwise_termsim(enrichment_result) 

  p <- emapplot(emap, color = "p.adjust") + 
    ggtitle("Enrichment Map") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  
  return(p)
}

plot_treeplot <- function(enrichment_result, showCategory = 10, title = "Treeplot") {
  
  if (is.null(enrichment_result) || nrow(as.data.frame(enrichment_result)) == 0) {
    stop("Enrichment result object is empty or invalid.")
  }
  
  # Calculate pairwise term similarity for treeplot
  emap <- pairwise_termsim(enrichment_result)
  
  # Generate treeplot
  p <- treeplot(emap, showCategory = showCategory, color = "pvalue") +
    ggtitle(title) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  
  return(p)
}

plot_pathview <- function(ranked_gene_list,
                          pathway_id,
                          cmp_name = NULL,
                          output_dir = NULL
                          ) {
  
  pathview_dir = paste0(output_dir, cmp_name, "_KEGG_pathview")
  dir.create(pathview_dir, recursive = TRUE)

  pathview(gene.data = ranked_gene_list,
           pathway.id = pathway_id,
           species = 'hsa',
           kegg.native = TRUE,
           out.suffix = 'colored'
           )

  files_to_copy <- list.files(path = ".", pattern = "^hsa", full.names = TRUE)
  file.copy(from = files_to_copy, to = pathview_dir, overwrite = TRUE, recursive = FALSE)

  return(pathview_dir)
}

plot_barplot <- function(enrichment_result, showCategory = 20, title = "Barplot") {
  
  if (is.null(enrichment_result) || nrow(as.data.frame(enrichment_result)) == 0) {
    stop("Enrichment result object is empty or invalid.")
  }
  
  p <- barplot(enrichment_result, showCategory = showCategory) +
    ggtitle(title) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  
  return(p)
}
