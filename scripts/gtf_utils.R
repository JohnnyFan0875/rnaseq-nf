build_transcript_gene_map_from_gtf <- function(gtf_file) {
  gtf_lines <- readLines(gtf_file, warn = FALSE)
  gtf_lines <- gtf_lines[!grepl("^#", gtf_lines)]

  if (length(gtf_lines) == 0) {
    stop("GTF file is empty or contains only comments: ", gtf_file)
  }

  fields <- strsplit(gtf_lines, "\t", fixed = TRUE)
  fields <- fields[lengths(fields) >= 9]

  if (length(fields) == 0) {
    stop("No valid feature rows found in GTF file: ", gtf_file)
  }

  attributes <- vapply(fields, `[`, character(1), 9)
  has_tx <- grepl('transcript_id "', attributes, fixed = TRUE)
  has_gene <- grepl('gene_id "', attributes, fixed = TRUE)
  attributes <- attributes[has_tx & has_gene]

  if (length(attributes) == 0) {
    stop("No transcript_id / gene_id pairs found in GTF file: ", gtf_file)
  }

  transcript_id <- sub('.*transcript_id "([^"]+)".*', "\\1", attributes)
  gene_id <- sub('.*gene_id "([^"]+)".*', "\\1", attributes)

  transcript_gene_map <- unique(data.frame(
    transcript_id = sub("\\..*$", "", transcript_id),
    gene_id = sub("\\..*$", "", gene_id),
    stringsAsFactors = FALSE
  ))

  transcript_gene_map <- transcript_gene_map[
    nzchar(transcript_gene_map$transcript_id) & nzchar(transcript_gene_map$gene_id),
  ]

  if (nrow(transcript_gene_map) == 0) {
    stop("Failed to derive transcript-to-gene mapping from GTF file: ", gtf_file)
  }

  transcript_gene_map
}
