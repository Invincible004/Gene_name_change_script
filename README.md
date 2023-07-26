# Gene_name_change_script
The R script program takes GMP file as a input and replace the gene name with the EntrezID from the file where data is available.

Rscript -e 'library(readr)

parse_gene_info <- function(gene_info_file) {
  gene_info <- read_tsv(gene_info_file, col_types = "cccc", na = "")
  symbol_to_entrez <- list()

  for (i in 1:nrow(gene_info)) {
    entrez_id <- gene_info[i, "GeneID"]
    symbol <- gene_info[i, "Symbol"]
    synonyms <- unlist(strsplit(gene_info[i, "Synonyms"], "\\\\|"))

    symbol_to_entrez[[symbol]] <- entrez_id
    for (synonym in synonyms) {
      symbol_to_entrez[[synonym]] <- entrez_id
    }
  }

  return(symbol_to_entrez)
}

replace_gene_names <- function(gmt_file, symbol_to_entrez) {
  gmt_data <- read_tsv(gmt_file, col_names = FALSE, col_types = cols(.default = "c"), na = "")
  new_gmt_data <- list()

  for (i in 1:nrow(gmt_data)) {
    pathway_name <- gmt_data[i, 1]
    pathway_description <- gmt_data[i, 2]
    genes <- gmt_data[i, -(1:2)]

    new_genes <- sapply(genes, function(gene) {
      if (gene %in% names(symbol_to_entrez)) {
        return(symbol_to_entrez[[gene]])
      } else {
        return(gene)
      }
    })

    new_row <- c(pathway_name, pathway_description, new_genes)
    new_gmt_data[[i]] <- new_row
  }

  return(do.call(rbind, new_gmt_data))
}

write_new_gmt <- function(new_gmt_data, output_file) {
  write_tsv(new_gmt_data, output_file, col_names = FALSE)
}

gene_info_file <- "Homo_sapiens.gene_info.gz"
gmt_file <- "h.all.v2023.1.Hs.symbols.gmt"
output_file <- "output.gmt"

symbol_to_entrez <- parse_gene_info(gene_info_file)
new_gmt_data <- replace_gene_names(gmt_file, symbol_to_entrez)
write_new_gmt(new_gmt_data, output_file)'
