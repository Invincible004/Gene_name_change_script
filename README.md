# Gene_name_change_script
The R script program takes GMP file as a input and replace the gene name with the EntrezID from the file where data is available.

# Load required libraries
library(data.table)
library(stringr)

# Function to create a mapping of Symbol to GeneID
create_gene_id_mapping <- function(gene_info_file) {
  symbol_to_gene_id <- list()
  
  chunk_size <- 10000  # Adjust the chunk size as per your system's memory capacity
  gene_info_dt <- fread(gene_info_file, sep = "\t", nrows = chunk_size, colClasses = c(GeneID = "character", Symbol = "character", Synonyms = "character"))
  
  while (!is.null(gene_info_dt)) {
    for (i in seq_len(nrow(gene_info_dt))) {
      symbols <- c(gene_info_dt$Symbol[i], unlist(strsplit(gene_info_dt$Synonyms[i], "\\|")))
      gene_id <- gene_info_dt$GeneID[i]
      
      for (symbol in symbols) {
        symbol_to_gene_id[[symbol]] <- gene_id
      }
    }
    
    gene_info_dt <- fread(gene_info_file, sep = "\t", nrows = chunk_size, colClasses = c(GeneID = "character", Symbol = "character", Synonyms = "character"))
  }
  
  return(symbol_to_gene_id)
}

# Function to replace symbols with Entrez IDs in GMT file
replace_symbols_with_entrez_ids <- function(gmt_file, gene_id_mapping) {
  gmt_lines <- readLines(gmt_file)
  new_gmt_lines <- character(length(gmt_lines))
  
  for (i in seq_len(length(gmt_lines))) {
    gmt_line <- gmt_lines[i]
    gene_symbols <- str_split(gmt_line, "\t")[[1]][-(1:2)]  # Skip first two fields (pathway name and description)
    gene_ids <- sapply(gene_symbols, function(symbol) {
      symbols <- str_trim(symbol)
      gene_id <- gene_id_mapping[[symbols]]
      if (!is.null(gene_id)) {
        return(gene_id)
      } else {
        return(NA)
      }
    })
    new_gmt_line <- paste(gmt_line[1:2], paste(gene_ids, collapse = "|"), sep = "\t")
    new_gmt_lines[i] <- new_gmt_line
  }
  
  return(new_gmt_lines)
}

# Input files
gene_info_file <- "Homo_sapiens.gene_info.gz"
gmt_file <- "h.all.v2023.1.Hs.symbols.gmt"

# Step 1: Create Symbol to GeneID mapping
symbol_to_gene_id_mapping <- create_gene_id_mapping(gene_info_file)

# Step 2: Replace symbols with Entrez IDs in the GMT file
new_gmt_lines <- replace_symbols_with_entrez_ids(gmt_file, symbol_to_gene_id_mapping)

# Output the new GMT file to the terminal
cat(new_gmt_lines, sep = "\n")
