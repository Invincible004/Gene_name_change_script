# Gene_name_change_script
The R script converts a GMT file (gene symbols) to a new GMT file with Entrez IDs. It requires 'gene_info' and 'gmt_file' as inputs, creates a mapping, replaces symbols with IDs, and writes output. Run as 'Rscript script.R gene_info gmt_file'.

# Here is the full command of the convert_gmt_to_entrez.R file which I have provided:
{
#!/usr/bin/env Rscript

# Load required libraries
library(R.utils)

# Function to create the gene symbol to Entrez ID mapping
create_gene_symbol_map <- function(gene_info_file) {
  gene_symbol_map <- list()

  # Read gene_info file
  gene_info <- read.delim(gene_info_file, header = FALSE, comment.char = "#", colClasses = "character")
  
  # Extract the required columns: GeneID, Symbol, Synonyms
  gene_info_selected <- gene_info[, c(2, 3, 5)]
  
  # Process each row and create the mapping
  for (i in 1:nrow(gene_info_selected)) {
    gene_id <- gene_info_selected[i, 1]
    symbol <- gene_info_selected[i, 2]
    synonyms <- strsplit(gene_info_selected[i, 3], "\\|")[[1]]
    
    # Add the main symbol and synonyms to the mapping
    gene_symbol_map[symbol] <- gene_id
    for (syn in synonyms) {
      gene_symbol_map[syn] <- gene_id
    }
  }

  return(gene_symbol_map)
}

# Function to replace gene symbols with Entrez IDs in GMT file
replace_symbols_with_entrez_ids <- function(gmt_file, gene_symbol_map, output_file) {
  output_lines <- character()

  # Read GMT file line by line
  con <- file(gmt_file, "rt")
  while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
    fields <- strsplit(line, "\t")[[1]]
    pathway_name <- fields[1]
    pathway_desc <- fields[2]
    genes <- fields[-c(1, 2)]
    
    # Replace gene symbols with Entrez IDs
    entrez_ids <- sapply(genes, function(gene) {
      if (gene %in% names(gene_symbol_map)) {
        return(gene_symbol_map[[gene]])
      } else {
        return(gene)
      }
    })
    
    # Construct the new GMT line with Entrez IDs
    new_line <- paste(c(pathway_name, pathway_desc, entrez_ids), collapse = "\t")
    output_lines <- c(output_lines, new_line)
  }
  close(con)
  
  # Write the new GMT file with Entrez IDs
  writeLines(output_lines, output_file)
}

# Main script
if (commandArgs(TRUE) < 2) {
  stop("Usage: Rscript convert_gmt_to_entrez.R Homo_sapiens.gene_info.gz h.all.v2023.1.Hs.symbols.gmt")
}

# Input files
gene_info_file <- commandArgs(TRUE)[1]
gmt_file <- commandArgs(TRUE)[2]

# Output file
output_file <- "output_entrez_ids.gmt"

# Create gene symbol to Entrez ID mapping
gene_symbol_map <- create_gene_symbol_map(gene_info_file)

# Replace gene symbols with Entrez IDs in GMT file and write the output
replace_symbols_with_entrez_ids(gmt_file, gene_symbol_map, output_file)

# Print completion message
cat("Gene symbols replaced with Entrez IDs and output written to", output_file, "\n")
}

# To run the script from the command line put the comand:
./convert_gmt_to_entrez.R Homo_sapiens.gene_info.gz h.all.v2023.1.Hs.symbols.gmt
