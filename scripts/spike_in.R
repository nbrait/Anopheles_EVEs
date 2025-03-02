########################################
# SCRIPT to spike in Anopheles specific viruses if not available in RVDB
########################################
# created by Nadja Brait 

library(dplyr)
library(tidyr)
library(readr)
library(stringr)

# Anopheles specific viruses from the Anopheles virus review
exogenous_viruses <- read_tsv("exogenous_viruses.tsv")

# RVDB input
data <- readLines("rvdb80_headers.txt")  
# extracting virus names from headers
references <- str_extract_all(data, "\\[([^]]+)\\]")
references_flat <- unlist(references)
references_flat <- gsub("\\[|\\]", "", references_flat)
unique_references <- unique(references_flat)
unique_references_df <- data.frame(virus_name = unique_references, stringsAsFactors = FALSE)

#matched_viruses <- merge(exogenous_viruses, unique_references_df, by = "virus_name")
# viruses for spike-in
non_matched_viruses <- anti_join(exogenous_viruses, unique_references_df, by = "virus_name")

###############
# download missing viruses for spike-in
###############

library(rentrez)

# function for taxid
get_taxid <- function(virus_name) {
  search_result <- entrez_search(db = "taxonomy", term = virus_name)
  #print(search_result)
  if (!is.null(search_result$ids) && length(search_result$ids) > 0) {
    taxid <- search_result$ids[1]  
    return(taxid)  
  } else {
    cat("No IDs found for:", virus_name, "\n")
    return(NA)
  }
}

# function for getting protein sequences from taxid
get_protein_sequences_by_taxid <- function(taxid) {
  search_result_noexp <- entrez_search(db = "protein", term = paste0("txid", taxid, "[Organism:noexp]"))
  if (search_result_noexp$count > 0) {
    protein_ids <- search_result_noexp$ids
    protein_sequences <- entrez_fetch(db = "protein", id = protein_ids, rettype = "fasta", retmode = "text")
    return(protein_sequences)
  }
  search_result_exp <- entrez_search(db = "protein", term = paste0("txid", taxid, "[Organism:exp]"))
  
  if (search_result_exp$count > 0) {
    protein_ids <- search_result_exp$ids
    protein_sequences <- entrez_fetch(db = "protein", id = protein_ids, rettype = "fasta", retmode = "text")
    return(protein_sequences)
  }
  
  cat("No protein sequences found for TaxID:", taxid, "\n")
  return(NULL)
}

all_protein_sequences <- list()

# get accessions
for (virus in non_matched_viruses$virus_name) {
  cat("Processing:", virus, "\n")
  taxid <- get_taxid(virus)
  if (!is.na(taxid)) {
    protein_sequences <- get_protein_sequences_by_taxid(taxid)
    if (!is.null(protein_sequences)) {
      all_protein_sequences[[virus]] <- protein_sequences  
    }
  } else {
    cat("No TaxID found for virus:", virus, "\n")
  }
}

fasta_output <- unlist(all_protein_sequences)
#writeLines(fasta_output, con = "spiked_sequences.fasta") 
