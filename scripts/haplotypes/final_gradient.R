
########################################
# Script to look at detailed EVE abundance between main and haplotype accessions
########################################
# created by Nadja Brait

library(dplyr)

# you will need list of accession pairs and similarity matrices of each accession pair 
# similarity matrices (csv) were generated in Geneious with a multiple alignment per pair

# paired list with accessions of main genomes and alternative haplotypes
pairings_df <- data.frame(
  main_genome = c("JABUIU01"
                  ,"GCF_943734705.1"
                  ,"GCF_943734745.1"
                  ,"GCF_943734845.2"
                  ,"GCA_964033645.1"
                  ,"GCF_943734735.2"
                  ,"GCF_943734695.1"
                  ,"GCF_943734725.1"
                  ,"GCF_943734755.1"
                  ,"GCF_943737925.1"
                  ,"GCF_943734765.1"
                  ,"GCF_943734685.1"
                  ,"CAJZCY01"
                  ,"GCA_943734665.2"
                  ,"GCA_963924065.2"),
  haplotype = c("JABUIV01"
                ,"CALSDW01"
                ,"CALSEC01"
                ,"CALSDO01"
                ,"CAWUOS01"
                ,"CALSDS01"
                ,"CALSEK01"
                ,"CALSEI01"
                ,"CALSED01"
                ,"CALSFZ01"
                ,"CALSEF01"
                ,"CALSEE01"
                ,"CALTRM01"
                ,"CALSEH02"
                ,"CAWUOX01")
)

all_detailed_results <- data.frame()
all_missing_results <- data.frame()


# this loop takes a long time (hours)
# this loop shows detailed differences of EVE abundance between main and haplotye genomes 
# each main accession is matched with the best haplotype match
# if a pair is formed they are no longer used for follow up matches
# loops until no matches are found between main and haplotype accessions anymore
# matches must be at least 95% identity
# unmatched EVEs are then put into a separate df
# a summary df with all EVEs (matched or unmatched) is created 

for (i in 1:nrow(pairings_df)) {
  main_genome <- pairings_df$main_genome[i]
  haplotype <- pairings_df$haplotype[i]
  filename <- paste0(main_genome, "_", haplotype, ".csv") # extract filenames for csv files
  similarity_matrix <- read.csv(filename, row.names = 1)
  
  # Set row and column names for the similarity matrix
  #rownames(similarity_matrix) <- rownames(similarity_matrix)
  #colnames(similarity_matrix) <- colnames(similarity_matrix)
  
  # keep only main genome (rows) and haplotype (columns)
  main_genome_rows <- grep(main_genome, rownames(similarity_matrix))
  haplotype_cols <- grep(haplotype, colnames(similarity_matrix))
  filtered_matrix <- similarity_matrix[main_genome_rows, haplotype_cols]
  # Skip if filtered_matrix is empty
  if (nrow(filtered_matrix) == 0 | ncol(filtered_matrix) == 0) {
    next
  }
  # df for summary
  detailed_df <- data.frame(
    Main_genome = character(),
    haplotype = character(),
    similarity_pair = numeric(),
    best_match = character(),
    similarity_best_match = numeric(),
    stringsAsFactors = FALSE
  )
  
  # df to track unmatched EVEs
  missing_df <- data.frame(
    Missing_EVEs = character(),
    similarity_best_match = numeric(),
    stringsAsFactors = FALSE
  )
  # rows and columns are only allowed to be used once
  used_rows <- character()
  used_cols <- character()
  
# potential_matches includes all potential matches (above 95%) before assigning them
  potential_matches <- data.frame(
    Main_genome = character(),
    haplotype = character(),
    similarity = numeric(),
    stringsAsFactors = FALSE
  )
# loop to get potential matches  
  for (j in 1:nrow(filtered_matrix)) {
    row_name <- rownames(filtered_matrix)[j]
    # Skip rows already used
    if (row_name %in% used_rows) {
      next
    }
    for (k in 1:ncol(filtered_matrix)) {
      col_name <- colnames(filtered_matrix)[k]
      # Skip columns already used
      if (col_name %in% used_cols) {
        next
      }
      # put potential matches into df
      similarity_value <- as.numeric(filtered_matrix[j, k])
      if (!is.na(similarity_value) && similarity_value >= 95) {
        potential_matches <- rbind(
          potential_matches,
          data.frame(
            Main_genome = row_name,
            haplotype = col_name,
            similarity = similarity_value,
            stringsAsFactors = FALSE
          )
        )
      }
    }
  }
  
  potential_matches <- potential_matches %>%
    arrange(desc(similarity))
  
  # Assign best matches in order, but skip rows and columns hat have  already been used
  for (i in 1:nrow(potential_matches)) {
    row_name <- potential_matches$Main_genome[i]
    col_name <- potential_matches$haplotype[i]
    similarity_value <- potential_matches$similarity[i]
    
    if (!(row_name %in% used_rows) && !(col_name %in% used_cols)) {
      detailed_df <- rbind(
        detailed_df,
        data.frame(
          Main_genome = row_name,
          haplotype = col_name,
          similarity_pair = similarity_value,
          best_match = col_name,
          similarity_best_match = similarity_value,
          stringsAsFactors = FALSE
        )
      )
      used_rows <- c(used_rows, row_name)
      used_cols <- c(used_cols, col_name)
    }
  }
  
  # need to add unmatched main genome accessions to summary df 
  for (i in rownames(filtered_matrix)) {
    if (!(i %in% used_rows)) {
      best_match_row <- colnames(filtered_matrix)[which.max(as.numeric(filtered_matrix[i, ]))]
      best_value <- max(as.numeric(filtered_matrix[i, ]), na.rm = TRUE)
      
      detailed_df <- rbind(
        detailed_df,
        data.frame(
          Main_genome = i,
          haplotype = NA,
          similarity_pair = NA,
          best_match = best_match_row,
          similarity_best_match = best_value,
          stringsAsFactors = FALSE
        )
      )
    }
  }
  
  # find unmatched haplotypes and add to summary df
  for (k in 1:ncol(filtered_matrix)) {
    col_name <- colnames(filtered_matrix)[k]
    if (col_name %in% used_cols) {
      next
    }
    row_best_match <- rownames(filtered_matrix)[which.max(as.numeric(filtered_matrix[, k]))]
    best_value <- max(as.numeric(filtered_matrix[, k]), na.rm = TRUE);
    detailed_df <- rbind(
      detailed_df,
      data.frame(
        Main_genome = NA,
        haplotype = col_name,
        similarity_pair = NA,
        best_match = row_best_match,
        similarity_best_match = best_value,
        stringsAsFactors = FALSE
      )
    )
    used_cols <- c(used_cols, col_name)
  }
  
  # missing_main <- detailed_df %>%
  #   filter(is.na(haplotype)) %>%
  #   select(Main_genome, similarity_best_match)
  # 
  # missing_haplotype <- detailed_df %>%
  #   filter(is.na(Main_genome)) %>%
  #   select(haplotype, similarity_best_match)
  # 
  # missing_haplotype <- missing_haplotype %>%
  #   rename(Main_genome = haplotype)
  # 
  # missing_df <- bind_rows(missing_main, missing_haplotype) %>%
  #   rename(Missing_EVEs = Main_genome)
  
  # unmatched main accessions
  missing_main <- detailed_df %>%
    filter(is.na(haplotype)) %>%
    select(Main_genome, similarity_best_match) %>%
    mutate(Source = "Main_genome")  
  
  # unmatched haplotype accessions
  missing_haplotype <- detailed_df %>%
    filter(is.na(Main_genome)) %>%
    select(haplotype, similarity_best_match) %>%
    rename(Main_genome = haplotype) %>%
    mutate(Source = "Haplotype")  
  
  # Combine all unmatched into missing_df
  missing_df <- bind_rows(missing_main, missing_haplotype) %>%
    rename(Missing_EVEs = Main_genome)
  
  # final output dfs
  all_detailed_results <- bind_rows(all_detailed_results, detailed_df)
  all_missing_results <- bind_rows(all_missing_results, missing_df)
}

write.table(all_detailed_results, file = "similarity_gradient_output.tsv", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE )
#write_tsv(all_detailed_results, file = "all_detailed_unmatched_EVEs.tsv")
#################

# summaries for main text regarding unmatched EVEs (also 95% matches)

main_genome_missing <- all_missing_results %>%
  filter(Source == "Main_genome")

haplotype_missing <- all_missing_results %>%
  filter(Source == "Haplotype")

main_over_95 <- main_genome_missing %>%
  filter(similarity_best_match >= 95) %>%
  nrow()

main_below_95 <- main_genome_missing %>%
  filter(similarity_best_match < 95) %>%
  nrow()

haplotype_over_95 <- haplotype_missing %>%
  filter(similarity_best_match >= 95) %>%
  nrow()

haplotype_below_95 <- haplotype_missing %>%
  filter(similarity_best_match < 95) %>%
  nrow()

cat("Main_genome missing EVEs with similarity >= 95%:", main_over_95, "\n")
cat("Main_genome missing EVEs with similarity < 95%:", main_below_95, "\n")

cat("Haplotype missing EVEs with similarity >= 95%:", haplotype_over_95, "\n")
cat("Haplotype missing EVEs with similarity < 95%:", haplotype_below_95, "\n")

####################

# missing EVEs per genome

main_genome_missing <- main_genome_missing %>%
  mutate(sample = sub("_[^_]+$", "", Missing_EVEs)) %>%  
  group_by(sample) %>%
  summarize(count = n()) %>%
  mutate(Source = "main_genome") 

haplotype_missing <- haplotype_missing %>%
  mutate(sample = sub("_[^_]+$", "", Missing_EVEs)) %>%  
  group_by(sample) %>%
  summarize(count = n()) %>%
  mutate(Source = "haplotype")  

combined_missing <- bind_rows(main_genome_missing, haplotype_missing)

##################

# Contig detection per missing EVE

main_genome_missing <- all_missing_results %>%
  filter(Source == "Main_genome")

haplotype_missing <- all_missing_results %>%
  filter(Source == "Haplotype")

combined_missing <- bind_rows(main_genome_missing, haplotype_missing)

# combined_missing <- combined_missing %>%
#   filter(similarity_best_match >= 95)
# 
# rest <- combined_missing %>%
#   filter(similarity_best_match < 95)

combined_missing <- combined_missing %>%
  rename(eve_id = Missing_EVEs) %>%  
  mutate(eve_id = gsub(" \\(reversed\\)", "", eve_id)) %>%
  mutate(eve_id = gsub("\\.\\.|reversed\\.", "", eve_id))

merged_results <- inner_join(combined_missing, metadata, by = "eve_id")  

summarized_counts <- merged_results %>%
  group_by(Contig,sample, Source) %>%  
  summarise(count = n(), .groups = 'drop')  

summarized_counts %>%
  summarize(
    average_count = mean(count, na.rm = TRUE),   
    sd_count = sd(count, na.rm = TRUE),           
    max_count = max(count, na.rm = TRUE),         
    min_count = min(count, na.rm = TRUE)          
  )

#write.table(summarized_counts, file = "missing_EVEs_per_contig.tsv", sep= "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)

###############

# visualization

library(ggplot2)
library(dplyr)

# contigs per sample
contig_counts <- summarized_counts %>%
  group_by(sample) %>%
  summarise(num_contigs = n(), .groups = 'drop')  

max_counts <- summarized_counts %>%
  group_by(sample) %>%
  summarize(max_count = max(count), .groups = 'drop')

combined_plot <- ggplot() +
  geom_bar(data = contig_counts, aes(x = sample, y = log1p(num_contigs)), 
           stat = "identity", fill = "steelblue", alpha = 0.4) +  
  geom_violin(data = summarized_counts, aes(x = sample, y = count, fill = Source), 
              trim = FALSE, alpha = 0.6, adjust = 0.5) +
              scale_fill_manual(values = c("forestgreen", "orange")) +
  geom_boxplot(data = summarized_counts, aes(x = sample, y = count, fill = Source), 
               width = 0.1, color = "black", outlier.shape = NA, alpha = 0.5) +  
  geom_point(data = max_counts, aes(x = sample, y = max_count), 
             color = "black", size = 2, shape = 21, fill = "red") +  
  geom_text(data = contig_counts, 
            aes(x = sample, y = -0.5, label = num_contigs),  
            vjust = 2, size = 3, color = "red") +  
  labs(title = "Distribution of EVEs Across Samples",
       x = "Sample",
       y = "Counts") +
  #scale_y_continuous(limits = c(0, max(summarized_counts$count) * 1.1)) +  
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

print(combined_plot)


