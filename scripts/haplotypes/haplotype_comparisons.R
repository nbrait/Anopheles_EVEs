#######################
# SCRIPT FOR HAPLOTYPE COMPARISONS - EVE counts and genome base pairs
#######################
# by Nadja Brait

setwd("C:/Users/nadja/Documents/LaptopAsus/PhD/Chapter_4/")

library(ggplot2)
library(gridExtra)
library(data.table)
library(patchwork)
library(tidyr)
library(dplyr)

accession_pairs <- read.csv("haplotype_pairs.csv", header = FALSE, stringsAsFactors = FALSE)
metadata <- read.table("input_haplotypes.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
metadata$sample <- as.character(metadata$sample)
viral_families <- unique(metadata$top_viral_family)

final_table <- data.frame(
  sample = character(),
  total_count = integer(),
  matrix(ncol = length(viral_families), nrow = 0)
)

colnames(final_table) <- c("sample", "total_count", viral_families)

calculate_counts <- function(accession) {
  accession_data <- metadata[metadata$sample == accession, ]
  total_eve_count <- nrow(accession_data)
  family_counts <- integer(length(viral_families))
  names(family_counts) <- viral_families
  if (total_eve_count > 0) {
    family_table <- table(accession_data$top_viral_family)
    family_counts[names(family_table)] <- as.numeric(family_table)
  }
  accession_row <- data.frame(
    sample = accession,
    total_count = total_eve_count,
    t(family_counts)  # Transpose to ensure correct row format
  )
  final_table <<- rbind(final_table, accession_row)
}

# Loop through each accession in the CSV file (main genome and haplotype accessions)
for (i in 1:nrow(accession_pairs)) {
  main_accession <- as.character(accession_pairs[i, 1])
  haplotype_accession <- as.character(accession_pairs[i, 2])
  calculate_counts(main_accession)
  calculate_counts(haplotype_accession)
}

final_table <- as.data.frame(final_table, stringsAsFactors = FALSE)
final_table$total_count <- as.numeric(final_table$total_count)
final_table[-(1:2)] <- lapply(final_table[-(1:2)], as.numeric)

# get rid of viral families with no matches
non_zero_families <- colSums(final_table[-(1:2)]) > 0
final_table <- final_table[, c(TRUE, TRUE, non_zero_families)]

#write.table(final_table, "accession_counts.tsv", sep = "\t", row.names = FALSE)

##############
# PLOTS
##############

# need to create single plot for each pair and then combine them

plot_data <- final_table %>%
  pivot_longer(cols = -c(sample, total_count), 
               names_to = "viral_family", 
               values_to = "count")
combined_plot <- NULL

for (i in 1:nrow(accession_pairs)) {
  pair <- accession_pairs[i, ]  
  pair_data <- plot_data %>%
    filter(sample %in% pair)
  
  total_counts <- final_table %>% 
    filter(sample %in% pair) %>%
    select(sample, total_count) 
  
  pair_data <- bind_rows(
    total_counts %>%
      mutate(viral_family = "Total Count", count = total_count),
    pair_data
  )
  
  pair_data$count <- log(pair_data$count + 1)
  
  # main genome should be first
  pair_data <- pair_data %>%
    mutate(sample = factor(sample, levels = c(pair[1], pair[2]))) 
  # order viral families
  pair_data$viral_family <- factor(pair_data$viral_family, 
                                   levels = c("Total Count", 
                                              colnames(final_table)[3:ncol(final_table)]))
  
  # BarPlot
  plot <- ggplot(pair_data, aes(x = viral_family, y = count, fill = sample)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          axis.title.x = element_blank()) +  
    labs(y = "EVE count (log)") +  
    scale_fill_manual(values = c("black", "grey"), name = NULL) +  
    scale_y_continuous(labels = function(x) round(exp(x) - 1))  
  
  combined_plot <- if (is.null(combined_plot)) {
    plot
  } else {
    combined_plot + plot
  }
}

combined_plot

##########
# make it dependent on genome basepairs
##########

basepair_info <- read.csv("basepair_info.tsv", sep="\t", stringsAsFactors=FALSE)
basepair_info$basepairs <- as.numeric(gsub(",", "", basepair_info$basepairs))

counts_per_bp <- final_table %>%
  left_join(basepair_info, by = "sample") %>%
  mutate(eve_percentage_per_genome = (total_count / basepairs )* 100 ) 

#write_tsv(counts_per_bp, file ="Sup_table_9.tsv")

haplotype_df <- accession_pairs %>%
  left_join(basepair_info, by = c("V1" = "sample")) %>%
  rename(main_basepairs = basepairs) %>%
  left_join(basepair_info, by = c("V2" = "sample")) %>%
  rename(haplotype_basepairs = basepairs) %>%
  mutate(haplotype_percentage_of_main = round((haplotype_basepairs / main_basepairs * 100),2)) %>%
  select(main_genome = V1, haplotype = V2, main_basepairs, haplotype_basepairs, haplotype_percentage_of_main)

# created df that includes EVEs per haplotype, basepair information, % of EVEs per basepairs
result_df <- data.frame(
  main_genome = character(),
  haplotype = character(),
  main_eve_count = numeric(),
  haplotype_eve_count = numeric(),
  eve_percentage = numeric(),
  main_basepairs = numeric(),
  haplotype_basepairs = numeric(),
  main_eve_percentage_per_genome = numeric(),
  haplotype_eve_percentage_per_genome = numeric(),
  EVE_basepair_percentage = numeric(),  # New column
  stringsAsFactors = FALSE
)

for (i in 1:nrow(accession_pairs)) {
  main_genome <- accession_pairs$V1[i]
  haplotype <- accession_pairs$V2[i]

  main_data <- counts_per_bp[counts_per_bp$sample == main_genome, ]
  haplotype_data <- counts_per_bp[counts_per_bp$sample == haplotype, ]

  if (nrow(main_data) == 1 && nrow(haplotype_data) == 1) {
    main_count <- main_data$total_count
    haplotype_count <- haplotype_data$total_count
    
    main_basepairs <- main_data$basepairs
    haplotype_basepairs <- haplotype_data$basepairs
    eve_percentage <- (haplotype_count / main_count) * 100
    main_eve_per_bp <- main_data$eve_percentage_per_genome
    haplotype_eve_per_bp <- haplotype_data$eve_percentage_per_genome
    
    EVE_basepair_percentage <- (haplotype_eve_per_bp / main_eve_per_bp) * 100
    
    result_df <- rbind(result_df, data.frame(
      main_genome = main_genome,
      haplotype = haplotype,
      main_eve_count = main_count,
      haplotype_eve_count = haplotype_count,
      eve_percentage = eve_percentage,
      main_basepairs = main_basepairs,
      haplotype_basepairs = haplotype_basepairs,
      main_eve_percentage_per_genome = main_eve_per_bp,
      haplotype_eve_percentage_per_genome = haplotype_eve_per_bp,
      EVE_basepair_percentage = EVE_basepair_percentage  # New column
    ))
  }
}


####
# STATISTICS
####

library(ggplot2)

data <- read.table("EVE_counts_assembly.tsv", header = TRUE, sep = "\t")
data$basepairs <- as.numeric(gsub(",", "", data$basepairs))
data$total_count <- as.numeric(gsub(",", "", data$total_count))
data$eve_percentage_per_genome <- as.numeric(data$eve_percentage_per_genome)
data$assembly <- factor(data$assembly)

ggplot(data, aes(x = basepairs, y = total_count, color = assembly)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Total EVE Count vs Genome Base Pairs",
       x = "Genome Base Pairs",
       y = "Total EVE Count")

shapiro.test(data$total_count)
shapiro.test(data$basepairs)

# Spearman correlation 
spearman_result <- cor.test(data$basepairs, data$total_count, method = "spearman")
print(spearman_result)

kendall_result <- cor.test(data$basepairs, data$total_count, method = "kendall")
print(kendall_result)

# Linear model with log-transformed total_count and assembly type as a factor
lm(log(total_count) ~ basepairs + assembly, data = data)


# # linearity and homoscedasticity
# plot(lm_model, which = 1) 
# 
# # Normality
# plot(lm_model, which = 2) 
