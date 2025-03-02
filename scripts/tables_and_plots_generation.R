### EVE ANALYSIS IN ANPHELES MOSQUITOES ###

# Analysis of detectEVE output and subsequent plot generation
# created by Nadja brait
# for input a combined tsv file (detectEVE output) is needed. It was pre-filtered for only high confidence hits.

# Path
setwd("C:/Users/nadja/Documents/LaptopAsus/PhD/Chapter_4")

# Libraries
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(tidyverse)

# colorbrewer 
palBrewerPlus <- c(
  "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#1B9E77", "#FB9A99", "#E31A1C",
  "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#B15928", "#1ff8ff","#FFFF99", "#D95F02", "#7570B3", "#E7298A","#66A61E", "#E6AB02", "#A6761D", "#666666",
  "#4b6a53","#b249d5","#7edc45","#5c47b8","#cfd251","#ff69b4","#69c86c","#cd3e50","#83d5af",
  "#da6130","#5e79b2","#c29545","#532a5a","#5f7b35","#c497cf","#773a27","#7cb9cb","#594e50", "#d3c4a8","#c17e7f"
)
barplot_colors <- c(
  "#A6CEE3","#5e79b2", "#1F78B4","#6A3D9A","#5c47b8","#7570B3","#c497cf","#CAB2D6","#b249d5","#E7298A","#ff69b4","#c17e7f","#cd3e50", "#E31A1C",
  "#33A02C","#66A61E","#7edc45", "#69c86c","#B2DF8A","#5f7b35","#4b6a53","#006d2c","#1B9E77","#83d5af",
  "#FFFF99","#da6130","#E6AB02","#FF7F00","#FDBF6F",
  "#d3c4a8","#594e50","#cfd251")

piechart_colors <- c(
  "#A6CEE3", "#1F78B4","#5c47b8","#CAB2D6","#b249d5",
  "#69c86c","#33A02C","#7edc45","#66A61E","#B2DF8A","#5f7b35","#4b6a53",
  "#FFFF99","#E6AB02","#FF7F00","#FDBF6F",
  "#d3c4a8","#cfd251","black")


########## analysis ###########

# reading df and updating it 
#df <- read.table('high_conf_hits_new.txt', sep = '\t', header = TRUE, stringsAsFactors = FALSE) # old
df <- read.table('high_conf_hits_final.txt', sep = '\t', header = TRUE, stringsAsFactors = FALSE)

# found out afterwards that certain accessions are either contaminated or not full assemblies and need to get rid of those rows

unwanted_accessions <- c("STHF01", "CALSDT01", "CALSDV02", "CALSDZ01", "CALSEJ01", 
                       "CAXITP01", "CALSDY02", "CALSDU01", "CALSDX01", "CALSEA01", 
                       "CALSGA01", "CALSEB02", "CALSDP01", "CALSEL02", "CALSDQ02", "CAWUOV02")

# Use filter_all to filter rows that do not contain any of the unwanted patterns
df <- df %>% filter_all(all_vars(!grepl(paste(unwanted_accessions, collapse="|"), .)))

#write.table(df, file = "Sup_Table_1.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# some rows have insufficient descriptions and taxid, which has to be fixed:
df <- df %>%
  mutate(top_viral_family = case_when(
    grepl("Baculoviral IAP repeat-containing protein 5|UBC core domain-containing protein", top_desc) ~ "Baculoviridae",
    TRUE ~ top_viral_family # Keep the original value if no conditions are met
  ))
df <- df %>%
  mutate(top_viral_desc = case_when(
    grepl("Baculoviral IAP repeat-containing protein 5", top_desc) ~ "Baculoviral IAP repeat-containing protein 5 n=33 Tax=Culicidae  RepID=A0A182RG31_ANOFN",
    grepl("UBC core domain-containing protein", top_desc) ~ "UBC core domain-containing protein n=4 Tax=Culicidae  RepID=A0A182SYF5_9DIPT",
    TRUE ~ top_viral_desc # Keep the original value if no conditions are met
  ))
df <- df %>%
  mutate(top_viral_family = case_when(
    grepl("Kaiowa virus|Croada virus|Cumbaru virus|Guato virus", top_viral_desc) ~ "Chuviridae",
    grepl("Byreska virus|orthomyxo-like virus", top_viral_desc) ~ "Orthomyxoviridae",
    grepl("reo-like virus|Hubei tetragnatha maxillosa virus 9", top_viral_desc) ~ "Sedoreoviridae",
    grepl("Bolahun virus", top_viral_desc) ~ "Xinmoviridae",
    grepl("virga-like virus", top_viral_desc) ~ "Virgaviridae",
    grepl("partiti-like virus", top_viral_desc) ~ "Partitiviridae",
    TRUE ~ top_viral_family # Keep the original value if no conditions are met
  ))
df <- df %>%
  mutate(top_viral_family = case_when(
    grepl("unclassified Reovirales family", top_viral_family) ~ "Sedoreoviridae",
    grepl("unclassified Viruses family", top_viral_family) ~ "unclassified",
    grepl("unclassified Bunyavirales family", top_viral_family) ~ "unclassified Bunyavirales",
    grepl("unclassified Mononegavirales family", top_viral_family) ~ "unclassified Mononegavirales",
    TRUE ~ top_viral_family # Keep the original value if no conditions are met
  ))

# make column for protein class
df <- df %>%
  mutate(Protein = sub("\\[.*$", "", top_viral_desc))

df <- df %>%
  mutate(Protein = case_when(
    grepl("glycoprotein|Glycoprotein|cell attachment protein|Collar head-to-tail connector protein", Protein) ~ "GP",
    grepl("PB1|PB2|PA|polymerase|polymerase PA|polymerase PB2|RNA-directed RNA polymerase|RNA-dependent RNA polymerase|replication-associated protein|RdRp|replicase|putative RNA dependent RNA polymerase|putative RdRP|VP1|RNA helicase|Protease|p125 protein", Protein) ~ "RdRp",
    grepl("nucleocapsid protein|nucleocapsid|nucleoprotein|Nucleoprotein|putative N protein", Protein) ~ "NP",
    grepl("IAP|Iap|iap|Iap-3|iap3|inhibitor of apoptosis 3|Inhibitor of apoptosis 3|inhibitor of apoptosis protein 3|inhibitor of apoptosis protein", Protein) ~ "iap3",
    grepl("Polyprotein|polyprotein", Protein) ~ "polyprotein",
    grepl("Putative NSs protein|putative non-structural protein NS1|putative NSs protein|ORF1|Nsp2|NSP2|E2|ORF3 protein n=1 Tax=Anopheles triannulatus orthophasmavirus", Protein) ~ "NS",
    grepl("WYL|Movement protein|D protein|E5|movement protein|overlapping protein/movement protein|functional/movement proteins", Protein) ~ "functional/movement proteins",
    grepl("Cap|capsid|Capsid protein|putative capsid|VP2|VP3|VP4|structural protein|ORF2|Coat protein|major core protein", Protein) ~ "Capsid",
    grepl("Uncharacterized protein|uncharacterized protein|hypothetical protein", Protein) ~ "uncharacterized protein",
    TRUE ~ Protein # Keep the original value if no conditions are met
  ))

# replication mechanism per viral family
replication <- read.table('replication_types.tsv', sep = '\t', header = TRUE, stringsAsFactors = FALSE)

df <- merge(df,replication, by="top_viral_family")
df <- df %>%
  mutate(replication = case_when(
    grepl("plant", replication) ~ "+veRNA",
    TRUE ~ replication))
    
# exclude retroviral hits, bacteriophages and dsDNA hits
#df <- df %>% filter(!str_detect(replication, "retrotransposon|bacteriophage|circular_dsDNA|dsDNA"))
df <- df[!grepl("retrotransposon|bacteriophage|circular_dsDNA|dsDNA", df$replication, ignore.case = TRUE), ]

# add Anopheles species per accession
Anopheles <- read.table('Anopheles_species.txt', sep = '\t', header = TRUE, stringsAsFactors = FALSE)
df <- merge(df,Anopheles, by="sample")

# final output file
#write.table(df, file = "updated_hits_final.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# calculate counts per protein class per viral family
segments <- df %>%
  group_by(top_viral_family, Protein) %>%
  summarise(count = n(), .groups = 'drop') %>%
  pivot_wider(names_from = Protein, values_from = count, values_fill = 0)


##########################################################
#CALCULATIONS
##########################################################

# # viral family counts for pie chart
# family_counts <- df %>%
#   mutate(top_viral_family = factor(top_viral_family, levels = replication_order)) %>%
#   group_by(top_viral_family, replication) %>%
#   summarize(count = n()) %>%
#   ungroup() %>%
#   arrange(match(top_viral_family, replication_order)) 
# 
# # Merge viral families with count less than 5 into "others"
# family_counts <- family_counts %>%
#   mutate(top_viral_family = ifelse(count < 5, "others", as.character(top_viral_family))) %>%
#   group_by(top_viral_family, replication) %>%
#   #summarize(count = n()) %>%
#   summarize(count = sum(count)) %>%
#   ungroup() #%>%
#   #arrange(desc(count))


# viral family counts for pie chart

# Order according to replication mechanism
replication_order <- df %>%
  distinct(replication, top_viral_family) %>%
  arrange(replication) %>%
  pull(top_viral_family) %>%
  unique()  

total_counts <- df %>%
  group_by(top_viral_family) %>%
  summarize(total_count = n(), .groups = "drop")

small_families <- total_counts %>%
  filter(total_count < 5) %>%
  pull(top_viral_family)

modified_df <- df %>%
  mutate(
    modified_family = ifelse(top_viral_family %in% small_families, "others", as.character(top_viral_family)),
    modified_replication = ifelse(modified_family == "others", NA, replication)
  )

# Calculate viral family counts for the pie chart
family_counts <- modified_df %>%
  mutate(modified_family = factor(modified_family, levels = c(replication_order, "others"))) %>%
  group_by(modified_replication, modified_family) %>%
  summarize(count = n(), .groups = "drop") %>%
  arrange(modified_replication, match(modified_family, c(replication_order, "others")))

# Recalculate fractions and positions for the pie chart
family_counts <- family_counts %>%
  mutate(
    fraction = count / sum(count),
    ymax = cumsum(fraction),
    ymin = lag(ymax, default = 0),
    labelPosition = (ymax + ymin) / 2,
    label = paste0("n: ", count)
  )


# viral protein counts for pie chart
protein_counts <- df %>%
  group_by(Protein) %>%
  summarize(count = n()) %>%
  arrange(desc(count))


protein_counts <- protein_counts %>%
mutate(
  fraction = count / sum(protein_counts$count),
  ymax = cumsum(fraction),
  ymin = lag(ymax, default = 0),
  labelPosition = (ymax + ymin) / 2,
  label = paste0("n: ", count))


# Viral families per sample for stacked barplot

organism_order <- readLines("species_order.txt")

family_counts_per_sample <- df %>%
  group_by(sample, top_viral_family, replication, organism_an) %>%
  summarize(count = n()) %>%
  arrange(sample, desc(count))

EVE_counts_per_sample <- df %>%
  group_by(sample) %>%
  summarize(count = n()) %>%
  arrange(sample, desc(count))

# Adjust the organism_an order to match the order from the text file
family_counts_per_sample <- family_counts_per_sample %>%
  mutate(organism_an = factor(organism_an, levels = rev(organism_order)))

# Order samples by organism
sample_order <- family_counts_per_sample %>%
  arrange(organism_an, sample) %>%
  pull(sample) %>%
  unique()

family_counts_per_sample <- family_counts_per_sample %>%
  mutate(sample = factor(sample, levels = sample_order))

family_counts_per_sample <- family_counts_per_sample %>%
  mutate(top_viral_family = factor(top_viral_family, levels = replication_order))


# # Replace families with count less than 5 with "Others"
# family_counts_per_sample <- family_counts_per_sample %>%
#   group_by(sample) %>%
#   mutate(top_viral_family = ifelse(count < 5, "Others", top_viral_family)) %>%
#   group_by(sample, top_viral_family) %>%
#   summarize(count = sum(count)) %>%
#   ungroup()


# Viral families percentages for stacked barplot

family_percentage_per_sample <- family_counts_per_sample %>%
  group_by(sample) %>%
  mutate(total_count = sum(count),
         percentage = (count / total_count) * 100) %>%
  ungroup()

#write.table(family_percentage_per_sample, file = "EVE_replication_counts.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

#########################################################
#PLOTS
#########################################################

# Pie chart for viral family counts
ggplot(family_counts, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=1, fill=modified_family)) +
  geom_rect() +
  geom_text( x=2, aes(y=labelPosition, label=label), color="black", size=6) + # x here controls label position (inner / outer)
  scale_fill_manual(values = piechart_colors) +  
  scale_color_manual(values = piechart_colors) +
  coord_polar(theta="y") +
  xlim(c(-1, 4)) +
  theme_void() +
  theme(legend.position = "right",
        legend.text = element_text(size = 14))

# Pie chart for viral proteins

ggplot(protein_counts, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=1, fill=Protein)) +
  geom_rect() +
  geom_text( x=2, aes(y=labelPosition, label=label), color="black", size=6) + # x here controls label position (inner / outer)
  #scale_fill_manual(values = palBrewerPlus) +  
  #scale_color_manual(values = palBrewerPlus) +
  scale_fill_viridis(discrete = TRUE, option = "C",na.value = "grey50")+
  scale_color_viridis(discrete = TRUE, option = "C",na.value = "grey50")+
  coord_polar(theta="y") +
  xlim(c(-1, 4)) +
  theme_void() +
  theme(legend.position = "right",
        legend.text = element_text(size = 14),
        )

# # Stacked barplot for viral families per sample
# ggplot(family_counts_per_sample, aes(x = sample, y = count, fill = top_viral_family)) +
#   geom_bar(stat = "identity", position = "stack") +
#   scale_fill_manual(values = barplot_colors) + 
#   labs(title = "Viral Family Counts per Sample",
#        x = "Sample",
#        y = "EVE count",
#        fill = "Viral Family") +
#   theme_minimal() +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1),  
#     legend.position = "bottom"  
#   ) +
#   coord_flip()  

# Calculate the total counts per sample
family_counts_per_sample$total_count <- ave(family_counts_per_sample$count, family_counts_per_sample$sample, FUN = sum)

# Stacked barplot with total count labels at the top of the bar
ggplot(family_counts_per_sample, aes(x = sample, y = count, fill = top_viral_family)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = barplot_colors) + 
  labs(title = "Viral Family Counts per Sample",
       x = "Sample",
       y = "EVE count",
       fill = "Viral Family") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  
    legend.position = "bottom"  
  ) +
  coord_flip() +
  geom_text(aes(x = sample, y = total_count, label = total_count), 
            size = 4, color = "darkgrey", vjust = 0.5, hjust = -1)




# Stacked barplot of viral family percentages per sample
ggplot(family_percentage_per_sample, aes(x = sample, y = percentage, fill = top_viral_family)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = barplot_colors) +  
  labs(title = "Percentage of Viral Families per Sample",
       x = "Genome",
       y = "Percentage",
       fill = "Viral Family") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"  
  ) +
  coord_flip()  


############

# # count check for paper
# 
# input_counts <- df %>%
#   group_by(top_viral_family,sample) %>%
#   summarize(total_count = n(), .groups = "drop")
# 
# input_counts <- df %>%
#   group_by(top_viral_family, sample) %>%
#   summarize(
#     total_count = n(),
#     organism_an = first(organism_an),
#     eve_is = first(eve_id),
#     .groups = "drop"
#   )

#########################

# calculate total EVE counts and family counts per sample for Supplementary table

accession_pairs <- unique(df$sample)
df$sample <- as.character(df$sample)
viral_families <- unique(df$top_viral_family)

final_table <- data.frame(
  sample = character(),
  total_count = integer(),
  matrix(0, ncol = length(viral_families), nrow = 0)
)
colnames(final_table) <- c("sample", "total_count", viral_families)

# Table with calculated total counts and counts per viral family for each accession
calculate_counts <- function(accession) {
  accession_data <- df[df$sample == accession, ]
  total_eve_count <- nrow(accession_data)
  family_counts <- setNames(rep(0, length(viral_families)), viral_families)
  if (total_eve_count > 0) {
    family_table <- table(accession_data$top_viral_family)
    family_counts[names(family_table)] <- as.numeric(family_table)
  }
  accession_row <- data.frame(
    sample = accession,
    total_count = total_eve_count,
    as.list(family_counts)  # Convert named vector to list
  )
  final_table <<- rbind(final_table, accession_row)
}

for (accession in accession_pairs) {
  calculate_counts(accession)
}

#final_table <- as.data.frame(final_table, stringsAsFactors = FALSE)

################

# Partitivirus specific calculations

# Count the number of rows that match the specific string in the 'top_viral_desc' column
count <- sum(df$top_viral_desc == "putative RNA dependent RNA polymerase [Hattula partiti-like virus] acc=UUV42351.1", na.rm = TRUE)
count

count <- sum(df$top_viral_desc == "putative capsid protein [Orestiada partiti-like virus] acc=QRD99858.1", na.rm = TRUE)
count

##########

# contig analysis for EVE haplotypes

##########
metadata <- read.table("input_haplotypes.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
metadata$sample <- as.character(metadata$sample)

# create contig column
metadata <- metadata %>%
  mutate(contig = sub("_[^_]*$", "", locus))

# Summarize contigs per value
contig_summary <- metadata %>%
  group_by(contig) %>%
  summarise(count = n(), .groups = 'drop')

#write_tsv(contig_summary, file = "contig_summary.tsv")

# Summarize contigs per sample
sample_contig_summary <- metadata %>%
  group_by(sample) %>%
  summarise(num_contigs = n_distinct(contig), .groups = 'drop')

# Count occurrences of contigs with only 1 EVE and more than 1 EVE
eves_count <- contig_summary %>%
  summarise(
    one_eve_count = sum(count == 1),  
    more_than_one_eve_count = sum(count > 1)  
  )
print(eves_count)


EVE_counts_per_sample <- df %>%
  group_by(organism_an,sample) %>%
  summarize(count = n()) %>%
  arrange(organism_an,sample, desc(count))

#write_tsv(EVE_counts_per_sample, file = "SuP-table5.tsv")


######
# statistical test identity versus replication type

kruskal.test(top_pident ~ replication, data = df)

# Perform pairwise - Wilcoxon tests
pairwise.wilcox.test(df$top_pident, df$replication, p.adjust.method = "BH")

mean_sd_identity <- df %>%
  group_by(replication) %>%
  summarise(
    mean_identity = mean(top_pident, na.rm = TRUE),
    sd_identity = sd(top_pident, na.rm = TRUE)
  )
print(mean_sd_identity)
#print(mean_identity)

# Perform one-way ANOVA
anova_result <- aov(top_pident ~ replication, data = df)
summary(anova_result)

# Perform Tukey's HSD test for pairwise comparisons
tukey_result <- TukeyHSD(anova_result)
print(tukey_result)

# Perform Kruskal-Wallis test for differences in identity between viral families
kruskal.test(top_pident ~ top_viral_family, data = df)

# Perform pairwise comparisons between viral families with p-value adjustment
pairwise.wilcox.test(df$top_pident, df$top_viral_family, p.adjust.method = "BH")

#######
# additional paper calculations

# average percentage identity per viral family
average_identity_per_family <- df %>%
  group_by(top_viral_family) %>%
  summarise(avg_identity = mean(top_pident, na.rm = TRUE))

# statistic percentage identity per viral family
kruskal.test(top_pident ~ top_viral_family, data = df)

# plot percentage identity per viral family
ggplot(df, aes(x = top_viral_family, y = top_pident)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Percentage Identity Across Viral Families",
       x = "Viral Family",
       y = "Percentage Identity")
