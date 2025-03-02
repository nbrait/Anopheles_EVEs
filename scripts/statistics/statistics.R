#####################
# SCRIPT FOR STATISTICS ON EVE COUNTS, SEQUENCING TECHNIQUES, BASEBAIR INFORMATION AND ANOPHELES SPECIES
#####################
# by Nadja Brait 

setwd("C:/Users/nadja/Documents/LaptopAsus/PhD/Chapter_4/PCA")

library(dplyr)
library(lme4)
#install.packages("performance")
library(performance)
#install.packages("DHARMa")
library(DHARMa)
library(ggplot2)
#install.packages("mgcv")
library(mgcv)

metadata <- read.delim("updated_hits_final.tsv", sep = "\t") 
sequencing <- read.delim("sequencing_technique.tsv", sep = "\t") 

merged_data <- merge(metadata, sequencing, by = "sample")
merged_data <- merged_data %>%
  filter(sequencing != "first") %>%
  mutate(sequencing = dplyr::recode(sequencing,
                                    "third" = "third_or_hybrid",
                                    "hybrid" = "third_or_hybrid"))

# EVE counts per sample
eve_counts <- merged_data %>%
  group_by(sample, sequencing, organism_an.x) %>%  
  summarise(EVE_counts = n()) 

if (!is.factor(eve_counts$sequencing)) {
  eve_counts$sequencing <- as.factor(eve_counts$sequencing)
}

# # Generalized linear mixed model with Anopheles species as random variable 
# glmm_model <- lmer(EVE_counts ~ sequencing + (1 | organism_an.x), data = eve_counts)
# summary(glmm_model)
# 
# # Diagnostics 
# fitted_values <- fitted(glmm_model)
# residuals <- residuals(glmm_model)
# plot(fitted_values, residuals, xlab = "Fitted Values", ylab = "Residuals")
# abline(h = 0, col = "red")
# check_model(glmm_model)
# qqnorm(residuals)
# qqline(residuals, col = "red")
# 
# # it needs to be logged

# Generalized linear mixed model with Anopheles species as random variable 
eve_counts$log_EVE_counts <- log(eve_counts$EVE_counts + 1) 
glmm_model_log <- lmer(log_EVE_counts ~ sequencing + (1 | organism_an.x), data = eve_counts)
summary(glmm_model_log)

# Diagnostics
fitted_values <- fitted(glmm_model_log)
residuals <- residuals(glmm_model_log)
plot(fitted_values, residuals, xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, col = "red")
check_model(glmm_model_log)
qqnorm(residuals)
qqline(residuals, col = "red")

#######################

# include basepair information to model
basepair_info <- read.csv("basepair_info.tsv", sep="\t", stringsAsFactors=FALSE)
basepair_info$basepairs <- as.numeric(gsub(",", "", basepair_info$basepairs))
EVE_combined <- merge(eve_counts, basepair_info, by = "sample", all.x = TRUE)

# Standardize basepairs and make GLMM
EVE_combined$basepairs_scaled <- scale(EVE_combined$basepairs)
glmm_model_log <- lmer(log_EVE_counts ~ sequencing + basepairs_scaled + (1 | organism_an.x), data = EVE_combined)
summary(glmm_model_log)

# vif_model <- lm(log_EVE_counts ~ sequencing + basepairs_scaled, data = EVE_combined)
# vif_results <- vif(vif_model)
# print(vif_results)

###########
# # ANOVA to compare the models and see if sequencing and basepairs don't interfere with each other

full_model <- lmer(log_EVE_counts ~ sequencing + basepairs_scaled + (1 | organism_an.x), data = EVE_combined)
reduced_model_seq <- lmer(log_EVE_counts ~ basepairs_scaled + (1 | organism_an.x), data = EVE_combined)
anova(reduced_model_seq, full_model)

reduced_model_bp <- lmer(log_EVE_counts ~ sequencing + (1 | organism_an.x), data = EVE_combined)
anova(reduced_model_bp, full_model)

###
# Partitiviruses might bias statistics - due to their abundance
### 

# same as aboth just without Partiti matches

non_partiti_metadata <- metadata %>%
  filter(top_viral_family != "Partitiviridae")

non_partiti <- merge(non_partiti_metadata, sequencing, by = "sample")

non_partiti <- non_partiti %>%
  filter(sequencing != "first") %>%
  mutate(sequencing = dplyr::recode(sequencing,
                                    "third" = "third_or_hybrid",
                                    "hybrid" = "third_or_hybrid"))
eve_counts <- non_partiti %>%
  group_by(sample, sequencing, organism_an.x) %>% 
  summarise(EVE_counts = n())  
if (!is.factor(eve_counts$sequencing)) {
  eve_counts$sequencing <- as.factor(eve_counts$sequencing)
}

eve_counts$log_EVE_counts <- log(eve_counts$EVE_counts + 1) 
EVE_combined <- merge(eve_counts, basepair_info, by = "sample", all.x = TRUE)
EVE_combined$basepairs_scaled <- scale(EVE_combined$basepairs)

# GLMM
glmm_model_log <- lmer(log_EVE_counts ~ sequencing + basepairs_scaled + (1 | organism_an.x), data = EVE_combined)
summary(glmm_model_log)

# ANOVA
full_model <- lmer(log_EVE_counts ~ sequencing + basepairs_scaled + (1 | organism_an.x), data = EVE_combined)
reduced_model_seq <- lmer(log_EVE_counts ~ basepairs_scaled + (1 | organism_an.x), data = EVE_combined)
anova(reduced_model_seq, full_model)

# Diagnostics
fitted_values <- fitted(glmm_model_log)
residuals <- residuals(glmm_model_log)
plot(fitted_values, residuals, xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, col = "red")
check_model(glmm_model_log)
qqnorm(residuals)
qqline(residuals, col = "red")


########################
## Poisson model
# poisson_model <- glmer(EVE_counts ~ sequencing + (1 | organism_an.x), family = poisson, data = eve_counts)
# 
# # Calculate Residual Deviance and Degrees of Freedom
# residual_deviance <- deviance(poisson_model)
# df_residual <- df.residual(poisson_model)
# 
# # Calculate Dispersion Statistic
# dispersion_stat <- residual_deviance / df_residual
# dispersion_stat

# ## negative binomial model
# nb_model <- glmer.nb(EVE_counts ~ sequencing + (1 | organism_an.x), data = eve_counts)
# # Summary of the model
# summary(nb_model)
# 
# check_model(nb_model)
# 
# # Residual diagnostics using DHARMa package
# # Simulate residuals for the model
# sim_res <- simulateResiduals(fittedModel = nb_model, n = 1000)
# 
# # Plotting the residuals
# plot(sim_res)
# 
# # Q-Q plot of residuals
# qqnorm(sim_res$scaledResiduals)
# qqline(sim_res$scaledResiduals, col = "red")
# 
# # Check for overdispersion
# testDispersion(sim_res)
# 
# # Posterior Predictive Check
# plot(sim_res, rank = TRUE)
# 
# # Check for zero inflation
# testZeroInflation(sim_res)
# 
# # Residuals vs Fitted Values
# plotResiduals(sim_res, form = fitted(nb_model))
# 
# # Homogeneity of variance check
# plot(sim_res, quantreg = TRUE)
# 
# # Linearity check (Optional)
# # Here, the assumption is that there's no clear pattern in residuals vs fitted values
# plotResiduals(sim_res, form = fitted(nb_model))
# 
# fitted_values <- predict(glmm_model_log, type = "response")
# residuals <- residuals(glmm_model_log, type = "pearson")
# 
# ggplot(data = eve_counts, aes(x = fitted_values, y = residuals)) +
#   geom_point() +
#   geom_smooth(method = "loess") +
#   labs(title = "Residual vs. Fitted Plot")
# 
# # Diagnostics
# fitted_values <- fitted(nb_model)
# residuals <- residuals(nb_model)
# # Plot
# plot(fitted_values, residuals, xlab = "Fitted Values", ylab = "Residuals")
# abline(h = 0, col = "red")

################

## checking only for Partitiviridae 

################

partitiviridae_metadata <- metadata %>%
  filter(top_viral_family == "Partitiviridae")

Partiti <- merge(partitiviridae_metadata, sequencing, by = "sample")
Partiti <- Partiti %>%
  filter(sequencing != "first") %>%
  mutate(sequencing = recode(sequencing,
                             "third" = "third_or_hybrid",
                             "hybrid" = "third_or_hybrid"))

eve_counts <- Partiti %>%
  group_by(sample, sequencing, organism_an.x) %>%  
  summarise(EVE_counts = n())  

if (!is.factor(eve_counts$sequencing)) {
  eve_counts$sequencing <- as.factor(eve_counts$sequencing)
}
eve_counts$log_EVE_counts <- log(eve_counts$EVE_counts + 1) 
glmm_model_log <- lmer(log_EVE_counts ~ sequencing + (1 | organism_an.x), data = eve_counts)
summary(glmm_model_log)

# Diagnostics
# fitted_values <- fitted(glmm_model_log)
# residuals <- residuals(glmm_model_log)
# plot(fitted_values, residuals, xlab = "Fitted Values", ylab = "Residuals")
# abline(h = 0, col = "red")
# check_model(glmm_model_log)
# qqnorm(residuals)
# qqline(residuals, col = "red")

########

non_partiti_metadata <- metadata %>%
  filter(top_viral_family != "Partitiviridae")

non_partiti <- merge(non_partiti_metadata, sequencing, by = "sample")

non_partiti <- non_partiti %>%
  filter(sequencing != "first") %>%
  mutate(sequencing = recode(sequencing,
                             "third" = "third_or_hybrid",
                             "hybrid" = "third_or_hybrid"))

# Calculate EVE counts per sample
eve_counts <- non_partiti %>%
  group_by(sample, sequencing, organism_an.x) %>%  
  summarise(EVE_counts = n()) 

# Convert 'sequencing' to a factor if not already
if (!is.factor(eve_counts$sequencing)) {
  eve_counts$sequencing <- as.factor(eve_counts$sequencing)
}

eve_counts$log_EVE_counts <- log(eve_counts$EVE_counts + 1) 
glmm_model_log <- lmer(log_EVE_counts ~ sequencing + (1 | organism_an.x), data = eve_counts)
summary(glmm_model_log)

# Diagnostics
fitted_values <- fitted(glmm_model_log)
residuals <- residuals(glmm_model_log)
# Plothttp://127.0.0.1:32039/graphics/cba0433f-bd60-452f-a764-d71c329b8193.png
plot(fitted_values, residuals, xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, col = "red")

check_model(glmm_model_log)

qqnorm(residuals)
qqline(residuals, col = "red")

####################
# MANN-WHITNEY-U TEST
####################

# Read in your data tables
metadata <- read.delim("updated_hits.tsv", sep = "\t")  
sequencing <- read.delim("sequencing_technique.tsv", sep = "\t") 

# Merge the tables on the 'sample' column
merged_data <- merge(metadata, sequencing, by = "sample")
merged_data <- merged_data %>%
  filter(sequencing != "first") %>%
  mutate(sequencing = recode(sequencing,
                             "third" = "third_or_hybrid",
                             "hybrid" = "third_or_hybrid"))

# Calculate EVE counts per sample
eve_counts <- merged_data %>%
  group_by(sample, sequencing, organism_an.x) %>%  
  summarise(EVE_counts = n())  

eve_counts <- eve_counts %>%
  mutate(sequencing = as.factor(sequencing),
         sample = as.factor(sample))

# Split by sequencing technique
group1 <- eve_counts %>% filter(sequencing == "second")
group2 <- eve_counts %>% filter(sequencing == "third_or_hybrid")

# Mann-Whitney U Test
test_result <- wilcox.test(group1$EVE_counts, group2$EVE_counts)
print(test_result)

#######
#### Using anova() for Likelihood Ratio Tests
#######
# will give a p-value based on a chi-square test comparing the two models

# Full model with sequencing
full_model <- lmer(log_EVE_counts ~ sequencing + (1 | organism_an.x), data = eve_counts)

# Reduced model without sequencing
reduced_model <- lmer(log_EVE_counts ~ (1 | organism_an.x), data = eve_counts)

# Perform likelihood ratio test
anova(reduced_model, full_model)

# try the same but whether organism_an adds to the model or not

reduced_model <- lmer(log_EVE_counts ~ sequencing + (1 | 1), data = eve_counts)

full_model <- lmer(log_EVE_counts ~ sequencing + (1 | organism_an.x), data = eve_counts)
anova(reduced_model, full_model)

# ##########################
# 
# # Load necessary packages
# library(dplyr)
# library(reshape2)
# library(vegan) # for Jaccard index
# library(ggplot2)
# 
# # Create a function to calculate pairwise Euclidean distance for EVE counts
# calculate_distance_matrix <- function(data) {
#   dist(as.matrix(data))
# }
# 
# # Aggregate total counts per sample
# input_counts_summarized <- input_counts %>%
#   group_by(sample) %>%
#   summarise(total_count = sum(total_count), .groups = 'drop')
# 
# # Extract sample names before spreading
# sample_names <- input_counts_summarized$sample
# 
# # Spread the data
# evecounts <- input_counts_summarized %>%
#   spread(sample, total_count)
# 
# # Convert to matrix
# evecounts_matrix <- as.matrix(evecounts[,-1])
# rownames(evecounts_matrix) <- evecounts$sample
# 
# # Example data preparation
# # Assume `df` is your dataframe with EVE counts and organism info
# evecounts <- input_counts %>% select(sample, total_count, organism_an) %>% spread(sample, total_count)
# evecounts_matrix <- as.matrix(evecounts[,-1])
# rownames(evecounts_matrix) <- evecounts$sample
# 
# # Calculate Euclidean distance for EVE counts
# distance_matrix <- calculate_distance_matrix(evecounts_matrix)
# 
# # Convert to data frame for easier manipulation
# distance_df <- as.data.frame(as.table(distance_matrix))
# distance_df <- distance_df %>%
#   rename(Sample1 = Var1, Sample2 = Var2, Distance = Freq) %>%
#   mutate(SameOrganism = ifelse(df$organism_an[match(Sample1, df$sample)] == df$organism_an[match(Sample2, df$sample)], "Same", "Different"))
# 
# # Perform Mann-Whitney U test
# wilcox.test(Distance ~ SameOrganism, data = distance_df)
# 
# # Optional: Jaccard Similarity for viral families
# viral_families <- df %>% select(sample, top_viral_family) %>% spread(sample, top_viral_family)
# jaccard_matrix <- vegdist(as.matrix(viral_families[,-1]), method = "jaccard")
# # Convert and analyze Jaccard matrix similarly
