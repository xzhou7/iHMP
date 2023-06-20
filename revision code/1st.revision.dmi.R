
#revision added analysis
#library(tidyverse)
library(dplyr)
library(phyloseq)
library(cowplot)
library(ggpubr)
library(patchwork)

setwd("~/Library/CloudStorage/Box-Box/")
getwd()
source("./human_microbiome_project/code/tools.R")

# diversity_all <- read.csv("./XinZhouFiles/Projects/ZXP2_HMP2_Core/Analysis/Diversity Table/All_Diversity.csv", header = T)
# richness_by_genus <- load("./XinZhouFiles/Projects/ZXP2_HMP2_Core/Analysis/Robject/Richness_Bygenus.RData")
# diversity_by_genus <- load("./XinZhouFiles/Projects/ZXP2_HMP2_Core/Analysis/Robject/Diversity_Datatable.RData")

load("./XinZhouFiles/Projects/ZXP2_HMP2_Core/Analysis/Robject/Diversity_Bygenus.RData")
load("./XinZhouFiles/Projects/ZXP2_HMP2_Core/Analysis/Robject/Prevalance.RData")

load("./human_microbiome_project/data_analysis/combine_microbiome/distance/stool/personalized_score_permutation")
load("./human_microbiome_project/data_analysis/combine_microbiome/distance/stool/personalized_score_permutation_trim")

stool.score <- personalized_score_permutation
stool.trim <- personalized_score_permutation_trim
stool.trim1 <- stool.trim

stool.score$trimmed <- ifelse(row.names(stool.score) %in% row.names(stool.trim), "no", "yes")
identical(filter(stool.score, trimmed=="no") %>% dplyr::select(-trimmed), stool.trim)

load("./human_microbiome_project/data_analysis/combine_microbiome/distance/skin/personalized_score_permutation")
load("./human_microbiome_project/data_analysis/combine_microbiome/distance/skin/personalized_score_permutation_trim")

skin.score <- personalized_score_permutation
skin.trim <- personalized_score_permutation_trim
skin.trim1 <- skin.trim

skin.score$trimmed <- ifelse(row.names(skin.score) %in% row.names(skin.trim), "no", "yes")
identical(filter(skin.score, trimmed=="no") %>% dplyr::select(-trimmed), skin.trim)

load("./human_microbiome_project/data_analysis/combine_microbiome/distance/oral/personalized_score_permutation")
load("./human_microbiome_project/data_analysis/combine_microbiome/distance/oral/personalized_score_permutation_trim")

oral.score <- personalized_score_permutation
oral.trim <- personalized_score_permutation_trim
oral.trim1 <- oral.trim

oral.score$trimmed <- ifelse(row.names(oral.score) %in% row.names(oral.trim), "no", "yes")
identical(filter(oral.score, trimmed=="no") %>% dplyr::select(-trimmed), oral.trim)

load("./human_microbiome_project/data_analysis/combine_microbiome/distance/nasal/personalized_score_permutation")
load("./human_microbiome_project/data_analysis/combine_microbiome/distance/nasal/personalized_score_permutation_trim")

nasal.score <- personalized_score_permutation
nasal.trim <- personalized_score_permutation_trim
nasal.trim1 <- nasal.trim

nasal.score$trimmed <- ifelse(row.names(nasal.score) %in% row.names(nasal.trim), "no", "yes")
identical(filter(nasal.score, trimmed=="no") %>% dplyr::select(-trimmed), nasal.trim)


dim(nasal.score)
dim(nasal.trim)

rm(personalized_score_permutation,personalized_score_permutation_trim)

DMI_score <- read.csv("./human_microbiome_project/Supplementary_data/Single tables/DMIscore.csv", header = T)

# create separate objects for each class
DMI_score_nasal <- DMI_score[DMI_score$class == "Nasal", ]
DMI_score_oral <- DMI_score[DMI_score$class == "Oral", ]
DMI_score_skin <- DMI_score[DMI_score$class == "Skin", ]
DMI_score_stool <- DMI_score[DMI_score$class == "Stool", ]

DMI_score_nasal$genus
oral.score$genus

stool_DMI_merged <- merge(DMI_score_stool, stool.score, by = "genus")
skin_DMI_merged <- merge(DMI_score_skin, skin.score, by = "genus")
oral_DMI_merged <- merge(DMI_score_oral,oral.score, by = "genus")
nasal_DMI_merged <- merge(DMI_score_nasal, nasal.score, by = "genus")

DMI_score

# create an empty vector to store the mean values
#Stool
mean_values <- numeric(length(colnames(richness_ST_byGenus)))

for (i in 1:length(colnames(richness_ST_byGenus))) {
  mean_values[i] <- mean(richness_ST_byGenus[,i], na.rm=TRUE)
}
names(mean_values) <- colnames(richness_ST_byGenus)

stool_ACE_bygenus <- mean_values %>%
  data.frame() %>%
  setNames("ACE") %>%
  filter(!is.na(ACE))

rm(mean_values)

#Skin
mean_values <- numeric(length(colnames(richness_SK_byGenus)))

for (i in 1:length(colnames(richness_SK_byGenus))) {
  mean_values[i] <- mean(richness_SK_byGenus[,i], na.rm=TRUE)
}

names(mean_values) <- colnames(richness_SK_byGenus)

skin_ACE_bygenus <- mean_values %>%
  data.frame() %>%
  setNames("ACE") %>%
  filter(!is.na(ACE))

rm(mean_values)

#Oral
mean_values <- numeric(length(colnames(richness_OR_byGenus)))

for (i in 1:length(colnames(richness_OR_byGenus))) {
  mean_values[i] <- mean(richness_OR_byGenus[,i], na.rm=TRUE)
}

names(mean_values) <- colnames(richness_OR_byGenus)

oral_ACE_bygenus <- mean_values %>%
  data.frame() %>%
  setNames("ACE") %>%
  filter(!is.na(ACE))

#Nasal
mean_values <- numeric(length(colnames(richness_NS_byGenus)))

for (i in 1:length(colnames(richness_NS_byGenus))) {
  mean_values[i] <- mean(richness_NS_byGenus[,i], na.rm=TRUE)
}

names(mean_values) <- colnames(richness_NS_byGenus)

nasal_ACE_bygenus <- mean_values %>%
  data.frame() %>%
  setNames("ACE") %>%
  filter(!is.na(ACE))

rm(mean_values)

stool_ACE_bygenus$bodysite <- "Stool"
skin_ACE_bygenus$bodysite <- "Skin"
oral_ACE_bygenus$bodysite <- "Oral"
nasal_ACE_bygenus$bodysite <- "Nasal"

#stool
mean_values <- numeric(length(colnames(diversity_ST_byGenus)))
for (i in 1:length(colnames(diversity_ST_byGenus))) {
  mean_values[i] <- mean(diversity_ST_byGenus[,i], na.rm=TRUE)
}
names(mean_values) <- colnames(diversity_ST_byGenus)

stool_Shannon_bygenus <- mean_values %>%
  data.frame() %>%
  setNames("Shannon") %>%
  filter(!is.na(Shannon)) %>%
  mutate(bodysite = "Stool")

rm(mean_values)

# Skin
mean_values <- numeric(length(colnames(diversity_SK_byGenus)))
for (i in 1:length(colnames(diversity_SK_byGenus))) {
  mean_values[i] <- mean(diversity_SK_byGenus[,i], na.rm=TRUE)
}
names(mean_values) <- colnames(diversity_SK_byGenus)

skin_Shannon_bygenus <- mean_values %>%
  data.frame() %>%
  setNames("Shannon") %>%
  filter(!is.na(Shannon)) %>%
  mutate(bodysite = "Skin")

rm(mean_values)

# Oral
mean_values <- numeric(length(colnames(diversity_OR_byGenus)))
for (i in 1:length(colnames(diversity_OR_byGenus))) {
  mean_values[i] <- mean(diversity_OR_byGenus[,i], na.rm=TRUE)
}
names(mean_values) <- colnames(diversity_OR_byGenus)

oral_Shannon_bygenus <- mean_values %>%
  data.frame() %>%
  setNames("Shannon") %>%
  filter(!is.na(Shannon)) %>%
  mutate(bodysite = "Oral")

rm(mean_values)

# Nasal
mean_values <- numeric(length(colnames(diversity_NS_byGenus)))
for (i in 1:length(colnames(diversity_NS_byGenus))) {
  mean_values[i] <- mean(diversity_NS_byGenus[,i], na.rm=TRUE)
}
names(mean_values) <- colnames(diversity_NS_byGenus)

nasal_Shannon_bygenus <- mean_values %>%
  data.frame() %>%
  setNames("Shannon") %>%
  filter(!is.na(Shannon)) %>%
  mutate(bodysite = "Nasal")
rm(mean_values)

stool_Shannon_bygenus$bodysite <- "Stool"
skin_Shannon_bygenus$bodysite <- "Skin"
oral_Shannon_bygenus$bodysite <- "Oral"
nasal_Shannon_bygenus$bodysite <- "Nasal"


stool_Pre_bygenus <- colMeans(ST.Pr) %>% data.frame() %>% setNames("Prevalance")
skin_Pre_bygenus <- colMeans(SK.Pr) %>% data.frame() %>% setNames("Prevalance")
oral_Pre_bygenus <- colMeans(OR.Pr) %>% data.frame() %>% setNames("Prevalance")
nasal_Pre_bygenus <- colMeans(NS.Pr) %>% data.frame() %>% setNames("Prevalance")


head(stool_Pre_bygenus)
head(stool_Shannon_bygenus)
head(stool_ACE_bygenus)

#stool
stool_Pre_bygenus <- rownames_to_column(stool_Pre_bygenus, var = "rowname")
stool_Shannon_bygenus <- rownames_to_column(stool_Shannon_bygenus, var = "rowname")
stool_ACE_bygenus <- rownames_to_column(stool_ACE_bygenus, var = "rowname")

merged_data <- inner_join(stool_Pre_bygenus, stool_Shannon_bygenus, by = "rowname") %>%
  dplyr::select(-bodysite)

merged_data <- inner_join(merged_data, stool_ACE_bygenus, by = "rowname")

rownames(merged_data) <- merged_data$rowname
merged_data$rowname <- NULL

Stool.merged <- merged_data
rm(merged_data)

# skin
skin_Pre_bygenus <- rownames_to_column(skin_Pre_bygenus, var = "rowname")
skin_Shannon_bygenus <- rownames_to_column(skin_Shannon_bygenus, var = "rowname")
skin_ACE_bygenus <- rownames_to_column(skin_ACE_bygenus, var = "rowname")

Skin.merged <- inner_join(skin_Pre_bygenus, skin_Shannon_bygenus, by = "rowname") %>%
  dplyr::select(-bodysite)

Skin.merged <- inner_join(Skin.merged, skin_ACE_bygenus, by = "rowname")

rownames(Skin.merged) <- Skin.merged$rowname
Skin.merged$rowname <- NULL

#oral
oral_Pre_bygenus <- rownames_to_column(oral_Pre_bygenus, var = "rowname")
oral_Shannon_bygenus <- rownames_to_column(oral_Shannon_bygenus, var = "rowname")
oral_ACE_bygenus <- rownames_to_column(oral_ACE_bygenus, var = "rowname")

Oral.merged <- inner_join(oral_Pre_bygenus, oral_Shannon_bygenus, by = "rowname") %>%
  dplyr::select(-bodysite)

Oral.merged <- inner_join(Oral.merged, oral_ACE_bygenus, by = "rowname") 

rownames(Oral.merged) <- Oral.merged$rowname
Oral.merged$rowname <- NULL

#Nasal
nasal_Pre_bygenus <- rownames_to_column(nasal_Pre_bygenus, var = "rowname")
nasal_Shannon_bygenus <- rownames_to_column(nasal_Shannon_bygenus, var = "rowname")
nasal_ACE_bygenus <- rownames_to_column(nasal_ACE_bygenus, var = "rowname")

Nasal.merged <- inner_join(nasal_Pre_bygenus, nasal_Shannon_bygenus, by = "rowname") %>%
  dplyr::select(-bodysite)

Nasal.merged <- inner_join(Nasal.merged, nasal_ACE_bygenus, by = "rowname") 

rownames(Nasal.merged) <- Nasal.merged$rowname
Nasal.merged$rowname <- NULL

# Remove unnecessary variables
rm(stool_Pre_bygenus, stool_Shannon_bygenus, stool_ACE_bygenus,
  skin_Pre_bygenus, skin_Shannon_bygenus, skin_ACE_bygenus,
   oral_Pre_bygenus, oral_Shannon_bygenus, oral_ACE_bygenus,
   nasal_Pre_bygenus, nasal_Shannon_bygenus, nasal_ACE_bygenus)

# Show the first few rows of the merged data for each body site
head(Stool.merged)
head(Skin.merged)
head(Oral.merged)
head(Nasal.merged)


# # select relevant columns and filter for non-missing values
# df <- stool_DMI_merged %>% select(trimmed, fc1, dmi) %>% na.omit()
# 
# # calculate correlation coefficient and p value
# corr <- cor.test(df$fc1, df$dmi)
# 
# # create plot
# pstool <- ggplot(df, aes(x = fc1, y = dmi, color = trimmed)) +  geom_point() + 
#   geom_smooth(method = "lm", se = FALSE) +
#   scale_color_manual(values = c("red", "blue")) +
#   annotate("text", x = Inf, y = Inf, hjust = 1, vjust = 1,
#            label = paste0("Corr: ", round(corr$estimate, 2), 
#                           ", p = ", format(corr$p.value, scientific = TRUE)))
# 
# pstool

phylum_order <- c("Actinobacteria", "Bacteroidetes", "Firmicutes", "Proteobacteria", "Other")

stool.data.frame <- merge(stool_DMI_merged, Stool.merged, by.x = "genus", by.y = "row.names", all = T)
skin.data.frame <- merge(skin_DMI_merged, Skin.merged, by.x = "genus", by.y = "row.names", all = T)
oral.data.frame <- merge(oral_DMI_merged, Oral.merged, by.x = "genus", by.y = "row.names", all = T)
nasal.data.frame <- merge(nasal_DMI_merged, Nasal.merged, by.x = "genus", by.y = "row.names", all = T)


combined_div_DMI_data_raw <- rbind(stool.data.frame, skin.data.frame,oral.data.frame,nasal.data.frame)
combined_div_DMI_data <- filter(combined_div_DMI_data_raw,!is.na(combined_div_DMI_data_raw$dmi))

table(combined_div_DMI_data$bodysite,combined_div_DMI_data$trimmed)

combined_div_DMI_data

model <- lm(dmi ~ bodysite * Shannon + Prevalance, data = combined_div_DMI_data)
summary(model)

# Get the residuals and fitted values from the model
residuals <- resid(model)
fitted_values <- fitted(model)

# Create a data frame with the residuals and fitted values
resid_df <- data.frame(Residuals = residuals, FittedValues = fitted_values)

# Plot the residuals against the fitted values
residual_plot_before <- ggplot(resid_df, aes(x = FittedValues, y = Residuals)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(x = "Fitted Values", y = "Residuals",
       title = "Residual Plot before Bootstrapping")
residual_plot_before

stool.data.frame.trimonly <- filter(stool.data.frame, trimmed == "no")
skin.data.frame.trimonly <- filter(skin.data.frame, trimmed == "no")
oral.data.frame.trimonly <- filter(oral.data.frame, trimmed == "no")
nasal.data.frame.trimonly <- filter(nasal.data.frame, trimmed == "no")

stool.data.frame.trimonly$Phylum <- factor(stool.data.frame.trimonly$Phylum, level = phylum_order)
skin.data.frame.trimonly$Phylum <- factor(skin.data.frame.trimonly$Phylum, level = phylum_order)
oral.data.frame.trimonly$Phylum <- factor(oral.data.frame.trimonly$Phylum, level = phylum_order)
nasal.data.frame.trimonly$Phylum <- factor(nasal.data.frame.trimonly$Phylum, level = phylum_order)

stool.compare <- stool.data.frame.trimonly %>% ggplot(aes(x=Phylum, y=dmi)) + 
  geom_jitter(aes(color=Phylum)) + 
  geom_boxplot(alpha = 0.4, outlier.alpha = 0) + base_theme +
  scale_color_manual(values = phylum_color) +
  stat_compare_means(comparisons = list(c("Actinobacteria", "Bacteroidetes"),
                                        c("Actinobacteria", "Firmicutes"),
                                        c("Actinobacteria", "Proteobacteria"),
                                        c("Actinobacteria", "Other"),
                                        c("Bacteroidetes", "Firmicutes"),
                                        c("Bacteroidetes", "Proteobacteria"),
                                        c("Bacteroidetes", "Other"),
                                        c("Firmicutes", "Proteobacteria"),
                                        c("Firmicutes", "Other"),
                                        c("Proteobacteria", "Other")),
                     method = "wilcox.test",
                     label = "p.signif",
                     p.adjust.method = "BH",
                     hide.ns = T) + ggtitle("The Stool Microbiome") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

skin.compare <- skin.data.frame.trimonly %>% ggplot(aes(x=Phylum, y=dmi)) +
  geom_jitter(aes(color=Phylum)) +
  geom_boxplot(alpha = 0.4, outlier.alpha = 0) + base_theme +
  scale_color_manual(values = phylum_color) +
  stat_compare_means(comparisons = list(c("Actinobacteria", "Bacteroidetes"),
                                        c("Actinobacteria", "Firmicutes"),
                                        c("Actinobacteria", "Proteobacteria"),
                                        c("Actinobacteria", "Other"),
                                        c("Bacteroidetes", "Firmicutes"),
                                        c("Bacteroidetes", "Proteobacteria"),
                                        c("Bacteroidetes", "Other"),
                                        c("Firmicutes", "Proteobacteria"),
                                        c("Firmicutes", "Other"),
                                        c("Proteobacteria", "Other")),
                     method = "wilcox.test",
                     label = "p.signif",
                     p.adjust.method = "BH",
                     hide.ns = T) + ggtitle("The Skin Microbiome") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
skin.compare

oral.compare <- oral.data.frame.trimonly %>% ggplot(aes(x=Phylum, y=dmi)) +
  geom_jitter(aes(color=Phylum)) +
  geom_boxplot(alpha = 0.4, outlier.alpha = 0) + base_theme +
  scale_color_manual(values = phylum_color) +
  stat_compare_means(comparisons = list(c("Actinobacteria", "Bacteroidetes"),
                                        c("Actinobacteria", "Firmicutes"),
                                        c("Actinobacteria", "Proteobacteria"),
                                        c("Actinobacteria", "Other"),
                                        c("Bacteroidetes", "Firmicutes"),
                                        c("Bacteroidetes", "Proteobacteria"),
                                        c("Bacteroidetes", "Other"),
                                        c("Firmicutes", "Proteobacteria"),
                                        c("Firmicutes", "Other"),
                                        c("Proteobacteria", "Other")),
                     method = "wilcox.test",
                     label = "p.signif",
                     p.adjust.method = "BH",
                     hide.ns = T) + ggtitle("The Oral Microbiome") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
oral.compare


nasal.compare <- nasal.data.frame.trimonly %>% ggplot(aes(x=Phylum, y=dmi)) +
  geom_jitter(aes(color=Phylum)) +
  geom_boxplot(alpha = 0.4, outlier.alpha = 0) + base_theme +
  scale_color_manual(values = phylum_color) +
  stat_compare_means(comparisons = list(c("Actinobacteria", "Bacteroidetes"),
                                        c("Actinobacteria", "Firmicutes"),
                                        c("Actinobacteria", "Proteobacteria"),
                                        c("Actinobacteria", "Other"),
                                        c("Bacteroidetes", "Firmicutes"),
                                        c("Bacteroidetes", "Proteobacteria"),
                                        c("Bacteroidetes", "Other"),
                                        c("Firmicutes", "Proteobacteria"),
                                        c("Firmicutes", "Other"),
                                        c("Proteobacteria", "Other")),
                     method = "wilcox.test",
                     label = "p.signif",
                     p.adjust.method = "BH",
                     hide.ns = T) + ggtitle("The Nasal Microbiome") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
nasal.compare 

stool.compare + skin.compare + oral.compare + nasal.compare + plot_layout(guides = "collect")

# ggsave("~/Desktop/1st.revision.CELL-S-22-04878/Figures/Figure2/Statistics/stool.stat.pdf",stool.compare, width = 4,height = 4, dpi=300)
# ggsave("~/Desktop/1st.revision.CELL-S-22-04878/Figures/Figure2/Statistics/skin.stat.pdf",skin.compare, width = 4,height = 4, dpi=300)
# ggsave("~/Desktop/1st.revision.CELL-S-22-04878/Figures/Figure2/Statistics/oral.stat.pdf",oral.compare, width = 4,height = 4, dpi=300)
# ggsave("~/Desktop/1st.revision.CELL-S-22-04878/Figures/Figure2/Statistics/nasal.stat.pdf",nasal.compare, width = 4,height = 4, dpi=300)


combined_div_DMI_trimmed_raw <- rbind(stool.data.frame.trimonly, skin.data.frame.trimonly,oral.data.frame.trimonly,nasal.data.frame.trimonly)
combined_div_DMI_trimmed_raw$dmi
model_trimmed <- lm(dmi ~ bodysite * Shannon + Prevalance, data = combined_div_DMI_trimmed_raw)
summary(model_trimmed)

# Get the residuals and fitted values from the model
residuals <- resid(model_trimmed)
fitted_values <- fitted(model_trimmed)

# Create a data frame with the residuals and fitted values
resid_df <- data.frame(Residuals = residuals, FittedValues = fitted_values)

# Plot the residuals against the fitted values
residual_plot <- ggplot(resid_df, aes(x = FittedValues, y = Residuals)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(x = "Fitted Values", y = "Residuals",
       title = "Residual Plot after Bootstrapping")


library(ggradar)

df <- data.frame(matrix(runif(30), ncol = 10))
df[, 1] <- paste0("G", 1:3)
colnames(df) <- c("Group", paste("Var", 1:9))

ggradar(df)

Stool.Radar <- stool.data.frame.trimonly %>%
  group_by(Phylum) %>%
  summarize(Stool_mean_dmi = mean(dmi), Stool_sd_dmi = sd(dmi)) %>% 
  mutate(Stool_uppersd = Stool_mean_dmi + Stool_sd_dmi) %>% 
  mutate(Stool_lowersd = Stool_mean_dmi - Stool_sd_dmi)

Skin.Radar <- skin.data.frame.trimonly %>%
  group_by(Phylum) %>%
  summarize(Skin_mean_dmi = mean(dmi), Skin_sd_dmi = sd(dmi)) %>%
  mutate(Skin_uppersd = Skin_mean_dmi + Skin_sd_dmi) %>%
  mutate(Skin_lowersd = Skin_mean_dmi - Skin_sd_dmi)

Oral.Radar <- oral.data.frame.trimonly %>%
  group_by(Phylum) %>%
  summarize(Oral_mean_dmi = mean(dmi), Oral_sd_dmi = sd(dmi)) %>% 
  mutate(Oral_uppersd = Oral_mean_dmi + Oral_sd_dmi) %>% 
  mutate(Oral_lowersd = Oral_mean_dmi - Oral_sd_dmi)

Nasal.Radar <- nasal.data.frame.trimonly %>%
  group_by(Phylum) %>%
  summarize(Nasal_mean_dmi = mean(dmi), Nasal_sd_dmi = sd(dmi)) %>% 
  mutate(Nasal_uppersd = Nasal_mean_dmi + Nasal_sd_dmi) %>% 
  mutate(Nasal_lowersd = Nasal_mean_dmi - Nasal_sd_dmi)

Nasal.Radar


radar.all <- cbind(Stool.Radar,Skin.Radar,Oral.Radar,Nasal.Radar)
radar.all

radar.all <- radar.all[, !duplicated(colnames(radar.all))]
radar.all <- dplyr::select(radar.all, select = -c(ends_with("_sd_dmi")))
ggradar(radar.all)

radar.mean <- dplyr::select(radar.all, Phylum, Stool_mean_dmi,Skin_mean_dmi,Oral_mean_dmi,Nasal_mean_dmi)
ggradar(radar.mean, grid.min = 0, grid.max =  0.4) + scale_color_manual(values = phylum_color)

radar.mean.t <- rbind(radar.mean$Stool_mean_dmi, radar.mean$Skin_mean_dmi, radar.mean$Oral_mean_dmi, radar.mean$Nasal_mean_dmi)
colnames(radar.mean.t) <- radar.mean$Phylum
rownames(radar.mean.t) <- c("Stool", "Skin", "Oral", "Nasal")
radar.mean.t <- as.data.frame(radar.mean.t)

radar.mean.t$BodySite <- c("Stool", "Skin", "Oral", "Nasal")
radar.mean.t <- radar.mean.t %>% dplyr::select(BodySite, everything())

radar.by.bodysite <- ggradar(radar.mean.t, grid.min = 0, grid.max = 0.4, values.radar = c("0", "0.4")) + scale_color_manual(values = body_site_color)
radar.by.bodysite
#ggsave2(filename = "~/Desktop/1st.revision.CELL-S-22-04878/Figures/Figure2/Radar.DMI.Bodysite.pdf",radar.by.bodysite, height = 9, width = 7, dpi = 300)

###
stool.trim <- stool.data.frame.trimonly[, c("Phylum", "dmi")]

stool.pairs <- combn(unique(stool.trim$Phylum), 2, simplify = FALSE)

stool.p_values <- lapply(stool.pairs, function(stool.pairs) {
  test_result <- wilcox.test(stool.trim$dmi[stool.trim$Phylum == stool.pairs[1]],
                             stool.trim$dmi[stool.trim$Phylum == stool.pairs[2]])
  return(data.frame(Comparison = paste(stool.pairs[1], stool.pairs[2], sep = "_"),
                    P_Value = test_result$p.value,
                    W = test_result$statistic))
})

stool.results_table <- do.call(rbind, stool.p_values)
stool.results_table$p.adj <- p.adjust(stool.results_table$P_Value,method = "BH", n = length(stool.results_table$P_Value))
stool.results_table$bodysite <- "Stool"
stool.results_table

#skin
skin.trim <- skin.data.frame.trimonly[, c("Phylum", "dmi")]

skin.pairs <- combn(unique(skin.trim$Phylum), 2, simplify = FALSE)

skin.p_values <- lapply(skin.pairs, function(skin.pairs) {
  test_result <- wilcox.test(skin.trim$dmi[skin.trim$Phylum == skin.pairs[1]],
                             skin.trim$dmi[skin.trim$Phylum == skin.pairs[2]])
  return(data.frame(Comparison = paste(skin.pairs[1], skin.pairs[2], sep = "_"),
                    P_Value = test_result$p.value,
                    W = test_result$statistic))
})

skin.results_table <- do.call(rbind, skin.p_values)
skin.results_table$p.adj <- p.adjust(skin.results_table$P_Value,method = "BH", n = length(skin.results_table$P_Value))
skin.results_table$bodysite <- "Skin"
skin.results_table

#oral
oral.trim <- oral.data.frame.trimonly[, c("Phylum", "dmi")]

oral.pairs <- combn(unique(oral.trim$Phylum), 2, simplify = FALSE)

oral.p_values <- lapply(oral.pairs, function(oral.pairs) {
  test_result <- wilcox.test(oral.trim$dmi[oral.trim$Phylum == oral.pairs[1]],
                             oral.trim$dmi[oral.trim$Phylum == oral.pairs[2]])
  return(data.frame(Comparison = paste(oral.pairs[1], oral.pairs[2], sep = "_"),
                    P_Value = test_result$p.value,
                    W = test_result$statistic))
})

oral.results_table <- do.call(rbind, oral.p_values)
oral.results_table$p.adj <- p.adjust(oral.results_table$P_Value,method = "BH", n = length(oral.results_table$P_Value))
oral.results_table$bodysite <- "Oral"
oral.results_table

#nasal
nasal.trim <- nasal.data.frame.trimonly[, c("Phylum", "dmi")]

nasal.pairs <- combn(unique(nasal.trim$Phylum), 2, simplify = FALSE)

nasal.p_values <- lapply(nasal.pairs, function(nasal.pairs) {
  test_result <- wilcox.test(nasal.trim$dmi[nasal.trim$Phylum == nasal.pairs[1]],
                             nasal.trim$dmi[nasal.trim$Phylum == nasal.pairs[2]])
  return(data.frame(Comparison = paste(nasal.pairs[1], nasal.pairs[2], sep = "_"),
                    P_Value = test_result$p.value,
                    W = test_result$statistic))
})

nasal.results_table <- do.call(rbind, nasal.p_values)
nasal.results_table$p.adj <- p.adjust(nasal.results_table$P_Value, method = "BH", n = length(nasal.results_table$P_Value))
nasal.results_table$bodysite <- "Nasal"
nasal.results_table

dmi_compare <- rbind(stool.results_table,skin.results_table,oral.results_table,nasal.results_table)
dmi_compare
#write.csv(file = "~/Desktop/1st.revision.CELL-S-22-04878/Supplementary_data/Supplimentaryfigure_added/Figure2.DMI_Phylum.csv",dmi_compare)


#https://r-charts.com/ranking/ggradar/ 
stool.radar <- stool.trim %>% dplyr::select(Phylum,dmi) %>% mutate(bodysite="Stool")
skin.radar <- skin.trim %>% dplyr::select(Phylum,dmi)%>% mutate(bodysite="Skin")
oral.radar <- oral.trim %>% dplyr::select(Phylum,dmi)%>% mutate(bodysite="Oral")
nasal.radar <- nasal.trim %>% dplyr::select(Phylum,dmi) %>% mutate(bodysite="Nasal")


radar.test.df <- rbind(stool.radar,skin.radar,oral.radar,nasal.radar)
radar.test.df

grouped_df_radar <- radar.test.df %>%
  group_by(Phylum)
grouped_df_radar

# run test
wli.result.byphyla <- grouped_df_radar %>%
  do(test = pairwise.wilcox.test(.$dmi, .$bodysite, p.adjust.method = "BH"))
wli.result.byphyla
wli.result.byphyla$test

# Define body sites
body_sites <- c("Nasal", "Oral", "Skin", "Stool")

# Manually extract p-values from each test
p_values_list <- list(
  actinobacteria = matrix(c(0.671, 1.000, 0.014, 0.671, 0.095, 0.107), nrow = 3, byrow = TRUE),
  bacteroidetes = matrix(c(0.00857, 0.34286, 0.00568, 0.14082, 0.00018, 0.00568), nrow = 3, byrow = TRUE),
  firmicutes = matrix(c(0.1312, 0.4121, 0.8975, 0.8975, 0.0017, 0.2162), nrow = 3, byrow = TRUE),
  proteobacteria = matrix(c(0.2777, 0.8726, 0.3164, 0.7812, 0.0096, 0.3952), nrow = 3, byrow = TRUE),
  other = matrix(c(0.59, 0.59, 0.59, 0.79, 0.59, 0.70), nrow = 3, byrow = TRUE)
)

# Combine p-values into a single data frame
p_values_df <- bind_rows(
  mutate(as.data.frame.table(p_values_list$actinobacteria), Phylum = "Actinobacteria"),
  mutate(as.data.frame.table(p_values_list$bacteroidetes), Phylum = "Bacteroidetes"),
  mutate(as.data.frame.table(p_values_list$firmicutes), Phylum = "Firmicutes"),
  mutate(as.data.frame.table(p_values_list$proteobacteria), Phylum = "Proteobacteria"),
  mutate(as.data.frame.table(p_values_list$other), Phylum = "Other")
)

# Add body site information
p_values_df$BodySite1 <- body_sites[p_values_df$Var1]
p_values_df$BodySite2 <- body_sites[p_values_df$Var2]

# Select and reorder columns
p_values_df <- select(p_values_df, Phylum, BodySite1, BodySite2, Freq)

# Rename columns
colnames(p_values_df) <- c("Phylum", "BodySite1", "BodySite2", "P_value")

# Print the table
print(p_values_df)
#write.csv(file = "~/Desktop/1st.revision.CELL-S-22-04878/Supplementary_data/Supplimentaryfigure_added/Figure2.Radar.byphylum.csv",p_values_df)

phyla <- unique(radar.test.df$Phylum)
results.f2c <- list()

for (phylum in phyla) {
  subset_data <- radar.test.df[radar.test.df$Phylum == phylum,]
  kruskal_result <- kruskal.test(dmi ~ bodysite, data = subset_data)
  results.f2c[[phylum]] <- kruskal_result
}

# Print results
results.f2c

stool.data.frame
#output_new_DMI table as supplement
sumpliment.dmi.table <-rbind(stool.data.frame %>% select(genus, dmi, Phylum, bodysite,fc1, fc1_sd,trimmed) %>% drop_na(),
      skin.data.frame %>% dplyr::select(genus, dmi, Phylum, bodysite,fc1, fc1_sd,trimmed) %>% drop_na(),
      oral.data.frame %>% dplyr::select(genus, dmi, Phylum, bodysite,fc1, fc1_sd,trimmed) %>% drop_na(),
      nasal.data.frame %>% dplyr::select(genus, dmi, Phylum, bodysite,fc1, fc1_sd,trimmed) %>% drop_na())
# write.csv(file = "~/Desktop/1st.revision.CELL-S-22-04878/Supplementary_data/Supplimentaryfigure_added/New_DMI_TABLE.csv", sumpliment.dmi.table)


#play with stool data
p.stool.dmi.Shannon <- filter(stool.data.frame,trimmed != "NA") %>%  ggplot(aes(x=dmi, y=Shannon)) + geom_point() + geom_smooth(method = "lm")
p.stool.dmi.Shannon <- p.stool.dmi.Shannon + facet_wrap(.~trimmed) + theme_cowplot() + ggtitle("Stool")
p.stool.dmi.Shannon

p.stool.dmi.pre <- filter(stool.data.frame,trimmed != "NA") %>%  ggplot(aes(x=dmi, y=Prevalance)) + geom_point() + geom_smooth(method = "lm")
p.stool.dmi.pre <- p.stool.dmi.pre + facet_wrap(.~trimmed) + theme_cowplot() + ggtitle("Stool")
p.stool.dmi.pre

p.stool.dmi.ACE <- filter(stool.data.frame,trimmed != "NA") %>%  ggplot(aes(x=dmi, y=ACE)) + geom_point() + geom_smooth(method = "lm")
p.stool.dmi.ACE <- p.stool.dmi.ACE + facet_wrap(.~trimmed) + theme_cowplot() + ggtitle("Stool")
p.stool.dmi.ACE

p.stool.Shannon.Pre <- filter(stool.data.frame,trimmed != "NA") %>%  ggplot(aes(x=Shannon, y=Prevalance, color=trimmed)) + geom_point() + geom_smooth(method = "lm")
p.stool.Shannon.Pre <- p.stool.Shannon.Pre + theme_cowplot() + ggtitle("Stool")
p.stool.Shannon.Pre

#Prevalence and DMI
p.stool.dmi.pre <- filter(stool.data.frame,trimmed != "NA") %>%  ggplot(aes(x=dmi, y=Prevalance)) + geom_point() + geom_smooth(method = "lm")
p.stool.dmi.pre <- p.stool.dmi.pre + facet_wrap(.~trimmed) + theme_cowplot() + ggtitle("Stool")
p.stool.dmi.pre

p.skin.dmi.pre <- filter(skin.data.frame,trimmed != "NA") %>%  ggplot(aes(x=dmi, y=Prevalance)) + geom_point() + geom_smooth(method = "lm")
p.skin.dmi.pre <- p.skin.dmi.pre + facet_wrap(.~trimmed) + theme_cowplot() + ggtitle("Skin")
p.skin.dmi.pre

p.oral.dmi.pre <- filter(oral.data.frame,trimmed != "NA") %>%  ggplot(aes(x=dmi, y=Prevalance)) + geom_point() + geom_smooth(method = "lm")
p.oral.dmi.pre <- p.oral.dmi.pre + facet_wrap(.~trimmed) + theme_cowplot() + ggtitle("Oral")
p.oral.dmi.pre

p.nasal.dmi.pre <- filter(nasal.data.frame,trimmed != "NA") %>%  ggplot(aes(x=dmi, y=Prevalance)) + geom_point() + geom_smooth(method = "lm")
p.nasal.dmi.pre <- p.nasal.dmi.pre + facet_wrap(.~trimmed) + theme_cowplot() + ggtitle("Nasal")
p.nasal.dmi.pre

#First figure
p.stool.pre.Shannon <- filter(stool.data.frame, trimmed != "NA") %>%  
  filter(Shannon != 0) %>% 
  ggplot(aes(x = Prevalance, y = Shannon)) + 
  geom_point() + 
  geom_smooth(method = "lm") +
  theme_cowplot() + 
  ggtitle("Stool \n Relationship between Prevalence and Shannon") + 
  xlab("Prevalence") + 
  ylab("Shannon") +
  stat_cor(label.x = 0.25, label.y = max(stool.data.frame$Shannon), label.sep = " ; ", size = 5, method = "pearson", 
           show.legend = TRUE)

p.stool.pre.Shannon

p.skin.pre.Shannon <- filter(skin.data.frame, trimmed != "NA") %>%  
  filter(Shannon != 0) %>% 
  ggplot(aes(x = Prevalance, y = Shannon)) + 
  geom_point() + 
  geom_smooth(method = "lm") +
  theme_cowplot() + 
  ggtitle("Skin \n Relationship between Prevalence and Shannon") + 
  xlab("Prevalence") + 
  ylab("Shannon") +
  stat_cor(label.x = 0.25, label.y = max(skin.data.frame$Shannon), label.sep = " ; ", size = 5, method = "pearson", 
           show.legend = TRUE)

p.skin.pre.Shannon

p.oral.pre.Shannon <- filter(oral.data.frame, trimmed != "NA") %>%  
  filter(Shannon != 0) %>% 
  ggplot(aes(x = Prevalance, y = Shannon)) + 
  geom_point() + 
  geom_smooth(method = "lm") +
  theme_cowplot() + 
  ggtitle("Oral \n Relationship between Prevalence and Shannon") + 
  xlab("Prevalence") + 
  ylab("Shannon") +
  stat_cor(label.x = 0.25, label.y = max(oral.data.frame$Shannon), label.sep = " ; ", size = 5, method = "pearson", 
           show.legend = TRUE)

p.oral.pre.Shannon

p.nasal.pre.Shannon <- filter(nasal.data.frame, trimmed != "NA") %>%  
  filter(Shannon != 0) %>% 
  ggplot(aes(x = Prevalance, y = Shannon)) + 
  geom_point() + 
  geom_smooth(method = "lm") +
  theme_cowplot() + 
  ggtitle("Nasal \n Relationship between Prevalence and Shannon") + 
  xlab("Prevalence") + 
  ylab("Shannon") +
  stat_cor(label.x = 0.25, label.y = max(nasal.data.frame$Shannon), label.sep = " ; ", size = 5, method = "pearson", 
           show.legend = TRUE)

p.nasal.pre.Shannon

p.stool.pre.Shannon + p.skin.pre.Shannon + p.oral.pre.Shannon + p.nasal.pre.Shannon


#second figure

p.stool.dmi.Shannon <- filter(stool.data.frame, trimmed == "no") %>%  
  ggplot(aes(x = dmi, y = Shannon)) + 
  geom_point() + 
  geom_smooth(method = "lm") +
  theme_cowplot() + 
  ggtitle("Stool \n Relationship between DMI and Shannon") + 
  stat_cor(label.x = 0.2, label.y = max(stool.data.frame$Shannon), label.sep = " ; ", size = 5, method = "pearson", 
           show.legend = TRUE)

p.stool.dmi.Shannon

p.skin.dmi.Shannon <- filter(skin.data.frame, trimmed == "no") %>%  
  ggplot(aes(x = dmi, y = Shannon)) + 
  geom_point() + 
  geom_smooth(method = "lm") +
  theme_cowplot() + 
  ggtitle("Skin \n Relationship between DMI and Shannon") + 
  stat_cor(label.x = 0.2, label.y = 0.25, label.sep = " ; ", size = 5, method = "pearson", 
           show.legend = TRUE)

p.skin.dmi.Shannon

p.oral.dmi.Shannon <- filter(oral.data.frame, trimmed == "no") %>%  
  ggplot(aes(x = dmi, y = Shannon)) + 
  geom_point() + 
  geom_smooth(method = "lm") +
  theme_cowplot() + 
  ggtitle("Oral \n Relationship between DMI and Shannon") + 
  stat_cor(label.x = 0.2, label.y = 1, label.sep = " ; ", size = 5, method = "pearson", 
           show.legend = TRUE)

p.oral.dmi.Shannon

p.nasal.dmi.Shannon <- filter(nasal.data.frame, trimmed == "no") %>%  
  ggplot(aes(x = dmi, y = Shannon)) + 
  geom_point() + 
  geom_smooth(method = "lm") +
  theme_cowplot() + 
  ggtitle("Nasal \n Relationship between DMI and Shannon") + 
  stat_cor(label.x = 0.2, label.y = max(nasal.data.frame$Shannon), label.sep = " ; ", size = 5, method = "pearson", 
           show.legend = TRUE)

p.nasal.dmi.Shannon

p.stool.dmi.Shannon + p.skin.dmi.Shannon + p.oral.dmi.Shannon + p.nasal.dmi.Shannon

###########
#stats for figure 2B
combined.dmi <- rbind(filter(stool.data.frame,trimmed == "no") %>% dplyr::select(class,dmi),
                  filter(skin.data.frame,trimmed == "no") %>% dplyr::select(class,dmi),
                  filter(oral.data.frame,trimmed == "no") %>% dplyr::select(class,dmi),
                  filter(nasal.data.frame,trimmed == "no") %>% dplyr::select(class,dmi))

# Perform pairwise Wilcoxon-Mann-Whitney tests
pairwise_wilcox_results.dmi <- pairwise.wilcox.test(combined.dmi$dmi, combined.dmi$class,
                                                p.adjust.method = "BH")

# Print the results
print(pairwise_wilcox_results.dmi)

#remake FS figure
FS_score <- read.csv("./human_microbiome_project/Supplementary_data/Single tables/Extended Data Table FS.csv", header = T)
FS_score

colnames(FS_score)[2] <- "bodysite"

colnames(FS_score)
colnames(sumpliment.dmi.table)
merged_table.figure2 <- left_join(FS_score, sumpliment.dmi.table, by = c("genus", "Phylum", "bodysite"))

p.hist.dmi <- ggplot(merged_table.figure2, aes(x = dmi, fill = bodysite)) +
  geom_histogram(alpha = 0.5, binwidth = 0.02) +
  geom_density(alpha = 0.5) +
  facet_wrap(~bodysite, ncol = 2) +
  xlab("dmi") +
  ylab("Frequency") +
  ggtitle("Distribution of dmi by bodysite") +
  scale_fill_manual(values = body_site_color) + base_theme + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p.hist.dmi

p.hist.FS <- ggplot(merged_table.figure2, aes(x = family_score, fill = bodysite)) +
  geom_histogram(alpha = 0.5, binwidth = 0.02) +
  geom_density(alpha = 0.5) +
  facet_wrap(~bodysite, ncol = 2, scales = "free") +
  xlab("Family Score") +
  ylab("Frequency") +
  ggtitle("Distribution of Family Score by bodysite") +
  scale_fill_manual(values = body_site_color) + base_theme + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p.hist.FS

p.dmi_fs <- p.hist.dmi + p.hist.FS
p.dmi_fs  


table(merged_table.figure2$bodysite)

result.fs <- kruskal.test(family_score ~ bodysite, data = merged_table.figure2)
print(result.fs)

result.fs.posthoc <- suppressWarnings(pairwise.wilcox.test(merged_table.figure2$family_score, merged_table.figure2$bodysite, p.adjust.method = "BH"))
result.fs.posthoc

result.dmi.posthoc <- suppressWarnings(pairwise.wilcox.test(merged_table.figure2$dmi, merged_table.figure2$bodysite, p.adjust.method = "BH"))
result.dmi.posthoc



#ggsave2(filename = "~/Desktop/1st.revision.CELL-S-22-04878/Supplementary_data/FS_DMI_HIST.pdf", p.dmi_fs, width = 10, height = 6, dpi=300) 
#write.csv(file = "~/Desktop/1st.revision.CELL-S-22-04878/Supplementary_data/Single tables/New_DMI_FS.csv",merged_table.figure2)

#see who affect cytokine more  
Cytokine.cor <- read.csv("./human_microbiome_project/Supplementary_data/Single tables/cytokinemicrobiomecor.csv", header = T)
Cytokine.cor

Cytokine.cor.stool <- filter(Cytokine.cor,bodysite == "Stool")
Cytokine.cor.skin <- filter(Cytokine.cor, bodysite == "Skin")
Cytokine.cor.oral <- filter(Cytokine.cor, bodysite == "Oral")
Cytokine.cor.nasal <- filter(Cytokine.cor, bodysite == "Nasal")

stool.data.frame <- stool.data.frame %>%
  mutate(cytokine.engagement = ifelse(genus %in% Cytokine.cor.stool$genus, "detected", "not_detected"))

stool.data.frame <- stool.data.frame %>% 
  mutate(cytokine.engagement = if_else(!is.na(trimmed) & trimmed == "yes", "not_detected", cytokine.engagement))

skin.data.frame <- skin.data.frame %>%
  mutate(cytokine.engagement = ifelse(genus %in% Cytokine.cor.skin$genus, "detected", "not_detected"))

skin.data.frame <- skin.data.frame %>%
  mutate(cytokine.engagement = if_else(!is.na(trimmed) & trimmed == "yes", "not_detected", cytokine.engagement))

oral.data.frame <- oral.data.frame %>%
  mutate(cytokine.engagement = ifelse(genus %in% Cytokine.cor.oral$genus, "detected", "not_detected"))

oral.data.frame <- oral.data.frame %>%
  mutate(cytokine.engagement = if_else(!is.na(trimmed) & trimmed == "yes", "not_detected", cytokine.engagement))

nasal.data.frame <- nasal.data.frame %>%
  mutate(cytokine.engagement = ifelse(genus %in% Cytokine.cor.nasal$genus, "detected", "not_detected"))

nasal.data.frame <- nasal.data.frame %>%
  mutate(cytokine.engagement = if_else(!is.na(trimmed) & trimmed == "yes", "not_detected", cytokine.engagement))


# Filter out rows with NA values and add order to the Phylum column
stool.filtered <- stool.data.frame %>%
  filter(!is.na(dmi)) %>%
  mutate(Phylum = factor(Phylum, levels = phylum_order)) %>%
  left_join(stool.data.frame %>%
              filter(!is.na(dmi)) %>%
              group_by(Phylum) %>%
              summarise(pvalue = t.test(dmi ~ cytokine.engagement, data = .)$p.value),
            by = "Phylum")

skin.filtered <- skin.data.frame %>%
  filter(!is.na(dmi)) %>%
  mutate(Phylum = factor(Phylum, levels = phylum_order)) %>%
  left_join(skin.data.frame %>%
              filter(!is.na(dmi)) %>%
              group_by(Phylum) %>%
              summarise(pvalue = t.test(dmi ~ cytokine.engagement, data = .)$p.value),
            by = "Phylum")

oral.filtered <- oral.data.frame %>%
  filter(!is.na(dmi)) %>%
  mutate(Phylum = factor(Phylum, levels = phylum_order)) %>%
  left_join(oral.data.frame %>%
              filter(!is.na(dmi)) %>%
              group_by(Phylum) %>%
              summarise(pvalue = t.test(dmi ~ cytokine.engagement, data = .)$p.value),
            by = "Phylum")

nasal.filtered <- nasal.data.frame %>%
  filter(!is.na(dmi)) %>%
  mutate(Phylum = factor(Phylum, levels = phylum_order)) %>%
  left_join(nasal.data.frame %>%
              filter(!is.na(dmi)) %>%
              group_by(Phylum) %>%
              summarise(pvalue = t.test(dmi ~ cytokine.engagement, data = .)$p.value),
            by = "Phylum")

# Calculate the maximum value of dmi in the stool.filtered dataframe
max_dmi_stool <- max(stool.filtered$dmi, na.rm = TRUE) + 0.1
max_dmi_skin <- max(skin.filtered$dmi, na.rm = TRUE) + 0.1
max_dmi_oral <- max(oral.filtered$dmi, na.rm = TRUE) + 0.1
max_dmi_nasal <- max(nasal.filtered$dmi, na.rm = TRUE) + 0.1

# Create a box plot comparing the DMI values based on cytokine engagement
stool_plot <- ggplot(stool.filtered, aes(x = cytokine.engagement, y = dmi, color = cytokine.engagement)) +
  geom_boxplot() +
  geom_jitter() +
  facet_wrap(~Phylum) +
  ggtitle("Comparison of DMI based on cytokine engagement by Phylum") +
  xlab("Cytokine engagement") +
  ylab("DMI")

skin_plot <- ggplot(skin.filtered, aes(x = cytokine.engagement, y = dmi, color = cytokine.engagement)) +
  geom_boxplot() +
  geom_jitter() +
  facet_wrap(~Phylum) +
  ggtitle("Comparison of DMI based on cytokine engagement by Phylum in skin dataset") +
  xlab("Cytokine engagement") +
  ylab("DMI")

oral_plot <- ggplot(oral.filtered, aes(x = cytokine.engagement, y = dmi, color = cytokine.engagement)) +
  geom_boxplot() +
  geom_jitter() +
  facet_wrap(~Phylum) +
  ggtitle("Comparison of DMI based on cytokine engagement by Phylum in oral dataset") +
  xlab("Cytokine engagement") +
  ylab("DMI")

nasal_plot <- ggplot(nasal.filtered, aes(x = cytokine.engagement, y = dmi, color = cytokine.engagement)) +
  geom_boxplot() +
  geom_jitter() +
  facet_wrap(~Phylum) +
  ggtitle("Comparison of DMI based on cytokine engagement by Phylum in nasal dataset") +
  xlab("Cytokine engagement") +
  ylab("DMI")

# Add the p-values to the plot

stool_plot.s <- stool_plot + stat_compare_means(aes(label = paste("p = ", round(pvalue, 3), sep="")), size = 4,
                                                label.y = max_dmi_stool, label.x = c(1,2),
                                                method = "t.test", paired = FALSE) + theme_pubr()
stool_plot.s

skin_plot.s <- skin_plot + stat_compare_means(aes(label = paste("p = ", round(pvalue, 3), sep="")), size = 4,
                                                label.y = max_dmi_stool, label.x = c(1,2),
                                                method = "t.test", paired = FALSE) + theme_pubr()
skin_plot.s

oral_plot.s <- oral_plot + stat_compare_means(aes(label = paste("p = ", round(pvalue, 3), sep="")), size = 4,
                                                label.y = max_dmi_stool, label.x = c(1,2),
                                                method = "t.test", paired = FALSE) + theme_pubr()
oral_plot.s

nasal_plot.s <- nasal_plot + stat_compare_means(aes(label = paste("p = ", round(pvalue, 3), sep="")), size = 4,
                                                label.y = max_dmi_stool, label.x = c(1,2),
                                                method = "t.test", paired = FALSE) + theme_pubr()
nasal_plot.s

# 
# ggsave2(filename = "~/Desktop/1st.revision.CELL-S-22-04878/Figures/Figure4/detailed/stool.dmi.cytokine.pdf", stool_plot.s, width = 8, height = 8, dpi = 300)
# ggsave2(filename = "~/Desktop/1st.revision.CELL-S-22-04878/Figures/Figure4/detailed/skin.dmi.cytokine.pdf", skin_plot.s, width = 8, height = 8, dpi = 300)
# ggsave2(filename = "~/Desktop/1st.revision.CELL-S-22-04878/Figures/Figure4/detailed/oral.dmi.cytokine.pdf", oral_plot.s, width = 8, height = 8, dpi = 300)
# ggsave2(filename = "~/Desktop/1st.revision.CELL-S-22-04878/Figures/Figure4/detailed/nasal.dmi.cytokine.pdf", nasal_plot.s, width = 8, height = 8, dpi = 300)

# Define the order of the Phylum column
phylum_order <- c("Actinobacteria", "Bacteroidetes", "Firmicutes", "Proteobacteria", "Other")

#Stool T test
ttest_actino <- t.test(dmi ~ cytokine.engagement, data = filter(stool.filtered, Phylum == "Actinobacteria"))
ttest_bacter <- t.test(dmi ~ cytokine.engagement, data = filter(stool.filtered, Phylum == "Bacteroidetes"))
ttest_firmic <- t.test(dmi ~ cytokine.engagement, data = filter(stool.filtered, Phylum == "Firmicutes"))
#ttest_other <- t.test(dmi ~ cytokine.engagement, data = filter(stool.filtered, Phylum == "Other"))
ttest_proteo <- t.test(dmi ~ cytokine.engagement, data = filter(stool.filtered, Phylum == "Proteobacteria"))
cat("Actinobacteria:\n")
print(ttest_actino$statistic)
print(ttest_actino$p.value)

cat("\nBacteroidetes:\n")
print(ttest_bacter$statistic)
print(ttest_bacter$p.value)

cat("\nFirmicutes:\n")
print(ttest_firmic$statistic)
print(ttest_firmic$p.value)

# cat("\nOther:\n")
# print(ttest_other$statistic)
# print(ttest_other$p.value)

cat("\nProteobacteria:\n")
print(ttest_proteo$statistic)
print(ttest_proteo$p.value)

#Skin T Test
ttest_actino <- t.test(dmi ~ cytokine.engagement, data = filter(skin.filtered, Phylum == "Actinobacteria"))
ttest_bacter <- t.test(dmi ~ cytokine.engagement, data = filter(skin.filtered, Phylum == "Bacteroidetes"))
ttest_firmic <- t.test(dmi ~ cytokine.engagement, data = filter(skin.filtered, Phylum == "Firmicutes"))
ttest_other <- t.test(dmi ~ cytokine.engagement, data = filter(skin.filtered, Phylum == "Other"))
ttest_proteo <- t.test(dmi ~ cytokine.engagement, data = filter(skin.filtered, Phylum == "Proteobacteria"))
cat("Actinobacteria:\n")
print(ttest_actino$statistic)
print(ttest_actino$p.value)

cat("\nBacteroidetes:\n")
print(ttest_bacter$statistic)
print(ttest_bacter$p.value)

cat("\nFirmicutes:\n")
print(ttest_firmic$statistic)
print(ttest_firmic$p.value)

cat("\nOther:\n")
print(ttest_other$statistic)
print(ttest_other$p.value)

cat("\nProteobacteria:\n")
print(ttest_proteo$statistic)
print(ttest_proteo$p.value)

#Oral
ttest_actino <- t.test(dmi ~ cytokine.engagement, data = filter(oral.filtered, Phylum == "Actinobacteria"))
ttest_bacter <- t.test(dmi ~ cytokine.engagement, data = filter(oral.filtered, Phylum == "Bacteroidetes"))
ttest_firmic <- t.test(dmi ~ cytokine.engagement, data = filter(oral.filtered, Phylum == "Firmicutes"))
ttest_other <- t.test(dmi ~ cytokine.engagement, data = filter(oral.filtered, Phylum == "Other"))
ttest_proteo <- t.test(dmi ~ cytokine.engagement, data = filter(oral.filtered, Phylum == "Proteobacteria"))
cat("Actinobacteria:\n")
print(ttest_actino$statistic)
print(ttest_actino$p.value)

cat("\nBacteroidetes:\n")
print(ttest_bacter$statistic)
print(ttest_bacter$p.value)

cat("\nFirmicutes:\n")
print(ttest_firmic$statistic)
print(ttest_firmic$p.value)

cat("\nOther:\n")
print(ttest_other$statistic)
print(ttest_other$p.value)

cat("\nProteobacteria:\n")
print(ttest_proteo$statistic)
print(ttest_proteo$p.value)

#Nasal
ttest_actino <- t.test(dmi ~ cytokine.engagement, data = filter(nasal.filtered, Phylum == "Actinobacteria"))
ttest_bacter <- t.test(dmi ~ cytokine.engagement, data = filter(nasal.filtered, Phylum == "Bacteroidetes"))
ttest_firmic <- t.test(dmi ~ cytokine.engagement, data = filter(nasal.filtered, Phylum == "Firmicutes"))
ttest_other <- t.test(dmi ~ cytokine.engagement, data = filter(nasal.filtered, Phylum == "Other"))
ttest_proteo <- t.test(dmi ~ cytokine.engagement, data = filter(nasal.filtered, Phylum == "Proteobacteria"))
cat("Actinobacteria:\n")
print(ttest_actino$statistic)
print(ttest_actino$p.value)

cat("\nBacteroidetes:\n")
print(ttest_bacter$statistic)
print(ttest_bacter$p.value)

cat("\nFirmicutes:\n")
print(ttest_firmic$statistic)
print(ttest_firmic$p.value)

cat("\nOther:\n")
print(ttest_other$statistic)
print(ttest_other$p.value)

cat("\nProteobacteria:\n")
print(ttest_proteo$statistic)
print(ttest_proteo$p.value)


summary_table <- data.frame(
  Phylum = c("Actinobacteria", "Bacteroidetes", "Firmicutes", "Other", "Proteobacteria"),
  Stool_t = c(-0.3931565, -0.6855807, -2.767941, NA, 2.108168),
  Stool_p = c(0.7077904, 0.5023091, 0.007271335, NA, 0.05470111),
  Skin_t = c(-3.548151, -2.461974, -3.11193, -3.655409, -5.615827),
  Skin_p = c(0.001183958, 0.02256107, 0.003092849, 0.003665233, 9.063625e-07),
  Oral_t = c(-1.59239, -1.73074, -3.052313, -0.9177608, -0.4434136),
  Oral_p = c(0.1637553, 0.1689681, 0.005966473, 0.3971465, 0.6677908),
  Nasal_t = c(-0.5277027, 0.4594002, -1.346606, -0.6564743, -0.216138),
  Nasal_p = c(0.6016988, 0.6498519, 0.1856893, 0.5230605, 0.830874)
)


summary_table


# Create a contingency table for figure 4A
contingency_table_cytokine <- matrix(c(477, 6572 - 477,
                              226, 4092 - 226,
                              318, 3968 - 318,
                              221, 4092 - 221),
                            nrow = 4, byrow = TRUE,
                            dimnames = list(BodySite = c("Stool", "Skin", "Oral", "Nasal"),
                                            Result = c("Significant", "Non-significant")))
contingency_table_cytokine

# Perform the Chi-squared test
chi_squared_test_cytokine4A <- chisq.test(contingency_table_cytokine)
chi_squared_test_cytokine4A

# Function to perform pairwise Chi-squared tests
pairwise_chi_squared <- function(contingency_table_cytokine, i, j) {
  pairwise_table <- contingency_table_cytokine[c(i, j), ]
  test_result <- chisq.test(pairwise_table)
  return(test_result)
}

# Perform pairwise comparisons
body_sites <- rownames(contingency_table_cytokine)
combinations <- combn(length(body_sites), 2)
pairwise_results <- list()

for (i in 1:ncol(combinations)) {
  site1 <- body_sites[combinations[1, i]]
  site2 <- body_sites[combinations[2, i]]
  pairwise_test <- pairwise_chi_squared(contingency_table_cytokine, combinations[1, i], combinations[2, i])
  pairwise_results[[paste(site1, "vs", site2)]] <- pairwise_test
}

# Print pairwise results
pairwise_results


#####################################################################################################################################
#Figure5
#####################################################################################################################################
metab.df <- read.csv("./human_microbiome_project/Supplementary_data/Single tables/Extended Data Table MSandMicro.csv", header = T)
metab.stool <- filter(metab.df, BodySite == "Stool")
metab.skin <- filter(metab.df, BodySite == "Skin")
metab.oral <- filter(metab.df, BodySite == "Oral")
metab.nasal <- filter(metab.df, BodySite == "Nasal")

#stool
metab.stool_relevant <- metab.stool %>%
  filter(Class == "Metabolome") %>%
  select(from_true_name = from_true_name, Class)

lipid.stool_relevant <- metab.stool %>%
  filter(Class == "Lipidome") %>%
  select(from_true_name = from_true_name, Class)

prot.stool_relevant <- metab.stool %>%
  filter(Class == "Proteome") %>%
  select(from_true_name = from_true_name, Class)

stool.data.frame.m <- stool.data.frame %>%
  mutate(metab.engagement = ifelse(genus %in% metab.stool_relevant$from_true_name, "detected", "not_detected")) %>%
  mutate(lipid.engagement = ifelse(genus %in% lipid.stool_relevant$from_true_name, "detected", "not_detected")) %>%
  mutate(prot.engagement = ifelse(genus %in% prot.stool_relevant$from_true_name, "detected", "not_detected")) %>%
  mutate(trim = ifelse(genus %in% stool.trim1$genus, "keep", "trimed"))

table(stool.data.frame.m$trim)

#skin
metab.skin_relevant <- metab.skin %>%
  filter(Class == "Metabolome") %>%
  select(from_true_name = from_true_name, Class)

lipid.skin_relevant <- metab.skin %>%
  filter(Class == "Lipidome") %>%
  select(from_true_name = from_true_name, Class)

prot.skin_relevant <- metab.skin %>%
  filter(Class == "Proteome") %>%
  select(from_true_name = from_true_name, Class)

skin.data.frame.m <- skin.data.frame %>%
  mutate(metab.engagement = ifelse(genus %in% metab.skin_relevant$from_true_name, "detected", "not_detected")) %>%
  mutate(lipid.engagement = ifelse(genus %in% lipid.skin_relevant$from_true_name, "detected", "not_detected")) %>%
  mutate(prot.engagement = ifelse(genus %in% prot.skin_relevant$from_true_name, "detected", "not_detected")) %>%
  mutate(trim = ifelse(genus %in% skin.trim1$genus, "keep", "trimmed"))

table(skin.data.frame.m$trim)

#oral
metab.oral_relevant <- metab.oral %>%
  filter(Class == "Metabolome") %>%
  select(from_true_name = from_true_name, Class)

lipid.oral_relevant <- metab.oral %>%
  filter(Class == "Lipidome") %>%
  select(from_true_name = from_true_name, Class)

prot.oral_relevant <- metab.oral %>%
  filter(Class == "Proteome") %>%
  select(from_true_name = from_true_name, Class)

oral.data.frame.m <- oral.data.frame %>%
  mutate(metab.engagement = ifelse(genus %in% metab.oral_relevant$from_true_name, "detected", "not_detected")) %>%
  mutate(lipid.engagement = ifelse(genus %in% lipid.oral_relevant$from_true_name, "detected", "not_detected")) %>%
  mutate(prot.engagement = ifelse(genus %in% prot.oral_relevant$from_true_name, "detected", "not_detected")) %>%
  mutate(trim = ifelse(genus %in% oral.trim1$genus, "keep", "trimed"))

#nasal
metab.nasal_relevant <- metab.nasal %>%
  filter(Class == "Metabolome") %>%
  select(from_true_name = from_true_name, Class)

lipid.nasal_relevant <- metab.nasal %>%
  filter(Class == "Lipidome") %>%
  select(from_true_name = from_true_name, Class)

prot.nasal_relevant <- metab.nasal %>%
  filter(Class == "Proteome") %>%
  select(from_true_name = from_true_name, Class)

nasal.data.frame.m <- nasal.data.frame %>%
  mutate(metab.engagement = ifelse(genus %in% metab.nasal_relevant$from_true_name, "detected", "not_detected"))%>%
  mutate(lipid.engagement = ifelse(genus %in% lipid.nasal_relevant$from_true_name, "detected", "not_detected")) %>%
  mutate(prot.engagement = ifelse(genus %in% prot.nasal_relevant$from_true_name, "detected", "not_detected")) %>%
  mutate(trim = ifelse(genus %in% nasal.trim1$genus, "keep", "trimed"))



stool.data.frame.m

phylum_order <- c("Actinobacteria", "Bacteroidetes", "Firmicutes", "Proteobacteria", "Other")

stool.data.frame.m <- stool.data.frame.m %>%
  filter(!is.na(Phylum))
stool.data.frame.m$Phylum <- factor(stool.data.frame.m$Phylum, levels = phylum_order)

skin.data.frame.m <- skin.data.frame.m %>%
  filter(!is.na(Phylum))
skin.data.frame.m$Phylum <- factor(skin.data.frame.m$Phylum, levels = phylum_order)

oral.data.frame.m <- oral.data.frame.m %>%
  filter(!is.na(Phylum))
oral.data.frame.m$Phylum <- factor(oral.data.frame.m$Phylum, levels = phylum_order)

nasal.data.frame.m <- nasal.data.frame.m %>%
  filter(!is.na(Phylum))
nasal.data.frame.m$Phylum <- factor(nasal.data.frame.m$Phylum, levels = phylum_order)

#remove low confident DMI
stool.data.frame.m <- stool.data.frame.m %>%
  mutate(cytokine.engagement = if_else(trim == "trimed", "not_detected", cytokine.engagement),
         metab.engagement = if_else(trim == "trimed", "not_detected", metab.engagement),
         lipid.engagement = if_else(trim == "trimed", "not_detected", lipid.engagement),
         prot.engagement = if_else(trim == "trimed", "not_detected", prot.engagement))

skin.data.frame.m <- skin.data.frame.m %>%
  mutate(cytokine.engagement = if_else(trim == "trimed", "not_detected", cytokine.engagement),
         metab.engagement = if_else(trim == "trimed", "not_detected", metab.engagement),
         lipid.engagement = if_else(trim == "trimed", "not_detected", lipid.engagement),
         prot.engagement = if_else(trim == "trimed", "not_detected", prot.engagement))

oral.data.frame.m <- oral.data.frame.m %>%
  mutate(cytokine.engagement = if_else(trim == "trimed", "not_detected", cytokine.engagement),
         metab.engagement = if_else(trim == "trimed", "not_detected", metab.engagement),
         lipid.engagement = if_else(trim == "trimed", "not_detected", lipid.engagement),
         prot.engagement = if_else(trim == "trimed", "not_detected", prot.engagement))

nasal.data.frame.m <- nasal.data.frame.m %>%
  mutate(cytokine.engagement = if_else(trim == "trimed", "not_detected", cytokine.engagement),
         metab.engagement = if_else(trim == "trimed", "not_detected", metab.engagement),
         lipid.engagement = if_else(trim == "trimed", "not_detected", lipid.engagement),
         prot.engagement = if_else(trim == "trimed", "not_detected", prot.engagement))

#plot stool
stool_metab_plot <- stool.data.frame.m %>%  ggplot(aes(x = metab.engagement, y = dmi, color = metab.engagement)) +
  geom_boxplot() +
  geom_jitter() +
  facet_wrap(~Phylum) +
  ggtitle("Comparison of DMI based on metab engagement by Phylum") +
  xlab("Metabolites engagement") +
  ylab("DMI") + stat_compare_means(method = "wilcox.test", label = "p.format")
stool_metab_plot

stool_lipid_plot <- stool.data.frame.m %>%  ggplot(aes(x = lipid.engagement, y = dmi, color = lipid.engagement)) +
  geom_boxplot() +
  geom_jitter() +
  facet_wrap(~ Phylum) +
  ggtitle("Comparison of DMI based on lipids engagement by Phylum") +
  xlab("Lipids engagement") +
  ylab("DMI") + stat_compare_means(method = "wilcox.test", label = "p.format")
stool_lipid_plot

stool_prot_plot <- stool.data.frame.m %>%  ggplot(aes(x = prot.engagement, y = dmi, color = prot.engagement)) +
  geom_boxplot() +
  geom_jitter() +
  facet_wrap(~ Phylum) +
  ggtitle("Comparison of DMI based on proteins engagement by Phylum") +
  xlab("Proteins engagement") +
  ylab("DMI") + stat_compare_means(method = "wilcox.test", label = "p.format")
stool_prot_plot

#skin
skin_metab_plot <- skin.data.frame.m %>% ggplot(aes(x = metab.engagement, y = dmi, color = metab.engagement)) +
  geom_boxplot() +
  geom_jitter() +
  facet_wrap(~Phylum) +
  ggtitle("Comparison of DMI based on metab engagement by Phylum in Skin") +
  xlab("Metabolites engagement") +
  ylab("DMI") + stat_compare_means(method = "wilcox.test", label = "p.format")
skin_metab_plot

skin_lipid_plot <- skin.data.frame.m %>% ggplot(aes(x = lipid.engagement, y = dmi, color = lipid.engagement)) +
  geom_boxplot() +
  geom_jitter() +
  facet_wrap(~ Phylum) +
  ggtitle("Comparison of DMI based on lipids engagement by Phylum in Skin") +
  xlab("Lipids engagement") +
  ylab("DMI") + stat_compare_means(method = "wilcox.test", label = "p.format")
skin_lipid_plot

skin_prot_plot <- skin.data.frame.m %>% ggplot(aes(x = prot.engagement, y = dmi, color = prot.engagement)) +
  geom_boxplot() +
  geom_jitter() +
  facet_wrap(~ Phylum) +
  ggtitle("Comparison of DMI based on proteins engagement by Phylum in Skin") +
  xlab("Proteins engagement") +
  ylab("DMI") + stat_compare_means(method = "wilcox.test", label = "p.format")
skin_prot_plot

#oral
oral_metab_plot <- oral.data.frame.m %>% ggplot(aes(x = metab.engagement, y = dmi, color = metab.engagement)) +
  geom_boxplot() +
  geom_jitter() +
  facet_wrap(~Phylum) +
  ggtitle("Comparison of DMI based on metab engagement by Phylum") +
  xlab("Metabolites engagement") +
  ylab("DMI") + stat_compare_means(method = "wilcox.test", label = "p.format")
oral_metab_plot

oral_lipid_plot <- oral.data.frame.m %>% ggplot(aes(x = lipid.engagement, y = dmi, color = lipid.engagement)) +
  geom_boxplot() +
  geom_jitter() +
  facet_wrap(~ Phylum) +
  ggtitle("Comparison of DMI based on lipids engagement by Phylum") +
  xlab("Lipids engagement") +
  ylab("DMI") + stat_compare_means(method = "wilcox.test", label = "p.format")
oral_lipid_plot

oral_prot_plot <- oral.data.frame.m %>% ggplot(aes(x = prot.engagement, y = dmi, color = prot.engagement)) +
  geom_boxplot() +
  geom_jitter() +
  facet_wrap(~ Phylum) +
  ggtitle("Comparison of DMI based on proteins engagement by Phylum") +
  xlab("Proteins engagement") +
  ylab("DMI") + stat_compare_means(method = "wilcox.test", label = "p.format")
oral_prot_plot

nasal_metab_plot <- nasal.data.frame.m %>%  ggplot(aes(x = metab.engagement, y = dmi, color = metab.engagement)) +
  geom_boxplot() +
  geom_jitter() +
  facet_wrap(~Phylum) +
  ggtitle("Comparison of DMI based on metab engagement by Phylum") +
  xlab("Metabolites engagement") +
  ylab("DMI") + stat_compare_means(method = "wilcox.test", label = "p.format")
nasal_metab_plot

nasal_lipid_plot <- nasal.data.frame.m %>%  ggplot(aes(x = lipid.engagement, y = dmi, color = lipid.engagement)) +
  geom_boxplot() +
  geom_jitter() +
  facet_wrap(~ Phylum) +
  ggtitle("Comparison of DMI based on lipids engagement by Phylum") +
  xlab("Lipids engagement") +
  ylab("DMI") + stat_compare_means(method = "wilcox.test", label = "p.format")
nasal_lipid_plot

nasal_prot_plot <- nasal.data.frame.m %>%  ggplot(aes(x = prot.engagement, y = dmi, color = prot.engagement)) +
  geom_boxplot() +
  geom_jitter() +
  facet_wrap(~ Phylum) +
  ggtitle("Comparison of DMI based on proteins engagement by Phylum") +
  xlab("Proteins engagement") +
  ylab("DMI") + stat_compare_means(method = "wilcox.test", label = "p.format")
nasal_prot_plot

#####################################################################################################################################
#Mediation analysis
#####################################################################################################################################
Mediation <- read.csv("./human_microbiome_project/Supplementary_data/Single tables/Mediation Result.csv", header = T)
clinical.class <- read.csv("~/Library/CloudStorage/Box-Box/human_microbiome_project/Figures/Figure6/clinicaldata/51labtestclass.csv", header = T)
Mediation$phenotype_class <- clinical.class$Class[match(Mediation$phenotype_true_name, clinical.class$Abbreviation)] 

Mediation

phenotype_class_color = c(
  "Kidney" = "black" ,
  "Hematologic"= "#DF536B",
  "Liver" = "#61D04F",
  "Immnue" = "#4DBBD599",
  "Electrolyte" = "#28E2E5",
  "Kidney/Hepatic" = "#CD0BBC",
  "Cholesterol/Lipid" = "#F5C710",
  "Glucose/Insulin" = "gray62")


#significance needed here

str(Mediation)

length(unique(Mediation$treat))
length(unique(Mediation$mediator))
length(unique(Mediation$phenotype))

Mediation_ALL <- subset(Mediation, Dataset == "All")
Mediation_IS <- subset(Mediation, Dataset == "IS")
Mediation_IR <- subset(Mediation, Dataset == "IR")

Mediation_ALL <- Mediation_ALL %>% mutate(unique.medi = paste(treat_class,mediator_class,phenotype_class, sep = "_"))
Mediation_IS <- Mediation_IS %>% mutate(unique.medi = paste(treat_class,mediator_class,phenotype_class, sep = "_"))
Mediation_IR <- Mediation_IR %>% mutate(unique.medi = paste(treat_class,mediator_class,phenotype_class, sep = "_"))

table(Mediation_ALL$unique.medi)
table(Mediation_IS$unique.medi)
table(Mediation_IR$unique.medi)

df_medi_all <- data.frame(variable = rownames(table(Mediation_ALL$unique.medi)), all = as.vector(table(Mediation_ALL$unique.medi)))
df_medi_IS <- data.frame(variable = rownames(table(Mediation_IS$unique.medi)), IS = as.vector(table(Mediation_IS$unique.medi)))
df_medi_IR <- data.frame(variable = rownames(table(Mediation_IR$unique.medi)), IR = as.vector(table(Mediation_IR$unique.medi)))

df_merged.0 <- merge(df_medi_all, df_medi_IS, by = "variable", all = TRUE)
df_merged_medi <- merge(df_merged.0, df_medi_IR, by = "variable", all = TRUE)

df_merged_medi[is.na(df_merged_medi)] <- 0

#compare between IR and IS
df_merged_medi <- df_merged_medi %>%
  filter(IS != 0 | IR != 0) %>%
  rowwise() %>%
  mutate(Xsq = chisq.test(matrix(c(IS, IR, 100000000, 100000000), nrow = 2))$statistic,
         p.value = chisq.test(matrix(c(IS, IR, 100000000, 100000000), nrow = 2))$p.value,
         fisher.p = fisher.test(matrix(c(IS, IR, 100000000, 100000000), nrow = 2))$p.value)
df_merged_medi

#write.csv(file = "~/Desktop/figure6_IRIS_Fisher.csv",df_merged_medi)

#compare between four body sites
df_medi_all.internal <- df_medi_all %>%
  separate(variable, into = c("group", "path"), sep = " ", remove = FALSE) %>%
  mutate(group = gsub("Stool |Skin |Oral |Nasal ", "", group))


df_medi_all.internal_wide <- df_medi_all.internal %>%
  select(-variable) %>%
  pivot_wider(names_from = group, values_from = all)

df_medi_all.internal_wide[is.na(df_medi_all.internal_wide)] <- 0
df_medi_all.internal_wide

df_medi_all.internal_wide.p <- df_medi_all.internal_wide %>% 
  rowwise() %>%
  mutate(Nasal_Oral_Fisher_P = fisher.test(matrix(c(Nasal, Oral, 100000000, 100000000), nrow = 2))$p.value,
         Nasal_Stool_Fisher_P = fisher.test(matrix(c(Nasal, Stool, 100000000, 100000000), nrow = 2))$p.value,
         Nasal_Skin_Fisher_P = fisher.test(matrix(c(Nasal, Skin, 100000000, 100000000), nrow = 2))$p.value,
         Oral_Stool_Fisher_P = fisher.test(matrix(c(Oral, Stool, 100000000, 100000000), nrow = 2))$p.value,
         Oral_Skin_Fisher_P = fisher.test(matrix(c(Oral, Skin, 100000000, 100000000), nrow = 2))$p.value,
         Stool_Skin_Fisher_P = fisher.test(matrix(c(Stool, Skin, 100000000, 100000000), nrow = 2))$p.value)

#IR
df_IR.internal <- df_medi_IR %>%
  separate(variable, into = c("group", "path"), sep = " ", remove = FALSE) %>%
  mutate(group = gsub("Stool |Skin |Oral |Nasal ", "", group))

df_IR.internal_wide <- df_IR.internal %>%
  select(-variable) %>%
  pivot_wider(names_from = group, values_from = IR)

df_IR.internal_wide[is.na(df_IR.internal_wide)] <- 0
df_IR.internal_wide

df_IR.internal_wide.p <- df_IR.internal_wide %>% 
  rowwise() %>%
  mutate(Nasal_Oral_Fisher_P = fisher.test(matrix(c(Nasal, Oral, 100000000, 100000000), nrow = 2))$p.value,
         Nasal_Stool_Fisher_P = fisher.test(matrix(c(Nasal, Stool, 100000000, 100000000), nrow = 2))$p.value,
         Nasal_Skin_Fisher_P = fisher.test(matrix(c(Nasal, Skin, 100000000, 100000000), nrow = 2))$p.value,
         Oral_Stool_Fisher_P = fisher.test(matrix(c(Oral, Stool, 100000000, 100000000), nrow = 2))$p.value,
         Oral_Skin_Fisher_P = fisher.test(matrix(c(Oral, Skin, 100000000, 100000000), nrow = 2))$p.value,
         Stool_Skin_Fisher_P = fisher.test(matrix(c(Stool, Skin, 100000000, 100000000), nrow = 2))$p.value)


#IS
df_IS.internal <- df_medi_IS %>%
  separate(variable, into = c("group", "path"), sep = " ", remove = FALSE) %>%
  mutate(group = gsub("Stool |Skin |Oral |Nasal ", "", group))

df_IS.internal_wide <- df_IS.internal %>%
  select(-variable) %>%
  pivot_wider(names_from = group, values_from = IS)

df_IS.internal_wide[is.na(df_IS.internal_wide)] <- 0
df_IS.internal_wide

df_IS.internal_wide.p <- df_IS.internal_wide %>%
  rowwise() %>%
  mutate(Nasal_Oral_Fisher_P = fisher.test(matrix(c(Nasal, Oral, 100000000, 100000000), nrow = 2))$p.value,
         Nasal_Stool_Fisher_P = fisher.test(matrix(c(Nasal, Stool, 100000000, 100000000), nrow = 2))$p.value,
         Nasal_Skin_Fisher_P = fisher.test(matrix(c(Nasal, Skin, 100000000, 100000000), nrow = 2))$p.value,
         Oral_Stool_Fisher_P = fisher.test(matrix(c(Oral, Stool, 100000000, 100000000), nrow = 2))$p.value,
         Oral_Skin_Fisher_P = fisher.test(matrix(c(Oral, Skin, 100000000, 100000000), nrow = 2))$p.value,
         Stool_Skin_Fisher_P = fisher.test(matrix(c(Stool, Skin, 100000000, 100000000), nrow = 2))$p.value)

colnames(df_medi_all.internal_wide.p)
colnames(df_IR.internal_wide.p)
colnames(df_IS.internal_wide.p)

df_medi_all.internal_wide.p$dataset <- "All"
df_IR.internal_wide.p$dataset <- "IR"
df_IS.internal_wide.p$dataset <- "IS"

internal.compare <- rbind(df_medi_all.internal_wide.p,df_IR.internal_wide.p) %>% rbind(df_IS.internal_wide.p)
internal.compare
#write.csv(file = "~/Desktop/1st.revision.CELL-S-22-04878/Supplementary_data/Supplimentaryfigure_added/Figure6.Comparewithindataset.csv",internal.compare)


IRIS.stool.compare <- rbind(stool.data.frame %>% filter(genus %in% (filter(Mediation_IS, treat_class== "Stool microbiome")$treat_true_name %>% unique())) %>% mutate(IRIS = "IS"),
                            stool.data.frame %>% filter(genus %in% (filter(Mediation_IR, treat_class== "Stool microbiome")$treat_true_name %>% unique())) %>% mutate(IRIS = "IR"))

IRIS.skin.compare <- rbind(skin.data.frame %>% filter(genus %in% (filter(Mediation_IS, treat_class== "Skin microbiome")$treat_true_name %>% unique())) %>% mutate(IRIS = "IS"),
                           skin.data.frame %>% filter(genus %in% (filter(Mediation_IR, treat_class== "Skin microbiome")$treat_true_name %>% unique())) %>% mutate(IRIS = "IR"))

IRIS.oral.compare <- rbind(oral.data.frame %>% filter(genus %in% (filter(Mediation_IS, treat_class== "Oral microbiome")$treat_true_name %>% unique())) %>% mutate(IRIS = "IS"),
                           oral.data.frame %>% filter(genus %in% (filter(Mediation_IR, treat_class== "Oral microbiome")$treat_true_name %>% unique())) %>% mutate(IRIS = "IR"))

IRIS.nasal.compare <- rbind(nasal.data.frame %>% filter(genus %in% (filter(Mediation_IS, treat_class== "Nasal microbiome")$treat_true_name %>% unique())) %>% mutate(IRIS = "IS"),
                            nasal.data.frame %>% filter(genus %in% (filter(Mediation_IR, treat_class== "Nasal microbiome")$treat_true_name %>% unique())) %>% mutate(IRIS = "IR"))

IRIS.stool.compare
iris_stool_plot <- ggplot(IRIS.stool.compare, aes(x = IRIS, y = dmi)) +
  geom_jitter() +
  geom_boxplot(fill = "#69b3a2", alpha = 0.4) +
  ggtitle("Comparison of DMI based on IRIS engagement") +
  xlab("IRIS engagement") +
  ylab("DMI") +
  stat_compare_means(method = "t.test", paired = FALSE, label = "p.format")

iris_stool_plot

iris_skin_plot <- ggplot(IRIS.skin.compare, aes(x = IRIS, y = dmi)) +
  geom_jitter() +
  geom_boxplot(fill = "#69b3a2", alpha = 0.4) +
  ggtitle("Comparison of DMI based on IRIS engagement") +
  xlab("IRIS engagement") +
  ylab("DMI") +
  stat_compare_means(method = "t.test", paired = FALSE, label = "p.format")

iris_skin_plot

iris_oral_plot <- ggplot(IRIS.oral.compare, aes(x = IRIS, y = dmi)) +
  geom_jitter() +
  geom_boxplot(fill = "#69b3a2", alpha = 0.4) +
  ggtitle("Comparison of DMI based on IRIS engagement") +
  xlab("IRIS engagement") +
  ylab("DMI") +
  stat_compare_means(method = "t.test", paired = FALSE, label = "p.format")

iris_oral_plot

iris_nasal_plot <- ggplot(IRIS.nasal.compare, aes(x = IRIS, y = dmi)) +
  geom_jitter() +
  geom_boxplot(fill = "#69b3a2", alpha = 0.4) +
  ggtitle("Comparison of DMI based on IRIS engagement") +
  xlab("IRIS engagement") +
  ylab("DMI") +
  stat_compare_means(method = "t.test", paired = FALSE, label = "p.format")

iris_nasal_plot


#add figure 3D data 
microbe_corre <- read.csv("~/Desktop/microbiome_correlation.csv", header = T)

table(microbe_corre$compare)

table(microbe_corre$compare, microbe_corre$from_class)
table(microbe_corre$compare, microbe_corre$to_class)

within.list.stool <- microbe_corre %>% filter(from_class == "stool") %>% filter(compare == "within")
within.list.skin <- microbe_corre %>% filter(from_class == "skin") %>% filter(compare == "within")
within.list.oral <- microbe_corre %>% filter(from_class == "oral") %>% filter(compare == "within")
within.list.nasal <- microbe_corre %>% filter(from_class == "nasal") %>% filter(compare == "within")

between.list.stool1 <- microbe_corre %>% filter(compare == "between") %>% filter(from_class == "stool") %>% select(from_Genus) %>% rename(Genus = from_Genus)
between.list.stool2 <- microbe_corre %>% filter(compare == "between") %>% filter(to_class == "stool") %>% select(to_Genus) %>% rename(Genus = to_Genus)
between.list.stool <- unique(rbind(between.list.stool1,between.list.stool2))

between.list.skin1 <- microbe_corre %>% filter(compare == "between") %>% filter(from_class == "skin") %>% select(from_Genus) %>% rename(Genus = from_Genus)
between.list.skin2 <- microbe_corre %>% filter(compare == "between") %>% filter(to_class == "skin") %>% select(to_Genus) %>% rename(Genus = to_Genus)
between.list.skin <- unique(rbind(between.list.skin1,between.list.skin2))

between.list.oral1 <- microbe_corre %>% filter(compare == "between") %>% filter(from_class == "oral") %>% select(from_Genus) %>% rename(Genus = from_Genus)
between.list.oral2 <- microbe_corre %>% filter(compare == "between") %>% filter(to_class == "oral") %>% select(to_Genus) %>% rename(Genus = to_Genus)
between.list.oral <- unique(rbind(between.list.oral1,between.list.oral2))

between.list.nasal1 <- microbe_corre %>% filter(compare == "between") %>% filter(from_class == "nasal") %>% select(from_Genus) %>% rename(Genus = from_Genus)
between.list.nasal2 <- microbe_corre %>% filter(compare == "between") %>% filter(to_class == "nasal") %>% select(to_Genus) %>% rename(Genus = to_Genus)
between.list.nasal <- unique(rbind(between.list.nasal1,between.list.nasal2))

#the within body site engagement
stool.data.frame.m2 <- stool.data.frame.m %>%
  mutate(within.engagement = ifelse(genus %in% within.list.stool$from_Genus, "detected", "not_detected"))

skin.data.frame.m2 <- skin.data.frame.m %>%
  mutate(within.engagement = ifelse(genus %in% within.list.skin$from_Genus, "detected", "not_detected"))

oral.data.frame.m2 <- oral.data.frame.m %>%
  mutate(within.engagement = ifelse(genus %in% within.list.oral$from_Genus, "detected", "not_detected"))

nasal.data.frame.m2 <- nasal.data.frame.m %>%
  mutate(within.engagement = ifelse(genus %in% within.list.nasal$from_Genus, "detected", "not_detected"))



#the between body site engagement
stool.data.frame.m3 <- stool.data.frame.m2 %>%
  mutate(between.engagement = ifelse(genus %in% between.list.stool$Genus, "detected", "not_detected"))

skin.data.frame.m3 <- skin.data.frame.m2 %>%
  mutate(between.engagement = ifelse(genus %in% between.list.skin$Genus, "detected", "not_detected"))

oral.data.frame.m3 <- oral.data.frame.m2 %>%
  mutate(between.engagement = ifelse(genus %in% between.list.oral$Genus, "detected", "not_detected"))

nasal.data.frame.m3 <- nasal.data.frame.m2 %>%
  mutate(between.engagement = ifelse(genus %in% between.list.nasal$Genus, "detected", "not_detected"))

#the between body site engagement
stool.data.frame.m4 <- stool.data.frame.m3 %>%
  mutate(medi.engagement_IS = ifelse(genus %in% filter(Mediation_IS, treat_class == "Stool microbiome")$Genus, "detected", "not_detected")) %>% 
  mutate(medi.engagement_IR = ifelse(genus %in% filter(Mediation_IR, treat_class == "Stool microbiome")$Genus, "detected", "not_detected")) %>% 
  mutate(medi.engagement_ALL = ifelse(genus %in% filter(Mediation_ALL, treat_class == "Stool microbiome")$Genus, "detected", "not_detected"))

skin.data.frame.m4 <- skin.data.frame.m3 %>%
  mutate(medi.engagement_IS = ifelse(genus %in% filter(Mediation_IS, treat_class == "Skin microbiome")$Genus, "detected", "not_detected")) %>% 
  mutate(medi.engagement_IR = ifelse(genus %in% filter(Mediation_IR, treat_class == "Skin microbiome")$Genus, "detected", "not_detected")) %>% 
  mutate(medi.engagement_ALL = ifelse(genus %in% filter(Mediation_ALL, treat_class == "Skin microbiome")$Genus, "detected", "not_detected"))

oral.data.frame.m4 <- oral.data.frame.m3 %>%
  mutate(medi.engagement_IS = ifelse(genus %in% filter(Mediation_IS, treat_class == "Oral microbiome")$Genus, "detected", "not_detected")) %>% 
  mutate(medi.engagement_IR = ifelse(genus %in% filter(Mediation_IR, treat_class == "Oral microbiome")$Genus, "detected", "not_detected")) %>% 
  mutate(medi.engagement_ALL = ifelse(genus %in% filter(Mediation_ALL, treat_class == "Oral microbiome")$Genus, "detected", "not_detected"))

nasal.data.frame.m4 <- nasal.data.frame.m3 %>%
  mutate(medi.engagement_IS = ifelse(genus %in% filter(Mediation_IS, treat_class == "Nasal microbiome")$Genus, "detected", "not_detected")) %>% 
  mutate(medi.engagement_IR = ifelse(genus %in% filter(Mediation_IR, treat_class == "Nasal microbiome")$Genus, "detected", "not_detected")) %>%
  mutate(medi.engagement_ALL = ifelse(genus %in% filter(Mediation_ALL, treat_class == "Nasal microbiome")$Genus, "detected", "not_detected"))

#within
p.within.stool.dmi <- stool.data.frame.m3 %>% filter(trim != "trimed") %>% ggplot(aes(x=within.engagement,y=dmi))  + geom_boxplot(fill=body_site_color[1],alpha=0.3, outlier.alpha = 0) + geom_jitter(size=0.7)
p.within.stool.dmi <- p.within.stool.dmi+ stat_compare_means(label="p.signif", color="red") + base_theme + theme(axis.title.x = element_blank()) + ggtitle("The Stool Microbiome")
p.within.stool.dmi

wilcox.test(dmi ~ within.engagement, data = stool.data.frame.m3 %>% filter(trim != "trimed"))

ggsave(filename = "~/Desktop/DMIengagement/p.within.stool.dmi.pdf",p.within.stool.dmi, width = 2, height = 2.5, dpi = 300)

p.within.skin.dmi <- skin.data.frame.m3 %>% filter(trim != "trimed") %>% ggplot(aes(x=within.engagement,y=dmi)) + geom_boxplot(fill=body_site_color[2],alpha=0.3,outlier.alpha = 0) + geom_jitter(size=0.7)
p.within.skin.dmi <- p.within.skin.dmi + stat_compare_means(label="p.signif",color="red") + base_theme  + theme(axis.title.x = element_blank()) + ggtitle("The Skin Microbiome")
p.within.skin.dmi

wilcox.test(dmi ~ within.engagement, data = skin.data.frame.m3 %>% filter(trim != "trimed"))

ggsave(filename = "~/Desktop/DMIengagement/p.within.skin.dmi.pdf",p.within.skin.dmi, width = 2, height = 2.5, dpi = 300)

p.within.oral.dmi <- oral.data.frame.m3 %>% filter(trim != "trimed") %>% ggplot(aes(x=within.engagement,y=dmi)) + geom_boxplot(fill=body_site_color[3],alpha=0.3,outlier.alpha = 0) + geom_jitter(size=0.7)
p.within.oral.dmi <- p.within.oral.dmi + stat_compare_means(label="p.signif",color="red") + base_theme  + theme(axis.title.x = element_blank()) + ggtitle("The Oral Microbiome")
p.within.oral.dmi

wilcox.test(dmi ~ within.engagement, data = oral.data.frame.m3 %>% filter(trim != "trimed"))

ggsave(filename = "~/Desktop/DMIengagement/p.within.oral.dmi.pdf",p.within.oral.dmi, width = 2, height = 2.5, dpi = 300)

p.within.nasal.dmi <- nasal.data.frame.m3 %>% filter(trim != "trimed") %>% ggplot(aes(x=within.engagement,y=dmi)) + geom_boxplot(fill=body_site_color[4],alpha=0.3,outlier.alpha = 0) + geom_jitter(size=0.7)
p.within.nasal.dmi <- p.within.nasal.dmi + stat_compare_means(label="p.signif", color="red") + base_theme  + theme(axis.title.x = element_blank()) + ggtitle("The Nasal Microbiome")
p.within.nasal.dmi

wilcox.test(dmi ~ within.engagement, data = nasal.data.frame.m3 %>% filter(trim != "trimed"))

ggsave(filename = "~/Desktop/DMIengagement/p.within.nasal.dmi.pdf",p.within.nasal.dmi, width = 2, height = 2.5, dpi = 300)


#between
p.between.stool.dmi <- stool.data.frame.m3 %>% filter(trim != "trimed") %>% ggplot(aes(x=between.engagement,y=dmi))  + geom_boxplot(fill=body_site_color[1],alpha=0.3,outlier.alpha = 0) + geom_jitter(size=0.5)
p.between.stool.dmi <- p.between.stool.dmi+ stat_compare_means(label="p.signif", color="red") + base_theme  + theme(axis.title.x = element_blank()) + ggtitle("The Stool Microbiome")
p.between.stool.dmi

wilcox.test(dmi ~ between.engagement, data = stool.data.frame.m3 %>% filter(trim != "trimed"))

ggsave(filename = "~/Desktop/DMIengagement/p.between.stool.dmi.pdf",p.between.stool.dmi, width = 3, height =4, dpi = 300)

p.between.skin.dmi <- skin.data.frame.m3 %>% filter(trim != "trimed") %>% ggplot(aes(x=between.engagement,y=dmi)) + geom_boxplot(fill=body_site_color[2],alpha=0.3,outlier.alpha = 0) + geom_jitter(size=0.5)
p.between.skin.dmi <- p.between.skin.dmi + stat_compare_means(label="p.signif", color="red") + base_theme  + theme(axis.title.x = element_blank()) + ggtitle("The Skin Microbiome")
p.between.skin.dmi

wilcox.test(dmi ~ between.engagement, data = skin.data.frame.m3 %>% filter(trim != "trimed"))

ggsave(filename = "~/Desktop/DMIengagement/p.between.skin.dmi.pdf",p.between.skin.dmi, width = 3, height =4, dpi = 300)

p.between.oral.dmi <- oral.data.frame.m3 %>% filter(trim != "trimed") %>% ggplot(aes(x=between.engagement,y=dmi)) + geom_boxplot(fill=body_site_color[3],alpha=0.3,outlier.alpha = 0) + geom_jitter(size=0.5)
p.between.oral.dmi <- p.between.oral.dmi + stat_compare_means(label="p.signif", color="red") + base_theme  + theme(axis.title.x = element_blank()) + ggtitle("The Oral Microbiome")
p.between.oral.dmi

wilcox.test(dmi ~ between.engagement, data = oral.data.frame.m3 %>% filter(trim != "trimed"))


ggsave(filename = "~/Desktop/DMIengagement/p.between.oral.dmi.pdf",p.between.oral.dmi, width = 3, height =4, dpi = 300)

p.between.nasal.dmi <- nasal.data.frame.m3 %>% filter(trim != "trimed") %>% ggplot(aes(x=between.engagement,y=dmi)) + geom_boxplot(fill=body_site_color[4],alpha=0.3, outlier.alpha = 0) + geom_jitter(size=0.5)
p.between.nasal.dmi <- p.between.nasal.dmi + stat_compare_means(label="p.signif", color="red") + base_theme  + theme(axis.title.x = element_blank()) + ggtitle("The Nasal Microbiome")
p.between.nasal.dmi

wilcox.test(dmi ~ between.engagement, data = nasal.data.frame.m3 %>% filter(trim != "trimed"))

ggsave(filename = "~/Desktop/DMIengagement/p.between.nasal.dmi.pdf", p.between.nasal.dmi, width = 3, height =4, dpi = 300)

#for cytokines 
p.cytokine.stool.dmi <- stool.data.frame.m4 %>% filter(trim != "trimed") %>% ggplot(aes(x=cytokine.engagement ,y=dmi))  + geom_boxplot(outlier.alpha = 0) + geom_jitter(size=0.5)
p.cytokine.stool.dmi <- p.cytokine.stool.dmi+ stat_compare_means(color="red") + theme_cowplot() + theme(axis.title.x = element_blank()) + ggtitle("The Stool Microbiome")
#p.cytokine.stool.dmi <- p.cytokine.stool.dmi + facet_grid(.~Phylum)
p.cytokine.stool.dmi

p.cytokine.skin.dmi <- skin.data.frame.m3 %>% filter(trim != "trimed") %>% ggplot(aes(x=cytokine.engagement ,y=dmi)) + 
  geom_boxplot(outlier.alpha = 0) + 
  geom_jitter(size=0.5) + 
  stat_compare_means(color="red") + 
  theme_cowplot() + 
  theme(axis.title.x = element_blank()) + 
  ggtitle("The Skin Microbiome") +
  facet_grid(.~Phylum)
p.cytokine.skin.dmi

p.cytokine.oral.dmi <- oral.data.frame.m3 %>% 
  filter(trim != "trimmed") %>% 
  ggplot(aes(x=cytokine.engagement ,y=dmi)) + 
  geom_boxplot(outlier.alpha = 0) + 
  geom_jitter(size=0.5) + 
  stat_compare_means(color="red") + 
  theme_cowplot() + 
  theme(axis.title.x = element_blank()) + 
  ggtitle("The Oral Microbiome") +
  facet_grid(.~Phylum)
p.cytokine.oral.dmi

p.cytokine.nasal.dmi <- nasal.data.frame.m3 %>% 
  filter(trim != "trimmed") %>% 
  ggplot(aes(x=cytokine.engagement ,y=dmi)) + 
  geom_boxplot(outlier.alpha = 0) + 
  geom_jitter(size=0.5) + 
  stat_compare_means(color="red") + 
  theme_cowplot() + 
  theme(axis.title.x = element_blank()) + 
  ggtitle("The Nasal Microbiome") +
  facet_grid(.~Phylum)
p.cytokine.nasal.dmi

# #for protein not finished yet
# p.protein.stool.dmi <- stool.data.frame.m3 %>% filter(trim != "trimmed") %>% ggplot(aes(x=prot.engagement ,y=dmi))  + geom_boxplot(outlier.alpha = 0) + geom_jitter(size=0.5)
# p.protein.stool.dmi <- p.protein.stool.dmi+ stat_compare_means(color="red") + theme_cowplot() + theme(axis.title.x = element_blank()) + ggtitle("The Stool Microbiome")
# p.protein.stool.dmi
# 
# p.within.skin.dmi <- skin.data.frame.m3 %>% filter(trim != "trimmed") %>% ggplot(aes(x=prot.engagement ,y=dmi)) + geom_boxplot(outlier.alpha = 0) + geom_jitter(size=0.5)
# p.within.skin.dmi <- p.within.skin.dmi + stat_compare_means(color="red") + theme_cowplot() + theme(axis.title.x = element_blank()) + ggtitle("The Skin Microbiome")
# p.within.skin.dmi
# 
# p.within.oral.dmi <- oral.data.frame.m3 %>% filter(trim != "trimmed") %>% ggplot(aes(x=prot.engagement ,y=dmi)) + geom_boxplot(outlier.alpha = 0) + geom_jitter(size=0.5)
# p.within.oral.dmi <- p.within.oral.dmi + stat_compare_means(color="red") + theme_cowplot() + theme(axis.title.x = element_blank()) + ggtitle("The Oral Microbiome")
# p.within.oral.dmi
# 
# p.within.nasal.dmi <- nasal.data.frame.m3 %>% filter(trim != "trimmed") %>% ggplot(aes(x=prot.engagement ,y=dmi)) + geom_boxplot(outlier.alpha = 0) + geom_jitter(size=0.5)
# p.within.nasal.dmi <- p.within.nasal.dmi + stat_compare_means(color="red") + theme_cowplot() + theme(axis.title.x = element_blank()) + ggtitle("The Nasal Microbiome")
# p.within.nasal.dmi
# 
# 
# #for lipid not finished yet
# p.within.stool.dmi <- stool.data.frame.m3 %>% filter(trim != "trimmed") %>% ggplot(aes(x=lipid.engagement ,y=dmi))  + geom_boxplot(outlier.alpha = 0) + geom_jitter(size=0.5)
# p.within.stool.dmi <- p.within.stool.dmi+ stat_compare_means(color="red") + theme_cowplot() + theme(axis.title.x = element_blank()) + ggtitle("The Stool Microbiome")
# p.within.stool.dmi
# 
# p.within.skin.dmi <- skin.data.frame.m3 %>% filter(trim != "trimmed") %>% ggplot(aes(x=lipid.engagement ,y=dmi)) + geom_boxplot(outlier.alpha = 0) + geom_jitter(size=0.5)
# p.within.skin.dmi <- p.within.skin.dmi + stat_compare_means(color="red") + theme_cowplot() + theme(axis.title.x = element_blank()) + ggtitle("The Skin Microbiome")
# p.within.skin.dmi
# 
# p.within.oral.dmi <- oral.data.frame.m3 %>% filter(trim != "trimmed") %>% ggplot(aes(x=lipid.engagement ,y=dmi)) + geom_boxplot(outlier.alpha = 0) + geom_jitter(size=0.5)
# p.within.oral.dmi <- p.within.oral.dmi + stat_compare_means(color="red") + theme_cowplot() + theme(axis.title.x = element_blank()) + ggtitle("The Oral Microbiome")
# p.within.oral.dmi
# 
# p.within.nasal.dmi <- nasal.data.frame.m3 %>% filter(trim != "trimmed") %>% ggplot(aes(x=lipid.engagement ,y=dmi)) + geom_boxplot(outlier.alpha = 0) + geom_jitter(size=0.5)
# p.within.nasal.dmi <- p.within.nasal.dmi + stat_compare_means(color="red") + theme_cowplot() + theme(axis.title.x = element_blank()) + ggtitle("The Nasal Microbiome")
# p.within.nasal.dmi

stool.data.frame.m3 %>% filter(trim != "trimed")

colnames(stool.data.frame.m3)

# Select relevant columns and reshape the data
#Perfrom Analysis for Stool
stool.data_long <- stool.data.frame.m4 %>% 
  filter(trim != "trimed") %>% 
  select(dmi, cytokine.engagement, metab.engagement, lipid.engagement, prot.engagement, within.engagement, between.engagement,medi.engagement_ALL) %>% 
  pivot_longer(-dmi, names_to = "engagement_type", values_to = "detection") %>% 
  mutate(engagement_type = factor(engagement_type, levels = c('within.engagement', 'between.engagement', 'cytokine.engagement', 'metab.engagement', 'lipid.engagement', 'prot.engagement', 'medi.engagement_ALL')))

stool.engagement_plots <- lapply(levels(stool.data_long$engagement_type), function(engagement) {
  p <- ggplot(filter(stool.data_long, engagement_type == engagement), aes(x = detection, y = dmi)) + 
    geom_boxplot() + 
    stat_compare_means(aes(group=detection), label="p.format", color= "red") +
    theme_minimal() + 
    labs(x = engagement, y = "DMI")
  return(p)
})

stool.engagement_plots[[7]]

stool_combined_plot <- plot_grid(plotlist = stool.engagement_plots, ncol = 7)
stool_combined_plot
#ggsave(filename = "~/Desktop/DMIengagement/stool_combined_plot.pdf",stool_combined_plot, width = 14, height = 5, dpi = 300)

#compare mediation analysis
stool.data_long_Medi <- stool.data.frame.m4 %>% 
  filter(trim != "trimed") %>% 
  select(dmi, medi.engagement_IS,medi.engagement_IR) %>% 
  pivot_longer(-dmi, names_to = "engagement_type", values_to = "detection") %>% 
  mutate(engagement_type = factor(engagement_type, levels = c('medi.engagement_IS', 'medi.engagement_IR')))

stool.engagement.medi_plots <- lapply(levels(stool.data_long_Medi$engagement_type), function(engagement) {
  p <- ggplot(filter(stool.data_long_Medi, engagement_type == engagement), aes(x = detection, y = dmi)) + 
    geom_boxplot() + 
    stat_compare_means(aes(group=detection), label="p.format", color= "red") +
    theme_minimal() + 
    labs(x = engagement, y = "DMI")
  return(p)
})

stool_combined_medi_plot <- plot_grid(plotlist = stool.engagement.medi_plots, ncol =2)
stool_combined_medi_plot
#ggsave(filename = "~/Desktop/DMIengagement/stool_combined_medi_plot.pdf",stool_combined_medi_plot, width = 5, height = 5, dpi = 300)

#Perfrom Analysis for Skin
skin.data_long <- skin.data.frame.m4 %>% 
  filter(trim != "trimed") %>% 
  select(dmi, cytokine.engagement, metab.engagement, lipid.engagement, prot.engagement, within.engagement, between.engagement,medi.engagement_ALL) %>% 
  pivot_longer(-dmi, names_to = "engagement_type", values_to = "detection") %>% 
  mutate(engagement_type = factor(engagement_type, levels = c('within.engagement', 'between.engagement', 'cytokine.engagement', 'metab.engagement', 'lipid.engagement', 'prot.engagement', 'medi.engagement_ALL')))

skin.engagement_plots <- lapply(levels(skin.data_long$engagement_type), function(engagement) {
  p <- ggplot(filter(skin.data_long, engagement_type == engagement), aes(x = detection, y = dmi)) + 
    geom_boxplot() + 
    stat_compare_means(aes(group=detection), label="p.format", color= "red") +
    theme_minimal() + 
    labs(x = engagement, y = "DMI")
  return(p)
})

skin.engagement_plots[[7]]

skin_combined_plot <- plot_grid(plotlist = skin.engagement_plots, ncol = 7)
skin_combined_plot
#ggsave(filename = "~/Desktop/DMIengagement/skin_combined_plot.pdf",skin_combined_plot, width = 14, height = 5, dpi = 300)

#compare mediation analysis
skin.data_long_Medi <- skin.data.frame.m4 %>% 
  filter(trim != "trimed") %>% 
  select(dmi, medi.engagement_IS,medi.engagement_IR) %>% 
  pivot_longer(-dmi, names_to = "engagement_type", values_to = "detection") %>% 
  mutate(engagement_type = factor(engagement_type, levels = c('medi.engagement_IS', 'medi.engagement_IR')))

skin.engagement.medi_plots <- lapply(levels(skin.data_long_Medi$engagement_type), function(engagement) {
  p <- ggplot(filter(skin.data_long_Medi, engagement_type == engagement), aes(x = detection, y = dmi)) + 
    geom_boxplot() + 
    stat_compare_means(aes(group=detection), label="p.format", color= "red") +
    theme_minimal() + 
    labs(x = engagement, y = "DMI")
  return(p)
})

skin_combined_medi_plot <- plot_grid(plotlist = skin.engagement.medi_plots, ncol =2)
skin_combined_medi_plot
#ggsave(filename = "~/Desktop/DMIengagement/skin_combined_medi_plot.pdf",skin_combined_medi_plot, width = 5, height = 5, dpi = 300)

# Perform Analysis for Oral
oral.data_long <- oral.data.frame.m4 %>% 
  filter(trim != "trimed") %>% 
  select(dmi, cytokine.engagement, metab.engagement, lipid.engagement, prot.engagement, within.engagement, between.engagement, medi.engagement_ALL) %>% 
  pivot_longer(-dmi, names_to = "engagement_type", values_to = "detection") %>% 
  mutate(engagement_type = factor(engagement_type, levels = c('within.engagement', 'between.engagement', 'cytokine.engagement', 'metab.engagement', 'lipid.engagement', 'prot.engagement', 'medi.engagement_ALL')))

oral.engagement_plots <- lapply(levels(oral.data_long$engagement_type), function(engagement) {
  p <- ggplot(filter(oral.data_long, engagement_type == engagement), aes(x = detection, y = dmi)) + 
    geom_boxplot() + 
    stat_compare_means(aes(group=detection), label="p.format", color= "red") +
    theme_minimal() + 
    labs(x = engagement, y = "DMI")
  return(p)
})

oral.engagement_plots[[7]]

oral_combined_plot <- plot_grid(plotlist = oral.engagement_plots, ncol = 7)
oral_combined_plot
#ggsave(filename = "~/Desktop/DMIengagement/oral_combined_plot.pdf", oral_combined_plot, width = 14, height = 5, dpi = 300)

# Compare mediation analysis
oral.data_long_Medi <- oral.data.frame.m4 %>% 
  filter(trim != "trimed") %>% 
  select(dmi, medi.engagement_IS, medi.engagement_IR) %>% 
  pivot_longer(-dmi, names_to = "engagement_type", values_to = "detection") %>% 
  mutate(engagement_type = factor(engagement_type, levels = c('medi.engagement_IS', 'medi.engagement_IR')))

oral.engagement.medi_plots <- lapply(levels(oral.data_long_Medi$engagement_type), function(engagement) {
  p <- ggplot(filter(oral.data_long_Medi, engagement_type == engagement), aes(x = detection, y = dmi)) + 
    geom_boxplot() + 
    stat_compare_means(aes(group=detection), label="p.format", color= "red") +
    theme_minimal() + 
    labs(x = engagement, y = "DMI")
  return(p)
})

oral_combined_medi_plot <- plot_grid(plotlist = oral.engagement.medi_plots, ncol = 2)
oral_combined_medi_plot
#ggsave(filename = "~/Desktop/DMIengagement/oral_combined_medi_plot.pdf", oral_combined_medi_plot, width = 5, height = 5, dpi = 300)

# Perform Analysis for Nasal
nasal.data_long <- nasal.data.frame.m4 %>% 
  filter(trim != "trimed") %>% 
  select(dmi, cytokine.engagement, metab.engagement, lipid.engagement, prot.engagement, within.engagement, between.engagement, medi.engagement_ALL) %>% 
  pivot_longer(-dmi, names_to = "engagement_type", values_to = "detection") %>% 
  mutate(engagement_type = factor(engagement_type, levels = c('within.engagement', 'between.engagement', 'cytokine.engagement', 'metab.engagement', 'lipid.engagement', 'prot.engagement', 'medi.engagement_ALL')))

nasal.engagement_plots <- lapply(levels(nasal.data_long$engagement_type), function(engagement) {
  p <- ggplot(filter(nasal.data_long, engagement_type == engagement), aes(x = detection, y = dmi)) + 
    geom_boxplot() + 
    stat_compare_means(aes(group=detection), label="p.format", color= "red") +
    theme_minimal() + 
    labs(x = engagement, y = "DMI")
  return(p)
})

nasal.engagement_plots[[7]]

nasal_combined_plot <- plot_grid(plotlist = nasal.engagement_plots, ncol = 7)
nasal_combined_plot
#ggsave(filename = "~/Desktop/DMIengagement/nasal_combined_plot.pdf", nasal_combined_plot, width = 14, height = 5, dpi = 300)

# Compare mediation analysis
nasal.data_long_Medi <- nasal.data.frame.m4 %>% 
  filter(trim != "trimed") %>% 
  select(dmi, medi.engagement_IS, medi.engagement_IR) %>% 
  pivot_longer(-dmi, names_to = "engagement_type", values_to = "detection") %>% 
  mutate(engagement_type = factor(engagement_type, levels = c('medi.engagement_IS', 'medi.engagement_IR')))

nasal.engagement.medi_plots <- lapply(levels(nasal.data_long_Medi$engagement_type), function(engagement) {
  p <- ggplot(filter(nasal.data_long_Medi, engagement_type == engagement), aes(x = detection, y = dmi)) + 
    geom_boxplot() + 
    stat_compare_means(aes(group=detection), label="p.format", color= "red") +
    theme_minimal() + 
    labs(x = engagement, y = "DMI")
  return(p)
})

nasal_combined_medi_plot <- plot_grid(plotlist = nasal.engagement.medi_plots, ncol = 2)
nasal_combined_medi_plot
#ggsave(filename = "~/Desktop/DMIengagement/nasal_combined_medi_plot.pdf", nasal_combined_medi_plot, width = 5, height = 5, dpi = 300)


