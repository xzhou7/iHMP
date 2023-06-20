#Season on single bacteria
library(readxl)
library(tidyverse)
library(ggpubr)
library(dplyr)
library(cowplot)

setwd("~/Library/CloudStorage/Box-Box/human_microbiome_project/")
source("./code/tools.R")

stool.asv.table <- read_xlsx("./data_analysis/stool_microbiome/season_analysis/asv_relative/season_p_value.xlsx")
skin.asv.table <- read_xlsx("./data_analysis/skin_microbiome/season_analysis/asv_relative/season_p_value.xlsx")
oral.asv.table <- read_xlsx("./data_analysis/oral_microbiome/season_analysis/asv_relative/season_p_value.xlsx")
nasal.asv.table <- read_xlsx("./data_analysis/nasal_microbiome/season_analysis/asv_relative/season_p_value.xlsx")

# stool.genus.table <- read_xlsx("./data_analysis/stool_microbiome/season_analysis/genus_relative/season_p_value.xlsx")
# skin.genus.table <- read_xlsx("./data_analysis/skin_microbiome/season_analysis/genus_relative/season_p_value.xlsx")
# oral.genus.table <- read_xlsx("./data_analysis/oral_microbiome/season_analysis/genus_relative/season_p_value.xlsx")
# nasal.genus.table <- read_xlsx("./data_analysis/nasal_microbiome/season_analysis/genus_relative/season_p_value.xlsx")

#1. using ASV
stool.table <- stool.asv.table
skin.table <- skin.asv.table
oral.table <- oral.asv.table
nasal.table <- nasal.asv.table

#2. using genus
# stool.table <- stool.genus.table
# skin.table <- skin.genus.table
# oral.table <- oral.genus.table
# nasal.table <- nasal.genus.table


stool.table$bodysite <- "Stool"
skin.table$bodysite <- "Skin"
oral.table$bodysite <- "Oral"
nasal.table$bodysite <- "Nasal"

season.p.value <- rbind(stool.table,skin.table,oral.table,nasal.table)
season.p.value$bodysite <- factor(season.p.value$bodysite, levels = c("Stool", "Skin", "Oral", "Nasal"))
table(season.p.value$fdr < 0.05)
table(season.p.value$fdr < 0.2)

filter(season.p.value, fdr < 0.05)

table(filter(season.p.value, p_value < 0.05)$bodysite)

table(filter(season.p.value, fdr < 0.05)$bodysite)

bdst.comparasion <- list(c("Stool", "Skin"), c("Stool", "Oral"), c("Stool", "Nasal"), c("Skin", "Oral"), c("Skin", "Nasal"), c("Oral", "Nasal"))

p1 <- filter(season.p.value, fdr < 0.2) %>% ggplot(aes(x=bodysite, y=-log(abs(winter_summer)), fill=bodysite )) + 
  geom_jitter() + geom_boxplot(alpha=0.4, outlier.alpha = 0) + stat_compare_means(comparisons = bdst.comparasion) + theme_cowplot() + scale_fill_manual(values = body_site_color)
p1



p2 <- filter(season.p.value, fdr < 0.2) %>% ggplot(aes(x=bodysite, y=-log(winter_mean + summer_mean))) + 
  geom_jitter() + geom_boxplot(alpha=0.4, outlier.alpha = 0)
p2

#write.csv(file = "~/Desktop/Season.ASV.pvalue.csv",season.p.value)
#write.csv(file = "~/Desktop/Season.genus.pvalue.csv",season.p.value)



filtered_data <- filter(season.p.value, fdr < 0.2)

mean_by_bodysite <- filtered_data %>%
  group_by(bodysite) %>%
  summarise(mean_winter_summer = mean(winter_summer, na.rm = TRUE))

print(mean_by_bodysite)







      