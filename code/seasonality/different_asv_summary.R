no_source()

# set work directory
setwd(masstools::get_project_wd())
library(tidyverse)
library(phyloseq)
rm(list = ls())

source("code/tools.R")

{
  ##load data
  load("data_analysis/nasal_microbiome/season_analysis/season_p_value")
  nasal_p_value <- season_p_value
  
  load("data_analysis/skin_microbiome/season_analysis/season_p_value")
  skin_p_value <- season_p_value
  
  load("data_analysis/stool_microbiome/season_analysis/season_p_value")
  stool_p_value <- season_p_value
  
  load("data_analysis/oral_microbiome/season_analysis/season_p_value")
  oral_p_value <- season_p_value
}

sum(skin_p_value$p_value < 0.05)
sum(skin_p_value$fdr < 0.05)
sum(skin_p_value$p_value < 0.05) / nrow(skin_p_value)
sum(skin_p_value$fdr < 0.05) / nrow(skin_p_value)

sum(nasal_p_value$p_value < 0.05)
sum(nasal_p_value$fdr < 0.05)
sum(nasal_p_value$p_value < 0.05) / nrow(nasal_p_value)
sum(nasal_p_value$fdr < 0.05) / nrow(nasal_p_value)

sum(oral_p_value$p_value < 0.05)
sum(oral_p_value$fdr < 0.05)
sum(oral_p_value$p_value < 0.05) / nrow(oral_p_value)
sum(oral_p_value$fdr < 0.05) / nrow(oral_p_value)

sum(stool_p_value$p_value < 0.05)
sum(stool_p_value$fdr < 0.05)
sum(stool_p_value$p_value < 0.05) / nrow(stool_p_value)
sum(stool_p_value$fdr < 0.05) / nrow(stool_p_value)

{
  ##load data
  load("data_analysis/nasal_microbiome/season_analysis/days_vs_microbiome/total_r2")
  nasal_total_r2 <- total_r2
  
  load("data_analysis/oral_microbiome/season_analysis/days_vs_microbiome/total_r2")
  oral_total_r2 <- total_r2
  
  load("data_analysis/stool_microbiome/season_analysis/days_vs_microbiome/total_r2")
  stool_total_r2 <- total_r2
  
  load("data_analysis/skin_microbiome/season_analysis/days_vs_microbiome/total_r2")
  skin_total_r2 <- total_r2
}

temp_data <-
  rbind(
    data.frame(body_site = "Stool",
               stool_total_r2),
    data.frame(body_site = "Skin",
               skin_total_r2),
    data.frame(body_site = "Nasal",
               nasal_total_r2),
    data.frame(body_site = "Oral",
               oral_total_r2)
  )

my_comparisons <- list(
  c("Skin", "Oral"),
  c("Skin", "Stool"),
  c("Skin", "Nasal"),
  c("Oral", "Stool"),
  c("Oral", "Nasal"),
  c("Stool", "Nasal")
)

library(ggpubr)

plot <- ggviolin(
  temp_data,
  x = "body_site",
  y = "r2",
  color = "body_site",
  palette = body_site_color[c("Skin", "Oral", "Stool", "Nasal")],
  # add = "jitter",
  add = c("boxplot", "jitter"),
  add.params = list(fill = "white"),
  legend = NULL
) +
  stat_compare_means(comparisons = my_comparisons) + # Add pairwise comparisons p-value
  stat_compare_means(label.y = 0.2)  +
  labs(x = "", y = "%") +
  base_theme

plot

ggsave(plot,
       filename = "all_asv_r2.pdf",
       width = 9,
       height = 7)





