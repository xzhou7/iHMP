no_source()

# set work directory
setwd(masstools::get_project_wd())
library(tidyverse)
library(phyloseq)
rm(list = ls())

source("code/tools.R")

{
  load("data_analysis/stool_microbiome/stool_microbiome_diet/total_r2")
  stool_microbiome_total_r2 = total_r2
  
  load("data_analysis/skin_microbiome/skin_microbiome_diet/total_r2")
  skin_microbiome_total_r2 = total_r2
  
  load("data_analysis/oral_microbiome/oral_microbiome_diet/total_r2")
  oral_microbiome_total_r2 = total_r2
  
  load("data_analysis/nasal_microbiome/nasal_microbiome_diet/total_r2")
  nasal_microbiome_total_r2 = total_r2
}

####plot to show
setwd(masstools::get_project_wd())
setwd("data_analysis/diet_affect_microbiome")

library(ggpubr)

temp_data =
  rbind(
    data.frame(body_site = "Skin",
               r2 = skin_microbiome_total_r2),
    data.frame(body_site = "Stool",
               r2 = stool_microbiome_total_r2),
    data.frame(body_site = "Nasal",
               r2 = nasal_microbiome_total_r2),
    data.frame(body_site = "Oral",
               r2 = oral_microbiome_total_r2)
  ) %>%
  dplyr::mutate(body_site = factor(body_site, levels = c("Skin", "Oral", "Stool", "Nasal")))

my_comparisons <- list(
  c("Skin", "Oral"),
  c("Skin", "Stool"),
  c("Skin", "Nasal"),
  c("Oral", "Stool"),
  c("Oral", "Nasal"),
  c("Stool", "Nasal")
)

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
  stat_compare_means(label.y = 1)  +
  labs(x = "", y = "R2") +
  base_theme

plot

# ggsave(plot,
#        filename = "body_site_r2.pdf",
#        width = 9,
#        height = 7)

# write.csv(temp_data, "data.csv", row.names = FALSE)
