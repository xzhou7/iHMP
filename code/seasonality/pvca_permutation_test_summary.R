no_source()

# set work directory
setwd(masstools::get_project_wd())
library(tidyverse)
library(phyloseq)
rm(list = ls())

source("code/tools.R")

{
  ##load data
  load(
    "data_analysis/nasal_microbiome/season_analysis/permutation_test/pvca_object_list"
  )
  nasal_object <- pvca_object_list
  
  load(
    "data_analysis/oral_microbiome/season_analysis/permutation_test/pvca_object_list"
  )
  oral_object <- pvca_object_list
  
  load(
    "data_analysis/skin_microbiome/season_analysis/permutation_test/pvca_object_list"
  )
  skin_object <- pvca_object_list
  
  load(
    "data_analysis/stool_microbiome/season_analysis/permutation_test/pvca_object_list"
  )
  stool_object <- pvca_object_list
}

####plot to show
setwd(masstools::get_project_wd())
dir.create("data_analysis/season_analysis/permutation_test")
setwd("data_analysis/season_analysis/permutation_test")
library(scatterpie)

temp_data <-
  rbind(
    data.frame(body_site = "Stool",
               days = stool_object$days),
    data.frame(body_site = "Skin",
               days = skin_object$days),
    data.frame(body_site = "Oral",
               days = oral_object$days),
    data.frame(body_site = "Nasal",
               days = nasal_object$days)
  )

temp_data$days <- temp_data$days * 100

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
  y = "days",
  color = "body_site",
  palette = body_site_color[c("Skin", "Oral", "Stool", "Nasal")],
  # add = "jitter",
  add = c("boxplot", "jitter"),
  add.params = list(fill = "white"),
  legend = NULL
) +
  stat_compare_means(comparisons = my_comparisons) + # Add pairwise comparisons p-value
  stat_compare_means(label.y = 7.5)  +
  labs(x = "", y = "%") +
  base_theme

plot

ggsave(plot,
       filename = "pvca_permutation_test.pdf",
       width = 9,
       height = 7)
