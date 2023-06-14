no_source()

# set work directory
setwd(masstools::get_project_wd())
library(tidyverse)
library(phyloseq)
rm(list = ls())

source("code/tools.R")

{
  load(
    "data_analysis/cytokine_affect_microbiome/IR/cytokine_nasal_microbiome/total_r2"
  )
  ir_nasal_microbiome_total_r2 = total_r2
  
  load(
    "data_analysis/cytokine_affect_microbiome/IR/cytokine_oral_microbiome/total_r2"
  )
  ir_oral_microbiome_total_r2 = total_r2
  
  load(
    "data_analysis/cytokine_affect_microbiome/IR/cytokine_skin_microbiome/total_r2"
  )
  ir_skin_microbiome_total_r2 = total_r2
  
  load(
    "data_analysis/cytokine_affect_microbiome/IR/cytokine_stool_microbiome/total_r2"
  )
  ir_stool_microbiome_total_r2 = total_r2
}

{
  load(
    "data_analysis/cytokine_affect_microbiome/IS/cytokine_nasal_microbiome/total_r2"
  )
  is_nasal_microbiome_total_r2 = total_r2
  
  load(
    "data_analysis/cytokine_affect_microbiome/IS/cytokine_oral_microbiome/total_r2"
  )
  is_oral_microbiome_total_r2 = total_r2
  
  load(
    "data_analysis/cytokine_affect_microbiome/IS/cytokine_skin_microbiome/total_r2"
  )
  is_skin_microbiome_total_r2 = total_r2
  
  load(
    "data_analysis/cytokine_affect_microbiome/IS/cytokine_stool_microbiome/total_r2"
  )
  is_stool_microbiome_total_r2 = total_r2
}

####plot to show
setwd(masstools::get_project_wd())
setwd("data_analysis/cytokine_affect_microbiome/")

library(ggpubr)

temp_data_ir =
  rbind(
    data.frame(body_site = "Skin",
               r2 = ir_skin_microbiome_total_r2),
    data.frame(body_site = "Stool",
               r2 = ir_stool_microbiome_total_r2),
    data.frame(body_site = "Nasal",
               r2 = ir_nasal_microbiome_total_r2),
    data.frame(body_site = "Oral",
               r2 = ir_oral_microbiome_total_r2)
  ) %>%
  dplyr::mutate(body_site = factor(body_site, levels = c("Stool", "Skin", "Oral", "Nasal")))

my_comparisons <- list(
  c("Skin", "Oral"),
  c("Skin", "Stool"),
  c("Skin", "Nasal"),
  c("Oral", "Stool"),
  c("Oral", "Nasal"),
  c("Stool", "Nasal")
)

plot_ir <- ggviolin(
  temp_data_ir,
  x = "body_site",
  y = "r2",
  color = "body_site",
  palette = body_site_color[c("Skin", "Oral", "Stool", "Nasal")],
  # add = "jitter",
  add = c("boxplot", "jitter"),
  add.params = list(fill = "white"),
  legend = NULL
) +
  # stat_compare_means(comparisons = my_comparisons) + # Add pairwise comparisons p-value
  # stat_compare_means(label.y = 1)  +
  labs(x = "", y = "R2") +
  scale_y_continuous(limits = c(0, 0.7)) +
  base_theme

plot_ir

temp_data_is =
  rbind(
    data.frame(body_site = "Skin",
               r2 = is_skin_microbiome_total_r2),
    data.frame(body_site = "Stool",
               r2 = is_stool_microbiome_total_r2),
    data.frame(body_site = "Nasal",
               r2 = is_nasal_microbiome_total_r2),
    data.frame(body_site = "Oral",
               r2 = is_oral_microbiome_total_r2)
  ) %>%
  dplyr::mutate(body_site = factor(body_site, levels = c("Stool", "Skin", "Oral", "Nasal")))

my_comparisons <- list(
  c("Skin", "Oral"),
  c("Skin", "Stool"),
  c("Skin", "Nasal"),
  c("Oral", "Stool"),
  c("Oral", "Nasal"),
  c("Stool", "Nasal")
)

plot_is <- ggviolin(
  temp_data_is,
  x = "body_site",
  y = "r2",
  color = "body_site",
  palette = body_site_color[c("Skin", "Oral", "Stool", "Nasal")],
  # add = "jitter",
  add = c("boxplot", "jitter"),
  add.params = list(fill = "white"),
  legend = NULL
) +
  # stat_compare_means(comparisons = my_comparisons) + # Add paiswise comparisons p-value
  # stat_compare_means(label.y = 1)  +
  labs(x = "", y = "R2") +
  scale_y_continuous(limits = c(0, 0.7)) +
  base_theme

plot_is

library(patchwork)

temp_data <-
  rbind(data.frame(temp_data_ir, IRIS = "IR"),
        data.frame(temp_data_is, IRIS = "IS")) %>%
  dplyr::filter(body_site %in% c("Stool", "Oral"))

plot <- 
temp_data %>%
  dplyr::mutate(body_site = as.character(body_site)) %>%
  dplyr::mutate(body_site = factor(body_site, levels = c("Stool", "Oral"))) %>%
  dplyr::mutate(IRIS = factor(IRIS, levels = c("IS", "IR"))) %>%
  ggplot(aes(x = IRIS, y = r2)) +
  geom_boxplot(outlier.shape = NA) +
  geom_dotplot(binaxis = "y", stackdir = "center") +
  facet_grid(cols = vars(body_site)) +
  theme_bw() +
  labs(x = "", y = "R2")


# ggsave(plot, filename = "summary_plot.pdf",
#        width = 9, height = 7)
# 
# write.csv(temp_data, file = "data_for_figures.csv", row.names = FALSE)
