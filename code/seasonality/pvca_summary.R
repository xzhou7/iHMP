no_source()

# set work directory
setwd(masstools::get_project_wd())
library(tidyverse)
library(phyloseq)
rm(list = ls())

source("code/tools.R")

{
  ####load data
  ##nasal
  load("data_analysis/nasal_microbiome/season_analysis/remain_idx")
  ns_remain_idx = remain_idx
  load("data_analysis/nasal_microbiome/season_analysis/p_value_adj_all")
  ns_p_value_adj_all = p_value_adj_all
  load("data_analysis/nasal_microbiome/season_analysis/p_value_all")
  ns_p_value_all = p_value_all
  
  ##stool
  load("data_analysis/stool_microbiome/season_analysis/remain_idx")
  st_remain_idx = remain_idx
  load("data_analysis/stool_microbiome/season_analysis/p_value_adj_all")
  st_p_value_adj_all = p_value_adj_all
  load("data_analysis/stool_microbiome/season_analysis/p_value_all")
  st_p_value_all = p_value_all
  
  ##skin
  load("data_analysis/skin_microbiome/season_analysis/remain_idx")
  skin_remain_idx = remain_idx
  load("data_analysis/skin_microbiome/season_analysis/p_value_adj_all")
  skin_p_value_adj_all = p_value_adj_all
  load("data_analysis/skin_microbiome/season_analysis/p_value_all")
  skin_p_value_all = p_value_all
  
  ##oral
  load("data_analysis/oral_microbiome/season_analysis/remain_idx")
  oral_remain_idx = remain_idx
  load("data_analysis/oral_microbiome/season_analysis/p_value_adj_all")
  oral_p_value_adj_all = p_value_adj_all
  load("data_analysis/oral_microbiome/season_analysis/p_value_all")
  oral_p_value_all = p_value_all
}

{
  ##load data
  load("data_analysis/nasal_microbiome/season_analysis/pvca_object")
  ns_object <- pvca_object
  
  load("data_analysis/skin_microbiome/season_analysis/pvca_object")
  skin_object <- pvca_object
  
  load("data_analysis/stool_microbiome/season_analysis/pvca_object")
  st_object <- pvca_object
  
  load("data_analysis/oral_microbiome/season_analysis/pvca_object")
  oral_object <- pvca_object
}

{
  load("data_analysis/nasal_microbiome/nasal_microbiome_exposome_chemical/total_r2")
  exposome_chemical_nasal_microbiome_total_r2 = total_r2
  
  load("data_analysis/stool_microbiome/stool_microbiome_exposome_chemical/total_r2")
  exposome_chemical_stool_microbiome_total_r2 = total_r2
  
  load("data_analysis/skin_microbiome/skin_microbiome_exposome_chemical/total_r2")
  exposome_chemical_skin_microbiome_total_r2 = total_r2
  
  load("data_analysis/oral_microbiome/oral_microbiome_exposome_chemical/total_r2")
  exposome_chemical_oral_microbiome_total_r2 = total_r2
}

length(ns_remain_idx)
length(st_remain_idx)
length(skin_remain_idx)
length(oral_remain_idx)

sum(ns_p_value_all < 0.05)
sum(ns_p_value_adj_all < 0.05)

sum(st_p_value_all < 0.05)
sum(st_p_value_adj_all < 0.05)

sum(skin_p_value_all < 0.05)
sum(skin_p_value_adj_all < 0.05)

sum(oral_p_value_all < 0.05)
sum(oral_p_value_adj_all < 0.05)

sum(ns_p_value_all < 0.05) / length(ns_remain_idx)
sum(ns_p_value_adj_all < 0.05) / length(ns_remain_idx)

sum(st_p_value_all < 0.05) / length(st_remain_idx)
sum(st_p_value_adj_all < 0.05) / length(st_remain_idx)

sum(skin_p_value_all < 0.05) / length(skin_remain_idx)
sum(skin_p_value_adj_all < 0.05) / length(skin_remain_idx)

sum(oral_p_value_all < 0.05) / length(oral_remain_idx)
sum(oral_p_value_adj_all < 0.05) / length(oral_remain_idx)

####plot to show
setwd(masstools::get_project_wd())
setwd("data_analysis/season_analysis/")
library(scatterpie)

temp_data =
  data.frame(
    body_site = c("Nasal", "Stool", "Skin", "Oral"),
    p = c(12.1, 5.7, 8.03, 34.4),
    p.adjust = c(1.97, 0.3, 1.41, 24.3),
    total_asv = log(c(52460, 52460, 8947, 8947), 10) / 5,
    remained_asv = c(916, 3569, 710, 584),
    removed_asv = c(52460, 52460, 8947, 8947) - c(916, 3569, 710, 584)
  )

plot =
  temp_data %>%
  ggplot(aes(x = p, p.adjust)) +
  geom_scatterpie(
    aes(
      x = p,
      y = p.adjust,
      group = body_site,
      r = total_asv
    ),
    data = temp_data,
    cols = c("remained_asv", "removed_asv"),
    color = NA,
    alpha = .8
  ) +
  ggrepel::geom_label_repel(aes(label = body_site),
                            size = 4) +
  scale_color_manual(values = body_site_color) +
  scale_fill_manual(values = c(
    "remained_asv" = "red",
    "removed_asv" = "grey"
  )) +
  labs(x = "Percent of ASV with p < 0.05 (%)",
       y = "Percent of ASV with p.adj < 0.05 (%)") +
  base_theme +
  coord_equal() +
  geom_scatterpie_legend(
    temp_data$total_asv,
    x = 10,
    y = 20,
    labeller = function(x)
      round((x * 5) ^ 10)
  )

plot

# ggsave(plot, filename = "season_stability.pdf", width = 9, height = 7)

temp_data <-
  rbind(ns_object,
        st_object,
        skin_object,
        oral_object) %>%
  as.data.frame() %>%
  dplyr::rename(Subject = subject_id,
                Days = days,
                Residual = resid)
rownames(temp_data) = c("Nasal", "Stool", "Skin", "Oral")
temp_data =
  temp_data %>%
  tibble::rownames_to_column(var = "class")

library(ggtern)

plot <-
  ggtern::ggtern(data = temp_data, aes(x = Days, y = Subject, z = Residual)) +
  ggtern::theme_bvbw() +
  geom_point(aes(fill = class),
             shape = 21,
             size = 6,
             alpha = 0.8) +
  geom_text(aes(label = class)) +
  scale_fill_manual(values = body_site_color) +
  labs(x = "Days",
       y = "Subject",
       z = "Residual",
       title = "") +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    plot.background = element_rect(fill = "transparent"),
    panel.background = element_rect(fill = "transparent")
  )

plot

# ggsave(plot, filename = "trenery_plot.pdf", width = 7, height = 7)


####exposome on microbiome
library(ggpubr)
wilcox.test(
  exposome_chemical_nasal_microbiome_total_r2,
  exposome_chemical_skin_microbiome_total_r2
)

temp_data =
  rbind(
    data.frame(body_site = "Skin",
               r2 = exposome_chemical_skin_microbiome_total_r2),
    data.frame(body_site = "Stool",
               r2 = exposome_chemical_stool_microbiome_total_r2),
    data.frame(body_site = "Nasal",
               r2 = exposome_chemical_nasal_microbiome_total_r2),
    data.frame(body_site = "Oral",
               r2 = exposome_chemical_oral_microbiome_total_r2)
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
  stat_compare_means(label.y = 1.3)  +
  labs(x = "", y = "R2") +
  base_theme

plot

# ggsave(plot,
#        filename = "body_site_r2.pdf",
#        width = 9,
#        height = 7)
