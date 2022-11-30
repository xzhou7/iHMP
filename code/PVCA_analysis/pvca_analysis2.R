##Principal Variance Component Analysis
###don't use this one, use the pvca_analysis2.R
no_source()

# set work directory
masstools::setwd_project()
library(tidyverse)
library(phyloseq)
rm(lins = ls())

source("code/tools.R")

####ggtern

######asv level
{
  ##load data 
  load(here::here("data_analysis/nasal_microbiome/PVCA_analysis/pvca_object_asv"))
  nasal_object <- pvca_object_asv
  
  load(here::here("data_analysis/skin_microbiome/PVCA_analysis/pvca_object_asv"))
  skin_object <- pvca_object_asv
  
  load(here::here("data_analysis/stool_microbiome/PVCA_analysis/pvca_object_asv"))
  stool_object <- pvca_object_asv
  
  load(here::here("data_analysis/oral_microbiome/PVCA_analysis/pvca_object_asv"))
  oral_object <- pvca_object_asv
}


temp_data <-
  rbind(nasal_object,
        stool_object,
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

ggsave(plot,
       filename = "trenery_plot_asv.pdf",
       width = 7,
       height = 7)










######class level
{
  ##load data 
  load(here::here("data_analysis/nasal_microbiome/PVCA_analysis/pvca_object_class"))
  nasal_object <- pvca_object_class
  
  load(here::here("data_analysis/skin_microbiome/PVCA_analysis/pvca_object_class"))
  skin_object <- pvca_object_class
  
  load(here::here("data_analysis/stool_microbiome/PVCA_analysis/pvca_object_class"))
  stool_object <- pvca_object_class
  
  load(here::here("data_analysis/oral_microbiome/PVCA_analysis/pvca_object_class"))
  oral_object <- pvca_object_class
}


temp_data <-
  rbind(nasal_object,
        stool_object,
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

ggsave(plot,
       filename = "trenery_plot_class.pdf",
       width = 7,
       height = 7)





######family level
{
  ##load data 
  load(here::here("data_analysis/nasal_microbiome/PVCA_analysis/pvca_object_family"))
  nasal_object <- pvca_object_family
  
  load(here::here("data_analysis/skin_microbiome/PVCA_analysis/pvca_object_family"))
  skin_object <- pvca_object_family
  
  load(here::here("data_analysis/stool_microbiome/PVCA_analysis/pvca_object_family"))
  stool_object <- pvca_object_family
  
  load(here::here("data_analysis/oral_microbiome/PVCA_analysis/pvca_object_family"))
  oral_object <- pvca_object_family
}


temp_data <-
  rbind(nasal_object,
        stool_object,
        skin_object,
        oral_object) %>%
  as.data.frame() %>%
  dplyr::rename(Subject = subject_id,
                Days = days,
                Residual = resid)
rownames(temp_data) = c("Nasal", "Stool", "Skin", "Oral")

temp_data = 
  temp_data %>% 
  tibble::rownames_to_column(var = "family")

library(ggtern)

plot <-
  ggtern::ggtern(data = temp_data, aes(x = Days, y = Subject, z = Residual)) +
  ggtern::theme_bvbw() +
  geom_point(aes(fill = family),
             shape = 21,
             size = 6,
             alpha = 0.8) +
  geom_text(aes(label = family)) +
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

ggsave(plot,
       filename = "trenery_plot_family.pdf",
       width = 7,
       height = 7)






######genus level
{
  ##load data 
  load(here::here("data_analysis/nasal_microbiome/PVCA_analysis/pvca_object_genus"))
  nasal_object <- pvca_object_genus
  
  load(here::here("data_analysis/skin_microbiome/PVCA_analysis/pvca_object_genus"))
  skin_object <- pvca_object_genus
  
  load(here::here("data_analysis/stool_microbiome/PVCA_analysis/pvca_object_genus"))
  stool_object <- pvca_object_genus
  
  load(here::here("data_analysis/oral_microbiome/PVCA_analysis/pvca_object_genus"))
  oral_object <- pvca_object_genus
}


temp_data <-
  rbind(nasal_object,
        stool_object,
        skin_object,
        oral_object) %>%
  as.data.frame() %>%
  dplyr::rename(Subject = subject_id,
                Days = days,
                Residual = resid)
rownames(temp_data) = c("Nasal", "Stool", "Skin", "Oral")

temp_data = 
  temp_data %>% 
  tibble::rownames_to_column(var = "genus")

library(ggtern)

plot <-
  ggtern::ggtern(data = temp_data, aes(x = Days, y = Subject, z = Residual)) +
  ggtern::theme_bvbw() +
  geom_point(aes(fill = genus),
             shape = 21,
             size = 6,
             alpha = 0.8) +
  geom_text(aes(label = genus)) +
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

ggsave(plot,
       filename = "trenery_plot_genus.pdf",
       width = 7,
       height = 7)









######order level
{
  ##load data 
  load(here::here("data_analysis/nasal_microbiome/PVCA_analysis/pvca_object_order"))
  nasal_object <- pvca_object_order
  
  load(here::here("data_analysis/skin_microbiome/PVCA_analysis/pvca_object_order"))
  skin_object <- pvca_object_order
  
  load(here::here("data_analysis/stool_microbiome/PVCA_analysis/pvca_object_order"))
  stool_object <- pvca_object_order
  
  load(here::here("data_analysis/oral_microbiome/PVCA_analysis/pvca_object_order"))
  oral_object <- pvca_object_order
}


temp_data <-
  rbind(nasal_object,
        stool_object,
        skin_object,
        oral_object) %>%
  as.data.frame() %>%
  dplyr::rename(Subject = subject_id,
                Days = days,
                Residual = resid)
rownames(temp_data) = c("Nasal", "Stool", "Skin", "Oral")

temp_data = 
  temp_data %>% 
  tibble::rownames_to_column(var = "order")

library(ggtern)

plot <-
  ggtern::ggtern(data = temp_data, aes(x = Days, y = Subject, z = Residual)) +
  ggtern::theme_bvbw() +
  geom_point(aes(fill = order),
             shape = 21,
             size = 6,
             alpha = 0.8) +
  geom_text(aes(label = order)) +
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

ggsave(plot,
       filename = "trenery_plot_order.pdf",
       width = 7,
       height = 7)






######phylum level
{
  ##load data 
  load(here::here("data_analysis/nasal_microbiome/PVCA_analysis/pvca_object_phylum"))
  nasal_object <- pvca_object_phylum
  
  load(here::here("data_analysis/skin_microbiome/PVCA_analysis/pvca_object_phylum"))
  skin_object <- pvca_object_phylum
  
  load(here::here("data_analysis/stool_microbiome/PVCA_analysis/pvca_object_phylum"))
  stool_object <- pvca_object_phylum
  
  load(here::here("data_analysis/oral_microbiome/PVCA_analysis/pvca_object_phylum"))
  oral_object <- pvca_object_phylum
}


temp_data <-
  rbind(nasal_object,
        stool_object,
        skin_object,
        oral_object) %>%
  as.data.frame() %>%
  dplyr::rename(Subject = subject_id,
                Days = days,
                Residual = resid)
rownames(temp_data) = c("Nasal", "Stool", "Skin", "Oral")

temp_data = 
  temp_data %>% 
  tibble::rownames_to_column(var = "phylum")

library(ggtern)

plot <-
  ggtern::ggtern(data = temp_data, aes(x = Days, y = Subject, z = Residual)) +
  ggtern::theme_bvbw() +
  geom_point(aes(fill = phylum),
             shape = 21,
             size = 6,
             alpha = 0.8) +
  geom_text(aes(label = phylum)) +
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

ggsave(plot,
       filename = "trenery_plot_phylum.pdf",
       width = 7,
       height = 7)
