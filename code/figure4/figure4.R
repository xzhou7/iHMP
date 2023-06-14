no_source()

setwd(masstools::get_project_wd())
rm(list = ls())
library(tidyverse)

source("code/tools.R")

stool <-
  readr::read_csv("Figures/Figure4/detailed/figure4E/stool.taxa.cytokine.csv")

skin <-
  readr::read_csv("Figures/Figure4/detailed/figure4E/skin.taxa.cytokine.csv")

oral <-
  readr::read_csv("Figures/Figure4/detailed/figure4E/oral.taxa.cytokine.csv")

nasal <-
  readr::read_csv("Figures/Figure4/detailed/figure4E/nasal.taxa.cytokine.csv")

setwd("data_analysis/figure4")

temp_data <-
  rbind(
    stool %>% dplyr::select(-`...1`) %>% mutate(class = "Stool"),
    skin %>% dplyr::select(-`...1`) %>% mutate(class = "Skin"),
    oral %>% dplyr::select(-`...1`) %>% mutate(class = "Oral"),
    nasal %>% dplyr::select(-`...1`) %>% mutate(class = "Nasal")
  )

library(tidyverse)
library(plyr)

colnames(temp_data)[2] <- "Cytokine"

plot <-
  temp_data %>%
  pivot_longer(
    cols = c(Cytokine, Total),
    names_to = "class2",
    values_to = "number"
  ) %>%
  plyr::dlply(.variables = .(class, class2)) %>%
  purrr::map(function(x) {
    x$number <-
      x$number * 100 / sum(x$number)
    x
  }) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(class = factor(class,
                               levels = c("Stool", "Skin", "Oral", "Nasal"))) %>%
  dplyr::mutate(class2 = factor(class2,
                                levels = c("Total", "Cytokine"))) %>%
  dplyr::mutate(phylum = factor(
    phylum,
    levels = c(
      "Actinobacteria",
      "Bacteroidetes",
      "Firmicutes",
      "Proteobacteria",
      "other"
    )
  )) %>%
  ggplot(aes(class2, number)) +
  geom_bar(aes(fill = phylum),
           color = "black",
           stat = "identity") +
  facet_grid(cols = vars(class)) +
  scale_fill_manual(values = c(phylum_color[c("Actinobacteria",
                                              "Bacteroidetes",
                                              "Firmicutes",
                                              "Proteobacteria")],
                               "other" = "grey")) +
  base_theme +
  labs(x = "", y = "Percentage (%)") +
  theme(axis.text.x = element_text(
    angle = 45,
    hjust = 1,
    vjust = 1
  )) +
  # scale_x_discrete(expand = expansion(mult = c(0.05,0.05))) +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.01)))

plot

# ggsave(plot,
#        filename = "figure4c.pdf",
#        width = 5,
#        height = 3)
