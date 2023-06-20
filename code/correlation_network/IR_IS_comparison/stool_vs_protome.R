#' ---
#' title: "Stool microbiome proteome correlation"
#' author:
#'   - name: "Xiaotao Shen"
#'     url: https://www.shenxt.info/
#'     affiliation: Stanford School of Medicine
#' date: "`r Sys.Date()`"
#' site: distill::distill_website
#' output:
#'   distill::distill_article:
#'     code_folding: false
#' ---

#+ r setup, echo=TRUE, eval = TRUE, include = TRUE

no_function()
# set work directory

setwd(masstools::get_project_wd())
library(tidyverse)
library(ggraph)
library(igraph)
library(tidygraph)

###example network
library(plyr)
rm(list = ls())

source("code/tools.R")

######work directory
setwd(masstools::get_project_wd())

data_IR <-
  purrr::map(1:20, function(i) {
    load(
      file.path(
        "data_analysis/correlation_network/whole_data_set_IR/stool_microbiome_vs_proteome_IR_permutation/",
        paste("edge_data", i, sep = "_")
      )
    )
    load(
      file.path(
        "data_analysis/correlation_network/whole_data_set_IR/stool_microbiome_vs_proteome_IR_permutation/",
        paste("node_data", i, sep = "_")
      )
    )
    data.frame(
      node_number = nrow(node_data),
      edge_number = nrow(edge_data),
      class = "IR"
    )
    
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()

data_IS <-
  purrr::map(1:20, function(i) {
    load(
      file.path(
        "data_analysis/correlation_network/whole_data_set_IS/stool_microbiome_vs_proteome_IS_permutation/",
        paste("edge_data", i, sep = "_")
      )
    )
    load(
      file.path(
        "data_analysis/correlation_network/whole_data_set_IS/stool_microbiome_vs_proteome_IS_permutation/",
        paste("node_data", i, sep = "_")
      )
    )
    data.frame(
      node_number = nrow(node_data),
      edge_number = nrow(edge_data),
      class = "IS"
    )
    
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()

dir.create("data_analysis/correlation_network/IR_vs_IS/stool_vs_proteome",
           recursive = TRUE)
setwd("data_analysis/correlation_network/IR_vs_IS/stool_vs_proteome")

temp_data <-
  rbind(data_IR,
        data_IS)

library(ggstatsplot)

plot1 <-
  ggbetweenstats(
    data  = temp_data,
    x     = class,
    y     = node_number,
    title = "Node number",
    p.adjust.method = "fdr",
    test.type =
  ) +
  ggplot2::scale_fill_manual(values = iris_color) +
  ggplot2::scale_color_manual(values = iris_color) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "", y = "Node number")

plot1

plot2 <-
  ggbetweenstats(
    data  = temp_data,
    x     = class,
    y     = edge_number,
    title = "Edge number",
    p.adjust.method = "fdr",
    test.type =
  ) +
  ggplot2::scale_fill_manual(values = iris_color) +
  ggplot2::scale_color_manual(values = iris_color) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "", y = "Edge number")

plot2

library(patchwork)

plot <- 
plot1 + plot2

ggsave(plot, filename = "stool_proteome.pdf", width = 14, height = 7)
