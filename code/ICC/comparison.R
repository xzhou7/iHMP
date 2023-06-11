#' ---
#' title: "nasal microbiome metabolome correlation"
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
rm(list = ls())
source("code/tools.R")

load("data_analysis/ICC_without_confounder/icc_data")
icc_data_without_confounder <- icc_data
load("data_analysis/ICC/icc_data")

dim(icc_data_without_confounder)
dim(icc_data)

temp_data <-
  icc_data_without_confounder %>%
  dplyr::inner_join(icc_data, by = c("dt", "Variables", "variable", "level"))

plot <- 
temp_data %>% 
  ggplot(aes(value.x, value.y)) +
  geom_point() +
  theme_bw() +
  labs(x = "ICC without adjust counfounders",
       x = "ICC") +
  geom_abline(color = "red")

dir.create("data_analysis/ICC_without_confounder/")
setwd("data_analysis/ICC_without_confounder/")

cor(temp_data$value.x,
    temp_data$value.y)
ggsave(plot, filename = "compaarison.pdf", width = 7, height = 7)

