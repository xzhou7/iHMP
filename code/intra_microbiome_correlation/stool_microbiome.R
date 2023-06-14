#' ---
#' title: "Stool microbiome metabolome correlation"
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

####load data
# ###load diversity
# load("data/from_xin/Diversity_Datatable.RData")
# 
# ##prevalence
# load("data/from_xin/Prevalance.RData")

###stool microbiome
{
  load("data_analysis/stool_microbiome/data_preparation/expression_data")
  load("data_analysis/stool_microbiome/data_preparation/sample_info")
  load("data_analysis/stool_microbiome/data_preparation/variable_info")  
}

stool_microbiome_expression_data = expression_data
stool_microbiome_sample_info = sample_info
stool_microbiome_variable_info = variable_info

###read genus table
expression_data =
readr::read_csv(here::here("data/from_xin/Genus Table/ST/Genus_ST.csv")) %>%
  tibble::column_to_rownames(var = "SampleID") %>%
  dplyr::select(-c(X1:batch)) %>%
  t() %>%
  as.data.frame()

stool_microbiome_variable_info =
  data.frame(variable_id = rownames(expression_data))

stool_microbiome_expression_data =
  expression_data

dim(stool_microbiome_sample_info)
length(unique(stool_microbiome_sample_info$subject_id))
dim(stool_microbiome_variable_info)

######work directory
setwd(masstools::get_project_wd())
setwd("data_analysis/intra_microbiome_correlation/stool_microbiome")

###only remain the subjects with at least >= 5
remian_subject_id = 
stool_microbiome_sample_info %>% 
  dplyr::group_by(subject_id) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(n > 5) %>% 
  dplyr::pull(subject_id)

stool_microbiome_sample_info = 
  stool_microbiome_sample_info %>% 
  dplyr::filter(subject_id %in% remian_subject_id)

stool_microbiome_expression_data = 
  stool_microbiome_expression_data[,stool_microbiome_sample_info$sample_id]

##only remain the genus at least in 10% subjects
remain_idx =
  which(rowSums(stool_microbiome_expression_data) > 0)

stool_microbiome_expression_data = stool_microbiome_expression_data[remain_idx,]
stool_microbiome_variable_info = stool_microbiome_variable_info[remain_idx,,drop = FALSE]

remain_idx =
  stool_microbiome_expression_data %>%
  apply(1, function(x){
    sum(as.numeric(x) == 0) / ncol(stool_microbiome_expression_data)
  }) %>%
  `<`(0.9) %>%
  which()

length(remain_idx)

stool_microbiome_expression_data = stool_microbiome_expression_data[remain_idx,]
stool_microbiome_variable_info = stool_microbiome_variable_info[remain_idx,,drop = FALSE]

##save data
{
  save(stool_microbiome_expression_data, file = "stool_microbiome_expression_data")
  save(stool_microbiome_variable_info, file = "stool_microbiome_variable_info")
  save(stool_microbiome_sample_info, file = "stool_microbiome_sample_info")
}

{
  load("stool_microbiome_expression_data")
  load("stool_microbiome_variable_info")
  load("stool_microbiome_sample_info")
}

