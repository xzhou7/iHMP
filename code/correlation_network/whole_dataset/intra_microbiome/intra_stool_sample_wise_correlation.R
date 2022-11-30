#' ---
#' title: "stool microbiome stool_microbiome correlation"
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

masstools::setwd_project()
library(tidyverse)
rm(list = ls())

source("code/tools.R")

######work directory
masstools::setwd_project()
dir.create("data_analysis/correlation_network/whole_data_set/intra_stool_microbiome")
setwd("data_analysis/correlation_network/whole_data_set/intra_stool_microbiome")

###load data
###stool microbiome
{
  load(here::here("data_analysis/stool_microbiome/data_preparation/expression_data"))
  load(here::here("data_analysis/stool_microbiome/data_preparation/sample_info"))
  load(here::here("data_analysis/stool_microbiome/data_preparation/variable_info"))
}

stool_microbiome_expression_data = expression_data
stool_microbiome_sample_info = sample_info
stool_microbiome_variable_info = variable_info

###read genus table
expression_data =
  data.table::fread(here::here("data/from_xin/Genus Table/ST/Genus_ST.csv")) %>%
  as.data.frame() %>%
  tibble::column_to_rownames(var = "SampleID") %>%
  dplyr::select(-c(V1:batch)) %>%
  t() %>%
  as.data.frame()

stool_microbiome_variable_info =
  stool_microbiome_variable_info[match(rownames(expression_data), stool_microbiome_variable_info$Genus),]

stool_microbiome_variable_info$Genus == rownames(expression_data)

###remove the variables which Genus are NA
remove_idx = which(is.na(stool_microbiome_variable_info$Genus))
remove_idx
if(length(remove_idx) > 0){
  stool_microbiome_variable_info = stool_microbiome_variable_info[-remove_idx,]
  expression_data = expression_data[-remove_idx,]
}

rownames(expression_data) = stool_microbiome_variable_info$variable_id

stool_microbiome_expression_data =
  expression_data

stool_microbiome_variable_info =
  stool_microbiome_variable_info %>%
  dplyr::filter(!stringr::str_detect(Genus, "Unclassified_Bacteria"))

stool_microbiome_expression_data =
  stool_microbiome_expression_data[stool_microbiome_variable_info$variable_id,]

dim(stool_microbiome_sample_info)
dim(stool_microbiome_variable_info)
dim(stool_microbiome_expression_data)

rownames(stool_microbiome_expression_data) == stool_microbiome_variable_info$variable_id
colnames(stool_microbiome_expression_data) == stool_microbiome_sample_info$sample_id

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

######--------------------------------------------------------------------------
library(plyr)

#####
dim(stool_microbiome_expression_data)

##https://rpubs.com/DKCH2020/578881
##https://ourcodingclub.github.io/tutorials/mixed-models/
##https://zhuanlan.zhihu.com/p/63092231
##https://www.linglab.cn/knowledge/10

library(compositions)
stool_microbiome_expression_data = 
  stool_microbiome_expression_data %>% 
  purrr::map(function(x){
    x = compositions::clr(x) %>% 
      as.numeric()
    x
  }) %>% 
  do.call(cbind, .) %>% 
  as.data.frame()

rownames(stool_microbiome_expression_data) = stool_microbiome_variable_info$variable_id

##step 1
###linear mixed model to adjust the subject ID random effect, and then use the partial correlation
###to get the correlation between microbiome and stool_microbiome
library(lme4)
library(rmcorr)

library(future)
library(furrr)

# intra_stool_microbiome_lm_adjusted_cor =
#   lm_adjusted_cor(
#     data_set1 = stool_microbiome_expression_data,
#     data_set2 = stool_microbiome_expression_data,
#     sample_info = stool_microbiome_sample_info,
#     method = "all",
#     threads = 8
#   )
# 
# intra_stool_microbiome_lm_adjusted_cor =
#   intra_stool_microbiome_lm_adjusted_cor[[1]]
# 
# ###remove duplicate correlation
# name =
#   intra_stool_microbiome_lm_adjusted_cor %>%
#   apply(1, function(x) {
#     paste(sort(as.character(x[c(1:2)])), collapse = "_")
#   })
# 
# intra_stool_microbiome_lm_adjusted_cor =
# intra_stool_microbiome_lm_adjusted_cor %>%
#   dplyr::mutate(name = name) %>%
#   dplyr::distinct(name, .keep_all = TRUE)
# 
# save(
#   intra_stool_microbiome_lm_adjusted_cor,
#   file = "intra_stool_microbiome_lm_adjusted_cor",
#   compress = "xz"
# )

load("intra_stool_microbiome_lm_adjusted_cor")


intra_stool_sample_wise_dim = 
  dim(stool_microbiome_expression_data)

save(intra_stool_sample_wise_dim, file = "intra_stool_sample_wise_dim")


###here we use the lm_adjusted_cor
sum(intra_stool_microbiome_lm_adjusted_cor$p_adjust < 0.2)
