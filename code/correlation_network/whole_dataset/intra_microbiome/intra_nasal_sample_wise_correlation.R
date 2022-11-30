#' ---
#' title: "nasal microbiome nasal_microbiome correlation"
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
dir.create("data_analysis/correlation_network/whole_data_set/intra_nasal_microbiome")
setwd("data_analysis/correlation_network/whole_data_set/intra_nasal_microbiome")

###load data
###nasal microbiome
{
  load(here::here("data_analysis/nasal_microbiome/data_preparation/expression_data"))
  load(here::here("data_analysis/nasal_microbiome/data_preparation/sample_info"))
  load(here::here("data_analysis/nasal_microbiome/data_preparation/variable_info"))
}

nasal_microbiome_expression_data = expression_data
nasal_microbiome_sample_info = sample_info
nasal_microbiome_variable_info = variable_info

###read genus table
expression_data =
  data.table::fread(here::here("data/from_xin/Genus Table/NS/Genus_NS.csv")) %>%
  as.data.frame() %>%
  tibble::column_to_rownames(var = "SampleID") %>%
  dplyr::select(-c(V1:batch)) %>%
  t() %>%
  as.data.frame()

nasal_microbiome_variable_info =
  nasal_microbiome_variable_info[match(rownames(expression_data), nasal_microbiome_variable_info$Genus),]

nasal_microbiome_variable_info$Genus == rownames(expression_data)

###remove the variables which Genus are NA
remove_idx = which(is.na(nasal_microbiome_variable_info$Genus))
remove_idx
if(length(remove_idx) > 0){
  nasal_microbiome_variable_info = nasal_microbiome_variable_info[-remove_idx,]
  expression_data = expression_data[-remove_idx,]
}

rownames(expression_data) = nasal_microbiome_variable_info$variable_id

nasal_microbiome_expression_data =
  expression_data

nasal_microbiome_variable_info =
  nasal_microbiome_variable_info %>%
  dplyr::filter(!stringr::str_detect(Genus, "Unclassified_Bacteria"))

nasal_microbiome_expression_data =
  nasal_microbiome_expression_data[nasal_microbiome_variable_info$variable_id,]

dim(nasal_microbiome_sample_info)
dim(nasal_microbiome_variable_info)
dim(nasal_microbiome_expression_data)

rownames(nasal_microbiome_expression_data) == nasal_microbiome_variable_info$variable_id
colnames(nasal_microbiome_expression_data) == nasal_microbiome_sample_info$sample_id

###only remain the subjects with at least >= 5
remian_subject_id =
  nasal_microbiome_sample_info %>%
  dplyr::group_by(subject_id) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::ungroup() %>%
  dplyr::filter(n > 5) %>%
  dplyr::pull(subject_id)

nasal_microbiome_sample_info =
  nasal_microbiome_sample_info %>%
  dplyr::filter(subject_id %in% remian_subject_id)

nasal_microbiome_expression_data =
  nasal_microbiome_expression_data[,nasal_microbiome_sample_info$sample_id]

##only remain the genus at least in 10% subjects
remain_idx =
  which(rowSums(nasal_microbiome_expression_data) > 0)

nasal_microbiome_expression_data = nasal_microbiome_expression_data[remain_idx,]
nasal_microbiome_variable_info = nasal_microbiome_variable_info[remain_idx,,drop = FALSE]

remain_idx =
  nasal_microbiome_expression_data %>%
  apply(1, function(x){
    sum(as.numeric(x) == 0) / ncol(nasal_microbiome_expression_data)
  }) %>%
  `<`(0.9) %>%
  which()

length(remain_idx)

nasal_microbiome_expression_data = nasal_microbiome_expression_data[remain_idx,]
nasal_microbiome_variable_info = nasal_microbiome_variable_info[remain_idx,,drop = FALSE]

######--------------------------------------------------------------------------
library(plyr)

#####
dim(nasal_microbiome_expression_data)

##https://rpubs.com/DKCH2020/578881
##https://ourcodingclub.github.io/tutorials/mixed-models/
##https://zhuanlan.zhihu.com/p/63092231
##https://www.linglab.cn/knowledge/10

library(compositions)
nasal_microbiome_expression_data = 
  nasal_microbiome_expression_data %>% 
  purrr::map(function(x){
    x = compositions::clr(x) %>% 
      as.numeric()
    x
  }) %>% 
  do.call(cbind, .) %>% 
  as.data.frame()

rownames(nasal_microbiome_expression_data) = nasal_microbiome_variable_info$variable_id

##step 1
###linear mixed model to adjust the subject ID random effect, and then use the partial correlation
###to get the correlation between microbiome and nasal_microbiome
library(lme4)
library(rmcorr)

library(future)
library(furrr)

# intra_nasal_microbiome_lm_adjusted_cor =
#   lm_adjusted_cor(
#     data_set1 = nasal_microbiome_expression_data,
#     data_set2 = nasal_microbiome_expression_data,
#     sample_info = nasal_microbiome_sample_info,
#     method = "all",
#     threads = 8
#   )
# 
# intra_nasal_microbiome_lm_adjusted_cor =
#   intra_nasal_microbiome_lm_adjusted_cor[[1]]
# 
# ###remove duplicate correlation
# name =
#   intra_nasal_microbiome_lm_adjusted_cor %>%
#   apply(1, function(x) {
#     paste(sort(as.character(x[c(1:2)])), collapse = "_")
#   })
# 
# intra_nasal_microbiome_lm_adjusted_cor =
# intra_nasal_microbiome_lm_adjusted_cor %>%
#   dplyr::mutate(name = name) %>%
#   dplyr::distinct(name, .keep_all = TRUE)
# 
# save(
#   intra_nasal_microbiome_lm_adjusted_cor,
#   file = "intra_nasal_microbiome_lm_adjusted_cor",
#   compress = "xz"
# )

intra_nasal_sample_wise_dim = 
dim(nasal_microbiome_expression_data)

save(intra_nasal_sample_wise_dim, file = "intra_nasal_sample_wise_dim")

load("intra_nasal_microbiome_lm_adjusted_cor")

###here we use the lm_adjusted_cor
sum(intra_nasal_microbiome_lm_adjusted_cor$p_adjust < 0.2)
