#' ---
#' title: "skin microbiome skin_microbiome correlation"
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
dir.create("data_analysis/correlation_network/whole_data_set/intra_skin_microbiome")
setwd("data_analysis/correlation_network/whole_data_set/intra_skin_microbiome")

###load data
###skin microbiome
{
  load(here::here("data_analysis/skin_microbiome/data_preparation/expression_data"))
  load(here::here("data_analysis/skin_microbiome/data_preparation/sample_info"))
  load(here::here("data_analysis/skin_microbiome/data_preparation/variable_info"))
}

skin_microbiome_expression_data = expression_data
skin_microbiome_sample_info = sample_info
skin_microbiome_variable_info = variable_info

##read genus table
expression_data =
  data.table::fread(here::here("data/from_xin/Genus Table/SK/Genus_SK.csv")) %>%
  as.data.frame() %>%
  tibble::column_to_rownames(var = "SampleID") %>%
  dplyr::select(-c(V1:SubjectID)) %>%
  t() %>%
  as.data.frame()

skin_microbiome_variable_info =
  skin_microbiome_variable_info[match(rownames(expression_data), skin_microbiome_variable_info$Genus),]

skin_microbiome_variable_info$Genus == rownames(expression_data)

###remove the variables which Genus are NA
remove_idx = which(is.na(skin_microbiome_variable_info$Genus))
remove_idx
if(length(remove_idx) > 0){
  skin_microbiome_variable_info = skin_microbiome_variable_info[-remove_idx,]
  expression_data = expression_data[-remove_idx,]
}

rownames(expression_data) = skin_microbiome_variable_info$variable_id

skin_microbiome_expression_data =
  expression_data

skin_microbiome_variable_info =
  skin_microbiome_variable_info %>%
  dplyr::filter(!stringr::str_detect(Genus, "Unclassified_Bacteria"))

skin_microbiome_expression_data =
  skin_microbiome_expression_data[skin_microbiome_variable_info$variable_id,]

dim(skin_microbiome_sample_info)
dim(skin_microbiome_variable_info)
dim(skin_microbiome_expression_data)

rownames(skin_microbiome_expression_data) == skin_microbiome_variable_info$variable_id
colnames(skin_microbiome_expression_data) == skin_microbiome_sample_info$sample_id

###only remain the subjects with at least >= 5
remian_subject_id =
  skin_microbiome_sample_info %>%
  dplyr::group_by(subject_id) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::ungroup() %>%
  dplyr::filter(n > 5) %>%
  dplyr::pull(subject_id)

skin_microbiome_sample_info =
  skin_microbiome_sample_info %>%
  dplyr::filter(subject_id %in% remian_subject_id)

skin_microbiome_expression_data =
  skin_microbiome_expression_data[,skin_microbiome_sample_info$sample_id]

##only remain the genus at least in 10% subjects
remain_idx =
  which(rowSums(skin_microbiome_expression_data) > 0)

skin_microbiome_expression_data = skin_microbiome_expression_data[remain_idx,]
skin_microbiome_variable_info = skin_microbiome_variable_info[remain_idx,,drop = FALSE]

remain_idx =
  skin_microbiome_expression_data %>%
  apply(1, function(x){
    sum(as.numeric(x) == 0) / ncol(skin_microbiome_expression_data)
  }) %>%
  `<`(0.9) %>%
  which()

length(remain_idx)

skin_microbiome_expression_data = skin_microbiome_expression_data[remain_idx,]
skin_microbiome_variable_info = skin_microbiome_variable_info[remain_idx,,drop = FALSE]

######--------------------------------------------------------------------------
library(plyr)

#####
dim(skin_microbiome_expression_data)

##https://rpubs.com/DKCH2020/578881
##https://ourcodingclub.github.io/tutorials/mixed-models/
##https://zhuanlan.zhihu.com/p/63092231
##https://www.linglab.cn/knowledge/10

library(compositions)
skin_microbiome_expression_data = 
  skin_microbiome_expression_data %>% 
  purrr::map(function(x){
    x = compositions::clr(x) %>% 
      as.numeric()
    x
  }) %>% 
  do.call(cbind, .) %>% 
  as.data.frame()

rownames(skin_microbiome_expression_data) = skin_microbiome_variable_info$variable_id

##step 1
###linear mixed model to adjust the subject ID random effect, and then use the partial correlation
###to get the correlation between microbiome and skin_microbiome
library(lme4)
library(rmcorr)

library(future)
library(furrr)

# intra_skin_microbiome_lm_adjusted_cor =
#   lm_adjusted_cor(
#     data_set1 = skin_microbiome_expression_data,
#     data_set2 = skin_microbiome_expression_data,
#     sample_info = skin_microbiome_sample_info,
#     method = "all",
#     threads = 8
#   )
# 
# intra_skin_microbiome_lm_adjusted_cor =
#   intra_skin_microbiome_lm_adjusted_cor[[1]]
# 
# ###remove duplicate correlation
# name =
#   intra_skin_microbiome_lm_adjusted_cor %>%
#   apply(1, function(x) {
#     paste(sort(as.character(x[c(1:2)])), collapse = "_")
#   })
# 
# intra_skin_microbiome_lm_adjusted_cor =
# intra_skin_microbiome_lm_adjusted_cor %>%
#   dplyr::mutate(name = name) %>%
#   dplyr::distinct(name, .keep_all = TRUE)
# 
# save(
#   intra_skin_microbiome_lm_adjusted_cor,
#   file = "intra_skin_microbiome_lm_adjusted_cor",
#   compress = "xz"
# )

load("intra_skin_microbiome_lm_adjusted_cor")


intra_skin_sample_wise_dim = 
  dim(skin_microbiome_expression_data)

save(intra_skin_sample_wise_dim, file = "intra_skin_sample_wise_dim")


###here we use the lm_adjusted_cor
sum(intra_skin_microbiome_lm_adjusted_cor$p_adjust < 0.2)
