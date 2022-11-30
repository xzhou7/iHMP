#' ---
#' title: "lipidome lipidome correlation"
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

load("data_analysis/stool_microbiome/data_preparation/sample_info")
stool_microbiome_sample_info <- 
  sample_info

######work directory
masstools::setwd_project()
dir.create("data_analysis/correlation_network/whole_data_set/lipidome_vs_phenotype_sample_wise")
setwd("data_analysis/correlation_network/whole_data_set/lipidome_vs_phenotype_sample_wise")
# 
# ####load data
# ###plasma lipidome
# {
#   load(here::here("data_analysis/lipidome/data_preparation/new_expression_data"))
#   load(here::here("data_analysis/lipidome/data_preparation/new_sample_info"))
#   load(here::here("data_analysis/lipidome/data_preparation/new_variable_info"))
# }
# 
# lipidome_expression_data = new_expression_data
# lipidome_sample_info = new_sample_info
# lipidome_variable_info = new_variable_info
# 
# ##remove QC samples
# lipidome_sample_info =
#   lipidome_sample_info %>%
#   dplyr::filter(subject_id != "QC")
# 
# lipidome_expression_data =
#   lipidome_expression_data[,lipidome_sample_info$sample_id]
# 
# dim(lipidome_expression_data)
# length(unique(lipidome_sample_info$subject_id))
# #
# ###phenotype
# {
#   load(here::here(
#     "data_analysis/phenotype_data_sample_wise/data_preparation/expression_data"
#   ))
#   load(here::here("data_analysis/phenotype_data_sample_wise/data_preparation/sample_info"))
#   load(here::here(
#     "data_analysis/phenotype_data_sample_wise/data_preparation/variable_info"
#   ))
# }
# 
# phenotype_expression_data = expression_data
# phenotype_sample_info = sample_info
# phenotype_variable_info = variable_info
# 
# phenotype_sample_info$colletion_data =
#   as.Date(phenotype_sample_info$colletion_data, "%m/%d/%y")
# 
# dim(phenotype_expression_data)
# length(unique(phenotype_sample_info$subject_id))
# 
# ###match samples
# dim(lipidome_sample_info)
# dim(phenotype_sample_info)
# 
# length(lipidome_sample_info$subject_id)
# length(unique(lipidome_sample_info$subject_id))
# 
# ###just matched samples according to sample id, only 1 missed
# intersect_sample_id =
#   intersect(lipidome_sample_info$sample_id,
#             phenotype_sample_info$sample_id)
# 
# length(intersect_sample_id)
# 
# lipidome_expression_data =
#   lipidome_expression_data[,intersect_sample_id]
# 
# phenotype_expression_data =
#   phenotype_expression_data[,intersect_sample_id]
# 
# lipidome_sample_info =
#   lipidome_sample_info[match(intersect_sample_id, lipidome_sample_info$sample_id),]
# 
# phenotype_sample_info =
#   phenotype_sample_info[match(intersect_sample_id, phenotype_sample_info$sample_id),]
# 
# length(unique(phenotype_sample_info$subject_id))
# 
# ###only remain the subjects with at least >= 5
# remian_subject_id =
#   lipidome_sample_info %>%
#   dplyr::group_by(subject_id) %>%
#   dplyr::summarise(n = n()) %>%
#   dplyr::ungroup() %>%
#   dplyr::filter(n > 5) %>%
#   dplyr::pull(subject_id)
# 
# phenotype_sample_info =
#   phenotype_sample_info %>%
#   dplyr::filter(subject_id %in% remian_subject_id)
# 
# phenotype_expression_data =
#   phenotype_expression_data[, phenotype_sample_info$sample_id]
# 
# lipidome_sample_info =
#   lipidome_sample_info %>%
#   dplyr::filter(subject_id %in% remian_subject_id)
# 
# lipidome_expression_data =
#   lipidome_expression_data[,lipidome_sample_info$sample_id]
# 
# ##save data
# {
#   save(lipidome_expression_data, file = "lipidome_expression_data")
#   save(lipidome_variable_info, file = "lipidome_variable_info")
#   save(lipidome_sample_info, file = "lipidome_sample_info")
# 
#   save(phenotype_expression_data, file = "phenotype_expression_data")
#   save(phenotype_variable_info, file = "phenotype_variable_info")
#   save(phenotype_sample_info, file = "phenotype_sample_info")
# }

{
  load("lipidome_expression_data")
  load("lipidome_variable_info")
  load("lipidome_sample_info")
  
  load("phenotype_expression_data")
  load("phenotype_variable_info")
  load("phenotype_sample_info")  
}

dim(lipidome_expression_data)
dim(phenotype_expression_data)

head(lipidome_sample_info)

#####remove phenotype who have a lot of zeros
na_percentage <-
  phenotype_expression_data %>% 
  apply(1, function(x){
    sum(is.na(x))/ncol(phenotype_expression_data)
  })

remain_idx <- unname(which(na_percentage < 0.5))

phenotype_expression_data <-
  phenotype_expression_data[remain_idx,]

phenotype_variable_info <-
  phenotype_variable_info[remain_idx,,drop = FALSE]

dim(lipidome_expression_data)
dim(phenotype_expression_data)
###finally, for lipidome, 302 protein, for phenotype, 52 variables

######--------------------------------------------------------------------------
##for raw data, we just log(x+1, 2)
library(plyr)

# phenotype_expression_data = 
#   log(phenotype_expression_data + 1, 2)

#####
dim(phenotype_expression_data)
dim(lipidome_expression_data)

lipidome_sample_info$subject_id == phenotype_sample_info$subject_id

##https://rpubs.com/DKCH2020/578881
##https://ourcodingclub.github.io/tutorials/mixed-models/
##https://zhuanlan.zhihu.com/p/63092231
##https://www.linglab.cn/knowledge/10

###step 1
###linear mixed model to adjust the subject ID random effect, and then use the partial correlation
###to get the correlation between microbiome and phenotype
library(lme4)
library(rmcorr)

library(future)
library(furrr)

dim(phenotype_sample_info)

match(unique(lipidome_sample_info$subject_id), unique(stool_microbiome_sample_info$subject_id))

lipidome_sample_info <- 
lipidome_sample_info %>% 
  dplyr::left_join(stool_microbiome_sample_info[,c("subject_id", "Gender", "Adj.age", "Ethnicity")] %>% 
                     dplyr::distinct(subject_id, .keep_all = TRUE),
                   by = "subject_id") %>% 
  dplyr::filter(!is.na(Gender))

lipidome_expression_data <-
  lipidome_expression_data[,lipidome_sample_info$sample_id]

phenotype_expression_data <-
  phenotype_expression_data[,lipidome_sample_info$sample_id]
# 
# lipidome_phenotype_lm_adjusted_cor <-
#   lm_adjusted_cor(
#     data_set1 = lipidome_expression_data,
#     data_set2 = phenotype_expression_data,
#     sample_info = lipidome_sample_info,
#     method = "all",
#     threads = 8
#   )
# 
# lipidome_phenotype_lm_adjusted_cor_spearman = lipidome_phenotype_lm_adjusted_cor[[1]]
# 
# save(
#   lipidome_phenotype_lm_adjusted_cor_spearman,
#   file = "lipidome_phenotype_lm_adjusted_cor_spearman",
#   compress = "xz"
# )

load("lipidome_phenotype_lm_adjusted_cor_spearman")

###here we use the lm_adjusted_cor
sum(lipidome_phenotype_lm_adjusted_cor_spearman$p_adjust < 0.2)

cor_data =
  lipidome_phenotype_lm_adjusted_cor_spearman %>%
  dplyr::filter(p_adjust < 0.2)

#

