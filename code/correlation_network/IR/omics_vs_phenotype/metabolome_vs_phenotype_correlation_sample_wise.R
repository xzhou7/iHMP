#' ---
#' title: "metabolome metabolome correlation"
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
dir.create("data_analysis/correlation_network/IR/metabolome_vs_phenotype_sample_wise")
setwd("data_analysis/correlation_network/IR/metabolome_vs_phenotype_sample_wise")

####load data
###metabolome
{
  load(here::here(
    "data_analysis/metabolome/data_preparation/expression_data"
  ))
  load(here::here(
    "data_analysis/metabolome/data_preparation/sample_info"
  ))
  load(here::here(
    "data_analysis/metabolome/data_preparation/variable_info"
  ))
}

metabolome_expression_data = expression_data
metabolome_sample_info = sample_info
metabolome_variable_info = variable_info

metabolome_sample_info$CollectionDate =
  as.Date(metabolome_sample_info$CollectionDate, "%m/%d/%y")

rownames(metabolome_expression_data) == metabolome_variable_info$variable_id
colnames(metabolome_expression_data) == metabolome_sample_info$sample_id

###phenotype
{
  load(here::here(
    "data_analysis/phenotype_data_sample_wise/data_preparation/expression_data"
  ))
  load(here::here("data_analysis/phenotype_data_sample_wise/data_preparation/sample_info"))
  load(here::here(
    "data_analysis/phenotype_data_sample_wise/data_preparation/variable_info"
  ))
}

phenotype_expression_data = expression_data
phenotype_sample_info = sample_info
phenotype_variable_info = variable_info

phenotype_sample_info$colletion_data =
  as.Date(phenotype_sample_info$colletion_data, "%m/%d/%y")

dim(phenotype_expression_data)
length(unique(phenotype_sample_info$subject_id))

#select IR sample first
subject_info <- 
  readr::read_csv("../../../../Figures/Figure6/metadata/metadata.subject.csv") %>% 
  dplyr::select(SubjectID, IRIS) %>% 
  dplyr::rename(subject_id = SubjectID) %>% 
  dplyr::filter(IRIS == "IR")

phenotype_sample_info <- 
  phenotype_sample_info %>% 
  dplyr::filter(subject_id %in% subject_info$subject_id)

metabolome_sample_info <- 
  metabolome_sample_info %>% 
  dplyr::filter(subject_id %in% subject_info$subject_id)

phenotype_expression_data <-
  phenotype_expression_data[,phenotype_sample_info$sample_id]  

metabolome_expression_data <-
  metabolome_expression_data[,metabolome_sample_info$sample_id]  


dim(phenotype_expression_data)
dim(metabolome_expression_data)

###match samples
dim(metabolome_sample_info)
dim(phenotype_sample_info)

length(metabolome_sample_info$subject_id)
length(unique(metabolome_sample_info$subject_id))

###just matched samples according to sample id, only 1 missed
intersect_sample_id =
  intersect(metabolome_sample_info$sample_id,
            phenotype_sample_info$sample_id)

length(intersect_sample_id)

metabolome_expression_data =
  metabolome_expression_data[,intersect_sample_id]

phenotype_expression_data =
  phenotype_expression_data[,intersect_sample_id]

metabolome_sample_info =
  metabolome_sample_info[match(intersect_sample_id, metabolome_sample_info$sample_id),]

phenotype_sample_info =
  phenotype_sample_info[match(intersect_sample_id, phenotype_sample_info$sample_id),]

length(unique(phenotype_sample_info$subject_id))

###only remain the subjects with at least >= 5
remian_subject_id =
  metabolome_sample_info %>%
  dplyr::group_by(subject_id) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::ungroup() %>%
  dplyr::filter(n > 5) %>%
  dplyr::pull(subject_id)

phenotype_sample_info =
  phenotype_sample_info %>%
  dplyr::filter(subject_id %in% remian_subject_id)

phenotype_expression_data =
  phenotype_expression_data[, phenotype_sample_info$sample_id]

metabolome_sample_info =
  metabolome_sample_info %>%
  dplyr::filter(subject_id %in% remian_subject_id)

metabolome_expression_data =
  metabolome_expression_data[,metabolome_sample_info$sample_id]

##save data
{
  save(metabolome_expression_data, file = "metabolome_expression_data")
  save(metabolome_variable_info, file = "metabolome_variable_info")
  save(metabolome_sample_info, file = "metabolome_sample_info")

  save(phenotype_expression_data, file = "phenotype_expression_data")
  save(phenotype_variable_info, file = "phenotype_variable_info")
  save(phenotype_sample_info, file = "phenotype_sample_info")
}

{
  load("metabolome_expression_data")
  load("metabolome_variable_info")
  load("metabolome_sample_info")
  
  load("phenotype_expression_data")
  load("phenotype_variable_info")
  load("phenotype_sample_info")  
}

dim(metabolome_expression_data)
dim(phenotype_expression_data)

head(metabolome_sample_info)

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

dim(metabolome_expression_data)
dim(phenotype_expression_data)
###finally, for metabolome, 302 protein, for phenotype, 52 variables

######--------------------------------------------------------------------------
##for raw data, we just log(x+1, 2)
library(plyr)

# phenotype_expression_data = 
#   log(phenotype_expression_data + 1, 2)

#####
dim(phenotype_expression_data)
dim(metabolome_expression_data)

metabolome_sample_info$subject_id == phenotype_sample_info$subject_id

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

metabolome_expression_data = 
  log(metabolome_expression_data + 1, 2)
 
metabolome_phenotype_lm_adjusted_cor <-
  lm_adjusted_cor(
    data_set1 = metabolome_expression_data,
    data_set2 = phenotype_expression_data,
    sample_info = metabolome_sample_info,
    method = "all",
    threads = 8
  )

metabolome_phenotype_lm_adjusted_cor_spearman = metabolome_phenotype_lm_adjusted_cor[[1]]

save(
  metabolome_phenotype_lm_adjusted_cor_spearman,
  file = "metabolome_phenotype_lm_adjusted_cor_spearman",
  compress = "xz"
)

load("metabolome_phenotype_lm_adjusted_cor_spearman")

###here we use the lm_adjusted_cor
sum(metabolome_phenotype_lm_adjusted_cor_spearman$p_adjust < 0.2)

cor_data =
  metabolome_phenotype_lm_adjusted_cor_spearman %>%
  dplyr::filter(p_adjust < 0.2)


       