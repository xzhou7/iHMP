#' ---
#' title: "cytokine cytokine correlation"
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
dir.create("data_analysis/correlation_network/IS/cytokine_vs_phenotype_sample_wise")
setwd("data_analysis/correlation_network/IS/cytokine_vs_phenotype_sample_wise")

####load data
###cytokine
{
  load(here::here(
    "data_analysis/cytokine/data_preparation/expression_data"
  ))
  load(here::here(
    "data_analysis/cytokine/data_preparation/sample_info"
  ))
  load(here::here(
    "data_analysis/cytokine/data_preparation/variable_info"
  ))
}

cytokine_expression_data = t(expression_data) %>%
  as.data.frame()
cytokine_sample_info = sample_info
cytokine_variable_info = variable_info

cytokine_sample_info$CollectionDate =
  as.Date(cytokine_sample_info$CollectionDate, "%m/%d/%y")

cytokine_sample_info <-
cytokine_sample_info %>%
  dplyr::rename(sample_id = SampleID,
                subject_id = SubjectID)

rownames(cytokine_expression_data) == cytokine_variable_info$variable_id
colnames(cytokine_expression_data) == cytokine_sample_info$sample_id

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

#select IS sample first
subject_info <- 
  readr::read_csv("../../../../Figures/Figure6/metadata/metadata.subject.csv") %>% 
  dplyr::select(SubjectID, IRIS) %>% 
  dplyr::rename(subject_id = SubjectID) %>% 
  dplyr::filter(IRIS == "IS")

phenotype_sample_info <- 
  phenotype_sample_info %>% 
  dplyr::filter(subject_id %in% subject_info$subject_id)

cytokine_sample_info <- 
  cytokine_sample_info %>% 
  dplyr::filter(subject_id %in% subject_info$subject_id)

phenotype_expression_data <-
  phenotype_expression_data[,phenotype_sample_info$sample_id]  

cytokine_expression_data <-
  cytokine_expression_data[,cytokine_sample_info$sample_id]  

cytokine_expression_data <- filter(cytokine_expression_data, 
                                  ! (rownames(cytokine_expression_data) %in% c("CHEX1", "CHEX2", "CHEX3", "CHEX4")))

dim(phenotype_expression_data)
dim(cytokine_expression_data)

###match samples
dim(cytokine_sample_info)
dim(phenotype_sample_info)

length(cytokine_sample_info$subject_id)
length(unique(cytokine_sample_info$subject_id))

###just matched samples according to sample id, only 1 missed
intersect_sample_id =
  intersect(cytokine_sample_info$sample_id,
            phenotype_sample_info$sample_id)

length(intersect_sample_id)

cytokine_expression_data =
  cytokine_expression_data[,intersect_sample_id]

phenotype_expression_data =
  phenotype_expression_data[,intersect_sample_id]

cytokine_sample_info =
  cytokine_sample_info[match(intersect_sample_id, cytokine_sample_info$sample_id),]

phenotype_sample_info =
  phenotype_sample_info[match(intersect_sample_id, phenotype_sample_info$sample_id),]

length(unique(phenotype_sample_info$subject_id))

###only remain the subjects with at least >= 5
remian_subject_id =
  cytokine_sample_info %>%
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

cytokine_sample_info =
  cytokine_sample_info %>%
  dplyr::filter(subject_id %in% remian_subject_id)

cytokine_expression_data =
  cytokine_expression_data[,cytokine_sample_info$sample_id]

##save data
{
  save(cytokine_expression_data, file = "cytokine_expression_data")
  save(cytokine_variable_info, file = "cytokine_variable_info")
  save(cytokine_sample_info, file = "cytokine_sample_info")

  save(phenotype_expression_data, file = "phenotype_expression_data")
  save(phenotype_variable_info, file = "phenotype_variable_info")
  save(phenotype_sample_info, file = "phenotype_sample_info")
}

{
  load("cytokine_expression_data")
  load("cytokine_variable_info")
  load("cytokine_sample_info")
  
  load("phenotype_expression_data")
  load("phenotype_variable_info")
  load("phenotype_sample_info")  
}

dim(cytokine_expression_data)
dim(phenotype_expression_data)

head(cytokine_sample_info)

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

dim(cytokine_expression_data)
dim(phenotype_expression_data)
###finally, for cytokine, 302 protein, for phenotype, 52 variables

######--------------------------------------------------------------------------
##for raw data, we just log(x+1, 2)
library(plyr)

plot(density(as.numeric(cytokine_expression_data[1,])))

cytokine_expression_data =
  log(cytokine_expression_data + 1, 2)

#####
dim(phenotype_expression_data)
dim(cytokine_expression_data)

cytokine_sample_info$subject_id == phenotype_sample_info$subject_id

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
 
cytokine_sample_info <- 
cytokine_sample_info %>% 
  dplyr::left_join(stool_microbiome_sample_info[,c("subject_id", "Gender", "Ethnicity", "Adj.age", "BMI")] %>% 
                     dplyr::distinct(subject_id, .keep_all = TRUE), by = "subject_id")

cytokine_phenotype_lm_adjusted_cor <-
  lm_adjusted_cor(
    data_set1 = cytokine_expression_data,
    data_set2 = phenotype_expression_data,
    sample_info = cytokine_sample_info,
    method = "all",
    threads = 8
  )

cytokine_phenotype_lm_adjusted_cor_spearman = cytokine_phenotype_lm_adjusted_cor[[1]]

save(
  cytokine_phenotype_lm_adjusted_cor_spearman,
  file = "cytokine_phenotype_lm_adjusted_cor_spearman",
  compress = "xz"
)

load("cytokine_phenotype_lm_adjusted_cor_spearman")

###here we use the lm_adjusted_cor
sum(cytokine_phenotype_lm_adjusted_cor_spearman$p_adjust < 0.2)

cor_data =
  cytokine_phenotype_lm_adjusted_cor_spearman %>%
  dplyr::filter(p_adjust < 0.2)

#

