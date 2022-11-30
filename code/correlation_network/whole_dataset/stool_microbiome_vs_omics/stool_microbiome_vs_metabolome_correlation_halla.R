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

masstools::setwd_project()
library(tidyverse)
rm(list = ls())

source("code/tools.R")

######work directory
masstools::setwd_project()
setwd("data_analysis/stool_microbiome_vs_metabolome")

####load data
{
  load("stool_microbiome_expression_data")
  load("stool_microbiome_variable_info")
  load("stool_microbiome_sample_info")
  
  load("metabolome_expression_data")
  load("metabolome_variable_info")
  load("metabolome_sample_info")  
}

dim(stool_microbiome_expression_data)
dim(metabolome_expression_data)

###finally, for stool microbiome, 107 genus, for metabolome, 169 metabolites

######--------------------------------------------------------------------------
##for raw data, we just log(x+1, 2)
library(plyr)

metabolome_expression_data = 
  log(metabolome_expression_data + 1, 2)

library(plyr)

###because our microbiome are percentage data, so here we use the CTL method
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

#####
dim(metabolome_expression_data)
dim(stool_microbiome_expression_data)

stool_microbiome_sample_info$subject_id == metabolome_sample_info$subject_id

#####here we normalize them first
dim(stool_microbiome_expression_data)

data1 = 
  stool_microbiome_expression_data %>% 
  t() %>% 
  as.data.frame() %>% 
  purrr::map(function(x){
    temp_data = 
      data.frame(stool_microbiome_sample_info, x)
    temp_data$Gender[temp_data$Gender == 'F'] = 0
    temp_data$Gender[temp_data$Gender == 'M'] = 1
    temp_data$Gender = as.numeric(temp_data$Gender)
    
    temp_data$Ethnicity[temp_data$Ethnicity == 'C'] = 0
    temp_data$Ethnicity[temp_data$Ethnicity == 'H'] = 1
    temp_data$Ethnicity[temp_data$Ethnicity == 'B'] = 2
    temp_data$Ethnicity[temp_data$Ethnicity == 'A'] = 3
    temp_data$Ethnicity = as.numeric(temp_data$Ethnicity)
    
    temp_data$SSPG = as.numeric(temp_data$SSPG)
    temp_data$FPG = as.numeric(temp_data$FPG)
    adjusted_x =
      lme4::lmer(formula = x ~ Gender + Adj.age + Ethnicity + (1 |
                                                                 subject_id),
                 data = temp_data) %>%
      residuals()
    adjusted_x
  }) %>% 
  do.call(rbind, .) %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "variable_id")

data2 = 
  metabolome_expression_data %>% 
  t() %>% 
  as.data.frame() %>% 
  purrr::map(function(x){
    temp_data = 
      data.frame(stool_microbiome_sample_info, x)
    temp_data$Gender[temp_data$Gender == 'F'] = 0
    temp_data$Gender[temp_data$Gender == 'M'] = 1
    temp_data$Gender = as.numeric(temp_data$Gender)
    
    temp_data$Ethnicity[temp_data$Ethnicity == 'C'] = 0
    temp_data$Ethnicity[temp_data$Ethnicity == 'H'] = 1
    temp_data$Ethnicity[temp_data$Ethnicity == 'B'] = 2
    temp_data$Ethnicity[temp_data$Ethnicity == 'A'] = 3
    temp_data$Ethnicity = as.numeric(temp_data$Ethnicity)
    
    temp_data$SSPG = as.numeric(temp_data$SSPG)
    temp_data$FPG = as.numeric(temp_data$FPG)
    adjusted_x =
      lme4::lmer(formula = x ~ Gender + Adj.age + Ethnicity + (1 |
                                                                 subject_id),
                 data = temp_data) %>%
      residuals()
    adjusted_x
  }) %>% 
  do.call(rbind, .) %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "variable_id")

write.table(
  data1,
  "data1.txt",
  quote = FALSE,
  sep = "\t",
  col.names = TRUE,
  row.names = FALSE
)

write.table(
  data2,
  "data2.txt",
  quote = FALSE,
  sep = "\t",
  col.names = TRUE,
  row.names = FALSE
)

###next code should be run in terminal
# cd /Users/xiaotaoshen/Box/Xiaotao Shen's Files/human_microbiome_project/data_analysis/stool_microbiome_vs_metabolome
# halla -x data1.txt -y data2.txt -m spearman -o halla_result

associations =
  read.table("halla_result/all_associations.txt", header = TRUE)

load("stool_microbiome_metabolome_lm_adjusted_cor_spearman")
head(stool_microbiome_metabolome_lm_adjusted_cor_spearman)

dim(stool_microbiome_metabolome_lm_adjusted_cor_spearman)
dim(associations)

temp = 
stool_microbiome_metabolome_lm_adjusted_cor_spearman[, c("microbiome", "metabolite", "cor", "p", "p_adjust2")] %>%
  dplyr::left_join(associations,
                   by = c("microbiome" = "X_features", "metabolite" = "Y_features"))


plot(temp$cor, temp$association)
abline(0,1)

plot(temp$p, temp$p.values)
abline(0,1)

plot(temp$p_adjust2, temp$q.values)
abline(0,1)

sum(temp$p_adjust2 < 0.05)
sum(temp$q.values < 0.05)



