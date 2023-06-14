
no_function()
# set work directory

setwd(masstools::get_project_wd())
library(tidyverse)
rm(list = ls())

source("code/tools.R")

######work directory
setwd(masstools::get_project_wd())
setwd("data_analysis/microbiome_stability")

####load data
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

rownames(stool_microbiome_expression_data) == stool_microbiome_variable_info$variable_id
colnames(stool_microbiome_expression_data) == stool_microbiome_sample_info$sample_id

####calculate stability
stool_stability = 
purrr::map(unique(stool_microbiome_sample_info$subject_id), function(x) {
  cat(x, " ")
  sample_id = 
    stool_microbiome_sample_info %>% 
    dplyr::filter(subject_id == x) %>% 
    pull(sample_id)
  
  if(length(sample_id) < 5){
   return(NULL)
  }
  
  stability = 
  apply(stool_microbiome_expression_data[,sample_id],1,function(x){
    mean(x)/sd(x)
  })
  
  stability
})

names(stool_stability) = 
  unique(stool_microbiome_sample_info$subject_id)

stool_stability = 
stool_stability %>% 
  dplyr::bind_cols() %>% 
  as.data.frame()

rownames(stool_stability) = 
  stool_microbiome_variable_info$Genus

###remove the variables which are all NA in all samples
remove_idx = 
stool_stability %>% 
  apply(1, function(x){
    sum(is.na(x))/ncol(stool_stability)
  }) %>% 
  `==`(1) %>% 
  which()

if(length(remove_idx) > 1){
  stool_stability = 
    stool_stability[-remove_idx,]
}

save(stool_stability, file = "stool_stability")


