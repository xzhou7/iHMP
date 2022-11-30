
no_function()
# set work directory

masstools::setwd_project()
library(tidyverse)
rm(list = ls())

source("code/tools.R")

######work directory
masstools::setwd_project()
setwd("data_analysis/microbiome_stability")


####load data
###oral microbiome
{
  load(here::here("data_analysis/oral_microbiome/data_preparation/expression_data"))
  load(here::here("data_analysis/oral_microbiome/data_preparation/sample_info"))
  load(here::here("data_analysis/oral_microbiome/data_preparation/variable_info"))
}

oral_microbiome_expression_data = expression_data
oral_microbiome_sample_info = sample_info
oral_microbiome_variable_info = variable_info

###read genus table
expression_data =
  data.table::fread(here::here("data/from_xin/Genus Table/OR/Genus_OR.csv")) %>%
  as.data.frame() %>%
  tibble::column_to_rownames(var = "SampleID") %>%
  dplyr::select(-c(V1:SubjectID)) %>%
  t() %>%
  as.data.frame()

oral_microbiome_variable_info =
  oral_microbiome_variable_info[match(rownames(expression_data), oral_microbiome_variable_info$Genus),]

oral_microbiome_variable_info$Genus == rownames(expression_data)

###remove the variables which Genus are NA
remove_idx = which(is.na(oral_microbiome_variable_info$Genus))
remove_idx
if(length(remove_idx) > 0){
  oral_microbiome_variable_info = oral_microbiome_variable_info[-remove_idx,]
  expression_data = expression_data[-remove_idx,]
}

rownames(expression_data) = oral_microbiome_variable_info$variable_id

oral_microbiome_expression_data =
  expression_data

oral_microbiome_variable_info =
  oral_microbiome_variable_info %>%
  dplyr::filter(!stringr::str_detect(Genus, "Unclassified_Bacteria"))

oral_microbiome_expression_data =
  oral_microbiome_expression_data[oral_microbiome_variable_info$variable_id,]

dim(oral_microbiome_sample_info)
dim(oral_microbiome_variable_info)
dim(oral_microbiome_expression_data)

rownames(oral_microbiome_expression_data) == oral_microbiome_variable_info$variable_id
colnames(oral_microbiome_expression_data) == oral_microbiome_sample_info$sample_id

####calculate stability
oral_stability = 
purrr::map(unique(oral_microbiome_sample_info$subject_id), function(x) {
  cat(x, " ")
  sample_id = 
    oral_microbiome_sample_info %>% 
    dplyr::filter(subject_id == x) %>% 
    pull(sample_id)
  
  if(length(sample_id) < 5){
   return(NULL)
  }
  
  stability = 
  apply(oral_microbiome_expression_data[,sample_id],1,function(x){
    mean(x)/sd(x)
  })
  
  stability
})

names(oral_stability) = 
  unique(oral_microbiome_sample_info$subject_id)

oral_stability = 
oral_stability %>% 
  dplyr::bind_cols() %>% 
  as.data.frame()

rownames(oral_stability) = 
  oral_microbiome_variable_info$Genus
  

###remove the variables which are all NA in all samples
remove_idx = 
oral_stability %>% 
  apply(1, function(x){
    sum(is.na(x))/ncol(oral_stability)
  }) %>% 
  `==`(1) %>% 
  which()

remove_idx

if(length(remove_idx) > 1){
  oral_stability = 
    oral_stability[-remove_idx,]
}

save(oral_stability, file = "oral_stability")


