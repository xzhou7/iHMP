
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

####calculate stability
nasal_stability = 
purrr::map(unique(nasal_microbiome_sample_info$subject_id), function(x) {
  cat(x, " ")
  sample_id = 
    nasal_microbiome_sample_info %>% 
    dplyr::filter(subject_id == x) %>% 
    pull(sample_id)
  
  if(length(sample_id) < 5){
   return(NULL)
  }
  
  stability = 
  apply(nasal_microbiome_expression_data[,sample_id],1,function(x){
    mean(x)/sd(x)
  })
  
  stability
})

names(nasal_stability) = 
  unique(nasal_microbiome_sample_info$subject_id)

nasal_stability = 
nasal_stability %>% 
  dplyr::bind_cols() %>% 
  as.data.frame()

rownames(nasal_stability) = 
  nasal_microbiome_variable_info$Genus

###remove the variables which are all NA in all samples
remove_idx = 
nasal_stability %>% 
  apply(1, function(x){
    sum(is.na(x))/ncol(nasal_stability)
  }) %>% 
  `==`(1) %>% 
  which()

remove_idx

if(length(remove_idx) > 1){
  nasal_stability = 
    nasal_stability[-remove_idx,]
}

save(nasal_stability, file = "nasal_stability")


