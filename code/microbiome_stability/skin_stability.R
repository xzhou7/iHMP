
no_function()
# set work directory

setwd(masstools::get_project_wd())
library(tidyverse)
rm(list = ls())

source("code/tools.R")

######work directory
setwd(masstools::get_project_wd())
setwd("data_analysis/microbiome_stability")

###skin microbiome
{
  load(here::here("data_analysis/skin_microbiome/data_preparation/expression_data"))
  load(here::here("data_analysis/skin_microbiome/data_preparation/sample_info"))
  load(here::here("data_analysis/skin_microbiome/data_preparation/variable_info"))
}

skin_microbiome_expression_data = expression_data
skin_microbiome_sample_info = sample_info
skin_microbiome_variable_info = variable_info

###read genus table
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

####calculate stability
skin_stability = 
purrr::map(unique(skin_microbiome_sample_info$subject_id), function(x) {
  cat(x, " ")
  sample_id = 
    skin_microbiome_sample_info %>% 
    dplyr::filter(subject_id == x) %>% 
    pull(sample_id)
  
  if(length(sample_id) < 5){
   return(NULL)
  }
  
  stability = 
  apply(skin_microbiome_expression_data[,sample_id],1,function(x){
    mean(x)/sd(x)
  })
  
  stability
})

names(skin_stability) = 
  unique(skin_microbiome_sample_info$subject_id)

skin_stability = 
skin_stability %>% 
  dplyr::bind_cols() %>% 
  as.data.frame()

rownames(skin_stability) = 
  skin_microbiome_variable_info$Genus
  

###remove the variables which are all NA in all samples
remove_idx = 
skin_stability %>% 
  apply(1, function(x){
    sum(is.na(x))/ncol(skin_stability)
  }) %>% 
  `==`(1) %>% 
  which()

remove_idx

if(length(remove_idx) > 1){
  skin_stability = 
    skin_stability[-remove_idx,]
}

save(skin_stability, file = "skin_stability")


