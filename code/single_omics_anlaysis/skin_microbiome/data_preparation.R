###
no_function()

masstools::setwd_project()
library(tidyverse)
library(phyloseq)
rm(list = ls())
####load raw data
load("data/from_xin/PhyloseqObject.RData")
load("data/from_xin/Revision_MultiOmes_0509.RData")
ls()

####UBM
physeq_UBM

expression_data = physeq_UBM@otu_table@.Data %>% 
  t() %>% 
  as.data.frame()

variable_info = as.data.frame(physeq_UBM@tax_table@.Data)

sample_info = get_variable(physeq = physeq_UBM)

sample_info = 
sample_info %>% 
  dplyr::select(SampleID, everything()) %>% 
  dplyr::rename(sample_id = SampleID, subject_id = SubjectID) %>% 
  dplyr::left_join(sc, by = c("subject_id" = "SubjectID"))

sum(sample_info$KitID == colnames(expression_data))

variable_info$variable_id = rownames(expression_data)

sample_info$SampleType %>% unique()

masstools::setwd_project()
setwd("data_analysis/skin_microbiome/data_preparation/")

sample_info = 
  sample_info %>% 
  dplyr::filter(SampleType == "Skin")

dim(sample_info)

expression_data = 
  expression_data[,sample_info$KitID]

colnames(expression_data) = sample_info$sample_id

###remove the duplicated samples
which(duplicated(sample_info$sample_id))
which(sample_info$sample_id == sample_info$sample_id[910])

sample_info[c(898,910),]
colnames(expression_data)[c(898,910)]

sample_info = 
  sample_info[-c(910),]

expression_data = 
  expression_data[,-910]

save(sample_info, file = "sample_info")
save(expression_data, file = "expression_data")
save(variable_info, file = "variable_info")



