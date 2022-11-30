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

####HMP
physeq_HMP

# expression_data = get_taxa(physeq_HMP)
# sample_info = get_variable(physeq = physeq_HMP)

expression_data = physeq_HMP@otu_table@.Data %>% 
  t() %>% 
  as.data.frame()

variable_info = as.data.frame(physeq_HMP@tax_table@.Data)

sample_info = get_variable(physeq = physeq_HMP)

sample_info = 
sample_info %>% 
  dplyr::select(SampleID, everything()) %>% 
  dplyr::rename(sample_id = SampleID, subject_id = SubjectID) %>% 
  dplyr::left_join(sc, by = c("subject_id" = "SubjectID"))

sum(sample_info$RandomID == colnames(expression_data))

variable_info$variable_id = rownames(expression_data)

sample_info$SampleType %>% unique()

masstools::setwd_project()
setwd("data_analysis/st_metabolome/data_preparation/")

sample_info = 
sample_info %>% 
  dplyr::filter(SampleType == "ST")

dim(sample_info)

expression_data = 
  expression_data[,sample_info$RandomID]

colnames(expression_data) = sample_info$sample_id

save(sample_info, file = "sample_info")
save(expression_data, file = "expression_data")
save(variable_info, file = "variable_info")




