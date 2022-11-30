###
no_function()

masstools::setwd_project()
library(tidyverse)
rm(list = ls())
####load raw data
load("data/from_xin/Revision_MultiOmes_0509.RData")

sample_info =
  ck.df %>%
  dplyr::select(SampleID:Plate, SubjectID:CL4)

expression_data =
  ck.df %>%
  dplyr::select(-c(SampleID:Plate, SubjectID:CL4))

rownames(expression_data) = sample_info$SampleID

variable_info = 
  data.frame(variable_id = colnames(expression_data))

dim(expression_data)
dim(variable_info)
dim(sample_info)

save(expression_data, file = here::here("data_analysis/cytokine/data_preparation/expression_data"))
save(variable_info, file = here::here("data_analysis/cytokine/data_preparation/variable_info"))
save(sample_info, file = here::here("data_analysis/cytokine/data_preparation/sample_info"))


