no_function()
# set work directory

setwd(masstools::get_project_wd())
library(tidyverse)
rm(list = ls())

source("code/tools.R")

setwd("Figures/Figure6/")

metadata_subject <- 
  readr::read_csv("clinicaldata/clinical.samples.csv")

metadata_subject <- 
  metadata_subject %>% 
  dplyr::select(SampleID:CL4) %>% 
  dplyr::rename(sample_id = SampleID,
                subject_id = SubjectID,
                colletion_data = CollectionDate) %>% 
  dplyr::select(sample_id:colletion_data, everything())

fpg <-
  readr::read_csv("clinicaldata/FPG_By_SampleID_XX.csv")

fpg <-
  fpg %>%
  dplyr::select(RandomID:FPG_Class) %>%
  dplyr::rename(subject_id = SubjectID,
                sample_id = SampleID,
                fpg_mg_ml = FPG..mg.ml.,
                fpg_class = FPG_Class) %>%
  dplyr::select(sample_id, subject_id, fpg_mg_ml, fpg_class, dplyr::everything())

sample_phenotype_data <-
  metadata_subject %>%
  dplyr::left_join(fpg, by = c("sample_id", "subject_id")) %>%
  dplyr::select(sample_id, subject_id:RandSampleID, everything())

expression_data <-
  sample_phenotype_data %>% 
  dplyr::select(A1C:WBC)

sample_info <- sample_phenotype_data %>% 
  dplyr::select(sample_id:CL4, fpg_class:RandSampleID)

variable_info <- 
  data.frame(variable_id = colnames(expression_data))

expression_data <-
  t(expression_data) %>% 
  as.data.frame()

colnames(expression_data) <- sample_info$sample_id
rownames(expression_data) <- variable_info$variable_id

setwd(masstools::get_project_wd())
save(expression_data,
     file = "data_analysis/phenotype_data_sample_wise/data_preparation/expression_data")

save(variable_info,
     file = "data_analysis/phenotype_data_sample_wise/data_preparation/variable_info")

save(sample_info,
     file = "data_analysis/phenotype_data_sample_wise/data_preparation/sample_info")

setwd(masstools::get_project_wd())
setwd("Figures/Figure6/")

subject_phenoty_data <-
  readr::read_csv("metadata/metadata.subject.csv")

subject_phenoty_data <- 
  subject_phenoty_data %>% 
  dplyr::rename(subject_id = SubjectID) %>% 
  dplyr::select(subject_id:Family)

expression_data <- 
subject_phenoty_data %>% 
  dplyr::select(FPG_Mean, SSPG, FPG, BMI, OGTT)

sample_info <-
  subject_phenoty_data %>% 
  dplyr::select(-one_of(colnames(expression_data))) %>% 
  dplyr::mutate(sample_id = subject_id)

variable_info <-
  data.frame(variable_id = colnames(expression_data))

expression_data <-
  t(expression_data) %>% 
  as.data.frame()

colnames(expression_data) <-
  sample_info$sample_id

rownames(expression_data) <-
  variable_info$variable_id

setwd(masstools::get_project_wd())
save(expression_data, 
     file = "data_analysis/phenotype_data_subject_wise/data_preparation/expression_data")
save(sample_info, 
     file = "data_analysis/phenotype_data_subject_wise/data_preparation/sample_info")
save(variable_info, 
     file = "data_analysis/phenotype_data_subject_wise/data_preparation/variable_info")


