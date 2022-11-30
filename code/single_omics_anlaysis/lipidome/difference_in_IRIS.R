###
no_function()

masstools::setwd_project()
library(tidyverse)
rm(list = ls())
library(tidyverse)
source("code/tools.R")

####load raw data
load("data_analysis/lipidome/data_preparation/new_expression_data")
load("data_analysis/lipidome/data_preparation/new_sample_info")
load("data_analysis/lipidome/data_preparation/new_variable_info")

metadata <-
  readr::read_csv("Figures/Figure6/metadata/metadata.subject.csv")

dir.create("data_analysis/lipidome/difference_in_IRIS/")
setwd("data_analysis/lipidome/difference_in_IRIS/")


####calculate the mean value

expression_data <-
  unique(new_sample_info$subject_id) %>%
  purrr::map(function(x) {
    idx <-
      which(new_sample_info$subject_id == x)
    apply(new_expression_data[, idx, drop = FALSE], 1, mean)
  }) %>%
  do.call(cbind, .) %>%
  as.data.frame()

colnames(expression_data) <- unique(new_sample_info$subject_id)


sample_info <-
  new_sample_info %>%
  dplyr::distinct(subject_id, .keep_all = TRUE)

sample_info$sample_id <-
  sample_info$subject_id

variable_info <-
  new_variable_info

sample_info <-
  sample_info %>%
  dplyr::left_join(metadata[, c("SubjectID", "IRIS", "SSPG")],
                   by = c("subject_id" = "SubjectID"))

sample_info <-
  sample_info %>%
  dplyr::filter(!is.na(IRIS)) %>%
  dplyr::filter(IRIS != "Unknown")

expression_data <-
  expression_data[, sample_info$sample_id]

rownames(expression_data)

lipid_iris_p_value <-
  seq_len(nrow(expression_data)) %>%
  purrr::map(function(i) {
    ir_idx <-
      which(sample_info$IRIS == "IR")
    is_idx <-
      which(sample_info$IRIS == "IS")
    
    x <-
      as.numeric(expression_data[i,])
    test <-
      wilcox.test(x[ir_idx], x[is_idx])
    
    data.frame(variable_id = rownames(expression_data)[i],
               p_value = test$p.value)
  }) %>%
  dplyr::bind_rows() %>%
  as.data.frame()

lipid_iris_p_value$p_value_BH <-
  p.adjust(lipid_iris_p_value$p_value, method = "BH")

library(openxlsx)

openxlsx::write.xlsx(
  lipid_iris_p_value,
  file = "lipid_iris_p_value.xlsx",
  asTable = TRUE,
  overwrite = TRUE
)

lipid_iris_p_value %>% 
  dplyr::filter(variable_id == "Lipid module_5_20")
