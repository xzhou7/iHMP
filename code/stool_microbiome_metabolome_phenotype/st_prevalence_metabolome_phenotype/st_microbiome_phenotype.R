##
no_function()

setwd(masstools::get_project_wd())
library(tidyverse)
rm(list = ls())

####load data
###load diversity
load("data/from_xin/Diversity_Datatable.RData")

##prevalence
load("data/from_xin/Prevalance.RData")
st_microbiome_expression_data = ST.Pr %>% 
  t() %>% 
  as.data.frame()

###stool microbiome
load("data_analysis/st_microbiome/data_preparation/sample_info")
load("data_analysis/st_microbiome/data_preparation/variable_info")

st_microbiome_sample_info = 
  sample_info %>% 
  dplyr::distinct(subject_id, .keep_all = TRUE) %>% 
  dplyr::filter(subject_id %in% colnames(st_microbiome_expression_data))

st_microbiome_expression_data = 
  st_microbiome_expression_data[,st_microbiome_sample_info$subject_id]

st_microbiome_variable_info = 
  data.frame(variable_id = rownames(st_microbiome_expression_data))

st_microbiome_sample_info$SSPG = as.numeric(st_microbiome_sample_info$SSPG)

##only remain the genus at least in 1% subjects
remain_idx = 
  which(rowSums(st_microbiome_expression_data) > 0)

st_microbiome_expression_data = st_microbiome_expression_data[remain_idx,]
st_microbiome_variable_info = st_microbiome_variable_info[remain_idx,,drop = FALSE]


st_microbiome_expression_data %>% 
  t() %>% 
  as.data.frame() %>% 
  purrr::map(function(x){
    temp_data =
      data.frame(x = x,
                 st_microbiome_sample_info) %>%
      dplyr::filter(!is.na(SSPG)) %>%
      dplyr::mutate(SSPG = as.numeric(SSPG),
                    x = as.numeric(x))
    
    temp_data$Gender[temp_data$Gender == 'F'] = 0
    temp_data$Gender[temp_data$Gender == 'M'] = 1
    temp_data$Gender = as.numeric(temp_data$Gender)
    
    temp_data$Ethnicity[temp_data$Ethnicity == 'C'] = 0
    temp_data$Ethnicity[temp_data$Ethnicity == 'H'] = 1
    temp_data$Ethnicity[temp_data$Ethnicity == 'B'] = 2
    temp_data$Ethnicity[temp_data$Ethnicity == 'A'] = 3
    temp_data$Ethnicity = as.numeric(temp_data$Ethnicity)
    
    glm_reg =
      glm(SSPG ~ x + Gender + Adj.age + Ethnicity,
          family = gaussian,
          temp_data)
    
    temp =
      summary(glm_reg)$coefficients %>%
      as.data.frame()
    
    c(estimate = temp$Estimate[2],
      glm_p = temp$`Pr(>|t|)`[2])
    
    library(ppcor)
    
    cor_value =
      pcor.test(
        x = temp_data$x,
        y = temp_data$SSPG,
        z = temp_data[, c("Gender", "Adj.age", "Ethnicity")],
        method = "spearman"
      )
    
    result =
      c(
        estimate = unname(temp$Estimate[2]),
        glm_p = unname(temp$`Pr(>|t|)`[2]),
        cor = unname(cor_value$estimate),
        cor_p = unname(cor_value$p.value)
      )
    
    if (is.na(result[3])) {
      result[3] = 0
    }
    
    if (is.na(result[4])) {
      result[4] = 1
    }
    
    result
  })







