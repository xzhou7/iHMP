##Principal Variance Component Analysis

no_source()

# set work directory
masstools::setwd_project()
library(tidyverse)
library(phyloseq)
rm(lioral = ls())

source("code/tools.R")

{
  load("data_analysis/oral_microbiome/data_preparation/sample_info")
  load("data_analysis/oral_microbiome/data_preparation/expression_data")
  load("data_analysis/oral_microbiome/data_preparation/variable_info")
  
  dim(variable_info)
  
  oral_microbiome_sample_info = sample_info
  oral_microbiome_expression_data = expression_data
  oral_microbiome_variable_info = variable_info
  
  oral_microbiome_sample_info = 
    oral_microbiome_sample_info %>% 
    dplyr::mutate(SSPG = as.numeric(SSPG)) %>% 
    dplyr::mutate(iris = case_when(
      SSPG > 125 ~ "IR",
      SSPG < 125 ~ "IS"
    )) 
  
  library(lubridate)
  
  oral_microbiome_sample_info$days = 
    oral_microbiome_sample_info$Date %>% 
    lubridate::as_date() %>% 
    yday() 
  
  oral_microbiome_sample_info = 
    oral_microbiome_sample_info %>% 
    dplyr::mutate(
      days =  yday(lubridate::as_date(Date)),
      months = lubridate::month(lubridate::as_date(oral_microbiome_sample_info$Date)),
      weeks = lubridate::week(lubridate::as_date(oral_microbiome_sample_info$Date))
    ) 
}

#######work directory
masstools::setwd_project()
setwd("data_analysis/oral_microbiome/season_analysis")

zero_percent = 
oral_microbiome_expression_data %>% 
  apply(1, function(x){
    sum(x == 0)/ncol(oral_microbiome_expression_data)
  })

sum(zero_percent > 0.99)/nrow(oral_microbiome_variable_info)

##here we remove the genus with 0 > 99%
remain_idx = which(zero_percent < 0.99)

length(remain_idx)

oral_microbiome_expression_data = 
  oral_microbiome_expression_data[remain_idx,]

oral_microbiome_variable_info = 
  oral_microbiome_variable_info[remain_idx,]

##remove duplicated sample
which(duplicated(oral_microbiome_sample_info$sample_id))
# oral_microbiome_sample_info$sample_id[910]
# oral_microbiome_sample_info = oral_microbiome_sample_info[-910,]
# oral_microbiome_expression_data = oral_microbiome_expression_data[,-910]

library(plyr)
library(lubridate)
library(mgcv)
library(gratia)
library(pvca)

temp_data <-
  apply(oral_microbiome_expression_data, 1, function(x){
    (x - mean(x))/sd(x)
  })

rownames(temp_data) == oral_microbiome_sample_info$sample_id

# ####pvca analysis
library(pvca)
assay_data <-
  as.matrix(t(temp_data))

pheno_data <-
  oral_microbiome_sample_info %>%
  dplyr::select(
    subject_id,
    sample_id,
    days
  ) %>%
  as.data.frame()

library(Biobase)

pheno_data <- AnnotatedDataFrame(data = pheno_data)

row.names(pheno_data) <- colnames(assay_data)

data_pvca <-
  Biobase::ExpressionSet(
    assayData = assay_data,
    phenoData = pheno_data
  )

pct_threshold <- 0.6

batch.factors <-
  c(
    "subject_id",#yes
    "days"
  )

pvcaObj <- pvcaBatchAssess(abatch = data_pvca,
                           batch.factors = batch.factors,
                           threshold = pct_threshold)

save(pvcaObj, file = "pvcaObj")
load("pvcaObj")

plot <-
  plot_pvca(object = pvcaObj)

plot

ggsave(plot, filename = "pvca_plot.pdf", height = 9, width = 7)

counts <- as.matrix(t(temp_data))

meta <- oral_microbiome_sample_info %>%
  dplyr::select(sample_id,
                subject_id,
                days) %>%
  as.data.frame() 
  
rownames(meta) = NULL

meta = 
meta %>% 
tibble::column_to_rownames(var = "sample_id")

pvca_object <-
  PVCA(counts = counts, meta = meta, threshold = 0.6, inter = FALSE)

pvca_object
save(pvca_object, file = "pvca_object")
load("pvca_object")
plot <- 
  PlotPVCA(pvca.res = pvca_object, title = "")
plot
ggsave(plot, filename = "pvca_plot.pdf", height = 7, width = 7)
