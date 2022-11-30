##Principal Variance Component Analysis


no_source()

# set work directory
masstools::setwd_project()
library(tidyverse)
library(phyloseq)
rm(list = ls())

source("code/tools.R")

{
  load("data_analysis/skin_microbiome/data_preparation/sample_info")
  load("data_analysis/skin_microbiome/data_preparation/expression_data")
  load("data_analysis/skin_microbiome/data_preparation/variable_info")
  
  dim(variable_info)
  
  skin_microbiome_sample_info = sample_info
  skin_microbiome_expression_data = expression_data
  skin_microbiome_variable_info = variable_info
  
  skin_microbiome_sample_info = 
    skin_microbiome_sample_info %>% 
    dplyr::mutate(SSPG = as.numeric(SSPG)) %>% 
    dplyr::mutate(iris = case_when(
      SSPG > 125 ~ "IR",
      SSPG < 125 ~ "IS"
    )) 
  
  library(lubridate)
  
  skin_microbiome_sample_info$days = 
    skin_microbiome_sample_info$Date %>% 
    lubridate::as_date() %>% 
    yday() 
  
  skin_microbiome_sample_info = 
    skin_microbiome_sample_info %>% 
    dplyr::mutate(
      days =  yday(lubridate::as_date(Date)),
      months = lubridate::month(lubridate::as_date(skin_microbiome_sample_info$Date)),
      weeks = lubridate::week(lubridate::as_date(skin_microbiome_sample_info$Date))
    ) 
}

#######work directory
masstools::setwd_project()
setwd("data_analysis/skin_microbiome/season_analysis")

zero_percent = 
skin_microbiome_expression_data %>% 
  apply(1, function(x){
    sum(x == 0)/ncol(skin_microbiome_expression_data)
  })

sum(zero_percent > 0.99)/nrow(skin_microbiome_variable_info)

##here we remove the genus with 0 > 99%
remain_idx = which(zero_percent < 0.99)

length(remain_idx)

skin_microbiome_expression_data = 
  skin_microbiome_expression_data[remain_idx,]

skin_microbiome_variable_info = 
  skin_microbiome_variable_info[remain_idx,]

##remove duplicated sample
which(duplicated(skin_microbiome_sample_info$sample_id))
skin_microbiome_sample_info$sample_id[910]
skin_microbiome_sample_info = skin_microbiome_sample_info[-910,]
skin_microbiome_expression_data = skin_microbiome_expression_data[,-910]

library(plyr)
library(lubridate)
library(mgcv)
library(gratia)
library(pvca)

temp_data <-
  apply(skin_microbiome_expression_data, 1, function(x){
    (x - mean(x))/sd(x)
  })

rownames(temp_data) == skin_microbiome_sample_info$sample_id

# ####pvca analysis
library(pvca)
assay_data <-
  as.matrix(t(temp_data))

pheno_data <-
  skin_microbiome_sample_info %>%
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

meta <- skin_microbiome_sample_info %>%
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
