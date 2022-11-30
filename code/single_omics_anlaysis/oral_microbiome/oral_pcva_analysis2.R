##Principal Variance Component Analysis
no_source()

# set work directory
masstools::setwd_project()
library(tidyverse)
library(phyloseq)
rm(lins = ls())

source("code/tools.R")

####load sample information
load("data_analysis/oral_microbiome/data_preparation/sample_info")
dim(variable_info)

oral_microbiome_sample_info = sample_info

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

######ASV level

{
  oral_microbiome_expression_data = 
    data.table::fread(here::here("data/from_xin/Genus Table/OR/ASV_OR.csv")) %>% 
    dplyr::select(SampleID, everything()) %>% 
    dplyr::select(-c(V1:SubjectID)) %>% 
    as.data.frame() %>% 
    tibble::column_to_rownames(var = "SampleID") %>% 
    t() %>% 
    as.data.frame()
  
  oral_microbiome_expression_data = oral_microbiome_expression_data[,oral_microbiome_sample_info$sample_id]
}

#######work directory
masstools::setwd_project()
setwd("data_analysis/oral_microbiome/PVCA_analysis/")

zero_percent =
  oral_microbiome_expression_data %>%
  apply(1, function(x) {
    sum(x == 0) / ncol(oral_microbiome_expression_data)
  })

sum(zero_percent > 0.99)/nrow(oral_microbiome_expression_data)

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

# # ####pvca analysis
# library(pvca)
# assay_data <-
#   as.matrix(t(temp_data))
# 
# pheno_data <-
#   oral_microbiome_sample_info %>%
#   dplyr::select(
#     subject_id,
#     sample_id,
#     days
#   ) %>%
#   as.data.frame()
# 
# library(Biobase)
# 
# pheno_data <- AnnotatedDataFrame(data = pheno_data)
# 
# row.names(pheno_data) <- colnames(assay_data)
# 
# data_pvca <-
#   Biobase::ExpressionSet(
#     assayData = assay_data,
#     phenoData = pheno_data
#   )
# 
# pct_threshold <- 0.6
# 
# batch.factors <-
#   c(
#     "subject_id",#yes
#     "days"
#   )
# 
# pvcaObj_asv <- pvcaBatchAssess(abatch = data_pvca,
#                            batch.factors = batch.factors,
#                            threshold = pct_threshold)
# 
# save(pvcaObj_asv, file = "pvcaObj_asv")
# load("pvcaObj_asv")
# 
# plot <-
#   plot_pvca(object = pvcaObj_asv)
# 
# plot
# 
# ggsave(plot, filename = "pvca_plot.pdf", height = 9, width = 7)

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

pvca_object_asv <-
  PVCA(counts = counts, meta = meta, threshold = 0.6, inter = FALSE)

pvca_object_asv
save(pvca_object_asv, file = "pvca_object_asv")
load("pvca_object_asv")
plot <- 
  PlotPVCA(pvca.res = pvca_object_asv, title = "")
plot
ggsave(plot, filename = "pvca_plot_asv.pdf", height = 7, width = 7)





######phylum level
{
  oral_microbiome_expression_data = 
    data.table::fread(here::here("data/from_xin/Genus Table/OR/Phylum_OR.csv")) %>% 
    dplyr::select(SampleID, everything()) %>% 
    dplyr::select(-c(V1:SubjectID)) %>% 
    as.data.frame() %>% 
    tibble::column_to_rownames(var = "SampleID") %>% 
    t() %>% 
    as.data.frame()
  
  oral_microbiome_expression_data = oral_microbiome_expression_data[,oral_microbiome_sample_info$sample_id]
}

#######work directory
masstools::setwd_project()
setwd("data_analysis/oral_microbiome/PVCA_analysis/")

zero_percent =
  oral_microbiome_expression_data %>%
  apply(1, function(x) {
    sum(x == 0) / ncol(oral_microbiome_expression_data)
  })

sum(zero_percent > 0.99)/nrow(oral_microbiome_expression_data)

##here we remove the genus with 0 > 99%
remain_idx = which(zero_percent < 0.99)

length(remain_idx)

oral_microbiome_expression_data = 
  oral_microbiome_expression_data[remain_idx,]

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

pvca_object_phylum <-
  PVCA(counts = counts, meta = meta, threshold = 0.6, inter = FALSE)

pvca_object_phylum
save(pvca_object_phylum, file = "pvca_object_phylum")
load("pvca_object_phylum")
plot <- 
  PlotPVCA(pvca.res = pvca_object_phylum, title = "")
plot
ggsave(plot, filename = "pvca_plot_phylum.pdf", height = 7, width = 7)



######class level
{
  oral_microbiome_expression_data = 
    data.table::fread(here::here("data/from_xin/Genus Table/OR/Class_OR.csv")) %>% 
    dplyr::select(SampleID, everything()) %>% 
    dplyr::select(-c(V1:SubjectID)) %>% 
    as.data.frame() %>% 
    tibble::column_to_rownames(var = "SampleID") %>% 
    t() %>% 
    as.data.frame()
  
  oral_microbiome_expression_data = oral_microbiome_expression_data[,oral_microbiome_sample_info$sample_id]
  
}

#######work directory
masstools::setwd_project()
setwd("data_analysis/oral_microbiome/PVCA_analysis/")

zero_percent =
  oral_microbiome_expression_data %>%
  apply(1, function(x) {
    sum(x == 0) / ncol(oral_microbiome_expression_data)
  })

sum(zero_percent > 0.99)/nrow(oral_microbiome_expression_data)

##here we remove the genus with 0 > 99%
remain_idx = which(zero_percent < 0.99)

length(remain_idx)

oral_microbiome_expression_data = 
  oral_microbiome_expression_data[remain_idx,]

##remove duplicated sample
which(duplicated(oral_microbiome_sample_info$sample_id))
# oral_microbiome_sample_info$sample_id[910]
# oral_microbiome_sample_info = oral_microbiome_sample_info[-910,]
# oral_microbiome_expression_data = oral_microbiome_expression_data[,-910]

temp_data <-
  apply(oral_microbiome_expression_data, 1, function(x){
    (x - mean(x))/sd(x)
  })

rownames(temp_data) == oral_microbiome_sample_info$sample_id

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

pvca_object_class <-
  PVCA(counts = counts, meta = meta, threshold = 0.6, inter = FALSE)

pvca_object_class
save(pvca_object_class, file = "pvca_object_class")
load("pvca_object_class")
plot <- 
  PlotPVCA(pvca.res = pvca_object_class, title = "")
plot
ggsave(plot, filename = "pvca_plot_class.pdf", height = 7, width = 7)




######family level
{
  oral_microbiome_expression_data = 
    data.table::fread(here::here("data/from_xin/Genus Table/OR/Family_OR.csv")) %>% 
    dplyr::select(SampleID, everything()) %>% 
    dplyr::select(-c(V1:SubjectID)) %>% 
    as.data.frame() %>% 
    tibble::column_to_rownames(var = "SampleID") %>% 
    t() %>% 
    as.data.frame()
  
  oral_microbiome_expression_data = oral_microbiome_expression_data[,oral_microbiome_sample_info$sample_id]
  
}

#######work directory
masstools::setwd_project()
setwd("data_analysis/oral_microbiome/PVCA_analysis/")

zero_percent =
  oral_microbiome_expression_data %>%
  apply(1, function(x) {
    sum(x == 0) / ncol(oral_microbiome_expression_data)
  })

sum(zero_percent > 0.99)/nrow(oral_microbiome_expression_data)

##here we remove the genus with 0 > 99%
remain_idx = which(zero_percent < 0.99)

length(remain_idx)

oral_microbiome_expression_data = 
  oral_microbiome_expression_data[remain_idx,]

##remove duplicated sample
which(duplicated(oral_microbiome_sample_info$sample_id))
# oral_microbiome_sample_info$sample_id[910]
# oral_microbiome_sample_info = oral_microbiome_sample_info[-910,]
# oral_microbiome_expression_data = oral_microbiome_expression_data[,-910]

temp_data <-
  apply(oral_microbiome_expression_data, 1, function(x){
    (x - mean(x))/sd(x)
  })

rownames(temp_data) == oral_microbiome_sample_info$sample_id

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

pvca_object_family <-
  PVCA(counts = counts, meta = meta, threshold = 0.6, inter = FALSE)

pvca_object_family
save(pvca_object_family, file = "pvca_object_family")
load("pvca_object_family")
plot <- 
  PlotPVCA(pvca.res = pvca_object_family, title = "")
plot
ggsave(plot, filename = "pvca_plot_family.pdf", height = 7, width = 7)




######genus level
{
  oral_microbiome_expression_data = 
    data.table::fread(here::here("data/from_xin/Genus Table/OR/Genus_OR.csv")) %>% 
    dplyr::select(SampleID, everything()) %>% 
    dplyr::select(-c(V1:SubjectID)) %>% 
    as.data.frame() %>% 
    tibble::column_to_rownames(var = "SampleID") %>% 
    t() %>% 
    as.data.frame()
  
  oral_microbiome_expression_data = oral_microbiome_expression_data[,oral_microbiome_sample_info$sample_id]
  
}

#######work directory
masstools::setwd_project()
setwd("data_analysis/oral_microbiome/PVCA_analysis/")

zero_percent =
  oral_microbiome_expression_data %>%
  apply(1, function(x) {
    sum(x == 0) / ncol(oral_microbiome_expression_data)
  })

sum(zero_percent > 0.99)/nrow(oral_microbiome_expression_data)

##here we remove the genus with 0 > 99%
remain_idx = which(zero_percent < 0.99)

length(remain_idx)

oral_microbiome_expression_data = 
  oral_microbiome_expression_data[remain_idx,]

##remove duplicated sample
which(duplicated(oral_microbiome_sample_info$sample_id))
# oral_microbiome_sample_info$sample_id[910]
# oral_microbiome_sample_info = oral_microbiome_sample_info[-910,]
# oral_microbiome_expression_data = oral_microbiome_expression_data[,-910]

temp_data <-
  apply(oral_microbiome_expression_data, 1, function(x){
    (x - mean(x))/sd(x)
  })

rownames(temp_data) == oral_microbiome_sample_info$sample_id

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

pvca_object_genus <-
  PVCA(counts = counts, meta = meta, threshold = 0.6, inter = FALSE)

pvca_object_genus
save(pvca_object_genus, file = "pvca_object_genus")
load("pvca_object_genus")
plot <- 
  PlotPVCA(pvca.res = pvca_object_genus, title = "")
plot
ggsave(plot, filename = "pvca_plot_genus.pdf", height = 7, width = 7)



######order level
{
  oral_microbiome_expression_data = 
    data.table::fread(here::here("data/from_xin/Genus Table/OR/Order_OR.csv")) %>% 
    dplyr::select(SampleID, everything()) %>% 
    dplyr::select(-c(V1:SubjectID)) %>% 
    as.data.frame() %>% 
    tibble::column_to_rownames(var = "SampleID") %>% 
    t() %>% 
    as.data.frame()
  
  oral_microbiome_expression_data = oral_microbiome_expression_data[,oral_microbiome_sample_info$sample_id]
  
}

#######work directory
masstools::setwd_project()
setwd("data_analysis/oral_microbiome/PVCA_analysis/")

zero_percent =
  oral_microbiome_expression_data %>%
  apply(1, function(x) {
    sum(x == 0) / ncol(oral_microbiome_expression_data)
  })

sum(zero_percent > 0.99)/nrow(oral_microbiome_expression_data)

##here we remove the order with 0 > 99%
remain_idx = which(zero_percent < 0.99)

length(remain_idx)

oral_microbiome_expression_data = 
  oral_microbiome_expression_data[remain_idx,]

##remove duplicated sample
which(duplicated(oral_microbiome_sample_info$sample_id))
# oral_microbiome_sample_info$sample_id[910]
# oral_microbiome_sample_info = oral_microbiome_sample_info[-910,]
# oral_microbiome_expression_data = oral_microbiome_expression_data[,-910]

temp_data <-
  apply(oral_microbiome_expression_data, 1, function(x){
    (x - mean(x))/sd(x)
  })

rownames(temp_data) == oral_microbiome_sample_info$sample_id

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

pvca_object_order <-
  PVCA(counts = counts, meta = meta, threshold = 0.6, inter = FALSE)

pvca_object_order
save(pvca_object_order, file = "pvca_object_order")
load("pvca_object_order")
plot <- 
  PlotPVCA(pvca.res = pvca_object_order, title = "")
plot
ggsave(plot, filename = "pvca_plot_order.pdf", height = 7, width = 7)

