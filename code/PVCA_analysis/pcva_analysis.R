##Principal Variance Component Analysis
###don't use this one, use the pvca_analysis2.R
no_source()

# set work directory
masstools::setwd_project()
library(tidyverse)
library(phyloseq)
rm(lins = ls())

source("code/tools.R")

####nasal microbiome
{
  load("data_analysis/nasal_microbiome/data_preparation/sample_info")
  load("data_analysis/nasal_microbiome/data_preparation/expression_data")
  load("data_analysis/nasal_microbiome/data_preparation/variable_info")
  
  dim(variable_info)
  
  nasal_microbiome_sample_info = sample_info
  nasal_microbiome_expression_data = expression_data
  nasal_microbiome_variable_info = variable_info
  
  nasal_microbiome_sample_info =
    nasal_microbiome_sample_info %>%
    dplyr::mutate(SSPG = as.numeric(SSPG)) %>%
    dplyr::mutate(iris = case_when(SSPG > 125 ~ "IR",
                                   SSPG < 125 ~ "IS"))
  
  library(lubridate)
  
  nasal_microbiome_sample_info$days =
    nasal_microbiome_sample_info$Date %>%
    lubridate::as_date() %>%
    yday()
  
  nasal_microbiome_sample_info =
    nasal_microbiome_sample_info %>%
    dplyr::mutate(
      days =  yday(lubridate::as_date(Date)),
      months = lubridate::month(lubridate::as_date(nasal_microbiome_sample_info$Date)),
      weeks = lubridate::week(lubridate::as_date(nasal_microbiome_sample_info$Date))
    )
}

zero_percent =
  nasal_microbiome_expression_data %>%
  apply(1, function(x) {
    sum(x == 0) / ncol(nasal_microbiome_expression_data)
  })

sum(zero_percent > 0.99) / nrow(nasal_microbiome_variable_info)

##here we remove the genus with 0 > 99%
remain_idx = which(zero_percent < 0.99)

length(remain_idx)

nasal_microbiome_expression_data =
  nasal_microbiome_expression_data[remain_idx, ]

nasal_microbiome_variable_info =
  nasal_microbiome_variable_info[remain_idx, ]

##remove duplicated sample
which(duplicated(nasal_microbiome_sample_info$sample_id))
# nasal_microbiome_sample_info$sample_id[910]
# nasal_microbiome_sample_info = nasal_microbiome_sample_info[-910,]
# nasal_microbiome_expression_data = nasal_microbiome_expression_data[,-910]


####oral microbiome
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
    dplyr::mutate(iris = case_when(SSPG > 125 ~ "IR",
                                   SSPG < 125 ~ "IS"))
  
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

zero_percent =
  oral_microbiome_expression_data %>%
  apply(1, function(x) {
    sum(x == 0) / ncol(oral_microbiome_expression_data)
  })

sum(zero_percent > 0.99) / nrow(oral_microbiome_variable_info)

##here we remove the genus with 0 > 99%
remain_idx = which(zero_percent < 0.99)

length(remain_idx)

oral_microbiome_expression_data =
  oral_microbiome_expression_data[remain_idx, ]

oral_microbiome_variable_info =
  oral_microbiome_variable_info[remain_idx, ]

##remove duplicated sample
which(duplicated(oral_microbiome_sample_info$sample_id))
# oral_microbiome_sample_info$sample_id[910]
# oral_microbiome_sample_info = oral_microbiome_sample_info[-910,]
# oral_microbiome_expression_data = oral_microbiome_expression_data[,-910]


###skin microbiome
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
    dplyr::mutate(iris = case_when(SSPG > 125 ~ "IR",
                                   SSPG < 125 ~ "IS"))
  
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

zero_percent =
  skin_microbiome_expression_data %>%
  apply(1, function(x) {
    sum(x == 0) / ncol(skin_microbiome_expression_data)
  })

sum(zero_percent > 0.99) / nrow(skin_microbiome_variable_info)

##here we remove the genus with 0 > 99%
remain_idx = which(zero_percent < 0.99)

length(remain_idx)

skin_microbiome_expression_data =
  skin_microbiome_expression_data[remain_idx, ]

skin_microbiome_variable_info =
  skin_microbiome_variable_info[remain_idx, ]

##remove duplicated sample
which(duplicated(skin_microbiome_sample_info$sample_id))
skin_microbiome_sample_info$sample_id[910]
skin_microbiome_sample_info = skin_microbiome_sample_info[-910, ]
skin_microbiome_expression_data = skin_microbiome_expression_data[, -910]


###stool microbiome
{
  load("data_analysis/stool_microbiome/data_preparation/sample_info")
  load("data_analysis/stool_microbiome/data_preparation/expression_data")
  load("data_analysis/stool_microbiome/data_preparation/variable_info")
  
  dim(variable_info)
  
  stool_microbiome_sample_info = sample_info
  stool_microbiome_expression_data = expression_data
  stool_microbiome_variable_info = variable_info
  
  stool_microbiome_sample_info =
    stool_microbiome_sample_info %>%
    dplyr::mutate(SSPG = as.numeric(SSPG)) %>%
    dplyr::mutate(iris = case_when(SSPG > 125 ~ "IR",
                                   SSPG < 125 ~ "IS"))
  
  library(lubridate)
  
  stool_microbiome_sample_info$days =
    stool_microbiome_sample_info$Date %>%
    lubridate::as_date() %>%
    yday()
  
  stool_microbiome_sample_info =
    stool_microbiome_sample_info %>%
    dplyr::mutate(
      days =  yday(lubridate::as_date(Date)),
      months = lubridate::month(lubridate::as_date(stool_microbiome_sample_info$Date)),
      weeks = lubridate::week(lubridate::as_date(stool_microbiome_sample_info$Date))
    )
}

zero_percent =
  stool_microbiome_expression_data %>%
  apply(1, function(x) {
    sum(x == 0) / ncol(stool_microbiome_expression_data)
  })

sum(zero_percent > 0.99) / nrow(stool_microbiome_variable_info)

##here we remove the genus with 0 > 99%
remain_idx = which(zero_percent < 0.99)

length(remain_idx)

stool_microbiome_expression_data =
  stool_microbiome_expression_data[remain_idx, ]

stool_microbiome_variable_info =
  stool_microbiome_variable_info[remain_idx, ]


#####set work directory
masstools::setwd_project()
setwd("data_analysis/pcva_analysis")

dim(oral_microbiome_expression_data)
dim(skin_microbiome_expression_data)
dim(stool_microbiome_expression_data)
dim(nasal_microbiome_expression_data)

colnames(stool_microbiome_expression_data)
colnames(skin_microbiome_expression_data)

intersect_sample_id =
  Reduce(intersect, x = list(
    colnames(oral_microbiome_expression_data),
    colnames(skin_microbiome_expression_data),
    colnames(stool_microbiome_expression_data),
    colnames(nasal_microbiome_expression_data)
  ))

library(plyr)
library(lubridate)
library(mgcv)
library(gratia)
library(pvca)

temp_data <-
  rbind(
    oral_microbiome_expression_data[, intersect_sample_id],
    skin_microbiome_expression_data[, intersect_sample_id],
    stool_microbiome_expression_data[, intersect_sample_id],
    nasal_microbiome_expression_data[, intersect_sample_id]
  )

zero_percent =
  temp_data %>%
  apply(1, function(x) {
    sum(x == 0) / ncol(temp_data)
  })

sum(zero_percent > 0.99) / nrow(temp_data)

##here we remove the genus with 0 > 99%
remain_idx = which(zero_percent < 0.99)

length(remain_idx)

temp_data =
  temp_data[remain_idx, ]

temp_data =
  temp_data %>%
  apply(1, function(x) {
    (x - mean(x)) / sd(x)
  })

sample_info =
  oral_microbiome_sample_info[match(intersect_sample_id, oral_microbiome_sample_info$sample_id), ]

rownames(temp_data) == sample_info$sample_id

# ####pvca analysis
library(pvca)
assay_data <-
  as.matrix(t(temp_data))

pheno_data <-
  sample_info %>%
  dplyr::select(subject_id,
                sample_id,
                days) %>%
  as.data.frame()

library(Biobase)

pheno_data <- AnnotatedDataFrame(data = pheno_data)

row.names(pheno_data) <- colnames(assay_data)

data_pvca <-
  Biobase::ExpressionSet(assayData = assay_data,
                         phenoData = pheno_data)

pct_threshold <- 0.6

batch.factors <-
  c("subject_id", #yes
    "days")

pvcaObj <- pvcaBatchAssess(abatch = data_pvca,
                           batch.factors = batch.factors,
                           threshold = pct_threshold)

save(pvcaObj, file = "pvcaObj")
load("pvcaObj")

plot <-
  plot_pvca(object = pvcaObj)

plot

ggsave(plot,
       filename = "pvca_plot.pdf",
       height = 9,
       width = 7)

counts <- as.matrix(t(temp_data))

meta <- sample_info %>%
  dplyr::select(sample_id,
                subject_id,
                days) %>%
  as.data.frame()

rownames(meta) = NULL

meta =
  meta %>%
  tibble::column_to_rownames(var = "sample_id")

pvca_object <-
  PVCA(
    counts = counts,
    meta = meta,
    threshold = 0.6,
    inter = FALSE
  )

pvca_object
save(pvca_object, file = "pvca_object")
load("pvca_object")
plot <-
  PlotPVCA(pvca.res = pvca_object, title = "")
plot
ggsave(plot,
       filename = "pvca_plot.pdf",
       height = 7,
       width = 7)



####ggtern

{
  ##load data
  load(here::here(
    "data_analysis/nasal_microbiome/season_analysis/pvca_object"
  ))
  nasal_object <- pvca_object
  
  load(here::here(
    "data_analysis/skin_microbiome/season_analysis/pvca_object"
  ))
  skin_object <- pvca_object
  
  load(here::here(
    "data_analysis/stool_microbiome/season_analysis/pvca_object"
  ))
  stool_object <- pvca_object
  
  load(here::here(
    "data_analysis/oral_microbiome/season_analysis/pvca_object"
  ))
  oral_object <- pvca_object
  
  load(here::here("data_analysis/pcva_analysis/pvca_object"))
  total_object <- pvca_object
}


temp_data <-
  rbind(nasal_object,
        stool_object,
        skin_object,
        oral_object,
        total_object) %>%
  as.data.frame() %>%
  dplyr::rename(Subject = subject_id,
                Days = days,
                Residual = resid)
rownames(temp_data) = c("Nasal", "Stool", "Skin", "Oral", "Total")

temp_data =
  temp_data %>%
  tibble::rownames_to_column(var = "class")

library(ggtern)

plot <-
  ggtern::ggtern(data = temp_data, aes(x = Days, y = Subject, z = Residual)) +
  ggtern::theme_bvbw() +
  geom_point(aes(fill = class),
             shape = 21,
             size = 6,
             alpha = 0.8) +
  geom_text(aes(label = class)) +
  scale_fill_manual(values = body_site_color) +
  labs(x = "Days",
       y = "Subject",
       z = "Residual",
       title = "") +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    plot.background = element_rect(fill = "transparent"),
    panel.background = element_rect(fill = "transparent")
  )

plot

ggsave(plot,
       filename = "trenery_plot.pdf",
       width = 7,
       height = 7)
