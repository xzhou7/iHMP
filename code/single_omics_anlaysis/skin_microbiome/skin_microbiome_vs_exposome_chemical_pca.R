###https://biofam.github.io/MOFA2/index.html
##https://genomebiology.biomedcentral.com/track/pdf/10.1186/s13059-020-02015-1.pdf
no_source()

# library(data.table)
# library(MOFA2)
# 
# data <- make_example_data(
#   n_views = 2, 
#   n_samples = 200, 
#   n_features = 1000, 
#   n_factors = 10
# )[[1]]
# 
# lapply(data,dim)
# 
# MOFAobject <- create_mofa(data)
# 
# N = ncol(data[[1]])
# groups = c(rep("A",N/2), rep("B",N/2))
# 
# MOFAobject <- create_mofa(data, groups=groups)
# 
# dt = fread("ftp://ftp.ebi.ac.uk/pub/databases/mofa/getting_started/data.txt.gz")
# head(dt)
# 
# dt[,group:=NULL]
# 
# MOFAobject <- create_mofa(dt)
# 
# print(MOFAobject)
# 
# 
# plot_data_overview(MOFAobject)
# 
# data_opts <- get_default_data_options(MOFAobject)
# head(data_opts)
# 
# model_opts <- get_default_model_options(MOFAobject)
# head(model_opts)
# 
# train_opts <- get_default_training_options(MOFAobject)
# head(train_opts)
# 
# MOFAobject <- prepare_mofa(
#   object = MOFAobject,
#   data_options = data_opts,
#   model_options = model_opts,
#   training_options = train_opts
# )
# 
# 
# outfile = file.path(getwd(),"model.hdf5")
# MOFAobject.trained <- run_mofa(MOFAobject, outfile)

# set work directory
masstools::setwd_project()
library(tidyverse)
library(phyloseq)
rm(list = ls())

source("code/tools.R")

###load skin_microbiome
{
  load("data_analysis/skin_microbiome/data_preparation/sample_info")
  load("data_analysis/skin_microbiome/data_preparation/expression_data")
  load("data_analysis/skin_microbiome/data_preparation/variable_info")
  
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
  
  ###only remain mike's samples
  skin_microbiome_sample_info =
    skin_microbiome_sample_info %>%
    dplyr::filter(subject_id == "69-001")
  
  skin_microbiome_expression_data = 
    skin_microbiome_expression_data[,skin_microbiome_sample_info$sample_id]
  

  zero_percent =
    skin_microbiome_expression_data %>%
    apply(1, function(x){
      sum(x == 0)/ncol(skin_microbiome_expression_data)
    })

  sum(zero_percent > 0.99)/nrow(skin_microbiome_variable_info)

  ##here we remove the genus with 0 > 99%
  remain_idx = which(zero_percent < 0.99)

  skin_microbiome_expression_data = skin_microbiome_expression_data[remain_idx,]
  skin_microbiome_variable_info = skin_microbiome_variable_info[remain_idx,]
  
  }

##load exposome_chemical
###exposome_chemical
{
  load("data_analysis/exposome_chemical/expression_data")
  exposome_chemical_expression_data <- expression_data
  load("data_analysis/exposome_chemical/sample_info")
  exposome_chemical_sample_info <- sample_info
  load("data_analysis/exposome_chemical/variable_info")
  exposome_chemical_variable_info <- variable_info  
  
  library(plyr)
  exposome_chemical_expression_data <-
    exposome_chemical_expression_data %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column(var = "Sample_ID") %>%
    dplyr::mutate(Sample_ID = stringr::str_replace(Sample_ID, "_[0-9]{1,2}", "")) %>%
    plyr::dlply(.variables = "Sample_ID") %>%
    lapply(function(x){
      apply(x[,-1], 2, mean)
    }) %>%
    do.call(rbind, .) %>%
    t() %>%
    as.data.frame()
}


#######work directory
masstools::setwd_project()
setwd("data_analysis/skin_microbiome/skin_microbiome_exposome_chemical")

#####just use the pca to do dimension reduction for microbiome
temp_expression_data <-
  log(skin_microbiome_expression_data + 1, 2)

temp_expression_data <- 
  temp_expression_data %>% 
  apply(1, function(x){
    (x - mean(x))/sd(x)
  }) %>% 
  t() %>% 
  as.data.frame()

rownames(temp_expression_data)

#PCA analysis
###PCA for date
pca_object <- 
  prcomp(x = t(temp_expression_data), center = FALSE, scale. = FALSE)

idx <-
  which(summary(pca_object)$importance[3, ] > 0.8)[1]

temp_data = 
summary(pca_object)$importance[,1:37] %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "class") %>% 
  dplyr::filter(class != "Standard deviation") %>% 
  tidyr::pivot_longer(cols = -class, names_to = "PC", values_to = "value") %>% 
  dplyr::mutate(PC = factor(PC, levels = unique(PC)))
  
plot = 
ggplot(data = temp_data %>% dplyr::filter(class == "Proportion of Variance"),
       aes(PC, value * 100)) +
  geom_bar(stat = "identity") +
  labs(x = "Principle component (PC)",
       y = "(Cumulative) Proportion of Variance (%)") +
  geom_line(
    aes(PC, value * 100, group = class),
    data = temp_data %>% dplyr::filter(class == "Cumulative Proportion")
  ) +
  geom_point(aes(PC, value * 100),
             size = 5,
             data = temp_data %>% dplyr::filter(class == "Cumulative Proportion"))+
  # geom_text(aes(PC, value * 100 + 2, 
  #               label = paste(round(value*100, 2), "%", sep = "")), 
  #           data = temp_data) +
  base_theme +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1))
plot
ggsave(plot, filename = "skin_microbiome_pca_pc.pdf",
       width = 9, height = 7)


x = pca_object$x[,c(1:idx)] %>% 
  t() %>% 
  as.data.frame()

colnames(x) == skin_microbiome_sample_info$sample_id

skin_microbiome_expression_data = x


#####just use the pca to do dimension reduction for exposome_chemical
temp_expression_data <-
  log(exposome_chemical_expression_data + 1, 2)

temp_expression_data <- 
  temp_expression_data %>% 
  apply(1, function(x){
    (x - mean(x))/sd(x)
  }) %>% 
  t() %>% 
  as.data.frame()

rownames(temp_expression_data)

#PCA analysis
pca_object <- 
  prcomp(x = t(temp_expression_data), center = FALSE, scale. = FALSE)

idx <-
  which(summary(pca_object)$importance[3, ] > 0.8)[1]

temp_data = 
  summary(pca_object)$importance[,1:10] %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "class") %>% 
  dplyr::filter(class != "Standard deviation") %>% 
  tidyr::pivot_longer(cols = -class, names_to = "PC", values_to = "value") %>% 
  dplyr::mutate(PC = factor(PC, levels = unique(PC)))


plot = 
  ggplot(data = temp_data %>% dplyr::filter(class == "Proportion of Variance"),
         aes(PC, value * 100)) +
  geom_bar(stat = "identity") +
  labs(x = "Principle component (PC)",
       y = "(Cumulative) Proportion of Variance (%)") +
  geom_line(
    aes(PC, value * 100, group = class),
    data = temp_data %>% dplyr::filter(class == "Cumulative Proportion")
  ) +
  geom_point(aes(PC, value * 100),
             size = 5,
             data = temp_data %>% dplyr::filter(class == "Cumulative Proportion"))+
  # geom_text(aes(PC, value * 100 + 2, 
  #               label = paste(round(value*100, 2), "%", sep = "")), 
  #           data = temp_data) +
  base_theme +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1))
plot
ggsave(plot, filename = "exposome_chemical_pca_pc.pdf",
       width = 9, height = 7)


x = pca_object$x[,c(1:idx)] %>% 
  t() %>% 
  as.data.frame()

colnames(x) == exposome_chemical_sample_info$sample_id

exposome_chemical_expression_data = x


#####match data
{
  sort(skin_microbiome_sample_info$Date)
  sort(exposome_chemical_sample_info$start_date)
  
  temp1 <- 
    exposome_chemical_sample_info$start_date %>% 
    as.Date() %>% 
    data.frame(date = ., 
               exposome_chemical = 1,
               stringsAsFactors = FALSE) %>% 
    dplyr::arrange(date)
  
  temp2 <- 
    skin_microbiome_sample_info$Date %>% 
    as.Date() %>% 
    data.frame(date = ., 
               skin_microbiome = 1,
               stringsAsFactors = FALSE) %>% 
    dplyr::arrange(date)
  
  temp <- 
    temp1 %>% 
    dplyr::full_join(temp2, by = "date") %>% 
    arrange(date)
  
  temp <- 
    temp %>% 
    tidyr::pivot_longer(cols = -date, names_to = "class", values_to = "value")
  
  temp <- temp %>% 
    dplyr::filter(!is.na(value))
  
  diff <- 
    lapply(skin_microbiome_sample_info$Date %>% as.Date(), function(x){
      x <- abs(x - as.Date(exposome_chemical_sample_info$start_date))
      as.Date(exposome_chemical_sample_info$start_date)[which(x >=0 & x <= 3)]
    }) 
  
  names(diff) <- skin_microbiome_sample_info$Date %>% as.character()
  
  diff =
    diff[unlist(lapply(diff, length)) != 0]
  
  diff = diff[sort(names(diff))]
  
  ####remove some duplicated samples
  diff$`2016-01-15` <- diff$`2016-01-15`[1]
  diff$`2016-01-19` <- diff$`2016-01-19`[2]
  diff$`2016-02-16` <- diff$`2016-02-16`[1]
  diff$`2016-02-24` <- diff$`2016-02-24`[2]
  
  
  diff <- 
    purrr::map2(diff, names(diff), .f = function(x,y){
      if(length(x) == 0){
        return(NULL)
      }
      z <- data.frame(y, x, stringsAsFactors = FALSE)
      colnames(z) <- c("skin_microbiome", "exposome_chemical")
      z
    }) 
  
  diff <- purrr::map2(.x = diff, .y = 1:length(diff), .f = function(x,y){
    data.frame(x, group = y, stringsAsFactors = FALSE)
  }) %>% 
    do.call(rbind, .)
  
  rownames(diff) <- NULL
  
  diff <- 
    diff %>% 
    mutate(exposome_chemical = as.character(exposome_chemical), 
           skin_microbiome = as.character(skin_microbiome)) %>% 
    tidyr::pivot_longer(cols = -group,
                        names_to = "class", 
                        values_to = "date")
  
  temp <-
    temp %>% 
    mutate(date = as.character(date)) %>% 
    left_join(diff, by = c("date","class"))
  
  temp$group[is.na(temp$group)] <- "No"
  
  diff
  
  library(plyr)
  
  diff <-
    diff %>% 
    plyr::dlply(.variables = .(group)) %>%
    lapply(function(x){
      x %>% 
        plyr::dlply(.variables = .(class)) %>% 
        do.call(cbind, .)
    }) %>% 
    do.call(rbind, .)
  
  diff$exposome_chemical.date <- as.Date(diff$exposome_chemical.date)
  
  diff$skin_microbiome.date <- as.Date(diff$skin_microbiome.date)
  
  diff$exposome_chemical.group <- as.character(diff$exposome_chemical.group)
  
  temp$class[temp$class == "exposome_chemical"] = "Exposome"
  temp$class[temp$class == "skin_microbiome"] = "Skin microbiome"
  diff$exposome_chemical.class[diff$exposome_chemical.class == "exposome_chemical"] = "Exposome"
  diff$skin_microbiome.class[diff$skin_microbiome.class == "skin_microbiome"] = "Skin microbiome"

}

temp_data = 
  temp %>%
  mutate(date = as.Date(date)) %>%
  dplyr::filter(date >= "2016-01-12" & date < "2016-04-06") %>%
  dplyr::mutate(class = factor(class,
                               levels = c("Exposome", "Skin microbiome"))) 

plot <-
  ggplot(aes(date, class), data = temp_data) +
  geom_segment(
    aes(
      x = exposome_chemical.date,
      xend = skin_microbiome.date,
      y = exposome_chemical.class,
      yend = skin_microbiome.class
    ),
    color = "black",
    data = diff,
    show.legend = FALSE
  ) +
  geom_point(
    aes(fill = class),
    color = "black",
    alpha = 1,
    size = 4,
    shape = 21,
    show.legend = FALSE,
    data = temp_data
  ) +
  scale_fill_manual(values = omics_color) +
  scale_x_continuous(
    trans = "date",
    breaks = c(as.Date(temp$date)),
    labels = as.character(temp$date)
  ) +
  labs(x = "", y = "") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 13),
    axis.text.x = element_text(
      size = 10,
      angle = 45,
      vjust = 1,
      hjust = 1
    ),
    axis.text = element_text(size = 10),
    panel.grid.minor = element_blank()
  )

plot

# ggsave(plot, filename = "exposome_chemical_skin_microbiome_match.pdf", width = 10, height = 7)


#-------------------------------------------------------------------------------
###prepare data
dim(skin_microbiome_expression_data)
dim(skin_microbiome_sample_info)

skin_microbiome_sample_info <-
  skin_microbiome_sample_info %>%
  dplyr::filter(Date %in% unique(diff$skin_microbiome.date))

skin_microbiome_sample_info

skin_microbiome_expression_data <-
  skin_microbiome_expression_data %>%
  dplyr::select(one_of(skin_microbiome_sample_info$sample_id))

skin_microbiome_sample_info <-
  skin_microbiome_sample_info %>%
  dplyr::left_join(diff[,c("skin_microbiome.date", "skin_microbiome.group")] %>% dplyr::distinct(skin_microbiome.date, skin_microbiome.group),
                   by = c("Date" = "skin_microbiome.date")) %>%
  dplyr::rename(group = skin_microbiome.group)

exposome_chemical_sample_info <-
  exposome_chemical_sample_info %>%
  dplyr::filter(as.Date(start_date) %in% diff$exposome_chemical.date)

###browser
exposome_chemical_expression_data <-
  exposome_chemical_expression_data %>%
  dplyr::select(one_of(exposome_chemical_sample_info$sample_id))

exposome_chemical_sample_info <-
  exposome_chemical_sample_info %>%
  mutate(start_date = as.Date(start_date)) %>%
  dplyr::left_join(diff[,c("exposome_chemical.date", "exposome_chemical.group")], by = c("start_date" = "exposome_chemical.date")) %>%
  dplyr::rename(group = exposome_chemical.group)

dim(skin_microbiome_sample_info)
dim(exposome_chemical_sample_info)

###sort date
skin_microbiome_sample_info$sample_id == colnames(skin_microbiome_expression_data)
exposome_chemical_sample_info$sample_id == colnames(exposome_chemical_expression_data)

skin_microbiome_sample_info =
skin_microbiome_sample_info %>%
  dplyr::arrange(Date)

skin_microbiome_expression_data =
  skin_microbiome_expression_data[,skin_microbiome_sample_info$sample_id]

exposome_chemical_sample_info =
  exposome_chemical_sample_info %>%
  dplyr::arrange(start_date)

exposome_chemical_expression_data =
  exposome_chemical_expression_data[,exposome_chemical_sample_info$sample_id]

exposome_chemical_sample_info =
exposome_chemical_sample_info %>%
  dplyr::mutate(date = start_date)

skin_microbiome_sample_info =
  skin_microbiome_sample_info %>%
  dplyr::mutate(date = Date)

save(exposome_chemical_expression_data, file = "exposome_chemical_expression_data")
save(skin_microbiome_expression_data, file = "skin_microbiome_expression_data")

save(exposome_chemical_sample_info, file = "exposome_chemical_sample_info")
save(skin_microbiome_sample_info, file = "skin_microbiome_sample_info")


####correlation analysis
load("exposome_chemical_expression_data")
load("skin_microbiome_expression_data")
load("exposome_chemical_sample_info")
load("skin_microbiome_sample_info")


####multiple linear regression
total_r2 <-
skin_microbiome_expression_data %>%
  t() %>%
  as.data.frame() %>%
  purrr::map(.f = function(x){
    temp_data <-
      rbind(y = x,
            exposome_chemical_expression_data) %>%
      t() %>%
      as.data.frame()

    lm_object <-
      lm(
        formula = y ~ .,
        data = temp_data
      )
    summary(lm_object)$r.squared
  }) %>%
  unlist()
save(total_r2, file = "total_r2")
load("total_r2")
