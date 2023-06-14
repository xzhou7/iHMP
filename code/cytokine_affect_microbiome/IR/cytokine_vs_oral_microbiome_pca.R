###https://biofam.github.io/MOFA2/index.html
##https://genomebiology.biomedcentral.com/track/pdf/10.1186/s13059-020-02015-1.pdf
no_source()

# set work directory
setwd(masstools::get_project_wd())
library(tidyverse)
library(phyloseq)
rm(list = ls())

source("code/tools.R")

######work directory
setwd(masstools::get_project_wd())

dir.create(
  "data_analysis/cytokine_affect_microbiome/IR/cytokine_oral_microbiome",
  recursive = TRUE
)

setwd("data_analysis/cytokine_affect_microbiome/IR/cytokine_oral_microbiome")

###load data
##oral microbiome
{
  load(here::here("data_analysis/oral_microbiome/data_preparation/expression_data"))
  load(here::here("data_analysis/oral_microbiome/data_preparation/sample_info"))
  load(here::here("data_analysis/oral_microbiome/data_preparation/variable_info"))
}

oral_microbiome_expression_data = expression_data
oral_microbiome_sample_info = sample_info
oral_microbiome_variable_info = variable_info

###read genus table
expression_data =
  data.table::fread(here::here("data/from_xin/Genus Table/OR/Genus_OR.csv")) %>%
  as.data.frame() %>%
  tibble::column_to_rownames(var = "SampleID") %>%
  dplyr::select(-c(V1:SubjectID)) %>%
  t() %>%
  as.data.frame()

oral_microbiome_variable_info =
  oral_microbiome_variable_info[match(rownames(expression_data), oral_microbiome_variable_info$Genus),]

oral_microbiome_variable_info$Genus == rownames(expression_data)

###remove the variables which Genus are NA
remove_idx = which(is.na(oral_microbiome_variable_info$Genus))
remove_idx
if (length(remove_idx) > 0) {
  oral_microbiome_variable_info = oral_microbiome_variable_info[-remove_idx,]
  expression_data = expression_data[-remove_idx,]
}

rownames(expression_data) = oral_microbiome_variable_info$variable_id

oral_microbiome_expression_data =
  expression_data

oral_microbiome_variable_info =
  oral_microbiome_variable_info %>%
  dplyr::filter(!stringr::str_detect(Genus, "Unclassified_Bacteria"))

oral_microbiome_expression_data =
  oral_microbiome_expression_data[oral_microbiome_variable_info$variable_id,]

dim(oral_microbiome_sample_info)
dim(oral_microbiome_variable_info)

rownames(oral_microbiome_expression_data) == oral_microbiome_variable_info$variable_id
colnames(oral_microbiome_expression_data) == oral_microbiome_sample_info$sample_id

oral_microbiome_sample_info =
  oral_microbiome_sample_info %>%
  dplyr::filter(IRIS == "IR")

oral_microbiome_expression_data =
  oral_microbiome_expression_data[, oral_microbiome_sample_info$sample_id]

###plasma cytokine
{
  load(here::here(
    "data_analysis/cytokine/data_preparation/expression_data"
  ))
  load(here::here("data_analysis/cytokine/data_preparation/sample_info"))
  load(here::here("data_analysis/cytokine/data_preparation/variable_info"))
}

cytokine_expression_data = t(expression_data) %>%
  as.data.frame()
cytokine_sample_info = sample_info %>%
  dplyr::rename(sample_id = SampleID,
                subject_id = SubjectID)

cytokine_variable_info = variable_info

cytokine_sample_info$CollectionDate =
  as.Date(cytokine_sample_info$CollectionDate, "%m/%d/%y")

dim(cytokine_expression_data)
length(unique(cytokine_sample_info$subject_id))

###match samples
dim(oral_microbiome_sample_info)
dim(cytokine_sample_info)

length(oral_microbiome_sample_info$subject_id)
length(unique(oral_microbiome_sample_info$subject_id))

###just matched samples according to sample id, only 1 missed
intersect_sample_id =
  intersect(oral_microbiome_sample_info$sample_id,
            cytokine_sample_info$sample_id)

length(intersect_sample_id)

oral_microbiome_expression_data =
  oral_microbiome_expression_data[, intersect_sample_id]

cytokine_expression_data =
  cytokine_expression_data[, intersect_sample_id]

oral_microbiome_sample_info =
  oral_microbiome_sample_info[match(intersect_sample_id, oral_microbiome_sample_info$sample_id),]

cytokine_sample_info =
  cytokine_sample_info[match(intersect_sample_id, cytokine_sample_info$sample_id),]

length(unique(cytokine_sample_info$subject_id))

###only remain the subjects with at least >= 5
remian_subject_id =
  oral_microbiome_sample_info %>%
  dplyr::group_by(subject_id) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::ungroup() %>%
  dplyr::filter(n > 5) %>%
  dplyr::pull(subject_id)

cytokine_sample_info =
  cytokine_sample_info %>%
  dplyr::filter(subject_id %in% remian_subject_id)

cytokine_expression_data =
  cytokine_expression_data[, cytokine_sample_info$sample_id]

oral_microbiome_sample_info =
  oral_microbiome_sample_info %>%
  dplyr::filter(subject_id %in% remian_subject_id)

oral_microbiome_expression_data =
  oral_microbiome_expression_data[, oral_microbiome_sample_info$sample_id]

##only remain the genus at least in 10% subjects
remain_idx =
  which(rowSums(oral_microbiome_expression_data) > 0)

oral_microbiome_expression_data = oral_microbiome_expression_data[remain_idx,]
oral_microbiome_variable_info = oral_microbiome_variable_info[remain_idx, , drop = FALSE]

remain_idx =
  oral_microbiome_expression_data %>%
  apply(1, function(x) {
    sum(as.numeric(x) == 0) / ncol(oral_microbiome_expression_data)
  }) %>%
  `<`(0.9) %>%
  which()

length(remain_idx)

oral_microbiome_expression_data = oral_microbiome_expression_data[remain_idx,]
oral_microbiome_variable_info = oral_microbiome_variable_info[remain_idx, , drop = FALSE]

##save data
  load(
    here::here(
      "data_analysis/mediation_analysis/sample_wise_IR/oral_microbiome_vs_cytokine/mediation_result"
    )
  )
  
  idx1 <-
    match(unique(mediation_result$treat), oral_microbiome_variable_info$variable_id)
  idx2 <-
    match(unique(mediation_result$mediator), cytokine_variable_info$variable_id)
  idx1
  idx2
  
  cytokine_expression_data <-
    cytokine_expression_data[idx2,]
  cytokine_variable_info <-
    cytokine_variable_info[idx2,]
  
  save(oral_microbiome_expression_data, file = "oral_microbiome_expression_data")
  save(oral_microbiome_variable_info, file = "oral_microbiome_variable_info")
  save(oral_microbiome_sample_info, file = "oral_microbiome_sample_info")
  
  save(cytokine_expression_data, file = "cytokine_expression_data")
  save(cytokine_variable_info, file = "cytokine_variable_info")
  save(cytokine_sample_info, file = "cytokine_sample_info")

#####just use the pca to do dimension reduction for microbiome
temp_expression_data <-
  log(oral_microbiome_expression_data + 1, 2)

temp_expression_data <-
  temp_expression_data %>%
  apply(1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()

rownames(temp_expression_data)

#PCA analysis
###PCA for date
pca_object <-
  prcomp(
    x = t(temp_expression_data),
    center = FALSE,
    scale. = FALSE
  )

idx <-
  which(summary(pca_object)$importance[3,] > 0.8)[1]

temp_data =
  summary(pca_object)$importance[, 1:31] %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "class") %>%
  dplyr::filter(class != "Standard deviation") %>%
  tidyr::pivot_longer(cols = -class,
                      names_to = "PC",
                      values_to = "value") %>%
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
  geom_point(
    aes(PC, value * 100),
    size = 5,
    data = temp_data %>% dplyr::filter(class == "Cumulative Proportion")
  ) +
  # geom_text(aes(PC, value * 100 + 2,
  #               label = panse(round(value*100, 2), "%", sep = "")),
  #           data = temp_data) +
  base_theme +
  theme(axis.text.x = element_text(
    size = 12,
    angle = 45,
    hjust = 1,
    vjust = 1
  ))
plot
# ggsave(plot, filename = "oral_microbiome_pca_pc.pdf",
#        width = 9, height = 7)

x = pca_object$x[, c(1:idx)] %>%
  t() %>%
  as.data.frame()

colnames(x) == oral_microbiome_sample_info$sample_id

oral_microbiome_expression_data = x

#####juns use the pca to do dimension reduction for cytokine
temp_expression_data <-
  log(cytokine_expression_data + 1, 2)

temp_expression_data <-
  temp_expression_data %>%
  apply(1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()

rownames(temp_expression_data)

#PCA analysis
pca_object <-
  prcomp(
    x = t(temp_expression_data),
    center = FALSE,
    scale. = FALSE
  )

idx <-
  which(summary(pca_object)$importance[3,] > 0.8)[1]

temp_data =
  summary(pca_object)$importance[, 1:10] %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "class") %>%
  dplyr::filter(class != "Standard deviation") %>%
  tidyr::pivot_longer(cols = -class,
                      names_to = "PC",
                      values_to = "value") %>%
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
  geom_point(
    aes(PC, value * 100),
    size = 5,
    data = temp_data %>% dplyr::filter(class == "Cumulative Proportion")
  ) +
  # geom_text(aes(PC, value * 100 + 2,
  #               label = panse(round(value*100, 2), "%", sep = "")),
  #           data = temp_data) +
  base_theme +
  theme(axis.text.x = element_text(
    size = 12,
    angle = 45,
    hjust = 1,
    vjust = 1
  ))
plot
# ggsave(plot, filename = "exposome_chemical_pca_pc.pdf",
#        width = 9, height = 7)

x = pca_object$x[, c(1:idx)] %>%
  t() %>%
  as.data.frame()

colnames(x) == cytokine_sample_info$sample_id

cytokine_expression_data = x

####multiple linear regression
total_r2 <-
  oral_microbiome_expression_data %>%
  t() %>%
  as.data.frame() %>%
  purrr::map(
    .f = function(x) {
      temp_data <-
        rbind(y = x,
              cytokine_expression_data) %>%
        t() %>%
        as.data.frame()
      
      lm_object <-
        lm(formula = y ~ .,
           data = temp_data)
      summary(lm_object)$r.squared
    }
  ) %>%
  unlist()
save(total_r2, file = "total_r2")
load("total_r2")
