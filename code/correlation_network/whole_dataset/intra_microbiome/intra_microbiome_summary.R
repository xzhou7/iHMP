#' ---
#' title: "oral microbiome nasal_microbiome correlation"
#' author:
#'   - name: "Xiaotao Shen"
#'     url: https://www.shenxt.info/
#'     affiliation: Stanford School of Medicine
#' date: "`r Sys.Date()`"
#' site: distill::distill_website
#' output:
#'   distill::distill_article:
#'     code_folding: false
#' ---

#+ r setup, echo=TRUE, eval = TRUE, include = TRUE

no_function()
# set work directory

masstools::setwd_project()
library(tidyverse)
rm(list = ls())

source("code/tools.R")

######work directory
masstools::setwd_project()

# ###load data
# ##stool vs skin
# {
#   load(
#     "data_analysis/correlation_network/whole_data_set/stool_microbiome_vs_skin_microbiome/stool_microbiome_skin_microbiome_lm_adjusted_cor_spearman"
#   )
#   stool_skin_sample_cor = stool_microbiome_skin_microbiome_lm_adjusted_cor_spearman
#
#   load(
#     "data_analysis/correlation_network/whole_data_set/stool_microbiome_vs_skin_microbiome/subject_wise/cor_data"
#   )
#
#   stool_skin_subject_cor = cor_data
#
#   ##stool vs oral
#   load(
#     "data_analysis/correlation_network/whole_data_set/stool_microbiome_vs_oral_microbiome/stool_microbiome_oral_microbiome_lm_adjusted_cor_spearman"
#   )
#   stool_oral_sample_cor = stool_microbiome_oral_microbiome_lm_adjusted_cor_spearman
#
#   load(
#     "data_analysis/correlation_network/whole_data_set/stool_microbiome_vs_oral_microbiome/subject_wise/cor_data"
#   )
#   stool_oral_subject_cor = cor_data
#
#   ##stool vs nasal
#   load(
#     "data_analysis/correlation_network/whole_data_set/stool_microbiome_vs_nasal_microbiome/stool_microbiome_nasal_microbiome_lm_adjusted_cor_spearman"
#   )
#   stool_nasal_sample_cor = stool_microbiome_nasal_microbiome_lm_adjusted_cor_spearman
#
#   load(
#     "data_analysis/correlation_network/whole_data_set/stool_microbiome_vs_nasal_microbiome/subject_wise/cor_data"
#   )
#   stool_nasal_subject_cor = cor_data
#
#
#   ##oral vs nasal
#   load(
#     "data_analysis/correlation_network/whole_data_set/oral_microbiome_vs_nasal_microbiome/oral_microbiome_nasal_microbiome_lm_adjusted_cor_spearman"
#   )
#   oral_nasal_sample_cor = oral_microbiome_nasal_microbiome_lm_adjusted_cor_spearman
#
#   load(
#     "data_analysis/correlation_network/whole_data_set/oral_microbiome_vs_nasal_microbiome/subject_wise/cor_data"
#   )
#   oral_nasal_subject_cor = cor_data
#
#   ##oral vs skin
#   load(
#     "data_analysis/correlation_network/whole_data_set/oral_microbiome_vs_skin_microbiome/oral_microbiome_skin_microbiome_lm_adjusted_cor_spearman"
#   )
#   oral_skin_sample_cor = oral_microbiome_skin_microbiome_lm_adjusted_cor_spearman
#
#   load(
#     "data_analysis/correlation_network/whole_data_set/oral_microbiome_vs_skin_microbiome/subject_wise/cor_data"
#   )
#   oral_skin_subject_cor = cor_data
#
#
#   ##skin vs nasal
#   load(
#     "data_analysis/correlation_network/whole_data_set/skin_microbiome_vs_nasal_microbiome/skin_microbiome_nasal_microbiome_lm_adjusted_cor_spearman"
#   )
#   skin_nasal_sample_cor = skin_microbiome_nasal_microbiome_lm_adjusted_cor_spearman
#
#   load(
#     "data_analysis/correlation_network/whole_data_set/skin_microbiome_vs_nasal_microbiome/subject_wise/cor_data"
#   )
#   skin_nasal_subject_cor = cor_data
# }
#
#
# ####intra-microbiome
# ####stool
# {
#   load(
#     "data_analysis/correlation_network/whole_data_set/intra_stool_microbiome/intra_stool_microbiome_lm_adjusted_cor"
#   )
#
#   stool_sample_cor <-
#     intra_stool_microbiome_lm_adjusted_cor
#
#   load(
#     "data_analysis/correlation_network/whole_data_set/intra_skin_microbiome/intra_skin_microbiome_lm_adjusted_cor"
#   )
#
#   skin_sample_cor <-
#     intra_skin_microbiome_lm_adjusted_cor
#
#   load(
#     "data_analysis/correlation_network/whole_data_set/intra_nasal_microbiome/intra_nasal_microbiome_lm_adjusted_cor"
#   )
#
#   nasal_sample_cor <-
#     intra_nasal_microbiome_lm_adjusted_cor
#
#   load(
#     "data_analysis/correlation_network/whole_data_set/intra_oral_microbiome/intra_oral_microbiome_lm_adjusted_cor"
#   )
#
#   oral_sample_cor <-
#     intra_oral_microbiome_lm_adjusted_cor
#
# }
#
#
####load the dims
####stool with skin
{
  load(
    "data_analysis/correlation_network/whole_data_set/stool_microbiome_vs_skin_microbiome/stool_skin_sample_wise_skin_dim"
  )
  load(
    "data_analysis/correlation_network/whole_data_set/stool_microbiome_vs_skin_microbiome/stool_skin_sample_wise_stool_dim"
  )
  load(
    "data_analysis/correlation_network/whole_data_set/stool_microbiome_vs_skin_microbiome/subject_wise/stool_skin_subject_wise_skin_dim"
  )
  load(
    "data_analysis/correlation_network/whole_data_set/stool_microbiome_vs_skin_microbiome/subject_wise/stool_skin_subject_wise_stool_dim"
  )
  
  ####stool with oral
  load(
    "data_analysis/correlation_network/whole_data_set/stool_microbiome_vs_oral_microbiome/stool_oral_sample_wise_stool_dim"
  )
  load(
    "data_analysis/correlation_network/whole_data_set/stool_microbiome_vs_oral_microbiome/stool_oral_sample_wise_oral_dim"
  )
  load(
    "data_analysis/correlation_network/whole_data_set/stool_microbiome_vs_oral_microbiome/subject_wise/stool_oral_subject_wise_oral_dim"
  )
  load(
    "data_analysis/correlation_network/whole_data_set/stool_microbiome_vs_oral_microbiome/subject_wise/stool_oral_subject_wise_stool_dim"
  )
  
  ####stool with nasal
  load(
    "data_analysis/correlation_network/whole_data_set/stool_microbiome_vs_nasal_microbiome/stool_nasal_sample_wise_stool_dim"
  )
  load(
    "data_analysis/correlation_network/whole_data_set/stool_microbiome_vs_nasal_microbiome/stool_nasal_sample_wise_nasal_dim"
  )
  load(
    "data_analysis/correlation_network/whole_data_set/stool_microbiome_vs_nasal_microbiome/subject_wise/stool_nasal_subject_wise_nasal_dim"
  )
  load(
    "data_analysis/correlation_network/whole_data_set/stool_microbiome_vs_nasal_microbiome/subject_wise/stool_nasal_subject_wise_stool_dim"
  )
  
  
  ####oral with nasal
  load(
    "data_analysis/correlation_network/whole_data_set/oral_microbiome_vs_nasal_microbiome/oral_nasal_sample_wise_oral_dim"
  )
  load(
    "data_analysis/correlation_network/whole_data_set/oral_microbiome_vs_nasal_microbiome/oral_nasal_sample_wise_nasal_dim"
  )
  load(
    "data_analysis/correlation_network/whole_data_set/oral_microbiome_vs_nasal_microbiome/subject_wise/oral_nasal_subject_wise_nasal_dim"
  )
  load(
    "data_analysis/correlation_network/whole_data_set/oral_microbiome_vs_nasal_microbiome/subject_wise/oral_nasal_subject_wise_oral_dim"
  )
  
  ####oral with skin
  load(
    "data_analysis/correlation_network/whole_data_set/oral_microbiome_vs_skin_microbiome/oral_skin_sample_wise_oral_dim"
  )
  load(
    "data_analysis/correlation_network/whole_data_set/oral_microbiome_vs_skin_microbiome/oral_skin_sample_wise_skin_dim"
  )
  load(
    "data_analysis/correlation_network/whole_data_set/oral_microbiome_vs_skin_microbiome/subject_wise/oral_skin_subject_wise_skin_dim"
  )
  load(
    "data_analysis/correlation_network/whole_data_set/oral_microbiome_vs_skin_microbiome/subject_wise/oral_skin_subject_wise_oral_dim"
  )
  
  ####oral with nasal
  load(
    "data_analysis/correlation_network/whole_data_set/oral_microbiome_vs_nasal_microbiome/oral_nasal_sample_wise_oral_dim"
  )
  load(
    "data_analysis/correlation_network/whole_data_set/oral_microbiome_vs_nasal_microbiome/oral_nasal_sample_wise_nasal_dim"
  )
  load(
    "data_analysis/correlation_network/whole_data_set/oral_microbiome_vs_nasal_microbiome/subject_wise/oral_nasal_subject_wise_nasal_dim"
  )
  load(
    "data_analysis/correlation_network/whole_data_set/oral_microbiome_vs_nasal_microbiome/subject_wise/oral_nasal_subject_wise_oral_dim"
  )
  
  ####skin with nasal
  load(
    "data_analysis/correlation_network/whole_data_set/skin_microbiome_vs_nasal_microbiome/skin_nasal_sample_wise_skin_dim"
  )
  load(
    "data_analysis/correlation_network/whole_data_set/skin_microbiome_vs_nasal_microbiome/skin_nasal_sample_wise_nasal_dim"
  )
  load(
    "data_analysis/correlation_network/whole_data_set/skin_microbiome_vs_nasal_microbiome/subject_wise/skin_nasal_subject_wise_nasal_dim"
  )
  load(
    "data_analysis/correlation_network/whole_data_set/skin_microbiome_vs_nasal_microbiome/subject_wise/skin_nasal_subject_wise_skin_dim"
  )
  
  
  #####intra-stool microbiome
  load(
    "data_analysis/correlation_network/whole_data_set/intra_stool_microbiome/intra_stool_sample_wise_dim"
  )
  load(
    "data_analysis/correlation_network/whole_data_set/intra_skin_microbiome/intra_skin_sample_wise_dim"
  )
  load(
    "data_analysis/correlation_network/whole_data_set/intra_nasal_microbiome/intra_nasal_sample_wise_dim"
  )
  load(
    "data_analysis/correlation_network/whole_data_set/intra_oral_microbiome/intra_oral_sample_wise_dim"
  )
  
  setwd(
    "data_analysis/correlation_network/whole_data_set/intra_microbiome_correlation"
  )
}
#
#
#
# ###---------------------------------------------------------------
# ######sample wise correlation
# stool_skin_sample_cor =
#   stool_skin_sample_cor %>%
#   dplyr::rename(from = microbiome, to = metabolite) %>%
#   dplyr::mutate(from_class = "stool", to_class = "skin") %>%
#   dplyr::mutate(from = paste("stool", from, sep = ""),
#                 to = paste("skin", to, sep = ""),)
#
# stool_nasal_sample_cor =
#   stool_nasal_sample_cor %>%
#   dplyr::rename(from = microbiome, to = metabolite) %>%
#   dplyr::mutate(from_class = "stool", to_class = "nasal") %>%
#   dplyr::mutate(from = paste("stool", from, sep = ""),
#                 to = paste("nasal", to, sep = ""),)
#
# stool_oral_sample_cor =
#   stool_oral_sample_cor %>%
#   dplyr::rename(from = microbiome, to = metabolite) %>%
#   dplyr::mutate(from_class = "stool", to_class = "oral") %>%
#   dplyr::mutate(from = paste("stool", from, sep = ""),
#                 to = paste("oral", to, sep = ""),)
#
# oral_nasal_sample_cor =
#   oral_nasal_sample_cor %>%
#   dplyr::rename(from = microbiome, to = metabolite) %>%
#   dplyr::mutate(from_class = "oral", to_class = "nasal") %>%
#   dplyr::mutate(from = paste("oral", from, sep = ""),
#                 to = paste("nasal", to, sep = ""),)
#
# oral_skin_sample_cor =
#   oral_skin_sample_cor %>%
#   dplyr::rename(from = microbiome, to = metabolite) %>%
#   dplyr::mutate(from_class = "oral", to_class = "skin") %>%
#   dplyr::mutate(from = paste("oral", from, sep = ""),
#                 to = paste("skin", to, sep = ""),)
#
# skin_nasal_sample_cor =
#   skin_nasal_sample_cor %>%
#   dplyr::rename(from = microbiome, to = metabolite) %>%
#   dplyr::mutate(from_class = "skin", to_class = "nasal") %>%
#   dplyr::mutate(from = paste("skin", from, sep = ""),
#                 to = paste("nasal", to, sep = ""),)
#
# stool_sample_cor <-
#   stool_sample_cor %>%
#   dplyr::rename(from = microbiome, to = metabolite) %>%
#   dplyr::mutate(from_class = "stool", to_class = "stool") %>%
#   dplyr::mutate(from = paste("stool", from, sep = ""),
#                 to = paste("stool", to, sep = ""),) %>%
#   dplyr::select(-name)
#
# skin_sample_cor <-
#   skin_sample_cor %>%
#   dplyr::rename(from = microbiome, to = metabolite) %>%
#   dplyr::mutate(from_class = "skin", to_class = "skin") %>%
#   dplyr::mutate(from = paste("skin", from, sep = ""),
#                 to = paste("skin", to, sep = ""),) %>%
#   dplyr::select(-name)
#
# nasal_sample_cor <-
#   nasal_sample_cor %>%
#   dplyr::rename(from = microbiome, to = metabolite) %>%
#   dplyr::mutate(from_class = "nasal", to_class = "nasal") %>%
#   dplyr::mutate(from = paste("nasal", from, sep = ""),
#                 to = paste("nasal", to, sep = ""),) %>%
#   dplyr::select(-name)
#
# oral_sample_cor <-
#   oral_sample_cor %>%
#   dplyr::rename(from = microbiome, to = metabolite) %>%
#   dplyr::mutate(from_class = "oral", to_class = "oral") %>%
#   dplyr::mutate(from = paste("oral", from, sep = ""),
#                 to = paste("oral", to, sep = ""),) %>%
#   dplyr::select(-name)

# #####sample wise
# sample_cor =
#   rbind(stool_skin_sample_cor,
#         stool_nasal_sample_cor,
#         stool_oral_sample_cor,
#         oral_nasal_sample_cor,
#         oral_skin_sample_cor,
#         skin_nasal_sample_cor,
#         stool_sample_cor,
#         skin_sample_cor,
#         nasal_sample_cor,
#         oral_sample_cor) %>%
#   dplyr::filter(p_adjust < 0.05)
#
# name = apply(sample_cor, 1, function(x){
#   paste(sort(as.character(x[1:2])), collapse = "_")
# })
#
# sample_cor =
#   sample_cor %>%
#   dplyr::mutate(name = name) %>%
#   dplyr::distinct(name, .keep_all = TRUE)
#
# save(sample_cor, file = "sample_cor")

load("sample_cor")

sum(sample_cor$from_class == sample_cor$to_class)
sum(sample_cor$from_class != sample_cor$to_class)

sample_cor$class_name =
  purrr::map(as.data.frame(t(sample_cor)), function(x) {
    paste(sort(as.character(x)[7:8]), collapse = "_")
  })  %>%
  unlist() %>%
  unname()

write.csv(sample_cor, file = "sample_cor.csv", row.names = FALSE)

sample_wise_number <-
  purrr::map(as.data.frame(t(sample_cor)), function(x) {
    c(sort(as.character(x)[7:8]), name = paste(sort(as.character(x)[7:8]), collapse = "_"))
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame() %>%
  dplyr::group_by(name) %>%
  dplyr::summarise(significant_n = n()) %>%
  dplyr::arrange(name) %>%
  dplyr::mutate(theoretical_n = NA)

sample_wise_number$theoretical_n[sample_wise_number$name == "nasal_nasal"] =
  (intra_nasal_sample_wise_dim[1] * (intra_nasal_sample_wise_dim[1] - 1)) /
  2

sample_wise_number$theoretical_n[sample_wise_number$name == "nasal_oral"] =
  oral_nasal_sample_wise_oral_dim[1] * oral_nasal_sample_wise_nasal_dim[1]

sample_wise_number$theoretical_n[sample_wise_number$name == "nasal_skin"] =
  skin_nasal_sample_wise_skin_dim[1] * skin_nasal_sample_wise_nasal_dim[1]

sample_wise_number$theoretical_n[sample_wise_number$name == "nasal_stool"] =
  stool_nasal_sample_wise_stool_dim[1] * stool_nasal_sample_wise_nasal_dim[1]

sample_wise_number$theoretical_n[sample_wise_number$name == "oral_oral"] =
  (intra_oral_sample_wise_dim[1] * (intra_oral_sample_wise_dim[1] - 1)) /
  2

sample_wise_number$theoretical_n[sample_wise_number$name == "oral_skin"] =
  oral_skin_sample_wise_oral_dim[1] * oral_skin_sample_wise_skin_dim[1]

sample_wise_number$theoretical_n[sample_wise_number$name == "oral_stool"] =
  stool_oral_sample_wise_oral_dim[1] * stool_oral_sample_wise_stool_dim[1]

sample_wise_number$theoretical_n[sample_wise_number$name == "skin_skin"] =
  (intra_skin_sample_wise_dim[1] * (intra_skin_sample_wise_dim[1] - 1)) /
  2

sample_wise_number$theoretical_n[sample_wise_number$name == "skin_stool"] =
  stool_skin_sample_wise_stool_dim[1] * stool_skin_sample_wise_skin_dim[1]

sample_wise_number$theoretical_n[sample_wise_number$name == "stool_stool"] =
  (intra_stool_sample_wise_dim[1] * (intra_stool_sample_wise_dim[1] - 1)) /
  2

sample_wise_number <-
  sample_wise_number %>%
  dplyr::mutate(percentage = significant_n * 100 / theoretical_n)

new_info =
  sample_wise_number$name %>%
  purrr::map(function(x) {
    significant_n_pos =
      sample_cor %>%
      dplyr::filter(class_name == x & cor > 0) %>%
      nrow()
    
    significant_n_neg =
      sample_cor %>%
      dplyr::filter(class_name == x & cor < 0) %>%
      nrow()
    
    data.frame(significant_n_pos = significant_n_pos,
               significant_n_neg = significant_n_neg)
    
  }) %>%
  dplyr::bind_rows()

sample_wise_number =
  cbind(sample_wise_number, new_info)

sample_wise_number =
  sample_wise_number %>%
  dplyr::mutate(
    percentage_pos = significant_n_pos * 100 / theoretical_n,
    percentage_neg = significant_n_neg * 100 / theoretical_n
  )

# ######subject wise
# subject_cor =
#   rbind(stool_skin_subject_cor,
#         stool_nasal_subject_cor,
#         stool_oral_subject_cor,
#         oral_nasal_subject_cor,
#         oral_skin_subject_cor,
#         skin_nasal_subject_cor) %>%
#   dplyr::filter(p_value_adjust < 0.05)
#
# name = apply(subject_cor, 1, function(x){
#   paste(sort(as.character(x[1:2])), collapse = "_")
# })
#
# subject_cor =
# subject_cor %>%
#   dplyr::mutate(name = name) %>%
#   dplyr::distinct(name, .keep_all = TRUE)
#
# save(subject_cor, file = "subject_cor")

load("subject_cor")

sum(subject_cor$from_class == subject_cor$to_class)
sum(subject_cor$from_class != subject_cor$to_class)

subject_cor$class_name =
  purrr::map(as.data.frame(t(subject_cor)), function(x) {
    paste(sort(as.character(x)[4:5]), collapse = "_")
  })  %>%
  unlist() %>%
  unname()

subject_wise_number =
  purrr::map(as.data.frame(t(subject_cor)), function(x) {
    c(sort(as.character(x)[4:5]), name = paste(sort(as.character(x)[4:5]), collapse = "_"))
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame() %>%
  dplyr::group_by(name) %>%
  dplyr::summarise(significant_n = n()) %>%
  dplyr::arrange(name) %>%
  dplyr::mutate(theoretical_n = NA)

subject_wise_number$theoretical_n[subject_wise_number$name == "nasal_nasal"] =
  (stool_nasal_subject_wise_nasal_dim[1] * (stool_nasal_subject_wise_nasal_dim[1] -
                                              1)) / 2

subject_wise_number$theoretical_n[subject_wise_number$name == "nasal_oral"] =
  oral_nasal_subject_wise_oral_dim[1] * oral_nasal_subject_wise_nasal_dim[1]

subject_wise_number$theoretical_n[subject_wise_number$name == "nasal_skin"] =
  skin_nasal_subject_wise_skin_dim[1] * skin_nasal_subject_wise_nasal_dim[1]

subject_wise_number$theoretical_n[subject_wise_number$name == "nasal_stool"] =
  stool_nasal_subject_wise_stool_dim[1] * stool_nasal_subject_wise_nasal_dim[1]

subject_wise_number$theoretical_n[subject_wise_number$name == "oral_oral"] =
  (stool_oral_subject_wise_oral_dim[1] * (stool_oral_subject_wise_oral_dim[1] -
                                            1)) / 2

subject_wise_number$theoretical_n[subject_wise_number$name == "oral_skin"] =
  oral_skin_subject_wise_oral_dim[1] * oral_skin_subject_wise_skin_dim[1]

subject_wise_number$theoretical_n[subject_wise_number$name == "oral_stool"] =
  stool_oral_subject_wise_oral_dim[1] * stool_oral_subject_wise_stool_dim[1]

subject_wise_number$theoretical_n[subject_wise_number$name == "skin_skin"] =
  (stool_skin_subject_wise_skin_dim[1] * (stool_skin_subject_wise_skin_dim[1] - 1)) /
  2

subject_wise_number$theoretical_n[subject_wise_number$name == "skin_stool"] =
  stool_skin_subject_wise_stool_dim[1] * stool_skin_subject_wise_skin_dim[1]

subject_wise_number$theoretical_n[subject_wise_number$name == "stool_stool"] =
  (stool_skin_subject_wise_stool_dim[1] * (stool_skin_subject_wise_stool_dim[1] - 1)) /
  2

subject_wise_number <-
  subject_wise_number %>%
  dplyr::mutate(percentage = significant_n * 100 / theoretical_n)


new_info =
  subject_wise_number$name %>%
  purrr::map(function(x) {
    significant_n_pos =
      sample_cor %>%
      dplyr::filter(class_name == x & cor > 0) %>%
      nrow()
    
    significant_n_neg =
      sample_cor %>%
      dplyr::filter(class_name == x & cor < 0) %>%
      nrow()
    
    data.frame(significant_n_pos = significant_n_pos,
               significant_n_neg = significant_n_neg)
    
  }) %>%
  dplyr::bind_rows()

subject_wise_number =
  cbind(subject_wise_number, new_info)

subject_wise_number =
  subject_wise_number %>%
  dplyr::mutate(
    percentage_pos = significant_n_pos * 100 / theoretical_n,
    percentage_neg = significant_n_neg * 100 / theoretical_n
  )

load(here::here(
  "data_analysis/stool_microbiome/data_preparation/variable_info"
))
stool_variable_info = variable_info

load(here::here(
  "data_analysis/skin_microbiome/data_preparation/variable_info"
))
skin_variable_info = variable_info

load(here::here(
  "data_analysis/nasal_microbiome/data_preparation/variable_info"
))
nasal_variable_info = variable_info

load(here::here(
  "data_analysis/oral_microbiome/data_preparation/variable_info"
))
oral_variable_info = variable_info


######phylum class distributation
###stool
dim(stool_variable_info)

from <-
  sample_cor %>%
  dplyr::filter(class_name == "stool_stool") %>%
  pull(from) %>%
  stringr::str_replace("stool", "") %>%
  unique()

to <-
  sample_cor %>%
  dplyr::filter(class_name == "stool_stool") %>%
  pull(to) %>%
  stringr::str_replace("stool", "") %>%
  unique()

temp <-
  stool_variable_info %>%
  dplyr::filter(variable_id %in% c(from, to))

plot <-
  temp %>%
  dplyr::count(Phylum) %>%
  ggplot(aes(x = "", y = n, fill = Phylum)) +
  geom_bar(stat = "identity",
           width = 1,
           color = "black") +
  scale_fill_manual(values = phylum_color) +
  coord_polar("y", start = 0) +
  theme_void()
plot
# ggsave(plot, filename = "stool_phylum_distributation.pdf", width = 7, height = 7)



###skin
dim(skin_variable_info)

from <-
  sample_cor %>%
  dplyr::filter(class_name == "skin_skin") %>%
  pull(from) %>%
  stringr::str_replace("skin", "") %>%
  unique()

to <-
  sample_cor %>%
  dplyr::filter(class_name == "skin_skin") %>%
  pull(to) %>%
  stringr::str_replace("skin", "") %>%
  unique()

temp <-
  skin_variable_info %>%
  dplyr::filter(variable_id %in% c(from, to))

plot <-
  temp %>%
  dplyr::count(Phylum) %>%
  ggplot(aes(x = "", y = n, fill = Phylum)) +
  geom_bar(stat = "identity",
           width = 1,
           color = "black") +
  scale_fill_manual(values = phylum_color) +
  coord_polar("y", start = 0) +
  theme_void()
plot
# ggsave(plot, filename = "skin_phylum_distributation.pdf", width = 7, height = 7)







###oral
dim(oral_variable_info)

from <-
  sample_cor %>%
  dplyr::filter(class_name == "oral_oral") %>%
  pull(from) %>%
  stringr::str_replace("oral", "") %>%
  unique()

to <-
  sample_cor %>%
  dplyr::filter(class_name == "oral_oral") %>%
  pull(to) %>%
  stringr::str_replace("oral", "") %>%
  unique()

temp <-
  oral_variable_info %>%
  dplyr::filter(variable_id %in% c(from, to))

plot <-
  temp %>%
  dplyr::count(Phylum) %>%
  ggplot(aes(x = "", y = n, fill = Phylum)) +
  geom_bar(stat = "identity",
           width = 1,
           color = "black") +
  scale_fill_manual(values = phylum_color) +
  coord_polar("y", start = 0) +
  theme_void()
plot
# ggsave(plot, filename = "oral_phylum_distributation.pdf", width = 7, height = 7)




###nasal
dim(nasal_variable_info)

from <-
  sample_cor %>%
  dplyr::filter(class_name == "nasal_nasal") %>%
  pull(from) %>%
  stringr::str_replace("nasal", "") %>%
  unique()

to <-
  sample_cor %>%
  dplyr::filter(class_name == "nasal_nasal") %>%
  pull(to) %>%
  stringr::str_replace("nasal", "") %>%
  unique()

temp <-
  nasal_variable_info %>%
  dplyr::filter(variable_id %in% c(from, to))

plot <-
  temp %>%
  dplyr::count(Phylum) %>%
  ggplot(aes(x = "", y = n, fill = Phylum)) +
  geom_bar(stat = "identity",
           width = 1,
           color = "black") +
  scale_fill_manual(values = phylum_color) +
  coord_polar("y", start = 0) +
  theme_void()
plot
# ggsave(plot, filename = "nasal_phylum_distributation.pdf", width = 7, height = 7)



######ggraph
#####sample wise
library(ggraph)
sample_wise_number

edge_data =
  sample_wise_number %>%
  dplyr::select(name, percentage) %>%
  tidyr::separate(col = name,
                  sep = "_",
                  into = c("from", "to")) %>%
  dplyr::mutate(from = stringr::str_to_title(from)) %>%
  dplyr::mutate(to = stringr::str_to_title(to))

node_data =
  data.frame(
    node = c("stool", "skin", "nasal", "oral"),
    true_name = c("stool", "skin", "nasal", "oral")
  ) %>%
  dplyr::mutate(node = stringr::str_to_title(node)) %>%
  dplyr::mutate(true_name = stringr::str_to_title(true_name))

#######within network
library(ggraph)
library(tidygraph)
library(igraph)

edge_data1 =
  edge_data %>%
  dplyr::filter(from == to) %>%
  dplyr::mutate(from = paste(from, "from", sep = '_')) %>%
  dplyr::mutate(to = paste(to, "to", sep = '_'))

node_data1 =
  rbind(node_data %>% dplyr::mutate(node = paste(node, "from", sep = "_")),
        node_data %>% dplyr::mutate(node = paste(node, "to", sep = "_")))

temp_graph <-
  tidygraph::tbl_graph(nodes = node_data1,
                       edges = edge_data1,
                       directed = FALSE)
#####up-down
g <- temp_graph
library(igraph)
library(ggraph)

V(g)$type <- bipartite_mapping(g)$type

coords <-
  create_layout(g, layout = "bipartite") %>%
  dplyr::select(node, true_name, x, y)

coords$x[grep("from$", coords$node)] = 0
coords$x[grep("to$", coords$node)] = 1

coords$y[coords$true_name == "Stool"] = 1
coords$y[coords$true_name == "Skin"] = 3
coords$y[coords$true_name == "Nasal"] = 2
coords$y[coords$true_name == "Oral"] = 4

my_graph <-
  create_layout(
    graph = g,
    layout = "manual",
    x = coords$x,
    y = coords$y
    # node.position = coords
  )

plot =
  ggraph(my_graph,
         layout = 'bipartite') +
  geom_edge_link(
    aes(edge_width = percentage,
        label = paste(round(percentage, 2), "%", sep = "")),
    angle_calc = 'along',
    label_dodge = unit(5, 'mm'),
    show.legend = TRUE
  ) +
  ggraph::geom_node_point(
    aes(fill = true_name),
    shape = 21,
    size = 8,
    show.legend = TRUE
  ) +
  guides(fill = guide_legend(title = "Body site"),
         size = FALSE) +
  # geom_node_text(aes(label = paste0(round(intra_number, 2), "%"))) +
  scale_fill_manual(values = body_site_color) +
  scale_edge_width(limits = c(1, 100), range = c(1, 20)) +
  ggraph::theme_graph() +
  theme(legend.position = "bottom")

plot

sample_cor2 =
  sample_cor %>%
  dplyr::mutate(
    from1 = to,
    to1 = from,
    from_class1 = to_class,
    to_class1 = from_class
  ) %>%
  dplyr::select(-c(from, to, from_class, to_class)) %>%
  dplyr::rename(
    from = from1,
    to = to1,
    from_class = from_class1,
    to_class = to_class1
  ) %>%
  dplyr::select(colnames(sample_cor))

sample_cor2$class_name =
  sample_cor2$class_name %>%
  stringr::str_split("_") %>%
  purrr::map(function(x) {
    paste(x[2:1], collapse = "_")
  }) %>%
  unlist()

temp_sample_cor = rbind(sample_cor,
                        sample_cor2)

result <-
  purrr::map(unique(temp_sample_cor$class_name), function(x) {
    temp_data =
      temp_sample_cor %>%
      dplyr::filter(class_name == x)
    
    from_class = unique(temp_data$from_class)
    to_class = unique(temp_data$to_class)
    
    from_id = unique(stringr::str_replace(temp_data$from, from_class, ""))
    to_id = unique(stringr::str_replace(temp_data$to, to_class, ""))
    
    if (from_class == "stool") {
      from_phylum = stool_variable_info$Phylum[match(from_id, stool_variable_info$variable_id)]
    }
    
    if (from_class == "skin") {
      from_phylum = skin_variable_info$Phylum[match(from_id, skin_variable_info$variable_id)]
    }
    
    if (from_class == "nasal") {
      from_phylum = nasal_variable_info$Phylum[match(from_id, nasal_variable_info$variable_id)]
    }
    
    if (from_class == "oral") {
      from_phylum = oral_variable_info$Phylum[match(from_id, oral_variable_info$variable_id)]
    }
    
    
    if (to_class == "stool") {
      to_phylum = stool_variable_info$Phylum[match(to_id, stool_variable_info$variable_id)]
    }
    
    if (to_class == "skin") {
      to_phylum = skin_variable_info$Phylum[match(to_id, skin_variable_info$variable_id)]
    }
    
    if (to_class == "nasal") {
      to_phylum = nasal_variable_info$Phylum[match(to_id, nasal_variable_info$variable_id)]
    }
    
    if (to_class == "oral") {
      to_phylum = oral_variable_info$Phylum[match(to_id, oral_variable_info$variable_id)]
    }
    
    
    rbind(
      data.frame(phylum = from_phylum) %>%
        dplyr::count(phylum) %>%
        dplyr::mutate(percentage = n * 100 / sum(n)) %>%
        dplyr::mutate(body_site = from_class) %>%
        dplyr::mutate(
          class = "from",
          from_class = from_class,
          to_class = to_class
        ),
      data.frame(phylum = to_phylum) %>%
        dplyr::count(phylum) %>%
        dplyr::mutate(percentage = n * 100 / sum(n)) %>%
        dplyr::mutate(body_site = to_class) %>%
        dplyr::mutate(
          class = "to",
          from_class = from_class,
          to_class = to_class
        )
    )
  }) %>%
  dplyr::bind_rows()

result =
  result %>%
  dplyr::mutate(
    phylum = case_when(
      phylum == "Actinobacteria" ~ "Actinobacteria",
      phylum == "Bacteroidetes" ~ "Bacteroidetes",
      phylum == "Firmicutes" ~ "Firmicutes",
      phylum == "Proteobacteria" ~ "Proteobacteria",
      TRUE ~ 'Other'
    )
  )

library(plyr)

result =
  result %>%
  plyr::dlply(.variables = .(phylum, body_site, class, from_class, to_class)) %>%
  purrr::map(function(x) {
    if (nrow(x) == 1) {
      return(x)
    } else{
      x$n = sum(x$n)
      x$percentage = sum(x$percentage)
      x =
        x %>%
        dplyr::distinct(phylum, .keep_all = TRUE)
      return(x)
    }
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()

rownames(result) = NULL

library(scatterpie)

phylum_name = unique(result$phylum)

temp_result =
  result %>%
  dplyr::filter(class == "from") %>%
  dplyr::filter(from_class == to_class)

temp_data =
  temp_result %>%
  tidyr::pivot_wider(names_from = phylum, values_from = percentage)

temp_data =
  temp_data %>%
  dplyr::mutate(x = 1) %>%
  dplyr::mutate(
    y = case_when(
      to_class == "stool" ~ 1,
      to_class == "nasal" ~ 2,
      to_class == "skin" ~ 3,
      to_class == "oral" ~ 4,
    )
  )

temp_data[is.na(temp_data)] = 0

plot1 =
  ggplot() +
  geom_scatterpie(
    aes(x = x, y = y, r = 0.3),
    cols = c(
      "Actinobacteria",
      "Bacteroidetes",
      "Firmicutes",
      "Proteobacteria",
      "Other"
    ),
    data = temp_data
  ) +
  scale_fill_manual(values = c(phylum_color, "Other" = "grey")) +
  coord_equal() +
  base_theme +
  labs(x = "", y = "") +
  scale_y_continuous(breaks = c(1, 2, 3, 4),
                     labels = c("Stool", "Nasal", "Skin", "Oral")) +
  theme(
    legend.position = "left",
    axis.ticks = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_blank()
  )

plot1

library(patchwork)

within_plot =
  plot1 + plot + plot_layout(ncol = 2)
within_plot
# ggsave(plot = within_plot, filename = "within_plot.pdf", width = 7, height = 7)

#######between network
library(ggraph)
library(tidygraph)
library(igraph)

edge_data1 =
  edge_data %>%
  dplyr::filter(from != to) %>%
  dplyr::mutate(from = paste(from, "from", sep = '_')) %>%
  dplyr::mutate(to = paste(to, "to", sep = '_')) %>%
  dplyr::arrange(desc(percentage))

edge_data1$from = masstools::name_duplicated(edge_data1$from)
edge_data1$to = masstools::name_duplicated(edge_data1$to)

node_data1 =
  data.frame(node = unique(c(edge_data1$from, edge_data1$to))) %>%
  dplyr::mutate(true_name = stringr::str_extract(node, "Stool|Skin|Nasal|Oral"))

temp_graph <-
  tidygraph::tbl_graph(nodes = node_data1,
                       edges = edge_data1,
                       directed = FALSE)
#####up-down
g <- temp_graph
library(igraph)
library(ggraph)

V(g)$type <- bipartite_mapping(g)$type

coords <-
  create_layout(g, layout = "bipartite") %>%
  dplyr::select(node, true_name, x, y)

coords$x[grep("from", coords$node)] = 0
coords$x[grep("to", coords$node)] = 1

coords$y[coords$node == "Nasal_from_1"] = 6
coords$y[coords$node == "Skin_to_1"] = 6

coords$y[coords$node == "Nasal_from_2"] = 5
coords$y[coords$node == "Stool_to_1"] = 5

coords$y[coords$node == "Skin_from"] = 4
coords$y[coords$node == "Stool_to_2"] = 4

coords$y[coords$node == "Nasal_from_3"] = 3
coords$y[coords$node == "Oral_to"] = 3

coords$y[coords$node == "Oral_from_1"] = 2
coords$y[coords$node == "Skin_to_2"] = 2

coords$y[coords$node == "Oral_from_2"] = 1
coords$y[coords$node == "Stool_to_3"] = 1

my_graph <-
  create_layout(
    graph = g,
    layout = "manual",
    x = coords$x,
    y = coords$y
    # node.position = coords
  )

plot =
  ggraph(my_graph,
         layout = 'bipartite') +
  geom_edge_link(
    aes(edge_width = percentage,
        label = paste(round(percentage, 2), "%", sep = "")),
    angle_calc = 'along',
    label_dodge = unit(5, 'mm'),
    show.legend = TRUE
  ) +
  ggraph::geom_node_point(
    aes(fill = true_name),
    shape = 21,
    size = 8,
    show.legend = TRUE
  ) +
  guides(fill = guide_legend(title = "Body site"),
         size = FALSE) +
  scale_fill_manual(values = body_site_color) +
  scale_edge_width(limits = c(1, 100), range = c(1, 20)) +
  ggraph::theme_graph() +
  theme(legend.position = "bottom")

plot

sample_cor2 =
  sample_cor %>%
  dplyr::mutate(
    from1 = to,
    to1 = from,
    from_class1 = to_class,
    to_class1 = from_class
  ) %>%
  dplyr::select(-c(from, to, from_class, to_class)) %>%
  dplyr::rename(
    from = from1,
    to = to1,
    from_class = from_class1,
    to_class = to_class1
  ) %>%
  dplyr::select(colnames(sample_cor))

sample_cor2$class_name =
  sample_cor2$class_name %>%
  stringr::str_split("_") %>%
  purrr::map(function(x) {
    paste(x[2:1], collapse = "_")
  }) %>%
  unlist()

temp_sample_cor = rbind(sample_cor,
                        sample_cor2)

result =
  purrr::map(unique(temp_sample_cor$class_name), function(x) {
    temp_data =
      temp_sample_cor %>%
      dplyr::filter(class_name == x)
    
    from_class = unique(temp_data$from_class)
    to_class = unique(temp_data$to_class)
    
    from_id = unique(stringr::str_replace(temp_data$from, from_class, ""))
    to_id = unique(stringr::str_replace(temp_data$to, to_class, ""))
    
    if (from_class == "stool") {
      from_phylum = stool_variable_info$Phylum[match(from_id, stool_variable_info$variable_id)]
    }
    
    if (from_class == "skin") {
      from_phylum = skin_variable_info$Phylum[match(from_id, skin_variable_info$variable_id)]
    }
    
    if (from_class == "nasal") {
      from_phylum = nasal_variable_info$Phylum[match(from_id, nasal_variable_info$variable_id)]
    }
    
    if (from_class == "oral") {
      from_phylum = oral_variable_info$Phylum[match(from_id, oral_variable_info$variable_id)]
    }
    
    
    if (to_class == "stool") {
      to_phylum = stool_variable_info$Phylum[match(to_id, stool_variable_info$variable_id)]
    }
    
    if (to_class == "skin") {
      to_phylum = skin_variable_info$Phylum[match(to_id, skin_variable_info$variable_id)]
    }
    
    if (to_class == "nasal") {
      to_phylum = nasal_variable_info$Phylum[match(to_id, nasal_variable_info$variable_id)]
    }
    
    if (to_class == "oral") {
      to_phylum = oral_variable_info$Phylum[match(to_id, oral_variable_info$variable_id)]
    }
    
    
    rbind(
      data.frame(phylum = from_phylum) %>%
        dplyr::count(phylum) %>%
        dplyr::mutate(percentage = n * 100 / sum(n)) %>%
        dplyr::mutate(body_site = from_class) %>%
        dplyr::mutate(
          class = "from",
          from_class = from_class,
          to_class = to_class
        ),
      data.frame(phylum = to_phylum) %>%
        dplyr::count(phylum) %>%
        dplyr::mutate(percentage = n * 100 / sum(n)) %>%
        dplyr::mutate(body_site = to_class) %>%
        dplyr::mutate(
          class = "to",
          from_class = from_class,
          to_class = to_class
        )
    )
  }) %>%
  dplyr::bind_rows()


result =
  result %>%
  dplyr::mutate(
    phylum = case_when(
      phylum == "Actinobacteria" ~ "Actinobacteria",
      phylum == "Bacteroidetes" ~ "Bacteroidetes",
      phylum == "Firmicutes" ~ "Firmicutes",
      phylum == "Proteobacteria" ~ "Proteobacteria",
      TRUE ~ 'Other'
    )
  )


library(plyr)

result =
  result %>%
  plyr::dlply(.variables = .(phylum, body_site, class, from_class, to_class)) %>%
  purrr::map(function(x) {
    if (nrow(x) == 1) {
      return(x)
    } else{
      x$n = sum(x$n)
      x$percentage = sum(x$percentage)
      x =
        x %>%
        dplyr::distinct(phylum, .keep_all = TRUE)
      return(x)
    }
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()

rownames(result) = NULL

library(scatterpie)

phylum_name = unique(result$phylum)

temp_result =
  result %>%
  dplyr::filter(class == "from") %>%
  dplyr::filter(
    from_class == "nasal" & to_class == "skin" |
      from_class == "nasal" & to_class == "stool" |
      from_class == "skin" & to_class == "stool" |
      from_class == "nasal" & to_class == "oral" |
      from_class == "oral" & to_class == "skin" |
      from_class == "oral" & to_class == "stool"
  )

temp_data =
  temp_result %>%
  tidyr::pivot_wider(names_from = phylum, values_from = percentage)

temp_data =
  temp_data %>%
  dplyr::mutate(x = 1) %>%
  dplyr::mutate(
    y = case_when(
      from_class == "nasal" & to_class == "skin" ~ 6,
      from_class == "nasal" & to_class == "stool" ~ 5,
      from_class == "skin" & to_class == "stool" ~ 4,
      from_class == "nasal" & to_class == "oral" ~ 3,
      from_class == "oral" & to_class == "skin" ~ 2,
      from_class == "oral" & to_class == "stool" ~ 1
    )
  )

temp_data[is.na(temp_data)] = 0

plot1 =
  ggplot() +
  geom_scatterpie(
    aes(x = x, y = y, r = 0.3),
    cols = c(
      "Actinobacteria",
      "Bacteroidetes",
      "Firmicutes",
      "Proteobacteria",
      "Other"
    ),
    data = temp_data
  ) +
  scale_fill_manual(values = c(phylum_color, "Other" = "grey")) +
  coord_equal() +
  base_theme +
  labs(x = "", y = "") +
  scale_y_continuous(
    breaks = c(1, 2, 3, 4, 5, 6),
    labels = c("Oral", "Oral", "Nasal", "Skin", "Nasal", "Nasal")
  ) +
  theme(
    legend.position = "left",
    axis.ticks = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_blank()
  )

plot1

temp_result =
  result %>%
  dplyr::filter(class == "to") %>%
  dplyr::filter(
    from_class == "nasal" & to_class == "skin" |
      from_class == "nasal" & to_class == "stool" |
      from_class == "skin" & to_class == "stool" |
      from_class == "nasal" & to_class == "oral" |
      from_class == "oral" & to_class == "skin" |
      from_class == "oral" & to_class == "stool"
  )

temp_data =
  temp_result %>%
  tidyr::pivot_wider(names_from = phylum, values_from = percentage)

temp_data =
  temp_data %>%
  dplyr::mutate(x = 1) %>%
  dplyr::mutate(
    y = case_when(
      from_class == "nasal" & to_class == "skin" ~ 6,
      from_class == "nasal" & to_class == "stool" ~ 5,
      from_class == "skin" & to_class == "stool" ~ 4,
      from_class == "nasal" & to_class == "oral" ~ 3,
      from_class == "oral" & to_class == "skin" ~ 2,
      from_class == "oral" & to_class == "stool" ~ 1
    )
  )

temp_data[is.na(temp_data)] = 0

plot2 =
  ggplot() +
  geom_scatterpie(
    aes(x = x, y = y, r = 0.3),
    cols = c(
      "Actinobacteria",
      "Bacteroidetes",
      "Firmicutes",
      "Proteobacteria",
      "Other"
    ),
    data = temp_data
  ) +
  scale_fill_manual(values = c(phylum_color, "Other" = "grey")) +
  coord_equal() +
  base_theme +
  labs(x = "", y = "") +
  scale_y_continuous(
    breaks = c(1, 2, 3, 4, 5, 6),
    labels = c("Stool", "Skin", "Oral", "Stool", "Stool", "Skin")
  ) +
  theme(
    legend.position = "right",
    axis.ticks = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_blank()
  )

plot2

library(patchwork)

between_plot =
  plot1 + plot + plot2 + plot_layout(ncol = 3)
between_plot
# ggsave(plot = between_plot, filename = "between_plot.pdf", width = 11, height = 10.5)





# #####high abundance geneara
# stool_id = c("Bacteroides", "Prevotella", "Phocaeicola", "Unclassified_Ruminococcaceae")
# load(here::here("data_analysis/stool_microbiome/data_preparation/variable_info"))
# stool_variable_info = variable_info
#
# stool_id =
#   paste0("stool",stool_variable_info$variable_id[match(stool_id, stool_variable_info$Genus)])
#
# skin_id = c("Cutibacterium", "Staphylococcus",
#             "Corynebacterium", "Anaerococcus")
# load(here::here("data_analysis/skin_microbiome/data_preparation/variable_info"))
# skin_variable_info = variable_info
#
# skin_id =
#   paste0("skin",skin_variable_info$variable_id[match(skin_id, skin_variable_info$Genus)])
#
#
# sample_cor %>%
#   dplyr::filter(from %in% stool_id | to %in% stool_id)
