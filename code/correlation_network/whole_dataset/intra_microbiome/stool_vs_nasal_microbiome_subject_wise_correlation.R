#' ---
#' title: "Stool microbiome nasal_microbiome correlation"
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
setwd("data_analysis/correlation_network/whole_data_set/stool_microbiome_vs_nasal_microbiome")


# ####load data
# ###stool microbiome
# {
#   load(here::here(
#     "data_analysis/stool_microbiome/data_preparation/expression_data"
#   ))
#   load(here::here(
#     "data_analysis/stool_microbiome/data_preparation/sample_info"
#   ))
#   load(here::here(
#     "data_analysis/stool_microbiome/data_preparation/variable_info"
#   ))
# }
# 
# stool_microbiome_expression_data = expression_data
# stool_microbiome_sample_info = sample_info
# stool_microbiome_variable_info = variable_info
# 
# ###read genus table
# expression_data =
#   data.table::fread(here::here("data/from_xin/Genus Table/ST/Genus_ST.csv")) %>%
#   as.data.frame() %>%
#   tibble::column_to_rownames(var = "SampleID") %>%
#   dplyr::select(-c(V1:batch)) %>%
#   t() %>%
#   as.data.frame()
# 
# stool_microbiome_variable_info =
#   stool_microbiome_variable_info[match(rownames(expression_data), stool_microbiome_variable_info$Genus),]
# 
# stool_microbiome_variable_info$Genus == rownames(expression_data)
# 
# rownames(expression_data) = stool_microbiome_variable_info$variable_id
# 
# stool_microbiome_expression_data =
#   expression_data
# 
# stool_microbiome_variable_info =
#   stool_microbiome_variable_info %>%
#   dplyr::filter(!stringr::str_detect(Genus, "Unclassified_Bacteria"))
# 
# stool_microbiome_expression_data =
#   stool_microbiome_expression_data[stool_microbiome_variable_info$variable_id,]
# 
# dim(stool_microbiome_sample_info)
# dim(stool_microbiome_variable_info)
# 
# rownames(stool_microbiome_expression_data) == stool_microbiome_variable_info$variable_id
# colnames(stool_microbiome_expression_data) == stool_microbiome_sample_info$sample_id
# 
# ##nasal_microbiome
# ###nasal microbiome
# {
#   load(here::here("data_analysis/nasal_microbiome/data_preparation/expression_data"))
#   load(here::here("data_analysis/nasal_microbiome/data_preparation/sample_info"))
#   load(here::here("data_analysis/nasal_microbiome/data_preparation/variable_info"))
# }
# 
# nasal_microbiome_expression_data = expression_data
# nasal_microbiome_sample_info = sample_info
# nasal_microbiome_variable_info = variable_info
# 
# ###read genus table
# expression_data =
#   data.table::fread(here::here("data/from_xin/Genus Table/NS/Genus_NS.csv")) %>%
#   as.data.frame() %>%
#   tibble::column_to_rownames(var = "SampleID") %>%
#   dplyr::select(-c(V1:batch)) %>%
#   t() %>%
#   as.data.frame()
# 
# nasal_microbiome_variable_info =
#   nasal_microbiome_variable_info[match(rownames(expression_data), nasal_microbiome_variable_info$Genus),]
# 
# nasal_microbiome_variable_info$Genus == rownames(expression_data)
# 
# ###remove the variables which Genus are NA
# remove_idx = which(is.na(nasal_microbiome_variable_info$Genus))
# remove_idx
# if(length(remove_idx) > 0){
#   nasal_microbiome_variable_info = nasal_microbiome_variable_info[-remove_idx,]
#   expression_data = expression_data[-remove_idx,]
# }
# 
# rownames(expression_data) = nasal_microbiome_variable_info$variable_id
# 
# nasal_microbiome_expression_data =
#   expression_data
# 
# nasal_microbiome_variable_info =
#   nasal_microbiome_variable_info %>%
#   dplyr::filter(!stringr::str_detect(Genus, "Unclassified_Bacteria"))
# 
# nasal_microbiome_expression_data =
#   nasal_microbiome_expression_data[nasal_microbiome_variable_info$variable_id,]
# 
# dim(nasal_microbiome_sample_info)
# dim(nasal_microbiome_variable_info)
# dim(nasal_microbiome_expression_data)
# 
# rownames(nasal_microbiome_expression_data) == nasal_microbiome_variable_info$variable_id
# colnames(nasal_microbiome_expression_data) == nasal_microbiome_sample_info$sample_id
# 
# ###match samples
# dim(stool_microbiome_sample_info)
# dim(nasal_microbiome_sample_info)
# 
# length(unique(stool_microbiome_sample_info$subject_id))
# 
# ###just matched samples according to sample id, only 1 missed
# intersect_sample_id =
#   intersect(stool_microbiome_sample_info$sample_id,
#             nasal_microbiome_sample_info$sample_id)
# 
# length(intersect_sample_id)
# 
# stool_microbiome_expression_data =
#   stool_microbiome_expression_data[,intersect_sample_id]
# 
# nasal_microbiome_expression_data =
#   nasal_microbiome_expression_data[,intersect_sample_id]
# 
# stool_microbiome_sample_info =
#   stool_microbiome_sample_info[match(intersect_sample_id, stool_microbiome_sample_info$sample_id),]
# 
# nasal_microbiome_sample_info =
#   nasal_microbiome_sample_info[match(intersect_sample_id, nasal_microbiome_sample_info$sample_id),]
# 
# length(unique(nasal_microbiome_sample_info$subject_id))
# 
# ###only remain the subjects with at least >= 5
# remian_subject_id =
#   stool_microbiome_sample_info %>%
#   dplyr::group_by(subject_id) %>%
#   dplyr::summarise(n = n()) %>%
#   dplyr::ungroup() %>%
#   dplyr::filter(n > 5) %>%
#   dplyr::pull(subject_id)
# 
# nasal_microbiome_sample_info =
#   nasal_microbiome_sample_info %>%
#   dplyr::filter(subject_id %in% remian_subject_id)
# 
# nasal_microbiome_expression_data =
#   nasal_microbiome_expression_data[, nasal_microbiome_sample_info$sample_id]
# 
# stool_microbiome_sample_info =
#   stool_microbiome_sample_info %>%
#   dplyr::filter(subject_id %in% remian_subject_id)
# 
# stool_microbiome_expression_data =
#   stool_microbiome_expression_data[,stool_microbiome_sample_info$sample_id]
# 
# ##only remain the genus at least in 10% subjects
# remain_idx =
#   which(rowSums(stool_microbiome_expression_data) > 0)
# 
# stool_microbiome_expression_data = stool_microbiome_expression_data[remain_idx,]
# stool_microbiome_variable_info = stool_microbiome_variable_info[remain_idx,,drop = FALSE]
# 
# remain_idx =
#   stool_microbiome_expression_data %>%
#   apply(1, function(x){
#     sum(as.numeric(x) == 0) / ncol(stool_microbiome_expression_data)
#   }) %>%
#   `<`(0.9) %>%
#   which()
# 
# length(remain_idx)
# 
# stool_microbiome_expression_data = stool_microbiome_expression_data[remain_idx,]
# stool_microbiome_variable_info = stool_microbiome_variable_info[remain_idx,,drop = FALSE]
# 
# ##only remain the genus at least in 10% subjects
# remain_idx =
#   which(rowSums(nasal_microbiome_expression_data) > 0)
# 
# nasal_microbiome_expression_data = nasal_microbiome_expression_data[remain_idx,]
# nasal_microbiome_variable_info = nasal_microbiome_variable_info[remain_idx,,drop = FALSE]
# 
# remain_idx =
#   nasal_microbiome_expression_data %>%
#   apply(1, function(x){
#     sum(as.numeric(x) == 0) / ncol(nasal_microbiome_expression_data)
#   }) %>%
#   `<`(0.9) %>%
#   which()
# 
# length(remain_idx)
# 
# nasal_microbiome_expression_data = nasal_microbiome_expression_data[remain_idx,]
# nasal_microbiome_variable_info = nasal_microbiome_variable_info[remain_idx,,drop = FALSE]
# 
# ##save data
# {
#   save(stool_microbiome_expression_data, file = "stool_microbiome_expression_data")
#   save(stool_microbiome_variable_info, file = "stool_microbiome_variable_info")
#   save(stool_microbiome_sample_info, file = "stool_microbiome_sample_info")
# 
#   save(nasal_microbiome_expression_data, file = "nasal_microbiome_expression_data")
#   save(nasal_microbiome_variable_info, file = "nasal_microbiome_variable_info")
#   save(nasal_microbiome_sample_info, file = "nasal_microbiome_sample_info")
# }
# 
{
  load("stool_microbiome_expression_data")
  load("stool_microbiome_variable_info")
  load("stool_microbiome_sample_info")

  load("nasal_microbiome_expression_data")
  load("nasal_microbiome_variable_info")
  load("nasal_microbiome_sample_info")
}

dim(stool_microbiome_expression_data)
dim(nasal_microbiome_expression_data)

#####set to other folder
setwd("subject_wise")

###finally, for stool microbiome, 106 genus, for nasal_microbiome, 76 genus

######--------------------------------------------------------------------------
library(plyr)

#####
dim(nasal_microbiome_expression_data)
dim(stool_microbiome_expression_data)

stool_microbiome_sample_info$subject_id == nasal_microbiome_sample_info$subject_id

##https://rpubs.com/DKCH2020/578881
##https://ourcodingclub.github.io/tutorials/mixed-models/
##https://zhuanlan.zhihu.com/p/63092231
##https://www.linglab.cn/knowledge/10

#######for each subject, just get the mean value for each genera
stool_microbiome_sample_info %>% 
  dplyr::count(subject_id)

library(plyr)

stool_data = 
  stool_microbiome_sample_info %>%
  plyr::dlply(.variables = .(subject_id)) %>% 
  purrr::map(function(x){
    stool_microbiome_expression_data[,x$sample_id] %>% 
      apply(1, mean)
  }) %>% 
  dplyr::bind_cols() %>% 
  as.data.frame()
  
rownames(stool_data) = rownames(stool_microbiome_expression_data)


nasal_data = 
  nasal_microbiome_sample_info %>%
  plyr::dlply(.variables = .(subject_id)) %>% 
  purrr::map(function(x){
    nasal_microbiome_expression_data[,x$sample_id] %>% 
      apply(1, mean)
  }) %>% 
  dplyr::bind_cols() %>% 
  as.data.frame()

rownames(nasal_data) = rownames(nasal_microbiome_expression_data)

colnames(stool_data) == colnames(nasal_data)


stool_nasal_subject_wise_stool_dim = 
  dim(stool_data)

stool_nasal_subject_wise_nasal_dim = 
  dim(nasal_data)

save(stool_nasal_subject_wise_stool_dim, file = "stool_nasal_subject_wise_stool_dim")
save(stool_nasal_subject_wise_nasal_dim, file = "stool_nasal_subject_wise_nasal_dim")


###Sparcc
library(discordant)

library(Biobase)

stool_data1 = round(as.matrix(stool_data) * 20000)
nasal_data1 = round(as.matrix(nasal_data) * 20000)

dim(stool_data1)
dim(nasal_data1)

106*76/2

sparcc_data = NULL

library(tictoc)
rownames(stool_data1) = paste("stool", rownames(stool_data1), sep = "")
rownames(nasal_data1) = paste("nasal", rownames(nasal_data1), sep = "")

# tic()
# cor_data =
# sparcc(data = t(rbind(stool_data1, nasal_data1)),
#        iter = 20,
#        inner_iter = 10,
#        th = 0.1)
# toc()
# 
# cor_data = cor_data$Cor
# 
# colnames(cor_data) = rownames(cor_data) = rownames(rbind(stool_data1, nasal_data1))
# 
# cor_data = as.data.frame(cor_data)
# 
# cor_data[upper.tri(cor_data)] = NA
# 
# library(tidyr)
# cor_data =
# cor_data %>%
#   tibble::rownames_to_column(var = "from") %>%
#   pivot_longer(cols = -from,
#                names_to = "to",
#                values_to = "cor") %>%
#   dplyr::mutate(from_class = case_when(
#     stringr::str_detect(from, "stool") ~ "stool",
#     stringr::str_detect(from, "nasal") ~ "nasal"
#   )) %>%
#   dplyr::mutate(to_class = case_when(
#     stringr::str_detect(to, "stool") ~ "stool",
#     stringr::str_detect(to, "nasal") ~ "nasal"
#   )) %>%
#   # dplyr::filter(from_class != to_class) %>%
#   dplyr::filter(!is.na(cor)) %>%
#   dplyr::filter(cor != 1)
# 
# ####calculate p value for each cor
# ###only calculate p values for abs(cor) > 0.15
# 
# permutation_cor_data =
#   purrr::map(seq_len(100), function(i){
#     cat(i, " ")
#     temp_data =
#       rbind(stool_data1, nasal_data1)
# 
#     temp_data =
#       temp_data %>%
#       apply(1, function(x){
#         x = sample(x)
#       }) %>%
#       t()
# 
#     colnames(temp_data) = colnames(stool_data1)
# 
#     temp_cor_data =
#       sparcc(data = t(temp_data),
#              iter = 20,
#              inner_iter = 10,
#              th = 0.1)
# 
#     temp_cor_data = temp_cor_data$Cor
# 
#     colnames(temp_cor_data) = rownames(temp_cor_data) =
#       rownames(rbind(stool_data1, nasal_data1))
# 
#     temp_cor_data = as.data.frame(temp_cor_data)
# 
#     temp_cor_data[upper.tri(temp_cor_data)] = NA
# 
#     temp_cor_data =
#       temp_cor_data %>%
#       tibble::rownames_to_column(var = "from") %>%
#       pivot_longer(cols = -from,
#                    names_to = "to",
#                    values_to = "cor") %>%
#       dplyr::mutate(from_class = case_when(
#         stringr::str_detect(from, "stool") ~ "stool",
#         stringr::str_detect(from, "nasal") ~ "nasal"
#       )) %>%
#       dplyr::mutate(to_class = case_when(
#         stringr::str_detect(to, "stool") ~ "stool",
#         stringr::str_detect(to, "nasal") ~ "nasal"
#       )) %>%
#       dplyr::filter(!is.na(cor)) %>%
#       dplyr::filter(cor != 1)
#   })
# 
# #####calculate p values for each pair
# 
# p_value =
# purrr::map(seq_len(nrow(cor_data)), function(i){
#   cat(i, " ")
# 
#   permutate_cor =
#   purrr::map(permutation_cor_data, function(x){
#     x$cor[which(x$from == cor_data$from[i] & x$to == cor_data$to[i])]
#   }) %>%
#     unlist()
# 
#   permutate_cor = permutate_cor[!is.na(permutate_cor)]
# 
#   if(cor_data$cor[i] > 0){
#     sum(permutate_cor > cor_data$cor[i])/length(permutate_cor)
#   }else{
#     sum(permutate_cor < cor_data$cor[i])/length(permutate_cor)
#   }
# })
# 
# 
# p_value = unlist(p_value)
# 
# dim(cor_data)
# length(p_value)
# 
# cor_data$p_value = p_value
# 
# head(cor_data)
# 
# cor_data =
#   cor_data %>%
#   dplyr::mutate(p_value_adjust = p.adjust(p_value, method = "BH"))
# 
# save(cor_data, file = "cor_data")

load("cor_data")

plot(cor_data$cor, -log(cor_data$p_value_adjust, 10))


cor_data %>% 
  dplyr::filter(from_class != to_class) %>% 
  dplyr::filter(abs(cor) > 0.25 & p_value_adjust < 0.05) %>% 
  dplyr::arrange(desc(abs(cor)))

plot(log(as.numeric(stool_data1["stoolASV117",])+1, 10),
     log(as.numeric(nasal_data1["nasalASV1",]) + 1, 10))

cor_data %>% 
  dplyr::filter(from_class == to_class) %>% 
  dplyr::filter(abs(cor) > 0.25 & p_value_adjust < 0.05) %>% 
  dplyr::arrange(desc(abs(cor)))

plot(log(as.numeric(stool_data1["stoolASV122",])+1, 10),
     log(as.numeric(stool_data1["stoolASV87",]) + 1, 10))








