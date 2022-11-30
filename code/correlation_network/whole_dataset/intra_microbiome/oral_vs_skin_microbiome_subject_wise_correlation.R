#' ---
#' title: "oral microbiome skin_microbiome correlation"
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
setwd("data_analysis/correlation_network/whole_data_set/oral_microbiome_vs_skin_microbiome")

# ####load data
# ###oral microbiome
# {
#   load(here::here(
#     "data_analysis/oral_microbiome/data_preparation/expression_data"
#   ))
#   load(here::here(
#     "data_analysis/oral_microbiome/data_preparation/sample_info"
#   ))
#   load(here::here(
#     "data_analysis/oral_microbiome/data_preparation/variable_info"
#   ))
# }
# 
# oral_microbiome_expression_data = expression_data
# oral_microbiome_sample_info = sample_info
# oral_microbiome_variable_info = variable_info
# 
# ###read genus table
# ###read genus table
# expression_data =
#   data.table::fread(here::here("data/from_xin/Genus Table/OR/Genus_OR.csv")) %>%
#   as.data.frame() %>%
#   tibble::column_to_rownames(var = "SampleID") %>%
#   dplyr::select(-c(V1:SubjectID)) %>%
#   t() %>%
#   as.data.frame()
# 
# oral_microbiome_variable_info =
#   oral_microbiome_variable_info[match(rownames(expression_data), oral_microbiome_variable_info$Genus),]
# 
# oral_microbiome_variable_info$Genus == rownames(expression_data)
# 
# rownames(expression_data) = oral_microbiome_variable_info$variable_id
# 
# oral_microbiome_expression_data =
#   expression_data
# 
# oral_microbiome_variable_info =
#   oral_microbiome_variable_info %>%
#   dplyr::filter(!stringr::str_detect(Genus, "Unclassified_Bacteria"))
# 
# oral_microbiome_expression_data =
#   oral_microbiome_expression_data[oral_microbiome_variable_info$variable_id,]
# 
# dim(oral_microbiome_sample_info)
# dim(oral_microbiome_variable_info)
# 
# rownames(oral_microbiome_expression_data) == oral_microbiome_variable_info$variable_id
# colnames(oral_microbiome_expression_data) == oral_microbiome_sample_info$sample_id
# 
# ##skin_microbiome
# ###skin microbiome
# {
#   load(here::here("data_analysis/skin_microbiome/data_preparation/expression_data"))
#   load(here::here("data_analysis/skin_microbiome/data_preparation/sample_info"))
#   load(here::here("data_analysis/skin_microbiome/data_preparation/variable_info"))
# }
# 
# skin_microbiome_expression_data = expression_data
# skin_microbiome_sample_info = sample_info
# skin_microbiome_variable_info = variable_info
# 
# ##read genus table
# expression_data =
#   data.table::fread(here::here("data/from_xin/Genus Table/SK/Genus_SK.csv")) %>%
#   as.data.frame() %>%
#   tibble::column_to_rownames(var = "SampleID") %>%
#   dplyr::select(-c(V1:SubjectID)) %>%
#   t() %>%
#   as.data.frame()
# 
# skin_microbiome_variable_info =
#   skin_microbiome_variable_info[match(rownames(expression_data), skin_microbiome_variable_info$Genus),]
# 
# skin_microbiome_variable_info$Genus == rownames(expression_data)
# 
# ###remove the variables which Genus are NA
# remove_idx = which(is.na(skin_microbiome_variable_info$Genus))
# remove_idx
# if(length(remove_idx) > 0){
#   skin_microbiome_variable_info = skin_microbiome_variable_info[-remove_idx,]
#   expression_data = expression_data[-remove_idx,]
# }
# 
# rownames(expression_data) = skin_microbiome_variable_info$variable_id
# 
# skin_microbiome_expression_data =
#   expression_data
# 
# skin_microbiome_variable_info =
#   skin_microbiome_variable_info %>%
#   dplyr::filter(!stringr::str_detect(Genus, "Unclassified_Bacteria"))
# 
# skin_microbiome_expression_data =
#   skin_microbiome_expression_data[skin_microbiome_variable_info$variable_id,]
# 
# dim(skin_microbiome_sample_info)
# dim(skin_microbiome_variable_info)
# dim(skin_microbiome_expression_data)
# 
# rownames(skin_microbiome_expression_data) == skin_microbiome_variable_info$variable_id
# colnames(skin_microbiome_expression_data) == skin_microbiome_sample_info$sample_id
# 
# ###match samples
# dim(oral_microbiome_sample_info)
# dim(skin_microbiome_sample_info)
# 
# length(unique(oral_microbiome_sample_info$subject_id))
# 
# ###just matched samples according to sample id, only 1 missed
# intersect_sample_id =
#   intersect(oral_microbiome_sample_info$sample_id,
#             skin_microbiome_sample_info$sample_id)
# 
# length(intersect_sample_id)
# 
# oral_microbiome_expression_data =
#   oral_microbiome_expression_data[,intersect_sample_id]
# 
# skin_microbiome_expression_data =
#   skin_microbiome_expression_data[,intersect_sample_id]
# 
# oral_microbiome_sample_info =
#   oral_microbiome_sample_info[match(intersect_sample_id, oral_microbiome_sample_info$sample_id),]
# 
# skin_microbiome_sample_info =
#   skin_microbiome_sample_info[match(intersect_sample_id, skin_microbiome_sample_info$sample_id),]
# 
# length(unique(skin_microbiome_sample_info$subject_id))
# 
# ###only remain the subjects with at least >= 5
# remian_subject_id =
#   oral_microbiome_sample_info %>%
#   dplyr::group_by(subject_id) %>%
#   dplyr::summarise(n = n()) %>%
#   dplyr::ungroup() %>%
#   dplyr::filter(n > 5) %>%
#   dplyr::pull(subject_id)
# 
# skin_microbiome_sample_info =
#   skin_microbiome_sample_info %>%
#   dplyr::filter(subject_id %in% remian_subject_id)
# 
# skin_microbiome_expression_data =
#   skin_microbiome_expression_data[, skin_microbiome_sample_info$sample_id]
# 
# oral_microbiome_sample_info =
#   oral_microbiome_sample_info %>%
#   dplyr::filter(subject_id %in% remian_subject_id)
# 
# oral_microbiome_expression_data =
#   oral_microbiome_expression_data[,oral_microbiome_sample_info$sample_id]
# 
# ##only remain the genus at least in 10% subjects
# remain_idx =
#   which(rowSums(oral_microbiome_expression_data) > 0)
# 
# oral_microbiome_expression_data = oral_microbiome_expression_data[remain_idx,]
# oral_microbiome_variable_info = oral_microbiome_variable_info[remain_idx,,drop = FALSE]
# 
# remain_idx =
#   oral_microbiome_expression_data %>%
#   apply(1, function(x){
#     sum(as.numeric(x) == 0) / ncol(oral_microbiome_expression_data)
#   }) %>%
#   `<`(0.9) %>%
#   which()
# 
# length(remain_idx)
# 
# oral_microbiome_expression_data = oral_microbiome_expression_data[remain_idx,]
# oral_microbiome_variable_info = oral_microbiome_variable_info[remain_idx,,drop = FALSE]
# 
# ##only remain the genus at least in 10% subjects
# remain_idx =
#   which(rowSums(skin_microbiome_expression_data) > 0)
# 
# skin_microbiome_expression_data = skin_microbiome_expression_data[remain_idx,]
# skin_microbiome_variable_info = skin_microbiome_variable_info[remain_idx,,drop = FALSE]
# 
# remain_idx =
#   skin_microbiome_expression_data %>%
#   apply(1, function(x){
#     sum(as.numeric(x) == 0) / ncol(skin_microbiome_expression_data)
#   }) %>%
#   `<`(0.9) %>%
#   which()
# 
# length(remain_idx)
# 
# skin_microbiome_expression_data = skin_microbiome_expression_data[remain_idx,]
# skin_microbiome_variable_info = skin_microbiome_variable_info[remain_idx,,drop = FALSE]
# 
# ##save data
# {
#   save(oral_microbiome_expression_data, file = "oral_microbiome_expression_data")
#   save(oral_microbiome_variable_info, file = "oral_microbiome_variable_info")
#   save(oral_microbiome_sample_info, file = "oral_microbiome_sample_info")
# 
#   save(skin_microbiome_expression_data, file = "skin_microbiome_expression_data")
#   save(skin_microbiome_variable_info, file = "skin_microbiome_variable_info")
#   save(skin_microbiome_sample_info, file = "skin_microbiome_sample_info")
# }

{
  load("oral_microbiome_expression_data")
  load("oral_microbiome_variable_info")
  load("oral_microbiome_sample_info")

  load("skin_microbiome_expression_data")
  load("skin_microbiome_variable_info")
  load("skin_microbiome_sample_info")
}

dim(oral_microbiome_expression_data)
dim(skin_microbiome_expression_data)

#####set to other folder
dir.create("subject_wise")
setwd("subject_wise")

###finally, for oral microbiome, 106 genus, for skin_microbiome, 76 genus

######--------------------------------------------------------------------------
library(plyr)

#####
dim(skin_microbiome_expression_data)
dim(oral_microbiome_expression_data)

oral_microbiome_sample_info$subject_id == skin_microbiome_sample_info$subject_id

##https://rpubs.com/DKCH2020/578881
##https://ourcodingclub.github.io/tutorials/mixed-models/
##https://zhuanlan.zhihu.com/p/63092231
##https://www.linglab.cn/knowledge/10

#######for each subject, just get the mean value for each genera
oral_microbiome_sample_info %>% 
  dplyr::count(subject_id)

library(plyr)

oral_data = 
  oral_microbiome_sample_info %>%
  plyr::dlply(.variables = .(subject_id)) %>% 
  purrr::map(function(x){
    oral_microbiome_expression_data[,x$sample_id] %>% 
      apply(1, mean)
  }) %>% 
  dplyr::bind_cols() %>% 
  as.data.frame()
  
rownames(oral_data) = rownames(oral_microbiome_expression_data)


skin_data = 
  skin_microbiome_sample_info %>%
  plyr::dlply(.variables = .(subject_id)) %>% 
  purrr::map(function(x){
    skin_microbiome_expression_data[,x$sample_id] %>% 
      apply(1, mean)
  }) %>% 
  dplyr::bind_cols() %>% 
  as.data.frame()

rownames(skin_data) = rownames(skin_microbiome_expression_data)


colnames(oral_data) == colnames(skin_data)


oral_skin_subject_wise_oral_dim = 
  dim(oral_data)

oral_skin_subject_wise_skin_dim = 
  dim(skin_data)

save(oral_skin_subject_wise_oral_dim, file = "oral_skin_subject_wise_oral_dim")
save(oral_skin_subject_wise_skin_dim, file = "oral_skin_subject_wise_skin_dim")


###Sparcc
library(discordant)

library(Biobase)

oral_data1 = round(as.matrix(oral_data) * 20000)
skin_data1 = round(as.matrix(skin_data) * 20000)

dim(oral_data1)
dim(skin_data1)

106*76/2

sparcc_data = NULL

library(tictoc)
rownames(oral_data1) = paste("oral", rownames(oral_data1), sep = "")
rownames(skin_data1) = paste("skin", rownames(skin_data1), sep = "")

# tic()
# cor_data =
# sparcc(data = t(rbind(oral_data1, skin_data1)),
#        iter = 20,
#        inner_iter = 10,
#        th = 0.1)
# toc()
# 
# cor_data = cor_data$Cor
# 
# colnames(cor_data) = rownames(cor_data) = rownames(rbind(oral_data1, skin_data1))
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
#     stringr::str_detect(from, "oral") ~ "oral",
#     stringr::str_detect(from, "skin") ~ "skin"
#   )) %>%
#   dplyr::mutate(to_class = case_when(
#     stringr::str_detect(to, "oral") ~ "oral",
#     stringr::str_detect(to, "skin") ~ "skin"
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
#       rbind(oral_data1, skin_data1)
# 
#     temp_data =
#       temp_data %>%
#       apply(1, function(x){
#         x = sample(x)
#       }) %>%
#       t()
# 
#     colnames(temp_data) = colnames(oral_data1)
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
#       rownames(rbind(oral_data1, skin_data1))
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
#         stringr::str_detect(from, "oral") ~ "oral",
#         stringr::str_detect(from, "skin") ~ "skin"
#       )) %>%
#       dplyr::mutate(to_class = case_when(
#         stringr::str_detect(to, "oral") ~ "oral",
#         stringr::str_detect(to, "skin") ~ "skin"
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

plot(log(as.numeric(oral_data1["oralASV117",])+1, 10),
     log(as.numeric(skin_data1["skinASV1",]) + 1, 10))

cor_data %>% 
  dplyr::filter(from_class == to_class) %>% 
  dplyr::filter(abs(cor) > 0.25 & p_value_adjust < 0.05) %>% 
  dplyr::arrange(desc(abs(cor)))

plot(log(as.numeric(oral_data1["oralASV122",])+1, 10),
     log(as.numeric(oral_data1["oralASV87",]) + 1, 10))








