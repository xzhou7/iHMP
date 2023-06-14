#' ---
#' title: "oral microbiome metabolome correlation"
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

setwd(masstools::get_project_wd())
library(tidyverse)
rm(list = ls())

source("code/tools.R")

####load data
###load diversity
load("data/from_xin/Diversity_Datatable.RData")

##prevalence
load("data/from_xin/Prevalance.RData")

###oral microbiome
load("data_analysis/oral_microbiome/data_preparation/expression_data")
load("data_analysis/oral_microbiome/data_preparation/sample_info")
load("data_analysis/oral_microbiome/data_preparation/variable_info")

oral_microbiome_expression_data = expression_data
oral_microbiome_sample_info = sample_info
oral_microbiome_variable_info = variable_info

expression_data = 
readr::read_csv(here::here("data/from_xin/Genus Table/OR/Genus_OR.csv")) %>% 
  tibble::column_to_rownames(var = "SampleID") %>% 
  dplyr::select(-c(X1:SubjectID)) %>% 
  t() %>% 
  as.data.frame()

oral_microbiome_variable_info = 
  data.frame(variable_id = rownames(expression_data))

oral_microbiome_expression_data = 
  expression_data

sum(colnames(oral_microbiome_expression_data) == oral_microbiome_sample_info$sample_id)

dim(oral_microbiome_sample_info)

#######work directory
setwd(masstools::get_project_wd())
setwd("data_analysis/intra_microbiome_correlation")

##only remain the genus at least in 10% subjects
remain_idx = 
  which(rowSums(oral_microbiome_expression_data) > 0)

oral_microbiome_expression_data = oral_microbiome_expression_data[remain_idx,]
oral_microbiome_variable_info = oral_microbiome_variable_info[remain_idx,,drop = FALSE]

remain_idx = 
  oral_microbiome_expression_data %>% 
  apply(1, function(x){
    sum(as.numeric(x) == 0) / ncol(oral_microbiome_expression_data)
  }) %>% 
  `<`(0.9) %>% 
  which()

length(remain_idx)

oral_microbiome_expression_data = oral_microbiome_expression_data[remain_idx,]
oral_microbiome_variable_info = oral_microbiome_variable_info[remain_idx,,drop = FALSE]

######--------------------------------------------------------------------------
##for microbiome data, we just log(x+1, 2)
library(plyr)

oral_microbiome_expression_data = 
  log(oral_microbiome_expression_data + 1, 2)

##scaled
oral_microbiome_expression_data =
  oral_microbiome_expression_data %>%
  apply(1, function(x){
    (x - mean(x))/sd(x)
  }) %>%
  t() %>%
  as.data.frame()


#####
dim(oral_microbiome_expression_data)

oral_microbiome_sample_info$subject_id == metabolome_sample_info$subject_id

my.rmc = 
rmcorr(participant = Subject, measure1 = PacO2, measure2 = pH, dataset = bland1995)
bland1995 %>% 
  ggplot(aes(x = PacO2, y = pH)) +
  geom_point(aes(color = as.character(Subject))) +
  geom_smooth(aes(color = as.character(Subject)), method = "lm", se = FALSE)




##https://rpubs.com/DKCH2020/578881
##https://ourcodingclub.github.io/tutorials/mixed-models/
##https://zhuanlan.zhihu.com/p/63092231
##https://www.linglab.cn/knowledge/10

###step 1
###liner mixed model to find the association between microbiome and metabolites
library(lme4)

# oral_microbiome_metabolome_cor =
#   rownames(oral_microbiome_expression_data) %>%
#   purrr::map(function(name) {
#     cat(name, " ")
#     x = as.numeric(oral_microbiome_expression_data[name,])
#     temp_partial_cor =
#     purrr::map(
#       as.data.frame(t(metabolome_expression_data)),
#       .f = function(y) {
#         temp_data =
#           data.frame(x = x,
#                      y = y,
#                      oral_microbiome_sample_info)
#         temp_data$Gender[temp_data$Gender == 'F'] = 0
#         temp_data$Gender[temp_data$Gender == 'M'] = 1
#         temp_data$Gender = as.numeric(temp_data$Gender)
# 
#         temp_data$Ethnicity[temp_data$Ethnicity == 'C'] = 0
#         temp_data$Ethnicity[temp_data$Ethnicity == 'H'] = 1
#         temp_data$Ethnicity[temp_data$Ethnicity == 'B'] = 2
#         temp_data$Ethnicity[temp_data$Ethnicity == 'A'] = 3
#         temp_data$Ethnicity = as.numeric(temp_data$Ethnicity)
# 
#         temp_data$SSPG = as.numeric(temp_data$SSPG)
#         temp_data$FPG = as.numeric(temp_data$FPG)
# 
#         ##remove samples without SSPG or FPG
#         temp_data =
#         temp_data %>%
#           dplyr::filter(!is.na(SSPG)) %>%
#           dplyr::filter(!is.na(FPG))
# 
#         ##only remain the subjects with less than 5 time points
#         remain_subject_id =
#           temp_data %>%
#           dplyr::select(subject_id) %>%
#           dplyr::group_by(subject_id) %>%
#           dplyr::summarise(n = n()) %>%
#           dplyr::filter(n >= 5) %>%
#           pull(subject_id)
# 
#         temp_data =
#           temp_data %>%
#           dplyr::filter(subject_id %in% remain_subject_id)
# 
#         ##partial correlation
#         cor_value =
#           ppcor::pcor.test(
#             x = temp_data$x,
#             y = temp_data$y,
#             z = temp_data[, c("Gender", "Adj.age", "Ethnicity", "SSPG", "FPG")],
#             method = "spearman"
#           )
# 
#         result =
#           c(
#             cor = unname(cor_value$estimate),
#             cor_p = unname(cor_value$p.value)
#           )
#         if (is.na(result[1])) {
#           result[1] = 0
#         }
# 
#         if (is.na(result[2])) {
#           result[2] = 1
#         }
#         result
#       }
#     ) %>%
#       do.call(rbind, .) %>%
#         as.data.frame()
# 
#     temp_partial_cor =
#     temp_partial_cor %>%
#       tibble::rownames_to_column(var = "metabolite") %>%
#       dplyr::mutate(microbiome = name) %>%
#       dplyr::select(microbiome, metabolite, dplyr::everything())
#     temp_partial_cor$cor_p_adjust = p.adjust(temp_partial_cor$cor_p, method = "BH")
#     temp_partial_cor
#   }) %>%
#   do.call(rbind, .) %>%
#   as.data.frame()
# 
# save(oral_microbiome_metabolome_cor, file = "oral_microbiome_metabolome_cor", compress = "xz")

load("oral_microbiome_metabolome_cor")

oral_microbiome_metabolome_cor$cor_p_adjust2 = 
  p.adjust(oral_microbiome_metabolome_cor$cor_p, method = "BH")

oral_microbiome_metabolome_cor =
  oral_microbiome_metabolome_cor %>%
  dplyr::filter(cor_p_adjust < 0.01)

which(abs(oral_microbiome_metabolome_cor$cor) > 0.4)

###remove the variables with low correlation
plot(oral_microbiome_metabolome_cor$cor)

###output plot
#####network to show the correlation between microbiome and metabolite
oral_microbiome_metabolome_cor

library(ggraph)
library(igraph)
library(tidygraph)

###example network
library(plyr)
edge_data =
  oral_microbiome_metabolome_cor %>%
  dplyr::mutate(p = -log(cor_p_adjust2, 10)) %>%
  dplyr::select(microbiome, metabolite, p, cor) %>%
  dplyr::rename(from = microbiome,
                to = metabolite,
                cor = cor) %>% 
  dplyr::filter(abs(cor) > 0.25)
  
node_data = 
  data.frame(node = unique(c(edge_data$from, edge_data$to))) %>% 
  dplyr::left_join(metabolome_variable_info[,c("variable_id", "Metabolite")],
                   by = c("node" = "variable_id")) %>% 
  dplyr::rename(true_name = Metabolite) %>% 
    dplyr::mutate(class = case_when(
      is.na(true_name) ~ "Oral microbiome",
      !is.na(true_name) ~ "Metabolite"
    )) 

node_data$true_name[is.na(node_data$true_name)] = 
  node_data$node[is.na(node_data$true_name)]

node_data = 
node_data %>% 
  dplyr::filter(!stringr::str_detect(true_name, "C[0-9]{1,2}H[0-9]{1,2}"))

node_data$true_name = stringr::str_replace_all(node_data$true_name, '\"', "")

edge_data = 
edge_data %>% 
  dplyr::filter(from %in% node_data$node & to %in% node_data$node)

node_data = 
node_data %>% 
  dplyr::filter(node %in% edge_data$from | node %in% edge_data$to)

temp_data <- 
  tidygraph::tbl_graph(nodes = node_data, 
                       edges = edge_data,
                       directed = FALSE) %>% 
  dplyr::mutate(Degree = centrality_degree(mode = 'all'))

angle <- 360 * (c(1:nrow(node_data)) - 0.5)/nrow(node_data)
hjust <- ifelse(angle > 180, 1, 0)
angle <- ifelse(angle > 180, 90 - angle + 180, 90 - angle)

plot <-
  ggraph(temp_data,
         layout = 'fr',
         circular = FALSE) +
  geom_edge_link(strength = 1,
                aes(color = cor, width = p),
                alpha = 1,
                show.legend = TRUE) +
  geom_node_point(aes(size = Degree,
                      color = class), 
                  shape = 16, 
                  alpha = 1, show.legend = TRUE) +
  # geom_node_text(aes(x = x ,
  #                    y = y ,
  #                    color = "black",
  #                    label = true_name),
  #                repel = TRUE, 
  #                size = 2) +
  shadowtext::geom_shadowtext(aes(x = x, y = y,
                                  label = ifelse(class == "Metabolite", NA, true_name),
                                  color = class), bg.colour = "white") +
  guides(edge_color = guide_edge_colorbar(title = "Correlation"),
         color = guide_legend(title = "Class")) +
  scale_edge_width_continuous(range = c(0.1,2)) +
  scale_size_continuous(range = c(3, 10)) +
  scale_color_manual(values = omics_color) +
  scale_edge_color_gradient2(low = viridis::cividis(n = 2)[1], 
                             mid = "white", 
                             high = "red") +
  ggraph::theme_graph() +
  theme(plot.background = element_rect(fill = "transparent", color = NA), 
        panel.background = element_rect(fill = "transparent", color = NA),
        legend.position = "right",
        legend.background = element_rect(fill = "transparent", color = NA))

plot

extrafont::loadfonts()

# ggsave(plot,
#        filename = "oral_microbiome_metabolome_cor_plot_example.pdf",
#        width = 8.5,
#        height = 7)








###complete network
library(plyr)
edge_data =
  oral_microbiome_metabolome_cor %>%
  dplyr::mutate(p = -log(cor_p_adjust2, 10)) %>%
  dplyr::select(microbiome, metabolite, p, cor) %>%
  dplyr::rename(from = microbiome,
                to = metabolite,
                cor = cor) 
  # dplyr::filter(abs(cor) > 0.25)

node_data = 
  data.frame(node = unique(c(edge_data$from, edge_data$to))) %>% 
  dplyr::left_join(metabolome_variable_info[,c("variable_id", "Metabolite")],
                   by = c("node" = "variable_id")) %>% 
  dplyr::rename(true_name = Metabolite) %>% 
  dplyr::mutate(class = case_when(
    is.na(true_name) ~ "Oral microbiome",
    !is.na(true_name) ~ "Metabolite"
  )) 

node_data$true_name[is.na(node_data$true_name)] = 
  node_data$node[is.na(node_data$true_name)]

node_data = 
  node_data %>% 
  dplyr::filter(!stringr::str_detect(true_name, "C[0-9]{1,2}H[0-9]{1,2}"))

node_data$true_name = stringr::str_replace_all(node_data$true_name, '\"', "")

edge_data = 
  edge_data %>% 
  dplyr::filter(from %in% node_data$node & to %in% node_data$node)

node_data = 
  node_data %>% 
  dplyr::filter(node %in% edge_data$from | node %in% edge_data$to)

temp_data <- 
  tidygraph::tbl_graph(nodes = node_data, 
                       edges = edge_data,
                       directed = FALSE) %>% 
  dplyr::mutate(Degree = centrality_degree(mode = 'all'))

angle <- 360 * (c(1:nrow(node_data)) - 0.5)/nrow(node_data)
hjust <- ifelse(angle > 180, 1, 0)
angle <- ifelse(angle > 180, 90 - angle + 180, 90 - angle)

plot <-
  ggraph(temp_data,
         layout = 'kk',
         circular = FALSE) +
  geom_edge_link(strength = 1,
                 aes(color = cor, width = p),
                 alpha = 1,
                 show.legend = TRUE) +
  geom_node_point(aes(size = Degree,
                      color = class), 
                  shape = 16, 
                  alpha = 1, show.legend = TRUE) +
  # geom_node_text(aes(x = x ,
  #                    y = y ,
  #                    color = "black",
  #                    label = true_name),
  #                repel = TRUE, 
  #                size = 2) +
  # shadowtext::geom_shadowtext(aes(x = x, y = y,
  #                                 label = ifelse(class == "Metabolite", NA, true_name),
  #                                 color = class), bg.colour = "white") +
  guides(edge_color = guide_edge_colorbar(title = "Correlation"),
         color = guide_legend(title = "Class")) +
  scale_edge_width_continuous(range = c(0.1,2)) +
  scale_size_continuous(range = c(3, 10)) +
  scale_color_manual(values = omics_color) +
  scale_edge_color_gradient2(low = viridis::cividis(n = 2)[1], 
                             mid = "white", 
                             high = "red") +
  ggraph::theme_graph() +
  theme(plot.background = element_rect(fill = "transparent", color = NA), 
        panel.background = element_rect(fill = "transparent", color = NA),
        legend.position = "right",
        legend.background = element_rect(fill = "transparent", color = NA))

plot

extrafont::loadfonts()

# ggsave(plot,
#        filename = "oral_microbiome_metabolome_cor_plot.pdf",
#        width = 8.5,
#        height = 7)



