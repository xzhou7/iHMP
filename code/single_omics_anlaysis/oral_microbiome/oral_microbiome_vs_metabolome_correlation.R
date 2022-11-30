#' ---
#' title: "nasal microbiome metabolome correlation"
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
setwd("data_analysis/nasal_microbiome_vs_metabolome")

####load data
###nasal microbiome
{
  load(here::here("data_analysis/nasal_microbiome/data_preparation/expression_data"))
  load(here::here("data_analysis/nasal_microbiome/data_preparation/sample_info"))
  load(here::here("data_analysis/nasal_microbiome/data_preparation/variable_info"))
}

nasal_microbiome_expression_data = expression_data
nasal_microbiome_sample_info = sample_info
nasal_microbiome_variable_info = variable_info

###read genus table
expression_data =
readr::read_csv(here::here("data/from_xin/Genus Table/NS/Genus_NS.csv")) %>%
  tibble::column_to_rownames(var = "SampleID") %>%
  dplyr::select(-c(X1:batch)) %>%
  t() %>%
  as.data.frame()

nasal_microbiome_variable_info =
  data.frame(variable_id = rownames(expression_data))

nasal_microbiome_expression_data =
  expression_data

dim(nasal_microbiome_sample_info)
length(unique(nasal_microbiome_sample_info$subject_id))
dim(nasal_microbiome_variable_info)

###plasma metabolome
{
  load(here::here("data_analysis/metabolome/data_preparation/expression_data"))
  load(here::here("data_analysis/metabolome/data_preparation/sample_info"))
  load(here::here("data_analysis/metabolome/data_preparation/variable_info"))  
}

metabolome_expression_data = expression_data
metabolome_sample_info = sample_info
metabolome_variable_info = variable_info

metabolome_sample_info$CollectionDate =
  as.Date(metabolome_sample_info$CollectionDate, "%m/%d/%y")

dim(metabolome_expression_data)
length(unique(metabolome_sample_info$subject_id))

# ###match samples
# dim(nasal_microbiome_sample_info)
# dim(metabolome_sample_info)
# 
# length(nasal_microbiome_sample_info$subject_id)
# length(unique(nasal_microbiome_sample_info$subject_id))
# 
# ###just matched samples according to sample id, only 1 missed
# intersect_sample_id =
#   intersect(nasal_microbiome_sample_info$sample_id,
#             metabolome_sample_info$sample_id)
# 
# length(intersect_sample_id)
# 
# nasal_microbiome_expression_data =
#   nasal_microbiome_expression_data[,intersect_sample_id]
# 
# metabolome_expression_data =
#   metabolome_expression_data[,intersect_sample_id]
# 
# nasal_microbiome_sample_info =
#   nasal_microbiome_sample_info[match(intersect_sample_id, nasal_microbiome_sample_info$sample_id),]
# 
# metabolome_sample_info =
#   metabolome_sample_info[match(intersect_sample_id, metabolome_sample_info$sample_id),]
# 
# length(unique(metabolome_sample_info$subject_id))
# 
# ###only remain the subjects with at least >= 5
# remian_subject_id =
# nasal_microbiome_sample_info %>%
#   dplyr::group_by(subject_id) %>%
#   dplyr::summarise(n = n()) %>%
#   dplyr::ungroup() %>%
#   dplyr::filter(n > 5) %>%
#   dplyr::pull(subject_id)
# 
# metabolome_sample_info =
# metabolome_sample_info %>%
#   dplyr::filter(subject_id %in% remian_subject_id)
# 
# metabolome_expression_data =
#   metabolome_expression_data[,metabolome_sample_info$sample_id]
# 
# nasal_microbiome_sample_info =
#   nasal_microbiome_sample_info %>%
#   dplyr::filter(subject_id %in% remian_subject_id)
# 
# nasal_microbiome_expression_data =
#   nasal_microbiome_expression_data[,nasal_microbiome_sample_info$sample_id]
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
#   save(nasal_microbiome_expression_data, file = "nasal_microbiome_expression_data")
#   save(nasal_microbiome_variable_info, file = "nasal_microbiome_variable_info")
#   save(nasal_microbiome_sample_info, file = "nasal_microbiome_sample_info")
# 
#   save(metabolome_expression_data, file = "metabolome_expression_data")
#   save(metabolome_variable_info, file = "metabolome_variable_info")
#   save(metabolome_sample_info, file = "metabolome_sample_info")
# }

{
  load("nasal_microbiome_expression_data")
  load("nasal_microbiome_variable_info")
  load("nasal_microbiome_sample_info")
  
  load("metabolome_expression_data")
  load("metabolome_variable_info")
  load("metabolome_sample_info")  
}

dim(nasal_microbiome_expression_data)

######--------------------------------------------------------------------------
##for raw data, we just log(x+1, 2)
library(plyr)

metabolome_expression_data = 
  log(metabolome_expression_data + 1, 2)

##for microbiome data, we just log(x+1, 2)
library(plyr)

###because our microbiome are percentage data, so here we use the CTL method
library(compositions)
nasal_microbiome_expression_data = 
nasal_microbiome_expression_data %>% 
  purrr::map(function(x){
    x = compositions::clr(x) %>% 
      as.numeric()
    x
  }) %>% 
  do.call(cbind, .) %>% 
    as.data.frame()

rownames(nasal_microbiome_expression_data) = nasal_microbiome_variable_info$variable_id

#####
dim(metabolome_expression_data)
dim(nasal_microbiome_expression_data)

nasal_microbiome_sample_info$subject_id == metabolome_sample_info$subject_id

##https://rpubs.com/DKCH2020/578881
##https://ourcodingclub.github.io/tutorials/mixed-models/
##https://zhuanlan.zhihu.com/p/63092231
##https://www.linglab.cn/knowledge/10

###step 1
###linear mixed model to adjust the subject ID random effect, and then use the partial correlation
###to get the correlation between microbiome and metabolome
library(lme4)
library(rmcorr)

library(future)
library(furrr)

plan(strategy = multisession(workers = 5))

# nasal_microbiome_metabolome_lm_adjusted_cor =
#   lm_adjusted_cor(
#     microbiome_data = nasal_microbiome_expression_data,
#     metabolome_data = metabolome_expression_data,
#     sample_info = nasal_microbiome_sample_info
#   )
# save(nasal_microbiome_metabolome_lm_adjusted_cor,
# file = "nasal_microbiome_metabolome_lm_adjusted_cor", compress = "xz")

load("nasal_microbiome_metabolome_lm_adjusted_cor")


# ##step 2
# ##partial correlation
# nasal_microbiome_metabolome_partial_cor =
#   partial_cor(microbiome_data = nasal_microbiome_expression_data,
#               metabolome_data = metabolome_expression_data,
#               sample_info = nasal_microbiome_sample_info)
# save(nasal_microbiome_metabolome_partial_cor, file = "nasal_microbiome_metabolome_partial_cor", compress = "xz")
load("nasal_microbiome_metabolome_partial_cor")


# ##step 3
# ##rmcorr
# nasal_microbiome_metabolome_rmcor =
#   rm_cor(
#     microbiome_data = nasal_microbiome_expression_data,
#     metabolome_data = metabolome_expression_data,
#     sample_info = nasal_microbiome_sample_info
#   )
# save(nasal_microbiome_metabolome_rmcor, file = "nasal_microbiome_metabolome_rmcor", compress = "xz")
load("nasal_microbiome_metabolome_rmcor")

sum(nasal_microbiome_metabolome_rmcor$p_adjust2 < 0.05)
sum(nasal_microbiome_metabolome_lm_adjusted_cor$p_adjust2 < 0.05)
sum(nasal_microbiome_metabolome_partial_cor$p_adjust2 < 0.05)

# ###calculate correlation for each subject
# nasal_microbiome_metabolome_individual_cor =
#   individual_cor(microbiome_data = nasal_microbiome_expression_data,
#                  metabolome_data = metabolome_expression_data,
#                  sample_info = nasal_microbiome_sample_info)
# save(nasal_microbiome_metabolome_individual_cor,
#      file =  "nasal_microbiome_metabolome_individual_cor")

load("nasal_microbiome_metabolome_individual_cor")


plot(nasal_microbiome_metabolome_rmcor$cor, 
     nasal_microbiome_metabolome_partial_cor$cor)

cor(nasal_microbiome_metabolome_rmcor$cor, 
     nasal_microbiome_metabolome_partial_cor$cor)

plot(nasal_microbiome_metabolome_rmcor$cor, 
     nasal_microbiome_metabolome_lm_adjusted_cor$cor)

cor(nasal_microbiome_metabolome_rmcor$cor, 
     nasal_microbiome_metabolome_lm_adjusted_cor$cor)

###here we use the lm_adjusted_cor
sum(nasal_microbiome_metabolome_lm_adjusted_cor$p_adjust2 < 0.05)

nasal_microbiome_metabolome_lm_adjusted_cor =
  nasal_microbiome_metabolome_lm_adjusted_cor %>%
  dplyr::filter(p_adjust2 < 0.05)


###output plot
#####network to show the correlation between microbiome and metabolite
nasal_microbiome_metabolome_lm_adjusted_cor

library(ggraph)
library(igraph)
library(tidygraph)


###example network
library(plyr)

edge_data =
  nasal_microbiome_metabolome_lm_adjusted_cor %>%
  dplyr::mutate(p = -log(p_adjust2, 10)) %>%
  dplyr::select(microbiome, metabolite, p, cor) %>%
  dplyr::rename(from = microbiome,
                to = metabolite,
                cor = cor) 
  
  node_data = 
  data.frame(node = unique(c(edge_data$from, edge_data$to))) %>% 
  dplyr::left_join(metabolome_variable_info[,c("variable_id", "Metabolite")],
                   by = c("node" = "variable_id")) %>% 
  dplyr::rename(true_name = Metabolite) %>% 
    dplyr::mutate(class = case_when(
      is.na(true_name) ~ "Nasal microbiome",
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
                                  # label = ifelse(class == "Metabolite", NA, true_name),
                                  label = true_name,
                                  color = class), bg.colour = "white") +
  guides(edge_color = guide_edge_colorbar(title = "Correlation"),
         color = guide_legend(title = "Class")) +
  scale_edge_width_continuous(range = c(0.1,2)) +
  scale_size_continuous(range = c(5, 10)) +
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

ggsave(plot,
       filename = "nasal_microbiome_metabolome_cor_plot_example.pdf",
       width = 8.5,
       height = 7)
