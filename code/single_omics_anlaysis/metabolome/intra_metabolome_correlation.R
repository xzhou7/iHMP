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

masstools::setwd_project()
library(tidyverse)
rm(list = ls())

source("code/tools.R")

######work directory
masstools::setwd_project()
setwd("data_analysis/metabolome/intra_metabolome_correlation")

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


###only remain the subjects with at least >= 5
remian_subject_id =
  metabolome_sample_info %>%
  dplyr::group_by(subject_id) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::ungroup() %>%
  dplyr::filter(n > 5) %>%
  dplyr::pull(subject_id)

metabolome_sample_info =
metabolome_sample_info %>%
  dplyr::filter(subject_id %in% remian_subject_id)

metabolome_expression_data =
  metabolome_expression_data[,metabolome_sample_info$sample_id]

##only remain the genus at least in 10% subjects
dim(metabolome_expression_data)

######--------------------------------------------------------------------------
##for raw data, we just log(x+1, 2)
library(plyr)

metabolome_expression_data = 
  log(metabolome_expression_data + 1, 2)

#####
dim(metabolome_expression_data)

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

# metabolome_lm_adjusted_cor =
#   intra_lm_adjusted_cor(expression_data = metabolome_expression_data,
#                         sample_info = metabolome_sample_info)
# 
# save(metabolome_lm_adjusted_cor,
#      file = "metabolome_lm_adjusted_cor",
#      compress = "xz")

load("metabolome_lm_adjusted_cor")


##step 2
##partial correlation
# metabolome_partial_cor =
#   intra_partial_cor(expression_data = metabolome_expression_data,
#                     sample_info = metabolome_sample_info)
# save(metabolome_partial_cor, file = "metabolome_partial_cor", compress = "xz")
load("metabolome_partial_cor")


# ##step 3
# ##rmcorr
# metabolome_rmcor =
#   intra_rm_cor(expression_data = metabolome_expression_data, 
#                sample_info = metabolome_sample_info)
# save(metabolome_rmcor, file = "metabolome_rmcor", compress = "xz")
load("metabolome_rmcor")

sum(metabolome_lm_adjusted_cor$p_adjust < 0.05)
sum(metabolome_partial_cor$p_adjust < 0.05)
sum(metabolome_rmcor$p_adjust < 0.05)




# ###calculate correlation for each subject
# oral_microbiome_metabolome_individual_cor =
#   individual_cor(microbiome_data = oral_microbiome_expression_data,
#                  metabolome_data = metabolome_expression_data,
#                  sample_info = oral_microbiome_sample_info)
# save(oral_microbiome_metabolome_individual_cor,
#      file =  "oral_microbiome_metabolome_individual_cor")

load("oral_microbiome_metabolome_individual_cor")


plot(oral_microbiome_metabolome_rmcor$cor, 
     oral_microbiome_metabolome_partial_cor$cor)

cor(oral_microbiome_metabolome_rmcor$cor, 
     oral_microbiome_metabolome_partial_cor$cor)

plot(oral_microbiome_metabolome_rmcor$cor, 
     oral_microbiome_metabolome_lm_adjusted_cor$cor)

cor(oral_microbiome_metabolome_rmcor$cor, 
     oral_microbiome_metabolome_lm_adjusted_cor$cor)

###here we use the lm_adjusted_cor
sum(oral_microbiome_metabolome_lm_adjusted_cor$p_adjust < 0.05)

oral_microbiome_metabolome_lm_adjusted_cor =
  oral_microbiome_metabolome_lm_adjusted_cor %>%
  dplyr::filter(p_adjust < 0.05)


###output plot
#####network to show the correlation between microbiome and metabolite
oral_microbiome_metabolome_lm_adjusted_cor

library(ggraph)
library(igraph)
library(tidygraph)


###example network
library(plyr)

edge_data =
  oral_microbiome_metabolome_lm_adjusted_cor %>%
  dplyr::mutate(p = -log(p_adjust, 10)) %>%
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
       filename = "oral_microbiome_metabolome_cor_plot_example.pdf",
       width = 8.5,
       height = 7)
