#' ---
#' title: "Nasal microbiome metabolome correlation"
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

######work directory
setwd(masstools::get_project_wd())
load("data/from_xin/Richness_Bygenus.RData")
load("data_analysis/nasal_microbiome/data_preparation/sample_info")
setwd("data_analysis/correlation_network/whole_data_set_enrichness/nasal_microbiome_vs_metabolome")

nasal_microbiome_expression_data =
  t(richness_NS_byGenus) %>%
  as.data.frame()

colnames(nasal_microbiome_expression_data) = 
  sample_info$sample_id[match(colnames(nasal_microbiome_expression_data), sample_info$RandomID)]

nasal_microbiome_sample_info = 
  sample_info[match(colnames(nasal_microbiome_expression_data), sample_info$sample_id),]

colnames(nasal_microbiome_expression_data) == nasal_microbiome_sample_info$sample_id

nasal_microbiome_variable_info = 
  data.frame(variable_id = rownames(nasal_microbiome_expression_data),
             Genus = rownames(nasal_microbiome_expression_data))

###plasma metabolome
{
  load(here::here(
    "data_analysis/metabolome/data_preparation/expression_data"
  ))
  load(here::here("data_analysis/metabolome/data_preparation/sample_info"))
  load(here::here(
    "data_analysis/metabolome/data_preparation/variable_info"
  ))
}

metabolome_expression_data = expression_data
metabolome_sample_info = sample_info
metabolome_variable_info = variable_info

metabolome_sample_info$CollectionDate =
  as.Date(metabolome_sample_info$CollectionDate, "%m/%d/%y")

dim(metabolome_expression_data)
length(unique(metabolome_sample_info$subject_id))

###match samples
dim(nasal_microbiome_sample_info)
dim(metabolome_sample_info)

length(nasal_microbiome_sample_info$subject_id)
length(unique(nasal_microbiome_sample_info$subject_id))

###just matched samples according to sample id, only 1 missed
intersect_sample_id =
  intersect(nasal_microbiome_sample_info$sample_id,
            metabolome_sample_info$sample_id)

length(intersect_sample_id)

nasal_microbiome_expression_data =
  nasal_microbiome_expression_data[,intersect_sample_id]

metabolome_expression_data =
  metabolome_expression_data[,intersect_sample_id]

nasal_microbiome_sample_info =
  nasal_microbiome_sample_info[match(intersect_sample_id, nasal_microbiome_sample_info$sample_id),]

metabolome_sample_info =
  metabolome_sample_info[match(intersect_sample_id, metabolome_sample_info$sample_id),]

length(unique(metabolome_sample_info$subject_id))

###only remain the subjects with at least >= 5
remian_subject_id =
  nasal_microbiome_sample_info %>%
  dplyr::group_by(subject_id) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::ungroup() %>%
  dplyr::filter(n > 5) %>%
  dplyr::pull(subject_id)

metabolome_sample_info =
  metabolome_sample_info %>%
  dplyr::filter(subject_id %in% remian_subject_id)

metabolome_expression_data =
  metabolome_expression_data[, metabolome_sample_info$sample_id]

nasal_microbiome_sample_info =
  nasal_microbiome_sample_info %>%
  dplyr::filter(subject_id %in% remian_subject_id)

nasal_microbiome_expression_data =
  nasal_microbiome_expression_data[,nasal_microbiome_sample_info$sample_id]

##only remain the genus at least in 10% subjects
remain_idx =
  which(rowSums(nasal_microbiome_expression_data) > 0)

nasal_microbiome_expression_data = nasal_microbiome_expression_data[remain_idx,]
nasal_microbiome_variable_info = nasal_microbiome_variable_info[remain_idx,,drop = FALSE]

remain_idx =
  nasal_microbiome_expression_data %>%
  apply(1, function(x){
    sum(as.numeric(x) == 0) / ncol(nasal_microbiome_expression_data)
  }) %>%
  `<`(0.9) %>%
  which()

length(remain_idx)

nasal_microbiome_expression_data = nasal_microbiome_expression_data[remain_idx,]
nasal_microbiome_variable_info = nasal_microbiome_variable_info[remain_idx,,drop = FALSE]

##save data
{
  save(nasal_microbiome_expression_data, file = "nasal_microbiome_expression_data")
  save(nasal_microbiome_variable_info, file = "nasal_microbiome_variable_info")
  save(nasal_microbiome_sample_info, file = "nasal_microbiome_sample_info")

  save(metabolome_expression_data, file = "metabolome_expression_data")
  save(metabolome_variable_info, file = "metabolome_variable_info")
  save(metabolome_sample_info, file = "metabolome_sample_info")
}

{
  load("nasal_microbiome_expression_data")
  load("nasal_microbiome_variable_info")
  load("nasal_microbiome_sample_info")
  
  load("metabolome_expression_data")
  load("metabolome_variable_info")
  load("metabolome_sample_info")  
}

dim(nasal_microbiome_expression_data)
dim(metabolome_expression_data)

###finally, for nasal microbiome, 107 genus, for metabolome, 169 metabolites

######--------------------------------------------------------------------------
##for raw data, we just log(x+1, 2)
library(plyr)

metabolome_expression_data = 
  log(metabolome_expression_data + 1, 2)

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

# nasal_microbiome_metabolome_lm_adjusted_cor =
#   lm_adjusted_cor(
#     data_set1 = nasal_microbiome_expression_data,
#     data_set2 = metabolome_expression_data,
#     sample_info = nasal_microbiome_sample_info,
#     method = "all",
#     threads = 8
#   )
# 
# nasal_microbiome_metabolome_lm_adjusted_cor_spearman = nasal_microbiome_metabolome_lm_adjusted_cor[[1]]
# 
# save(
#   nasal_microbiome_metabolome_lm_adjusted_cor_spearman,
#   file = "nasal_microbiome_metabolome_lm_adjusted_cor_spearman",
#   compress = "xz"
# )

load("nasal_microbiome_metabolome_lm_adjusted_cor_spearman")

###here we use the lm_adjusted_cor
sum(nasal_microbiome_metabolome_lm_adjusted_cor_spearman$p_adjust < 0.05)
sum(nasal_microbiome_metabolome_lm_adjusted_cor_spearman$p_adjust < 0.2)

cor_data =
  nasal_microbiome_metabolome_lm_adjusted_cor_spearman %>%
  dplyr::filter(p_adjust < 0.2)

###output plot
#####network to show the correlation between microbiome and metabolite
dim(cor_data)

library(ggraph)
library(igraph)
library(tidygraph)

###example network
library(plyr)

edge_data =
  cor_data %>%
  dplyr::mutate(p = -log(p_adjust, 10)) %>%
  dplyr::select(microbiome, metabolite, p, p_adjust, cor) %>%
  dplyr::rename(from = microbiome,
                to = metabolite,
                cor = cor) 
  
node_data =
  data.frame(node = unique(c(edge_data$from, edge_data$to))) %>%
  dplyr::left_join(metabolome_variable_info[, c("variable_id", "Metabolite")],
                   by = c("node" = "variable_id")) %>%
  dplyr::rename(true_name = Metabolite) %>%
  dplyr::mutate(class = case_when(
    is.na(true_name) ~ "Nasal microbiome",
    !is.na(true_name) ~ "Metabolite"
  )) %>%
  dplyr::left_join(nasal_microbiome_variable_info[, c("variable_id", "Genus")],
                   by = c("node" = "variable_id")) %>%
  dplyr::mutate(true_name =
                  case_when(is.na(true_name) ~ Genus,
                            !is.na(true_name) ~ true_name)) %>%
  dplyr::select(-Genus)


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

####output data
edge_data$from_true_name = 
  node_data$true_name[match(edge_data$from, node_data$node)]

edge_data$to_true_name = 
  node_data$true_name[match(edge_data$to, node_data$node)]


edge_data = 
  edge_data %>% 
  dplyr::mutate(significance = case_when(
    p_adjust < 0.05 ~ "p.adj<0.05",
    p_adjust >= 0.05 ~ "0.05<p.adj<0.2",
  ))


library(openxlsx)
wb <- createWorkbook()
modifyBaseFont(wb, fontSize = 12, fontName = "Time New Roma")
addWorksheet(wb, sheetName = "Node data", gridLines = TRUE)
addWorksheet(wb, sheetName = "Edge data", gridLines = TRUE)
freezePane(wb, sheet = 1, firstRow = TRUE, firstCol = TRUE)
freezePane(wb, sheet = 2, firstRow = TRUE, firstCol = TRUE)
writeDataTable(wb, sheet = 1, x = node_data,
               colNames = TRUE, rowNames = FALSE)
writeDataTable(wb, sheet = 2, x = edge_data,
               colNames = TRUE, rowNames = FALSE)
saveWorkbook(wb = wb, file = "nasal_microbiome_vs_metabolome.xlsx", overwrite = TRUE)

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
  geom_edge_link(
    strength = 1,
    aes(
      color = cor,
      width = abs(cor),
      linetype = significance
    ),
    alpha = 1,
    show.legend = TRUE
  ) +
  scale_edge_linetype_manual(values = c("p.adj<0.05" = 1,
                                        "0.05<p.adj<0.2" = 2))+
  geom_node_point(aes(size = Degree,
                      color = class), 
                  shape = 16, 
                  alpha = 1, show.legend = TRUE) +
  geom_node_text(aes(x = x ,
                     y = y ,
                     color = "black",
                     label = true_name),
                 repel = TRUE,
                 size = 2) +
  # shadowtext::geom_shadowtext(aes(x = x, y = y, 
  #                                 # label = ifelse(class == "Metabolite", NA, true_name),
  #                                 label = true_name,
  #                                 color = class), bg.colour = "white") +
  guides(edge_color = guide_edge_colorbar(title = "Correlation"),
         color = guide_legend(title = "Class")) +
  scale_edge_width_continuous(range = c(0.1,2)) +
  scale_size_continuous(range = c(5, 10)) +
  scale_color_manual(values = omics_color[c("Metabolite", "Nasal microbiome")]) +
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
