#' ---
#' title: "Stool microbiome cytokine correlation"
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
dir.create("data_analysis/correlation_network/whole_data_set/stool_microbiome_vs_cytokine")
setwd("data_analysis/correlation_network/whole_data_set/stool_microbiome_vs_cytokine")

####load data
###stool microbiome
{
  load(here::here(
    "data_analysis/stool_microbiome/data_preparation/expression_data"
  ))
  load(here::here(
    "data_analysis/stool_microbiome/data_preparation/sample_info"
  ))
  load(here::here(
    "data_analysis/stool_microbiome/data_preparation/variable_info"
  ))
}

stool_microbiome_expression_data = expression_data
stool_microbiome_sample_info = sample_info
stool_microbiome_variable_info = variable_info

###read genus table
expression_data =
  data.table::fread(here::here("data/from_xin/Genus Table/ST/Genus_ST.csv")) %>%
  as.data.frame() %>%
  tibble::column_to_rownames(var = "SampleID") %>%
  dplyr::select(-c(V1:batch)) %>%
  t() %>%
  as.data.frame()

stool_microbiome_variable_info =
  stool_microbiome_variable_info[match(rownames(expression_data), stool_microbiome_variable_info$Genus),]

stool_microbiome_variable_info$Genus == rownames(expression_data)

rownames(expression_data) = stool_microbiome_variable_info$variable_id

stool_microbiome_expression_data =
  expression_data

stool_microbiome_variable_info =
  stool_microbiome_variable_info %>%
  dplyr::filter(!stringr::str_detect(Genus, "Unclassified_Bacteria"))

stool_microbiome_expression_data =
  stool_microbiome_expression_data[stool_microbiome_variable_info$variable_id,]

dim(stool_microbiome_sample_info)
dim(stool_microbiome_variable_info)

rownames(stool_microbiome_expression_data) == stool_microbiome_variable_info$variable_id
colnames(stool_microbiome_expression_data) == stool_microbiome_sample_info$sample_id

###plasma cytokine
{
  load(here::here(
    "data_analysis/cytokine/data_preparation/expression_data"
  ))
  load(here::here("data_analysis/cytokine/data_preparation/sample_info"))
  load(here::here(
    "data_analysis/cytokine/data_preparation/variable_info"
  ))
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
dim(stool_microbiome_sample_info)
dim(cytokine_sample_info)

length(stool_microbiome_sample_info$subject_id)
length(unique(stool_microbiome_sample_info$subject_id))

###remove CHEX
cytokine_expression_data <- filter(cytokine_expression_data, ! rownames(cytokine_expression_data) %in% c("CHEX1", "CHEX2", "CHEX3", "CHEX4"))
dim(cytokine_expression_data)

###just matched samples according to sample id, only 1 missed
intersect_sample_id =
  intersect(stool_microbiome_sample_info$sample_id,
            cytokine_sample_info$sample_id)

length(intersect_sample_id)

stool_microbiome_expression_data =
  stool_microbiome_expression_data[,intersect_sample_id]

cytokine_expression_data =
  cytokine_expression_data[,intersect_sample_id]

stool_microbiome_sample_info =
  stool_microbiome_sample_info[match(intersect_sample_id, stool_microbiome_sample_info$sample_id),]

cytokine_sample_info =
  cytokine_sample_info[match(intersect_sample_id, cytokine_sample_info$sample_id),]

length(unique(cytokine_sample_info$subject_id))

###only remain the subjects with at least >= 5
remian_subject_id =
  stool_microbiome_sample_info %>%
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

stool_microbiome_sample_info =
  stool_microbiome_sample_info %>%
  dplyr::filter(subject_id %in% remian_subject_id)

stool_microbiome_expression_data =
  stool_microbiome_expression_data[,stool_microbiome_sample_info$sample_id]

##only remain the genus at least in 10% subjects
remain_idx =
  which(rowSums(stool_microbiome_expression_data) > 0)

stool_microbiome_expression_data = stool_microbiome_expression_data[remain_idx,]
stool_microbiome_variable_info = stool_microbiome_variable_info[remain_idx,,drop = FALSE]

remain_idx =
  stool_microbiome_expression_data %>%
  apply(1, function(x){
    sum(as.numeric(x) == 0) / ncol(stool_microbiome_expression_data)
  }) %>%
  `<`(0.9) %>%
  which()

length(remain_idx)

stool_microbiome_expression_data = stool_microbiome_expression_data[remain_idx,]
stool_microbiome_variable_info = stool_microbiome_variable_info[remain_idx,,drop = FALSE]

##save data
{
  save(stool_microbiome_expression_data, file = "stool_microbiome_expression_data")
  save(stool_microbiome_variable_info, file = "stool_microbiome_variable_info")
  save(stool_microbiome_sample_info, file = "stool_microbiome_sample_info")

  save(cytokine_expression_data, file = "cytokine_expression_data")
  save(cytokine_variable_info, file = "cytokine_variable_info")
  save(cytokine_sample_info, file = "cytokine_sample_info")
}

{
  load("stool_microbiome_expression_data")
  load("stool_microbiome_variable_info")
  load("stool_microbiome_sample_info")
  
  load("cytokine_expression_data")
  load("cytokine_variable_info")
  load("cytokine_sample_info")  
}

dim(stool_microbiome_expression_data)
dim(cytokine_expression_data)

###finally, for stool microbiome, 107 genus, for cytokine, 169 metabolites

######--------------------------------------------------------------------------
##for raw data, we just log(x+1, 2)
library(plyr)

plot(density(as.numeric(cytokine_expression_data[1,])))

cytokine_expression_data =
  log(cytokine_expression_data + 1, 2)

library(plyr)

###because our microbiome are percentage data, so here we use the CTL method
library(compositions)
stool_microbiome_expression_data = 
stool_microbiome_expression_data %>% 
  purrr::map(function(x){
    x = compositions::clr(x) %>% 
      as.numeric()
    x
  }) %>% 
  do.call(cbind, .) %>% 
    as.data.frame()

rownames(stool_microbiome_expression_data) = stool_microbiome_variable_info$variable_id

#####
dim(cytokine_expression_data)
dim(stool_microbiome_expression_data)

stool_microbiome_sample_info$subject_id == cytokine_sample_info$subject_id

##https://rpubs.com/DKCH2020/578881
##https://ourcodingclub.github.io/tutorials/mixed-models/
##https://zhuanlan.zhihu.com/p/63092231
##https://www.linglab.cn/knowledge/10

###step 1
###linear mixed model to adjust the subject ID random effect, and then use the partial correlation
###to get the correlation between microbiome and cytokine
library(lme4)
library(rmcorr)

library(future)
library(furrr)

stool_microbiome_cytokine_lm_adjusted_cor =
  lm_adjusted_cor(
    data_set1 = stool_microbiome_expression_data,
    data_set2 = cytokine_expression_data,
    sample_info = stool_microbiome_sample_info,
    method = "all",
    threads = 8
  )

stool_microbiome_cytokine_lm_adjusted_cor_spearman = stool_microbiome_cytokine_lm_adjusted_cor[[1]]

save(
  stool_microbiome_cytokine_lm_adjusted_cor_spearman,
  file = "stool_microbiome_cytokine_lm_adjusted_cor_spearman",
  compress = "xz"
)

load("stool_microbiome_cytokine_lm_adjusted_cor_spearman")


###here we use the lm_adjusted_cor
sum(stool_microbiome_cytokine_lm_adjusted_cor_spearman$p_adjust < 0.2)

cor_data =
  stool_microbiome_cytokine_lm_adjusted_cor_spearman %>%
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
  dplyr::left_join(cytokine_variable_info[, c("variable_id", "variable_id")] %>% 
                     dplyr::rename(protein_name = variable_id.1),
                   by = c("node" = "variable_id")) %>%
  dplyr::rename(true_name = protein_name) %>%
  dplyr::mutate(class = case_when(
    is.na(true_name) ~ "Stool microbiome",
    !is.na(true_name) ~ "cytokine"
  )) %>%
  dplyr::left_join(stool_microbiome_variable_info[, c("variable_id", "Genus")],
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
saveWorkbook(wb = wb, file = "stool_microbiome_vs_cytokine.xlsx", overwrite = TRUE)

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
  #                                 # label = ifelse(class == "cytokine", NA, true_name),
  #                                 label = true_name,
  #                                 color = class), bg.colour = "white") +
  guides(edge_color = guide_edge_colorbar(title = "Correlation"),
         color = guide_legend(title = "Class")) +
  scale_edge_width_continuous(range = c(0.1,2)) +
  scale_size_continuous(range = c(5, 10)) +
  scale_color_manual(values = omics_color[c("Cytokine", "Stool microbiome")]) +
  scale_edge_color_gradient2(low = viridis::cividis(n = 2)[1], 
                             mid = "white", 
                             high = "red") +
  ggraph::theme_graph() +
  theme(plot.background = element_rect(fill = "transparent", color = NA), 
        panel.background = element_rect(fill = "transparent", color = NA),
        legend.position = "right",
        legend.background = element_rect(fill = "transparent", color = NA))

plot

library(extrafont)
extrafont::loadfonts()

ggsave(plot,
       filename = "stool_microbiome_cytokine_cor_plot_example.pdf",
       width = 8.5,
       height = 7)



