#' ---
#' title: "Nasal microbiome proteome correlation"
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

setwd("~/Library/CloudStorage/Box-Box/human_microbiome_project/")
library(tidyverse)
rm(list = ls())

source("code/tools.R")

######work directory
masstools::setwd_project()
dir.create("data_analysis/correlation_network/IS/nasal_microbiome_vs_phenotype_sample_wise")
setwd("data_analysis/correlation_network/IS/nasal_microbiome_vs_phenotype_sample_wise")

###load data
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
  data.table::fread(here::here("data/from_xin/Genus Table/NS/Genus_NS.csv")) %>%
  as.data.frame() %>%
  tibble::column_to_rownames(var = "SampleID") %>%
  dplyr::select(-c(V1:batch)) %>%
  t() %>%
  as.data.frame()

nasal_microbiome_variable_info =
  nasal_microbiome_variable_info[match(rownames(expression_data), nasal_microbiome_variable_info$Genus),]

nasal_microbiome_variable_info$Genus == rownames(expression_data)

###remove the variables which Genus are NA
remove_idx = which(is.na(nasal_microbiome_variable_info$Genus))
remove_idx
if(length(remove_idx) > 0){
  nasal_microbiome_variable_info = nasal_microbiome_variable_info[-remove_idx,]
  expression_data = expression_data[-remove_idx,]
}

rownames(expression_data) = nasal_microbiome_variable_info$variable_id

nasal_microbiome_expression_data =
  expression_data

nasal_microbiome_variable_info =
nasal_microbiome_variable_info %>%
  dplyr::filter(!stringr::str_detect(Genus, "Unclassified_Bacteria"))

nasal_microbiome_expression_data =
  nasal_microbiome_expression_data[nasal_microbiome_variable_info$variable_id,]

dim(nasal_microbiome_sample_info)
dim(nasal_microbiome_variable_info)

rownames(nasal_microbiome_expression_data) == nasal_microbiome_variable_info$variable_id
colnames(nasal_microbiome_expression_data) == nasal_microbiome_sample_info$sample_id
###phenotype
{
  load(here::here(
    "data_analysis/phenotype_data_sample_wise/data_preparation/expression_data"
  ))
  load(here::here("data_analysis/phenotype_data_sample_wise/data_preparation/sample_info"))
  load(here::here(
    "data_analysis/phenotype_data_sample_wise/data_preparation/variable_info"
  ))
}

phenotype_expression_data = expression_data
phenotype_sample_info = sample_info
phenotype_variable_info = variable_info

phenotype_sample_info$colletion_data =
  as.Date(phenotype_sample_info$colletion_data, "%m/%d/%y")

dim(phenotype_expression_data)
length(unique(phenotype_sample_info$subject_id))

#select IS sample first
subject_info <- 
  readr::read_csv("../../../../Figures/Figure6/metadata/metadata.subject.csv") %>% 
  dplyr::select(SubjectID, IRIS) %>% 
  dplyr::rename(subject_id = SubjectID) %>% 
  dplyr::filter(IRIS == "IS")

phenotype_sample_info <- 
  phenotype_sample_info %>% 
  dplyr::filter(subject_id %in% subject_info$subject_id)

nasal_microbiome_sample_info <- 
  nasal_microbiome_sample_info %>% 
  dplyr::filter(subject_id %in% subject_info$subject_id)

phenotype_expression_data <-
  phenotype_expression_data[,phenotype_sample_info$sample_id]  

nasal_microbiome_expression_data <-
  nasal_microbiome_expression_data[,nasal_microbiome_sample_info$sample_id]  

dim(phenotype_expression_data)
dim(nasal_microbiome_expression_data)

###match samples
dim(nasal_microbiome_sample_info)
dim(phenotype_sample_info)

length(nasal_microbiome_sample_info$subject_id)
length(unique(nasal_microbiome_sample_info$subject_id))

###just matched samples according to sample id, only 1 missed
intersect_sample_id =
  intersect(nasal_microbiome_sample_info$sample_id,
            phenotype_sample_info$sample_id)

length(intersect_sample_id)

nasal_microbiome_expression_data =
  nasal_microbiome_expression_data[,intersect_sample_id]

phenotype_expression_data =
  phenotype_expression_data[,intersect_sample_id]

nasal_microbiome_sample_info =
  nasal_microbiome_sample_info[match(intersect_sample_id, nasal_microbiome_sample_info$sample_id),]

phenotype_sample_info =
  phenotype_sample_info[match(intersect_sample_id, phenotype_sample_info$sample_id),]

length(unique(phenotype_sample_info$subject_id))

###only remain the subjects with at least >= 5
remian_subject_id =
  nasal_microbiome_sample_info %>%
  dplyr::group_by(subject_id) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::ungroup() %>%
  dplyr::filter(n > 5) %>%
  dplyr::pull(subject_id)

phenotype_sample_info =
  phenotype_sample_info %>%
  dplyr::filter(subject_id %in% remian_subject_id)

phenotype_expression_data =
  phenotype_expression_data[, phenotype_sample_info$sample_id]

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

  save(phenotype_expression_data, file = "phenotype_expression_data")
  save(phenotype_variable_info, file = "phenotype_variable_info")
  save(phenotype_sample_info, file = "phenotype_sample_info")
}

{
  load("nasal_microbiome_expression_data")
  load("nasal_microbiome_variable_info")
  load("nasal_microbiome_sample_info")
  
  load("phenotype_expression_data")
  load("phenotype_variable_info")
  load("phenotype_sample_info")  
}

dim(nasal_microbiome_expression_data)
dim(phenotype_expression_data)

#####remove phenotype who have a lot of zeros
na_percentage <-
  phenotype_expression_data %>% 
  apply(1, function(x){
    sum(is.na(x))/ncol(phenotype_expression_data)
  })

remain_idx <- unname(which(na_percentage < 0.5))

phenotype_expression_data <-
  phenotype_expression_data[remain_idx,]

phenotype_variable_info <-
  phenotype_variable_info[remain_idx,,drop = FALSE]

###finally, for nasal microbiome, 108 genus, for phenotype, 52 variables

######--------------------------------------------------------------------------
##for raw data, we just log(x+1, 2)
library(plyr)

# phenotype_expression_data = 
#   log(phenotype_expression_data + 1, 2)

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
dim(phenotype_expression_data)
dim(nasal_microbiome_expression_data)

nasal_microbiome_sample_info$subject_id == phenotype_sample_info$subject_id

##https://rpubs.com/DKCH2020/578881
##https://ourcodingclub.github.io/tutorials/mixed-models/
##https://zhuanlan.zhihu.com/p/63092231
##https://www.linglab.cn/knowledge/10

###step 1
###linear mixed model to adjust the subject ID random effect, and then use the partial correlation
###to get the correlation between microbiome and phenotype
library(lme4)
library(rmcorr)

library(future)
library(furrr)
 
nasal_microbiome_phenotype_lm_adjusted_cor <-
  lm_adjusted_cor(
    data_set1 = nasal_microbiome_expression_data,
    data_set2 = phenotype_expression_data,
    sample_info = nasal_microbiome_sample_info,
    method = "all",
    threads = 8
  )

nasal_microbiome_phenotype_lm_adjusted_cor_spearman = nasal_microbiome_phenotype_lm_adjusted_cor[[1]]

save(
  nasal_microbiome_phenotype_lm_adjusted_cor_spearman,
  file = "nasal_microbiome_phenotype_lm_adjusted_cor_spearman",
  compress = "xz"
)

load("nasal_microbiome_phenotype_lm_adjusted_cor_spearman")

###here we use the lm_adjusted_cor
sum(nasal_microbiome_phenotype_lm_adjusted_cor_spearman$p_adjust < 0.2)

cor_data =
  nasal_microbiome_phenotype_lm_adjusted_cor_spearman %>%
  dplyr::filter(p_adjust < 0.2)

###output plot
#####network to show the correlation between microbiome and metabolite
dim(cor_data)

library(ggraph)
library(igraph)
library(tidygraph)

###example network
library(plyr)

edge_data <-
  cor_data %>%
  dplyr::mutate(p = -log(p_adjust, 10)) %>%
  dplyr::select(microbiome, metabolite, p, p_adjust, cor) %>%
  dplyr::rename(from = microbiome,
                to = metabolite,
                cor = cor) 
  
node_data =
  data.frame(node = unique(c(edge_data$from, edge_data$to))) %>%
  dplyr::left_join(phenotype_variable_info[,c("variable_id", "variable_id")] %>% 
                     dplyr::rename(name2 = variable_id.1),
                   by = c("node" = "variable_id")) %>%
  dplyr::rename(true_name = name2) %>%
  dplyr::mutate(class = case_when(
    is.na(true_name) ~ "Nasal microbiome",
    !is.na(true_name) ~ "Phenotype"
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
saveWorkbook(wb = wb, file = "nasal_microbiome_vs_phenotype.xlsx", overwrite = TRUE)

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
  #                                 # label = ifelse(class == "phenotype", NA, true_name),
  #                                 label = true_name,
  #                                 color = class), bg.colour = "white") +
  guides(edge_color = guide_edge_colorbar(title = "Correlation"),
         color = guide_legend(title = "Class")) +
  scale_edge_width_continuous(range = c(0.1,2)) +
  scale_size_continuous(range = c(5, 10)) +
  scale_color_manual(values = c(omics_color[c("Nasal microbiome")], 
                                phenotype_color)) +
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
       filename = "nasal_microbiome_phenotype_cor_plot_example.pdf",
       width = 8.5,
       height = 7)



