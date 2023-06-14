#' ---
#' title: "Stool microbiome metabolome correlation"
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

####load data
###load diversity
load("data/from_xin/Diversity_Datatable.RData")

##prevalence
load("data/from_xin/Prevalance.RData")

###stool microbiome
load("data_analysis/st_microbiome/data_preparation/expression_data")
load("data_analysis/st_microbiome/data_preparation/sample_info")
load("data_analysis/st_microbiome/data_preparation/variable_info")

st_microbiome_expression_data = expression_data
st_microbiome_sample_info = sample_info
st_microbiome_variable_info = variable_info

##remove subjects without sspg
remain_idx = 
  which(!is.na(st_microbiome_sample_info$SSPG))

st_microbiome_expression_data = 
  st_microbiome_expression_data[,remain_idx]

st_microbiome_sample_info = st_microbiome_sample_info[remain_idx,]

###plasma metabolome
load("data_analysis/metabolome/data_preparation/expression_data")
load("data_analysis/metabolome/data_preparation/sample_info")
load("data_analysis/metabolome/data_preparation/variable_info")

metabolome_expression_data = expression_data
metabolome_sample_info = sample_info
metabolome_variable_info = variable_info

metabolome_sample_info$CollectionDate = 
  as.Date(metabolome_sample_info$CollectionDate, "%m/%d/%y")

##remove subjects without sspg
remain_idx = 
  which(!is.na(metabolome_sample_info$SSPG))

metabolome_expression_data = 
  metabolome_expression_data[,remain_idx]

metabolome_sample_info = metabolome_sample_info[remain_idx,]

###remove subjects which have less than 5 samples
remain_subject_id =
  metabolome_sample_info %>%
  dplyr::group_by(subject_id) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::ungroup() %>%
  dplyr::filter(n > 5) %>%
  pull(subject_id)

metabolome_sample_info =
  metabolome_sample_info %>%
  dplyr::filter(subject_id %in% remain_subject_id)

metabolome_expression_data =
  metabolome_expression_data[, metabolome_sample_info$sample_id]

dim(metabolome_expression_data)
colnames(metabolome_expression_data)

length(unique(metabolome_sample_info$subject_id))

#######work directory
setwd(masstools::get_project_wd())
setwd("data_analysis/st_microbiome_metabolome")

plot = 
  metabolome_sample_info %>% 
  dplyr::group_by(subject_id) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::arrange(desc(n)) %>% 
  dplyr::mutate(subject_id = factor(subject_id, levels = subject_id)) %>% 
  ggplot(aes(subject_id, n)) +
  geom_point() +
  theme_bw() +
  labs(x = "", y = "Sample number") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,
                                   size = 5),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 13),
        axis.text.y = element_text(size = 12),
        panel.grid.minor = element_blank())

plot
# ggsave(plot, filename = "metabolome_subject_sample_number.pdf", width = 10, height = 7)

###match samples
dim(st_microbiome_sample_info)
dim(metabolome_sample_info)

###just matched samples according to sample id, only 1 missed
intersect_sample_id = 
  intersect(st_microbiome_sample_info$sample_id,
            metabolome_sample_info$sample_id)

length(intersect_sample_id)

st_microbiome_expression_data = 
  st_microbiome_expression_data[,intersect_sample_id]

metabolome_expression_data = 
  metabolome_expression_data[,intersect_sample_id]

st_microbiome_sample_info = 
  st_microbiome_sample_info[match(intersect_sample_id, st_microbiome_sample_info$sample_id),]  

metabolome_sample_info = 
  metabolome_sample_info[match(intersect_sample_id, metabolome_sample_info$sample_id),]  

plot = 
  metabolome_sample_info %>% 
  dplyr::group_by(subject_id) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::arrange(desc(n)) %>% 
  dplyr::mutate(subject_id = factor(subject_id, levels = subject_id)) %>% 
  ggplot(aes(subject_id, n)) +
  geom_point() +
  theme_bw() +
  labs(x = "", y = "Sample number") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,
                                   size = 10),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 13),
        axis.text.y = element_text(size = 12),
        panel.grid.minor = element_blank())
plot
# ggsave(plot, filename = "metabolome_subject_sample_number2.pdf", width = 10, height = 7)

metabolome_sample_info %>% 
  dplyr::group_by(subject_id) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::ungroup() %>% 
  pull(n) %>% 
  range()

metabolome_sample_info %>% 
  dplyr::group_by(subject_id) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::ungroup() %>% 
  pull(n) %>% 
  mean()

###use the prevalence as the stool microbiome expression data
st_microbiome_expression_data = ST.Pr %>% 
  t() %>% 
  as.data.frame()

st_microbiome_sample_info = 
  sample_info %>% 
  dplyr::distinct(subject_id, .keep_all = TRUE) %>% 
  dplyr::filter(subject_id %in% colnames(st_microbiome_expression_data))

st_microbiome_expression_data = 
  st_microbiome_expression_data[,st_microbiome_sample_info$subject_id]

st_microbiome_variable_info = 
  data.frame(variable_id = rownames(st_microbiome_expression_data))

# ###remove the unclassified variables
# grep("class", st_microbiome_variable_info$variable_id, value = TRUE)
# 
# st_microbiome_variable_info =
#   st_microbiome_variable_info %>% 
#   dplyr::filter(!stringr::str_detect(variable_id, "Unclassified"))
# 
# st_microbiome_expression_data = 
#   st_microbiome_expression_data[st_microbiome_variable_info$variable_id,]

###remove the samples which have no SSPG
st_microbiome_sample_info =
  st_microbiome_sample_info %>%
  dplyr::filter(!is.na(SSPG)) %>%
  dplyr::filter(subject_id %in% metabolome_sample_info$subject_id)

st_microbiome_expression_data = 
  st_microbiome_expression_data[,st_microbiome_sample_info$subject_id]

dim(st_microbiome_sample_info)
dim(st_microbiome_expression_data)

sum(rowSums(st_microbiome_expression_data) == 0)
##only remain the genus at least in 1% subjects
remain_idx = 
  which(rowSums(st_microbiome_expression_data) > 0)

st_microbiome_expression_data = st_microbiome_expression_data[remain_idx,]
st_microbiome_variable_info = st_microbiome_variable_info[remain_idx,,drop = FALSE]

remain_idx = 
  st_microbiome_expression_data %>% 
  apply(1, function(x){
    sum(as.numeric(x) == 0) / ncol(st_microbiome_expression_data)
  }) %>% 
  `<`(0.5) %>% 
  which()

st_microbiome_expression_data = st_microbiome_expression_data[remain_idx,]
st_microbiome_variable_info = st_microbiome_variable_info[remain_idx,,drop = FALSE]

# ##remove the genus which only > 0 in one sample
# dim(st_microbiome_expression_data)
# 
# remain_idx = 
# apply(st_microbiome_expression_data, 1, function(x){
#   sum(x == 0)
# }) %>% 
#   `<`(34) %>% 
#   which()
# length(remain_idx)
# st_microbiome_expression_data = st_microbiome_expression_data[remain_idx,]
# st_microbiome_variable_info = st_microbiome_variable_info[remain_idx,,drop = FALSE]

###calculate the rsd
rsd =
  apply(st_microbiome_expression_data, 1, function(x) {
    sd(x) * 100 / mean(x)
  })
plot(rsd)
which(rsd < 30)
plot(as.numeric(st_microbiome_expression_data[38,]))

library(ComplexHeatmap)
library(circlize)

col_fun = circlize::colorRamp2(breaks = c(0, 1), colors = c("white", "red"))

plot = 
Heatmap(
  t(st_microbiome_expression_data),
  show_row_names = TRUE,
  show_column_names = FALSE, 
  col = col_fun, border = TRUE, 
  name = "Prevalence", 
  row_names_gp = grid::gpar(fontsize = 8),
  column_names_gp = grid::gpar(fontsize = 8, angle = 45)
)

plot = ggplotify::as.ggplot(plot)
plot
# ggsave(plot, filename = "st_microbiome_prevalence.pdf", width = 14, height = 5)


######--------------------------------------------------------------------------
###for metabolome data, we just use the mean value for each metabolite across all sample for each object
##for raw data, we just log(x+1, 2)
library(plyr)

metabolome_expression_data = 
  log(metabolome_expression_data + 1, 2)

# ##scaled
# metabolome_expression_data = 
#   metabolome_expression_data %>% 
#   apply(1, function(x){
#     (x - mean(x))/sd(x)
#   }) %>% 
#   t() %>% 
#   as.data.frame()

metabolome_expression_data =
  metabolome_expression_data %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "sample_id") %>%
  dplyr::left_join(metabolome_sample_info[, c("subject_id", "sample_id")],
                   by = "sample_id") %>%
  plyr::dlply(.variables = .(subject_id)) %>%
  purrr::map(function(x) {
    x = x %>%
      dplyr::select(-c(subject_id, sample_id))
    apply(x, 2, function(y) {
      mean(y)
    })
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()

dim(metabolome_expression_data)  

rownames(metabolome_expression_data)

metabolome_expression_data 

##match st microbiome and plasma metabolome
intersect_subject_id = 
intersect(st_microbiome_sample_info$subject_id, 
          metabolome_sample_info$subject_id)

length(st_microbiome_sample_info$subject_id)
length(metabolome_sample_info$subject_id)

st_microbiome_sample_info = 
st_microbiome_sample_info %>% 
  dplyr::filter(subject_id %in% intersect_subject_id)

st_microbiome_expression_data = 
  st_microbiome_expression_data[,st_microbiome_sample_info$subject_id]

colnames(st_microbiome_expression_data) == rownames(metabolome_expression_data)

metabolome_expression_data =
  metabolome_expression_data[st_microbiome_sample_info$subject_id,] %>%
  t() %>%
  as.data.frame()

metabolome_sample_info =
  metabolome_sample_info[match(colnames(metabolome_expression_data), metabolome_sample_info$subject_id),]

colnames(st_microbiome_expression_data) == colnames(metabolome_expression_data)

#####
dim(metabolome_expression_data)
dim(st_microbiome_expression_data)

rsd =
  metabolome_expression_data %>%
  apply(1, function(x) {
    sd(x) * 100 / mean(x)
  }) 

plot(rsd)

rsd

which(rsd < 3)

plot(as.numeric(metabolome_expression_data[154,]))

st_microbiome_sample_info$subject_id == metabolome_sample_info$subject_id

st_microbiome_sample_info$SSPG
metabolome_sample_info$SSPG



# ###step 1
# ###liner model to find the association between microbiome core microbiome and SSPG
# st_microbiome_sspg_cor =
#   st_microbiome_expression_data %>%
#   t() %>%
#   as.data.frame() %>%
#   purrr::map(function(x) {
#     temp_data =
#       data.frame(x = x,
#                  st_microbiome_sample_info) %>%
#       dplyr::filter(!is.na(SSPG)) %>%
#       dplyr::mutate(SSPG = as.numeric(SSPG),
#                     x = as.numeric(x))
# 
#     temp_data$Gender[temp_data$Gender == 'F'] = 0
#     temp_data$Gender[temp_data$Gender == 'M'] = 1
#     temp_data$Gender = as.numeric(temp_data$Gender)
# 
#     temp_data$Ethnicity[temp_data$Ethnicity == 'C'] = 0
#     temp_data$Ethnicity[temp_data$Ethnicity == 'H'] = 1
#     temp_data$Ethnicity[temp_data$Ethnicity == 'B'] = 2
#     temp_data$Ethnicity[temp_data$Ethnicity == 'A'] = 3
#     temp_data$Ethnicity = as.numeric(temp_data$Ethnicity)
# 
#     glm_reg =
#       glm(SSPG ~ x + Gender + Adj.age + Ethnicity,
#           family = gaussian,
#           temp_data)
# 
#     temp =
#       summary(glm_reg)$coefficients %>%
#       as.data.frame()
# 
#     c(estimate = temp$Estimate[2],
#       glm_p = temp$`Pr(>|t|)`[2])
# 
#     library(ppcor)
# 
#     # cor_value =
#     #   cor.test(temp_data$x, temp_data$SSPG)
#     cor_value =
#       pcor.test(
#         x = temp_data$x,
#         y = temp_data$SSPG,
#         z = temp_data[, c("Gender", "Adj.age", "Ethnicity")],
#         method = "spearman"
#       )
# 
#     result =
#       c(
#         estimate = unname(temp$Estimate[2]),
#         glm_p = unname(temp$`Pr(>|t|)`[2]),
#         cor = unname(cor_value$estimate),
#         cor_p = unname(cor_value$p.value)
#       )
# 
#     if (is.na(result[3])) {
#       result[3] = 0
#     }
# 
#     if (is.na(result[4])) {
#       result[4] = 1
#     }
# 
#     result
# 
#   }) %>%
#   do.call(rbind, .) %>%
#   as.data.frame() %>%
#   tibble::rownames_to_column(var = "microbiome") %>%
#   dplyr::mutate(phenotype = "SSPG") %>%
#   dplyr::select(microbiome, phenotype, dplyr::everything())
# 
# ###remove the variables with low correlation
# plot(st_microbiome_sspg_cor$cor)
# 
# st_microbiome_sspg_cor =
# st_microbiome_sspg_cor %>%
#   dplyr::filter(abs(cor) > 0.2)
# 
# st_microbiome_sspg_cor$glm_p_adjust = p.adjust(st_microbiome_sspg_cor$glm_p, method = "fdr")
# st_microbiome_sspg_cor$cor_p_adjust = p.adjust(st_microbiome_sspg_cor$cor_p, method = "fdr")
# 
# save(st_microbiome_sspg_cor, file = "st_microbiome_sspg_cor")
load("st_microbiome_sspg_cor")

which(st_microbiome_sspg_cor$glm_p_adjust < 0.05)
which(st_microbiome_sspg_cor$cor_p_adjust < 0.05)
which(st_microbiome_sspg_cor$glm_p < 0.05)
which(st_microbiome_sspg_cor$cor_p < 0.05)

plot(st_microbiome_sspg_cor$cor)

idx =
  which(st_microbiome_sspg_cor$glm_p_adjust < 0.05 
          # st_microbiome_sspg_cor$cor_p_adjust < 0.05
          )
st_microbiome_sspg_cor[idx, ]


###output plot
dir("st_microbiome_sspg")
# for(i in idx){
#   cat(i, " ")
#   id1 = st_microbiome_sspg_cor$microbiome[i]
#   id2 = st_microbiome_sspg_cor$phenotype[i]
#   cor = st_microbiome_sspg_cor$cor[i]
#   cor_p_adjust = st_microbiome_sspg_cor$cor_p_adjust[i]
#   title = paste("Correlation: ", round(cor,2), "/p.adjust: ",round(cor_p_adjust, 4),
#                 sep = "")
#   plot =
#   data.frame(
#     x = as.numeric(st_microbiome_expression_data[id1,]),
#     sspg = as.numeric(st_microbiome_sample_info$SSPG)
#   ) %>%
#     ggplot(aes(x*100, sspg)) +
#     geom_point(size = 3) +
#     geom_smooth(method = "glm", color = "red") +
#     theme_bw() +
#     labs(x = paste("Prevalence % (",id1, ")", sep = ""),
#          y = "SSPG",
#          title = title) +
#     base_theme
#   ggsave(plot, filename = file.path("st_microbiome_sspg", paste(id1, ".pdf", sep = "")),
#          width = 9, height = 7)
# }


# ###--------------------------------------------------------------------------
# ####liner model to find the association between st_microbiome and metabolite
# st_microbiome_metabolome_cor =
#   st_microbiome_expression_data %>%
#   t() %>%
#   as.data.frame() %>%
#   purrr::map(function(x) {
#     temp_cor =
#     metabolome_expression_data %>%
#       t() %>%
#       as.data.frame() %>%
#       purrr::map(function(y) {
#         temp_data =
#           data.frame(x = x,
#                      y = y,
#                      metabolome_sample_info)
# 
#             temp_data$Gender[temp_data$Gender == 'F'] = 0
#             temp_data$Gender[temp_data$Gender == 'M'] = 1
#             temp_data$Gender = as.numeric(temp_data$Gender)
# 
#             temp_data$Ethnicity[temp_data$Ethnicity == 'C'] = 0
#             temp_data$Ethnicity[temp_data$Ethnicity == 'H'] = 1
#             temp_data$Ethnicity[temp_data$Ethnicity == 'B'] = 2
#             temp_data$Ethnicity[temp_data$Ethnicity == 'A'] = 3
#             temp_data$Ethnicity = as.numeric(temp_data$Ethnicity)
# 
#             ###y is metabolite and x is microbiome
#             glm_reg =
#               glm(y ~ x + Gender + Adj.age + Ethnicity,
#                   family = gaussian,
#                   temp_data)
# 
#         temp =
#           summary(glm_reg)$coefficients %>%
#           as.data.frame()
# 
#         # cor_value =
#         #   cor.test(temp_data$x, temp_data$y)
# 
#         cor_value =
#           pcor.test(
#             x = temp_data$x,
#             y = temp_data$y,
#             z = temp_data[, c("Gender", "Adj.age", "Ethnicity")],
#             method = "spearman"
#           )
# 
#         result =
#           c(
#             estimate = unname(temp$Estimate[2]),
#             glm_p = unname(temp$`Pr(>|t|)`[2]),
#             cor = unname(cor_value$estimate),
#             cor_p = unname(cor_value$p.value)
#           )
#         if (is.na(result[3])) {
#           result[3] = 0
#         }
#         if (is.na(result[4])) {
#           result[4] = 1
#         }
#         result
#       }) %>%
#       do.call(rbind, .) %>%
#       as.data.frame() %>%
#       tibble::rownames_to_column(var = "metabolite")
# 
#     temp_cor =
#     temp_cor %>%
#       dplyr::filter(abs(cor) > 0.2)
# 
#     temp_cor$glm_p_adjust = p.adjust(temp_cor$glm_p, method = "fdr")
#     temp_cor$cor_p_adjust = p.adjust(temp_cor$cor_p, method = "fdr")
#     temp_cor
#   })
# 
# st_microbiome_metabolome_cor =
#   purrr::map2(
#     .x = st_microbiome_metabolome_cor,
#     .y = names(st_microbiome_metabolome_cor),
#     .f = function(x, y) {
#       data.frame(microbiome = y, x, stringsAsFactors = FALSE)
#     }
#   ) %>%
#   do.call(rbind, .) %>%
#   as.data.frame()
# 
# save(st_microbiome_metabolome_cor, file = "st_microbiome_metabolome_cor")

load("st_microbiome_metabolome_cor")

st_microbiome_metabolome_cor$cor
sum(st_microbiome_metabolome_cor$cor_p_adjust < 0.05)
sum(st_microbiome_metabolome_cor$glm_p_adjust < 0.05)

idx = 
which(
  st_microbiome_metabolome_cor$cor_p_adjust < 0.05 
    # st_microbiome_metabolome_cor$glm_p_adjust < 0.05
)


dir.create("st_microbiome_metabolome")
# for(i in idx){
#   cat(i, " ")
#   id1 = st_microbiome_metabolome_cor$microbiome[i]
#   id2 = st_microbiome_metabolome_cor$metabolite[i]
#   cor = st_microbiome_metabolome_cor$cor[i]
#   cor_p_adjust = st_microbiome_metabolome_cor$cor_p_adjust[i]
#   title = paste("Correlation: ", round(cor,2), "/p.adjust: ",round(cor_p_adjust, 4),
#                 sep = "")
#   plot =
#     data.frame(
#       x = as.numeric(st_microbiome_expression_data[id1,]),
#       y = as.numeric(metabolome_expression_data[id2,])
#     ) %>%
#     ggplot(aes(x*100, y)) +
#     geom_point(size = 3) +
#     geom_smooth(method = "glm", color = "red") +
#     theme_bw() +
#     labs(x = paste("Prevalence % (",id1, ")", sep = ""),
#          y = metabolome_variable_info$Metabolite[match(id2, metabolome_variable_info$variable_id)],
#          title = title) +
#     base_theme
#   
#   name = metabolome_variable_info$Metabolite[match(id2, metabolome_variable_info$variable_id)] %>% 
#     stringr::str_replace("/", "_")
#   
#   ggsave(plot, filename = file.path("st_microbiome_metabolome",
#                                     paste(id1, "_", name,
#                                           ".pdf", sep = "")),
#          width = 9, height = 7)
# }


#####network to show the correlation between microbiome and metabolite
st_microbiome_metabolome_cor

library(ggraph)
library(igraph)
library(tidygraph)

edge_data = 
  st_microbiome_metabolome_cor %>% 
  dplyr::filter(glm_p_adjust < 0.05) %>% 
  dplyr::mutate(p = -log(glm_p_adjust, 10)) %>% 
  dplyr::select(microbiome, metabolite, p,cor) %>% 
  dplyr::rename(from = microbiome, 
                to = metabolite, 
                cor = cor)

node_data = 
  data.frame(node = unique(c(edge_data$from, edge_data$to))) %>% 
  dplyr::left_join(metabolome_variable_info[,c("variable_id", "Metabolite")],
                   by = c("node" = "variable_id")) %>% 
  dplyr::rename(true_name = Metabolite) %>% 
    dplyr::mutate(class = case_when(
      is.na(true_name) ~ "Stool microbiome",
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
         layout = 'linear',
         circular = TRUE) +
  geom_edge_arc(strength = 1,
                aes(color = cor),
                alpha = 1,
                show.legend = TRUE) +
  geom_node_point(aes(size = Degree,
                      color = class), 
                  shape = 16, 
                  alpha = 0.5, show.legend = TRUE) +
  geom_node_text(aes(x = x * 1.05,
                     y = y * 1.05,
                     color = class,
                     label = true_name, 
                     angle = angle, 
                     hjust = hjust),
                 repel = FALSE, size = 2) +
  guides(edge_color = guide_edge_colorbar(title = "Correlation"),
         color = guide_legend(title = "Class")) +
  scale_size_continuous(range = c(3, 10)) +
  scale_color_manual(values = omics_color) +
  scale_edge_color_gradient2(low = viridis::cividis(n = 2)[1], 
                             mid = "white", 
                             high = viridis::cividis(n = 2)[2]) +
  ggraph::theme_graph() +
  theme(plot.background = element_rect(fill = "transparent", color = NA), 
        panel.background = element_rect(fill = "transparent", color = NA),
        legend.position = "right",
        legend.background = element_rect(fill = "transparent", color = NA)) +
  coord_cartesian(xlim=c(-1.4,1.4), ylim=c(-1.4,1.4))

plot

ggsave(plot, filename = "st_microbiome_metabolome_cor_plot.pdf", width = 8.5, height = 7)


#####---------------------------------------------------------------------------
###mediation analysis
dim(st_microbiome_expression_data)
dim(metabolome_expression_data)

colnames(st_microbiome_expression_data) == colnames(metabolome_expression_data)

library(mediation)

st_microbiome_sspg_cor2 = 
st_microbiome_sspg_cor %>% 
  dplyr::filter(glm_p_adjust < 0.05 & cor_p_adjust < 0.05)

st_microbiome_metabolome_cor2 = 
  st_microbiome_metabolome_cor %>% 
  dplyr::filter(glm_p_adjust < 0.05 & cor_p_adjust < 0.05) %>% 
  dplyr::filter(microbiome %in% st_microbiome_sspg_cor2$microbiome)

# mediation_result = vector(mode = "list", length = nrow(st_microbiome_metabolome_cor2))
# 
# for(i in 1:nrow(st_microbiome_metabolome_cor2)) {
#   cat(i, "\n")
#   id1 = st_microbiome_metabolome_cor2$microbiome[i]
#   id2 = st_microbiome_metabolome_cor2$metabolite[i]
#   temp_data =
#   data.frame(x = as.numeric(st_microbiome_expression_data[id1,]),
#              y = as.numeric(metabolome_expression_data[id2,]),
#              st_microbiome_sample_info
#              )
# 
#   temp_data$Gender[temp_data$Gender == 'F'] = 0
#   temp_data$Gender[temp_data$Gender == 'M'] = 1
#   temp_data$Gender = as.numeric(temp_data$Gender)
# 
#   temp_data$Ethnicity[temp_data$Ethnicity == 'C'] = 0
#   temp_data$Ethnicity[temp_data$Ethnicity == 'H'] = 1
#   temp_data$Ethnicity[temp_data$Ethnicity == 'B'] = 2
#   temp_data$Ethnicity[temp_data$Ethnicity == 'A'] = 3
#   temp_data$Ethnicity = as.numeric(temp_data$Ethnicity)
#   temp_data$SSPG = as.numeric(temp_data$SSPG)
# 
#   ###y is metabolite and x is microbiome
#   glm_reg1 =
#     glm(y ~ x + Gender + Adj.age + Ethnicity,
#         family = gaussian,
#         temp_data)
# 
#   ###SSPG is SSPG and x is microbiome, y is metabolome
#   glm_reg2 =
#     glm(SSPG ~ x + y + Gender + Adj.age + Ethnicity,
#         family = gaussian,
#         temp_data)
# 
#     result = mediate(
#       glm_reg1 ,
#       glm_reg2,
#       treat = "x",
#       mediator = "y",
#       boot = FALSE,
#       sims = 1000
#     )
# 
#   return_result = list(id1 = id1,
#                        id2 = id2,
#                        id3 = "SSPG",
#                        result = result)
#   mediation_result[[i]] = return_result
# }
# 
# save(mediation_result, file = "mediation_result")
load("mediation_result")
length(mediation_result)

mediate =
  purrr::map(mediation_result, function(x) {
    phenotype = x$id3
    mediator = x$id2
    treat = x$id1
    
    results = summary(x$result)
    acme = results$d1
    acme_ci_lower = results$d1.ci[1]
    acme_ci_upper = results$d1.ci[2]
    acme_p = results$d1.p
    c(phenotype = phenotype,
      mediator = mediator,
      treat = treat,
      acme = acme,
      acme_ci_lower = acme_ci_lower,
      acme_ci_upper = acme_ci_upper,
      acme_p = acme_p)
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()

mediate$acme = as.numeric(mediate$acme)
mediate$`acme_ci_lower.2.5%` = as.numeric(mediate$`acme_ci_lower.2.5%`)
mediate$`acme_ci_upper.97.5%` = as.numeric(mediate$`acme_ci_upper.97.5%`)
.$acme_p = as.numeric(mediate$acme_p)

# mediate =
# mediate %>%
#   dplyr::filter(acme_p < 0.05 & acme > 0 & acme < 1 &
#                   `acme_ci_upper.97.5%` < 1 & `acme_ci_lower.2.5%` > 0)

save(mediate, file = "mediate")

load("mediate")

###get the correlation between treat with mediator, mediator and phenotype and treat with phenotype
cor_value =
t(mediate) %>%
  as.data.frame() %>%
purrr::map(function(x){
 ##treat vs mediator
   cor_treat_mediator_test =
  cor.test(x = as.numeric(exposome_expression_data[x[3],]),
           y = as.numeric(internal_ome_expression_data[x[2],]),
           method = "spearman")
  cor_treat_mediator = cor_treat_mediator_test$estimate
  cor_treat_mediator_p = cor_treat_mediator_test$p.value

  ##mediator vs phenotype
  cor_mediator_phenotype_test =
    cor.test(x = as.numeric(exposome_sample_info[,x[1]]),
             y = as.numeric(internal_ome_expression_data[x[2],]),
             method = "spearman")
  cor_mediator_phenotype = cor_mediator_phenotype_test$estimate
  cor_mediator_phenotype_p = cor_mediator_phenotype_test$p.value

  ##treat vs phenotype
  cor_treat_phenotype_test =
    cor.test(x = as.numeric(exposome_sample_info[,x[1]]),
             y = as.numeric(exposome_expression_data[x[3],]),
             method = "spearman")
  cor_treat_phenotype = cor_treat_phenotype_test$estimate
  cor_treat_phenotype_p = cor_treat_phenotype_test$p.value

  c(cor_treat_mediator = cor_treat_mediator,
    cor_treat_mediator_p = cor_treat_mediator_p,
    cor_mediator_phenotype = cor_mediator_phenotype,
    cor_mediator_phenotype_p = cor_mediator_phenotype_p,
    cor_treat_phenotype = cor_treat_phenotype,
    cor_treat_phenotype_p = cor_treat_phenotype_p
    )
}) %>%
  do.call(rbind, .) %>%
  as.data.frame()

mediate_result =
  data.frame(mediate, cor_value)

##only remain the relations with significant cor
mediate_result =
mediate_result %>%
  dplyr::filter(cor_treat_mediator_p < 0.05 &
                  cor_mediator_phenotype_p < 0.05 &
                  cor_treat_phenotype_p < 0.05) %>%
  dplyr::filter(((cor_treat_mediator.rho * cor_mediator_phenotype.rho) > 0 & cor_treat_phenotype.rho > 0) |
                  ((cor_treat_mediator.rho * cor_mediator_phenotype.rho) < 0 & cor_treat_phenotype.rho < 0)) %>%
  dplyr::arrange(desc(acme))

save(mediate_result, file = "mediate_result")
