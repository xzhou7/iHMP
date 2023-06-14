#' ---
#' title: "skin microbiome metabolome correlation"
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
load(here::here("data_analysis/st_microbiome_vs_metabolome/st_microbiome_metabolome_cor"))
load(here::here("data_analysis/skin_microbiome_vs_metabolome/skin_microbiome_metabolome_cor"))
load(here::here("data_analysis/oral_microbiome_vs_metabolome/oral_microbiome_metabolome_cor"))
load(here::here("data_analysis/nasal_microbiome_vs_metabolome/nasal_microbiome_metabolome_cor"))

st_microbiome_metabolome_cor$cor_p_adjust2 = 
  p.adjust(st_microbiome_metabolome_cor$cor_p, method = "BH")

skin_microbiome_metabolome_cor$cor_p_adjust2 = 
  p.adjust(skin_microbiome_metabolome_cor$cor_p, method = "BH")

oral_microbiome_metabolome_cor$cor_p_adjust2 = 
  p.adjust(oral_microbiome_metabolome_cor$cor_p, method = "BH")

nasal_microbiome_metabolome_cor$cor_p_adjust2 = 
  p.adjust(nasal_microbiome_metabolome_cor$cor_p, method = "BH")

plot(density(abs(st_microbiome_metabolome_cor$cor)), col = "black")
lines(density(abs(skin_microbiome_metabolome_cor$cor)), col = "red")

dim(st_microbiome_metabolome_cor)
dim(skin_microbiome_metabolome_cor)

###this is the number of microbiome
st_number1 = length(unique(st_microbiome_metabolome_cor$microbiome))
skin_number1 = length(unique(skin_microbiome_metabolome_cor$microbiome))
nasal_number1 = length(unique(nasal_microbiome_metabolome_cor$microbiome))
oral_number1 = length(unique(oral_microbiome_metabolome_cor$microbiome))

unique(skin_microbiome_metabolome_cor$microbiome)
unique(nasal_microbiome_metabolome_cor$microbiome)
unique(oral_microbiome_metabolome_cor$microbiome)

library(ComplexUpset)

temp_data = 
rbind(
  data.frame(genus = unique(st_microbiome_metabolome_cor$microbiome), class = "Stool"),
  data.frame(genus = unique(skin_microbiome_metabolome_cor$microbiome), class = "Skin"),
  data.frame(genus = unique(nasal_microbiome_metabolome_cor$microbiome), class = "Nasal"),
  data.frame(genus = unique(oral_microbiome_metabolome_cor$microbiome), class = "Oral")
) %>% 
  tidyr::pivot_wider(names_from = class, values_from = class) %>% 
  dplyr::mutate(Stool = case_when(
    is.na(Stool) ~ FALSE,
    !is.na(Stool) ~ TRUE
  )) %>% 
  dplyr::mutate(Skin = case_when(
    is.na(Skin) ~ FALSE,
    !is.na(Skin) ~ TRUE
  )) %>% 
  dplyr::mutate(Nasal = case_when(
    is.na(Nasal) ~ FALSE,
    !is.na(Nasal) ~ TRUE
  )) %>% 
  dplyr::mutate(Oral = case_when(
    is.na(Oral) ~ FALSE,
    !is.na(Oral) ~ TRUE
  ))

plot = 
upset(data = temp_data, intersect = c("Stool", "Skin", "Nasal", "Oral"))
plot
setwd("data_analysis/microbiome_vs_metabolome")

# ggsave(plot, filename = "body_site_microbe_upset.pdf", width = 9, height = 7)














#####we only remain the correlation with adjusted p < 0.05
st_microbiome_metabolome_cor =
  st_microbiome_metabolome_cor %>%
  dplyr::filter(cor_p_adjust < 0.01)

skin_microbiome_metabolome_cor =
  skin_microbiome_metabolome_cor %>%
  dplyr::filter(cor_p_adjust < 0.01)

oral_microbiome_metabolome_cor =
  oral_microbiome_metabolome_cor %>%
  dplyr::filter(cor_p_adjust < 0.01)

nasal_microbiome_metabolome_cor =
  nasal_microbiome_metabolome_cor %>%
  dplyr::filter(cor_p_adjust < 0.01)

temp_data = 
  rbind(
    data.frame(genus = unique(st_microbiome_metabolome_cor$microbiome), class = "Stool"),
    data.frame(genus = unique(skin_microbiome_metabolome_cor$microbiome), class = "Skin"),
    data.frame(genus = unique(nasal_microbiome_metabolome_cor$microbiome), class = "Nasal"),
    data.frame(genus = unique(oral_microbiome_metabolome_cor$microbiome), class = "Oral")
  ) %>% 
  tidyr::pivot_wider(names_from = class, values_from = class) %>% 
  dplyr::mutate(Stool = case_when(
    is.na(Stool) ~ FALSE,
    !is.na(Stool) ~ TRUE
  )) %>% 
  dplyr::mutate(Skin = case_when(
    is.na(Skin) ~ FALSE,
    !is.na(Skin) ~ TRUE
  )) %>% 
  dplyr::mutate(Nasal = case_when(
    is.na(Nasal) ~ FALSE,
    !is.na(Nasal) ~ TRUE
  )) %>% 
  dplyr::mutate(Oral = case_when(
    is.na(Oral) ~ FALSE,
    !is.na(Oral) ~ TRUE
  ))

plot = 
  upset(data = temp_data, intersect = c("Stool", "Skin", "Nasal", "Oral"))
plot

# ggsave(plot, filename = "final_network_body_site_microbe_upset.pdf", width = 9, height = 7)






temp_data = 
  rbind(
    data.frame(genus = unique(st_microbiome_metabolome_cor$metabolite), class = "Stool"),
    data.frame(genus = unique(skin_microbiome_metabolome_cor$metabolite), class = "Skin"),
    data.frame(genus = unique(nasal_microbiome_metabolome_cor$metabolite), class = "Nasal"),
    data.frame(genus = unique(oral_microbiome_metabolome_cor$metabolite), class = "Oral")
  ) %>% 
  tidyr::pivot_wider(names_from = class, values_from = class) %>% 
  dplyr::mutate(Stool = case_when(
    is.na(Stool) ~ FALSE,
    !is.na(Stool) ~ TRUE
  )) %>% 
  dplyr::mutate(Skin = case_when(
    is.na(Skin) ~ FALSE,
    !is.na(Skin) ~ TRUE
  )) %>% 
  dplyr::mutate(Nasal = case_when(
    is.na(Nasal) ~ FALSE,
    !is.na(Nasal) ~ TRUE
  )) %>% 
  dplyr::mutate(Oral = case_when(
    is.na(Oral) ~ FALSE,
    !is.na(Oral) ~ TRUE
  ))

plot = 
  upset(data = temp_data, intersect = c("Stool", "Skin", "Nasal", "Oral"))
plot

# ggsave(plot, filename = "final_network_body_site_metabolites_upset.pdf", width = 9, height = 7)


st_number2 = length(unique(st_microbiome_metabolome_cor$microbiome))
skin_number2 = length(unique(skin_microbiome_metabolome_cor$microbiome))
nasal_number2 = length(unique(nasal_microbiome_metabolome_cor$microbiome))
oral_number2 = length(unique(oral_microbiome_metabolome_cor$microbiome))

st_number2*100/st_number1
skin_number2*100/skin_number1
nasal_number2*100/nasal_number1
oral_number2*100/oral_number1

nrow(st_microbiome_metabolome_cor)/st_number2
nrow(skin_microbiome_metabolome_cor)/skin_number2
nrow(nasal_microbiome_metabolome_cor)/nasal_number2
nrow(oral_microbiome_metabolome_cor)/oral_number2

nrow(st_microbiome_metabolome_cor)
nrow(skin_microbiome_metabolome_cor)
nrow(nasal_microbiome_metabolome_cor)
nrow(oral_microbiome_metabolome_cor)


length(unique(st_microbiome_metabolome_cor$microbiome)) + length(unique(st_microbiome_metabolome_cor$metabolite))
length(unique(skin_microbiome_metabolome_cor$microbiome)) + length(unique(skin_microbiome_metabolome_cor$metabolite))
length(unique(nasal_microbiome_metabolome_cor$microbiome)) + length(unique(nasal_microbiome_metabolome_cor$metabolite))
length(unique(oral_microbiome_metabolome_cor$microbiome)) + length(unique(oral_microbiome_metabolome_cor$metabolite))


########integrative network

###plasma metabolome
load(here::here("data_analysis/metabolome/data_preparation/variable_info"))
metabolome_variable_info = variable_info

dim(st_microbiome_metabolome_cor)
dim(skin_microbiome_metabolome_cor)
dim(nasal_microbiome_metabolome_cor)
dim(oral_microbiome_metabolome_cor)

edge_data_stool =
  st_microbiome_metabolome_cor %>%
  dplyr::mutate(p = -log(cor_p_adjust2, 10)) %>%
  dplyr::select(microbiome, metabolite, p, cor) %>%
  dplyr::rename(from = microbiome,
                to = metabolite,
                cor = cor) %>% 
  dplyr::mutate(from = paste("stool", from, sep = "_")) %>% 
  dplyr::filter(abs(cor) > 0.25)

node_data_stool = 
  data.frame(node = unique(c(edge_data_stool$from, edge_data_stool$to))) %>% 
  dplyr::left_join(metabolome_variable_info[,c("variable_id", "Metabolite")],
                   by = c("node" = "variable_id")) %>% 
  dplyr::rename(true_name = Metabolite) %>% 
  dplyr::mutate(class = case_when(
    is.na(true_name) ~ "Stool microbiome",
    !is.na(true_name) ~ "Metabolite"
  )) 

edge_data_skin =
  skin_microbiome_metabolome_cor %>%
  dplyr::mutate(p = -log(cor_p_adjust2, 10)) %>%
  dplyr::select(microbiome, metabolite, p, cor) %>%
  dplyr::rename(from = microbiome,
                to = metabolite,
                cor = cor) %>% 
  dplyr::mutate(from = paste("skin", from, sep = "_")) %>% 
  dplyr::filter(abs(cor) > 0.25)

node_data_skin = 
  data.frame(node = unique(c(edge_data_skin$from, edge_data_skin$to))) %>% 
  dplyr::left_join(metabolome_variable_info[,c("variable_id", "Metabolite")],
                   by = c("node" = "variable_id")) %>% 
  dplyr::rename(true_name = Metabolite) %>% 
  dplyr::mutate(class = case_when(
    is.na(true_name) ~ "Skin microbiome",
    !is.na(true_name) ~ "Metabolite"
  )) 




edge_data_oral =
  oral_microbiome_metabolome_cor %>%
  dplyr::mutate(p = -log(cor_p_adjust2, 10)) %>%
  dplyr::select(microbiome, metabolite, p, cor) %>%
  dplyr::rename(from = microbiome,
                to = metabolite,
                cor = cor) %>% 
  dplyr::mutate(from = paste("oral", from, sep = "_")) %>% 
  dplyr::filter(abs(cor) > 0.25)

node_data_oral = 
  data.frame(node = unique(c(edge_data_oral$from, edge_data_oral$to))) %>% 
  dplyr::left_join(metabolome_variable_info[,c("variable_id", "Metabolite")],
                   by = c("node" = "variable_id")) %>% 
  dplyr::rename(true_name = Metabolite) %>% 
  dplyr::mutate(class = case_when(
    is.na(true_name) ~ "Oral microbiome",
    !is.na(true_name) ~ "Metabolite"
  )) 



edge_data_nasal =
  nasal_microbiome_metabolome_cor %>%
  dplyr::mutate(p = -log(cor_p_adjust2, 10)) %>%
  dplyr::select(microbiome, metabolite, p, cor) %>%
  dplyr::rename(from = microbiome,
                to = metabolite,
                cor = cor) %>% 
  dplyr::mutate(from = paste("nasal", from, sep = "_")) %>% 
  dplyr::filter(abs(cor) > 0.25)

node_data_nasal = 
  data.frame(node = unique(c(edge_data_nasal$from, edge_data_nasal$to))) %>% 
  dplyr::left_join(metabolome_variable_info[,c("variable_id", "Metabolite")],
                   by = c("node" = "variable_id")) %>% 
  dplyr::rename(true_name = Metabolite) %>% 
  dplyr::mutate(class = case_when(
    is.na(true_name) ~ "Nasal microbiome",
    !is.na(true_name) ~ "Metabolite"
  )) 



node_data = 
  rbind(
    node_data_stool,
    node_data_skin,
    node_data_nasal,
    node_data_oral
  ) %>% 
  dplyr::distinct(node, .keep_all = TRUE)

edge_data = 
  rbind(edge_data_stool,
        edge_data_nasal,
        edge_data_skin,
        edge_data_oral)

node_data$true_name[is.na(node_data$true_name)] = 
  node_data$node[is.na(node_data$true_name)] %>% 
  stringr::str_replace("stool_|skin_|nasal_|oral_", "")

node_data = 
  node_data %>% 
  dplyr::filter(!stringr::str_detect(true_name, "C[0-9]{1,2}H[0-9]{1,2}"))

node_data$true_name = stringr::str_replace_all(node_data$true_name, '\"', "")

edge_data = 
  edge_data %>% 
  dplyr::filter(from %in% node_data$node & to %in% node_data$node)

node_data = 
  node_data %>% 
  dplyr::filter(node %in% edge_data$from | node %in% edge_data$to) %>% 
  dplyr::arrange(class) %>% 
  dplyr::mutate(node = factor(node, levels = node))

edge_data = 
edge_data %>% 
  dplyr::mutate(edge_class = case_when(
    stringr::str_detect(from, "stool") ~ "Stool",
    stringr::str_detect(from, "skin") ~ "Skin",
    stringr::str_detect(from, "nasal") ~ "Nasal",
    stringr::str_detect(from, "oral") ~ "Oral"
  ))


temp_data <- 
  tidygraph::tbl_graph(nodes = node_data, 
                       edges = edge_data,
                       directed = FALSE) %>% 
  dplyr::mutate(Degree = centrality_degree(mode = 'all'))

angle <- 360 * (c(1:nrow(node_data)) - 0.5)/nrow(node_data)
hjust <- ifelse(angle > 180, 1, 0)
angle <- ifelse(angle > 180, 90 - angle + 180, 90 - angle)

# plot <-
#   ggraph(temp_data,
#          layout = 'fr',
#          circular = FALSE) +
#   geom_edge_link(strength = 1,
#                  aes(color = cor, width = p),
#                  alpha = 1,
#                  show.legend = TRUE) +
#   geom_node_point(aes(size = Degree,
#                       color = class), 
#                   shape = 16, 
#                   alpha = 1, show.legend = TRUE) +
#   geom_node_text(aes(x = x ,
#                      y = y ,
#                      color = "black",
#                      label = true_name),
#                  repel = TRUE,
#                  size = 2) +
#   # shadowtext::geom_shadowtext(aes(x = x, y = y,
#   #                                 label = ifelse(class == "Metabolite", NA, true_name),
#   #                                 color = class), bg.colour = "white") +
#   guides(edge_color = guide_edge_colorbar(title = "Correlation"),
#          color = guide_legend(title = "Class")) +
#   scale_edge_width_continuous(range = c(0.1,2)) +
#   scale_size_continuous(range = c(3, 10)) +
#   scale_color_manual(values = omics_color) +
#   scale_edge_color_gradient2(low = viridis::cividis(n = 2)[1], 
#                              mid = "white", 
#                              high = "red") +
#   ggraph::theme_graph() +
#   theme(plot.background = element_rect(fill = "transparent", color = NA), 
#         panel.background = element_rect(fill = "transparent", color = NA),
#         legend.position = "right",
#         legend.background = element_rect(fill = "transparent", color = NA))
# 
# plot



plot <-
  ggraph(temp_data,
         layout = 'linear',
         circular = TRUE) +
  geom_edge_arc(strength = 1,
                 aes(color = edge_class, 
                     width = p),
                 alpha = 0.7,
                 show.legend = TRUE) +
  geom_node_point(aes(size = Degree,
                      color = class), 
                  shape = 16, 
                  alpha = 0.7, 
                  show.legend = TRUE) +
  geom_node_text(
    aes(
      x = x ,
      y = y ,
      color = class,
      label = true_name
    ),
    angle = angle,
    hjust = hjust,
    repel = FALSE,
    size = 2
  ) +
  guides(edge_color = guide_edge_colorbar(title = "Correlation"),
         color = guide_legend(title = "Class")) +
  scale_edge_width_continuous(range = c(0.1,2)) +
  scale_size_continuous(range = c(3, 10)) +
  scale_color_manual(values = omics_color) +
  # scale_edge_color_gradient2(low = viridis::cividis(n = 2)[1], 
  #                            mid = "white", 
  #                            high = "red") +
  scale_edge_color_manual(values = body_site_color) +
  ggraph::theme_graph() +
  theme(plot.background = element_rect(fill = "transparent", color = NA), 
        panel.background = element_rect(fill = "transparent", color = NA),
        legend.position = "right",
        legend.background = element_rect(fill = "transparent", color = NA))

plot

extrafont::loadfonts()

ggsave(plot, filename = "microbiome_vs_metabolome.pdf", width = 14, height = 12)































