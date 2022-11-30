###
no_function()

masstools::setwd_project()
library(tidyverse)
rm(list = ls())
library(tidyverse)
source("code/tools.R")

####load raw data
load("data/LipidData2020/Step5_wiff_assembled_9_2020-Apr-20.RData")

load(here::here("data_analysis/metabolome/data_preparation/sample_info"))
metabolome_sample_info = sample_info

variable_info =
  collapsed_df_filtered_impute_final %>%
  dplyr::select(
    -c(
      CollapsedSampleName,
      log_sample_nmol_per_g_concentration,
      IsImputed,
      ClassLogSumConcentration,
      ClassLogMedian,
      CollapsedTAGs,
      TAGsum,
      TAGmedian
    )
  ) %>%
  dplyr::distinct(LipidSpecies, .keep_all = TRUE) %>%
  dplyr::rename(variable_id = LipidSpecies)

sample_info =
  collapsed_df_filtered_impute_final %>%
  dplyr::select(CollapsedSampleName) %>%
  dplyr::rename(sample_id = CollapsedSampleName) %>%
  dplyr::distinct(sample_id, .keep_all = TRUE)

expression_data = collapsed_df_filtered_impute_final %>%
  dplyr::select(CollapsedSampleName,
                LipidSpecies,
                log_sample_nmol_per_g_concentration) %>%
  tidyr::pivot_wider(names_from = CollapsedSampleName, values_from = log_sample_nmol_per_g_concentration) %>%
  tibble::column_to_rownames(var = "LipidSpecies")

colnames(expression_data) == sample_info$sample_id
rownames(expression_data) == variable_info$variable_id

lipid_info <-
  read.table(
    here::here("data/LipidData2020/Lipomat05.txt"),
    header = TRUE,
    sep = "\t"
  )

variable_info =
  variable_info %>%
  dplyr::left_join(
    lipid_info %>% dplyr::distinct(Lipid_ID, .keep_all = TRUE),
    by = c("variable_id" = "Lipid_ID")
  )

batch =
  stringr::str_extract(sample_info$sample_id,
                       "batch[0-9]{1,3}")

sample_id =
  stringr::str_extract(sample_info$sample_id,
                       "[0-9]{2}\\-[0-9]{3,4}\\-[0-9]{2,6}|QC_[0-9]{1,3}")

date =
  stringr::str_extract(sample_info$sample_id,
                       "[0-9]{8}")

date =
  as.Date(date, "%m%d%Y")

sample_info$sample_id[which(is.na(sample_id))]

sample_id[!is.na(sample_id)]

sample_info =
  data.frame(sample_info, true_sample_id = sample_id, date, batch)

sample_info =
  sample_info %>%
  dplyr::filter(!is.na(true_sample_id))

sample_info$true_sample_id[grep("QC", sample_info$true_sample_id)]

sample_info$true_sample_id[grep("QC", sample_info$true_sample_id)] =
  paste(sample_info$batch[grep("QC", sample_info$true_sample_id)],
        sample_info$true_sample_id[grep("QC", sample_info$true_sample_id)],
        sep = "_")

sample_info =
  sample_info %>%
  dplyr::distinct(true_sample_id,
                  .keep_all = TRUE)

expression_data =
  expression_data[, sample_info$sample_id]

sample_info =
  sample_info %>%
  dplyr::rename(sample_id = true_sample_id,
                original_sample_id = sample_id) %>%
  dplyr::select(sample_id, dplyr::everything())

colnames(expression_data) = sample_info$sample_id

sample_info$subject_id =
  stringr::str_extract(string = sample_info$sample_id, pattern = "[0-9]{2}\\-[0-9]{3,4}")

sample_info$subject_id[is.na(sample_info$subject_id)] = "QC"

temp_sample_id =
  sample_info$sample_id[grep("69-028", sample_info$sample_id)]

metabolome_sample_info[match(temp_sample_id, metabolome_sample_info$sample_id), c("sample_id", "subject_id")]

sample_info$subject_id[grep("69-028", sample_info$sample_id)] = "70-1014"

dim(sample_info)
dim(variable_info)
dim(expression_data)

rownames(expression_data) == variable_info$variable_id
colnames(expression_data) == sample_info$sample_id

setwd("data_analysis/lipidome/data_preparation/")

######get the intra-lipidomics correlation network and see if some lipids have high correlations
which(variable_info$variable_id == "TAG58.6.FA20.4")
which(variable_info$variable_id == "TAG58.8.FA18.2")

plot(as.numeric(expression_data[830,]),
     as.numeric(expression_data[832,]))

#####k-means clustering
library(Mfuzz)
temp_data <-
  expression_data

rownames(temp_data)

time <- c(1:ncol(temp_data))

temp_data <- rbind(time, temp_data)

temp_data2 =
  temp_data %>%
  apply(1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()

row.names(temp_data)[1] <- "time"

# write.table(
#   temp_data,
#   file = "temp_data.txt",
#   sep = '\t',
#   quote = FALSE,
#   col.names = NA
# )

#read it back in as an expression set
library(Mfuzz)
data <- table2eset(filename = "temp_data.txt")

data.s <- standardise(data)
m1 <- mestimate(data.s)
m1

# plot <-
# Dmin(
#   data.s,
#   m = m1,
#   crange = seq(2, 40, 1),
#   repeats = 3,
#   visu = TRUE
# )
#
# plot <-
# plot %>%
#   data.frame(distance = plot,
#              k = seq(2,40,1)) %>%
#   ggplot(aes(k, distance)) +
#   geom_point(shape = 21, size = 4, fill = "black") +
#   geom_smooth() +
#   geom_segment(aes(x = k, y = 0, xend = k, yend = distance)) +
#   theme_bw() +
#   theme(
#     # legend.position = c(0, 1),
#     # legend.justification = c(0, 1),
#     panel.grid = element_blank(),
#     axis.title = element_text(size = 13),
#     axis.text = element_text(size = 12),
#     panel.background = element_rect(fill = "transparent", color = NA),
#     plot.background = element_rect(fill = "transparent", color = NA),
#     legend.background = element_rect(fill = "transparent", color = NA)
#   ) +
#   labs(
#     x = "Cluster number",
#     y = "Min. centroid distance"
#   ) +
#   scale_y_continuous(expand = expansion(mult = c(0,0.1)))
#
# plot
#
# ggsave(plot, filename = "distance_k_number.pdf", width = 7, height = 7)

clust = 6

# c <- mfuzz(data.s, c = clust, m = m1)
#
# mfuzz.plot(eset = data.s,
#            min.mem = 0.8,
#            cl = c,
#            mfrow=c(4,3),
#            time.labels = time,
#            new.window = FALSE)
#
# names(c$cluster) <- rownames(temp_data2)[-1]
# rownames(c$membership) <- rownames(temp_data2)[-1]
# save(c, file = "c")
load("c")

center <- c$centers
rownames(center) <- paste("Cluster", rownames(center), sep = ' ')
corrplot::corrplot(
  corr = cor(t(center)),
  type = "full",
  diag = TRUE,
  order = "hclust",
  hclust.method = "ward.D",
  # addrect = 5,
  col = colorRampPalette(colors = rev(
    RColorBrewer::brewer.pal(n = 11, name = "Spectral")
  ))(n = 100),
  number.cex = .7,
  addCoef.col = "black"
)

acore <- acore(data.s, c, min.acore = 0)
acore

centers <- c$centers
names(c$cluster) == rownames(c$membership)

cluster_info <-
  data.frame(
    variable_id = names(c$cluster),
    c$membership,
    cluster = c$cluster,
    stringsAsFactors = FALSE
  ) %>%
  arrange(cluster)

# openxlsx::write.xlsx(x = cluster_info,
#                      file = "cluster_info.xlsx", asTable = TRUE)

#####output the expression data of different clusters

###plot for each cluster
# for (cluster_idx in 1:clust) {
#   cat(cluster_idx, " ")
#   dir.create(paste("cluster", cluster_idx, sep = "_"))
#   cluster_data <-
#     cluster_info %>%
#     dplyr::filter(cluster_idx == cluster_idx) %>%
#     dplyr::select(1, 1 + cluster_idx)
#
#   colnames(cluster_data) <- c("variable_id", "membership")
#
#   cluster_data <-
#     cluster_data %>%
#     dplyr::filter(membership > 0.9)
#
#   openxlsx::write.xlsx(
#     x = cluster_data,
#     file = file.path(
#       paste("cluster", cluster_idx, sep = "_"),
#       paste("cluster", cluster_idx, ".xlsx", sep = "")
#     ),
#     asTable = TRUE,
#     overwrite = TRUE
#   )
#
#   ###cluster plot
#
#   temp =
#     temp_data2[cluster_data$variable_id,] %>%
#     data.frame(
#       membership = cluster_data$membership,
#       .,
#       stringsAsFactors = FALSE,
#       check.names = FALSE
#     ) %>%
#     tibble::rownames_to_column(var = "variable_id") %>%
#     tidyr::pivot_longer(
#       cols = -c(variable_id, membership),
#       names_to = "sample_id",
#       values_to = "value"
#     )
#   # dplyr::left_join(sample_info[, c("sample_id", "accurate_time")], by = "sample_id")
#
#   plot <-
#     ggplot() +
#     geom_hline(yintercept = 0) +
#     geom_line(aes(sample_id, value, group = variable_id, color = membership),
#               data = temp) +
#     scale_color_gradientn(colours = c(RColorBrewer::brewer.pal(n = 9, name = "OrRd")[c(1:7)])) +
#     theme_bw() +
#     theme(
#       panel.grid = element_blank(),
#       legend.position = c(0, 1),
#       legend.justification = c(0, 1),
#       panel.grid.minor = element_blank(),
#       axis.title = element_text(size = 13),
#       axis.text = element_text(size = 12),
#       axis.text.x = element_blank(),
#       axis.ticks.x = element_blank(),
#       panel.background = element_rect(fill = "transparent", color = NA),
#       plot.background = element_rect(fill = "transparent", color = NA),
#       legend.background = element_rect(fill = "transparent", color = NA)
#     ) +
#     labs(
#       x = "",
#       y = "Z-score",
#       title = paste(
#         "Cluster ",
#         cluster_idx,
#         " (",
#         nrow(cluster_data),
#         " lipids)",
#         sep = ""
#       )
#     )
#
#   plot
#
#   ggsave(
#     plot,
#     filename = file.path(
#       paste("cluster", cluster_idx, sep = "_"),
#       paste("cluster", cluster_idx, ".pdf", sep = "")
#     ),
#     width = 20,
#     height = 7
#   )
# }

## (4) feature number
cluster1 <- readxl::read_xlsx("cluster_1/cluster1.xlsx")
cluster2 <- readxl::read_xlsx("cluster_2/cluster2.xlsx")
cluster3 <- readxl::read_xlsx("cluster_3/cluster3.xlsx")
cluster4 <- readxl::read_xlsx("cluster_4/cluster4.xlsx")
cluster5 <- readxl::read_xlsx("cluster_5/cluster5.xlsx")
cluster6 <- readxl::read_xlsx("cluster_6/cluster6.xlsx")

###annotation for each cluster
cluster1 = data.frame(cluster1, cluster = "1")
cluster2 = data.frame(cluster2, cluster = "2")
cluster3 = data.frame(cluster3, cluster = "3")
cluster4 = data.frame(cluster4, cluster = "4")
cluster5 = data.frame(cluster5, cluster = "5")
cluster6 = data.frame(cluster6, cluster = "6")

cluster =
  rbind(cluster1,
        cluster2,
        cluster3,
        cluster4,
        cluster5,
        cluster6)

variable_info =
  variable_info %>%
  dplyr::left_join(cluster, by = "variable_id")

variable_info$cluster[is.na(variable_info$cluster)] = "Other"

variable_info$mol_name[!is.na(variable_info$Lipid_Name)] =
  variable_info$Lipid_Name[!is.na(variable_info$Lipid_Name)]

####get global correlation
dim(variable_info)
dim(expression_data)
#
# cor_data =
#   purrr::map(1:(nrow(expression_data) - 1), function(idx1) {
#     cat(idx1, " ")
#     purrr::map((idx1 + 1):nrow(expression_data), function(idx2) {
#       x = as.numeric(expression_data[idx1, ])
#       y = as.numeric(expression_data[idx2, ])
#       temp =
#         cor.test(x, y, method = "spearman")
#       data.frame(
#         name1 = variable_info$variable_id[idx1],
#         name2 = variable_info$variable_id[idx2],
#         cor = unname(temp$estimate),
#         p = unname(temp$p.value)
#       )
#     }) %>%
#       do.call(rbind, .) %>%
#       as.data.frame()
#   }) %>%
#   do.call(rbind, .) %>%
#   as.data.frame()
# save(cor_data, file = "cor_data")

load("cor_data")

cor_data$p_adjust = p.adjust(cor_data$p, method = "BH")

plot <- 
cor_data %>% 
  ggplot(aes(cor)) +
  geom_histogram(color = "black", binwidth = 0.05) +
  theme_bw() +
  labs(x = "Spearman correlation",
       y = "Count")

ggsave(plot, filename = "correlation distributation.pdf", width = 9, height = 7)

cor_data =
  cor_data %>%
  dplyr::filter(p_adjust < 0.05)

cor_data$p_adjust[cor_data$p_adjust == 0] = min(cor_data$p_adjust[cor_data$p_adjust != 0])

dim(cor_data)

cor_data %>%
  dplyr::arrange(desc(cor)) %>%
  head()

plot <-
  data.frame(x = as.numeric(expression_data["TAG52.3.FA16.0",]),
             y = as.numeric(expression_data["TAG52.3.FA18.2",])) %>%
  ggplot(aes(x, y)) +
  geom_point() +
  theme_bw() +
  labs(x = "TAG52.3.FA16.0", y = "TAG52.3.FA18.2") +
  geom_abline(slope = 1,
              intercept = 0,
              color = "red") +
  theme(panel.grid.minor = element_blank()) +
  annotate(geom = "text",
           x = 2,
           y = 3.0,
           label = "Spearman cor: 0.9994382") +
  annotate(geom = "text",
           x = 2,
           y = 2.8,
           label = "adjuste p value: 1.43279e-322")

plot

# ggsave(plot, filename = "example_cor.pdf", width = 7, height = 7)

plot(as.numeric(expression_data["TAG52.3.FA16.0",]),
     as.numeric(expression_data["TAG52.3.FA18.2",]))

cor.test(as.numeric(expression_data["TAG52.3.FA16.0",]),
         as.numeric(expression_data["TAG52.3.FA18.2",]),
         method = "spearman")

##### network for each cluster
# for (cluster_idx in 1:clust) {
#   cat(cluster_idx, " ")
#   edge_data =
#     cor_data %>%
#     dplyr::filter(abs(cor) > 0.8) %>%
#     dplyr::rename(from = name1, to = name2)
#
#   node_data =
#     variable_info %>%
#     dplyr::rename(node = variable_id) %>%
#     dplyr::filter(node %in% c(edge_data$from, edge_data$to)) %>%
#     dplyr::filter(cluster == cluster_idx)
#
#   edge_data =
#     edge_data %>%
#     dplyr::filter(from %in% node_data$node & to %in% node_data$node)
#
#   node_data =
#     node_data %>%
#     dplyr::filter(node %in% c(edge_data$from, edge_data$to))
#
#   if (nrow(node_data) > 0) {
#     library(ggraph)
#     library(igraph)
#     library(tidygraph)
#
#     graph <-
#       tidygraph::tbl_graph(nodes = node_data,
#                            edges = edge_data,
#                            directed = FALSE) %>%
#       dplyr::mutate(Degree = centrality_degree(mode = 'all'))
#
#     # subnetworks <-
#     #   igraph::cluster_edge_betweenness(graph = graph,
#     #                                    weights = abs(edge_attr(graph,
#     #                                                            "cor")))
#
#     subnetworks <-
#       igraph::cluster_fast_greedy(graph = graph,
#                                   weights = abs(edge_attr(graph,
#                                                           "cor")))
#
#     plot =
#       modularity_plot(subnetworks = subnetworks)
#
#     ggsave(
#       plot,
#       filename = file.path(
#         paste("cluster", cluster_idx, sep = "_"),
#         "all_modularity.pdf"
#       ),
#       width = 9,
#       height = 7
#     )
#
#     # table(igraph::membership(communities = subnetworks))
#
#     # membership =
#     #   igraph::cut_at(communities = subnetworks, no = 12)
#
#     # if (max(subnetworks$modularity) < 0.4) {
#     #   module = paste(cluster_idx, rep(1, length(
#     #     igraph::membership(communities = subnetworks)
#     #   )), sep = "_")
#     # } else{
#       module = paste(cluster_idx, igraph::membership(communities = subnetworks), sep = "_")
#     # }
#
#     remove_module_name =
#     data.frame(module) %>%
#       dplyr::group_by(module) %>%
#       dplyr::summarise(n = n()) %>%
#       dplyr::filter(n == 1) %>%
#       dplyr::pull(module)
#
#     module[module %in% remove_module_name] = "Other"
#
#     # module[module != "1_2"] = NA
#
#     graph =
#       graph %>%
#       tidygraph::activate(what = "nodes") %>%
#       dplyr::mutate(module = module)
#
#     node =
#       vertex_attr(graph) %>%
#       dplyr::bind_cols() %>%
#       as.data.frame()
#
#     save(node, file = file.path(paste("cluster", cluster_idx, sep = "_"), "node"))
#
#     plot <-
#       ggraph(graph, layout = "fr",
#              circular = FALSE) +
#       geom_edge_link(aes(color = cor,
#                          width = -log(p_adjust, 10)),
#                      alpha = 1,
#                      show.legend = TRUE) +
#       geom_node_point(
#         aes(size = Degree,
#             fill = module),
#         shape = 21,
#         alpha = 0.7,
#         # fill = class_color["lipidomics"],
#         show.legend = FALSE
#       ) +
#       ggforce::geom_mark_hull(
#         aes(
#           x = x,
#           y = y,
#           group = module,
#           fill = module,
#           filter = module != 'Other'
#         ),
#         concavity = 30,
#         show.legend = FALSE,
#         alpha = 0.15
#       ) +
#       scale_size_continuous(range = c(2, 8)) +
#       guides(
#         linetype = "none",
#         color = guide_colorbar(
#           title = "Lagged correlation",
#           override.aes = list(linetype = "none")
#         ),
#         size = guide_legend(
#           title = "Degree",
#           override.aes = list(
#             linetype = NA,
#             fill = "transparent",
#             shape = 21,
#             color = "black"
#           )
#         ),
#         fill = guide_legend(
#           title = "Class",
#           override.aes = list(
#             shape = 21,
#             size = 3,
#             alpha = 1
#           )
#         )
#       ) +
#       scale_edge_width_continuous(range = c(0.3, 2)) +
#       ggraph::scale_edge_color_gradient2(
#         low = alpha("#3B4992FF", 0.7),
#         mid = "white",
#         high = alpha("#EE0000FF", 0.7)
#       ) +
#       ggraph::theme_graph() +
#       theme(
#         plot.background = element_rect(fill = "transparent", color = NA),
#         panel.background = element_rect(fill = "transparent", color = NA),
#         legend.position = "right",
#         legend.background = element_rect(fill = "transparent", color = NA)
#       )
#
#     # if (nrow(node_data) < 150) {
#       plot =
#         plot +
#         geom_node_text(aes(x = x,
#                            y = y,
#                            label = mol_name),
#                        size = 2,
#                        check_overlap = TRUE)
#     # }
#
#     # extrafont::loadfonts()
#
#     # plot
#     ggsave(
#       plot,
#       filename = file.path(paste("cluster", cluster_idx, sep = "_"), "network.pdf"),
#       width = 9,
#       height = 7
#     )
#   } else{
#     node = node_data
#     save(node, file = file.path(paste("cluster", cluster_idx, sep = "_"), "node"))
#   }
# }


#####load the new node information
load("cluster_1/node")
node1 = node
load("cluster_2/node")
node2 = node
load("cluster_3/node")
node3 = node
load("cluster_4/node")
node4 = node
load("cluster_5/node")
node5 = node
load("cluster_6/node")
node6 = node

node_info =
  rbind(node1,
        node2,
        node3,
        node4,
        node5,
        node6)

#####plot for each module
# for (module_idx in 1:length(unique(node_info$module))) {
#   cat(module_idx, " ")
#   module = unique(node_info$module)[module_idx]
#   cluster_idx = stringr::str_split(module, pattern = "_")[[1]][1]
#   ###cluster plot
#
#   temp =
#     temp_data2[node_info$node[node_info$module == module],] %>%
#     tibble::rownames_to_column(var = "variable_id") %>%
#     tidyr::pivot_longer(
#       cols = -c(variable_id),
#       names_to = "sample_id",
#       values_to = "value"
#     )
#
#   plot <-
#     ggplot() +
#     geom_hline(yintercept = 0) +
#     geom_line(aes(sample_id, value, group = variable_id),
#               data = temp) +
#     scale_color_gradientn(colours = c(RColorBrewer::brewer.pal(n = 9, name = "OrRd")[c(1:7)])) +
#     theme_bw() +
#     theme(
#       panel.grid = element_blank(),
#       legend.position = c(0, 1),
#       legend.justification = c(0, 1),
#       panel.grid.minor = element_blank(),
#       axis.title = element_text(size = 13),
#       axis.text = element_text(size = 12),
#       axis.text.x = element_blank(),
#       axis.ticks.x = element_blank(),
#       panel.background = element_rect(fill = "transparent", color = NA),
#       plot.background = element_rect(fill = "transparent", color = NA),
#       legend.background = element_rect(fill = "transparent", color = NA)
#     ) +
#     labs(
#       x = "",
#       y = "Z-score",
#       title = paste("Module ",
#                     module,
#                     " (",
#                     length(unique(temp$variable_id)),
#                     " lipids)",
#                     sep = "")
#     )
#
#   plot
#
#   ggsave(
#     plot,
#     filename = file.path(
#       paste("cluster", cluster_idx, sep = "_"),
#       paste("module", module, ".pdf", sep = "")
#     ),
#     width = 20,
#     height = 7
#   )
# }


#####output the combined lipidomics data
dim(node_info)
dim(variable_info)

final_info =
  variable_info %>%
  dplyr::left_join(node_info[, c("node", "module")],
                   by = c("variable_id" = "node"))

final_info$module[is.na(final_info$module)] = "Other"

final_info %>%
  dplyr::filter(module == "1_1") %>%
  dplyr::pull(LipidClass)

final_info %>%
  dplyr::filter(module == "1_2") %>%
  dplyr::pull(LipidClass)

library(plyr)

final_info %>%
  plyr::dlply(.(module)) %>%
  purrr::map(function(x) {
    table(x$LipidClass)
  })

###output
library(openxlsx)

final_info =
  final_info %>%
  dplyr::arrange(module, LipidClass)

# openxlsx::write.xlsx(final_info,
#                      file = "final_info.xlsx",
#                      asTable = TRUE,
#                      overwrite = TRUE)

table(final_info$module)

#####annotation for each module
library(plyr)

final_info %>%
  plyr::dlply(.variables = .(module)) %>%
  purrr::map(function(x) {
    x[, c("variable_id", "LipidClass")]
  })

dim(expression_data)
dim(final_info)
dim(sample_info)

rownames(expression_data) == final_info$variable_id


variable_info$Lipid_Name[is.na(variable_info$Lipid_Name)] =
  variable_info$variable_id[is.na(variable_info$Lipid_Name)]

final_info$Lipid_Name[is.na(final_info$Lipid_Name)] =
  final_info$variable_id[is.na(final_info$Lipid_Name)]

new_variable_info =
  final_info %>%
  plyr::dlply(.variables = .(module)) %>%
  purrr::map(function(x) {
    cat(x$variable_id[1], " ")
    if (x$module[1] == "Other") {
      x$annotation = x$Lipid_Name
      return(x)
    }
    if (nrow(x) == 1) {
      x$annotation = x$Lipid_Name
    } else{
      x$variable_id = paste(x$variable_id, collapse = ";")
      x$annotation = paste(c("Lipid module", unique(x$module), unique(x$subclass)), collapse = "_")
      x$LipidClass = paste(x$LipidClass, collapse = ";")
      x$Lipid_Name = paste(x$Lipid_Name, collapse = ";")
      x$Lipid_Name_Tot = paste(x$Lipid_Name_Tot, collapse = ";")
      x$Lipid_Name_LION = paste(x$Lipid_Name_LION, collapse = ";")
      x$Lipid_Class_Detailed = paste(x$Lipid_Class_Detailed, collapse = ";")
      x$Lipid_Class = paste(x$Lipid_Class, collapse = ";")
      x$FA_1 = paste(x$FA_1, collapse = ";")
      x$FA_2 = paste(x$FA_2, collapse = ";")
      x$FA_All = paste(x$FA_All, collapse = ";")
      x$Sat_status_FA1 = paste(x$Sat_status_FA1, collapse = ";")
      x$Sat_status_FA2 = paste(x$Sat_status_FA2, collapse = ";")
      x$Saturation_All = paste(x$Saturation_All, collapse = ";")
      x$OddEven_FA1 = paste(x$OddEven_FA1, collapse = ";")
      x$OddEven_FA2 = paste(x$OddEven_FA2, collapse = ";")
      x$OddEven_All = paste(x$OddEven_All, collapse = ";")
      x$Omega_FA1 = paste(x$Omega_FA1, collapse = ";")
      x$Omega_FA2 = paste(x$Omega_FA2, collapse = ";")
      x$Omega_All = paste(x$Omega_All, collapse = ";")
      x$HMDB_ID = paste(x$HMDB_ID, collapse = ";")
      x$KEGG_ID = paste(x$KEGG_ID, collapse = ";")
      x$LION_ID = paste(x$LION_ID, collapse = ";")
      x$cellular_component = paste(x$cellular_component, collapse = ";")
      x$Function = paste(x$Function, collapse = ";")
      x$Total_Carb = paste(x$Total_Carb, collapse = ";")
      x$Total_Unsat = paste(x$Total_Unsat, collapse = ";")
      x$Carbon_numb_FA1 = paste(x$Carbon_numb_FA1, collapse = ";")
      x$Carbon_numb_FA2 = paste(x$Carbon_numb_FA2, collapse = ";")
      x$CarbonFA_All = paste(x$CarbonFA_All, collapse = ";")
      x$Unsat_numb_FA1 = paste(x$Unsat_numb_FA1, collapse = ";")
      x$Unsat_numb_FA2 = paste(x$Unsat_numb_FA2, collapse = ";")
      x$UnsatFA_All = paste(x$UnsatFA_All, collapse = ";")
      x$PhysicoChemical = paste(x$PhysicoChemical, collapse = ";")
      x$PathwaysKegg = paste(x$PathwaysKegg, collapse = ";")
      x$membership = paste(x$membership, collapse = ";")
      x$cluster = paste(x$cluster, collapse = ";")
      x$mol_name = paste(x$mol_name, collapse = ";")
      x$module = paste(x$module, collapse = ";")
      x =
        x %>%
        dplyr::distinct(variable_id, .keep_all = TRUE)
    }
    x
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()

new_expression_data =
  new_variable_info$variable_id %>%
  purrr::map(function(x) {
    x = stringr::str_split(x, pattern = ";")[[1]]
    expression_data[x,] %>%
      colSums()
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()

rownames(new_expression_data) = new_variable_info$annotation

new_variable_info$old_variable_id = new_variable_info$variable_id
new_variable_info$variable_id = new_variable_info$annotation

new_sample_info = sample_info

colnames(new_expression_data) == new_sample_info$sample_id
rownames(new_expression_data) == new_variable_info$variable_id

new_variable_info$variable_id

# save(new_expression_data, file = "new_expression_data")
# save(new_sample_info, file = "new_sample_info")
# save(new_variable_info, file = "new_variable_info")
#
# save(sample_info,
#      file = here::here("data_analysis/lipidome/data_preparation/sample_info"))
# save(
#   expression_data,
#   file = here::here("data_analysis/lipidome/data_preparation/expression_data")
# )
# save(
#   variable_info,
#   file = here::here("data_analysis/lipidome/data_preparation/variable_info")
# )


openxlsx::write.xlsx(
  new_variable_info,
  file = "new_variable_info.xlsx",
  asTable = TRUE,
  overwrite = TRUE
)
