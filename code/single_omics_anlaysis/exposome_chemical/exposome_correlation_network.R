##avoid source
no_exist_function()

##load data
masstools::setwd_project()
library(tidyverse)
rm(list = ls())

load("data_analysis/exposomeChemical_metabolome/exp_expression_data")
load("data_analysis/exposomeChemical_metabolome/met_expression_data")
load("data_analysis/exposomeChemical_metabolome/exp_sample_info")
load("data_analysis/exposomeChemical_metabolome/met_sample_info")

# sample_info <-
#   sample_info %>% 
#   dplyr::filter(start_date >= "2016-01-15" & start_date <= "2016-03-03")
# 
# sample_id <-
#   sample_info$sample_id
# 
# sample_id <-
#   c(paste(sample_id, c(1), sep = "_"),
#     paste(sample_id, c(2), sep = "_"),
#     paste(sample_id, c(3), sep = "_"))
# 
# expression_data <-
#   expression_data %>% 
#   dplyr::select(one_of(sample_id))

dim(exp_expression_data)
dim(met_expression_data)

colnames(exp_expression_data)
colnames(met_expression_data)

####load metabolomics data
load("data_20200511/metabolome/variable_info")
met_variable_info <- variable_info
met_annotation_info <- readr::read_csv("data_20200511/metabolome/annotation_table.csv")
met_annotation_info %>% dim()

variable_info <-
  variable_info %>% 
  left_join(met_annotation_info, by = c("peak_name" = "name"))

load("data_20200511/exposome/variable_info")
exp_variable_info <- variable_info

setwd("data_analysis/exposomeChemical_metabolome/exp_correlation_network")

#######correlation analysis
exp_expression_data <- log(exp_expression_data + 1, 2)

cor_value <- 
cor(t(exp_expression_data), method = "spearman")

cor_value[lower.tri(cor_value)] <- 0

cor_value <- 
cor_value %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "from") %>% 
  tidyr::pivot_longer(-from, names_to = "to", values_to = "cor") %>% 
  dplyr::filter(abs(cor) > 0.9) %>% 
  dplyr::filter(cor != 1)

p_value <- 
  purrr::map(as.data.frame(t(cor_value)), .f = function(x){
    x <- as.character(x)
    cor.test(as.numeric(exp_expression_data[x[1],]),
             as.numeric(exp_expression_data[x[2],]), method = "spearman")$p.value
  }) %>% 
  unlist()

p_adjust <- p.adjust(p_value, method = "fdr")

cor_value <- data.frame(cor_value, p_value, p_adjust, stringsAsFactors = FALSE)

cor_value <- 
  cor_value %>% 
  dplyr::filter(abs(cor) > 0.9 & p_adjust < 0.05)

cor_value <- 
  cor_value %>% 
  dplyr::left_join(exp_variable_info[,c(1:2)], by = c("from" = "peak_ID")) %>% 
  dplyr::rename(from_id = MetabID) %>% 
  dplyr::left_join(exp_variable_info[,c(1,2)], by = c("to" = "peak_ID")) %>% 
  dplyr::rename(to_id = MetabID) 

library(igraph)
library(ggraph)
library(tidygraph)

edge_data <-  
  cor_value %>% 
  dplyr::select(-c(from, to)) %>% 
  dplyr::rename(from = from_id, 
                to = to_id, 
                Correlation = cor) %>% 
  dplyr::mutate(p.value = -log(p_adjust, 10))

node_data <- 
  cor_value %>% 
  dplyr::select(-c(from, to)) %>% 
  dplyr::rename(from = from_id, to = to_id) %>% 
  dplyr::select(from, to) %>% 
  tidyr::pivot_longer(cols = c(from, to), 
                      names_to = "class", values_to = "node") %>% 
  dplyr::select(node) %>% 
  dplyr::left_join(exp_variable_info, by = c("node" = "MetabID")) %>% 
  dplyr::distinct(node, .keep_all = TRUE)

node_data <- 
  node_data[node_data$Chemical_class != "Unknown",] %>% 
  dplyr::arrange(Chemical_class)

edge_data <- 
  edge_data %>% 
  dplyr::filter((from %in% node_data$node) & (to %in% node_data$node))

###output node data and edge data
library(openxlsx)
wb = createWorkbook()
modifyBaseFont(wb, fontSize = 12, fontName = "Arial Narrow")
addWorksheet(wb, sheetName = "Node information", gridLines = TRUE)
addWorksheet(wb, sheetName = "Edge information", gridLines = TRUE)
freezePane(wb, sheet = 1, firstRow = TRUE, firstCol = TRUE) 
freezePane(wb, sheet = 2, firstRow = TRUE, firstCol = TRUE) 
writeDataTable(wb, sheet = 1, x = node_data,
               colNames = TRUE, rowNames = FALSE)
writeDataTable(wb, sheet = 2, x = edge_data %>% dplyr::select(from, to, everything()),
               colNames = TRUE, rowNames = FALSE)

saveWorkbook(wb, "exposomeChemical_correlation.xlsx", overwrite = TRUE)




temp_data <- 
  tidygraph::tbl_graph(nodes = node_data, 
                       edges = edge_data,
                       directed = FALSE) %>% 
  dplyr::mutate(Degree = centrality_degree(mode = 'all'))


chemical_class <- 
  node_data$Chemical_class %>% unique() %>% 
  sort()

library(ggsci)

col <-
  ggsci::pal_igv()(15)

names(col) <- chemical_class

library(wesanderson)
names(wes_palettes)
wes_palette(name = "Zissou1", n = 100, type = "continuous")

pal <-
  wesanderson::wes_palette(name = "Zissou1", n = 100, type = "continuous")

plot <- 
  ggraph(temp_data,
         layout = 'kk',
         circular = FALSE) +
  geom_edge_link(
    aes(color = Correlation),
    width = 0.3,
    show.legend = TRUE
    # arrow = arrow(length = unit(2, 'mm'))
  ) +
  geom_node_point(aes(fill = Chemical_class, 
                      size = Degree), 
                  shape = 21,
                  show.legend = TRUE) +
  ggraph::geom_node_text(aes(label = node), repel = TRUE, size = 2) +
  ggraph::scale_edge_color_gradientn(colours = pal) +
  ggraph::scale_edge_width(range = c(0.05, 0.5)) +
  scale_size_continuous(range = c(0.3, 7)) +
  scale_fill_manual(values = col) +
  guides(fill = guide_legend(title = "Chemical class", override.aes = list(size = 4))) +
  theme_void() +
  theme(plot.background = element_rect(fill = "transparent", color = NA), 
        panel.background = element_rect(fill = "transparent", color = NA))
plot

ggsave(plot, filename = "exp_correlation_network.png", 
       width = 11, height = 7, bg = "transparent")

ggsave(plot, filename = "exp_correlation_network.pdf", 
       width = 11, height = 7, bg = "transparent")

save(cor_value, file = "cor_value")
write.csv(cor_value, "cor_value.csv", row.names = FALSE)


cor_value

core_met <- "Pymetrozin"

grep(core_met, cor_value$from_id, value = TRUE)
grep(core_met, cor_value$to_id, value = TRUE)

idx1 <- grep(core_met, cor_value$from_id, value = FALSE)
idx2 <- grep(core_met, cor_value$to_id, value = FALSE)

id <- 
unique(c(
  cor_value[idx1,]$from,
  cor_value[idx1,]$to,
  cor_value[idx2,]$from
)
)

temp_data <- 
  exp_expression_data[id,] %>% 
  apply(1, function(x){
    (x - mean(x))/sd(x)
  }) %>% 
  t() %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "peak_id") %>% 
  dplyr::left_join(node_data[,c("node", "peak_ID")], by = c("peak_id" = "peak_ID")) 


plot <- 
temp_data %>% 
  tidyr::pivot_longer(cols = -c(peak_id, node), 
                      names_to = "time", values_to = "value") %>% 
  ggplot(aes(time, value)) +
  geom_line(aes(color = node, group = node)) +
  geom_point(aes(fill = node, group = node), shape = 21, 
             size = 2) +
  ggsci::scale_color_futurama() +
  ggsci::scale_fill_futurama() +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13),
    panel.grid.minor = element_blank(),
    # legend.position = c(0, 1),
    # legend.justification = c(0, 1),
    legend.background = element_rect(fill = "transparent"),
    plot.background = element_rect(fill = "transparent"),
    panel.background = element_rect(fill = "transparent")
  ) +
  labs(
    x = "Time point",
    y = "Scaled intensity"
  )

plot

ggsave(plot, filename = paste(core_met, "_line_plot.pdf", sep = ""), width = 10, height = 7)
ggsave(plot, filename = paste(core_met, "_line_plot.png", sep = ""), width = 10, height = 7)


