##avoid source
no_function()

masstools::setwd_project()
library(tidyverse)
setwd("data_analysis/exposome_biological/")
rm(list = ls())

##load data
load("dna_expression_data")
load("dna_sample_info")
load("dna_variable_info")

load("rna_expression_data")
load("rna_sample_info")
load("rna_variable_info")

dim(dna_variable_info)

table(dna_variable_info$level)

plot <-
  dna_sample_info %>%
  dplyr::mutate(subject_id = "69-001") %>%
  ggplot(aes(as.Date(date.start), y = subject_id)) +
  geom_point(
    size = 5,
    shape = 21,
    color = "black",
    fill = ggsci::pal_d3()(10)[5]
  ) +
  theme_bw() +
  scale_x_continuous(
    trans = "date",
    breaks = c(as.Date(dna_sample_info$date.start)),
    labels = as.character(dna_sample_info$date.start)
  ) +
  labs(x = "", y = "") +
  theme(
    axis.title = element_text(size = 10),
    axis.text.x = element_text(
      size = 10,
      angle = 45,
      vjust = 1,
      hjust = 1
    ),
    axis.text.y = element_text(size = 10),
    plot.margin = unit(c(0, 0, 0, 0), "pt")
  )

plot



# ggsave(plot,
#        file = "microbiome_sample_collection.pdf",
#        width = 9,
#        height = 7)

temp_data = 
dna_variable_info$level %>%
  table() %>% 
  data.frame()

colnames(temp_data) = c("Level", "Freq")


plot = 
temp_data %>%
  dplyr::mutate(Level = factor(Level, levels = c("phylum", "family", "genus", "species"))) %>%
  ggplot(aes(x = 2, y = Freq)) +
  geom_bar(
    stat = "identity",
    position = "stack",
    aes(fill = Level),
    color = "black",
    show.legend = FALSE
  ) +
  ggsci::scale_fill_npg() +
  theme_void() +
  xlim(0.5, 2.5) +
  ggplot2::coord_polar(theta = "y")  

# ggsave(plot, filename = "dna_level.pdf", width = 7, height = 7)

##kingdom 界 phylum 门 class 纲
###DNA
level <- "genus"

temp_expression_data <-
  dna_expression_data[rownames(dna_expression_data) %in%
                        dna_variable_info$variable_id[dna_variable_info$level ==
                                                        level],]

colnames(temp_expression_data) <-
  as.character(dna_sample_info$date.start)

temp_expression_data <- 
  temp_expression_data %>% 
  apply(2, function(x){
    x/sum(x)
  }) %>% 
  as.data.frame()

idx <- 
apply(temp_expression_data, 2, function(x){
  names(x) <- rownames(temp_expression_data)
  x <- sort(x, decreasing = TRUE)
  names(x)[1:3]
}) %>% 
  as.vector() %>% 
  unique()

temp_expression_data1 <- 
  temp_expression_data[idx,]

temp_expression_data2 <- 
  temp_expression_data[-match(idx, rownames(temp_expression_data)),] %>% 
  apply(2, sum) %>% 
  matrix(nrow = 1) %>% 
  as.data.frame()

colnames(temp_expression_data2) <- colnames(temp_expression_data1)
rownames(temp_expression_data2) <- paste(level, "Other", sep = "_")

temp_expression_data <- rbind(temp_expression_data1,
                              temp_expression_data2)

rownames(temp_expression_data)

library(Polychrome)

level_color <-
  createPalette(N = nrow(temp_expression_data),  
                seedcolors = ggsci::pal_aaas()(n=8))
  # colorRampPalette(colors = ggsci::pal_aaas()(n=8))(n=nrow(temp_expression_data))
  # ggsci::pal_aaas()(n = nrow(temp_expression_data))

names(level_color) <-
  stringr::str_split(rownames(temp_expression_data), "_", n = 2) %>%
  purrr::map(
    .f = function(x) {
      x[2]
    }
  ) %>%
  unlist()

level_color["Other"] <- "grey"

temp_expression_data %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "variable_id") %>%
  dplyr::mutate(mean_int = rowMeans(dplyr::select(., -variable_id))) %>%
  tidyr::pivot_longer(
    cols = -c(variable_id, mean_int),
    names_to = "sample",
    values_to = "value"
  ) %>%
  dplyr::left_join(dna_variable_info[, c("variable_id", "short_name")], by = "variable_id") %>%
  dplyr::mutate(short_name = 
                  case_when(
                    is.na(short_name) ~ "Other",
                    TRUE ~ short_name
                  )) %>% 
  dplyr::arrange(mean_int) %>%
  dplyr::mutate(short_name = factor(short_name, levels = unique(short_name))) %>%
  ggplot(aes(x = sample, y = value)) +
  geom_bar(aes(fill = short_name),
           stat = "identity",
           show.legend = TRUE) +
  theme_bw() +
  scale_fill_manual(values = level_color) +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  scale_x_discrete(expand = expansion(mult = c(0, 0))) +
  labs(x = "", y = "Proportion") +
  theme(
    axis.text.x = element_text(
      size = 10,
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    panel.grid = element_blank()
  )

###sankey
library(ggalluvial)

colnames(temp_expression_data) <-
  stringr::str_replace(colnames(temp_expression_data), "2016-", "")

plot <-
  temp_expression_data %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "variable_id") %>%
  dplyr::mutate(mean_int = rowMeans(dplyr::select(., -variable_id))) %>%
  tidyr::pivot_longer(
    cols = -c(variable_id, mean_int),
    names_to = "sample",
    values_to = "value"
  ) %>%
  dplyr::left_join(dna_variable_info[, c("variable_id", "short_name")], by = "variable_id") %>%
  dplyr::mutate(short_name = 
                  case_when(
                    is.na(short_name) ~ "Other",
                    TRUE ~ short_name
                  )) %>% 
  dplyr::arrange(mean_int) %>%
  dplyr::mutate(short_name = factor(short_name, levels = unique(short_name))) %>% 
  ggplot(aes(
    y = value,
    x = sample,
    stratum = short_name,
    alluvium = short_name,
    fill = short_name
  )) +
  geom_alluvium(width = 0.3, alpha = 0.5) +
  geom_stratum(width = 0.5) +
  labs(x = "", y = "Proportion") +
  theme_bw() +
  scale_fill_manual(values = level_color) +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  scale_x_discrete(expand = expansion(mult = c(0, 0))) +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13),
    axis.text.x = element_text(
      size = 10,
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "bottom"
  )

plot

# ggsave(
#   plot,
#   filename = paste("top_3_", level, ".pdf", sep = ""),
#   width = 9,
#   height = 7
# )
# 
# ggsave(
#   plot,
#   filename = paste("top_3_", level, ".png", sep = ""),
#   width = 9,
#   height = 7
# )


####PCA analysis
##kingdom 界 phylum 门 class 纲
###DNA
level <- "genus"
temp_expression_data <-
  dna_expression_data[rownames(dna_expression_data) %in%
                        dna_variable_info$variable_id[dna_variable_info$level ==
                                                        level],]

colnames(temp_expression_data) <-
  as.character(dna_sample_info$date.start)

##remove rows with 50% zero
remain_idx <- 
apply(temp_expression_data, 1, function(x){
  sum(x == 0)/ncol(temp_expression_data)
}) %>% 
  `<`(0.8) %>% 
  which()

temp_expression_data <-
  temp_expression_data[remain_idx,]

temp_expression_data <-
  log(temp_expression_data + 1, 2)

temp_expression_data <- 
  temp_expression_data %>% 
  apply(1, function(x){
    (x - mean(x))/sd(x)
  }) %>% 
  t() %>% 
  as.data.frame()

rownames(temp_expression_data)

#PCA analysis
###correlate date
###PCA for date
pca_object <- 
  prcomp(x = t(temp_expression_data), center = FALSE, scale. = FALSE)

library(wesanderson)
names(wes_palettes)
wes_palette(name = "Zissou1", n = 100, type = "continuous")

pal <-
  wesanderson::wes_palette(name = "Zissou1", n = 100, type = "continuous")

as.Date_origin <- function(x) {
  as.Date(x, origin = '1970-01-01')
}

plot <-
  pca_object$x %>%
  as.data.frame() %>%
  rownames_to_column(var = "date") %>%
  left_join(dna_sample_info %>% dplyr::mutate(date = as.character(date.start)), 
            by = c("date")) %>%
  mutate(date2 = as.character(date)) %>%
  mutate(date = as.integer(as.Date(date))) %>%
  ggplot(aes(PC1, PC2, colour = date)) +
  geom_vline(xintercept = 0,
             linetype = 2,
             color = "black") +
  geom_hline(yintercept = 0,
             linetype = 2,
             color = "black") +
  geom_point(shape = 21,
             size = 4,
             aes(fill = date),
             color = "black") +
  guides(colour = guide_colourbar(title = "Date")) +
  scale_colour_gradientn(colours = pal,
                         labels = as.Date_origin) +
  scale_fill_gradientn(colours = pal,
                       labels = as.Date_origin) +
  ggrepel::geom_text_repel(aes(PC1, PC2, label = paste(date2, location, sep = "_")), 
                            show.legend = FALSE) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13),
    panel.grid.minor = element_blank(),
    legend.position = c(1, 1),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = "transparent"),
    plot.background = element_rect(fill = "transparent"),
    panel.background = element_rect(fill = "transparent")
  ) +
  labs(
    x = paste("PC1 (", round(summary(pca_object)$importance[2, 1], 4) * 100, "%)", sep = ""),
    y = paste("PC2 (", round(summary(pca_object)$importance[2, 2], 4) *
                100, "%)", sep = "")
  )

plot

# ggsave(plot,
#        file = "pca_date.pdf",
#        width = 7,
#        height = 7)
# 
# ggsave(plot,
#        file = "pca_data.png",
#        width = 7,
#        height = 7)








