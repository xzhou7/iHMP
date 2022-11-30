##avoid source
no_exist_function()
masstools::setwd_project()
rm(list = ls())

source("R/20200511/tools.R")

library(tidyverse)

##load exposome data
load("data_20200511/exposome/expression_data")
load("data_20200511/exposome/sample_info")
load("data_20200511/exposome/variable_info")

sample_info <-
  sample_info %>% 
  dplyr::filter(start_date >= as.Date("2016-01-15") & start_date <= as.Date("2016-03-03"))

sample_id <-
  sample_info$sample_id

sample_id <-
  c(paste(sample_id, c(1), sep = "_"),
    paste(sample_id, c(2), sep = "_"),
    paste(sample_id, c(3), sep = "_"))

expression_data <-
  expression_data %>% 
  dplyr::select(one_of(sample_id))

setwd("data_analysis/exposome/data_overview")

##Distance-based redundancy analysis (db-RDA) analysis
library(vegan)

##pvca analysis
# library(pvca)
# data <-
#   expression_data %>%
#   select(-contains("Blank"))
# 
# data <- log(data + 1, 2)
# 
# data <-
#   data %>%
#   t() %>% 
#   as.data.frame() %>%
#   rownames_to_column(var = "Sample_ID") %>%
#   mutate(Sample_ID = stringr::str_replace(Sample_ID, "_[0-9]{1,2}", "")) %>%
#   plyr::dlply(.variables = "Sample_ID") %>%
#   lapply(function(x) {
#     apply(x[, -1], 2, mean)
#   }) %>%
#   do.call(rbind, .)
# 
# data <-
#   apply(data, 2, function(x) {
#     (x - mean(x)) / sd(x)
#   }) 
# 
# 
# assay_data <- 
#   as.matrix(t(data))
# 
# pheno_data <-
#   sample_info %>%
#   dplyr::select(
#     sample_id,
#     start_date,
#     location
#   ) %>%
#   as.data.frame()
# 
# pheno_data$start_date <- 
#   as.numeric(as.Date(pheno_data$start_date) - as.Date("2016-01-15"))
# 
# pheno_data$location <- factor(pheno_data$location)
# 
# library(Biobase)
# 
# pheno_data <- AnnotatedDataFrame(data = pheno_data)
# 
# row.names(pheno_data) <- colnames(assay_data)
# 
# data_pvca <- 
#   Biobase::ExpressionSet(
#     assayData = assay_data,
#     phenoData = pheno_data
#   )
# 
# pct_threshold <- 0.6
# 
# batch.factors <-
#   c(
#     "start_date",
#     "location"
#   )
# 
# pvcaObj <- pvcaBatchAssess(abatch = data_pvca,
#                            batch.factors = batch.factors,
#                            threshold = pct_threshold)
# 
# save(pvcaObj, file = "pvcaObj")
# load("pvcaObj")
# 
# plot <- 
#   plot_pvca(object = pvcaObj)
# 
# plot
# ggsave(plot, filename = "pvca_plot.pdf", height = 9, width = 7)





#PCA analysis
data <-
  expression_data %>%
  dplyr::select(-contains("Blank"))

data <- log(data + 1, 2)

data <-
  data %>%
  t() %>% 
  as.data.frame() %>%
  rownames_to_column(var = "Sample_ID") %>%
  mutate(Sample_ID = stringr::str_replace(Sample_ID, "_[0-9]{1,2}", "")) %>%
  plyr::dlply(.variables = "Sample_ID") %>%
  lapply(function(x) {
    apply(x[, -1], 2, mean)
  }) %>%
  do.call(rbind, .)

data <-
  apply(data, 2, function(x) {
    (x - mean(x)) / sd(x)
  }) 

###correlate date
data_date_corrected <- 
data %>% 
  apply(2, function(x){
    x <- as.numeric(x)
    temp_data <- 
      data.frame(date = as.numeric(as.Date(sample_info$start_date) - as.Date("2016-01-15")), 
                 x, stringsAsFactors = FALSE)
    lm_result <- lm(formula = x ~ date, data = temp_data)
    lm_result$residuals
  })

rownames(data_date_corrected) <- rownames(data)

###correlate location
data_location_corrected <- 
  data %>% 
  apply(2, function(x){
    x <- as.numeric(x)
    temp_data <- 
      data.frame(location = sample_info$location, 
                 x, stringsAsFactors = FALSE)
    lm_result <- lm(formula = x ~ location, data = temp_data)
    lm_result$residuals
  })

rownames(data_location_corrected) <- rownames(data)

pca_object <-
  prcomp(x = data)

pca_object_date_corrected <-
  prcomp(x = data_date_corrected)

pca_object_location_corrected <-
  prcomp(x = data_location_corrected)

as.Date_origin <- function(x) {
  as.Date(x, origin = '1970-01-01')
}

###PCA for date
library(wesanderson)
names(wes_palettes)
wes_palette(name = "Zissou1", n = 100, type = "continuous")

pal <-
  wesanderson::wes_palette(name = "Zissou1", n = 100, type = "continuous")

plot <-
  pca_object$x %>%
  as.data.frame() %>%
  rownames_to_column(var = "sample_id") %>%
  left_join(sample_info, by = "sample_id") %>%
  mutate(date = as.integer(as.Date(start_date))) %>%
  mutate(date2 = as.character(start_date)) %>%
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
    legend.position = c(0, 1),
    legend.justification = c(0, 1),
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

plot <-
  pca_object$x %>%
  as.data.frame() %>%
  rownames_to_column(var = "sample_id") %>%
  left_join(sample_info, by = "sample_id") %>%
  mutate(date = as.character(start_date)) %>%
  mutate(location = stringr::str_replace(location, "\\(.{1,20}\\)", "") %>% stringr::str_trim()) %>%
  mutate(location = stringr::str_replace(location, "weekdays", "weekday")) %>%
  ggplot(aes(PC1, PC2, colour = location)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_point(shape = 21,
             size = 4,
             aes(fill = location),
             color = "black") +
  ggsci::scale_color_aaas() +
  ggrepel::geom_label_repel(aes(PC1, PC2, label = location)) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13),
    panel.grid.minor = element_blank(),
    legend.position = c(0, 0),
    legend.justification = c(0, 0),
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
#        file = "pca_location.pdf",
#        width = 7,
#        height = 7)


####heatmap
library(ComplexHeatmap)

data <-
  expression_data %>%
  select(-contains("Blank"))

data <- log(data + 1, 2)

data <-
  data %>%
  t() %>% 
  as.data.frame() %>%
  rownames_to_column(var = "Sample_ID") %>%
  mutate(Sample_ID = stringr::str_replace(Sample_ID, "_[0-9]{1,2}", "")) %>%
  plyr::dlply(.variables = "Sample_ID") %>%
  lapply(function(x) {
    apply(x[, -1], 2, mean)
  }) %>%
  do.call(rbind, .)

data <-
  apply(data, 2, function(x) {
    (x - mean(x)) / sd(x)
  }) 

data[, 1]

rownames(data) <- as.character(sample_info$start_date)

# sample_info <-
#   sample_info %>%
#   dplyr::filter(sample_id %in% rownames(data))

rownames(data) == sample_info$sample_id

range(data)

data[data > 2.56] <- 2.56

date <- as.numeric(as.Date(sample_info$start_date) - as.Date("2016-01-15"))

location_color <-
  ggsci::pal_aaas()(n = 10)[1:4]

names(location_color) <-
  c("Campus", "Houston", "Montana", "SF and Novato, CA")

ha1 = HeatmapAnnotation(
  location = sample_info$location,
  date = date,
  col = list(
    location = location_color,
    date = circlize::colorRamp2(
      breaks = c(
        min(date),
        mean(date),
        max(date)
      ),
      colors = c(wes_palettes$Zissou1[c(1, 3, 5)])
    )
    
  ),
  gp = gpar(col = "black"),
  annotation_name_side = c("left")
)

library(circlize)
wesanderson::wes_palette(name = names(wes_palettes)[19], 
                         type = "discrete")

scales::show_col(viridis::magma(n = 100))

col_fun = colorRamp2(
  breaks = seq(min(data), max(data), length.out = 90),
  colors =
    viridis::magma(n = 100)[-c(1:10)],
  transparency = 0
)

rownames(data) <-
  stringr::str_replace(rownames(data), "2016-", "")

plot <-
  Heatmap(
    matrix = t(data),
    show_row_names = FALSE,
    clustering_method_columns = "complete",
    clustering_method_rows = "complete",
    clustering_distance_rows = "euclidean",
    clustering_distance_columns = "euclidean",
    top_annotation = ha1,
    border = TRUE,
    col = col_fun,
    km = 2,
    # gap = gpar(col = "black"),
    rect_gp = gpar(col = "white", lwd = 0.05),
    column_names_gp = gpar(angle = 45), 
    column_names_rot = 45
  )

plot

library(ggfortify)
library(ggforce)
library(ggplotify)

plot <- as.ggplot(plot)

# ggsave(plot, filename = "heatmap.pdf", width = 7, height = 7)

###fuzzy c means clustering
library(Mfuzz)
#first get the time point data together:
# bind that to the data frame
temp_data <- exp_expression_data

##scale
temp_data <-
  t(data)
  
time <- c(1:11)

temp_data <- rbind(time, temp_data)

row.names(temp_data)[1] <- "time"

#save it to a temp file so ti doesnt clutter up my blog directory
tmp <- tempfile(tmpdir = "exposomeChemical_k_mean_clustering")

# write.table(
#   temp_data,
#   file = tmp,
#   sep = '\t',
#   quote = F,
#   col.names = NA
# )

#read it back in as an expression set
data <- table2eset(file = tmp)

# data.s <- standardise(data)
data.s <- data

m1 <- mestimate(data.s)
m1

Dmin(data.s, m=m1, crange=seq(2,22,1), repeats=3, visu=TRUE)

clust = 3

c <- mfuzz(data.s,c=clust,m=m1)
# save(c, file = file.path("exposomeChemical_k_mean_clustering","c"))
# load("c")
load("exposomeChemical_k_mean_clustering/c")

mfuzz.plot(
  eset = data.s,
  # min.mem = 0.6,
  cl = c,
  mfrow = c(2, 2),
  time.labels = time,
  new.window = FALSE
)

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

# xlsx::write.xlsx(cluster_info,
#                  file.path("exposomeChemical_k_mean_clustering/","cluster_info.xlsx"),
#                  row.names = FALSE)

decrease_class <- 
cluster_info %>%
  dplyr::filter(cluster == 2) %>%
  dplyr::select(variable_id) %>%
  dplyr::left_join(variable_info, by = c("variable_id" = "peak_ID")) %>% 
  dplyr::pull(Chemical_class)

increase_class <- 
  cluster_info %>%
  dplyr::filter(cluster != 2) %>%
  dplyr::select(variable_id) %>%
  dplyr::left_join(variable_info, by = c("variable_id" = "peak_ID")) %>% 
  dplyr::pull(Chemical_class)

chemical_class <- 
  variable_info$Chemical_class %>% unique() %>% 
  sort()

chemical_class <-
  c(chemical_class[-which(chemical_class == "Unknown")],
    "Unknown")

library(ggsci)

col <-
  c(ggsci::pal_igv()(15), "black")

names(col) <- chemical_class

plot1 <- 
  data.frame(class = increase_class, stringsAsFactors = FALSE) %>% 
  dplyr::group_by(class) %>% 
  dplyr::summarise(num = n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::arrange(num) %>% 
  dplyr::mutate(class = factor(class, levels = class)) %>% 
  dplyr::mutate(prop = num * 100/sum(num)) %>% 
  dplyr::mutate(lab.ypos = cumsum(prop) - 0.5*prop) %>%
  dplyr::mutate(class = factor(class, levels = rev(class))) %>% 
  ggplot(aes(x = 2, y = prop)) +
  geom_bar(
    stat = "identity",
    position = "stack",
    aes(fill = class),
    color = "black",
    show.legend = FALSE
  ) +
  geom_text(aes(y = lab.ypos, label = round(prop,1)), color = "white") + 
  scale_fill_manual(values = col) +
  theme_void() +
  xlim(0.5, 2.5) +
  ggplot2::coord_polar(theta = "y")



plot2 <- 
  data.frame(class = decrease_class, stringsAsFactors = FALSE) %>% 
  dplyr::group_by(class) %>% 
  dplyr::summarise(num = n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::arrange(num) %>% 
  dplyr::mutate(class = factor(class, levels = class)) %>% 
  dplyr::mutate(prop = num * 100/sum(num)) %>% 
  dplyr::mutate(lab.ypos = cumsum(prop) - 0.5*prop) %>%
  dplyr::mutate(class = factor(class, levels = rev(class))) %>% 
  ggplot(aes(x = 2, y = prop)) +
  geom_bar(
    stat = "identity",
    position = "stack",
    aes(fill = class),
    color = "black",
    show.legend = FALSE
  ) +
  geom_text(aes(y = lab.ypos, label = round(prop,1)), color = "white") + 
  scale_fill_manual(values = col) +
  theme_void() +
  xlim(0.5, 2.5) +
  ggplot2::coord_polar(theta = "y")


# ggsave(plot1,
#        filename = "plot_increase_plot.pdf",
#        width = 7,
#        height = 7)

# ggsave(plot2,
#        filename = "plot_decrease_plot.pdf",
#        width = 7,
#        height = 7)
