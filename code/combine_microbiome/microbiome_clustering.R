###https://biofam.github.io/MOFA2/index.html
##https://genomebiology.biomedcentral.com/track/pdf/10.1186/s13059-020-02015-1.pdf
no_source()

# set work directory
setwd(masstools::get_project_wd())
library(tidyverse)
library(phyloseq)
rm(list = ls())

source("code/tools.R")

###load stool_microbiome
###stool
{
  load("data_analysis/stool_microbiome/clustering/phylum/new_expression_data")
  load("data_analysis/stool_microbiome/clustering/phylum/new_sample_info")
  load("data_analysis/stool_microbiome/clustering/phylum/new_variable_info")
  
  temp_data <-
    unique(new_sample_info$class) %>%
    purrr::map(function(x) {
      new_expression_data[, which(new_sample_info$class == x)] %>%
        apply(1, mean)
    }) %>%
    do.call(cbind, .) %>%
    as.data.frame()
  
  colnames(temp_data) <- unique(new_sample_info$class)
  
  ###remove the tax that all 0 in all the samples
  idx <-
    apply(temp_data, 1, function(x) {
      sum(x == 0) / ncol(temp_data)
    }) %>%
    `<`(0.5) %>%
    which()
  
  temp_data <-
    temp_data[idx, ]
  
  new_variable_info <-
    new_variable_info[idx, ]
  
  rownames(temp_data) == new_variable_info$variable_id
  
  temp_data_stool <- temp_data
  variable_info_stool <- new_variable_info
}

###skin
{
  load("data_analysis/skin_microbiome/clustering/phylum/new_expression_data")
  load("data_analysis/skin_microbiome/clustering/phylum/new_sample_info")
  load("data_analysis/skin_microbiome/clustering/phylum/new_variable_info")
  
  temp_data <-
    unique(new_sample_info$class) %>%
    purrr::map(function(x) {
      new_expression_data[, which(new_sample_info$class == x)] %>%
        apply(1, mean)
    }) %>%
    do.call(cbind, .) %>%
    as.data.frame()
  
  colnames(temp_data) <- unique(new_sample_info$class)
  
  ###remove the tax that all 0 in all the samples
  idx <-
    apply(temp_data, 1, function(x) {
      sum(x == 0) / ncol(temp_data)
    }) %>%
    `<`(0.5) %>%
    which()
  
  temp_data <-
    temp_data[idx, ]
  
  new_variable_info <-
    new_variable_info[idx, ]
  
  rownames(temp_data) == new_variable_info$variable_id
  
  temp_data_skin <- temp_data
  variable_info_skin <- new_variable_info
}



###nasal
{
  load("data_analysis/nasal_microbiome/clustering/phylum/new_expression_data")
  load("data_analysis/nasal_microbiome/clustering/phylum/new_sample_info")
  load("data_analysis/nasal_microbiome/clustering/phylum/new_variable_info")
  
  temp_data <-
    unique(new_sample_info$class) %>%
    purrr::map(function(x) {
      new_expression_data[, which(new_sample_info$class == x)] %>%
        apply(1, mean)
    }) %>%
    do.call(cbind, .) %>%
    as.data.frame()
  
  colnames(temp_data) <- unique(new_sample_info$class)
  
  ###remove the tax that all 0 in all the samples
  idx <-
    apply(temp_data, 1, function(x) {
      sum(x == 0) / ncol(temp_data)
    }) %>%
    `<`(0.5) %>%
    which()
  
  temp_data <-
    temp_data[idx, ]
  
  new_variable_info <-
    new_variable_info[idx, ]
  
  rownames(temp_data) == new_variable_info$variable_id
  
  temp_data_nasal <- temp_data
  variable_info_nasal <- new_variable_info
}


###oral
{
  load("data_analysis/oral_microbiome/clustering/phylum/new_expression_data")
  load("data_analysis/oral_microbiome/clustering/phylum/new_sample_info")
  load("data_analysis/oral_microbiome/clustering/phylum/new_variable_info")
  
  temp_data <-
    unique(new_sample_info$class) %>%
    purrr::map(function(x) {
      new_expression_data[, which(new_sample_info$class == x)] %>%
        apply(1, mean)
    }) %>%
    do.call(cbind, .) %>%
    as.data.frame()
  
  colnames(temp_data) <- unique(new_sample_info$class)
  
  ###remove the tax that all 0 in all the samples
  idx <-
    apply(temp_data, 1, function(x) {
      sum(x == 0) / ncol(temp_data)
    }) %>%
    `<`(0.5) %>%
    which()
  
  temp_data <-
    temp_data[idx, ]
  
  new_variable_info <-
    new_variable_info[idx, ]
  
  rownames(temp_data) == new_variable_info$variable_id
  
  temp_data_oral <- temp_data
  variable_info_oral <- new_variable_info
}

variable_info_nasal <-
  variable_info_nasal %>%
  dplyr::mutate(variable_id = paste(variable_id, "nasal", sep = "_"))
rownames(temp_data_nasal) <-
  variable_info_nasal$variable_id

variable_info_oral <-
  variable_info_oral %>%
  dplyr::mutate(variable_id = paste(variable_id, "oral", sep = "_"))
rownames(temp_data_oral) <-
  variable_info_oral$variable_id

variable_info_skin <-
  variable_info_skin %>%
  dplyr::mutate(variable_id = paste(variable_id, "skin", sep = "_"))
rownames(temp_data_skin) <-
  variable_info_skin$variable_id

variable_info_stool <-
  variable_info_stool %>%
  dplyr::mutate(variable_id = paste(variable_id, "stool", sep = "_"))
rownames(temp_data_stool) <-
  variable_info_stool$variable_id

temp_data <-
  rbind(temp_data_nasal,
        temp_data_oral,
        temp_data_skin,
        temp_data_stool)

variable_info <-
  rbind(variable_info_nasal,
        variable_info_oral,
        variable_info_skin,
        variable_info_stool)

rownames(temp_data) == variable_info$variable_id

#######work directory
setwd(masstools::get_project_wd())
dir.create("data_analysis/combine_microbiome/clustering/phylum",
           recursive = TRUE)
setwd("data_analysis/combine_microbiome/clustering/phylum")

###clustering
library(Mfuzz)
#first get the time point data together:
time <- c(1:5)

temp_data2 <-
  apply(temp_data, 1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()

temp_data <- rbind(time, temp_data)

row.names(temp_data)[1] <- "time"

#save it to a temp file so ti doesnt clutter up my blog directory
write.table(
  temp_data,
  file = "temp_data.txt",
  sep = '\t',
  quote = FALSE,
  col.names = NA
)

data <- table2eset(filename = "temp_data.txt")

data.s <- standardise(data)

m1 <- mestimate(data.s)
m1

plot <-
  Dmin(
    data.s,
    m = m1,
    crange = seq(2, 20, 1),
    repeats = 3,
    visu = TRUE
  )

plot <-
  plot %>%
  data.frame(distance = plot,
             k = seq(2, 20, 1)) %>%
  ggplot(aes(k, distance)) +
  geom_point(shape = 21,
             size = 4,
             fill = "black") +
  geom_smooth() +
  geom_segment(aes(
    x = k,
    y = 0,
    xend = k,
    yend = distance
  )) +
  theme_bw() +
  theme(
    # legend.position = c(0, 1),
    # legend.justification = c(0, 1),
    panel.grid = element_blank(),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA)
  ) +
  labs(x = "Cluster number",
       y = "Min. centroid distance") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

plot

ggsave(plot,
       filename = "distance_k_number.pdf",
       width = 7,
       height = 7)

clust = 5

c <- mfuzz(data.s, c = clust, m = m1)

mfuzz.plot(
  eset = data.s,
  min.mem = 0.5,
  cl = c,
  mfrow = c(4, 4),
  time.labels = time,
  new.window = FALSE
)

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

openxlsx::write.xlsx(x = cluster_info,
                     file = "cluster_info.xlsx", asTable = TRUE)

#####output the expression data of different clusters

##plot for each cluster
for (cluster_idx in 1:clust) {
  cat(cluster_idx, " ")
  dir.create(paste("cluster", cluster_idx, sep = "_"))
  cluster_data <-
    cluster_info %>%
    dplyr::filter(cluster_idx == cluster_idx) %>%
    dplyr::select(1, 1 + cluster_idx)
  
  colnames(cluster_data) <- c("variable_id", "membership")
  
  cluster_data <-
    cluster_data %>%
    dplyr::filter(membership > 0.5) %>%
    dplyr::left_join(variable_info, by = "variable_id")
  
  openxlsx::write.xlsx(
    x = cluster_data,
    file = file.path(
      paste("cluster", cluster_idx, sep = "_"),
      paste("cluster", cluster_idx, ".xlsx", sep = "")
    ),
    asTable = TRUE,
    overwrite = TRUE
  )
  
  ###cluster plot
  
  temp =
    temp_data2[cluster_data$variable_id,] %>%
    data.frame(
      membership = cluster_data$membership,
      .,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ) %>%
    tibble::rownames_to_column(var = "variable_id") %>%
    tidyr::pivot_longer(
      cols = -c(variable_id, membership),
      names_to = "sample_id",
      values_to = "value"
    ) %>%
    dplyr::mutate(sample_id = factor(
      sample_id,
      levels = c("negtaive_h", "ee", "el", "re", "positive_h")
    ))
  # dplyr::left_join(sample_info[, c("sample_id", "accurate_time")], by = "sample_id")
  
  plot <-
    ggplot() +
    geom_hline(yintercept = 0) +
    geom_line(aes(sample_id, value, group = variable_id, color = membership),
              data = temp) +
    scale_color_gradientn(colours = c(RColorBrewer::brewer.pal(n = 9, name = "OrRd")[c(1:7)])) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      legend.position = c(0, 1),
      legend.justification = c(0, 1),
      panel.grid.minor = element_blank(),
      axis.title = element_text(size = 13),
      axis.text = element_text(size = 12),
      # axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.background = element_rect(fill = "transparent", color = NA),
      plot.background = element_rect(fill = "transparent", color = NA),
      legend.background = element_rect(fill = "transparent", color = NA)
    ) +
    labs(
      x = "",
      y = "Z-score",
      title = paste(
        "Cluster ",
        cluster_idx,
        " (",
        nrow(cluster_data),
        " taxa)",
        sep = ""
      )
    )
  
  plot
  
  ggsave(
    plot,
    filename = file.path(
      paste("cluster", cluster_idx, sep = "_"),
      paste("cluster", cluster_idx, ".pdf", sep = "")
    ),
    width = 8,
    height = 6
  )
}
