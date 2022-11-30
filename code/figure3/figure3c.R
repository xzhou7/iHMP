masstools::setwd_project()
rm(list = ls())
library(tidyverse)

source("code/tools.R")

data <-
  readr::read_csv("Figures/Figure3/detailed/figure3c_data/beta.estimation.fourbodysite.csv")
setwd("Figures/Figure3/detailed/figure3c_data")
remain_idx <-
  apply(data, 1, function(x) {
    sum(is.na(x))
  }) %>%
  `<`(1) %>%
  which()

data <-
  data[remain_idx, ]

library(ComplexHeatmap)

temp_data <-
  data %>%
  dplyr::select(-IRIS.x) %>%
  tibble::column_to_rownames(var = "subject_id1")

ha1 <- HeatmapAnnotation(
  IRIS = data$IRIS.x,
  col = list(IRIS = c(iris_color, "Unknown" = "grey")),
  gp = gpar(col = "black")
)

library(circlize)

col_fun <-
  colorRamp2(breaks = c(-0.5, 0, 1),
             colors = c("darkblue", "white", "red"))

for (i in c("euclidean",
            "maximum",
            "manhattan",
            "canberra",
            "binary",
            "minkowski")) {
  for (j in c(
    "ward.D",
    "ward.D2",
    "single",
    "complete",
    "average",
    "mcquitty",
    "median",
    "centroid"
  )) {
    cat(paste(i, j, sep = "_"), " ")
    plot <-
      Heatmap(
        t(temp_data),
        name = paste(i, j, sep = "_"),
        top_annotation = ha1,
        clustering_distance_rows = i,
        clustering_distance_columns = i,
        clustering_method_rows = j,
        clustering_method_columns = j,
        border = TRUE,
        col = col_fun
      )
    
    plot <- ggplotify::as.ggplot(plot)
    ggsave(
      plot,
      height = 3,
      width = 7,
      filename = paste0(i, "_", j, ".pdf")
    )
    
  }
}

dim(temp_data)

cor(as.matrix(temp_data))


dist(t(as.matrix(temp_data)), method = "canberra")
