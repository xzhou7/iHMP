###https://biofam.github.io/MOFA2/index.html
##https://genomebiology.biomedcentral.com/track/pdf/10.1186/s13059-020-02015-1.pdf
no_source()

# set work directory
setwd(masstools::get_project_wd())
library(tidyverse)
library(phyloseq)
rm(list = ls())

source("code/tools.R")

####permutation
# for (i in 1:100) {
#   cat(i, " ")
#
#   ###load stool_microbiome
#   ###stool
#   {
#     load(
#       file.path(
#         "data_analysis/stool_microbiome/clustering/phylum_IS_permutation/",
#         paste("new_expression_data", i, sep = "_")
#       )
#     )
#     load(
#       file.path(
#         "data_analysis/stool_microbiome/clustering/phylum_IS_permutation",
#         paste("new_sample_info", i, sep = "_")
#       )
#     )
#     load(
#       file.path(
#         "data_analysis/stool_microbiome/clustering/phylum_IS_permutation",
#         paste("new_variable_info", i, sep = "_")
#       )
#     )
#
#     temp_data <-
#       unique(new_sample_info$class) %>%
#       purrr::map(function(x) {
#         new_expression_data[, which(new_sample_info$class == x)] %>%
#           apply(1, mean)
#       }) %>%
#       do.call(cbind, .) %>%
#       as.data.frame()
#
#     colnames(temp_data) <- unique(new_sample_info$class)
#
#     ###remove the tax that all 0 in all the samples
#     idx <-
#       apply(temp_data, 1, function(x) {
#         sum(x == 0) / ncol(temp_data)
#       }) %>%
#       `<`(0.5) %>%
#       which()
#
#     temp_data <-
#       temp_data[idx, ]
#
#     new_variable_info <-
#       new_variable_info[idx, ]
#
#     rownames(temp_data) == new_variable_info$variable_id
#
#     temp_data_stool <- temp_data
#     variable_info_stool <- new_variable_info
#   }
#
#   ###skin
#   {
#     load(
#       file.path(
#         "data_analysis/skin_microbiome/clustering/phylum_IS_permutation",
#         paste("new_expression_data", i, sep = "_")
#       )
#     )
#     load(
#       file.path(
#         "data_analysis/skin_microbiome/clustering/phylum_IS_permutation",
#         paste("new_sample_info", i, sep = "_")
#       )
#     )
#     load(
#       file.path(
#         "data_analysis/skin_microbiome/clustering/phylum_IS_permutation",
#         paste("new_variable_info", i, sep = "_")
#       )
#     )
#
#     temp_data <-
#       unique(new_sample_info$class) %>%
#       purrr::map(function(x) {
#         new_expression_data[, which(new_sample_info$class == x)] %>%
#           apply(1, mean)
#       }) %>%
#       do.call(cbind, .) %>%
#       as.data.frame()
#
#     colnames(temp_data) <- unique(new_sample_info$class)
#
#     ###remove the tax that all 0 in all the samples
#     idx <-
#       apply(temp_data, 1, function(x) {
#         sum(x == 0) / ncol(temp_data)
#       }) %>%
#       `<`(0.5) %>%
#       which()
#
#     temp_data <-
#       temp_data[idx, ]
#
#     new_variable_info <-
#       new_variable_info[idx, ]
#
#     rownames(temp_data) == new_variable_info$variable_id
#
#     temp_data_skin <- temp_data
#     variable_info_skin <- new_variable_info
#   }
#
#
#
#   ###nasal
#   {
#     load(
#       file.path(
#         "data_analysis/nasal_microbiome/clustering/phylum_IS_permutation",
#         paste("new_expression_data", i, sep = "_")
#       )
#     )
#     load(
#       file.path(
#         "data_analysis/nasal_microbiome/clustering/phylum_IS_permutation",
#         paste("new_sample_info", i, sep = "_")
#       )
#     )
#     load(
#       file.path(
#         "data_analysis/nasal_microbiome/clustering/phylum_IS_permutation",
#         paste("new_variable_info", i, sep = "_")
#       )
#     )
#
#     temp_data <-
#       unique(new_sample_info$class) %>%
#       purrr::map(function(x) {
#         new_expression_data[, which(new_sample_info$class == x)] %>%
#           apply(1, mean)
#       }) %>%
#       do.call(cbind, .) %>%
#       as.data.frame()
#
#     colnames(temp_data) <- unique(new_sample_info$class)
#
#     ###remove the tax that all 0 in all the samples
#     idx <-
#       apply(temp_data, 1, function(x) {
#         sum(x == 0) / ncol(temp_data)
#       }) %>%
#       `<`(0.5) %>%
#       which()
#
#     temp_data <-
#       temp_data[idx, ]
#
#     new_variable_info <-
#       new_variable_info[idx, ]
#
#     rownames(temp_data) == new_variable_info$variable_id
#
#     temp_data_nasal <- temp_data
#     variable_info_nasal <- new_variable_info
#   }
#
#
#   ###oral
#   {
#     load(
#       file.path(
#         "data_analysis/oral_microbiome/clustering/phylum_IS_permutation",
#         paste("new_expression_data", i, sep = "_")
#       )
#     )
#     load(
#       file.path(
#         "data_analysis/oral_microbiome/clustering/phylum_IS_permutation",
#         paste("new_sample_info", i, sep = "_")
#       )
#     )
#     load(
#       file.path(
#         "data_analysis/oral_microbiome/clustering/phylum_IS_permutation",
#         paste("new_variable_info", i, sep = "_")
#       )
#     )
#
#     temp_data <-
#       unique(new_sample_info$class) %>%
#       purrr::map(function(x) {
#         new_expression_data[, which(new_sample_info$class == x)] %>%
#           apply(1, mean)
#       }) %>%
#       do.call(cbind, .) %>%
#       as.data.frame()
#
#     colnames(temp_data) <- unique(new_sample_info$class)
#
#     ###remove the tax that all 0 in all the samples
#     idx <-
#       apply(temp_data, 1, function(x) {
#         sum(x == 0) / ncol(temp_data)
#       }) %>%
#       `<`(0.5) %>%
#       which()
#
#     temp_data <-
#       temp_data[idx, ]
#
#     new_variable_info <-
#       new_variable_info[idx, ]
#
#     rownames(temp_data) == new_variable_info$variable_id
#
#     temp_data_oral <- temp_data
#     variable_info_oral <- new_variable_info
#   }
#
#
#   ###proteome
#   {
#     load(file.path(
#       "data_analysis/proteome/clustering_IS_permutation",
#       paste("new_expression_data", i, sep = "_")
#     ))
#     load(file.path(
#       "data_analysis/proteome/clustering_IS_permutation",
#       paste("new_sample_info", i, sep = "_")
#     ))
#     load(file.path(
#       "data_analysis/proteome/clustering_IS_permutation",
#       paste("new_variable_info", i, sep = "_")
#     ))
#     temp_data <-
#       unique(new_sample_info$class) %>%
#       purrr::map(function(x) {
#         new_expression_data[, which(new_sample_info$class == x)] %>%
#           apply(1, mean)
#       }) %>%
#       do.call(cbind, .) %>%
#       as.data.frame()
#
#     colnames(temp_data) <- unique(new_sample_info$class)
#
#     ###remove the tax that all 0 in all the samples
#     idx <-
#       apply(temp_data, 1, function(x) {
#         sum(x == 0) / ncol(temp_data)
#       }) %>%
#       `<`(0.5) %>%
#       which()
#
#     temp_data <-
#       temp_data[idx, ]
#
#     new_variable_info <-
#       new_variable_info[idx, ]
#
#     rownames(temp_data) == new_variable_info$variable_id
#
#     temp_data_proteome <- temp_data
#     variable_info_proteome <- new_variable_info
#   }
#
#
#   ###metabolome
#   {
#     load(file.path(
#       "data_analysis/metabolome/clustering_IS_permutation",
#       paste("new_expression_data", i, sep = "_")
#     ))
#     load(file.path(
#       "data_analysis/metabolome/clustering_IS_permutation",
#       paste("new_sample_info", i, sep = "_")
#     ))
#     load(file.path(
#       "data_analysis/metabolome/clustering_IS_permutation",
#       paste("new_variable_info", i, sep = "_")
#     ))
#
#     temp_data <-
#       unique(new_sample_info$class) %>%
#       purrr::map(function(x) {
#         new_expression_data[, which(new_sample_info$class == x)] %>%
#           apply(1, mean)
#       }) %>%
#       do.call(cbind, .) %>%
#       as.data.frame()
#
#     colnames(temp_data) <- unique(new_sample_info$class)
#
#     ###remove the tax that all 0 in all the samples
#     idx <-
#       apply(temp_data, 1, function(x) {
#         sum(x == 0) / ncol(temp_data)
#       }) %>%
#       `<`(0.5) %>%
#       which()
#
#     temp_data <-
#       temp_data[idx, ]
#
#     new_variable_info <-
#       new_variable_info[idx, ]
#
#     rownames(temp_data) == new_variable_info$variable_id
#
#     temp_data_metabolome <- temp_data
#     variable_info_metabolome <- new_variable_info
#   }
#
#
#   ###lipidome
#   {
#     load(file.path(
#       "data_analysis/lipidome/clustering_IS_permutation",
#       paste("new_expression_data", i, sep = "_")
#     ))
#     load(file.path(
#       "data_analysis/lipidome/clustering_IS_permutation",
#       paste("new_sample_info", i, sep = "_")
#     ))
#     load(file.path(
#       "data_analysis/lipidome/clustering_IS_permutation",
#       paste("new_variable_info", i, sep = "_")
#     ))
#
#     temp_data <-
#       unique(new_sample_info$class) %>%
#       purrr::map(function(x) {
#         new_expression_data[, which(new_sample_info$class == x)] %>%
#           apply(1, mean)
#       }) %>%
#       do.call(cbind, .) %>%
#       as.data.frame()
#
#     colnames(temp_data) <- unique(new_sample_info$class)
#
#     ###remove the tax that all 0 in all the samples
#     idx <-
#       apply(temp_data, 1, function(x) {
#         sum(x == 0) / ncol(temp_data)
#       }) %>%
#       `<`(0.5) %>%
#       which()
#
#     temp_data <-
#       temp_data[idx, ]
#
#     new_variable_info <-
#       new_variable_info[idx, ]
#
#     rownames(temp_data) == new_variable_info$variable_id
#
#     temp_data_lipidome <- temp_data
#     variable_info_lipidome <- new_variable_info
#   }
#
#
#
#   variable_info_nasal <-
#     variable_info_nasal %>%
#     dplyr::mutate(variable_id = paste(variable_id, "nasal", sep = "_"))
#   rownames(temp_data_nasal) <-
#     variable_info_nasal$variable_id
#
#   variable_info_oral <-
#     variable_info_oral %>%
#     dplyr::mutate(variable_id = paste(variable_id, "oral", sep = "_"))
#   rownames(temp_data_oral) <-
#     variable_info_oral$variable_id
#
#   variable_info_skin <-
#     variable_info_skin %>%
#     dplyr::mutate(variable_id = paste(variable_id, "skin", sep = "_"))
#   rownames(temp_data_skin) <-
#     variable_info_skin$variable_id
#
#   variable_info_stool <-
#     variable_info_stool %>%
#     dplyr::mutate(variable_id = paste(variable_id, "stool", sep = "_"))
#   rownames(temp_data_stool) <-
#     variable_info_stool$variable_id
#
#
#
#   temp_data <-
#     rbind(
#       temp_data_nasal,
#       temp_data_oral,
#       temp_data_skin,
#       temp_data_stool,
#       temp_data_proteome,
#       temp_data_metabolome,
#       temp_data_lipidome
#     )
#
#   variable_info_lipidome <-
#     variable_info_lipidome %>%
#     dplyr::select(-membership)
#
#   variable_info <-
#     variable_info_nasal %>%
#     dplyr::full_join(variable_info_nasal, by = colnames(variable_info_nasal)) %>%
#     dplyr::full_join(variable_info_oral, by = colnames(variable_info_nasal)) %>%
#     dplyr::full_join(variable_info_skin, by = colnames(variable_info_nasal)) %>%
#     dplyr::full_join(variable_info_stool, by = colnames(variable_info_nasal)) %>%
#     dplyr::full_join(variable_info_proteome, by = "variable_id") %>%
#     dplyr::full_join(variable_info_metabolome, by = "variable_id") %>%
#     dplyr::full_join(variable_info_lipidome, by = "variable_id")
#
#   rownames(temp_data) == variable_info$variable_id
#   ###clustering
#   library(Mfuzz)
#   #first get the time point data together:
#   time <- c(1:5)
#
#   temp_data2 <-
#     apply(temp_data, 1, function(x) {
#       (x - mean(x)) / sd(x)
#     }) %>%
#     t() %>%
#     as.data.frame()
#
#   temp_data <- rbind(time, temp_data)
#
#   row.names(temp_data)[1] <- "time"
#
#   #save it to a temp file so ti doesnt clutter up my blog directory
#   dir.create("data_analysis/combine_omics/clustering_IS_permutation",
#              recursive = TRUE)
#   write.table(
#     temp_data,
#     file = "data_analysis/combine_omics/clustering_IS_permutation/temp_data.txt",
#     sep = '\t',
#     quote = FALSE,
#     col.names = NA
#   )
#
#   data <-
#     table2eset(filename = "data_analysis/combine_omics/clustering_IS_permutation/temp_data.txt")
#
#   data.s <- standardise(data)
#
#   m1 <- mestimate(data.s)
#   m1
#
#   clust = 3
#
#   c <- mfuzz(data.s, c = clust, m = m1)
#
#   cluster_info <-
#     data.frame(
#       variable_id = names(c$cluster),
#       c$membership,
#       cluster = c$cluster,
#       stringsAsFactors = FALSE
#     ) %>%
#     arrange(cluster)
#
#   openxlsx::write.xlsx(
#     x = cluster_info,
#     file = file.path(
#       "data_analysis/combine_omics/clustering_IS_permutation",
#       paste0("cluster_info_", i, ".xlsx")
#     ),
#     asTable = TRUE
#   )
#
# }

original_info <-
  readxl::read_xlsx("data_analysis/combine_omics/clustering_IS/cluster_info.xlsx") %>%
  as.data.frame()

info_permutation <- vector(mode = "list", length = 100)
for (i in 1:100) {
  cat(i, " ")
  info_permutation[[i]] <-
    readxl::read_xlsx(file.path(
      "data_analysis/combine_omics/clustering_IS_permutation/",
      paste0("cluster_info_", i, ".xlsx")
    )) %>%
    as.data.frame()
}

#####cluster1
cluster1 <-
  readxl::read_xlsx("data_analysis/combine_omics/clustering_IS/cluster_1/cluster1.xlsx")

matched_cluster <-
  lapply(info_permutation, function(x) {
    x %>%
      dplyr::filter(variable_id %in% cluster1$variable_id) %>%
      dplyr::count(cluster) %>%
      dplyr::filter(n == max(n)) %>%
      pull(cluster)
  }) %>%
  unlist()

temp_data <-
  purrr::map2(info_permutation,
              matched_cluster, function(x, y) {
                temp <-
                  cluster1[, c("variable_id")] %>%
                  dplyr::left_join(x[, c(1, y + 1, 5)],
                                   by = "variable_id") %>%
                  as.data.frame()
                colnames(temp)[2] <- "membership"
                temp
              })

random_cluster <-
  cluster1$variable_id %>%
  purrr::map(function(x) {
    cat(x, " ")
    temp <-
      lapply(temp_data, function(y) {
        y %>%
          dplyr::filter(variable_id == x) %>%
          pull(cluster)
      }) %>%
      unlist()
    temp
    # temp[!is.na(temp)]
  })

p_value <-
  purrr::map(random_cluster, function(x) {
    idx <- sample(1:length(matched_cluster), 100000, replace = TRUE)
    (1 - sum(x[idx] == matched_cluster[idx], na.rm = TRUE) / length(x[idx][!is.na(x[idx])]))
  }) %>%
  unlist()

p_value[p_value == 0] <-
  min(p_value[p_value != 0])

plot(cluster1$membership, -log(p_value, 10))
sum(p_value < 0.05)
cluster1$p_value <- p_value
openxlsx::write.xlsx(
  cluster1,
  "data_analysis/combine_omics/clustering_IS/cluster_1/cluster1_with_p_value.xlsx"
)






#####cluster2
cluster2 <-
  readxl::read_xlsx("data_analysis/combine_omics/clustering_IS/cluster_2/cluster2.xlsx")

matched_cluster <-
  lapply(info_permutation, function(x) {
    x %>%
      dplyr::filter(variable_id %in% cluster2$variable_id) %>%
      dplyr::count(cluster) %>%
      dplyr::filter(n == max(n)) %>%
      pull(cluster)
  }) %>%
  unlist()

temp_data <-
  purrr::map2(info_permutation,
              matched_cluster, function(x, y) {
                temp <-
                  cluster2[, c("variable_id")] %>%
                  dplyr::left_join(x[, c(1, y + 1, 5)],
                                   by = "variable_id") %>%
                  as.data.frame()
                colnames(temp)[2] <- "membership"
                temp
              })

random_cluster <-
  cluster2$variable_id %>%
  purrr::map(function(x) {
    cat(x, " ")
    temp <-
      lapply(temp_data, function(y) {
        y %>%
          dplyr::filter(variable_id == x) %>%
          pull(cluster)
      }) %>%
      unlist()
    temp
    # temp[!is.na(temp)]
  })

p_value <-
  purrr::map(random_cluster, function(x) {
    idx <- sample(1:length(matched_cluster), 100000, replace = TRUE)
    (1 - sum(x[idx] == matched_cluster[idx], na.rm = TRUE) / length(x[idx][!is.na(x[idx])]))
  }) %>%
  unlist()

p_value[p_value == 0] <-
  min(p_value[p_value != 0])

plot(cluster2$membership, -log(p_value, 10))
sum(p_value < 0.05)

cluster2$p_value <- p_value

openxlsx::write.xlsx(
  cluster2,
  "data_analysis/combine_omics/clustering_IS/cluster_2/cluster2_with_p_value.xlsx"
)


#####cluster3
cluster3 <-
  readxl::read_xlsx("data_analysis/combine_omics/clustering_IS/cluster_3/cluster3.xlsx")

matched_cluster <-
  lapply(info_permutation, function(x) {
    x %>%
      dplyr::filter(variable_id %in% cluster3$variable_id) %>%
      dplyr::count(cluster) %>%
      dplyr::filter(n == max(n)) %>%
      pull(cluster)
  }) %>%
  unlist()

temp_data <-
  purrr::map2(info_permutation,
              matched_cluster, function(x, y) {
                temp <-
                  cluster3[, c("variable_id")] %>%
                  dplyr::left_join(x[, c(1, y + 1, 5)],
                                   by = "variable_id") %>%
                  as.data.frame()
                colnames(temp)[2] <- "membership"
                temp
              })

random_cluster <-
  cluster3$variable_id %>%
  purrr::map(function(x) {
    cat(x, " ")
    temp <-
      lapply(temp_data, function(y) {
        y %>%
          dplyr::filter(variable_id == x) %>%
          pull(cluster)
      }) %>%
      unlist()
    temp
    # temp[!is.na(temp)]
  })

p_value <-
  purrr::map(random_cluster, function(x) {
    idx <- sample(1:length(matched_cluster), 100000, replace = TRUE)
    (1 - sum(x[idx] == matched_cluster[idx], na.rm = TRUE) / length(x[idx][!is.na(x[idx])]))
  }) %>%
  unlist()

p_value[p_value == 0] <-
  min(p_value[p_value != 0])

plot(cluster3$membership, -log(p_value, 10))
cor(cluster3$membership, -log(p_value, 10))
sum(p_value < 0.05)

cluster3$p_value <- p_value

openxlsx::write.xlsx(
  cluster3,
  "data_analysis/combine_omics/clustering_IS/cluster_3/cluster3_with_p_value.xlsx"
)
