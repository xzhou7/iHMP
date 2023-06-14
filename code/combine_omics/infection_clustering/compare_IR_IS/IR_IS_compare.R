###https://biofam.github.io/MOFA2/index.html
##https://genomebiology.biomedcentral.com/track/pdf/10.1186/s13059-020-02015-1.pdf
no_source()

# set work directory
setwd(masstools::get_project_wd())
library(tidyverse)
library(phyloseq)
rm(list = ls())

source("code/tools.R")

###IS
original_info_IS <-
  readxl::read_xlsx("data_analysis/combine_omics/clustering_IS/cluster_info.xlsx") %>%
  as.data.frame()

info_permutation_IS <- vector(mode = "list", length = 100)
for (i in 1:100) {
  cat(i, " ")
  info_permutation_IS[[i]] <-
    readxl::read_xlsx(file.path(
      "data_analysis/combine_omics/clustering_IS_permutation/",
      paste0("cluster_info_", i, ".xlsx")
    )) %>%
    as.data.frame()
}


###IR
original_info_IR <-
  readxl::read_xlsx("data_analysis/combine_omics/clustering_IR/cluster_info.xlsx") %>%
  as.data.frame()

info_permutation_IR <- vector(mode = "list", length = 100)
for (i in 1:100) {
  cat(i, " ")
  info_permutation_IR[[i]] <-
    readxl::read_xlsx(file.path(
      "data_analysis/combine_omics/clustering_IR_permutation/",
      paste0("cluster_info_", i, ".xlsx")
    )) %>%
    as.data.frame()
}


#####IS
##cluster 1 down
##cluster 2 up
##
##IR
##cluster 1 up
##cluster 3 down


#####Down
#####IS
#####cluster1
cluster1_IS <-
  readxl::read_xlsx("data_analysis/combine_omics/clustering_IS/cluster_1/cluster1.xlsx")

matched_cluster1_IS <-
  lapply(info_permutation_IS, function(x) {
    x %>%
      dplyr::filter(variable_id %in% cluster1_IS$variable_id) %>%
      dplyr::count(cluster) %>%
      dplyr::filter(n == max(n)) %>%
      pull(cluster)
  }) %>%
  unlist()

temp_data1_IS <-
  purrr::map2(info_permutation_IS,
              matched_cluster1_IS, function(x, y) {
                temp <-
                  cluster1_IS[, c("variable_id")] %>%
                  dplyr::left_join(x[, c(1, y + 1, 5)],
                                   by = "variable_id") %>%
                  as.data.frame()
                colnames(temp)[2] <- "membership"
                temp
              })


random_cluster1_IS <-
  cluster1_IS$variable_id %>%
  purrr::map(function(x) {
    cat(x, " ")
    temp <-
      lapply(temp_data1_IS, function(y) {
        y %>%
          dplyr::filter(variable_id == x) %>%
          pull(cluster)
      }) %>%
      unlist()
    temp
    # temp[!is.na(temp)]
  })


number_value1_IS <-
  purrr::map(random_cluster1_IS, function(x) {
    sum(x == matched_cluster1_IS, na.rm = TRUE)
  }) %>%
  unlist()

names(number_value1_IS) <- 
  cluster1_IS$variable_id


#####cluster3
cluster3_IR <-
  readxl::read_xlsx("data_analysis/combine_omics/clustering_IR/cluster_3/cluster3.xlsx")

matched_cluster3_IR <-
  lapply(info_permutation_IR, function(x) {
    x %>%
      dplyr::filter(variable_id %in% cluster3_IR$variable_id) %>%
      dplyr::count(cluster) %>%
      dplyr::filter(n == max(n)) %>%
      pull(cluster)
  }) %>%
  unlist()

temp_data3_IR <-
  purrr::map2(info_permutation_IR,
              matched_cluster3_IR, function(x, y) {
                temp <-
                  cluster3_IR[, c("variable_id")] %>%
                  dplyr::left_join(x[, c(1, y + 1, 5)],
                                   by = "variable_id") %>%
                  as.data.frame()
                colnames(temp)[2] <- "membership"
                temp
              })

random_cluster3_IR <-
  cluster3_IR$variable_id %>%
  purrr::map(function(x) {
    cat(x, " ")
    temp <-
      lapply(temp_data3_IR, function(y) {
        y %>%
          dplyr::filter(variable_id == x) %>%
          pull(cluster)
      }) %>%
      unlist()
    temp
    # temp[!is.na(temp)]
  })

number_value3_IR <-
  purrr::map(random_cluster3_IR, function(x) {
    sum(x == matched_cluster3_IR, na.rm = TRUE)
  }) %>%
  unlist()

names(number_value3_IR) <- 
  cluster3_IR$variable_id


####compare
length(number_value1_IS)
length(number_value3_IR)

intersect(names(number_value1_IS),
          names(number_value3_IR))

temp_data_down <- 
data.frame(variable_id = names(number_value1_IS),
           number = number_value1_IS) %>% 
  dplyr::full_join(
    data.frame(variable_id = names(number_value3_IR),
               number = number_value3_IR),
    by = "variable_id"
  )

p_value <- 
apply(temp_data_down, 1, function(x){
  x = as.numeric(x[2:3])
  x[is.na(x)] <- 0
  temp1 = c(as.numeric(x[1]), 100 - as.numeric(x[1]))
  temp2 = c(as.numeric(x[2]), 100 - as.numeric(x[2]))
  temp <- rbind(temp1, temp2)
  chisq.test(temp)$p.value
})

temp_data_down$p_value <- 
  p_value

temp_data_down$number.x[is.na(temp_data_down$number.x)] <- 0
temp_data_down$number.y[is.na(temp_data_down$number.y)] <- 0

colnames(temp_data_down)[2:3] <-c("number_IS", "number_IR")

dir.create("data_analysis/combine_omics/infection_clustering_IR_IS_comparison")

write.csv(temp_data_down, "data_analysis/combine_omics/infection_clustering_IR_IS_comparison/down.csv", row.names = FALSE)












#####UP
#####IS
#####cluster2
cluster2_IS <-
  readxl::read_xlsx("data_analysis/combine_omics/clustering_IS/cluster_2/cluster2.xlsx")

matched_cluster2_IS <-
  lapply(info_permutation_IS, function(x) {
    x %>%
      dplyr::filter(variable_id %in% cluster2_IS$variable_id) %>%
      dplyr::count(cluster) %>%
      dplyr::filter(n == max(n)) %>%
      pull(cluster)
  }) %>%
  unlist()

temp_data2_IS <-
  purrr::map2(info_permutation_IS,
              matched_cluster2_IS, function(x, y) {
                temp <-
                  cluster2_IS[, c("variable_id")] %>%
                  dplyr::left_join(x[, c(1, y + 1, 5)],
                                   by = "variable_id") %>%
                  as.data.frame()
                colnames(temp)[2] <- "membership"
                temp
              })

random_cluster2_IS <-
  cluster2_IS$variable_id %>%
  purrr::map(function(x) {
    cat(x, " ")
    temp <-
      lapply(temp_data2_IS, function(y) {
        y %>%
          dplyr::filter(variable_id == x) %>%
          pull(cluster)
      }) %>%
      unlist()
    temp
    # temp[!is.na(temp)]
  })


number_value2_IS <-
  purrr::map(random_cluster2_IS, function(x) {
    sum(x == matched_cluster2_IS, na.rm = TRUE)
  }) %>%
  unlist()

names(number_value2_IS) <- 
  cluster2_IS$variable_id


#####cluster1
cluster1_IR <-
  readxl::read_xlsx("data_analysis/combine_omics/clustering_IR/cluster_1/cluster1.xlsx")

matched_cluster1_IR <-
  lapply(info_permutation_IR, function(x) {
    x %>%
      dplyr::filter(variable_id %in% cluster1_IR$variable_id) %>%
      dplyr::count(cluster) %>%
      dplyr::filter(n == max(n)) %>%
      pull(cluster)
  }) %>%
  unlist()

temp_data1_IR <-
  purrr::map2(info_permutation_IR,
              matched_cluster1_IR, function(x, y) {
                temp <-
                  cluster1_IR[, c("variable_id")] %>%
                  dplyr::left_join(x[, c(1, y + 1, 5)],
                                   by = "variable_id") %>%
                  as.data.frame()
                colnames(temp)[2] <- "membership"
                temp
              })

random_cluster1_IR <-
  cluster1_IR$variable_id %>%
  purrr::map(function(x) {
    cat(x, " ")
    temp <-
      lapply(temp_data1_IR, function(y) {
        y %>%
          dplyr::filter(variable_id == x) %>%
          pull(cluster)
      }) %>%
      unlist()
    temp
    # temp[!is.na(temp)]
  })

number_value1_IR <-
  purrr::map(random_cluster1_IR, function(x) {
    sum(x == matched_cluster1_IR, na.rm = TRUE)
  }) %>%
  unlist()

names(number_value1_IR) <- 
  cluster1_IR$variable_id

####compare
length(number_value2_IS)
length(number_value1_IR)

intersect(names(number_value2_IS),
          names(number_value1_IR))

temp_data_up <- 
  data.frame(variable_id = names(number_value2_IS),
             number = number_value2_IS) %>% 
  dplyr::full_join(
    data.frame(variable_id = names(number_value1_IR),
               number = number_value1_IR),
    by = "variable_id"
  )

p_value <- 
  apply(temp_data_up, 1, function(x){
    x = as.numeric(x[2:3])
    x[is.na(x)] <- 0
    temp1 = c(as.numeric(x[1]), 100 - as.numeric(x[1]))
    temp2 = c(as.numeric(x[2]), 100 - as.numeric(x[2]))
    temp <- rbind(temp1, temp2)
    chisq.test(temp)$p.value
  })

temp_data_up$p_value <- 
  p_value

temp_data_up$number.x[is.na(temp_data_up$number.x)] <- 0
temp_data_up$number.y[is.na(temp_data_up$number.y)] <- 0

colnames(temp_data_up)[2:3] <-c("number_IS", "number_IR")

dir.create("data_analysis/combine_omics/infection_clustering_IR_IS_comparison")

write.csv(temp_data_up, "data_analysis/combine_omics/infection_clustering_IR_IS_comparison/up.csv", row.names = FALSE)


