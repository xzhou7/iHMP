###no source
setwd(masstools::get_project_wd())
rm(list = ls())
library(tidyverse)

source("code/tools.R")
load("data/from_xin/physeq_clean.rda")
load("data_analysis/oral_microbiome/data_preparation/sample_info")

family_info <-
  readxl::read_xlsx("data/from_xin/Related People In iPOP_from.sophia.xlsx",
                    col_names = FALSE)
dir.create("data_analysis/combine_microbiome/distance/oral/permutation")
setwd("data_analysis/combine_microbiome/distance/oral/permutation")

colnames(family_info) <- 
  c("family", "subject_id", "family_role")

family_info =
  family_info %>%
  dplyr::mutate(subject_id = stringr::str_replace(subject_id, "1636\\-", ""))

subject_info =
  sample_info %>%
  dplyr::select(subject_id,
                IRIS,
                SSPG,
                FPG,
                SSPG.Date,
                Class,
                Gender,
                Ethnicity,
                Adj.age)

library(phyloseq)
library(tidyverse)
library(cowplot)
library(parallel)

ls()

physeq_oral = physeq_OR

variable_info =
  phyloseq::tax_table(physeq_oral) %>%
  as.data.frame()

sample_info =
  suppressWarnings(as_tibble(sample_data(physeq_oral))) %>%
  dplyr::rename(RandomID = KitID) %>%
  dplyr::select(RandomID:SampleID, batch:SubjectID) %>%
  dplyr::left_join(family_info, by = c("SubjectID" = "subject_id"))

allGenus = as.vector(tax_table(physeq_oral)[, "Genus"])

uniqueGenus = unique(allGenus)

message(length(uniqueGenus), " genus to calculate")

####remain the distance == 1 which are not 0 in both two samples

library(future)
library(tictoc)

load(
  here::here(
    "data_analysis/combine_microbiome/distance/oral/oral_braydist_by_genus"
  )
)

####remove the genus with within or between < 5
remain_genus =
  oral_braydist_by_genus %>%
  dplyr::group_by(genus, type) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::filter(n >= 5) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(genus) %>%
  dplyr::summarise(n = sum(type == "between") + sum(type == "within")) %>%
  dplyr::filter(n == 2) %>%
  dplyr::pull(genus)

remove_name =
  c("Zea",
    variable_info %>%
      dplyr::filter(Kingdom == "Archaea") %>%
      dplyr::pull(Genus)) %>%
  unique()

remain_genus = remain_genus[!remain_genus %in% remove_name]

oral_braydist_by_genus =
  oral_braydist_by_genus %>%
  dplyr::filter(genus %in% remain_genus)

#######
library(plyr)
library(tictoc)

genus_list <- unique(oral_braydist_by_genus$genus)

for (j in 31:50) {
  cat(j)
  x <- oral_braydist_by_genus %>%
    dplyr::filter(genus == genus_list[j])
  
  temp <-
    purrr::map(1:20, function(i) {
      cat(i, " ")
      idx <- sample(1:nrow(x), nrow(x), replace = TRUE) %>%
        unique() %>%
        sort()
      
      x <- x[idx, ]
      ####fold change
      between_mean1 = mean(x$dist[x$type == "between"])
      within_mean1 = mean(x$dist[x$type == "within"])
      family_mean1 = mean(x$dist[x$type == "family"])
      family_mean1 = ifelse(is.nan(family_mean1), NA, family_mean1)
      ####dist == 1 number
      one_info =
        x %>%
        dplyr::group_by(type) %>%
        dplyr::summarise(
          total_number = n(),
          number_1 = sum(dist == 1),
          freq = number_1 / total_number
        ) %>%
        dplyr::arrange(type) %>%
        tidyr::pivot_wider(names_from = type,
                           values_from = c(total_number, number_1, freq))
      
      if (ncol(one_info) < 9) {
        one_info =
          one_info %>%
          dplyr::mutate(
            total_number_family = NA,
            number_1_family = NA,
            freq_family = NA
          )
      }
      
      one_info =
        one_info %>%
        dplyr::select(
          total_number_between,
          total_number_family,
          total_number_within,
          number_1_between,
          number_1_family,
          number_1_within,
          freq_between,
          freq_family,
          freq_within
        )
      
      ###between within and family
      dist1 = x$dist[x$type == "between"]
      dist2 = x$dist[x$type == "within"]
      dist3 = x$dist[x$type == "family"]
      all_dist = list(dist1, dist2, dist3)
      number1 = length(dist1)
      number2 = length(dist2)
      number3 = length(dist3)
      all_number = list(number1, number2, number3)
      
      ###calcualte fc1_p
      if (number1 > number2) {
        ####between should larger than within
        if (number1 > 3 * number2) {
          temp =
            purrr::map(1:1000, function(i) {
              wilcox.test(
                x = sample(dist1, number2),
                y = dist2,
                alternative = "greater"
              )$p.value
            }) %>%
            unlist()
          fc1_p = median(temp)
        } else{
          fc1_p = wilcox.test(x = dist1,
                              y = dist2,
                              alternative = "greater")$p.value
        }
      } else{
        if (number1 * 3 < number2) {
          temp =
            purrr::map(1:1000, function(i) {
              wilcox.test(
                x = dist1,
                y = sample(dist2, number1),
                alternative = "greater"
              )$p.value
            }) %>%
            unlist()
          fc1_p = median(temp)
        } else{
          fc1_p = wilcox.test(x = dist1,
                              y = dist2,
                              alternative = "greater")$p.value
        }
      }
      
      #####calculate fc2 (family vs between)
      if (number3 > 0) {
        if (number1 > number3) {
          if (number1 > 3 * number3) {
            temp =
              purrr::map(1:1000, function(i) {
                wilcox.test(
                  x = sample(dist1, number3),
                  y = dist3,
                  alternative = "greater"
                )$p.value
              }) %>%
              unlist()
            ##family should less than between
            fc2_p = median(temp)
          } else{
            fc2_p = wilcox.test(x = dist1,
                                y = dist3,
                                alternative = "greater")$p.value
          }
        } else{
          if (number1 * 3 < number3) {
            temp =
              purrr::map(1:1000, function(i) {
                wilcox.test(
                  x = dist1,
                  y = sample(dist3, number1),
                  alternative = "greater"
                )$p.value
              }) %>%
              unlist()
            ##family should less than between
            fc2_p = median(temp)
          } else{
            fc2_p =  wilcox.test(x = dist1,
                                 y = dist3,
                                 alternative = "greater")$p.value
          }
        }
      } else{
        fc2_p = NA
      }
      
      #####calculate fc3 (family vs within)
      if (number3 > 0) {
        if (number2 > number3) {
          if (number2 > 3 * number3) {
            temp =
              purrr::map(1:1000, function(i) {
                mean(sample(dist2, number3))
                wilcox.test(
                  x = sample(dist2, number3),
                  y = dist3,
                  alternative = "less"
                )$p.value
              }) %>%
              unlist()
            ##family should higher than within
            fc3_p = median(temp)
          } else{
            fc3_p = wilcox.test(x = dist2,
                                y = dist3,
                                alternative = "less")$p.value
          }
          
        } else{
          if (number2 * 3 < number3) {
            temp =
              purrr::map(1:1000, function(i) {
                wilcox.test(
                  x = dist2,
                  y = sample(dist3, number2),
                  alternative = "less"
                )$p.value
              }) %>%
              unlist()
            ##family should higher than within
            ##temp is family
            fc3_p = median(temp)
          } else{
            fc3_p = wilcox.test(x = dist2,
                                y = dist3,
                                alternative = "less")$p.value
          }
        }
      } else{
        fc3_p = NA
      }
      
      data.frame(
        genus = x$genus[1],
        dataset = x$dataset[1],
        fc1_p,
        fc2_p,
        fc3_p,
        between_mean1,
        within_mean1,
        family_mean1,
        one_info
      )
    }) %>%
    dplyr::bind_rows()
  
  save(temp, file = paste(j))
  
}
