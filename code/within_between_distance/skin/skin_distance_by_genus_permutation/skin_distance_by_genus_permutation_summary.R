###no source
setwd(masstools::get_project_wd())
rm(list = ls())
library(tidyverse)
library(tidyr)
library(plyr)
source("code/tools.R")

dir.create("data_analysis/combine_microbiome/distance/skin/permutation")
setwd("data_analysis/combine_microbiome/distance/skin/permutation")

load("../personalized_score")

personalized_score$fc1 = personalized_score$between_mean1 - personalized_score$within_mean1

library(phyloseq)
library(tidyverse)
library(cowplot)
library(parallel)

ls()

file_name <- dir()

personalized_score_permutation <-
  purrr::map(1:length(file_name), function(i) {
    load(paste(i))
    temp
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()

personalized_score_permutation <-
  personalized_score_permutation %>%
  plyr::dlply(.variables = .(genus)) %>%
  purrr::map(function(x) {
    x$fc1 = x$between_mean1 - x$within_mean1
    x$fc1_sd = sd(x$fc1)
    x %>%
      dplyr::distinct(genus, .keep_all = TRUE)
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()

personalized_score_permutation$genus == personalized_score$genus

plot(personalized_score_permutation$fc1,
     personalized_score$fc1)

cor.test(personalized_score_permutation$fc1,
    personalized_score$fc1)

plot(personalized_score_permutation$fc1,
     personalized_score_permutation$fc1_sd)
abline(0, 0.15)

plot(personalized_score_permutation$fc1,
     personalized_score_permutation$fc1_sd)

personalized_score_permutation_trim <- 
  filter(personalized_score_permutation, fc1 > 15 * fc1_sd)

plot(personalized_score_permutation_trim$fc1,
     personalized_score_permutation_trim$fc1_sd)

head(personalized_score_permutation$fc1)

save(personalized_score_permutation, file = "../personalized_score_permutation")
save(personalized_score_permutation_trim, file = "../personalized_score_permutation_trim")



