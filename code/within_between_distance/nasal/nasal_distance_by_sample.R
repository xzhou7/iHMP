###no source
masstools::setwd_project()
rm(list = ls())
library(tidyverse)

source("code/tools.R")
load("data/from_xin/physeq_clean.rda")
load("data_analysis/nasal_microbiome/data_preparation/sample_info")
load("data/from_xin/Diversity_Datatable.RData")
# family_info <-
#   readxl::read_xlsx("data/from_xin/Related People In iPOP_from.sophia.xlsx",
#                                 col_names = FALSE)
dir.create("data_analysis/combine_microbiome/distance/nasal_by_sample/")
setwd("data_analysis/combine_microbiome/distance/nasal_by_sample/")

# colnames(family_info) = c("family", "subject_id", "family_role")
#
# family_info =
#   family_info %>%
#   dplyr::mutate(subject_id = stringr::str_replace(subject_id, "1636\\-", ""))

subject_info <-
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

physeq_nasal = physeq_NS

variable_info =
  phyloseq::tax_table(physeq_nasal) %>%
  as.data.frame()

sample_info =
  suppressWarnings(as_tibble(sample_data(physeq_nasal))) %>%
  dplyr::select(RandomID:Date, batch)
# dplyr::left_join(family_info, by = c("SubjectID" = "subject_id"))


####remain the distance == 1 which are not 0 in both two samples

subdat2 = prune_taxa(taxa_sums(physeq_nasal) > 0, physeq_nasal)
subdat3 = prune_samples(sample_sums(subdat2) > 0, subdat2)

# cormat =
#   as.matrix(distance(subdat3, "bray"))
#
# dist =
#   cormat %>%
#   as.data.frame() %>%
#   tibble::rownames_to_column(var = "sample_id1") %>%
#   tidyr::pivot_longer(cols = -sample_id1,
#                       names_to = "sample_id2",
#                       values_to = "dist") %>%
#   dplyr::filter(!is.na(dist))
#
# save(dist, file = "dist")
load("dist")

dist <-
  dist %>%
  dplyr::left_join(
    stool.nasal.ASV.diversity %>%
      tibble::rownames_to_column(var = "sample_id") %>%
      dplyr::select(sample_id, depth),
    by = c("sample_id1" = "sample_id")
  ) %>%
  dplyr::rename(depth1 = depth) %>%
  dplyr::left_join(
    stool.nasal.ASV.diversity %>%
      tibble::rownames_to_column(var = "sample_id") %>%
      dplyr::select(sample_id, depth),
    by = c("sample_id2" = "sample_id")
  ) %>%
  dplyr::rename(depth2 = depth)

dist <-
  dist %>%
  dplyr::filter(sample_id1 != sample_id2)


dist %>%
  ggplot(aes(depth1, dist)) +
  geom_hex()

dist %>%
  ggplot(aes(depth2, dist)) +
  geom_hex()

dist %>%
  ggplot(aes(depth1 + depth2, dist)) +
  geom_hex()

cor.test(dist$depth1, dist$dist)
cor.test(dist$depth2, dist$dist)
cor.test(dist$depth1 + dist$depth2,
         dist$dist, method = "spearman")
