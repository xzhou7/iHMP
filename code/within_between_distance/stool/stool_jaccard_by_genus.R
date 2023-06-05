###no source
setwd(masstools::get_project_wd())
rm(list = ls())
library(tidyverse)

source("code/tools.R")
load("data/from_xin/physeq_clean.rda")
load("data_analysis/stool_microbiome/data_preparation/sample_info")
family_info = readxl::read_xlsx("data/from_xin/Related People In iPOP_from.sophia.xlsx", 
                                col_names = FALSE) 
setwd("data_analysis/combine_microbiome/distance/stool/")

colnames(family_info) = c("family", "subject_id", "family_role")

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

physeq_stool = physeq_ST

variable_info =
    phyloseq::tax_table(physeq_stool) %>%
    as.data.frame()

sample_info =
    suppressWarnings(as_tibble(sample_data(physeq_stool))) %>%
    dplyr::select(RandomID:Date, batch) %>% 
    dplyr::left_join(family_info, by = c("SubjectID" = "subject_id"))

allGenus = as.vector(tax_table(physeq_stool)[, "Genus"])

uniqueGenus = unique(allGenus)

message(length(uniqueGenus), " genus to calculate")

####remain the distance == 1 which are not 0 in both two samples

library(future)
library(tictoc)
# Run sequentially
future::plan(multisession)

stool_jaccard_by_genus =
    1:length(uniqueGenus) %>%
    furrr::future_map(
        .f = function(i) {
            message(uniqueGenus[i], " ", i, " of ", length(uniqueGenus))

            temp_genus = uniqueGenus[i]
            subdat = prune_taxa(allGenus == temp_genus, physeq_stool)
            subdat2 = prune_taxa(taxa_sums(subdat) > 0, subdat)
            subdat3 = prune_samples(sample_sums(subdat2) > 0, subdat2)

            cormat =
                as.matrix(distance(subdat3, "jaccard"))

            cormat[lower.tri(cormat, diag = TRUE)] = as.numeric(NA)

            dist =
                cormat %>%
                as.data.frame() %>%
                tibble::rownames_to_column(var = "sample_id1") %>%
                tidyr::pivot_longer(
                    cols = -sample_id1,
                    names_to = "sample_id2",
                    values_to = "dist"
                ) %>%
                dplyr::filter(!is.na(dist))

            if (nrow(dist) == 0) {
                return(NULL)
            }

            dist =
                dist %>%
                dplyr::left_join(sample_info, by = c("sample_id1" = "RandomID")) %>%
                dplyr::rename(
                    SampleID1 = SampleID,
                    subject_id1 = SubjectID,
                    date1 = Date,
                    batch1 = batch,
                    family1 = family,
                    family_role1 = family_role
                ) %>%
                dplyr::left_join(sample_info, by = c("sample_id2" = "RandomID")) %>%
                dplyr::rename(
                    SampleID2 = SampleID,
                    subject_id2 = SubjectID,
                    date2 = Date,
                    batch2 = batch,
                    family2 = family,
                    family_role2 = family_role
                ) %>%
                dplyr::mutate(
                    type = case_when(
                        subject_id1 == subject_id2 ~ "within",
                        subject_id1 != subject_id2 ~ "between"
                    )
                ) %>%
                dplyr::mutate(type = case_when(
                    (type == "between") & (family1 == family2) ~ "family",
                    TRUE ~ type
                )) %>%
                dplyr::select(-c(sample_id1, sample_id2)) %>%
                dplyr::rename(sample_id1 = SampleID1,
                              sample_id2 = SampleID2) %>%
                dplyr::mutate(dataset = "stool",
                              genus = uniqueGenus[i]) %>%
                dplyr::mutate(
                    batch_type = ifelse(batch1 == batch2, "same_batch", "diff_batch"),
                    diffdays = abs(as.integer(date1 - date2)),
                    season1 = quarters(date1),
                    season2 = quarters(date2),
                    season_type = ifelse(season1 == season2, "in_season", "out_of_season")
                )

            dist =
                dist %>%
                dplyr::select(
                    subject_id1,
                    subject_id2,
                    sample_id1,
                    sample_id2,
                    date1,
                    date2,
                    dist,
                    dplyr::everything()
                )
        },
        .progress = TRUE
    ) %>%
    dplyr::bind_rows()

####remove the genus with within or between < 5
remain_genus = 
    stool_jaccard_by_genus %>% 
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
          dplyr::pull(Genus)
    ) %>% 
    unique()

remain_genus = remain_genus[!remain_genus %in% remove_name]

stool_jaccard_by_genus =
  stool_jaccard_by_genus %>%
  dplyr::filter(genus %in% remain_genus)

save(stool_jaccard_by_genus,
     file = "stool_jaccard_by_genus",
     compress = "xz")
