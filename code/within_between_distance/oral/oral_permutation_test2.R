##
# no_function()

setwd(masstools::get_project_wd())
rm(list = ls())
load("data/from_xin/physeq_clean.rda")
setwd("data_analysis/combine_microbiome/distance/oral/")
load("all_genus")
library(tidyverse)
library(phyloseq)
dir.create("permutation_test")
setwd("permutation_test/")

family_info = readxl::read_xlsx(here::here("data/from_xin/Related People In iPOP_from.sophia.xlsx"),
                                col_names = FALSE) 

colnames(family_info) = c("family", "subject_id", "family_role")

family_info =
  family_info %>%
  dplyr::mutate(subject_id = stringr::str_replace(subject_id, "1636\\-", ""))

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

# uniqueGenus = unique(allGenus)

#####-------------------------------------------------------------------------
####permutatation test analysis

temp.fun2 = function(idx,
                     permutation_result,
                     iteration) {
  x = permutation_result[[idx]]
  x2 =
    x %>%
    dplyr::filter(dist != 1)

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

      ####for x2
      ###wilcox test
      if (sum(x2$type == "between") >= 5 &
          sum(x2$type == "within") >= 5) {
        ####fold change
        between_mean2 = mean(x2$dist[x2$type == "between"])
        within_mean2 = mean(x2$dist[x2$type == "within"])
        family_mean2 = mean(x2$dist[x2$type == "family"])
        family_mean2 = ifelse(is.nan(family_mean2), NA, family_mean2)
      } else{
        between_mean2 = NA
        within_mean2 = NA
        family_mean2 = NA
      }

      data.frame(
        genus = x$genus[1],
        dataset = x$dataset[1],
        between_mean1,
        within_mean1,
        family_mean1,
        between_mean2,
        within_mean2,
        family_mean2,
        one_info,
        iteration = iteration
      )
}


library(BiocParallel)
if (masstools::get_os() == "windows") {
  bpparam = BiocParallel::SnowParam(workers = 4,
                                    progressbar = TRUE)
} else{
  bpparam = BiocParallel::MulticoreParam(workers = 4,
                                         progressbar = TRUE)
}

library(plyr)

for(iteration in 101:200){
  cat(iteration, "\n")
  library(dtplyr)
  library(data.table)
  library(future)
  
  final_dist =
    1:length(all_genus) %>%
    furrr::future_map(
      .f = function(i) {
        temp_genus = all_genus[i]
        subdat = prune_taxa(allGenus == temp_genus, physeq_oral)
        subdat2 = prune_taxa(taxa_sums(subdat) > 0, subdat)
        subdat3 = prune_samples(sample_sums(subdat2) > 0, subdat2)
        
        cormat =
          as.matrix(distance(subdat3, "bray"))
        
        cormat[lower.tri(cormat, diag = TRUE)] = as.numeric(NA)
        
        random_name = colnames(cormat)
        random_name = sample(random_name)
        colnames(cormat) =
          rownames(cormat) = random_name
        
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
          dplyr::mutate(type = case_when(
            subject_id1 == subject_id2 ~ "within",
            subject_id1 != subject_id2 ~ "between"
          )) %>%
          dplyr::mutate(type = case_when((type == "between") &
                                           (family1 == family2) ~ "family",
                                         TRUE ~ type)) %>%
          dplyr::select(-c(sample_id1, sample_id2)) %>%
          dplyr::rename(sample_id1 = SampleID1,
                        sample_id2 = SampleID2) %>%
          dplyr::mutate(dataset = "oral",
                        genus = all_genus[i]) %>%
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
  
  ###remove genus without > 5
  remain_genus = 
    final_dist %>% 
    dplyr::group_by(genus, type) %>% 
    dplyr::summarise(n = n()) %>% 
    dplyr::filter(n >= 5) %>% 
    dplyr::ungroup() %>% 
    dplyr::group_by(genus) %>% 
    dplyr::summarise(n = sum(type == "between") + sum(type == "within")) %>% 
    dplyr::filter(n == 2) %>% 
    dplyr::pull(genus)
  
  final_dist =
    final_dist %>%
    dplyr::filter(genus %in% remain_genus) 
  
  permutation_result =
    final_dist %>% 
    plyr::dlply(.variables = .(genus))
  
  permutation_result = 
    BiocParallel::bplapply(
      X = 1:length(unique(final_dist$genus)),
      FUN = temp.fun2,
      BPPARAM = bpparam,
      permutation_result = permutation_result,
      iteration = iteration
    ) %>% 
    dplyr::bind_rows()

  save(permutation_result, 
       file = file.path(".", paste("permutation_result", iteration, sep = "_")), 
       compress = "xz")
}

