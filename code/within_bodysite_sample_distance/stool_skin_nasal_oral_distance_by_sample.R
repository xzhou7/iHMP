###no source
setwd(masstools::get_project_wd())
rm(list = ls())
library(tidyverse)
library(phyloseq)
library(cowplot)
library(parallel)

setwd(masstools::get_project_wd())

source("code/tools.R")
load("data/from_xin/physeq_clean.rda")
load("data_analysis/stool_microbiome/data_preparation/sample_info")
load("data/from_xin/ls")

physeq_stool = physeq_ST
physeq_skin = physeq_SK
physeq_oral = physeq_OR
physeq_nasal = physeq_NS

dir.create("data_analysis/combine_microbiome/within_bodysite_sample_distance/",
           recursive = TRUE)
setwd("data_analysis/combine_microbiome/within_bodysite_sample_distance/")

# #####stool
# stool_subject_info =
#   sample_info %>%
#   dplyr::select(subject_id,
#                 IRIS,
#                 SSPG,
#                 FPG,
#                 SSPG.Date,
#                 Class,
#                 Gender,
#                 Ethnicity,
#                 Adj.age)
#
# stool_variable_info =
#   phyloseq::tax_table(physeq_stool) %>%
#   as.data.frame()
#
# stool_sample_info =
#   suppressWarnings(as_tibble(sample_data(physeq_stool))) %>%
#   dplyr::select(RandomID:Date, batch)
#
# stool_tax_table = as.data.frame(tax_table(physeq_stool))
# stool_otu_table = otu_table(physeq_stool)
# stool_sample_data = as.data.frame(sample_data(physeq_stool))
#
# stool_otu_table <-
#   unique(stool_tax_table$Genus) %>%
#   purrr::map(function(x){
#     idx = which(stool_tax_table$Genus == x)
#     apply(stool_otu_table[,idx], 1, sum)
#   }) %>%
#   do.call(cbind, .) %>%
#   as.data.frame()
#
# colnames(stool_otu_table) = unique(stool_tax_table$Genus)
#
# stool_sample_info$RandomID == rownames(stool_otu_table)
#
#
# #####skin
# skin_subject_info =
#   sample_info %>%
#   dplyr::select(subject_id,
#                 IRIS,
#                 SSPG,
#                 FPG,
#                 SSPG.Date,
#                 Class,
#                 Gender,
#                 Ethnicity,
#                 Adj.age)
#
# skin_variable_info =
#   phyloseq::tax_table(physeq_skin) %>%
#   as.data.frame()
#
# skin_sample_info =
#     suppressWarnings(as_tibble(sample_data(physeq_skin))) %>%
#     dplyr::rename(RandomID = KitID) %>%
#     dplyr::select(RandomID:SampleID, batch:SubjectID)
#
# skin_tax_table = as.data.frame(tax_table(physeq_skin))
# skin_otu_table = otu_table(physeq_skin)
# skin_sample_data = as.data.frame(sample_data(physeq_skin))
#
# skin_otu_table <-
#   unique(skin_tax_table$Genus) %>%
#   purrr::map(function(x){
#     idx = which(skin_tax_table$Genus == x)
#     apply(skin_otu_table[,idx], 1, sum)
#   }) %>%
#   do.call(cbind, .) %>%
#   as.data.frame()
#
# colnames(skin_otu_table) = unique(skin_tax_table$Genus)
#
# skin_sample_info$RandomID == rownames(skin_otu_table)
#
# #####nasal
# nasal_subject_info =
#   sample_info %>%
#   dplyr::select(subject_id,
#                 IRIS,
#                 SSPG,
#                 FPG,
#                 SSPG.Date,
#                 Class,
#                 Gender,
#                 Ethnicity,
#                 Adj.age)
#
# nasal_variable_info =
#   phyloseq::tax_table(physeq_nasal) %>%
#   as.data.frame()
#
# nasal_sample_info =
#       suppressWarnings(as_tibble(sample_data(physeq_nasal))) %>%
#       dplyr::select(RandomID:Date, batch)
#
# nasal_tax_table = as.data.frame(tax_table(physeq_nasal))
# nasal_otu_table = otu_table(physeq_nasal)
# nasal_sample_data = as.data.frame(sample_data(physeq_nasal))
#
# nasal_otu_table <-
#   unique(nasal_tax_table$Genus) %>%
#   purrr::map(function(x){
#     idx = which(nasal_tax_table$Genus == x)
#     apply(nasal_otu_table[,idx], 1, sum)
#   }) %>%
#   do.call(cbind, .) %>%
#   as.data.frame()
#
# colnames(nasal_otu_table) = unique(nasal_tax_table$Genus)
#
# nasal_sample_info$RandomID == rownames(nasal_otu_table)
#
#
# #####oral
# oral_subject_info =
#   sample_info %>%
#   dplyr::select(subject_id,
#                 IRIS,
#                 SSPG,
#                 FPG,
#                 SSPG.Date,
#                 Class,
#                 Gender,
#                 Ethnicity,
#                 Adj.age)
#
# oral_variable_info =
#   phyloseq::tax_table(physeq_oral) %>%
#   as.data.frame()
#
# oral_sample_info =
#   suppressWarnings(as_tibble(sample_data(physeq_oral))) %>%
#   dplyr::rename(RandomID = KitID) %>%
#   dplyr::select(RandomID:SampleID, batch:SubjectID)
#
# oral_tax_table = as.data.frame(tax_table(physeq_oral))
# oral_otu_table = otu_table(physeq_oral)
# oral_sample_data = as.data.frame(sample_data(physeq_oral))
#
# oral_otu_table <-
#   unique(oral_tax_table$Genus) %>%
#   purrr::map(function(x){
#     idx = which(oral_tax_table$Genus == x)
#     apply(oral_otu_table[,idx], 1, sum)
#   }) %>%
#   do.call(cbind, .) %>%
#   as.data.frame()
#
# colnames(oral_otu_table) = unique(oral_tax_table$Genus)
#
# oral_sample_info$RandomID == rownames(oral_otu_table)
#
# ######match samples
# stool_sample_info <-
#   stool_sample_info %>%
#   dplyr::select(SampleID, SubjectID, RandomID, Date)
#
# skin_sample_info <-
#   skin_sample_info %>%
#   dplyr::select(SampleID, SubjectID, RandomID, Date)
#
# oral_sample_info <-
#   oral_sample_info %>%
#   dplyr::select(SampleID, SubjectID, RandomID, Date)
#
# nasal_sample_info <-
#   nasal_sample_info %>%
#   dplyr::select(SampleID, SubjectID, RandomID, Date)
#
# sample_info =
#  rbind(stool_sample_info %>%
#          dplyr::rename(sample_id = SampleID) %>%
#          dplyr::mutate(old_sample_id = sample_id) %>%
#          dplyr::mutate(sample_id = paste("stool", sample_id, sep = "_")),
#        skin_sample_info %>%
#          dplyr::rename(sample_id = SampleID) %>%
#          dplyr::mutate(old_sample_id = sample_id) %>%
#          dplyr::mutate(sample_id = paste("skin", sample_id, sep = "_")),
#        nasal_sample_info %>%
#          dplyr::rename(sample_id = SampleID) %>%
#          dplyr::mutate(old_sample_id = sample_id) %>%
#          dplyr::mutate(sample_id = paste("nasal", sample_id, sep = "_")),
#        oral_sample_info %>%
#          dplyr::rename(sample_id = SampleID) %>%
#          dplyr::mutate(old_sample_id = sample_id) %>%
#          dplyr::mutate(sample_id = paste("oral", sample_id, sep = "_")))
#
# stool_otu_table <-
# as.data.frame(t(stool_otu_table))
#
# colnames(stool_otu_table) =
#   paste("stool",stool_sample_info$SampleID, sep = "_")
#
# skin_otu_table <-
#   as.data.frame(t(skin_otu_table))
#
# colnames(skin_otu_table) =
#   paste("skin",skin_sample_info$SampleID, sep = "_")
#
# nasal_otu_table <-
#   as.data.frame(t(nasal_otu_table))
#
# colnames(nasal_otu_table) =
#   paste("nasal",nasal_sample_info$SampleID, sep = "_")
#
# oral_otu_table <-
#   as.data.frame(t(oral_otu_table))
#
# colnames(oral_otu_table) =
#   paste("oral",oral_sample_info$SampleID, sep = "_")
#
# union_variable_id =
#   Reduce(f = union, x = list(
#     rownames(stool_otu_table),
#     rownames(skin_otu_table),
#     rownames(nasal_otu_table),
#     rownames(oral_otu_table)
#   ))
#
# expression_data =
#   cbind(stool_otu_table[union_variable_id,],
#         skin_otu_table[union_variable_id,],
#         nasal_otu_table[union_variable_id,],
#         oral_otu_table[union_variable_id,])
#
# dim(expression_data)
# dim(sample_info)
#
# expression_data <-
# expression_data[,sample_info$sample_id]
#
# colnames(expression_data) == sample_info$sample_id
#
# expression_data[is.na(expression_data)] = 0
#
#
# dist =
#   vegan::vegdist(x = as.matrix(t(expression_data))) %>%
#   as.matrix()
#
# dist[lower.tri(dist, diag = TRUE)] = as.numeric(NA)
#
# dist =
#   dist %>%
#   as.data.frame() %>%
#   tibble::rownames_to_column(var = "sample_id1") %>%
#   tidyr::pivot_longer(
#     cols = -sample_id1,
#     names_to = "sample_id2",
#     values_to = "dist"
#   ) %>%
#   dplyr::filter(!is.na(dist))
#
# dist <-
#   dist %>%
#   dplyr::left_join(sample_info, by = c("sample_id1" = "sample_id")) %>%
#   dplyr::rename(
#     subject_id1 = SubjectID,
#     RandomID1 = RandomID,
#     date1 = Date,
#     old_sample_id1 = old_sample_id
#   ) %>%
#   dplyr::left_join(sample_info, by = c("sample_id2" = "sample_id")) %>%
#   dplyr::rename(
#     subject_id2 = SubjectID,
#     RandomID2 = RandomID,
#     date2 = Date,
#     old_sample_id2 = old_sample_id
#   )
#
#
# dist =
#   dist %>%
#   dplyr::filter(old_sample_id1 == old_sample_id2)
#
# dist <-
#   dist %>%
#   dplyr::mutate(class1 = stringr::str_extract(sample_id1, "stool|skin|nasal|oral")) %>%
#   dplyr::mutate(class2 = stringr::str_extract(sample_id2, "stool|skin|nasal|oral"))
#
# dist =
# dist %>%
#   dplyr::arrange(old_sample_id1)
#
#
# save(dist,
#      file = "dist",
#      compress = "xz")

load("dist")

dist <-
  dist %>%
  dplyr::select(subject_id1, old_sample_id1, dist, date1, class1, class2) %>%
  dplyr::rename(subject_id = subject_id1,
                sample_id = old_sample_id1,
                date = date1)

dist$name =
  apply(dist, 1, function(x) {
    x = paste(sort(as.character(x[5:6])), collapse = "_")
  })

dist

dist %>%
  dplyr::count(subject_id) %>%
  dplyr::arrange(desc(n))

library(plyr)

temp_data =
  dist %>%
  dplyr::left_join(ls[, c("SampleID", "CL4")], by = c("sample_id" = "SampleID")) %>%
  dplyr::filter(CL4 == "Healthy" |
                  CL4 == "Infection" | CL4 == "Infection_L") %>%
  dplyr::mutate(CL4 = case_when(CL4 == "Infection_L" ~ "Infection",
                                TRUE ~ CL4)) %>%
  dplyr::mutate(date = as.character(date)) %>%
  dplyr::mutate(dist = round(dist, 2)) %>%
  plyr::dlply(.variables = .(subject_id, name)) %>%
  purrr::map(function(x) {
    if (sum(x$CL4 == "Healthy") < 5 | sum(x$CL4 == "Infection") < 5) {
      return(NULL)
    } else{
      return(x)
    }
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()

library(plyr)

temp_data %>%
  plyr::dlply(.variables = .(name)) %>%
  purrr::map(function(x) {
    plot =
      x %>%
      ggplot(aes(CL4, dist)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(aes(color = CL4), show.legend = FALSE) +
      base_theme +
      scale_color_manual(values = infection_color) +
      facet_wrap(facets = vars(subject_id), scales = "free_y") +
      theme(
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        strip.text = element_text(size = 8)
      ) +
      labs(x = "")
    
    ggsave(
      plot,
      filename = file.path("box_plot/", paste0(x$name[1], ".pdf")),
      width = 11,
      height = 9
    )
    
  })
