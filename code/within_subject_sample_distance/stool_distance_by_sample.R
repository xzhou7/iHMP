###no source
masstools::setwd_project()
rm(list = ls())
library(tidyverse)

source("code/tools.R")
load("data/from_xin/physeq_clean.rda")
load("data_analysis/stool_microbiome/data_preparation/sample_info")

dir.create(
  "data_analysis/combine_microbiome/within_subject_distance_by_sample/stool/",
  recursive = TRUE
)
setwd("data_analysis/combine_microbiome/within_subject_distance_by_sample/stool/")

# subject_info =
#     sample_info %>%
#     dplyr::select(subject_id,
#                   IRIS,
#                   SSPG,
#                   FPG,
#                   SSPG.Date,
#                   Class,
#                   Gender,
#                   Ethnicity,
#                   Adj.age)
# 
# library(phyloseq)
# library(tidyverse)
# library(cowplot)
# library(parallel)
# 
# ls()
# 
# physeq_stool = physeq_ST
# 
# variable_info =
#     phyloseq::tax_table(physeq_stool) %>%
#     as.data.frame()
# 
# sample_info =
#     suppressWarnings(as_tibble(sample_data(physeq_stool))) %>%
#     dplyr::select(RandomID:Date, batch)
# 
# ####remain the distance == 1 which are not 0 in both two samples
# 
# tax_table = as.data.frame(tax_table(physeq_stool))
# otu_table = otu_table(physeq_stool)
# sample_data = as.data.frame(sample_data(physeq_stool))
# 
# new_otu_table <-
# unique(tax_table$Genus) %>%
#   purrr::map(function(x){
#     idx = which(tax_table$Genus == x)
#     apply(otu_table[,idx], 1, sum)
#   }) %>%
#   do.call(cbind, .) %>%
#   as.data.frame()
# 
# colnames(new_otu_table) = unique(tax_table$Genus)
# 
# sample_info$RandomID == rownames(new_otu_table)
# 
# remain_subject_id = 
#   sample_info %>% 
#   dplyr::count(SubjectID) %>% 
#   dplyr::filter(n >= 5) %>% 
#   pull(SubjectID)
# 
# sample_info = 
#   sample_info %>% 
#   dplyr::filter(SubjectID %in% remain_subject_id)
# 
# new_otu_table = new_otu_table[sample_info$RandomID,]
# 
# dist =
# vegan::vegdist(x = as.matrix(new_otu_table)) %>%
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
# dist =
#   dist %>%
#   dplyr::left_join(sample_info, by = c("sample_id1" = "RandomID")) %>%
#   dplyr::rename(
#     SampleID1 = SampleID,
#     subject_id1 = SubjectID,
#     date1 = Date,
#     batch1 = batch
#   ) %>%
#   dplyr::left_join(sample_info, by = c("sample_id2" = "RandomID")) %>%
#   dplyr::rename(
#     SampleID2 = SampleID,
#     subject_id2 = SubjectID,
#     date2 = Date,
#     batch2 = batch
#   ) %>%
#   dplyr::mutate(
#     type = case_when(
#       subject_id1 == subject_id2 ~ "within",
#       subject_id1 != subject_id2 ~ "between"
#     )
#   ) %>%
#   dplyr::select(-c(sample_id1, sample_id2)) %>%
#   dplyr::rename(sample_id1 = SampleID1,
#                 sample_id2 = SampleID2) %>%
#   dplyr::mutate(dataset = "stool") %>%
#   dplyr::mutate(
#     batch_type = ifelse(batch1 == batch2, "same_batch", "diff_batch"),
#     diffdays = abs(as.integer(date1 - date2)),
#     season1 = quarters(date1),
#     season2 = quarters(date2),
#     season_type = ifelse(season1 == season2, "in_season", "out_of_season")
#   )
# 
# dist =
#   dist %>%
#   dplyr::select(
#     subject_id1,
#     subject_id2,
#     sample_id1,
#     sample_id2,
#     date1,
#     date2,
#     dist,
#     dplyr::everything()
#   ) %>%
#   dplyr::filter(type == "within")
# 
# stool_braydist_by_sample = dist
# 
# remain_subject_id =
#   stool_braydist_by_sample %>%
#   dplyr::count(subject_id1) %>%
#   dplyr::filter(n >= 5) %>%
#   pull(subject_id1)
# 
# stool_braydist_by_sample <-
#   stool_braydist_by_sample %>%
#   dplyr::filter(subject_id1 %in% remain_subject_id)
# 
# save(stool_braydist_by_sample,
#      file = "stool_braydist_by_sample",
#      compress = "xz")

load("stool_braydist_by_sample")

library(plyr)
summary_info = 
stool_braydist_by_sample %>% 
  plyr::dlply(.variables = .(subject_id1)) %>% 
  purrr::map(function(x){
    cor = cor.test(x$diffdays, x$dist, method = "spearman")$estimate
    p_value = cor.test(x$diffdays, x$dist, method = "spearman")$p.value
    
    data.frame(
      subject_id = x$subject_id1[1],
      number = nrow(x),
      cor = unname(cor),
      p_value = unname(p_value)
    )
    
  }) %>% 
  dplyr::bind_rows()
  
summary_info$p_value_adjust = 
  p.adjust(summary_info$p_value, method = "BH") 

summary_info <-
  summary_info %>%
  dplyr::filter(p_value_adjust < 0.05)

temp_data = 
stool_braydist_by_sample %>% 
  dplyr::left_join(summary_info, by = c("subject_id1" = "subject_id")) %>% 
  dplyr::filter(subject_id1 %in% summary_info$subject_id)

temp_data$name =
  paste(
    temp_data$subject_id1,
    temp_data$number,
    round(temp_data$cor, 2),
    round(temp_data$p_value_adjust, 4),
    sep = ","
  )

plot1 = 
temp_data %>% 
  ggplot(aes(diffdays, dist)) +
  # geom_point(aes(color = subject_id1), show.legend = FALSE) +
  geom_smooth(
    aes(color = name),
    method = "lm",
    show.legend = FALSE,
    se = FALSE
  ) +
  base_theme

plot1


stool_subject_id = 
unique(temp_data$subject_id1)
save(stool_subject_id, file = "stool_subject_id")

plot2 = 
temp_data %>% 
  ggplot(aes(diffdays, dist)) +
  geom_point(aes(color = name), 
             show.legend = FALSE) +
  geom_smooth(
   color = "black",
    method = "lm",
    show.legend = FALSE,
    se = FALSE
  ) +
  base_theme +
  facet_wrap(facets = vars(name), 
             scales = "free_x") +
  theme(axis.text = element_text(size = 10))
plot2

ggsave(plot1, filename = "total_plot.pdf", width = 9, height = 7)
ggsave(plot2, filename = "each_plot.pdf", width = 9, height = 7)


########IR vs IS
load(here::here("data/from_xin/ls"))
load(here::here("data_analysis/metabolome/data_preparation/sample_info"))

sample_info <-
  sample_info %>% 
  dplyr::select(subject_id, IRIS) %>% 
  dplyr::distinct(subject_id, .keep_all = TRUE) %>% 
  dplyr::filter(!is.na(IRIS)) %>% 
  dplyr::filter(IRIS != "Unknown")

########IR
summary_info_ir = 
summary_info %>%
  dplyr::left_join(sample_info, by = "subject_id") %>% 
  dplyr::filter(!is.na(IRIS)) %>% 
  dplyr::filter(IRIS == "IR")

plot1 = 
  temp_data %>% 
  dplyr::filter(subject_id1 %in% summary_info_ir$subject_id) %>% 
  ggplot(aes(diffdays, dist)) +
  geom_smooth(
    aes(color = name),
    method = "lm",
    show.legend = FALSE,
    se = FALSE
  ) +
  base_theme

plot1

plot2 = 
  temp_data %>% 
  dplyr::filter(subject_id1 %in% summary_info_ir$subject_id) %>% 
  ggplot(aes(diffdays, dist)) +
  geom_point(aes(color = name), 
             show.legend = FALSE) +
  geom_smooth(
    color = "black",
    method = "lm",
    show.legend = FALSE,
    se = FALSE
  ) +
  base_theme +
  facet_wrap(facets = vars(name), 
             scales = "free_x") +
  theme(axis.text = element_text(size = 10))
plot2

ggsave(plot1, filename = "total_plot_ir.pdf", width = 9, height = 7)
ggsave(plot2, filename = "each_plot_ir.pdf", width = 9, height = 7)





########IS
summary_info_is = 
  summary_info %>%
  dplyr::left_join(sample_info, by = "subject_id") %>% 
  dplyr::filter(!is.na(IRIS)) %>% 
  dplyr::filter(IRIS == "IS")

plot1 = 
  temp_data %>% 
  dplyr::filter(subject_id1 %in% summary_info_is$subject_id) %>% 
  ggplot(aes(diffdays, dist)) +
  geom_smooth(
    aes(color = name),
    method = "lm",
    show.legend = FALSE,
    se = FALSE
  ) +
  base_theme

plot1

plot2 = 
  temp_data %>% 
  dplyr::filter(subject_id1 %in% summary_info_is$subject_id) %>% 
  ggplot(aes(diffdays, dist)) +
  geom_point(aes(color = name), 
             show.legend = FALSE) +
  geom_smooth(
    color = "black",
    method = "lm",
    show.legend = FALSE,
    se = FALSE
  ) +
  base_theme +
  facet_wrap(facets = vars(name), 
             scales = "free_x") +
  theme(axis.text = element_text(size = 10))
plot2

ggsave(plot1, filename = "total_plot_is.pdf", width = 9, height = 7)
ggsave(plot2, filename = "each_plot_is.pdf", width = 9, height = 7)

plot = 
  rbind(summary_info_ir,
        summary_info_is) %>% 
  ggplot(aes(IRIS, cor)) +
  geom_boxplot(aes(color = IRIS), show.legend = FALSE) +
  geom_jitter(aes(color = IRIS, size = -log(p_value_adjust, 10)),
              show.legend = FALSE) +
  scale_fill_manual(values = iris_color) +
  scale_color_manual(values = iris_color) +
  base_theme +
  labs(x = "", y = "Correlation") +
  scale_size_continuous(range = c(1,10))

plot

wilcox.test(summary_info_ir$cor,
            summary_info_is$cor)

ggsave(plot, filename = "IRIS_cor.pdf", width = 7, height = 7)



