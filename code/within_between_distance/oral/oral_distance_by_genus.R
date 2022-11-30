###no source
masstools::setwd_project()
rm(list = ls())
library(tidyverse)

source("code/tools.R")
load("data/from_xin/physeq_clean.rda")
load("data_analysis/oral_microbiome/data_preparation/sample_info")
family_info = readxl::read_xlsx("data/from_xin/Related People In iPOP_from.sophia.xlsx", 
                                col_names = FALSE) 
setwd("data_analysis/combine_microbiome/distance/oral/")

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
# # Run sequentially
# future::plan(multisession)
# 
# oral_braydist_by_genus =
#   1:length(uniqueGenus) %>%
#   furrr::future_map(
#     .f = function(i) {
#       message(uniqueGenus[i], " ", i, " of ", length(uniqueGenus))
# 
#       temp_genus = uniqueGenus[i]
#       subdat = prune_taxa(allGenus == temp_genus, physeq_oral)
#       subdat2 = prune_taxa(taxa_sums(subdat) > 0, subdat)
#       subdat3 = prune_samples(sample_sums(subdat2) > 0, subdat2)
# 
#       cormat =
#         as.matrix(distance(subdat3, "bray"))
# 
#       cormat[lower.tri(cormat, diag = TRUE)] = as.numeric(NA)
# 
#       dist =
#         cormat %>%
#         as.data.frame() %>%
#         tibble::rownames_to_column(var = "sample_id1") %>%
#         tidyr::pivot_longer(cols = -sample_id1,
#                             names_to = "sample_id2",
#                             values_to = "dist") %>%
#         dplyr::filter(!is.na(dist))
# 
#       if (nrow(dist) == 0) {
#         return(NULL)
#       }
# 
#       dist =
#         dist %>%
#         dplyr::left_join(sample_info, by = c("sample_id1" = "RandomID")) %>%
#         dplyr::rename(
#           SampleID1 = SampleID,
#           subject_id1 = SubjectID,
#           date1 = Date,
#           batch1 = batch,
#           family1 = family,
#           family_role1 = family_role
#         ) %>%
#         dplyr::left_join(sample_info, by = c("sample_id2" = "RandomID")) %>%
#         dplyr::rename(
#           SampleID2 = SampleID,
#           subject_id2 = SubjectID,
#           date2 = Date,
#           batch2 = batch,
#           family2 = family,
#           family_role2 = family_role
#         ) %>%
#         dplyr::mutate(type = case_when(
#           subject_id1 == subject_id2 ~ "within",
#           subject_id1 != subject_id2 ~ "between"
#         )) %>%
#         dplyr::mutate(type = case_when((type == "between") &
#                                          (family1 == family2) ~ "family",
#                                        TRUE ~ type)) %>%
#         dplyr::select(-c(sample_id1, sample_id2)) %>%
#         dplyr::rename(sample_id1 = SampleID1,
#                       sample_id2 = SampleID2) %>%
#         dplyr::mutate(dataset = "oral",
#                       genus = uniqueGenus[i]) %>%
#         dplyr::mutate(
#           batch_type = ifelse(batch1 == batch2, "same_batch", "diff_batch"),
#           diffdays = abs(as.integer(date1 - date2)),
#           season1 = quarters(date1),
#           season2 = quarters(date2),
#           season_type = ifelse(season1 == season2, "in_season", "out_of_season")
#         )
# 
#       dist =
#         dist %>%
#         dplyr::select(
#           subject_id1,
#           subject_id2,
#           sample_id1,
#           sample_id2,
#           date1,
#           date2,
#           dist,
#           dplyr::everything()
#         )
#     },
#     .progress = TRUE
#   ) %>%
#   dplyr::bind_rows()
# 
# save(oral_braydist_by_genus,
#      file = "oral_braydist_by_genus",
#      compress = "xz")

load(here::here("data_analysis/combine_microbiome/distance/oral/oral_braydist_by_genus"))

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
          dplyr::pull(Genus)
    ) %>% 
    unique()

remain_genus = remain_genus[!remain_genus %in% remove_name]

oral_braydist_by_genus =
  oral_braydist_by_genus %>%
  dplyr::filter(genus %in% remain_genus)

# #######
# library(plyr)
# library(tictoc)
# tic()
# personalized_score =
#   oral_braydist_by_genus %>%
#   plyr::dlply(.variables = .(genus)) %>%
#   furrr::future_map(function(x) {
#     cat(x$genus[1], " ")
#     ####fold change
#     between_mean1 = mean(x$dist[x$type == "between"])
#     within_mean1 = mean(x$dist[x$type == "within"])
#     family_mean1 = mean(x$dist[x$type == "family"])
#     family_mean1 = ifelse(is.nan(family_mean1), NA, family_mean1)
#     ####dist == 1 number
#     one_info =
#       x %>%
#       dplyr::group_by(type) %>%
#       dplyr::summarise(
#         total_number = n(),
#         number_1 = sum(dist == 1),
#         freq = number_1 / total_number
#       ) %>%
#       dplyr::arrange(type) %>%
#       tidyr::pivot_wider(names_from = type,
#                          values_from = c(total_number, number_1, freq))
# 
#     if (ncol(one_info) < 9) {
#       one_info =
#         one_info %>%
#         dplyr::mutate(
#           total_number_family = NA,
#           number_1_family = NA,
#           freq_family = NA
#         )
#     }
# 
#     one_info =
#       one_info %>%
#       dplyr::select(
#         total_number_between,
#         total_number_family,
#         total_number_within,
#         number_1_between,
#         number_1_family,
#         number_1_within,
#         freq_between,
#         freq_family,
#         freq_within
#       )
# 
#     ###between within and family
#     dist1 = x$dist[x$type == "between"]
#     dist2 = x$dist[x$type == "within"]
#     dist3 = x$dist[x$type == "family"]
#     all_dist = list(dist1, dist2, dist3)
#     number1 = length(dist1)
#     number2 = length(dist2)
#     number3 = length(dist3)
#     all_number = list(number1, number2, number3)
# 
#     ###calcualte fc1_p
#     if(number1 > number2){
#       ####between should larger than within
#       if(number1 > 3 * number2){
#         temp =
#           purrr::map(1:1000, function(i){
#             wilcox.test(x = sample(dist1, number2), y = dist2, alternative = "greater")$p.value
#           }) %>%
#           unlist()
#         fc1_p = median(temp)
#       }else{
#         fc1_p = wilcox.test(x = dist1, y = dist2, alternative = "greater")$p.value
#       }
#     }else{
#       if(number1 * 3 < number2) {
#         temp =
#           purrr::map(1:1000, function(i){
#             wilcox.test(x = dist1, y = sample(dist2, number1), alternative = "greater")$p.value
#           }) %>%
#           unlist()
#         fc1_p = median(temp)
#       }else{
#         fc1_p = wilcox.test(x = dist1, y = dist2, alternative = "greater")$p.value
#       }
#     }
# 
#     #####calculate fc2 (family vs between)
#     if (number3 > 0) {
#       if(number1 > number3){
#         if(number1 > 3 * number3){
#           temp =
#             purrr::map(1:1000, function(i){
#               wilcox.test(x = sample(dist1, number3), y = dist3, alternative = "greater")$p.value
#             }) %>%
#             unlist()
#           ##family should less than between
#           fc2_p = median(temp)
#         }else{
#           fc2_p = wilcox.test(x = dist1, y = dist3, alternative = "greater")$p.value
#         }
#       }else{
#         if(number1 * 3 < number3){
#           temp =
#             purrr::map(1:1000, function(i){
#               wilcox.test(x = dist1, y = sample(dist3, number1), alternative = "greater")$p.value
#             }) %>%
#             unlist()
#           ##family should less than between
#           fc2_p = median(temp)
#         }else{
#           fc2_p =  wilcox.test(x = dist1, y = dist3, alternative = "greater")$p.value
#         }
#       }
#     } else{
#       fc2_p = NA
#     }
# 
#     #####calculate fc3 (family vs within)
#     if (number3 > 0) {
#       if(number2 > number3){
#         if(number2 > 3 * number3){
#           temp =
#             purrr::map(1:1000, function(i){
#               mean(sample(dist2, number3))
#               wilcox.test(x = sample(dist2, number3),
#                           y = dist3,
#                           alternative = "less")$p.value
#             }) %>%
#             unlist()
#           ##family should higher than within
#           fc3_p = median(temp)
#         }else{
#           fc3_p = wilcox.test(x = dist2,
#                               y = dist3,
#                               alternative = "less")$p.value
#         }
# 
#       }else{
#         if(number2 * 3 < number3){
#           temp =
#             purrr::map(1:1000, function(i){
#               wilcox.test(x = dist2,
#                           y = sample(dist3, number2),
#                           alternative = "less")$p.value
#             }) %>%
#             unlist()
#           ##family should higher than within
#           ##temp is family
#           fc3_p = median(temp)
#         }else{
#           fc3_p = wilcox.test(x = dist2,
#                               y = dist3,
#                               alternative = "less")$p.value
#         }
#       }
#     } else{
#       fc3_p = NA
#     }
# 
#     data.frame(
#       genus = x$genus[1],
#       dataset = x$dataset[1],
#       fc1_p,
#       fc2_p,
#       fc3_p,
#       between_mean1,
#       within_mean1,
#       family_mean1,
#       one_info
#     )
#   }) %>%
#   dplyr::bind_rows()
# 
# toc()
# save(personalized_score, file = "personalized_score")

load("personalized_score")

personalized_score$fc1_p_adjust = p.adjust(personalized_score$fc1_p,
                                           method = "BH")

personalized_score$fc2_p_adjust = p.adjust(personalized_score$fc2_p,
                                           method = "BH")

personalized_score$fc3_p_adjust = p.adjust(personalized_score$fc3_p,
                                           method = "BH")

personalized_score$fc_mean1 =
  personalized_score$between_mean1 - personalized_score$within_mean1

plot = 
personalized_score %>% 
  dplyr::filter(!is.na(fc_mean1)) %>% 
    ggplot(aes(freq_between - freq_within, fc_mean1)) +
    geom_point() +
    geom_smooth(se = FALSE, method = "lm") +
    base_theme +
    labs(x = "Frequency of distance = 1 (Between - Withinin)",
         y = "Mean distance (Between - Within)")
plot
# ggsave(plot,
#        file = "personalization_score_comparason.pdf",
#        width = 7, height = 7)

plot = 
    personalized_score %>%
    dplyr::mutate(diff = between_mean1 - within_mean1) %>% 
    dplyr::mutate(within_mean1 = within_mean1 * -1) %>% 
    dplyr::arrange(desc(diff)) %>% 
    dplyr::mutate(genus = factor(genus, levels = unique(genus))) %>% 
    ggplot() +
    geom_hline(yintercept = 0) +
    geom_point(aes(x = genus, y = between_mean1),
               color = between_within_color["between"]) +
    geom_point(aes(x = genus, y = within_mean1), color = between_within_color["within"]) +
    geom_point(aes(x = genus, y = family_mean1), color = between_within_color["family"]) +
    geom_point(aes(x = genus, y = -1*family_mean1), color = between_within_color["family"]) +
    geom_segment(aes(
        x = genus,
        xend = genus,
        y = 0,
        yend = between_mean1
    ),color = between_within_color["between"]
    ) +
    geom_segment(aes(
        x = genus,
        xend = genus,
        y = 0,
        yend = within_mean1
    ),color = between_within_color["within"]
    ) +
    geom_point(aes(genus, diff)) +
    base_theme +
    theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(
            angle = 45,
            hjust = 1,
            vjust = 1,
            size = 5
        )
    ) 
plot

# ggsave(plot, filename = "persionzied_score_distributation.pdf", width = 20, height = 7)
personalized_score %>% 
  dplyr::filter(!is.na(number_1_family)) %>% 
  dplyr::select(number_1_between, number_1_within, number_1_family)

####how many subjects in the final dataset
head(oral_braydist_by_genus)
table(oral_braydist_by_genus$subject_id1)

####if the diff days affect the within and between distance
masstools::setwd_project()
setwd("data_analysis/combine_microbiome/distance/oral/")

all_genus = oral_braydist_by_genus$genus %>% 
    unique()

# save(all_genus, file = "all_genus")

# final_result =
#   purrr::map(1:length(unique(oral_braydist_by_genus$genus)), function(idx) {
#     cat(idx, " ")
#     temp_genus = unique(oral_braydist_by_genus$genus)[idx]
#     temp_data_between =
#       oral_braydist_by_genus %>%
#       dplyr::filter(genus == unique(oral_braydist_by_genus$genus)[idx]) %>%
#       dplyr::filter(dist != 1) %>%
#       dplyr::filter(type == "between")
# 
#     temp_data_within =
#       oral_braydist_by_genus %>%
#       dplyr::filter(genus == unique(oral_braydist_by_genus$genus)[idx]) %>%
#       dplyr::filter(dist != 1) %>%
#       dplyr::filter(type == "within")
# 
#     if (nrow(temp_data_between) < 5 |
#         nrow(temp_data_within) < 5) {
#       return(NULL)
#     }
# 
#     test =
#       cor.test(temp_data_between$dist,
#                temp_data_between$diffdays,
#                method = "spearman")
#     cor_between =
#       as.numeric(test$estimate)
# 
#     p_between_cor = test$p.value
# 
#     scatter_plot_between =
#       temp_data_between %>%
#       ggplot(aes(diffdays, dist)) +
#       geom_point(aes(diffdays, dist, color = subject_id1)) +
#       base_theme +
#       geom_smooth(method = "lm",
#                   se = FALSE,
#                   color = "black") +
#       labs(title = "Between")
# 
#     scatter_plot_between
# 
#     test =
#       cor.test(temp_data_within$dist,
#                temp_data_within$diffdays,
#                method = "spearman")
#     cor_within =
#       as.numeric(test$estimate)
# 
#     p_within_cor = test$p.value
# 
#     scatter_plot_within =
#       temp_data_within %>%
#       ggplot(aes(diffdays, dist)) +
#       geom_point(aes(diffdays, dist, color = subject_id1)) +
#       base_theme +
#       geom_smooth(method = "lm",
#                   se = FALSE,
#                   color = "black") +
#       labs(title = "Within")
# 
#     scatter_plot_within
# 
#     result =
#       data.frame(
#         genus = temp_genus,
#         cor_within = cor_within,
#         number_within = nrow(temp_data_within),
#         p_within_cor = p_within_cor,
#         cor_between = cor_between,
#         p_between_cor = p_between_cor,
#         number_between = nrow(temp_data_between)
#       )
# 
#     library(patchwork)
#     if (any(c(result$p_within_cor,
#               result$p_between_cor) < 0.05)) {
#       plot =
#         scatter_plot_within + scatter_plot_between
#       if (nrow(temp_data_between) < 1000 &
#           nrow(temp_data_within) < 1000) {
#         ggsave(
#           plot,
#           filename = file.path("plot", paste(temp_genus, "pdf", sep = ".")),
#           width = 15,
#           height = 7
#         )
#       }
#     }
#     result
#   })
# 
# final_result =
#   final_result %>%
#   dplyr::bind_rows()
# 
# rownames(final_result) = NULL
# 
# save(final_result, file = "final_result")

load("final_result")

final_result$p_adjust_within_cor = p.adjust(final_result$p_within_cor, method = "BH")
final_result$p_adjust_between_cor = p.adjust(final_result$p_between_cor, method = "BH")

final_result

head(final_result)

plot = 
final_result %>% 
  ggplot() +
  geom_hline(yintercept = 0, color = "red", linetype = 2) +
  geom_vline(xintercept = 0, color = "red", linetype = 2) +
  geom_point(aes(x = cor_within, y = -log(p_adjust_within_cor, 10),
                 size = number_within),
             alpha = 0.8) +
  base_theme +
  labs("Correlation (diff days vs distance)", y = "Density") +
  ggrepel::geom_text_repel(aes(cor_within, 
                               -log(p_adjust_within_cor, 10), 
                               label = ifelse(abs(cor_within) > 0.3 & p_adjust_within_cor < 0.05, genus, NA)))
plot
# ggsave(plot, filename = "within_cor_p.pdf", width = 9, height = 7)

plot = 
  final_result %>% 
  ggplot() +
  geom_hline(yintercept = 0, color = "red", linetype = 2) +
  geom_vline(xintercept = 0, color = "red", linetype = 2) +
  geom_point(aes(x = cor_between, y = -log(p_adjust_between_cor, 10),
                 size = number_between),
             alpha = 0.8) +
  base_theme +
  labs("Correlation (diff days vs distance)", y = "Density") +
  ggrepel::geom_text_repel(aes(cor_between, 
                               -log(p_adjust_between_cor, 10), 
                               label = ifelse(abs(cor_between) > 0.3 & p_adjust_between_cor < 0.05, genus, NA)))
plot
# ggsave(plot, filename = "between_cor_p.pdf", width = 9, height = 7)


# #####after run the permutataion test, then use this to calculate the p values for fc1 and fc2
# final_permutation_result =
#   purrr::map(
#     1:length(dir("permutation_test/")),
#     .f = function(i) {
#       cat(i, " ")
#       load(dir("permutation_test", full.names = TRUE)[i])
#       permutation_result %>%
#         dplyr::mutate(fc1 = between_mean1 - within_mean1) %>%
#         dplyr::mutate(fc2 = (between_mean1 - family_mean1) / (between_mean1 - within_mean1)) %>%
#         dplyr::select(genus, fc1, fc2)
#     }
#   ) %>%
#   dplyr::bind_rows() %>%
#   dplyr::arrange(genus)
# 
# unique_genus = unique(personalized_score$genus)
# 
# p_values =
#   purrr::map(unique_genus, function(temp_genus) {
#     cat(temp_genus, " ")
#     
#     between_mean1 = personalized_score$between_mean1[which(personalized_score$genus == temp_genus)]
#     within_mean1 = personalized_score$within_mean1[which(personalized_score$genus == temp_genus)]
#     family_mean1 = personalized_score$family_mean1[which(personalized_score$genus == temp_genus)]
#     
#     fc1 = between_mean1 - within_mean1
#     fc2 = (between_mean1 - family_mean1) / (between_mean1 - within_mean1)
#     
#     permutation_fc1 = final_permutation_result$fc1[final_permutation_result$genus == temp_genus]
#     permutation_fc2 = final_permutation_result$fc2[final_permutation_result$genus == temp_genus]
#     permutation_fc1 = permutation_fc1[!is.na(permutation_fc1)]
#     permutation_fc2 = permutation_fc2[!is.na(permutation_fc2)]
#     
#     fc1_p = sum(permutation_fc1 > fc1) / length(permutation_fc1)
#     fc2_p = sum(permutation_fc2 > fc2) / length(permutation_fc2)
#     data.frame(genus = temp_genus, fc1_p, fc2_p)
#   }) %>%
#   dplyr::bind_rows()
# 
# 
# p_values$fc1_p_adjust = p.adjust(p_values$fc1_p, method = "BH")
# p_values$fc2_p_adjust = p.adjust(p_values$fc2_p, method = "BH")
# 
# permutation_p_values = p_values
# 
# save(permutation_p_values, file = "permutation_p_values")
load("permutation_p_values")

########output the box plot for each genus
dir.create("boxplot")

# for (temp_genus in all_genus) {
#   cat(temp_genus, " ")
#   temp_data =
#   oral_braydist_by_genus %>%
#     dplyr::filter(genus == temp_genus) %>%
#     dplyr::mutate(type = factor(type, levels = c("within", "family", "between")))
# 
#   nrow(temp_data)
# 
#   p1 = permutation_p_values$fc1_p_adjust[permutation_p_values$genus == temp_genus]
#   p2 = permutation_p_values$fc2_p_adjust[permutation_p_values$genus == temp_genus]
# 
#   p3 =
#   personalized_score %>%
#     dplyr::filter(genus == temp_genus) %>%
#     dplyr::select(fc1_p_adjust)
# 
#   p4 =
#     personalized_score %>%
#     dplyr::filter(genus == temp_genus) %>%
#     dplyr::select(fc2_p_adjust)
# 
#   p5 =
#     personalized_score %>%
#     dplyr::filter(genus == temp_genus) %>%
#     dplyr::select(fc3_p_adjust)
# 
#   text =
#     paste(
#       paste("permutation_test (Between - within), ", p1),
#       paste("permutation_test (family), ", p2),
#       paste("wilcox test (between vs within)", p3),
#       paste("wilcox test (between vs family)", p4),
#       paste("wilcox test (within vs family)", p5),
#       sep = "\n"
#     )
# 
#   plot =
#     temp_data %>%
#     ggplot() +
#     geom_boxplot(aes(type, dist, color = type),
#                  outlier.shape = NA,
#                  show.legend = FALSE) +
#     scale_color_manual(values = between_within_color) +
#     theme_bw() +
#     annotate(
#       geom = "text",
#       x = -Inf,
#       y = Inf,
#       label = text,
#       hjust = 0,
#       vjust = 1
#     ) +
#     labs(x = "", y = "Bray-Curtis distance")
# 
#   if(all(as.character(unique(temp_data$type)) != "family")){
#     temp_data$x[temp_data$type == "within"] =
#       1 + sample(seq(-0.3, 0.3, length.out = 1000),
#                  sum(temp_data$type == "within"),
#                  replace = TRUE)
# 
#     if(sum(temp_data$type == "within") > 2000){
#       plot =
#       plot +
#         geom_hex(aes(x = x, y = dist, fill = stat(count)),
#                  alpha = 0.8,
#                  data = temp_data %>% dplyr::filter(type == "within"),
#                  show.legend = TRUE) +
#         scale_fill_gradient(low = "white", high = between_within_color["within"])
#     }else{
#       plot =
#         plot +
#         geom_point(aes(x = x, y = dist),
#                    alpha = 0.8,
#                    data = temp_data %>% dplyr::filter(type == "within"),
#                    show.legend = FALSE,
#                    color = between_within_color["within"])
# 
#     }
# 
#     temp_data$x[temp_data$type == "between"] =
#       2 + sample(seq(-0.3, 0.3, length.out = 1000),
#                  sum(temp_data$type == "between"),
#                  replace = TRUE)
# 
# 
#     if(sum(temp_data$type == "between") > 2000){
#       plot =
#         plot +
#         ggnewscale::new_scale_fill() +
#         geom_hex(aes(x = x, y = dist, fill = stat(count)),
#                  alpha = 0.8,
#                  data = temp_data %>% dplyr::filter(type == "between"),
#                  show.legend = TRUE) +
#         scale_fill_gradient(low = "white", high = between_within_color["between"]) +
#         ggnewscale::new_scale_fill()
#     }else{
#       plot =
#         plot +
#         geom_point(aes(x = x, y = dist),
#                    alpha = 0.8,
#                    data = temp_data %>% dplyr::filter(type == "between"),
#                    show.legend = FALSE,
#                    color = between_within_color["between"])
# 
#     }
# 
# 
#   }else{
#     temp_data$x[temp_data$type == "within"] =
#       1 + sample(seq(-0.3, 0.3, length.out = 1000),
#                  sum(temp_data$type == "within"),
#                  replace = TRUE)
# 
#     if(sum(temp_data$type == "within") > 2000){
#       plot =
#         plot +
#         geom_hex(
#           aes(x = x, y = dist, fill = stat(count)),
#           alpha = 0.8,
#           data = temp_data %>% dplyr::filter(type == "within"),
#           show.legend = TRUE
#         ) +
#         scale_fill_gradient(low = "white", high = between_within_color["within"])
#     }else{
#       plot =
#         plot +
#         geom_point(aes(x = x, y = dist),
#                    alpha = 0.8,
#                    data = temp_data %>% dplyr::filter(type == "within"),
#                    show.legend = FALSE,
#                    color = between_within_color["within"])
# 
#     }
# 
# 
#     temp_data$x[temp_data$type == "family"] =
#       2 + sample(seq(-0.3, 0.3, length.out = 1000),
#                  sum(temp_data$type == "family"),
#                  replace = TRUE)
# 
#     if(sum(temp_data$type == "family") > 2000){
#       plot =
#         plot +
#         ggnewscale::new_scale_fill() +
#         geom_hex(
#           aes(x = x, y = dist, fill = stat(count)),
#           alpha = 0.8,
#           data = temp_data %>% dplyr::filter(type == "family"),
#           show.legend = TRUE
#         ) +
#         scale_fill_gradient(low = "white", high = between_within_color["family"])
#     }else{
#       plot =
#         plot +
#         geom_point(aes(x = x, y = dist),
#                    alpha = 0.8,
#                    data = temp_data %>% dplyr::filter(type == "family"),
#                    show.legend = FALSE,
#                    color = between_within_color["family"])
# 
#     }
# 
# 
#     temp_data$x[temp_data$type == "between"] =
#       3 + sample(seq(-0.3, 0.3, length.out = 1000),
#                  sum(temp_data$type == "between"),
#                  replace = TRUE)
# 
#     if(sum(temp_data$type == "between") > 2000){
#       plot =
#         plot +
#         ggnewscale::new_scale_fill() +
#         geom_hex(
#           aes(x = x, y = dist, fill = stat(count)),
#           alpha = 0.8,
#           data = temp_data %>% dplyr::filter(type == "between"),
#           show.legend = TRUE
#         ) +
#         scale_fill_gradient(low = "white", high = between_within_color["between"])
#     }else{
#       plot =
#         plot +
#         geom_point(aes(x = x, y = dist),
#                    alpha = 0.8,
#                    data = temp_data %>% dplyr::filter(type == "between"),
#                    show.legend = FALSE,
#                    color = between_within_color["between"])
# 
#     }
#   }
# 
#   ggsave(
#     plot,
#     filename = file.path(
#       "boxplot",
#       paste(temp_genus %>%
#               stringr::str_replace("\\/", "_"), ".pdf", sep = "")
#     ),
#     width = 7,
#     height = 7
#   )
# }

