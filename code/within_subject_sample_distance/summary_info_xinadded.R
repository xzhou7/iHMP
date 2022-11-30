###no source
no_function()
masstools::setwd_project()
rm(list = ls())
library(tidyverse)

source("code/tools.R")
load("data/from_xin/physeq_clean.rda")
load("data_analysis/oral_microbiome/data_preparation/sample_info")

dir.create(
  "data_analysis/combine_microbiome/within_subject_distance_by_sample/summary/",
  recursive = TRUE
)

setwd("data_analysis/combine_microbiome/within_subject_distance_by_sample/summary/")

load("../oral/oral_braydist_by_sample")
load("../stool/stool_braydist_by_sample")
load("../nasal/nasal_braydist_by_sample")
load("../skin/skin_braydist_by_sample")

library(plyr)

save(stool_braydist_by_sample, stool_summary_info,
     oral_braydist_by_sample, oral_summary_info, 
     skin_braydist_by_sample, skin_summary_info, 
     nasal_braydist_by_sample, nasal_summary_info,
     file = "~/Library/CloudStorage/Box-Box/human_microbiome_project/Figures/Figure3/Figure3b_data/rawdata/raw.dist.RData")

####oral
oral_summary_info =
  oral_braydist_by_sample %>%
  plyr::dlply(.variables = .(subject_id1)) %>%
  purrr::map(function(x) {
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

oral_summary_info$p_value_adjust =
  p.adjust(oral_summary_info$p_value, method = "BH")

oral_summary_info <-
  oral_summary_info %>%
  dplyr::filter(p_value_adjust < 0.05)

####nasal
nasal_summary_info =
  nasal_braydist_by_sample %>%
  plyr::dlply(.variables = .(subject_id1)) %>%
  purrr::map(function(x) {
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

nasal_summary_info$p_value_adjust =
  p.adjust(nasal_summary_info$p_value, method = "BH")

nasal_summary_info <-
  nasal_summary_info %>%
  dplyr::filter(p_value_adjust < 0.05)

####stool
stool_summary_info =
  stool_braydist_by_sample %>%
  plyr::dlply(.variables = .(subject_id1)) %>%
  purrr::map(function(x) {
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

stool_summary_info$p_value_adjust =
  p.adjust(stool_summary_info$p_value, method = "BH")

stool_summary_info <-
  stool_summary_info %>%
  dplyr::filter(p_value_adjust < 0.05)

####skin
skin_summary_info =
  skin_braydist_by_sample %>%
  plyr::dlply(.variables = .(subject_id1)) %>%
  purrr::map(function(x) {
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

skin_summary_info$p_value_adjust =
  p.adjust(skin_summary_info$p_value, method = "BH")

skin_summary_info <-
  skin_summary_info %>%
  dplyr::filter(p_value_adjust < 0.05)

stool_summary_info <-
  stool_summary_info %>%
  dplyr::mutate(class = "stool")

skin_summary_info <-
  skin_summary_info %>%
  dplyr::mutate(class = "skin")

nasal_summary_info <-
  nasal_summary_info %>%
  dplyr::mutate(class = "nasal")

oral_summary_info <-
  oral_summary_info %>%
  dplyr::mutate(class = "oral")

temp_data =
  rbind(stool_summary_info,
        skin_summary_info,
        nasal_summary_info,
        oral_summary_info)

plot =
  temp_data %>%
  ggplot(aes(class, cor)) +
  geom_boxplot(aes(color = stringr::str_to_title(class)),
               outlier.shape = NA) +
  scale_color_manual(values = body_site_color) +
  guides(color = guide_none()) +
  ggnewscale::new_scale_color() +
  geom_point(aes(color = subject_id,
                 size = -log(p_value_adjust, 10)),
             alpha = 0.7) +
  guides(color = guide_none(),
         size = guide_legend(title = "-log(p_adjust, 10)")) +
  scale_size_continuous(range = c(1, 10)) +
  base_theme +
  labs(x = "", y = "Correlation")

plot

intersect_subject_id =
  Reduce(
    f = intersect,
    x = list(
      stool_summary_info$subject_id,
      skin_summary_info$subject_id,
      oral_summary_info$subject_id,
      nasal_summary_info$subject_id
    )
  )

plot

# ggsave(plot, filename = "summary_plot.pdf", width = 9, height = 7)



###########oral
oral_braydist_by_sample
library(plyr)
summary_info =
  oral_braydist_by_sample %>%
  plyr::dlply(.variables = .(subject_id1)) %>%
  purrr::map(function(x) {
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
  oral_braydist_by_sample %>%
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

temp_data_oral = temp_data

library(plyr)
oral_coefficients =
  temp_data_oral %>%
  plyr::dlply(.variables = .(name)) %>%
  purrr::map(function(x) {
    coefficients(lm(dist ~ diffdays, data = x))
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame() %>%
  colMeans()

###########nasal
nasal_braydist_by_sample
library(plyr)
summary_info =
  nasal_braydist_by_sample %>%
  plyr::dlply(.variables = .(subject_id1)) %>%
  purrr::map(function(x) {
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
  nasal_braydist_by_sample %>%
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

temp_data_nasal = temp_data

library(plyr)
nasal_coefficients =
  temp_data_nasal %>%
  plyr::dlply(.variables = .(name)) %>%
  purrr::map(function(x) {
    coefficients(lm(dist ~ diffdays, data = x))
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame() %>%
  colMeans()



###########stool
stool_braydist_by_sample
library(plyr)
summary_info =
  stool_braydist_by_sample %>%
  plyr::dlply(.variables = .(subject_id1)) %>%
  purrr::map(function(x) {
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

temp_data_stool = temp_data

library(plyr)
stool_coefficients =
  temp_data_stool %>%
  plyr::dlply(.variables = .(name)) %>%
  purrr::map(function(x) {
    coefficients(lm(dist ~ diffdays, data = x))
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame() %>%
  colMeans()


###########skin
skin_braydist_by_sample
library(plyr)
summary_info =
  skin_braydist_by_sample %>%
  plyr::dlply(.variables = .(subject_id1)) %>%
  purrr::map(function(x) {
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
  skin_braydist_by_sample %>%
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

temp_data_skin = temp_data

library(plyr)
skin_coefficients =
  temp_data_skin %>%
  plyr::dlply(.variables = .(name)) %>%
  purrr::map(function(x) {
    coefficients(lm(dist ~ diffdays, data = x))
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame() %>%
  colMeans()

temp =
  data.frame(x = 1000, y = 1)

temp_slope =
  rbind(stool_coefficients,
        skin_coefficients,
        nasal_coefficients,
        oral_coefficients) %>%
  as.data.frame()
colnames(temp_slope) = c("Intercept", "slope")
temp_slope$class = c("Stool", "Skin", "Nasal", "Oral")

plot =
  ggplot(aes(x, y), data = temp) +
  geom_point(shape = NA) +
  scale_x_continuous(limits = c(0, 1000)) +
  scale_y_continuous(limits = c(0.25, 1)) +
  base_theme +
  geom_abline(aes(
    slope = slope,
    intercept = Intercept,
    color = class
  ),
  data = temp_slope) +
  scale_color_manual(values = body_site_color) +
  labs(x = "Different days", y = "Bray-curtis Distance")
plot
# ggsave(plot,
#        filename = "different_days_vs_distance.pdf",
#        width = 8,
#        height = 7)
# save(temp_data_stool, file = "temp_data_stool")
# save(temp_data_skin, file = "temp_data_skin")
# save(temp_data_nasal, file = "temp_data_nasal")
# save(temp_data_oral, file = "temp_data_oral")






