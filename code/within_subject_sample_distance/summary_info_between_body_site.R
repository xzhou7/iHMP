###no source
no_function()
masstools::setwd_project()
rm(list = ls())
library(tidyverse)

source("code/tools.R")gg
load("data/from_xin/physeq_clean.rda")
load("data_analysis/oral_microbiome/data_preparation/sample_info")

load(
  "data_analysis/combine_microbiome/within_subject_distance_by_sample/stool/stool_subject_id"
)
load(
  "data_analysis/combine_microbiome/within_subject_distance_by_sample/skin/skin_subject_id"
)
load(
  "data_analysis/combine_microbiome/within_subject_distance_by_sample/nasal/nasal_subject_id"
)
load(
  "data_analysis/combine_microbiome/within_subject_distance_by_sample/oral/oral_subject_id"
)

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

stool_braydist_by_sample$name =
  stool_braydist_by_sample %>%
  apply(1, function(x) {
    paste(sort(x[c(3, 4)]), collapse = "_")
  })

skin_braydist_by_sample$name =
  skin_braydist_by_sample %>%
  apply(1, function(x) {
    paste(sort(x[c(3, 4)]), collapse = "_")
  })

oral_braydist_by_sample$name =
  oral_braydist_by_sample %>%
  apply(1, function(x) {
    paste(sort(x[c(3, 4)]), collapse = "_")
  })

nasal_braydist_by_sample$name =
  nasal_braydist_by_sample %>%
  apply(1, function(x) {
    paste(sort(x[c(3, 4)]), collapse = "_")
  })

######stool with skin
intersect_name =
  intersect(stool_braydist_by_sample$name,
            skin_braydist_by_sample$name)

temp_data =
  data.frame(
    stool = stool_braydist_by_sample$dist[match(intersect_name, stool_braydist_by_sample$name)],
    skin = skin_braydist_by_sample$dist[match(intersect_name, skin_braydist_by_sample$name)],
    diff_days = stool_braydist_by_sample$diffdays[match(intersect_name, stool_braydist_by_sample$name)],
    subject_id = stool_braydist_by_sample$subject_id1[match(intersect_name, stool_braydist_by_sample$name)]
  )

###remove subjects with samples < 5
remain_subject_id =
  temp_data %>%
  dplyr::count(subject_id) %>%
  dplyr::filter(n >= 10) %>%
  pull(subject_id)

plot =
  temp_data %>%
  dplyr::filter(subject_id %in% unique(c(stool_subject_id, skin_subject_id))) %>%
  dplyr::filter(subject_id %in% remain_subject_id) %>%
  ggplot(aes(stool, skin)) +
  geom_point(aes(color = diff_days)) +
  # geom_hex() +
  geom_smooth(method = "lm") +
  labs(x = "Stool", y = "Skin") +
  base_theme +
  scale_color_gradient(low = "darkblue", high = "red") +
  facet_wrap(facets = vars(subject_id), scales = "free") +
  theme(axis.text = element_text(size = 10))
plot
ggsave(plot,
       filename = "stool_vs_skin_sample_distance.pdf",
       width = 7,
       height = 7)

plot =
  temp_data %>%
  dplyr::filter(subject_id %in% unique(c(stool_subject_id, skin_subject_id))) %>%
  dplyr::filter(subject_id %in% remain_subject_id) %>%
  ggplot(aes(stool, skin)) +
  geom_point(aes(color = diff_days)) +
  # geom_hex() +
  geom_smooth(method = "lm") +
  labs(x = "Stool", y = "Skin") +
  base_theme +
  scale_color_gradient(low = "darkblue", high = "red") +
  # facet_wrap(facets = vars(subject_id), scales = "free") +
  theme(axis.text = element_text(size = 10))
plot
ggsave(plot,
       filename = "stool_vs_skin_sample_distance_all.pdf",
       width = 7,
       height = 7)





######stool with nasal
intersect_name =
  intersect(stool_braydist_by_sample$name,
            nasal_braydist_by_sample$name)

temp_data =
  data.frame(
    stool = stool_braydist_by_sample$dist[match(intersect_name, stool_braydist_by_sample$name)],
    nasal = nasal_braydist_by_sample$dist[match(intersect_name, nasal_braydist_by_sample$name)],
    diff_days = stool_braydist_by_sample$diffdays[match(intersect_name, stool_braydist_by_sample$name)],
    subject_id = stool_braydist_by_sample$subject_id1[match(intersect_name, stool_braydist_by_sample$name)]
  )


###remove subjects with samples < 5
remain_subject_id =
  temp_data %>%
  dplyr::count(subject_id) %>%
  dplyr::filter(n >= 10) %>%
  pull(subject_id)


plot =
  temp_data %>%
  dplyr::filter(subject_id %in% unique(c(stool_subject_id, nasal_subject_id))) %>%
  dplyr::filter(subject_id %in% remain_subject_id) %>%
  ggplot(aes(stool, nasal)) +
  geom_point(aes(color = diff_days)) +
  # geom_hex() +
  geom_smooth(method = "lm") +
  labs(x = "Stool", y = "nasal") +
  base_theme +
  scale_color_gradient(low = "darkblue", high = "red") +
  facet_wrap(facets = vars(subject_id), scales = "free") +
  theme(axis.text = element_text(size = 10))
plot
ggsave(plot,
       filename = "stool_vs_nasal_sample_distance.pdf",
       width = 7,
       height = 7)

plot =
  temp_data %>%
  dplyr::filter(subject_id %in% unique(c(stool_subject_id, nasal_subject_id))) %>%
  dplyr::filter(subject_id %in% remain_subject_id) %>%
  ggplot(aes(stool, nasal)) +
  geom_point(aes(color = diff_days)) +
  # geom_hex() +
  geom_smooth(method = "lm") +
  labs(x = "Stool", y = "nasal") +
  base_theme +
  scale_color_gradient(low = "darkblue", high = "red") +
  # facet_wrap(facets = vars(subject_id), scales = "free") +
  theme(axis.text = element_text(size = 10))
plot
ggsave(plot,
       filename = "stool_vs_nasal_sample_distance_all.pdf",
       width = 7,
       height = 7)





######stool with oral
intersect_name =
  intersect(stool_braydist_by_sample$name,
            oral_braydist_by_sample$name)

temp_data =
  data.frame(
    stool = stool_braydist_by_sample$dist[match(intersect_name, stool_braydist_by_sample$name)],
    oral = oral_braydist_by_sample$dist[match(intersect_name, oral_braydist_by_sample$name)],
    diff_days = stool_braydist_by_sample$diffdays[match(intersect_name, stool_braydist_by_sample$name)],
    subject_id = stool_braydist_by_sample$subject_id1[match(intersect_name, stool_braydist_by_sample$name)]
  )


###remove subjects with samples < 5
remain_subject_id =
  temp_data %>%
  dplyr::count(subject_id) %>%
  dplyr::filter(n >= 10) %>%
  pull(subject_id)


plot =
  temp_data %>%
  dplyr::filter(subject_id %in% unique(c(stool_subject_id, oral_subject_id))) %>%
  dplyr::filter(subject_id %in% remain_subject_id) %>%
  ggplot(aes(stool, oral)) +
  geom_point(aes(color = diff_days)) +
  # geom_hex() +
  geom_smooth(method = "lm") +
  labs(x = "Stool", y = "oral") +
  base_theme +
  scale_color_gradient(low = "darkblue", high = "red") +
  facet_wrap(facets = vars(subject_id), scales = "free") +
  theme(axis.text = element_text(size = 10))
plot
ggsave(plot,
       filename = "stool_vs_oral_sample_distance.pdf",
       width = 7,
       height = 7)

plot =
  temp_data %>%
  dplyr::filter(subject_id %in% unique(c(stool_subject_id, oral_subject_id))) %>%
  dplyr::filter(subject_id %in% remain_subject_id) %>%
  ggplot(aes(stool, oral)) +
  geom_point(aes(color = diff_days)) +
  # geom_hex() +
  geom_smooth(method = "lm") +
  labs(x = "Stool", y = "oral") +
  base_theme +
  scale_color_gradient(low = "darkblue", high = "red") +
  # facet_wrap(facets = vars(subject_id), scales = "free") +
  theme(axis.text = element_text(size = 10))
plot
ggsave(plot,
       filename = "stool_vs_oral_sample_distance_all.pdf",
       width = 7,
       height = 7)









######skin with nasal
intersect_name =
  intersect(skin_braydist_by_sample$name,
            nasal_braydist_by_sample$name)

temp_data =
  data.frame(
    skin = skin_braydist_by_sample$dist[match(intersect_name, skin_braydist_by_sample$name)],
    nasal = nasal_braydist_by_sample$dist[match(intersect_name, nasal_braydist_by_sample$name)],
    diff_days = skin_braydist_by_sample$diffdays[match(intersect_name, skin_braydist_by_sample$name)],
    subject_id = skin_braydist_by_sample$subject_id1[match(intersect_name, skin_braydist_by_sample$name)]
  )

###remove subjects with samples < 5
remain_subject_id =
  temp_data %>%
  dplyr::count(subject_id) %>%
  dplyr::filter(n >= 10) %>%
  pull(subject_id)


plot =
  temp_data %>%
  dplyr::filter(subject_id %in% unique(c(skin_subject_id, nasal_subject_id))) %>%
  dplyr::filter(subject_id %in% remain_subject_id) %>%
  ggplot(aes(skin, nasal)) +
  geom_point(aes(color = diff_days)) +
  # geom_hex() +
  geom_smooth(method = "lm") +
  labs(x = "skin", y = "nasal") +
  base_theme +
  scale_color_gradient(low = "darkblue", high = "red") +
  facet_wrap(facets = vars(subject_id), scales = "free") +
  theme(axis.text = element_text(size = 10))
plot
ggsave(plot,
       filename = "skin_vs_nasal_sample_distance.pdf",
       width = 7,
       height = 7)

plot =
  temp_data %>%
  dplyr::filter(subject_id %in% unique(c(skin_subject_id, nasal_subject_id))) %>%
  dplyr::filter(subject_id %in% remain_subject_id) %>%
  ggplot(aes(skin, nasal)) +
  geom_point(aes(color = diff_days)) +
  # geom_hex() +
  geom_smooth(method = "lm") +
  labs(x = "skin", y = "nasal") +
  base_theme +
  scale_color_gradient(low = "darkblue", high = "red") +
  # facet_wrap(facets = vars(subject_id), scales = "free") +
  theme(axis.text = element_text(size = 10))
plot
ggsave(plot,
       filename = "skin_vs_nasal_sample_distance_all.pdf",
       width = 7,
       height = 7)


######skin with oral
intersect_name =
  intersect(skin_braydist_by_sample$name,
            oral_braydist_by_sample$name)

temp_data =
  data.frame(
    skin = skin_braydist_by_sample$dist[match(intersect_name, skin_braydist_by_sample$name)],
    oral = oral_braydist_by_sample$dist[match(intersect_name, oral_braydist_by_sample$name)],
    diff_days = skin_braydist_by_sample$diffdays[match(intersect_name, skin_braydist_by_sample$name)],
    subject_id = skin_braydist_by_sample$subject_id1[match(intersect_name, skin_braydist_by_sample$name)]
  )


###remove subjects with samples < 5
remain_subject_id =
  temp_data %>%
  dplyr::count(subject_id) %>%
  dplyr::filter(n >= 10) %>%
  pull(subject_id)


plot =
  temp_data %>%
  dplyr::filter(subject_id %in% unique(c(skin_subject_id, oral_subject_id))) %>%
  dplyr::filter(subject_id %in% remain_subject_id) %>%
  ggplot(aes(skin, oral)) +
  geom_point(aes(color = diff_days)) +
  # geom_hex() +
  geom_smooth(method = "lm") +
  labs(x = "skin", y = "oral") +
  base_theme +
  scale_color_gradient(low = "darkblue", high = "red") +
  facet_wrap(facets = vars(subject_id), scales = "free") +
  theme(axis.text = element_text(size = 10))
plot
ggsave(plot,
       filename = "skin_vs_oral_sample_distance.pdf",
       width = 7,
       height = 7)

plot =
  temp_data %>%
  dplyr::filter(subject_id %in% unique(c(skin_subject_id, oral_subject_id))) %>%
  dplyr::filter(subject_id %in% remain_subject_id) %>%
  ggplot(aes(skin, oral)) +
  geom_point(aes(color = diff_days)) +
  # geom_hex() +
  geom_smooth(method = "lm") +
  labs(x = "skin", y = "oral") +
  base_theme +
  scale_color_gradient(low = "darkblue", high = "red") +
  # facet_wrap(facets = vars(subject_id), scales = "free") +
  theme(axis.text = element_text(size = 10))
plot
ggsave(plot,
       filename = "skin_vs_oral_sample_distance_all.pdf",
       width = 7,
       height = 7)









######nasal with oral
intersect_name =
  intersect(nasal_braydist_by_sample$name,
            oral_braydist_by_sample$name)

temp_data =
  data.frame(
    nasal = nasal_braydist_by_sample$dist[match(intersect_name, nasal_braydist_by_sample$name)],
    oral = oral_braydist_by_sample$dist[match(intersect_name, oral_braydist_by_sample$name)],
    diff_days = nasal_braydist_by_sample$diffdays[match(intersect_name, nasal_braydist_by_sample$name)],
    subject_id = nasal_braydist_by_sample$subject_id1[match(intersect_name, nasal_braydist_by_sample$name)]
  )


###remove subjects with samples < 5
remain_subject_id =
  temp_data %>%
  dplyr::count(subject_id) %>%
  dplyr::filter(n >= 10) %>%
  pull(subject_id)


plot =
  temp_data %>%
  dplyr::filter(subject_id %in% unique(c(nasal_subject_id, oral_subject_id))) %>%
  dplyr::filter(subject_id %in% remain_subject_id) %>%
  ggplot(aes(nasal, oral)) +
  geom_point(aes(color = diff_days)) +
  # geom_hex() +
  geom_smooth(method = "lm") +
  labs(x = "nasal", y = "oral") +
  base_theme +
  scale_color_gradient(low = "darkblue", high = "red") +
  facet_wrap(facets = vars(subject_id), scales = "free") +
  theme(axis.text = element_text(size = 10))
plot
ggsave(plot,
       filename = "nasal_vs_oral_sample_distance.pdf",
       width = 7,
       height = 7)

plot =
  temp_data %>%
  dplyr::filter(subject_id %in% unique(c(nasal_subject_id, oral_subject_id))) %>%
  dplyr::filter(subject_id %in% remain_subject_id) %>%
  ggplot(aes(nasal, oral)) +
  geom_point(aes(color = diff_days)) +
  # geom_hex() +
  geom_smooth(method = "lm") +
  labs(x = "nasal", y = "oral") +
  base_theme +
  scale_color_gradient(low = "darkblue", high = "red") +
  # facet_wrap(facets = vars(subject_id), scales = "free") +
  theme(axis.text = element_text(size = 10))
plot
ggsave(plot,
       filename = "nasal_vs_oral_sample_distance_all.pdf",
       width = 7,
       height = 7)
