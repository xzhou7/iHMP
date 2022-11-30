




###
no_function()

masstools::setwd_project()
library(tidyverse)
library(phyloseq)
rm(list = ls())

####load raw data
source(here::here("code/tools.R"))

###stool
load("data_analysis/combine_microbiome/distance/stool/stool_braydist_by_genus")

####remove the genus with within or between < 5
remain_genus =
  stool_braydist_by_genus %>%
  dplyr::group_by(genus, type) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::filter(n >= 5) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(genus) %>%
  dplyr::summarise(n = sum(type == "between") + sum(type == "within")) %>%
  dplyr::filter(n == 2) %>%
  dplyr::pull(genus)

stool_braydist_by_genus =
  stool_braydist_by_genus %>%
  dplyr::filter(genus %in% remain_genus)

###skin
load("data_analysis/combine_microbiome/distance/skin/skin_braydist_by_genus")

####remove the genus with within or between < 5
remain_genus =
  skin_braydist_by_genus %>%
  dplyr::group_by(genus, type) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::filter(n >= 5) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(genus) %>%
  dplyr::summarise(n = sum(type == "between") + sum(type == "within")) %>%
  dplyr::filter(n == 2) %>%
  dplyr::pull(genus)

skin_braydist_by_genus =
  skin_braydist_by_genus %>%
  dplyr::filter(genus %in% remain_genus)

###nasal
load("data_analysis/combine_microbiome/distance/nasal/nasal_braydist_by_genus")

####remove the genus with within or between < 5
remain_genus =
  nasal_braydist_by_genus %>%
  dplyr::group_by(genus, type) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::filter(n >= 5) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(genus) %>%
  dplyr::summarise(n = sum(type == "between") + sum(type == "within")) %>%
  dplyr::filter(n == 2) %>%
  dplyr::pull(genus)

nasal_braydist_by_genus =
  nasal_braydist_by_genus %>%
  dplyr::filter(genus %in% remain_genus)

###oral
load("data_analysis/combine_microbiome/distance/oral/oral_braydist_by_genus")

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

oral_braydist_by_genus =
  oral_braydist_by_genus %>%
  dplyr::filter(genus %in% remain_genus)


oral_braydist_by_genus$type

library(gghalves)

setwd("data_analysis/combine_microbiome/distance/")

temp_data1 =
  stool_braydist_by_genus %>%
  dplyr::mutate(class = "Stool") %>%
  dplyr::mutate(type = factor(type, levels = c("within", "family", "between")))

temp_data1 %>%
  dplyr::mutate(is1 = dist == 1) %>%
  dplyr::count(type, is1)

temp_data1 %>%
  dplyr::filter(dist == 1 & type == "within") %>%
  pull(subject_id1) %>%
  table()

temp_data2 =
  skin_braydist_by_genus %>%
  dplyr::mutate(class = "Skin") %>%
  dplyr::mutate(type = factor(type, levels = c("within", "family", "between")))

temp_data2 %>%
  dplyr::mutate(is1 = dist == 1) %>%
  dplyr::count(type, is1)

temp_data3 =
  oral_braydist_by_genus %>%
  dplyr::mutate(class = "Oral") %>%
  dplyr::mutate(type = factor(type, levels = c("within", "family", "between")))

temp_data3 %>%
  dplyr::mutate(is1 = dist == 1) %>%
  dplyr::count(type, is1)

temp_data4 =
  nasal_braydist_by_genus %>%
  dplyr::mutate(class = "Nasal") %>%
  dplyr::mutate(type = factor(type, levels = c("within", "family", "between")))

temp_data4 %>%
  dplyr::mutate(is1 = dist == 1) %>%
  dplyr::count(type, is1)

###statistical test
stool_within_family_test <-
  wilcox.test(temp_data1$dist[which(temp_data1$type == "within")],
              temp_data1$dist[which(temp_data1$type == "family")])

stool_within_between_test <-
  wilcox.test(temp_data1$dist[which(temp_data1$type == "within")],
              temp_data1$dist[which(temp_data1$type == "between")])

stool_family_between_test <-
  wilcox.test(temp_data1$dist[which(temp_data1$type == "family")],
              temp_data1$dist[which(temp_data1$type == "between")])


skin_within_family_test <-
  wilcox.test(temp_data2$dist[which(temp_data2$type == "within")],
              temp_data2$dist[which(temp_data2$type == "family")])

skin_within_between_test <-
  wilcox.test(temp_data2$dist[which(temp_data2$type == "within")],
              temp_data2$dist[which(temp_data2$type == "between")])

skin_family_between_test <-
  wilcox.test(temp_data2$dist[which(temp_data2$type == "family")],
              temp_data2$dist[which(temp_data2$type == "between")])




oral_within_family_test <-
  wilcox.test(temp_data3$dist[which(temp_data3$type == "within")],
              temp_data3$dist[which(temp_data3$type == "family")])

oral_within_between_test <-
  wilcox.test(temp_data3$dist[which(temp_data3$type == "within")],
              temp_data3$dist[which(temp_data3$type == "between")])

oral_family_between_test <-
  wilcox.test(temp_data3$dist[which(temp_data3$type == "family")],
              temp_data3$dist[which(temp_data3$type == "between")])



nasal_within_family_test <-
  wilcox.test(temp_data4$dist[which(temp_data4$type == "within")],
              temp_data4$dist[which(temp_data4$type == "family")])

nasal_within_between_test <-
  wilcox.test(temp_data4$dist[which(temp_data4$type == "within")],
              temp_data4$dist[which(temp_data4$type == "between")])

nasal_family_between_test <-
  wilcox.test(temp_data4$dist[which(temp_data4$type == "family")],
              temp_data4$dist[which(temp_data4$type == "between")])


#####write test result
sink(file = "test_result.txt")
cat("Stool within vs family")
stool_within_family_test
cat("Stool within vs between")
stool_within_between_test
cat("Stool family vs between")
stool_family_between_test

cat("Skin within vs family")
skin_within_family_test
cat("Skin within vs between")
skin_within_between_test
cat("Skin family vs between")
skin_family_between_test

cat("Oral within vs family")
oral_within_family_test
cat("Oral within vs between")
oral_within_between_test
cat("Oral family vs between")
oral_family_between_test

cat("Nasal within vs family")
nasal_within_family_test
cat("Nasal within vs between")
nasal_within_between_test
cat("Nasal family vs between")
nasal_family_between_test

sink()
# unlink('test_result.txt')

stool_within_dist_quantile <-
  quantile(temp_data1$dist[temp_data1$type == "within"])
stool_within_dist_mean <-
  mean(temp_data1$dist[temp_data1$type == "within"])
stool_family_dist_quantile <-
  quantile(temp_data1$dist[temp_data1$type == "family"])
stool_family_dist_mean <-
  mean(temp_data1$dist[temp_data1$type == "family"])
stool_between_dist_quantile <-
  quantile(temp_data1$dist[temp_data1$type == "between"])
stool_between_dist_mean <-
  mean(temp_data1$dist[temp_data1$type == "between"])

skin_within_dist_quantile <-
  quantile(temp_data2$dist[temp_data2$type == "within"])
skin_within_dist_mean <-
  mean(temp_data2$dist[temp_data2$type == "within"])
skin_family_dist_quantile <-
  quantile(temp_data2$dist[temp_data2$type == "family"])
skin_family_dist_mean <-
  mean(temp_data2$dist[temp_data2$type == "family"])
skin_between_dist_quantile <-
  quantile(temp_data2$dist[temp_data2$type == "between"])
skin_between_dist_mean <-
  mean(temp_data2$dist[temp_data2$type == "between"])


oral_within_dist_quantile <-
  quantile(temp_data3$dist[temp_data3$type == "within"])
oral_within_dist_mean <-
  mean(temp_data3$dist[temp_data3$type == "within"])
oral_family_dist_quantile <-
  quantile(temp_data3$dist[temp_data3$type == "family"])
oral_family_dist_mean <-
  mean(temp_data3$dist[temp_data3$type == "family"])
oral_between_dist_quantile <-
  quantile(temp_data3$dist[temp_data3$type == "between"])
oral_between_dist_mean <-
  mean(temp_data3$dist[temp_data3$type == "between"])


nasal_within_dist_quantile <-
  quantile(temp_data4$dist[temp_data4$type == "within"])
nasal_within_dist_mean <-
  mean(temp_data4$dist[temp_data4$type == "within"])
nasal_family_dist_quantile <-
  quantile(temp_data4$dist[temp_data4$type == "family"])
nasal_family_dist_mean <-
  mean(temp_data4$dist[temp_data4$type == "family"])
nasal_between_dist_quantile <-
  quantile(temp_data4$dist[temp_data4$type == "between"])
nasal_between_dist_mean <-
  mean(temp_data4$dist[temp_data4$type == "between"])


sink(file = "median_mean_value.txt")
cat("Stool within mean\n")
stool_within_dist_mean
cat("Stool within quantile\n")
stool_within_dist_quantile
cat("Stool family mean\n")
stool_family_dist_mean
cat("Stool family quantile\n")
stool_family_dist_quantile
cat("Stool between mean\n")
stool_between_dist_mean
cat("Stool between quantile\n")
stool_between_dist_quantile

cat("skin within mean\n")
skin_within_dist_mean
cat("skin within quantile\n")
skin_within_dist_quantile
cat("skin family mean\n")
skin_family_dist_mean
cat("skin family quantile\n")
skin_family_dist_quantile
cat("skin between mean\n")
skin_between_dist_mean
cat("skin between quantile\n")
skin_between_dist_quantile

cat("oral within mean\n")
oral_within_dist_mean
cat("oral within quantile\n")
oral_within_dist_quantile
cat("oral family mean\n")
oral_family_dist_mean
cat("oral family quantile\n")
oral_family_dist_quantile
cat("oral between mean\n")
oral_between_dist_mean
cat("oral between quantile\n")
oral_between_dist_quantile

cat("nasal within mean\n")
nasal_within_dist_mean
cat("nasal within quantile\n")
nasal_within_dist_quantile
cat("nasal family mean\n")
nasal_family_dist_mean
cat("nasal family quantile\n")
nasal_family_dist_quantile
cat("nasal between mean\n")
nasal_between_dist_mean
cat("nasal between quantile\n")
nasal_between_dist_quantile

sink()


# save(temp_data1, file = "temp_data1")
# save(temp_data2, file = "temp_data2")
# save(temp_data3, file = "temp_data3")
# save(temp_data4, file = "temp_data4")

plot1 =
  temp_data1 %>%
  ggplot(aes(type, dist)) +
  geom_boxplot(aes(color = type),
               show.legend = FALSE,
               outlier.shape = NA) +
  scale_color_manual(values = between_within_color) +
  facet_wrap(facets = vars(class),
             scales = "free",
             nrow = 1) +
  base_theme +
  labs(y = "Distance", x = "")
plot1
# ggsave(plot1, filename = "stool_distance.pdf", width = 7, height = 6)



plot2 =
  temp_data2 %>%
  ggplot(aes(type, dist)) +
  geom_boxplot(aes(color = type),
               show.legend = FALSE,
               outlier.shape = NA) +
  scale_color_manual(values = between_within_color) +
  facet_wrap(facets = vars(class),
             scales = "free",
             nrow = 1) +
  base_theme +
  labs(y = "Distance", x = "")
# plot2
ggsave(plot2,
       filename = "skin_distance.pdf",
       width = 7,
       height = 6)



plot3 =
  temp_data3 %>%
  ggplot(aes(type, dist)) +
  geom_boxplot(aes(color = type),
               show.legend = FALSE,
               outlier.shape = NA) +
  scale_color_manual(values = between_within_color) +
  facet_wrap(facets = vars(class),
             scales = "free",
             nrow = 1) +
  base_theme +
  labs(y = "Distance", x = "")
# plot3
ggsave(plot3,
       filename = "oral_distance.pdf",
       width = 7,
       height = 6)


plot4 =
  temp_data4 %>%
  ggplot(aes(type, dist)) +
  geom_boxplot(aes(color = type),
               show.legend = FALSE,
               outlier.shape = NA) +
  scale_color_manual(values = between_within_color) +
  facet_wrap(facets = vars(class),
             scales = "free",
             nrow = 1) +
  base_theme +
  labs(y = "Distance", x = "")
# plot4
ggsave(plot4,
       filename = "nasal_distance.pdf",
       width = 7,
       height = 6)
