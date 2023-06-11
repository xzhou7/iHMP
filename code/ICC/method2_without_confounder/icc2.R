#' ---
#' title: "nasal microbiome metabolome correlation"
#' author: 
#'   - name: "Xiaotao Shen" 
#'     url: https://www.shenxt.info/
#'     affiliation: Stanford School of Medicine
#' date: "`r Sys.Date()`"
#' site: distill::distill_website
#' output: 
#'   distill::distill_article:
#'     code_folding: false
#' ---

#+ r setup, echo=TRUE, eval = TRUE, include = TRUE

no_function()
# set work directory

setwd(masstools::get_project_wd())
library(tidyverse)
rm(list = ls())
source("code/tools.R")

###ASV
####load data asv level
load(here::here("data_analysis/nasal_microbiome/ICC_without_confounder/vd_nasal_asv2"))
load(here::here("data_analysis/oral_microbiome/ICC_without_confounder/vd_oral_asv2"))
load(here::here("data_analysis/skin_microbiome/ICC_without_confounder/vd_skin_asv2"))
load(here::here("data_analysis/stool_microbiome/ICC_without_confounder/vd_stool_asv2"))

#Combine all vd tables
vd.comb = 
  rbind(
    data.frame(vd_nasal_asv2, dt = "Nasal"),
    data.frame(vd_oral_asv2, dt = "Oral"),
    data.frame(vd_skin_asv2, dt = "Skin"),
    data.frame(vd_stool_asv2, dt = "Stool")    
  ) %>% 
  dplyr::filter(!is.na(random_Subject) & !is.na(random_Residual))

vd.comb$ICC = vd.comb$random_Subject*100/(vd.comb$random_Subject+vd.comb$random_Residual)
vd.comb = 
vd.comb %>% 
  dplyr::filter(!is.na(ICC))

library(reshape2)
pd = vd.comb[, 4:5]
pd$Variables = rownames(pd)
pd <- melt(pd)

pd = 
  pd %>% 
  dplyr::filter(!is.na(value)) %>% 
    dplyr::mutate(dt = factor(dt, levels = c("Skin", "Oral", "Nasal", "Stool")))

library(gghalves)

pd %>% 
  dplyr::group_by(dt) %>% 
  dplyr::summarise(mean = mean(value),
                   median = median(value),
                   sd = sd(value))

library(ggsignif)

plot = 
ggplot(data = pd, aes(x = dt, y = value)) +
  geom_jitter(aes(color = dt), show.legend = FALSE, alpha = 1, size = 4) +
  geom_boxplot(color = "black",
               fill = "transparent",
               show.legend = FALSE) +
  # geom_half_violin(aes(fill = dt), show.legend = FALSE) +
  # scale_y_continuous(limits = c(0,1)) +
  # geom_half_boxplot(aes(fill = dt), outlier.shape = NA, side = "r") + 
  # geom_dotplot(binaxis = "y", method="histodot", stackdir="up",
  #              aes(fill = dt), dotsize = 0.5,
  #              , show.legend = FALSE) +
  scale_color_manual(values = body_site_color) +
  scale_fill_manual(values = body_site_color) +
  base_theme +
  labs(x = "", y = "ICC of Normalized Variables") +
  geom_signif(
    comparisons = list(c("Skin", "Oral"),
                       c("Skin", "Nasal"),
                       c("Skin", "Stool"),
                       c("Oral", "Nasal"),
                       c("Oral", "Stool"),
                       c("Nasal", "Stool")),
    step_increase = 0.05,
    test = "wilcox.test",
    map_signif_level = TRUE
  )
plot

setwd(masstools::get_project_wd())
dir.create("data_analysis/ICC_without_confounder")
setwd("data_analysis/ICC_without_confounder")
plot
# ggsave(plot, filename = "icc_of_different_body_site_microbiome_asv2.pdf", width = 9, height = 7)

pd_asv = pd

# plot = 
#   ggplot(data=pd, aes(x = dt, y = value)) +
#   geom_boxplot(aes(fill = dt), outlier.shape = NA, show.legend = FALSE) + 
#   scale_color_manual(values = body_site_color) +
#   scale_fill_manual(values = body_site_color) +
#   base_theme +
#   labs(x = "", y = "ICC of Normalized Variables")
# plot
# ggsave(plot, filename = "icc_of_different_body_site_microbiome2_asv.pdf", width = 9, height = 7)




###Class
####load data class level
load(here::here("data_analysis/nasal_microbiome/ICC_without_confounder/vd_nasal_class2"))
load(here::here("data_analysis/oral_microbiome/ICC_without_confounder/vd_oral_class2"))
load(here::here("data_analysis/skin_microbiome/ICC_without_confounder/vd_skin_class2"))
load(here::here("data_analysis/stool_microbiome/ICC_without_confounder/vd_stool_class2"))

#Combine all vd tables
vd.comb = 
  rbind(
    data.frame(vd_nasal_class2, dt = "Nasal"),
    data.frame(vd_oral_class2, dt = "Oral"),
    data.frame(vd_skin_class2, dt = "Skin"),
    data.frame(vd_stool_class2, dt = "Stool")    
  ) %>% 
  dplyr::filter(!is.na(random_Subject) & !is.na(random_Residual))

vd.comb$ICC = vd.comb$random_Subject*100/(vd.comb$random_Subject+vd.comb$random_Residual)
vd.comb = 
  vd.comb %>% 
  dplyr::filter(!is.na(ICC))

library(reshape2)
pd = vd.comb[, 4:5]
pd$Variables = rownames(pd)
pd <- melt(pd)

pd = 
  pd %>% 
  dplyr::filter(!is.na(value)) %>% 
  dplyr::mutate(dt = factor(dt, levels = c("Skin", "Oral", "Nasal", "Stool")))

library(gghalves)

pd %>% 
  dplyr::group_by(dt) %>% 
  dplyr::summarise(mean = mean(value),
                   median = median(value),
                   sd = sd(value))


plot = 
  ggplot(data = pd, aes(x = dt, y = value)) +
  geom_jitter(aes(color = dt), show.legend = FALSE, alpha = 1,
              size = 4) +
  geom_boxplot(color = "black",
               fill = "transparent",
               show.legend = FALSE, outlier.shape = NA) +
  # geom_half_violin(aes(fill = dt), show.legend = FALSE) +
  # scale_y_continuous(limits = c(0,1)) +
  # geom_half_boxplot(aes(fill = dt), outlier.shape = NA, side = "r") + 
  # geom_dotplot(binaxis = "y", method="histodot", stackdir="up",
  #              aes(fill = dt), dotsize = 0.5,
  #              , show.legend = FALSE) +
  scale_color_manual(values = body_site_color) +
  scale_fill_manual(values = body_site_color) +
  base_theme +
  labs(x = "", y = "ICC of Normalized Variables") +
  geom_signif(
    comparisons = list(c("Skin", "Oral"),
                       c("Skin", "Nasal"),
                       c("Skin", "Stool"),
                       c("Oral", "Nasal"),
                       c("Oral", "Stool"),
                       c("Nasal", "Stool")),
    step_increase = 0.05,
    test = "wilcox.test",
    map_signif_level = TRUE
  )
plot

# ggsave(plot, filename = "icc_of_different_body_site_microbiome_class2.pdf", width = 9, height = 7)

pd_class = pd

# plot = 
#   ggplot(data=pd, aes(x = dt, y = value)) +
#   geom_boxplot(aes(fill = dt), outlier.shape = NA, show.legend = FALSE) + 
#   scale_color_manual(values = body_site_color) +
#   scale_fill_manual(values = body_site_color) +
#   base_theme +
#   labs(x = "", y = "ICC of Normalized Variables")
# plot
# ggsave(plot, filename = "icc_of_different_body_site_microbiome2_class.pdf", width = 9, height = 7)






###family
####load data family level
load(here::here("data_analysis/nasal_microbiome/ICC_without_confounder/vd_nasal_family2"))
load(here::here("data_analysis/oral_microbiome/ICC_without_confounder/vd_oral_family2"))
load(here::here("data_analysis/skin_microbiome/ICC_without_confounder/vd_skin_family2"))
load(here::here("data_analysis/stool_microbiome/ICC_without_confounder/vd_stool_family2"))

#Combine all vd tables
vd.comb = 
  rbind(
    data.frame(vd_nasal_family2, dt = "Nasal"),
    data.frame(vd_oral_family2, dt = "Oral"),
    data.frame(vd_skin_family2, dt = "Skin"),
    data.frame(vd_stool_family2, dt = "Stool")    
  ) %>% 
  dplyr::filter(!is.na(random_Subject) & !is.na(random_Residual))

vd.comb$ICC = vd.comb$random_Subject*100/(vd.comb$random_Subject+vd.comb$random_Residual)
vd.comb = 
  vd.comb %>% 
  dplyr::filter(!is.na(ICC))

library(reshape2)
pd = vd.comb[, 4:5]
pd$Variables = rownames(pd)
pd <- melt(pd)

pd = 
  pd %>% 
  dplyr::filter(!is.na(value)) %>% 
  dplyr::mutate(dt = factor(dt, levels = c("Skin", "Oral", "Nasal", "Stool")))

library(gghalves)

pd %>% 
  dplyr::group_by(dt) %>% 
  dplyr::summarise(mean = mean(value),
                   median = median(value),
                   sd = sd(value))


plot = 
  ggplot(data = pd, aes(x = dt, y = value)) +
  geom_jitter(aes(color = dt), show.legend = FALSE, alpha = 1,
              size = 4) +
  geom_boxplot(color = "black",
               fill = "transparent",
               show.legend = FALSE, outlier.shape = NA) +
  # geom_half_violin(aes(fill = dt), show.legend = FALSE) +
  # scale_y_continuous(limits = c(0,1)) +
  # geom_half_boxplot(aes(fill = dt), outlier.shape = NA, side = "r") + 
  # geom_dotplot(binaxis = "y", method="histodot", stackdir="up",
  #              aes(fill = dt), dotsize = 0.5,
  #              , show.legend = FALSE) +
  scale_color_manual(values = body_site_color) +
  scale_fill_manual(values = body_site_color) +
  base_theme +
  labs(x = "", y = "ICC of Normalized Variables") +
  geom_signif(
    comparisons = list(c("Skin", "Oral"),
                       c("Skin", "Nasal"),
                       c("Skin", "Stool"),
                       c("Oral", "Nasal"),
                       c("Oral", "Stool"),
                       c("Nasal", "Stool")),
    step_increase = 0.05,
    test = "wilcox.test",
    map_signif_level = TRUE
  )
plot

# ggsave(plot, filename = "icc_of_different_body_site_microbiome_family2.pdf", width = 9, height = 7)

pd_family = pd

# plot = 
#   ggplot(data=pd, aes(x = dt, y = value)) +
#   geom_boxplot(aes(fill = dt), outlier.shape = NA, show.legend = FALSE) + 
#   scale_color_manual(values = body_site_color) +
#   scale_fill_manual(values = body_site_color) +
#   base_theme +
#   labs(x = "", y = "ICC of Normalized Variables")
# plot
# ggsave(plot, filename = "icc_of_different_body_site_microbiome2_family.pdf", width = 9, height = 7)



###genus
####load data genus level
load(here::here("data_analysis/nasal_microbiome/ICC_without_confounder/vd_nasal_genus2"))
load(here::here("data_analysis/oral_microbiome/ICC_without_confounder/vd_oral_genus2"))
load(here::here("data_analysis/skin_microbiome/ICC_without_confounder/vd_skin_genus2"))
load(here::here("data_analysis/stool_microbiome/ICC_without_confounder/vd_stool_genus2"))

#Combine all vd tables
vd.comb = 
  rbind(
    data.frame(vd_nasal_genus2, dt = "Nasal"),
    data.frame(vd_oral_genus2, dt = "Oral"),
    data.frame(vd_skin_genus2, dt = "Skin"),
    data.frame(vd_stool_genus2, dt = "Stool")    
  ) %>% 
  dplyr::filter(!is.na(random_Subject) & !is.na(random_Residual))

vd.comb$ICC = vd.comb$random_Subject*100/(vd.comb$random_Subject+vd.comb$random_Residual)
vd.comb = 
  vd.comb %>% 
  dplyr::filter(!is.na(ICC))

library(reshape2)
pd = vd.comb[, 4:5]
pd$Variables = rownames(pd)
pd <- melt(pd)

pd = 
  pd %>% 
  dplyr::filter(!is.na(value)) %>% 
  dplyr::mutate(dt = factor(dt, levels = c("Skin", "Oral", "Nasal", "Stool")))

library(gghalves)

pd %>% 
  dplyr::group_by(dt) %>% 
  dplyr::summarise(mean = mean(value),
                   median = median(value),
                   sd = sd(value))


plot = 
  ggplot(data = pd, aes(x = dt, y = value)) +
  geom_jitter(aes(color = dt), show.legend = FALSE, alpha = 1,
              size = 4) +
  geom_boxplot(color = "black",
               fill = "transparent",
               show.legend = FALSE, outlier.shape = NA) +
  # geom_half_violin(aes(fill = dt), show.legend = FALSE) +
  # scale_y_continuous(limits = c(0,1)) +
  # geom_half_boxplot(aes(fill = dt), outlier.shape = NA, side = "r") + 
  # geom_dotplot(binaxis = "y", method="histodot", stackdir="up",
  #              aes(fill = dt), dotsize = 0.5,
  #              , show.legend = FALSE) +
  scale_color_manual(values = body_site_color) +
  scale_fill_manual(values = body_site_color) +
  base_theme +
  labs(x = "", y = "ICC of Normalized Variables") +
  geom_signif(
    comparisons = list(c("Skin", "Oral"),
                       c("Skin", "Nasal"),
                       c("Skin", "Stool"),
                       c("Oral", "Nasal"),
                       c("Oral", "Stool"),
                       c("Nasal", "Stool")),
    step_increase = 0.05,
    test = "wilcox.test",
    map_signif_level = TRUE
  )
plot

# ggsave(plot, filename = "icc_of_different_body_site_microbiome_genus2.pdf", width = 9, height = 7)

pd_genus = pd

# plot = 
#   ggplot(data=pd, aes(x = dt, y = value)) +
#   geom_boxplot(aes(fill = dt), outlier.shape = NA, show.legend = FALSE) + 
#   scale_color_manual(values = body_site_color) +
#   scale_fill_manual(values = body_site_color) +
#   base_theme +
#   labs(x = "", y = "ICC of Normalized Variables")
# plot
# ggsave(plot, filename = "icc_of_different_body_site_microbiome2_genus.pdf", width = 9, height = 7)



###order
####load data order level
load(here::here("data_analysis/nasal_microbiome/ICC_without_confounder/vd_nasal_order2"))
load(here::here("data_analysis/oral_microbiome/ICC_without_confounder/vd_oral_order2"))
load(here::here("data_analysis/skin_microbiome/ICC_without_confounder/vd_skin_order2"))
load(here::here("data_analysis/stool_microbiome/ICC_without_confounder/vd_stool_order2"))

#Combine all vd tables
vd.comb = 
  rbind(
    data.frame(vd_nasal_order2, dt = "Nasal"),
    data.frame(vd_oral_order2, dt = "Oral"),
    data.frame(vd_skin_order2, dt = "Skin"),
    data.frame(vd_stool_order2, dt = "Stool")    
  ) %>% 
  dplyr::filter(!is.na(random_Subject) & !is.na(random_Residual))

vd.comb$ICC = vd.comb$random_Subject*100/(vd.comb$random_Subject+vd.comb$random_Residual)
vd.comb = 
  vd.comb %>% 
  dplyr::filter(!is.na(ICC))

library(reshape2)
pd = vd.comb[, 4:5]
pd$Variables = rownames(pd)
pd <- melt(pd)

pd = 
  pd %>% 
  dplyr::filter(!is.na(value)) %>% 
  dplyr::mutate(dt = factor(dt, levels = c("Skin", "Oral", "Nasal", "Stool")))

library(gghalves)

pd %>% 
  dplyr::group_by(dt) %>% 
  dplyr::summarise(mean = mean(value),
                   median = median(value),
                   sd = sd(value))


plot = 
  ggplot(data = pd, aes(x = dt, y = value)) +
  geom_jitter(aes(color = dt), show.legend = FALSE, alpha = 1,
              size = 4) +
  geom_boxplot(color = "black",
               fill = "transparent",
               show.legend = FALSE, outlier.shape = NA) +
  # geom_half_violin(aes(fill = dt), show.legend = FALSE) +
  # scale_y_continuous(limits = c(0,1)) +
  # geom_half_boxplot(aes(fill = dt), outlier.shape = NA, side = "r") + 
  # geom_dotplot(binaxis = "y", method="histodot", stackdir="up",
  #              aes(fill = dt), dotsize = 0.5,
  #              , show.legend = FALSE) +
  scale_color_manual(values = body_site_color) +
  scale_fill_manual(values = body_site_color) +
  base_theme +
  labs(x = "", y = "ICC of Normalized Variables") +
  geom_signif(
    comparisons = list(c("Skin", "Oral"),
                       c("Skin", "Nasal"),
                       c("Skin", "Stool"),
                       c("Oral", "Nasal"),
                       c("Oral", "Stool"),
                       c("Nasal", "Stool")),
    step_increase = 0.05,
    test = "wilcox.test",
    map_signif_level = TRUE
  )
plot

# ggsave(plot, filename = "icc_of_different_body_site_microbiome_order2.pdf", width = 9, height = 7)

pd_order = pd

# plot = 
#   ggplot(data=pd, aes(x = dt, y = value)) +
#   geom_boxplot(aes(fill = dt), outlier.shape = NA, show.legend = FALSE) + 
#   scale_color_manual(values = body_site_color) +
#   scale_fill_manual(values = body_site_color) +
#   base_theme +
#   labs(x = "", y = "ICC of Normalized Variables")
# plot
# ggsave(plot, filename = "icc_of_different_body_site_microbiome2_order.pdf", width = 9, height = 7)



###phylum
####load data phylum level
load(here::here("data_analysis/nasal_microbiome/ICC_without_confounder/vd_nasal_phylum2"))
load(here::here("data_analysis/oral_microbiome/ICC_without_confounder/vd_oral_phylum2"))
load(here::here("data_analysis/skin_microbiome/ICC_without_confounder/vd_skin_phylum2"))
load(here::here("data_analysis/stool_microbiome/ICC_without_confounder/vd_stool_phylum2"))

#Combine all vd tables
vd.comb = 
  rbind(
    data.frame(vd_nasal_phylum2, dt = "Nasal"),
    data.frame(vd_oral_phylum2, dt = "Oral"),
    data.frame(vd_skin_phylum2, dt = "Skin"),
    data.frame(vd_stool_phylum2, dt = "Stool")    
  ) %>% 
  dplyr::filter(!is.na(random_Subject) & !is.na(random_Residual))

vd.comb$ICC = vd.comb$random_Subject*100/(vd.comb$random_Subject+vd.comb$random_Residual)
vd.comb = 
  vd.comb %>% 
  dplyr::filter(!is.na(ICC))

library(reshape2)
pd = vd.comb[, 4:5]
pd$Variables = rownames(pd)
pd <- melt(pd)

pd = 
  pd %>% 
  dplyr::filter(!is.na(value)) %>% 
  dplyr::mutate(dt = factor(dt, levels = c("Skin", "Oral", "Nasal", "Stool")))

library(gghalves)

pd %>% 
  dplyr::group_by(dt) %>% 
  dplyr::summarise(mean = mean(value),
                   median = median(value),
                   sd = sd(value))


plot = 
  ggplot(data = pd, aes(x = dt, y = value)) +
  geom_jitter(aes(color = dt), show.legend = FALSE, alpha = 1,
              size = 4) +
  geom_boxplot(color = "black",
               fill = "transparent",
               show.legend = FALSE, outlier.shape = NA) +
  # geom_half_violin(aes(fill = dt), show.legend = FALSE) +
  # scale_y_continuous(limits = c(0,1)) +
  # geom_half_boxplot(aes(fill = dt), outlier.shape = NA, side = "r") + 
  # geom_dotplot(binaxis = "y", method="histodot", stackdir="up",
  #              aes(fill = dt), dotsize = 0.5,
  #              , show.legend = FALSE) +
  scale_color_manual(values = body_site_color) +
  scale_fill_manual(values = body_site_color) +
  base_theme +
  labs(x = "", y = "ICC of Normalized Variables") +
  geom_signif(
    comparisons = list(c("Skin", "Oral"),
                       c("Skin", "Nasal"),
                       c("Skin", "Stool"),
                       c("Oral", "Nasal"),
                       c("Oral", "Stool"),
                       c("Nasal", "Stool")),
    step_increase = 0.05,
    test = "wilcox.test",
    map_signif_level = TRUE
  )
plot

# ggsave(plot, filename = "icc_of_different_body_site_microbiome_phylum2.pdf", width = 9, height = 7)

pd_phylum = pd

# plot = 
  # ggplot(data=pd, aes(x = dt, y = value)) +
  # geom_boxplot(aes(fill = dt), outlier.shape = NA, show.legend = FALSE) +
  # scale_color_manual(values = body_site_color) +
  # scale_fill_manual(values = body_site_color) +
  # base_theme +
  # labs(x = "", y = "ICC of Normalized Variables")
# plot
# ggsave(plot, filename = "icc_of_different_body_site_microbiome2_phylum.pdf", width = 9, height = 7)



###combine them together
temp_data =
  rbind(
    data.frame(pd_phylum, level = "Phylum"),
    data.frame(pd_class, level = "Class"),
    data.frame(pd_order, level = "Order"),
    data.frame(pd_family, level = "Family"),
    data.frame(pd_genus, level = "Genus"),
    data.frame(pd_asv, level = "ASV")
  ) %>% 
  dplyr::mutate(level = factor(level,
                               levels = c("Phylum", "Class", "Order", "Family", "Genus", "ASV")))


head(temp_data)

ggplot(data = temp_data,
       x = level, 
       y = value,
       color = dt) +
  geom_point(aes(x = level, 
                 y = value,
                 fill = dt),
             alpha = 1,
             shape = 16,
             size = 1,
             position=position_jitterdodge(dodge.width = 0.75),
             show.legend = TRUE) +
  geom_boxplot(
    aes(x = level, 
        y = value,
        fill = dt),
    position = position_dodge(width = 0.75),
    alpha = 0.7,
    show.legend = TRUE,
    outlier.shape = NA
  ) +
  stat_summary(fun = median,
               geom = "line",
               aes(
                 x = level,
                 y = value,
                 group = dt,
                 color = dt
               ),
               position = position_dodge(width = 0.75)) +
  scale_color_manual(values = body_site_color[-5]) +
  scale_fill_manual(values = body_site_color[-5]) +
  base_theme +
  labs(x = "", y = "ICC of Normalized Variables")


library(ggpubr)
# Box plots
temp_data$dt <- 
  factor(x = temp_data$dt, levels = c("Stool", "Skin", "Oral", "Nasal"))

save(temp_data, file = "temp_data")

bxp <- ggboxplot(
  temp_data,
  x = "level",
  y = "value",
  fill = "dt",
  outlier.shape = NA,
  width = 0.75
) +
  geom_point(
    aes(x = level,
        y = value,
        fill = dt),
    alpha = 0.7,
    shape = 16,
    size = 1,
    position = position_jitterdodge(dodge.width = 0.75),
    show.legend = TRUE
  ) +
  stat_summary(
    fun = median,
    geom = "line",
    aes(
      x = level,
      y = value,
      group = dt,
      color = dt
    ),
    position = position_dodge(width = 0.75)
  ) +
  scale_color_manual(values = body_site_color[-5]) +
  scale_fill_manual(values = body_site_color[-5]) +
  base_theme +
  labs(x = "", y = "ICC of Normalized Variables")
bxp

ggsave(bxp, filename = "icc_plot.pdf", width = 12, height = 7)

library(rstatix)

stat.test <- temp_data %>%
  dplyr::group_by(level) %>%
  rstatix::wilcox_test(value ~ dt, p.adjust.method = "BH") %>% 
  dplyr::filter(p.adj < 0.05)

stat.test 

stat.test <- stat.test %>%
  add_xy_position(x = "level", dodge = 0.75)

bxp2 = 
bxp +
  stat_pvalue_manual(
    data = stat.test,
    label = "p.adj.signif",
    tip.length = 0.01,
    bracket.nudge.y = -2
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))

bxp2

pwc <- 
  temp_data %>%
  group_by(dt) %>%
  rstatix::wilcox_test(value ~ level, p.adjust.method = "BH") %>% 
  # dplyr::filter(p.adj < 0.05) %>% 
  dplyr::mutate(name = paste(group1, group2, sep = "_")) %>% 
  dplyr::filter(name %in% c("Genus_ASV", 
                            "Family_Genus", 
                            "Order_Family",
                            "Class_Order", 
                            "Phylum_Class",
                            "ASV_Genus", 
                            "Genus_Family", 
                            "Family_Order",
                            "Order_Class", 
                            "Class_Phylum"
                            )) %>% 
  dplyr::select(-name)

pwc <- pwc %>% add_xy_position(x = "level")

bxp3 = 
bxp2 +
  stat_pvalue_manual(
    pwc,
    color = "dt",
    step.group.by = "dt",
    tip.length = 0.1,
    step.increase = 0,
    label = "p"
  )

bxp3
bxp2

# ggsave(bxp3, filename = "icc_plot2.pdf", width = 12, height = 7)

temp_data

icc_data <- temp_data
# save(icc_data, file = "icc_data")

# 
# #######only for nasal
# ###nasal variable information
# load(here::here("data_analysis/nasal_microbiome/data_preparation/variable_info"))
# 
# nasal_microbiome_variable_info = variable_info
# 
# nasal_microbiome_variable_info
# 
# temp_data_nasal =
#   temp_data %>%
#   dplyr::filter(dt == "Nasal")
# 
# segment_data1 = 
#   nasal_microbiome_variable_info[,c("Genus", "variable_id")] %>% 
#   dplyr::rename(ASV = variable_id) %>% 
#   dplyr::left_join(temp_data_nasal %>% 
#                      dplyr::filter(level == "Genus") %>% 
#                      dplyr::select(Variables, value), by = c("Genus" = "Variables")) %>% 
#   dplyr::rename(value1 = value) %>% 
#   dplyr::left_join(temp_data_nasal %>% 
#                      dplyr::filter(level == "ASV") %>% 
#                      dplyr::select(Variables, value), by = c("ASV" = "Variables")) %>% 
#   dplyr::rename(value2 = value) %>% 
#   dplyr::filter(!is.na(value1) & !is.na(value2)) %>% 
#   dplyr::distinct(.keep_all = TRUE)
# 
# 
# segment_data2 =
#   nasal_microbiome_variable_info[, c("Family", "Genus")] %>%
#   dplyr::left_join(
#     temp_data_nasal %>%
#       dplyr::filter(level == "Family") %>%
#       dplyr::select(Variables, value),
#     by = c("Family" = "Variables")
#   ) %>%
#   dplyr::rename(value1 = value) %>%
#   dplyr::left_join(
#     temp_data_nasal %>%
#       dplyr::filter(level == "Genus") %>%
#       dplyr::select(Variables, value),
#     by = c("Genus" = "Variables")
#   ) %>%
#   dplyr::rename(value2 = value) %>%
#   dplyr::filter(!is.na(value1) & !is.na(value2)) %>% 
#   dplyr::distinct(.keep_all = TRUE)
# 
# 
# segment_data3 =
#   nasal_microbiome_variable_info[, c("Order", "Family")] %>%
#   dplyr::left_join(
#     temp_data_nasal %>%
#       dplyr::filter(level == "Order") %>%
#       dplyr::select(Variables, value),
#     by = c("Order" = "Variables")
#   ) %>%
#   dplyr::rename(value1 = value) %>%
#   dplyr::left_join(
#     temp_data_nasal %>%
#       dplyr::filter(level == "Family") %>%
#       dplyr::select(Variables, value),
#     by = c("Family" = "Variables")
#   ) %>%
#   dplyr::rename(value2 = value) %>%
#   dplyr::filter(!is.na(value1) & !is.na(value2)) %>% 
#   dplyr::distinct(.keep_all = TRUE)
# 
# 
# segment_data4 =
#   nasal_microbiome_variable_info[, c("Class", "Order")] %>%
#   dplyr::left_join(
#     temp_data_nasal %>%
#       dplyr::filter(level == "Class") %>%
#       dplyr::select(Variables, value),
#     by = c("Class" = "Variables")
#   ) %>%
#   dplyr::rename(value1 = value) %>%
#   dplyr::left_join(
#     temp_data_nasal %>%
#       dplyr::filter(level == "Order") %>%
#       dplyr::select(Variables, value),
#     by = c("Order" = "Variables")
#   ) %>%
#   dplyr::rename(value2 = value) %>%
#   dplyr::filter(!is.na(value1) & !is.na(value2)) %>% 
#   dplyr::distinct(.keep_all = TRUE)
# 
# 
# segment_data5 =
#   nasal_microbiome_variable_info[, c("Phylum", "Class")] %>%
#   dplyr::left_join(
#     temp_data_nasal %>%
#       dplyr::filter(level == "Phylum") %>%
#       dplyr::select(Variables, value),
#     by = c("Phylum" = "Variables")
#   ) %>%
#   dplyr::rename(value1 = value) %>%
#   dplyr::left_join(
#     temp_data_nasal %>%
#       dplyr::filter(level == "Class") %>%
#       dplyr::select(Variables, value),
#     by = c("Class" = "Variables")
#   ) %>%
#   dplyr::rename(value2 = value) %>%
#   dplyr::filter(!is.na(value1) & !is.na(value2)) %>% 
#   dplyr::distinct(.keep_all = TRUE)
#   
# 
# temp_data_nasal %>%
#   ggplot(x = level,
#          y = value,
#          color = dt) +
#   geom_boxplot(
#     aes(x = level,
#         y = value,
#         color = dt),
#     alpha = 0.7,
#     show.legend = FALSE,
#     outlier.shape = NA
#   ) +
#   geom_segment(data = segment_data1,
#                aes(
#                  x = "Genus",
#                  y = value1,
#                  xend = "ASV",
#                  yend = value2
#                ),
#                color = "grey") +
#   geom_segment(data = segment_data2,
#                aes(
#                  x = "Family",
#                  y = value1,
#                  xend = "Genus",
#                  yend = value2
#                ),
#                color = "grey") +
#   geom_segment(data = segment_data3,
#                aes(
#                  x = "Order",
#                  y = value1,
#                  xend = "Family",
#                  yend = value2
#                ),
#                color = "grey") +
#   geom_segment(data = segment_data4,
#                aes(
#                  x = "Class",
#                  y = value1,
#                  xend = "Order",
#                  yend = value2
#                ),
#                color = "grey") +
#   geom_segment(data = segment_data5,
#                aes(
#                  x = "Phylum",
#                  y = value1,
#                  xend = "Class",
#                  yend = value2
#                ),
#                color = "grey") +
#   geom_point(
#     aes(x = level,
#         y = value),
#     alpha = 1,
#     shape = 16,
#     size = 2,
#     show.legend = FALSE
#   ) +
#   scale_color_manual(values = body_site_color[-5]) +
#   base_theme +
#   labs(x = "", y = "ICC of Normalized Variables")
# 
# 
# 
# 
# 
# 
# 
# 
# #######only for skin
# ###skin variable information
# load(here::here("data_analysis/skin_microbiome/data_preparation/variable_info"))
# 
# skin_microbiome_variable_info = variable_info
# 
# skin_microbiome_variable_info
# 
# temp_data_skin =
#   temp_data %>%
#   dplyr::filter(dt == "Skin")
# 
# segment_data1 = 
#   skin_microbiome_variable_info[,c("Genus", "variable_id")] %>% 
#   dplyr::rename(ASV = variable_id) %>% 
#   dplyr::left_join(temp_data_skin %>% 
#                      dplyr::filter(level == "Genus") %>% 
#                      dplyr::select(Variables, value), by = c("Genus" = "Variables")) %>% 
#   dplyr::rename(value1 = value) %>% 
#   dplyr::left_join(temp_data_skin %>% 
#                      dplyr::filter(level == "ASV") %>% 
#                      dplyr::select(Variables, value), by = c("ASV" = "Variables")) %>% 
#   dplyr::rename(value2 = value) %>% 
#   dplyr::filter(!is.na(value1) & !is.na(value2)) %>% 
#   dplyr::distinct(.keep_all = TRUE)
# 
# 
# segment_data2 =
#   skin_microbiome_variable_info[, c("Family", "Genus")] %>%
#   dplyr::left_join(
#     temp_data_skin %>%
#       dplyr::filter(level == "Family") %>%
#       dplyr::select(Variables, value),
#     by = c("Family" = "Variables")
#   ) %>%
#   dplyr::rename(value1 = value) %>%
#   dplyr::left_join(
#     temp_data_skin %>%
#       dplyr::filter(level == "Genus") %>%
#       dplyr::select(Variables, value),
#     by = c("Genus" = "Variables")
#   ) %>%
#   dplyr::rename(value2 = value) %>%
#   dplyr::filter(!is.na(value1) & !is.na(value2)) %>% 
#   dplyr::distinct(.keep_all = TRUE)
# 
# 
# segment_data3 =
#   skin_microbiome_variable_info[, c("Order", "Family")] %>%
#   dplyr::left_join(
#     temp_data_skin %>%
#       dplyr::filter(level == "Order") %>%
#       dplyr::select(Variables, value),
#     by = c("Order" = "Variables")
#   ) %>%
#   dplyr::rename(value1 = value) %>%
#   dplyr::left_join(
#     temp_data_skin %>%
#       dplyr::filter(level == "Family") %>%
#       dplyr::select(Variables, value),
#     by = c("Family" = "Variables")
#   ) %>%
#   dplyr::rename(value2 = value) %>%
#   dplyr::filter(!is.na(value1) & !is.na(value2)) %>% 
#   dplyr::distinct(.keep_all = TRUE)
# 
# 
# segment_data4 =
#   skin_microbiome_variable_info[, c("Class", "Order")] %>%
#   dplyr::left_join(
#     temp_data_skin %>%
#       dplyr::filter(level == "Class") %>%
#       dplyr::select(Variables, value),
#     by = c("Class" = "Variables")
#   ) %>%
#   dplyr::rename(value1 = value) %>%
#   dplyr::left_join(
#     temp_data_skin %>%
#       dplyr::filter(level == "Order") %>%
#       dplyr::select(Variables, value),
#     by = c("Order" = "Variables")
#   ) %>%
#   dplyr::rename(value2 = value) %>%
#   dplyr::filter(!is.na(value1) & !is.na(value2)) %>% 
#   dplyr::distinct(.keep_all = TRUE)
# 
# 
# segment_data5 =
#   skin_microbiome_variable_info[, c("Phylum", "Class")] %>%
#   dplyr::distinct(.keep_all = TRUE) %>% 
#   dplyr::left_join(
#     temp_data_skin %>%
#       dplyr::filter(level == "Phylum") %>%
#       dplyr::select(Variables, value),
#     by = c("Phylum" = "Variables")
#   ) %>%
#   dplyr::rename(value1 = value) %>%
#   dplyr::left_join(
#     temp_data_skin %>%
#       dplyr::filter(level == "Class") %>%
#       dplyr::select(Variables, value),
#     by = c("Class" = "Variables")
#   ) %>%
#   dplyr::rename(value2 = value) %>%
#   dplyr::filter(!is.na(value1) & !is.na(value2)) %>% 
#   dplyr::distinct(.keep_all = TRUE)
# 
# 
# temp_data_skin %>%
#   ggplot(x = level,
#          y = value,
#          color = dt) +
#   geom_boxplot(
#     aes(x = level,
#         y = value,
#         color = dt),
#     alpha = 0.7,
#     show.legend = FALSE,
#     outlier.shape = NA
#   ) +
#   geom_segment(data = segment_data1,
#                aes(
#                  x = "Genus",
#                  y = value1,
#                  xend = "ASV",
#                  yend = value2
#                ),
#                color = "grey") +
#   geom_segment(data = segment_data2,
#                aes(
#                  x = "Family",
#                  y = value1,
#                  xend = "Genus",
#                  yend = value2
#                ),
#                color = "grey") +
#   geom_segment(data = segment_data3,
#                aes(
#                  x = "Order",
#                  y = value1,
#                  xend = "Family",
#                  yend = value2
#                ),
#                color = "grey") +
#   geom_segment(data = segment_data4,
#                aes(
#                  x = "Class",
#                  y = value1,
#                  xend = "Order",
#                  yend = value2
#                ),
#                color = "grey") +
#   geom_segment(data = segment_data5,
#                aes(
#                  x = "Phylum",
#                  y = value1,
#                  xend = "Class",
#                  yend = value2
#                ),
#                color = "grey") +
#   geom_point(
#     aes(x = level,
#         y = value),
#     alpha = 1,
#     shape = 16,
#     size = 2,
#     show.legend = FALSE
#   ) +
#   scale_color_manual(values = body_site_color[-5]) +
#   base_theme +
#   labs(x = "", y = "ICC of Normalized Variables")
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
