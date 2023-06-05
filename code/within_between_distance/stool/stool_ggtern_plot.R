
###
no_function()

setwd(masstools::get_project_wd())
library(tidyverse)
library(phyloseq)
rm(list = ls())

####load raw data
ls()
source(here::here("code/tools.R"))

###stool
load("data_analysis/combine_microbiome/distance/stool/personalized_score")
load("data_analysis/combine_microbiome/distance/stool/permutation_p_values")

personalized_score$fc1_p_adjust = p.adjust(personalized_score$fc1_p,
                                           method = "BH")

personalized_score$fc2_p_adjust = p.adjust(personalized_score$fc2_p,
                                           method = "BH")

personalized_score$fc3_p_adjust = p.adjust(personalized_score$fc3_p,
                                           method = "BH")

####remove the genus whose permutation test > 0.05
remain_genus = 
  permutation_p_values %>% 
  dplyr::filter(fc1_p_adjust < 0.05) %>% 
  dplyr::pull(genus)

personalized_score = 
  personalized_score %>% 
  dplyr::filter(genus %in% remain_genus)

personalized_score$fc1 = personalized_score$between_mean1 - personalized_score$within_mean1 

personalized_score$fc2 = 
  (personalized_score$between_mean1 - personalized_score$family_mean1)/(personalized_score$between_mean1 - personalized_score$within_mean1)

personalized_score$fc2[which(personalized_score$fc2_p_adjust > 0.05 |
                               personalized_score$fc3_p_adjust > 0.05)] <- NA

personalized_score$fc2[which(personalized_score$family_mean1 > personalized_score$between_mean1)] = NA

setwd(masstools::get_project_wd())
setwd("data_analysis/combine_microbiome/distance/stool/")

library(plotly)

library(scatter3D)

plot_ly(
  x = personalized_score,
  y = personalized_score$family_mean1,
  z = personalized_score$between_mean1,
  type = "scatter3d",
  mode = "markers"
)


library(plot3D)

scatter3D(x = personalized_score$within_mean1,
          y = personalized_score$family_mean1,
          z = personalized_score$between_mean1)
