###no source
setwd(masstools::get_project_wd())
rm(list = ls())
library(tidyverse)

source("code/tools.R")

#####nasal
data = readr::read_csv("Figures/Figure3/strainreplacement/nasal.replacement.csv")

load(here::here("data_analysis/combine_microbiome/distance/nasal/nasal_braydist_by_genus"))
load(here::here("data_analysis/combine_microbiome/distance/nasal/nasal_jaccard_by_genus"))

nasal_braydist_by_genus = 
nasal_braydist_by_genus %>% 
  dplyr::filter(type == "within") %>% 
  dplyr::select("sample_id1", "sample_id2", "dist", "genus")

nasal_braydist_by_genus$name = 
  apply(nasal_braydist_by_genus, 1, function(x){
    x = as.character(x)
    paste(sort(x[1:2]), collapse = "_")
  })

data$name = 
  apply(data, 1, function(x){
    x = as.character(x)
    paste(sort(x[2:3]), collapse = "_")
  })

match(data$name, nasal_braydist_by_genus$name)

data =
  data %>%
  dplyr::left_join(nasal_braydist_by_genus, by = c("name", "taxa" = "genus"))


data = 
data  %>% 
  dplyr::left_join(nasal_jaccard_by_genus, by = c("Start", "End", "taxa", "SubjectID"))

readr::write_csv(x = data, file = "Figures/Figure3/strainreplacement/nasal.replacement_new.csv")







###no source
setwd(masstools::get_project_wd())
rm(list = ls())
library(tidyverse)

source("code/tools.R")

#####oral
data = readr::read_csv("Figures/Figure3/strainreplacement/oral.replacement.csv")

load(here::here("data_analysis/combine_microbiome/distance/oral/oral_braydist_by_genus"))
load(here::here("data_analysis/combine_microbiome/distance/oral/oral_jaccard_by_genus"))

oral_braydist_by_genus = 
  oral_braydist_by_genus %>% 
  dplyr::filter(type == "within") %>% 
  dplyr::select("sample_id1", "sample_id2", "dist", "genus")

oral_braydist_by_genus$name = 
  apply(oral_braydist_by_genus, 1, function(x){
    x = as.character(x)
    paste(sort(x[1:2]), collapse = "_")
  })

data$name = 
  apply(data, 1, function(x){
    x = as.character(x)
    paste(sort(x[2:3]), collapse = "_")
  })

match(data$name, oral_braydist_by_genus$name)

data =
  data %>%
  dplyr::left_join(oral_braydist_by_genus, by = c("name", "taxa" = "genus"))


data = 
  data  %>% 
  dplyr::left_join(oral_jaccard_by_genus, by = c("Start", "End", "taxa", "SubjectID"))

readr::write_csv(x = data, file = "Figures/Figure3/strainreplacement/oral.replacement_new.csv")








###no source
setwd(masstools::get_project_wd())
rm(list = ls())
library(tidyverse)

source("code/tools.R")

#####stool
data = readr::read_csv("Figures/Figure3/strainreplacement/stool.replacement.csv")

load(here::here("data_analysis/combine_microbiome/distance/stool/stool_braydist_by_genus"))
load(here::here("data_analysis/combine_microbiome/distance/stool/stool_jaccard_by_genus"))

stool_braydist_by_genus = 
  stool_braydist_by_genus %>% 
  dplyr::filter(type == "within") %>% 
  dplyr::select("sample_id1", "sample_id2", "dist", "genus")

stool_braydist_by_genus$name = 
  apply(stool_braydist_by_genus, 1, function(x){
    x = as.character(x)
    paste(sort(x[1:2]), collapse = "_")
  })

data$name = 
  apply(data, 1, function(x){
    x = as.character(x)
    paste(sort(x[2:3]), collapse = "_")
  })

match(data$name, stool_braydist_by_genus$name)

data =
  data %>%
  dplyr::left_join(stool_braydist_by_genus, by = c("name", "taxa" = "genus"))


data = 
  data  %>% 
  dplyr::left_join(stool_jaccard_by_genus, by = c("Start", "End", "taxa", "SubjectID"))

readr::write_csv(x = data, file = "Figures/Figure3/strainreplacement/stool.replacement_new.csv")






###no source
setwd(masstools::get_project_wd())
rm(list = ls())
library(tidyverse)

source("code/tools.R")

#####skin
data = readr::read_csv("Figures/Figure3/strainreplacement/skin.replacement.csv")

load(here::here("data_analysis/combine_microbiome/distance/skin/skin_braydist_by_genus"))
load(here::here("data_analysis/combine_microbiome/distance/skin/skin_jaccard_by_genus"))

skin_braydist_by_genus = 
  skin_braydist_by_genus %>% 
  dplyr::filter(type == "within") %>% 
  dplyr::select("sample_id1", "sample_id2", "dist", "genus")

skin_braydist_by_genus$name = 
  apply(skin_braydist_by_genus, 1, function(x){
    x = as.character(x)
    paste(sort(x[1:2]), collapse = "_")
  })

data$name = 
  apply(data, 1, function(x){
    x = as.character(x)
    paste(sort(x[2:3]), collapse = "_")
  })

match(data$name, skin_braydist_by_genus$name)

data =
  data %>%
  dplyr::left_join(skin_braydist_by_genus, by = c("name", "taxa" = "genus"))


data = 
  data  %>% 
  dplyr::left_join(skin_jaccard_by_genus, by = c("Start", "End", "taxa", "SubjectID"))

readr::write_csv(x = data, file = "Figures/Figure3/strainreplacement/skin.replacement_new.csv")

