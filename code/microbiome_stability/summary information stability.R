
no_function()
# set work directory

masstools::setwd_project()
library(tidyverse)
rm(list = ls())

source("code/tools.R")

######work directory
masstools::setwd_project()
setwd("data_analysis/microbiome_stability")

####load data
load("nasal_stability")
load("stool_stability")
load("skin_stability")
load("oral_stability")

####
load(here::here("data/from_xin/Revision_MultiOmes_0509.RData"))

temp = 
t(stool_stability) %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "subject_id") %>%
  dplyr::left_join(sc[,c("SubjectID", "IRIS")], by = c("subject_id" = "SubjectID"))

library(plyr)
temp2 = 
temp %>% 
plyr::dlply(.variables = .(IRIS))  %>% 
  purrr::map(function(x){
    x %>% 
      dplyr::select(-c(subject_id, IRIS)) %>% 
      colMeans(na.rm = TRUE)
  })
  

boxplot(temp2)





dim(stool_stability)
dim(skin_stability)
dim(oral_stability)
dim(nasal_stability)

rownames(stool_stability)
rownames(skin_stability)
rownames(oral_stability)
rownames(nasal_stability)

######load dmi
load(here::here("data_analysis/combine_microbiome/distance/stool/personalized_score"))
stool_personalized_score = personalized_score
stool_personalized_score$fc1 = stool_personalized_score$between_mean1 - stool_personalized_score$within_mean1

load(here::here("data_analysis/combine_microbiome/distance/skin/personalized_score"))
skin_personalized_score = personalized_score
skin_personalized_score$fc1 = skin_personalized_score$between_mean1 - skin_personalized_score$within_mean1

load(here::here("data_analysis/combine_microbiome/distance/oral/personalized_score"))
oral_personalized_score = personalized_score
oral_personalized_score$fc1 = oral_personalized_score$between_mean1 - oral_personalized_score$within_mean1

load(here::here("data_analysis/combine_microbiome/distance/nasal/personalized_score"))
nasal_personalized_score = personalized_score
nasal_personalized_score$fc1 = nasal_personalized_score$between_mean1 - nasal_personalized_score$within_mean1


stool_id = intersect(stool_personalized_score$genus, rownames(stool_stability))
temp_data_stool = 
  data.frame(stability = apply(stool_stability, 1, function(x){mean(x, na.rm = TRUE)})[stool_id],
             dmi = stool_personalized_score$fc1[match(stool_id, stool_personalized_score$genus)],
             class = "Stool")




skin_id = intersect(skin_personalized_score$genus, rownames(skin_stability))
temp_data_skin = 
  data.frame(stability = apply(skin_stability, 1, function(x){mean(x, na.rm = TRUE)})[skin_id],
             dmi = skin_personalized_score$fc1[match(skin_id, skin_personalized_score$genus)],
             class = "Skin")


nasal_id = intersect(nasal_personalized_score$genus, rownames(nasal_stability))
temp_data_nasal = 
  data.frame(stability = apply(nasal_stability, 1, function(x){mean(x, na.rm = TRUE)})[nasal_id],
             dmi = nasal_personalized_score$fc1[match(nasal_id, nasal_personalized_score$genus)],
             class = "Nasal")


oral_id = intersect(oral_personalized_score$genus, rownames(oral_stability))
temp_data_oral = 
  data.frame(stability = apply(oral_stability, 1, function(x){mean(x, na.rm = TRUE)})[oral_id],
             dmi = oral_personalized_score$fc1[match(oral_id, oral_personalized_score$genus)],
             class = "Oral")


temp_data =
  rbind(temp_data_stool,
        temp_data_skin,
        temp_data_oral,
        temp_data_nasal)

plot = 
temp_data %>% 
  ggplot(aes(dmi, stability)) +
  geom_point(aes(color = class), show.legend = FALSE) +
  geom_smooth(method = "lm", color = "black") +
  facet_wrap(facets = vars(class), scales = "free") +
  base_theme +
  scale_color_manual(values = body_site_color) +
  labs(x = "DMI", y = "Stability")

plot
ggsave(plot, filename = "stability_plot.pdf", width = 8, height = 7)

####add phylum to temp_data
load(here::here("data_analysis/stool_microbiome/data_preparation/variable_info"))
stool_variable_info =
  variable_info

load(here::here("data_analysis/skin_microbiome/data_preparation/variable_info"))
skin_variable_info =
  variable_info

load(here::here("data_analysis/oral_microbiome/data_preparation/variable_info"))
oral_variable_info =
  variable_info

load(here::here("data_analysis/nasal_microbiome/data_preparation/variable_info"))
nasal_variable_info =
  variable_info


library(plyr)


temp_data = 
temp_data %>%
  tibble::rownames_to_column(var = "Genus") %>%
  dplyr::mutate(Genus = stringr::str_replace(Genus, "[1,2]{1}", "")) %>% 
  plyr::dlply(.variables = .(class)) %>%
  purrr::map(function(x) {
    if(x$class[1] == "Nasal"){
      x = 
      x %>% 
        dplyr::left_join(nasal_variable_info[,c("Phylum", "Genus")] %>% 
                           dplyr::distinct(.keep_all = TRUE), by = "Genus") %>% 
        dplyr::filter(!is.na(Phylum))
    }
    
    if(x$class[1] == "Stool"){
      x = 
        x %>% 
        dplyr::left_join(stool_variable_info[,c("Phylum", "Genus")] %>% 
                           dplyr::distinct(.keep_all = TRUE), 
                         by = "Genus") %>% 
        dplyr::filter(!is.na(Phylum))
    }
    
    if(x$class[1] == "Skin"){
      x = 
        x %>% 
        dplyr::left_join(skin_variable_info[,c("Phylum", "Genus")] %>% 
                           dplyr::distinct(.keep_all = TRUE), 
                         by = "Genus") %>% 
        dplyr::filter(!is.na(Phylum))
    }
    
    if(x$class[1] == "Oral"){
      x = 
        x %>% 
        dplyr::left_join(oral_variable_info[,c("Phylum", "Genus")] %>% 
                           dplyr::distinct(.keep_all = TRUE), 
                         by = "Genus") %>% 
        dplyr::filter(!is.na(Phylum))
    }
    
    x
    
  }) %>% 
  dplyr::bind_rows()


######stool stability plot
###only remain the phylum at least 5 points
temp_data_stool = 
temp_data %>% 
  dplyr::filter(class == "Stool")

remain_phylum = 
temp_data_stool %>% 
  dplyr::count(Phylum)  %>% 
  dplyr::filter(n >= 5) %>% 
  pull(Phylum)

stool_stability_plot = 
temp_data_stool %>% 
  dplyr::filter(Phylum %in% remain_phylum) %>%
  ggplot(aes(dmi, stability)) +
  geom_point(aes(color = Phylum), show.legend = FALSE) +
  geom_smooth(method = "lm", color = "black") +
  facet_wrap(facets = vars(Phylum), scales = "free") +
  base_theme +
  scale_color_manual(values = phylum_color) +
  labs(x = "DMI", y = "Stability")
stool_stability_plot
ggsave(stool_stability_plot,
       filename = "stool_stability_plot.pdf", 
       width = 8, height = 7)




######skin stability plot
###only remain the phylum at least 5 points
temp_data_skin = 
  temp_data %>% 
  dplyr::filter(class == "Skin")

remain_phylum = 
  temp_data_skin %>% 
  dplyr::count(Phylum)  %>% 
  dplyr::filter(n >= 5) %>% 
  pull(Phylum)

skin_stability_plot = 
  temp_data_skin %>% 
  dplyr::filter(Phylum %in% remain_phylum) %>%
  ggplot(aes(dmi, stability)) +
  geom_point(aes(color = Phylum), show.legend = FALSE) +
  geom_smooth(method = "lm", color = "black") +
  facet_wrap(facets = vars(Phylum), scales = "free") +
  base_theme +
  scale_color_manual(values = phylum_color) +
  labs(x = "DMI", y = "Stability")
skin_stability_plot
ggsave(skin_stability_plot,
       filename = "skin_stability_plot.pdf", 
       width = 8, height = 7)


######nasal stability plot
###only remain the phylum at least 5 points
temp_data_nasal = 
  temp_data %>% 
  dplyr::filter(class == "Nasal")

remain_phylum = 
  temp_data_nasal %>% 
  dplyr::count(Phylum)  %>% 
  dplyr::filter(n >= 5) %>% 
  pull(Phylum)

nasal_stability_plot = 
  temp_data_nasal %>% 
  dplyr::filter(Phylum %in% remain_phylum) %>%
  ggplot(aes(dmi, stability)) +
  geom_point(aes(color = Phylum), show.legend = FALSE) +
  geom_smooth(method = "lm", color = "black") +
  facet_wrap(facets = vars(Phylum), scales = "free") +
  base_theme +
  scale_color_manual(values = phylum_color) +
  labs(x = "DMI", y = "Stability")
nasal_stability_plot
ggsave(nasal_stability_plot,
       filename = "nasal_stability_plot.pdf", 
       width = 12, height = 4)




######oral stability plot
###only remain the phylum at least 5 points
temp_data_oral = 
  temp_data %>% 
  dplyr::filter(class == "Oral")

remain_phylum = 
  temp_data_oral %>% 
  dplyr::count(Phylum)  %>% 
  dplyr::filter(n >= 5) %>% 
  pull(Phylum)

oral_stability_plot = 
  temp_data_oral %>% 
  dplyr::filter(Phylum %in% remain_phylum) %>%
  ggplot(aes(dmi, stability)) +
  geom_point(aes(color = Phylum), show.legend = FALSE) +
  geom_smooth(method = "lm", color = "black") +
  facet_wrap(facets = vars(Phylum), scales = "free") +
  base_theme +
  scale_color_manual(values = phylum_color) +
  labs(x = "DMI", y = "Stability")
oral_stability_plot
ggsave(oral_stability_plot,
       filename = "oral_stability_plot.pdf", 
       width = 4, height = 4)









