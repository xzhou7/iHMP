
###
no_function()

masstools::setwd_project()
library(tidyverse)
library(phyloseq)
rm(list = ls())

####load raw data
load("data/from_xin/DetailedPhyloseq.RData")
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

####remove the genus whose permutation test > 0.05, if the permutation p value > 0
####means that the between and within distance are same, so the personization score
####is 0, so here we just remove the genus whose personization score is 0
####remove Unidentified_ASV
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

# personalized_score$fc2[which(personalized_score$fc2_p_adjust > 0.05 &
#                                personalized_score$fc3_p_adjust > 0.05)] <- NA

personalized_score$fc2[which(personalized_score$family_mean1 > personalized_score$between_mean1)] = 0

stool_personalized_score = personalized_score

#####skin
load("data_analysis/combine_microbiome/distance/skin/personalized_score")
load("data_analysis/combine_microbiome/distance/skin/permutation_p_values")

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

# personalized_score$fc2[which(personalized_score$fc2_p_adjust > 0.2 &
#                                personalized_score$fc3_p_adjust > 0.2)] <- NA

personalized_score$fc2[which(personalized_score$family_mean1 > personalized_score$between_mean1)] = 0

skin_personalized_score = personalized_score

#####oral
load("data_analysis/combine_microbiome/distance/oral/personalized_score")
load("data_analysis/combine_microbiome/distance/oral/permutation_p_values")

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

# personalized_score$fc2[which(personalized_score$fc2_p_adjust > 0.05 &
#                                personalized_score$fc3_p_adjust > 0.05)] <- NA

personalized_score$fc2[which(personalized_score$family_mean1 > personalized_score$between_mean1)] = NA

oral_personalized_score = personalized_score

#####nasal
load("data_analysis/combine_microbiome/distance/nasal/personalized_score")
load("data_analysis/combine_microbiome/distance/nasal/permutation_p_values")

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

# personalized_score$fc2[which(personalized_score$fc2_p_adjust > 0.05 &
#                                personalized_score$fc3_p_adjust > 0.05)] <- NA

personalized_score$fc2[which(personalized_score$family_mean1 > personalized_score$between_mean1)] = NA

nasal_personalized_score = personalized_score

dir.create("data_analysis/combine_microbiome")
dir.create("data_analysis/combine_microbiome/cladogram/")
setwd("data_analysis/combine_microbiome/cladogram/")

temp_data = 
  rbind(stool_personalized_score,
        skin_personalized_score,
        oral_personalized_score,
        nasal_personalized_score) %>% 
  dplyr::rename(class = dataset) %>% 
  dplyr::mutate(class = stringr::str_to_sentence(class))

personalized_score = temp_data

library(ggvenn)
library(ggVennDiagram)

plot = 
ggVennDiagram(x = 
  list(
    Stool = unique(stool_personalized_score$genus),
    Skin = unique(stool_personalized_score$genus),
    Oral = unique(oral_personalized_score$genus),
    Nasal = unique(nasal_personalized_score$genus)
  ), label_color = "white", 
  label_geom = "text"
) +
  scale_colour_manual(values = body_site_color)
plot
# ggsave(plot, filename = "4_body_site_ovrelap.pdf", width = 9, height = 7)

####HMP
physeqGenus_ST
physeqGenus_SK
physeqGenus_NS
physeqGenus_OR

physeqGenus_ST =
  subset_taxa(physeq = physeqGenus_ST, Genus %in% temp_data$genus)
physeqGenus_SK =
  subset_taxa(physeq = physeqGenus_SK, Genus %in% temp_data$genus)
physeqGenus_OR =
  subset_taxa(physeq = physeqGenus_OR, Genus %in% temp_data$genus)
physeqGenus_NS =
  subset_taxa(physeq = physeqGenus_NS, Genus %in% temp_data$genus)

#####combine them together as a new phyloseq class
expression_data_stool =
  physeqGenus_ST@otu_table@.Data %>%
  as.data.frame()

variable_info_stool = as.data.frame(physeqGenus_ST@tax_table@.Data)
sample_info_stool = get_variable(physeq = physeqGenus_ST)

expression_data_skin =
  physeqGenus_SK@otu_table@.Data %>%
  as.data.frame()
variable_info_skin = as.data.frame(physeqGenus_SK@tax_table@.Data)
sample_info_skin = get_variable(physeq = physeqGenus_SK)

expression_data_nasal =
  physeqGenus_NS@otu_table@.Data %>%
  as.data.frame()
variable_info_nasal = as.data.frame(physeqGenus_NS@tax_table@.Data)
sample_info_nasal = get_variable(physeq = physeqGenus_NS)

expression_data_oral =
  physeqGenus_OR@otu_table@.Data %>%
  as.data.frame()
variable_info_oral = as.data.frame(physeqGenus_OR@tax_table@.Data)
sample_info_oral = get_variable(physeq = physeqGenus_OR)

dim(sample_info_stool)

variable_info_stool

variable_info = 
  dplyr::full_join(variable_info_stool, 
                   variable_info_skin, 
                   by = colnames(variable_info_stool)) %>% 
  dplyr::full_join(variable_info_nasal, 
                   by = colnames(variable_info_stool)) %>% 
  dplyr::full_join(variable_info_oral, 
                   by = colnames(variable_info_stool)) %>% 
  dplyr::distinct(Genus, .keep_all = TRUE)

rownames(variable_info) = 
  paste("asv", 1:nrow(variable_info), sep = "_")

expression_data = 
  matrix(1:(nrow(variable_info) * 2), ncol = 2)

colnames(expression_data) = paste("sample", 1:ncol(expression_data), sep = "_")
rownames(expression_data) = rownames(variable_info)

sample_info = 
  data.frame(sample_id = colnames(expression_data))
rownames(sample_info) = sample_info$sample_id

##remove Archaea
remove_name = 
  c("Zea",
    variable_info %>% 
      dplyr::filter(Kingdom == "Archaea") %>% 
      dplyr::pull(Genus)
  ) %>% 
  unique()

if(length(remove_name) > 0) {
  variable_info = 
    variable_info %>% 
    dplyr::filter(!Genus %in% remove_name)
  expression_data = 
    expression_data[rownames(variable_info),]
}

expression_data = otu_table(expression_data, taxa_are_rows = TRUE)
variable_info = tax_table(as.matrix(variable_info))
sample_info = sample_data(sample_info)

physeqGenus = phyloseq(expression_data, variable_info, sample_info)

physeqGenus2 =
  subset_taxa(physeq = physeqGenus, Genus %in% temp_data$genus)

# ##remove some genus which have no root
# ##
# physeqGenus2 =
#   subset_taxa(physeq = physeqGenus2, !Genus %in% "Unclassified_Actinobacteria")

library(microbiomeViz)
library(ggtreeExtra)
library(ggtree)
library(phyloseq)
library(dplyr)

GP <- physeqGenus2
GP
GP = fix_duplicate_tax(GP)
GP

###create tree
tree = 
  microbiomeViz::parsePhyloseq(
  physeq = GP,
  use_abundance = FALSE,
  node.size.scale = 0,
  node.size.offset = 0
)

raw_p <-
  tree.backbone(
  tree = tree,
  size = 0.3,
  shape = 16,
  layout = "circular",
  fill = "black",
  color = "black"
)

raw_p

raw_p$data$nodeClass2 = 
  as.character(raw_p$data$nodeClass)

raw_p$data$nodeClass2[is.na(raw_p$data$nodeClass2)] = "Root"
raw_p$data$nodeClass2[raw_p$data$nodeClass2 == "f"] = "Family"
raw_p$data$nodeClass2[raw_p$data$nodeClass2 == "c"] = "Class"
raw_p$data$nodeClass2[raw_p$data$nodeClass2 == "o"] = "Order"
raw_p$data$nodeClass2[raw_p$data$nodeClass2 == "p"] = "Phylum"
raw_p$data$nodeClass2[raw_p$data$nodeClass2 == "g"] = "Genus"
raw_p$data$nodeClass2[raw_p$data$nodeClass2 == "k"] = "Kingdom"

raw_p$data$nodeSize2 = raw_p$data$nodeSize
raw_p$data$nodeSize2[raw_p$data$nodeClass2 == "Root"] = 3.5
raw_p$data$nodeSize2[raw_p$data$nodeClass2 == "Kingdom"] = 3
raw_p$data$nodeSize2[raw_p$data$nodeClass2 == "Phylum"] =2.5
raw_p$data$nodeSize2[raw_p$data$nodeClass2 == "Class"] = 2
raw_p$data$nodeSize2[raw_p$data$nodeClass2 == "Order"] = 1.5
raw_p$data$nodeSize2[raw_p$data$nodeClass2 == "Family"] = 1
raw_p$data$nodeSize2[raw_p$data$nodeClass2 == "Genus"] = 0.5

# #####add node point
raw_p =
raw_p +
  geom_point2(aes(color = nodeClass2,
                  size = nodeSize2),
             show.legend = FALSE) +
  ggsci::scale_color_jama() +
  ggnewscale::new_scale(new_aes = "fill")

raw_p

###label 2 is the name of taxa
raw_p$data$label2 = 
  raw_p$data$label %>% 
  stringr::str_split("__") %>% 
  purrr::map(function(x){
    x[2]
  }) %>% 
  unlist()

raw_p$data$label2[as.character(raw_p$data$nodeClass) != "g"] = NA

#####add new information
temp_data_fc1 =
  temp_data %>% 
  dplyr::select(genus, class, fc1) %>% 
  tidyr::pivot_wider(names_from = class, values_from = "fc1")

colnames(temp_data_fc1)[-1] = paste(colnames(temp_data_fc1)[-1], "fc1", sep = "_")

temp_data_fc2 =
  temp_data %>% 
  dplyr::select(genus, class, fc2) %>% 
  tidyr::pivot_wider(names_from = class, values_from = "fc2")

colnames(temp_data_fc2)[-1] = paste(colnames(temp_data_fc2)[-1], "fc2", sep = "_")

temp_data_fc1_p_adjust =
  temp_data %>% 
  dplyr::select(genus, class, fc1_p_adjust) %>% 
  tidyr::pivot_wider(names_from = class, values_from = "fc1_p_adjust")

colnames(temp_data_fc1_p_adjust)[-1] = 
  paste(colnames(temp_data_fc1_p_adjust)[-1], "fc1_p_adjust", sep = "_")
  
new_info =
  data.frame(Genus = raw_p$data$label2) %>%
  dplyr::left_join(as.data.frame(variable_info)[, c("Genus", "Kingdom", "Phylum", "Class", "Order", "Family")],
                   by = "Genus") %>%
  dplyr::left_join(temp_data_fc1, by = c("Genus" = "genus")) %>% 
  dplyr::left_join(temp_data_fc1_p_adjust, by = c("Genus" = "genus")) %>% 
  dplyr::left_join(temp_data_fc2, by = c("Genus" = "genus")) %>% 
  dplyr::mutate(Stool = case_when(
    !is.na(Stool_fc1) ~ "Stool",
    TRUE ~ "no"
  )) %>% 
  dplyr::mutate(Skin = case_when(
    !is.na(Skin_fc1) ~ "Skin",
    TRUE ~ "no"
  )) %>% 
  dplyr::mutate(Oral = case_when(
    !is.na(Oral_fc1) ~ "Oral",
    TRUE ~ "no"
  )) %>% 
  dplyr::mutate(Nasal = case_when(
    !is.na(Nasal_fc1) ~ "Nasal",
    TRUE ~ "no"
  )) %>% 
  dplyr::mutate(Stool_fc1_star = case_when(
    !is.na(Stool_fc1_p_adjust) & Stool_fc1_p_adjust < 0.05 ~ "*",
    TRUE ~ ""
  )) %>% 
  dplyr::mutate(Skin_fc1_star = case_when(
    !is.na(Skin_fc1_p_adjust) & Skin_fc1_p_adjust < 0.05 ~ "*",
    TRUE ~ ""
  )) %>% 
  dplyr::mutate(Oral_fc1_star = case_when(
    !is.na(Oral_fc1_p_adjust) & Oral_fc1_p_adjust < 0.05 ~ "*",
    TRUE ~ ""
  )) %>% 
  dplyr::mutate(Nasal_fc1_star = case_when(
    !is.na(Nasal_fc1_p_adjust) & Nasal_fc1_p_adjust < 0.05 ~ "*",
    TRUE ~ ""
  ))

raw_p$data =
  cbind(raw_p$data, new_info)

###add tip label (genus label)
# #####add point to show this body is here or not
# raw_p = 
# raw_p +
#   geom_tippoint(
#     aes(fill = Stool),
#     shape = 21,
#     x = 6.15,
#     size = 1,
#     show.legend = FALSE
#   ) +
#   scale_fill_manual(values = c(body_site_color, "no" = "white")) +
#   ggnewscale::new_scale(new_aes = "fill") +
#   geom_tippoint(
#     aes(fill = Skin),
#     shape = 21,
#     x = 6.35,
#     size = 1,
#     show.legend = FALSE
#   ) +
#   scale_fill_manual(values = c(body_site_color, "no" = "white")) +
#   ggnewscale::new_scale(new_aes = "fill") +
#   geom_tippoint(
#     aes(fill = Oral),
#     shape = 21,
#     x = 6.55,
#     size = 1,
#     show.legend = FALSE
#   ) +
#   scale_fill_manual(values = c(body_site_color, "no" = "white")) +
#   ggnewscale::new_scale(new_aes = "fill") +
#   geom_tippoint(
#     aes(fill = Nasal),
#     shape = 21,
#     x = 6.75,
#     size = 1,
#     show.legend = FALSE
#   ) +
#   scale_fill_manual(values = c(body_site_color, "no" = "white"))
# 
#  
# raw_p


######add highlight
hight_data = 
  raw_p$data %>% 
  dplyr::filter(stringr::str_detect(label, "p__")) %>% 
  dplyr::select(node, label) %>% 
  dplyr::mutate(label = stringr::str_replace_all(label, "p__", "")) %>% 
  dplyr::rename(id = node, type = label)

p1 = 
raw_p +
  geom_hilight(
    data = hight_data,
    mapping = aes(node = id,
                  fill = type),
    alpha = .4
  ) +
  scale_fill_manual(values = phylum_color)

p1

######add heatmap
##add fc1 information
p1$data$Stool_fc1[which(p1$data$Stool_fc1 < 0)] = 0
p1$data$Skin_fc1[which(p1$data$Skin_fc1 < 0)] = 0
p1$data$Oral_fc1[which(p1$data$Oral_fc1 < 0)] = 0
p1$data$Nasal_fc1[which(p1$data$Nasal_fc1 < 0)] = 0

##add multiple tip information
stool_fc1_info =
  p1$data[, c("Stool_fc1", "isTip")]
rownames(stool_fc1_info) = p1$data$label
stool_fc1_info = stool_fc1_info %>%
  dplyr::filter(isTip) %>%
  dplyr::select(-isTip)

skin_fc1_info =
  p1$data[, c("Skin_fc1", "isTip")]
rownames(skin_fc1_info) = p1$data$label
skin_fc1_info = skin_fc1_info %>%
  dplyr::filter(isTip) %>%
  dplyr::select(-isTip)

oral_fc1_info =
  p1$data[, c("Oral_fc1", "isTip")]
rownames(oral_fc1_info) = p1$data$label
oral_fc1_info = oral_fc1_info %>%
  dplyr::filter(isTip) %>%
  dplyr::select(-isTip)

nasal_fc1_info =
  p1$data[, c("Nasal_fc1", "isTip")]
rownames(nasal_fc1_info) = p1$data$label
nasal_fc1_info = nasal_fc1_info %>%
  dplyr::filter(isTip) %>%
  dplyr::select(-isTip)

range(
  c(
    p1$data$Stool_fc1,
    p1$data$Skin_fc1,
    p1$data$Nasal_fc1,
    p1$data$Oral_fc1
  ),
  na.rm = TRUE
)

p1 = p1 +
  ggnewscale::new_scale_fill() 

p2_1 =
  gheatmap(
    p = p1,
    data = stool_fc1_info,
    offset = -0.1,
    width = .08,
    colnames_angle = 95,
    colnames_offset_y = .25, 
    colnames = FALSE,
    color = alpha(body_site_color["Stool"], 1),
    legend_title = "Index1"
  ) +
  scale_fill_gradientn(colours = c(RColorBrewer::brewer.pal(n = 11, name = "BrBG"))[-c(5:7)], 
                       na.value = "white", limits = c(0,0.8)) +
  ggnewscale::new_scale(new_aes = "fill") +
    geom_tippoint(
      mapping = aes(shape = Stool_fc1_star),
      x = 6.35,
      size = 2,
      show.legend = FALSE
    ) +
    scale_shape_manual(values = c("*" = "*")) +
  ggnewscale::new_scale(new_aes = "shape")
  
p2_1

p2_2 =   
gheatmap(p = p2_1,
    data = skin_fc1_info,
    offset = 0.5,
    width = .08,
    colnames_angle = 95,
    colnames_offset_y = .25, 
    colnames = FALSE,
    color = alpha(body_site_color["Skin"], 1),
    legend_title = "Index1"
  ) +
  scale_fill_gradientn(colours = c(RColorBrewer::brewer.pal(n = 11, name = "BrBG"))[-c(5:7)], 
                       na.value = "white",
                       limits = c(0,0.8)) +
  ggnewscale::new_scale(new_aes = "fill") +
  geom_tippoint(
    mapping = aes(shape = Skin_fc1_star),
    x = 6.95,
    size = 2,
    show.legend = FALSE
  ) +
  scale_shape_manual(values = c("*" = "*")) +
  ggnewscale::new_scale(new_aes = "shape")
  
p2_2

p2_3 =   
  gheatmap(p = p2_2,
           data = oral_fc1_info,
           offset = 1.1,
           width = .08,
           colnames_angle = 95,
           colnames_offset_y = .25, 
           colnames = FALSE,
           color = alpha(body_site_color["Oral"], 1),
           legend_title = "Index1"
  ) +
  scale_fill_gradientn(colours = c(RColorBrewer::brewer.pal(n = 11, name = "BrBG"))[-c(5:7)], 
                       na.value = "white",
                       limits = c(0,0.8)) +
  ggnewscale::new_scale(new_aes = "fill") +
  geom_tippoint(
    mapping = aes(shape = Oral_fc1_star),
    x = 7.55,
    size = 2,
    show.legend = FALSE
  ) +
  scale_shape_manual(values = c("*" = "*")) +
  ggnewscale::new_scale(new_aes = "shape")

p2_3

p2_4 =   
  gheatmap(p = p2_3,
           data = nasal_fc1_info,
           offset = 1.7,
           width = .08,
           colnames_angle = 95,
           colnames_offset_y = .25, 
           colnames = FALSE,
           color = alpha(body_site_color["Nasal"], 1),
           legend_title = "Index1"
  ) +
  scale_fill_gradientn(colours = c(RColorBrewer::brewer.pal(n = 11, name = "BrBG"))[-c(5:7)], 
                       na.value = "white",
                       limits = c(0,0.8)) +
  ggnewscale::new_scale(new_aes = "fill") +
  geom_tippoint(
    mapping = aes(shape = Nasal_fc1_star),
    x = 8.15,
    size = 2,
    show.legend = FALSE
  ) +
  scale_shape_manual(values = c("*" = "*")) +
  ggnewscale::new_scale(new_aes = "shape")

p2_4

##add tip lab
##only add some tip points
# idx1 = which(!is.na(raw_p$data$Stool_fc1_p_adjust) & raw_p$data$Stool_fc1_p_adjust < 0.001)
# idx1 = sample(idx1, 10)
# raw_p$data$label2[-idx1] = NA
p3 =
  p2_4 +
  ggnewscale::new_scale(new_aes = "color") +
  geom_tiplab(
    aes(label = label2,
        color = Phylum),
    offset = 2.5,
    size = 2,
    show.legend = FALSE
  )

p3

# ggsave(p3, filename = "fc1_cladogram.pdf", width = 14, height = 14)

###output results
temp_data =
  p3$data

temp_data =
  temp_data %>%
  dplyr::filter(isTip)
  
# openxlsx::write.xlsx(temp_data,
#                      file = "fc1_tree_data.xlsx",
#                      asTable = TRUE,
#                      overwrite = TRUE)

#####summary information
plot <- 
personalized_score %>%
  dplyr::mutate(significant = case_when(
    fc1_p_adjust < 0.05 ~ class,
    TRUE ~ "no"
  )) %>%
  dplyr::mutate(class = factor(class, levels = c("Stool", "Skin", "Oral", "Nasal"))) %>%
  ggplot(aes(class, fc1)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(size = -log(fc1_p_adjust, 10),
                  fill = significant),
              alpha = 0.9,
              shape = 21) +
  scale_fill_manual(values = c(body_site_color, "no" = "white")) +
  guides(size = guide_legend(title = "-log10(p-adjust)")) +
  base_theme +
  labs(x = "", y = "Personized score")

plot

# ggsave(plot, filename = "Personized_score_boxplot.pdf", width = 10, height = 7)

temp <-
  personalized_score %>%
  dplyr::mutate(significant = case_when(
    fc1_p_adjust < 0.05 ~ class,
    TRUE ~ "no"
  )) %>%
  dplyr::mutate(class = factor(class, levels = c("Stool", "Skin", "Oral", "Nasal")))

stool_vs_skin_test <-
  wilcox.test(temp$fc1[temp$class == "Stool"],
              temp$fc1[temp$class == "Skin"])

stool_vs_oral_test <-
  wilcox.test(temp$fc1[temp$class == "Stool"],
              temp$fc1[temp$class == "Oral"])

stool_vs_nasal_test <-
  wilcox.test(temp$fc1[temp$class == "Stool"],
              temp$fc1[temp$class == "Nasal"])

skin_vs_oral_test <-
  wilcox.test(temp$fc1[temp$class == "Skin"],
              temp$fc1[temp$class == "Oral"])

skin_vs_nasal_test <-
  wilcox.test(temp$fc1[temp$class == "Skin"],
              temp$fc1[temp$class == "Nasal"])

oral_vs_nasal_test <-
  wilcox.test(temp$fc1[temp$class == "Oral"],
              temp$fc1[temp$class == "Nasal"])


sink(file = "test_result.txt")
cat("stool vs skin")
stool_vs_skin_test
cat("stool vs oral")
stool_vs_oral_test
cat("stool vs nasal")
stool_vs_nasal_test
cat("skin vs oral")
skin_vs_oral_test
cat("skin vs nasal")
skin_vs_nasal_test
cat("oral vs nasal")
oral_vs_nasal_test
sink()


stool_mean <- 
  mean(temp$fc1[temp$class == "Stool"])
stool_quantile <- 
  quantile(temp$fc1[temp$class == "Stool"])

skin_mean <- 
  mean(temp$fc1[temp$class == "Skin"])
skin_quantile <- 
  quantile(temp$fc1[temp$class == "Skin"])

oral_mean <- 
  mean(temp$fc1[temp$class == "Oral"])
oral_quantile <- 
  quantile(temp$fc1[temp$class == "Oral"])

nasal_mean <- 
  mean(temp$fc1[temp$class == "Nasal"])
nasal_quantile <- 
  quantile(temp$fc1[temp$class == "Nasal"])


sink(file = "median_mean_value.txt")
cat("stool mean\n")
stool_mean
cat("stool quantile\n")
stool_quantile

cat("skin mean\n")
skin_mean
cat("skin quantile\n")
skin_quantile

cat("oral mean\n")
oral_mean
cat("oral quantile\n")
oral_quantile

cat("nasal mean\n")
nasal_mean
cat("nasal quantile\n")
nasal_quantile

sink()





library(plyr)

variable_info =
  rbind(variable_info_stool,
        variable_info_skin,
        variable_info_oral,
        variable_info_nasal) %>% 
  dplyr::select(Phylum, Genus) %>% 
  dplyr::distinct(.keep_all = TRUE)

remain_phylum = 
personalized_score %>%
  dplyr::left_join(variable_info, by = c("genus" = "Genus")) %>%
  dplyr::group_by(Phylum, class) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::filter(n >= 5) %>% 
  dplyr::group_by(Phylum) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::filter(n > 1) %>% 
  pull(Phylum)

remain_phylum = remain_phylum[remain_phylum != "Cyanobacteria/Chloroplast"]

plot = 
  personalized_score %>%
  dplyr::mutate(significant = case_when(
    fc1_p_adjust < 0.05 ~ class,
    TRUE ~ "no"
  )) %>%
  dplyr::left_join(variable_info, by = c("genus" = "Genus")) %>% 
  dplyr::mutate(class = factor(class, levels = c("Stool", "Skin", "Oral", "Nasal"))) %>%
  dplyr::mutate(fc1 = round(fc1, 2)) %>% 
  dplyr::filter(Phylum %in% remain_phylum) %>% 
  ggplot(aes(class, fc1)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(size = -log(fc1_p_adjust, 10),
                  fill = significant),
              alpha = 0.9,
              shape = 21) +
  scale_fill_manual(values = c(body_site_color, "no" = "white")) +
  guides(size = guide_legend(title = "-log10(p-adjust)")) +
  base_theme +
  scale_size_continuous(range = c(0.01, 1)) +
  labs(x = "", y = "Personized score") +
  facet_wrap(facets = vars(Phylum), scales = "free_y")

plot

# ggsave(plot,
#        filename = "Personized_score_boxplot_based_on_phylum.pdf",
#        width = 16,
#        height = 10)


temp_data <- 
personalized_score %>%
  dplyr::mutate(significant = case_when(
    fc1_p_adjust < 0.05 ~ class,
    TRUE ~ "no"
  )) %>%
  dplyr::left_join(variable_info, by = c("genus" = "Genus")) %>% 
  dplyr::mutate(class = factor(class, levels = c("Stool", "Skin", "Oral", "Nasal"))) %>%
  dplyr::mutate(fc1 = round(fc1, 2)) %>% 
  # dplyr::filter(Phylum %in% remain_phylum) %>% 
  dplyr::mutate(Phylum = case_when(
    Phylum == "Actinobacteria" | Phylum == "Bacteroidetes" | Phylum == "Firmicutes" |  Phylum == "Proteobacteria" ~ Phylum,
    TRUE ~ "Other"
  )) %>% 
  dplyr::mutate(Phylum = factor(Phylum, levels = c("Actinobacteria", 
                                                   "Bacteroidetes",
                                                   "Firmicutes",
                                                   "Proteobacteria",
                                                   "Other")))


plot = 
  temp_data %>% 
  ggplot(aes(Phylum, fc1)) +
  geom_jitter(aes(fill = Phylum),
              size = 3.5,
              shape = 21,
              alpha = 0.7,
              show.legend = FALSE) +
  geom_boxplot(outlier.shape = NA, fill = "transparent") +
  scale_fill_manual(values = c(phylum_color, "Other" = "grey")) +
  # guides(size = guide_legend(title = "-log10(p-adjust)")) +
  base_theme +
  # scale_size_continuous(range = c(0.05, 2)) +
  labs(x = "", y = "Personized score") +
  facet_wrap(facets = vars(class),
             # scales = "free_y",
             nrow = 1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
        strip.text = element_text(size = 10))

plot

# ggsave(plot,
#        filename = "Personized_score_boxplot_based_on_body_site.pdf",
#        width = 10,
#        height = 4)


temp_data <- 
  temp_data %>% 
  dplyr::select(-c(fc1_p, fc2_p, fc3_p, fc1_p_adjust, fc2_p_adjust, fc3_p_adjust, fc2)) %>% 
  dplyr::rename(dmi = fc1)

write.csv(temp_data, "dmi_based_on_phylum.csv", row.names = FALSE)





