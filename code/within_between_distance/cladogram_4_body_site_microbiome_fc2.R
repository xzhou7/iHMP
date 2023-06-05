


###
no_function()

setwd(masstools::get_project_wd())
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
  (personalized_score$between_mean1 - personalized_score$family_mean1) /
  (personalized_score$between_mean1 - personalized_score$within_mean1)

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
  (personalized_score$between_mean1 - personalized_score$family_mean1) /
  (personalized_score$between_mean1 - personalized_score$within_mean1)

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
  (personalized_score$between_mean1 - personalized_score$family_mean1) /
  (personalized_score$between_mean1 - personalized_score$within_mean1)

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
  (personalized_score$between_mean1 - personalized_score$family_mean1) /
  (personalized_score$between_mean1 - personalized_score$within_mean1)

# personalized_score$fc2[which(personalized_score$fc2_p_adjust > 0.05 &
#                                personalized_score$fc3_p_adjust > 0.05)] <- NA

personalized_score$fc2[which(personalized_score$family_mean1 > personalized_score$between_mean1)] = NA

nasal_personalized_score = personalized_score

dir.create("data_analysis/combine_microbiome")
dir.create("data_analysis/combine_microbiome/cladogram/")
setwd("data_analysis/combine_microbiome/cladogram/")

temp_data =
  rbind(
    stool_personalized_score,
    skin_personalized_score,
    oral_personalized_score,
    nasal_personalized_score
  ) %>%
  dplyr::rename(class = dataset) %>%
  dplyr::mutate(class = stringr::str_to_sentence(class)) %>%
  dplyr::filter(!is.na(fc2))

personalized_score <-
  temp_data

library(ggvenn)
library(ggVennDiagram)

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
      dplyr::pull(Genus)) %>%
  unique()

if (length(remove_name) > 0) {
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

# raw_p

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
raw_p$data$nodeSize2[raw_p$data$nodeClass2 == "Phylum"] = 2.5
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

# raw_p

###label 2 is the name of taxa
raw_p$data$label2 =
  raw_p$data$label %>%
  stringr::str_split("__") %>%
  purrr::map(function(x) {
    x[2]
  }) %>%
  unlist()

raw_p$data$label2[as.character(raw_p$data$nodeClass) != "g"] = NA

#####add new information
temp_data_fc2 =
  temp_data %>%
  dplyr::select(genus, class, fc2) %>%
  tidyr::pivot_wider(names_from = class, values_from = "fc2")

colnames(temp_data_fc2)[-1] = paste(colnames(temp_data_fc2)[-1], "fc2", sep = "_")

new_info =
  data.frame(Genus = raw_p$data$label2) %>%
  dplyr::left_join(as.data.frame(variable_info)[, c("Genus", "Kingdom", "Phylum", "Class", "Order", "Family")],
                   by = "Genus") %>%
  dplyr::left_join(temp_data_fc2, by = c("Genus" = "genus")) %>%
  dplyr::mutate(Stool = case_when(!is.na(Stool_fc2) ~ "Stool",
                                  TRUE ~ "no")) %>%
  dplyr::mutate(Skin = case_when(!is.na(Skin_fc2) ~ "Skin",
                                 TRUE ~ "no")) %>%
  dplyr::mutate(Oral = case_when(!is.na(Oral_fc2) ~ "Oral",
                                 TRUE ~ "no")) %>%
  dplyr::mutate(Nasal = case_when(!is.na(Nasal_fc2) ~ "Nasal",
                                  TRUE ~ "no"))

raw_p$data =
  cbind(raw_p$data, new_info)

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

# p1

######add heatmap
##add fc2 information
p1$data$Stool_fc2[which(p1$data$Stool_fc2 < 0)] = 0
p1$data$Skin_fc2[which(p1$data$Skin_fc2 < 0)] = 0
p1$data$Oral_fc2[which(p1$data$Oral_fc2 < 0)] = 0
p1$data$Nasal_fc2[which(p1$data$Nasal_fc2 < 0)] = 0

p1$data$Stool_fc2[which(p1$data$Stool_fc2 > 1)] = 1
p1$data$Skin_fc2[which(p1$data$Skin_fc2 > 1)] = 1
p1$data$Oral_fc2[which(p1$data$Oral_fc2 > 1)] = 1
p1$data$Nasal_fc2[which(p1$data$Nasal_fc2 > 1)] = 1

##add multiple tip information
stool_fc2_info =
  p1$data[, c("Stool_fc2", "isTip")]
rownames(stool_fc2_info) = p1$data$label
stool_fc2_info = stool_fc2_info %>%
  dplyr::filter(isTip) %>%
  dplyr::select(-isTip)

skin_fc2_info =
  p1$data[, c("Skin_fc2", "isTip")]
rownames(skin_fc2_info) = p1$data$label
skin_fc2_info = skin_fc2_info %>%
  dplyr::filter(isTip) %>%
  dplyr::select(-isTip)

oral_fc2_info =
  p1$data[, c("Oral_fc2", "isTip")]
rownames(oral_fc2_info) = p1$data$label
oral_fc2_info = oral_fc2_info %>%
  dplyr::filter(isTip) %>%
  dplyr::select(-isTip)

nasal_fc2_info =
  p1$data[, c("Nasal_fc2", "isTip")]
rownames(nasal_fc2_info) = p1$data$label
nasal_fc2_info = nasal_fc2_info %>%
  dplyr::filter(isTip) %>%
  dplyr::select(-isTip)

range(
  c(
    p1$data$Stool_fc2,
    p1$data$Skin_fc2,
    p1$data$Nasal_fc2,
    p1$data$Oral_fc2
  ),
  na.rm = TRUE
)

p1 = p1 +
  ggnewscale::new_scale_fill()

p2_1 =
  gheatmap(
    p = p1,
    data = stool_fc2_info,
    offset = -0.1,
    width = .08,
    colnames_angle = 95,
    colnames_offset_y = .25,
    colnames = FALSE,
    color = alpha(body_site_color["Stool"], 1),
    legend_title = "Index1"
  ) +
  scale_fill_gradientn(colours = rev(c(
    RColorBrewer::brewer.pal(n = 11, name = "Spectral")
  )[-c(5:7)]),
  na.value = "white",
  limits = c(0, 1)) +
  ggnewscale::new_scale(new_aes = "fill")

p2_1

p2_2 =
  gheatmap(
    p = p2_1,
    data = skin_fc2_info,
    offset = 0.5,
    width = .08,
    colnames_angle = 95,
    colnames_offset_y = .25,
    colnames = FALSE,
    color = alpha(body_site_color["Skin"], 1),
    legend_title = "Index1"
  ) +
  scale_fill_gradientn(colours = rev(c(
    RColorBrewer::brewer.pal(n = 11, name = "Spectral")
  )[-c(5:7)]),
  na.value = "white",
  limits = c(0, 1)) +
  ggnewscale::new_scale(new_aes = "fill")

# p2_2

p2_3 =
  gheatmap(
    p = p2_2,
    data = oral_fc2_info,
    offset = 1.1,
    width = .08,
    colnames_angle = 95,
    colnames_offset_y = .25,
    colnames = FALSE,
    color = alpha(body_site_color["Oral"], 1),
    legend_title = "Index1"
  ) +
  scale_fill_gradientn(colours = rev(c(
    RColorBrewer::brewer.pal(n = 11, name = "Spectral")
  )[-c(5:7)]),
  na.value = "white",
  limits = c(0, 1)) +
  ggnewscale::new_scale(new_aes = "fill")

# p2_3

p2_4 =
  gheatmap(
    p = p2_3,
    data = nasal_fc2_info,
    offset = 1.7,
    width = .08,
    colnames_angle = 95,
    colnames_offset_y = .25,
    colnames = FALSE,
    color = alpha(body_site_color["Nasal"], 1),
    legend_title = "Index1"
  ) +
  scale_fill_gradientn(colours = rev(c(
    RColorBrewer::brewer.pal(n = 11, name = "Spectral")
  )[-c(5:7)]),
  na.value = "white",
  limits = c(0, 1)) +
  ggnewscale::new_scale(new_aes = "fill")

# p2_4

##add tip lab
##only add some tip points
# idx1 = which(!is.na(raw_p$data$Stool_fc2_p_adjust) & raw_p$data$Stool_fc2_p_adjust < 0.001)
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

# ggsave(p3, filename = "fc2_cladogram.pdf", width = 14, height = 14)


###output results
temp_data =
  p3$data

temp_data =
  temp_data %>%
  dplyr::filter(isTip)

# openxlsx::write.xlsx(temp_data,
#                      file = "fc2_tree_data.xlsx",
#                      asTable = TRUE,
#                      overwrite = TRUE)

#####summary information
plot =
  personalized_score %>%
  dplyr::mutate(class = factor(class, levels = c("Stool", "Skin", "Oral", "Nasal"))) %>%
  dplyr::mutate(fc2 = case_when(fc2 > 1 ~ 1,
                                fc2 < 0 ~ 0,
                                TRUE ~ fc2)) %>%
  ggplot(aes(class, fc2)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(fill = class),
              alpha = 0.9,
              size = 6,
              shape = 21) +
  guides(fill = guide_legend(title = "-log10(p-adjust)")) +
  scale_fill_manual(values = c(body_site_color, "no" = "white")) +
  base_theme +
  labs(x = "", y = "Family score")

plot

# ggsave(plot, filename = "family_score_boxplot.pdf", width = 10, height = 7)

temp <-
  personalized_score %>%
  dplyr::mutate(class = factor(class, levels = c("Stool", "Skin", "Oral", "Nasal"))) %>%
  dplyr::mutate(fc2 = case_when(fc2 > 1 ~ 1,
                                fc2 < 0 ~ 0,
                                TRUE ~ fc2))

stool_vs_skin_test <-
  wilcox.test(temp$fc2[temp$class == "Stool"],
              temp$fc2[temp$class == "Skin"])

stool_vs_oral_test <-
  wilcox.test(temp$fc2[temp$class == "Stool"],
              temp$fc2[temp$class == "Oral"])

stool_vs_nasal_test <-
  wilcox.test(temp$fc2[temp$class == "Stool"],
              temp$fc2[temp$class == "Nasal"])

skin_vs_oral_test <-
  wilcox.test(temp$fc2[temp$class == "Skin"],
              temp$fc2[temp$class == "Oral"])

skin_vs_nasal_test <-
  wilcox.test(temp$fc2[temp$class == "Skin"],
              temp$fc2[temp$class == "Nasal"])

oral_vs_nasal_test <-
  wilcox.test(temp$fc2[temp$class == "Oral"],
              temp$fc2[temp$class == "Nasal"])


sink(file = "family_score_test_result.txt")
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
  mean(temp$fc2[temp$class == "Stool"])
stool_quantile <-
  quantile(temp$fc2[temp$class == "Stool"])

skin_mean <-
  mean(temp$fc2[temp$class == "Skin"])
skin_quantile <-
  quantile(temp$fc2[temp$class == "Skin"])

oral_mean <-
  mean(temp$fc2[temp$class == "Oral"])
oral_quantile <-
  quantile(temp$fc2[temp$class == "Oral"])

nasal_mean <-
  mean(temp$fc2[temp$class == "Nasal"])
nasal_quantile <-
  quantile(temp$fc2[temp$class == "Nasal"])


sink(file = "family_score_median_mean_value.txt")
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

personalized_score$genus

library(plyr)

variable_info <-
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

remain_phylum <-
  remain_phylum[remain_phylum != "Cyanobacteria/Chloroplast"]


temp_data <-
  personalized_score %>%
  dplyr::left_join(variable_info, by = c("genus" = "Genus")) %>%
  dplyr::mutate(class = factor(class, levels = c("Stool", "Skin", "Oral", "Nasal"))) %>%
  dplyr::mutate(fc2 = round(fc2, 2)) %>%
  dplyr::mutate(fc2 = case_when(fc2 > 1 ~ 1,
                                fc2 < 0 ~ 0,
                                TRUE ~ fc2)) %>%
  dplyr::mutate(
    Phylum = case_when(
      Phylum == "Actinobacteria" |
        Phylum == "Bacteroidetes" |
        Phylum == "Firmicutes" |
        Phylum == "Proteobacteria" ~ Phylum,
      TRUE ~ "Other"
    )
  ) %>%
  dplyr::mutate(Phylum = factor(
    Phylum,
    levels = c(
      "Actinobacteria",
      "Bacteroidetes",
      "Firmicutes",
      "Fusobacteria",
      "Proteobacteria",
      "Other"
    )
  ))

plot <-
  temp_data %>%
  ggplot(aes(Phylum, fc2)) +
  geom_jitter(
    aes(fill = Phylum),
    size = 3.5,
    shape = 21,
    alpha = 0.7,
    show.legend = FALSE
  ) +
  geom_boxplot(outlier.shape = NA, fill = "transparent") +
  scale_fill_manual(values = c(phylum_color, "Other" = "grey")) +
  # guides(size = guide_legend(title = "-log10(p-adjust)")) +
  base_theme +
  # scale_size_continuous(range = c(0.05, 2)) +
  labs(x = "", y = "Family score") +
  facet_wrap(facets = vars(class),
             # scales = "free_y",
             nrow = 1) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 10
    ),
    strip.text = element_text(size = 10)
  )

plot

ggsave(plot,
       filename = "family_score_boxplot_based_on_phylum.pdf",
       width = 10,
       height = 4)


colnames(temp_data)

temp_data <-
  temp_data %>%
  dplyr::select(-c(
    fc1_p,
    fc2_p,
    fc3_p,
    fc1_p_adjust,
    fc2_p_adjust,
    fc3_p_adjust,
    fc1
  )) %>%
  dplyr::rename(family_score = fc2)

write.csv(temp_data, "family_score_based_on_phylum.csv", row.names = FALSE)
