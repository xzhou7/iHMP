#Inter_bodysite_microbiome correlation (Fig3d)

library(xxx)

library(phyloseq)
library(ggplot2)
library(reshape2)
library(dplyr)
library(cowplot)
library(stringr)
library(tidyverse)
library(broom)
library(lme4)
library(coefplot2)
library(reghelper)
library(jtools)
library(lmerTest)
library(ggpubr)
library(survival)
library(survminer)
library(broom)
library(gridExtra)
library(patchwork)
library(igraph)
library(ggraph)
library(tidygraph)
library(corrr)

base_theme = theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 13),
        panel.grid.minor = element_blank())

body_site_color = c( "Stool" = ggsci::pal_jama()(n=7)[2],
                     "Skin" = ggsci::pal_jama()(n=7)[3],
                     "Oral" = ggsci::pal_jama()(n=7)[4],
                     "Nasal" = ggsci::pal_jama()(n=7)[5])

body_site_color2= c("stool" = ggsci::pal_jama()(n=7)[2],
                    "skin" = ggsci::pal_jama()(n=7)[3],
                    "oral" = ggsci::pal_jama()(n=7)[4],
                    "nasal" = ggsci::pal_jama()(n=7)[5])

setwd("~/Library/CloudStorage/Box-Box/XinZhouFiles/Projects/ZXP2_HMP2_Core/")
source("./Analysis/Analysis Script/sxt.tools.R")
load("../../../human_microbiome_project/data_analysis/correlation_network/whole_data_set/intra_microbiome_correlation/sample_cor")
load("../../../human_microbiome_project/data_analysis/correlation_network/whole_data_set/intra_microbiome_correlation/subject_cor")
load("./Analysis/Robject/DetailedPhyloseq.RData")

ASV_table <- tax_table(physeq_ST) %>% data.frame()
ASV_table$ASV <- row.names(ASV_table)

OTU_table <- tax_table(physeq_SK) %>% data.frame()
OTU_table$ASV <- row.names(OTU_table)

taxa_table <- rbind(ASV_table,OTU_table)

sample_cor$from_taxa <- sample_cor$from
sample_cor$to_taxa <- sample_cor$to

sample_cor$from_taxa <- str_remove(sample_cor$from_taxa,"stool")
sample_cor$from_taxa <- str_remove(sample_cor$from_taxa,"skin")
sample_cor$from_taxa <- str_remove(sample_cor$from_taxa,"oral")
sample_cor$from_taxa <- str_remove(sample_cor$from_taxa,"nasal")

sample_cor$to_taxa <- str_remove(sample_cor$to_taxa,"stool")
sample_cor$to_taxa <- str_remove(sample_cor$to_taxa,"skin")
sample_cor$to_taxa <- str_remove(sample_cor$to_taxa,"oral")
sample_cor$to_taxa <- str_remove(sample_cor$to_taxa,"nasal")

subject_cor$from_taxa <- subject_cor$from
subject_cor$to_taxa <- subject_cor$to

subject_cor$from_taxa <- str_remove(subject_cor$from_taxa,"stool")
subject_cor$from_taxa <- str_remove(subject_cor$from_taxa,"skin")
subject_cor$from_taxa <- str_remove(subject_cor$from_taxa,"oral")
subject_cor$from_taxa <- str_remove(subject_cor$from_taxa,"nasal")

subject_cor$to_taxa <- str_remove(subject_cor$to_taxa,"stool")
subject_cor$to_taxa <- str_remove(subject_cor$to_taxa,"skin")
subject_cor$to_taxa <- str_remove(subject_cor$to_taxa,"oral")
subject_cor$to_taxa <- str_remove(subject_cor$to_taxa,"nasal")


sample_cor_taxa0 <- merge(sample_cor,taxa_table, by.x = "from_taxa", by.y = "ASV", all.x=T)
colnames(sample_cor_taxa0)[12:18] <- paste("from",colnames(sample_cor_taxa0)[12:18], sep="_")

sample_cor_taxa <- merge(sample_cor_taxa0,taxa_table, by.x = "to_taxa", by.y = "ASV", all.x=T)
colnames(sample_cor_taxa)[19:25] <- paste("to",colnames(sample_cor_taxa)[19:25], sep="_")

subject_cor_taxa0 <-  merge(subject_cor,taxa_table, by.x = "from_taxa", by.y = "ASV", all.x=T)
colnames(subject_cor_taxa0)[11:17] <- paste("from",colnames(subject_cor_taxa0)[11:17], sep="_")

subject_cor_taxa <- merge(subject_cor_taxa0,taxa_table, by.x = "to_taxa", by.y = "ASV", all.x=T)
colnames(subject_cor_taxa)[18:24] <- paste("to",colnames(subject_cor_taxa)[18:24], sep="_")

sample_cor_taxa$compare <- "between"
sample_cor_taxa$compare[which(sample_cor_taxa$from_class == sample_cor_taxa$to_class)] <- "within"

subject_cor_taxa$compare <- "between"
subject_cor_taxa$compare[which(subject_cor_taxa$from_class == subject_cor_taxa$to_class)] <- "within"

#############
stool_table <- tax_table(filter_taxa(physeqGenus_ST, function(x) sum(x) > 0, T)) %>% data.frame()
skin_table <- tax_table(filter_taxa(physeqGenus_SK, function(x) sum(x) > 0, T)) %>% data.frame()
oral_table <- tax_table(filter_taxa(physeqGenus_OR, function(x) sum(x) > 0, T)) %>% data.frame()
nasal_table <- tax_table(filter_taxa(physeqGenus_NS, function(x) sum(x) > 0, T)) %>% data.frame()

ST_phy_freq <- table(stool_table$Phylum) %>% data.frame()
SK_phy_freq <- table(skin_table$Phylum) %>% data.frame()
OR_phy_freq <- table(oral_table$Phylum) %>% data.frame()
NS_phy_freq <- table(nasal_table$Phylum) %>% data.frame()

sample_cor_taxa$compareclass <- paste(sample_cor_taxa$from_class,sample_cor_taxa$to_class, sep="_")

sample_cor_taxa <- sample_cor_taxa %>% filter(cor < 1)
subject_cor_taxa <- subject_cor_taxa %>% filter(cor < 1)
# dim(filter(sample_cor_taxa,compare == "within") %>% filter(from_class == "stool"))
# dim(filter(sample_cor_taxa,compare == "within" & from_class == "skin"))
# dim(filter(sample_cor_taxa,compare == "within" & from_class == "oral"))
# dim(filter(sample_cor_taxa,compare == "within" & from_class == "nasal"))
# 
# dim(filter(subject_cor_taxa,compare == "within" & from_class == "stool"))
# dim(filter(subject_cor_taxa,compare == "within" & from_class == "skin"))
# dim(filter(subject_cor_taxa,compare == "within" & from_class == "oral"))
# dim(filter(subject_cor_taxa,compare == "within" & from_class == "nasal"))
# 
# dim(filter(sample_cor_taxa, compare == "between" & from_class == "stool"))
# dim(filter(sample_cor_taxa,compare == "between" & to_class == "stool"))
# 
# dim(filter(sample_cor_taxa,compare == "within" & from_class == "stool"))[1]/dim(stool_table)[1]^2
# dim(filter(sample_cor_taxa,compare == "within" & from_class == "skin"))[1]/dim(skin_table)[1]^2
# dim(filter(sample_cor_taxa,compare == "within" & from_class == "oral"))[1]/dim(oral_table)[1]^2
# dim(filter(sample_cor_taxa,compare == "within" & from_class == "nasal"))[1]/dim(nasal_table)[1]^2
# 
# write.csv(file = "../../../human_microbiome_project/Figures/Figure3/detailed/inter-bodysite-corre/sample_wise_correlation.csv",sample_cor_taxa)
# write.csv(file = "../../../human_microbiome_project/Figures/Figure3/detailed/inter-bodysite-corre/sample_wise_within_correlation.csv",filter(sample_cor_taxa,compare=="within"))
# write.csv(file = "../../../human_microbiome_project/Figures/Figure3/detailed/inter-bodysite-corre/sample_wise_between_correlation.csv",filter(sample_cor_taxa,compare=="between"))
# 
# write.csv(file = "../../../human_microbiome_project/Figures/Figure3/detailed/inter-bodysite-corre/subject_wise_correlation.csv",subject_cor_taxa)
# write.csv(file = "../../../human_microbiome_project/Figures/Figure3/detailed/inter-bodysite-corre/subject_wise_within_correlation.csv",filter(subject_cor_taxa,compare=="within"))
# write.csv(file = "../../../human_microbiome_project/Figures/Figure3/detailed/inter-bodysite-corre/subject_wise_between_correlation.csv",filter(subject_cor_taxa,compare=="between"))

######################################################################################
#set several examples for the correlation 


dim(filter(sample_cor_taxa,compare == "within" & from_class == "stool"))

table(sample_cor_taxa$from_class,sample_cor_taxa$compare)

ggplot(sample_cor_taxa, aes(x=cor)) + geom_density() + facet_wrap(.~compareclass) + ggtitle("Sample.Wise.Corr.Coe.Density") + base_theme
ggplot(subject_cor_taxa,aes(x=cor)) + geom_density() + facet_grid(.~compare) + ggtitle("Subject.Wise.Corr.Coe.Density")+ base_theme
 
# p.stooltostool <- filter(sample_cor_taxa, cor > 0.2 | cor < -0.2) %>%  filter(compareclass == "stool_stool") %>% 
#   select(from_Genus, to_Genus, cor) %>%
#   as_tbl_graph(directed = FALSE) %>%
#   ggraph(layout = "fr") + 
#   geom_edge_link(colour = "lightgray") + 
#   geom_node_point() +
#   geom_node_text(aes(label = name), size = 5, repel = TRUE, color="red") +
#   theme_graph() + ggtitle("Stool")
# p.stooltostool 
# #ggsave(filename = "../../../human_microbiome_project/Figures/Figure3/detailed/inter-bodysite-corre/stool.to.stool.sample.pdf", p.stooltostool, width = 10, height = 10, dpi = 300)
# 
# p.skintoskin <- filter(sample_cor_taxa, cor > 0.3 | cor < -0.3) %>%  filter(compareclass == "skin_skin") %>% 
#   select(from_Genus, to_Genus, cor) %>%
#   as_tbl_graph(directed = FALSE) %>%
#   ggraph(layout = "fr") + 
#   geom_edge_link(colour = "lightgray") + 
#   geom_node_point() +
#   geom_node_text(aes(label = name), size = 5, repel = TRUE, color="red") +
#   theme_graph() + ggtitle("Skin")
# p.skintoskin 
# #ggsave(filename = "../../../human_microbiome_project/Figures/Figure3/detailed/inter-bodysite-corre/skin.to.skin.sample.pdf", p.skintoskin, width = 10, height = 10, dpi = 300)
# 
# p.oraltooral <- filter(sample_cor_taxa, cor > 0.3 | cor < -0.3) %>%  filter(compareclass == "oral_oral") %>% 
#   select(from_Genus, to_Genus, cor) %>%
#   as_tbl_graph(directed = FALSE) %>%
#   ggraph(layout = "fr") + 
#   geom_edge_link(colour = "lightgray") + 
#   geom_node_point() +
#   geom_node_text(aes(label = name), size = 5, repel = TRUE, color="red") +
#   theme_graph() + ggtitle("oral")
# p.oraltooral 
# #ggsave(filename = "../../../human_microbiome_project/Figures/Figure3/detailed/inter-bodysite-corre/oral.to.oral.sample.pdf", p.oraltooral, width = 10, height = 10, dpi = 300)
# 
# p.nasaltonasal <- filter(sample_cor_taxa, cor > 0.2 | cor < -0.2) %>%  filter(compareclass == "nasal_nasal") %>% 
#   select(from_Genus, to_Genus, cor) %>%
#   as_tbl_graph(directed = FALSE) %>%
#   ggraph(layout = "fr") + 
#   geom_edge_link(colour = "lightgray") + 
#   geom_node_point() +
#   geom_node_text(aes(label = name), size = 5, repel = TRUE, color="red") +
#   theme_graph() + ggtitle("Nasal")
# p.nasaltonasal
# #ggsave(filename = "../../../human_microbiome_project/Figures/Figure3/detailed/inter-bodysite-corre/nasal.to.nasal.sample.pdf", p.nasaltonasal, width = 10, height = 10, dpi = 300)

sample_cor_taxa[str_detect(sample_cor_taxa$from_Genus, pattern = "Staphylococcus") & str_detect(sample_cor_taxa$to_Genus, pattern = "Staphylococcus"),]

subject_cor_taxa[str_detect(subject_cor_taxa$from_Genus, pattern = "Staphylococcus") & str_detect(subject_cor_taxa$to_Genus, pattern = "Staphylococcus"),]

table(sample_cor_taxa$from_Genus) %>% sort

filter(sample_cor_taxa, from_class!=to_class) %>% filter(from_Genus == to_Genus)
filter(sample_cor_taxa, from_class!=to_class) %>% filter(from_Genus == "Staphylococcus")
filter(sample_cor_taxa, from_class!=to_class) %>% filter(to_Genus == "Staphylococcus")

filter(sample_cor_taxa, from_class!=to_class) %>% filter(from_Genus == "Cutibacterium")
filter(sample_cor_taxa, from_class!=to_class) %>% filter(to_Genus == "Cutibacterium")

filter(sample_cor_taxa, from_class!=to_class) %>% filter(from_Genus == "Corynebacterium")
filter(sample_cor_taxa, from_class!=to_class) %>% filter(to_Genus == "Corynebacterium")

filter(sample_cor_taxa, from_class!=to_class) %>% filter(from_Genus == "Streptococcus")
filter(sample_cor_taxa, from_class!=to_class) %>% filter(to_Genus == "Streptococcus")

filter(sample_cor_taxa, from_class!=to_class) %>% filter(from_Genus == "Tenericutes")
filter(sample_cor_taxa, from_class!=to_class) %>% filter(to_Genus == "Tenericutes")


filter(subject_cor_taxa, from_class!=to_class) %>% filter(from_Genus == to_Genus)

filter(subject_cor_taxa, from_class!=to_class) %>% filter(from_Genus == "Bacteroides")
filter(subject_cor_taxa, from_class!=to_class) %>% filter(to_Genus == "Bacteroides")

filter(subject_cor_taxa, from_class!=to_class) %>% filter(from_Genus == "Unclassified_Ruminococcaceae")
filter(subject_cor_taxa, from_class!=to_class) %>% filter(to_Genus == "Unclassified_Ruminococcaceae")

filter(subject_cor_taxa, from_class!=to_class) %>% filter(from_Genus == "Prevotella")
filter(subject_cor_taxa, from_class!=to_class) %>% filter(to_Genus == "Prevotella")

filter(subject_cor_taxa, from_class!=to_class) %>% filter(from_Genus == "Leptotrichia")
filter(subject_cor_taxa, from_class!=to_class) %>% filter(to_Genus == "Leptotrichia")

filter(subject_cor_taxa, from_class!=to_class) %>% filter(from_Genus == "Cutibacterium")
filter(subject_cor_taxa, from_class!=to_class) %>% filter(to_Genus == "Cutibacterium")

filter(subject_cor_taxa, from_class!=to_class) %>% filter(from_Genus == "Corynebacterium")
filter(subject_cor_taxa, from_class!=to_class) %>% filter(to_Genus == "Corynebacterium")

filter(subject_cor_taxa, from_class!=to_class) %>% filter(from_Genus == "Staphylococcus")
filter(subject_cor_taxa, from_class!=to_class) %>% filter(to_Genus == "Staphylococcus")


filter(subject_cor_taxa, from_class!=to_class) %>% filter(from_Genus == "Tenericutes")

write.csv(file = "../../../human_microbiome_project/Supplementary_data/Extended Data Table subjectcor.csv",subject_cor_taxa )


