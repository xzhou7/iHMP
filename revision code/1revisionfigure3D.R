#orgnize tables for figure 3
library(phyloseq)
library(dplyr)

#set directory
setwd("~/Library/CloudStorage/Box-Box/human_microbiome_project/data_analysis/correlation_network/whole_data_set/")

# old version
load("./intra_stool_microbiome/intra_stool_microbiome_lm_adjusted_cor")
# load("./intra_skin_microbiome/intra_skin_microbiome_lm_adjusted_cor")
load("./intra_oral_microbiome/intra_oral_microbiome_lm_adjusted_cor")
# load("./intra_nasal_microbiome/intra_nasal_microbiome_lm_adjusted_cor")
# 
# load("./stool_microbiome_vs_skin_microbiome/stool_microbiome_skin_microbiome_lm_adjusted_cor_spearman")
# load("./stool_microbiome_vs_oral_microbiome/stool_microbiome_oral_microbiome_lm_adjusted_cor_spearman")
# load("./stool_microbiome_vs_nasal_microbiome/stool_microbiome_nasal_microbiome_lm_adjusted_cor_spearman")
# load("./skin_microbiome_vs_nasal_microbiome/skin_microbiome_nasal_microbiome_lm_adjusted_cor_spearman")
# load("./oral_microbiome_vs_nasal_microbiome/oral_microbiome_nasal_microbiome_lm_adjusted_cor_spearman")
# load("./oral_microbiome_vs_skin_microbiome/oral_microbiome_skin_microbiome_lm_adjusted_cor_spearman")

#new version
sample_wise_corre_betweenbodysite <- read.csv("./intra_microbiome_correlation/sample_cor_output.csv", header = T)
sample_wise_corre_withinbodysite <- read.csv("./intra_microbiome_correlation/intra.csv", header = T)

sample_wise_corre_withinbodysite$comparison <- "within_bodysite"
sample_wise_corre_betweenbodysite$comparison <- "between_bodysite"

load("./intra_microbiome_correlation/subject_cor")

load("../../../../XinZhouFiles/Projects/ZXP2_HMP2_Core/Analysis/Robject/Prevalance.RData")
load("../../../../XinZhouFiles/Projects/ZXP2_HMP2_Core/Analysis/Robject/DetailedPhyloseq.RData")

table(intra_stool_microbiome_lm_adjusted_cor$p_adjust < 0.05)
table(intra_skin_microbiome_lm_adjusted_cor$p_adjust < 0.05)
table(intra_oral_microbiome_lm_adjusted_cor$p_adjust < 0.05)
table(intra_nasal_microbiome_lm_adjusted_cor$p_adjust < 0.05)

ASV_mapping <- tax_table(physeq_ST) %>% data.frame()
ASV_mapping$ASV <- rownames(ASV_mapping)
OTU_mapping <- tax_table(physeq_OR) %>% data.frame()
OTU_mapping$OTU <- rownames(OTU_mapping)

length(table(intra_stool_microbiome_lm_adjusted_cor$microbiome))

#intra-individual
#within985.17 
intra_stool_microbiome_lm_adjusted_cor$from_genus <- ASV_mapping$Genus[match(intra_stool_microbiome_lm_adjusted_cor$microbiome, ASV_mapping$ASV)]
intra_stool_microbiome_lm_adjusted_cor$to_genus <- ASV_mapping$Genus[match(intra_stool_microbiome_lm_adjusted_cor$metabolite, ASV_mapping$ASV)]

intra_skin_microbiome_lm_adjusted_cor$from_genus <- OTU_mapping$Genus[match(intra_skin_microbiome_lm_adjusted_cor$microbiome, OTU_mapping$OTU)]
intra_skin_microbiome_lm_adjusted_cor$to_genus <- OTU_mapping$Genus[match(intra_skin_microbiome_lm_adjusted_cor$metabolite, OTU_mapping$OTU)]

intra_oral_microbiome_lm_adjusted_cor$from_genus <- OTU_mapping$Genus[match(intra_oral_microbiome_lm_adjusted_cor$microbiome, OTU_mapping$OTU)]
intra_oral_microbiome_lm_adjusted_cor$to_genus <- OTU_mapping$Genus[match(intra_oral_microbiome_lm_adjusted_cor$metabolite, OTU_mapping$OTU)]

intra_nasal_microbiome_lm_adjusted_cor$from_genus <- ASV_mapping$Genus[match(intra_nasal_microbiome_lm_adjusted_cor$microbiome, ASV_mapping$ASV)]
intra_nasal_microbiome_lm_adjusted_cor$to_genus <- ASV_mapping$Genus[match(intra_nasal_microbiome_lm_adjusted_cor$metabolite, ASV_mapping$ASV)]

#between
stool_microbiome_nasal_microbiome_lm_adjusted_cor_spearman$from_genus <-  ASV_mapping$Genus[match(stool_microbiome_nasal_microbiome_lm_adjusted_cor_spearman$microbiome, ASV_mapping$ASV)]
stool_microbiome_nasal_microbiome_lm_adjusted_cor_spearman$to_genus <- ASV_mapping$Genus[match(stool_microbiome_nasal_microbiome_lm_adjusted_cor_spearman$metabolite, ASV_mapping$ASV)]

stool_microbiome_skin_microbiome_lm_adjusted_cor_spearman$from_genus <- ASV_mapping$Genus[match(stool_microbiome_skin_microbiome_lm_adjusted_cor_spearman$microbiome, ASV_mapping$ASV)]
stool_microbiome_skin_microbiome_lm_adjusted_cor_spearman$to_genus <- OTU_mapping$Genus[match(stool_microbiome_skin_microbiome_lm_adjusted_cor_spearman$metabolite, OTU_mapping$OTU)]

stool_microbiome_oral_microbiome_lm_adjusted_cor_spearman$from_genus <- ASV_mapping$Genus[match(stool_microbiome_oral_microbiome_lm_adjusted_cor_spearman$microbiome, ASV_mapping$ASV)]
stool_microbiome_oral_microbiome_lm_adjusted_cor_spearman$to_genus <- OTU_mapping$Genus[match(stool_microbiome_oral_microbiome_lm_adjusted_cor_spearman$metabolite, OTU_mapping$OTU)]

skin_microbiome_nasal_microbiome_lm_adjusted_cor_spearman$from_genus <- OTU_mapping$Genus[match(skin_microbiome_nasal_microbiome_lm_adjusted_cor_spearman$microbiome, OTU_mapping$OTU)]
skin_microbiome_nasal_microbiome_lm_adjusted_cor_spearman$to_genus <- ASV_mapping$Genus[match(skin_microbiome_nasal_microbiome_lm_adjusted_cor_spearman$metabolite, ASV_mapping$ASV)]

oral_microbiome_skin_microbiome_lm_adjusted_cor_spearman$from_genus <- OTU_mapping$Genus[match(oral_microbiome_skin_microbiome_lm_adjusted_cor_spearman$microbiome, OTU_mapping$OTU)]
oral_microbiome_skin_microbiome_lm_adjusted_cor_spearman$to_genus <- OTU_mapping$Genus[match(oral_microbiome_skin_microbiome_lm_adjusted_cor_spearman$metabolite, OTU_mapping$OTU)]

oral_microbiome_nasal_microbiome_lm_adjusted_cor_spearman$from_genus <- OTU_mapping$Genus[match(oral_microbiome_nasal_microbiome_lm_adjusted_cor_spearman$microbiome, OTU_mapping$OTU)]
oral_microbiome_nasal_microbiome_lm_adjusted_cor_spearman$to_genus <- ASV_mapping$Genus[match(oral_microbiome_nasal_microbiome_lm_adjusted_cor_spearman$metabolite, ASV_mapping$ASV)]



table(intra_stool_microbiome_lm_adjusted_cor[intra_stool_microbiome_lm_adjusted_cor$from_genus == intra_stool_microbiome_lm_adjusted_cor$to_genus,]$cor) 

table(intra_stool_microbiome_lm_adjusted_cor$cor == 1) 


table(intra_stool_microbiome_lm_adjusted_cor$from_genus == intra_stool_microbiome_lm_adjusted_cor$to_genus)
table(intra_stool_microbiome_lm_adjusted_cor$cor == 1)

unique(intra_stool_microbiome_lm_adjusted_cor$microbiome) %>% length()




table(intra_stool_microbiome_lm_adjusted_cor$p == 0)

table(intra_stool_microbiome_lm_adjusted_cor$p_adjust < 0.05)


table(stool_microbiome_skin_microbiome_lm_adjusted_cor_spearman$p_adjust < 0.05)

table(intra_stool_microbiome_lm_adjusted_cor$cor < 1)

intra_stool_microbiome_lm_adjusted_cor$from_genus 


sample_wise_within <- read.csv("../../../Figures/Figure3/detailed/inter-bodysite-corre/sample_wise_within_correlation.csv", header = T)

table(sample_wise_within$compareclass)


subject_wise_within <- read.csv("../../../Figures/Figure3/detailed/inter-bodysite-corre/subject_wise_within_correlation.csv", header = T)
table(subject_wise_within$to_taxa)
table(subject_wise_within$from_taxa)

table(subject_wise_within$from_Genus == subject_wise_within$to_Genus)




