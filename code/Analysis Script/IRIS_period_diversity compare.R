#IRIS Diversity comparasion 
#this is to compare is IR and IS show different diversity in four body sites

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
library(pracma)

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
load("../../../human_microbiome_project/data_analysis/combine_microbiome/within_bodysite_sample_distance/dist")
load("./Analysis/Robject/Diversity_Datatable.RData")
meta.data <- read.csv("./Analysis/ori_meta_table/metadata.subject.csv")
meta.data2 <- read.csv("./Analysis/ori_meta_table/hmp2_infection_metadata_0323.csv", header = T)
Diversity.Bysample <- read.csv("./Analysis/Diversity Table/All_Diversity.csv", header = T, row.names = 1)
# 
ASV_diversity <- rbind(select(stool.nasal.ASV.diversity, camargo:Fisher,"SampleID","SubjectID","Date","SampleType","depth"),
                       select(skin.oral.ASV.diversity,camargo:Fisher,"SampleID","SubjectID","Date","SampleType","depth"))
# 
# Genus.diversity <- rbind(select(stool.nasal.Genus.diversity, camargo:Fisher,"SampleID","SubjectID","Date","SampleType","depth"),
#                          select(skin.oral.Genus.diversity,camargo:Fisher,"SampleID","SubjectID","Date","SampleType","depth"))

colnames(ASV_diversity)[1:14] <- paste(colnames(ASV_diversity)[1:14], "untrimed", sep="_")

p.diver <- ggplot(Diversity.Bysample, aes(x=ACE, y=pielou_e, group=Bodysite, color=Bodysite)) + stat_density_2d(adjust=1.5)
p.diver <- p.diver + base_theme + scale_color_manual(values = body_site_color) + labs(title = "Ecological Characters of Microbiome", x="Richness(ACE)", y="Evenness(Pielou)")
p.diver <- p.diver + theme(
  legend.position = c(.95, .05),
  legend.justification = c("right", "bottom"),
  legend.box.just = "right",
  legend.margin = margin(6, 6, 6, 6),
  legend.box.background = element_rect(color="black", size=0.5),)
p.diver

Diversity.Bysample.meta <- merge(Diversity.Bysample, select(ASV_diversity, -(1:14)), by = "row.names")
Diversity.Bysample.meta.subject <- merge(Diversity.Bysample.meta,meta.data, by= "SubjectID")

Diversity_metatable <- merge(Diversity.Bysample.meta.subject, meta.data2,by = "SampleID", all.x = T)
Diversity_metatable$Bodysite <- factor(Diversity_metatable$Bodysite, levels = c("Stool","Skin","Oral","Nasal"))

Diversity_metatable$CL4[is.na(Diversity_metatable$CL4)] <- "unknown"

table(Diversity_metatable$CL4)

#Calculate on how many individuals are IR and how many are IS
length(unique(filter(Diversity_metatable, IRIS == "IS")$SubjectID.x))
length(unique(filter(Diversity_metatable, IRIS == "IR")$SubjectID.x))
length(unique(filter(Diversity_metatable, IRIS == "Unknown")$SubjectID.x))

#Calculate the samples belong to IR or IS
dim(filter(Diversity_metatable, IRIS == "IS"))
dim(filter(Diversity_metatable, IRIS == "IR"))
dim(filter(Diversity_metatable, IRIS == "Unknown"))

#Diversity Metatable Date format correction
Diversity_metatable$CollectionDate <- as.Date(Diversity_metatable$CollectionDate, format = "%m/%d/%y")

p.diver.IS <- ggplot(filter(Diversity_metatable, IRIS != "Unknown"), aes(x=ACE, y=pielou_e, group=IRIS, color=IRIS)) + stat_density_2d(adjust=1.5)
p.diver.IS <- p.diver.IS + base_theme  + scale_color_manual(values = iris_color) +labs(title = "Ecological Characters of Microbiome IS/IR", x="Richness(ACE)", y="Evenness(Pielou)")
p.diver.IS <- p.diver.IS + theme(
  legend.position = c(.95, .05),
  legend.justification = c("right", "bottom"),
  legend.box.just = "right",
  legend.margin = margin(6, 6, 6, 6),
  legend.box.background = element_rect(color="black", size=0.5),)
p.diver.IS <- p.diver.IS + facet_wrap(.~Bodysite, scales = "free")
p.diver.IS

#Estimate the center for each body site
contours <- layer_data(p.diver.IS)
cen <- data.frame()
for (i in 1:4){
  object <- filter(contours, colour=="#BC3C29FF" & PANEL== i & piece == 2)
  pracma::poly_center(object$x, object$y)
  cen[i,1:2] <- pracma::poly_center(object$x, object$y)
}

for (i in 1:4){
  object <- filter(contours, colour=="#0072B5FF" & PANEL== i & piece == 2)
  pracma::poly_center(object$x, object$y)
  cen[i+4,1:2] <- pracma::poly_center(object$x, object$y)
}
cen$bodysite <- rep(rev(c("Stool","Skin","Oral","Nasal")),2)
cen$IRIS <- c(rep("IS", 4), rep("IR", 4))
cen

PdiversityIRIS <- ggplot(cen, aes(x=V1, y=V2, group=bodysite, color=IRIS)) + geom_point(size=1.5) + xlab("Richness") + ylab("Evenness") + base_theme  + scale_color_manual(values = iris_color)
PdiversityIRIS <- PdiversityIRIS + geom_line(aes(group = bodysite), color="black", lty=2)
PdiversityIRIS

#ggsave(filename = "./Analysis/Diversity Table/Diversity_IRIS.pdf", PdiversityIRIS,width = 4, height = 3, dpi = 300 )


#summarize 1/cv for 
table(Diversity_metatable$Bodysite)
Diversity_metatable_test <- Diversity_metatable %>% group_by(SubjectID.x,Bodysite) %>% 
  mutate(AEC_CV = 1/cv(ACE)) %>% 
  mutate(AEC_mean = mean(ACE)) %>% 
  mutate(Shannon_CV = 1/cv(Shannon)) %>% 
  mutate(Shannon_mean = mean(Shannon)) %>% 
  mutate(Evenness = mean(pielou_e))
  ungroup()


Diversity_metatable_test$AEC_CV
#to validate the cv function:
#sd(filter(Diversity_metatable, SubjectID.x == "69-001" & Bodysite == "Stool")$ACE)/mean(filter(Diversity_metatable, SubjectID.x == "69-001" & Bodysite == "Stool")$ACE)

ACE_CV_table <- unique(Diversity_metatable_test %>% dplyr::select(SubjectID.x,IRIS, AEC_CV,AEC_mean,Shannon_CV,Evenness,Shannon_mean,Bodysite,SSPG,BMI))
ACE_CV_table <- ACE_CV_table[!is.na(ACE_CV_table$AEC_CV),]

p <- ggplot(filter(ACE_CV_table, IRIS != "Unknown"), aes(x=factor(IRIS, levels=c("IS", "IR")), y=AEC_mean, fill=Bodysite)) + geom_jitter(size=0.5) +  geom_boxplot(alpha=0.4,outlier.alpha = 0)  + facet_wrap(.~Bodysite, scale= "free") + scale_y_continuous(trans='log10')
p <- p + scale_fill_manual(values = body_site_color) + base_theme + ggtitle("Shannon Diversity between IS/IR individuals") + ylab("Shannon Diversity") + xlab("")
p
#ggsave(filename = "./Analysis/Diversity Table/DiversityIRIS.pdf",p, width = 5, height = 4, dpi = 300)

p <- ggplot(filter(ACE_CV_table, IRIS != "Unknown"), aes(x=factor(IRIS, levels=c("IS", "IR")), y=Evenness, fill=Bodysite)) + geom_jitter(size=0.5) +  geom_boxplot(alpha=0.4,outlier.alpha = 0)  + facet_wrap(.~Bodysite, scale= "free") + scale_y_continuous(trans='log10')
p <- p + scale_fill_manual(values = body_site_color) + base_theme + ggtitle("Evenness between IS/IR individuals") + ylab("Pielou Evenness") + xlab("")
p
ggsave(filename = "./Analysis/Diversity Table/EVENNESS_IRIS.pdf",p, width = 5, height = 4, dpi = 300)

p <- ggplot(filter(ACE_CV_table, IRIS != "Unknown"), aes(x=factor(IRIS, levels=c("IS", "IR")), y=AEC_mean, fill=Bodysite)) + geom_jitter(size=0.5) +  geom_boxplot(alpha=0.4,outlier.alpha = 0)  + facet_wrap(.~Bodysite, scale= "free") + scale_y_continuous(trans='log10')
p <- p + scale_fill_manual(values = body_site_color) + base_theme + ggtitle("ACE Richness between IS/IR individuals") + ylab("Shannon Diversity") + xlab("")
p
ggsave(filename = "./Analysis/Diversity Table/ACE_IRIS.pdf",p, width = 5, height = 4, dpi = 300)

p <- ggplot(filter(ACE_CV_table, IRIS != "Unknown"), aes(x=IRIS, y=AEC_CV)) + geom_point() + geom_boxplot() + facet_wrap(.~Bodysite,scale= "free") +  scale_y_continuous(trans='log10')
p

#test the shannon diversity difference
t.test(filter(ACE_CV_table, IRIS=="IS" & Bodysite=="Stool")$Shannon_mean, filter(ACE_CV_table, IRIS=="IR" & Bodysite=="Stool")$Shannon_mean)
t.test(filter(ACE_CV_table, IRIS=="IS" & Bodysite=="Skin")$Shannon_mean, filter(ACE_CV_table, IRIS=="IR" & Bodysite=="Skin")$Shannon_mean)
t.test(filter(ACE_CV_table, IRIS=="IS" & Bodysite=="Oral")$Shannon_mean, filter(ACE_CV_table, IRIS=="IR" & Bodysite=="Oral")$Shannon_mean) 
t.test(filter(ACE_CV_table, IRIS=="IS" & Bodysite=="Nasal")$Shannon_mean, filter(ACE_CV_table, IRIS=="IR" & Bodysite=="Nasal")$Shannon_mean) 
#test the shannon diversity difference (no significance)
t.test(filter(ACE_CV_table, IRIS=="IS" & Bodysite=="Stool")$Shannon_CV, filter(ACE_CV_table, IRIS=="IR" & Bodysite=="Stool")$Shannon_CV)
t.test(filter(ACE_CV_table, IRIS=="IS" & Bodysite=="Skin")$Shannon_CV, filter(ACE_CV_table, IRIS=="IR" & Bodysite=="Skin")$Shannon_CV)
t.test(filter(ACE_CV_table, IRIS=="IS" & Bodysite=="Oral")$Shannon_CV, filter(ACE_CV_table, IRIS=="IR" & Bodysite=="Oral")$Shannon_CV)
t.test(filter(ACE_CV_table, IRIS=="IS" & Bodysite=="Nasal")$Shannon_CV, filter(ACE_CV_table, IRIS=="IR" & Bodysite=="Nasal")$Shannon_CV)

#test the shannon diversity difference
t.test(filter(ACE_CV_table, IRIS=="IS" & Bodysite=="Stool")$AEC_mean, filter(ACE_CV_table, IRIS=="IR" & Bodysite=="Stool")$AEC_mean)
t.test(filter(ACE_CV_table, IRIS=="IS" & Bodysite=="Skin")$AEC_mean, filter(ACE_CV_table, IRIS=="IR" & Bodysite=="Skin")$AEC_mean)
t.test(filter(ACE_CV_table, IRIS=="IS" & Bodysite=="Oral")$AEC_mean, filter(ACE_CV_table, IRIS=="IR" & Bodysite=="Oral")$AEC_mean) 
t.test(filter(ACE_CV_table, IRIS=="IS" & Bodysite=="Nasal")$AEC_mean, filter(ACE_CV_table, IRIS=="IR" & Bodysite=="Nasal")$AEC_mean) 
#test the shannon diversity difference
t.test(filter(ACE_CV_table, IRIS=="IS" & Bodysite=="Stool")$Evenness, filter(ACE_CV_table, IRIS=="IR" & Bodysite=="Stool")$Evenness)
t.test(filter(ACE_CV_table, IRIS=="IS" & Bodysite=="Skin")$Evenness, filter(ACE_CV_table, IRIS=="IR" & Bodysite=="Skin")$Evenness)
t.test(filter(ACE_CV_table, IRIS=="IS" & Bodysite=="Oral")$Evenness, filter(ACE_CV_table, IRIS=="IR" & Bodysite=="Oral")$Evenness) 
t.test(filter(ACE_CV_table, IRIS=="IS" & Bodysite=="Nasal")$Evenness, filter(ACE_CV_table, IRIS=="IR" & Bodysite=="Nasal")$Evenness) 


filter(ACE_CV_table, IRIS=="IS" & Bodysite=="Stool")$Shannon_mean

lastingtime_table <- Diversity_metatable %>% dplyr::filter(!is.na(Diversity_metatable$CollectionDate))%>%  dplyr::select(SubjectID.x, CollectionDate, IRIS) %>% unique()  %>% group_by(SubjectID.x) %>% mutate(participatelength = max(CollectionDate) - min(CollectionDate)) %>% unique()
#make sure there are 86 individuals
unique(lastingtime_table$SubjectID.x)

lastingtime_table <- lastingtime_table %>% dplyr::select(-CollectionDate) %>% unique()
p.lasting <- ggplot(lastingtime_table, aes(x=IRIS, y= as.numeric(participatelength), fill=IRIS)) + geom_jitter()+ geom_boxplot(alpha=0.4, outlier.alpha = 0)
p.lasting <- p.lasting +scale_fill_manual(values = iris_color) + base_theme + ylab("Duration of the Observation (Days)") + xlab("")
p.lasting

t.test(filter(lastingtime_table, IRIS=="IS")$participatelength,filter(lastingtime_table, IRIS=="IR")$participatelength)

mean(lastingtime_table$participatelength)
sd(lastingtime_table$participatelength)

mean(filter(lastingtime_table, IRIS=="IS")$participatelength)
sd(filter(lastingtime_table, IRIS=="IS")$participatelength)

mean(filter(lastingtime_table, IRIS=="IR")$participatelength)
sd(filter(lastingtime_table, IRIS=="IR")$participatelength)


table(lastingtime_table$participatelength)

Diversity_metatable[,]


