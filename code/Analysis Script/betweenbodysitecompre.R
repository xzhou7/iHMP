#between microbiome dissimilarity 

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
meta.data <- read.csv("./Analysis/ori_meta_table/metadata.subject.csv")

#starting here
individual.list <- unique(dist$subject_id1)
distance_table <- data.frame()
for (i in 1:length(individual.list)){
  #get subject ID
  subject.id <- individual.list[i]
  print(subject.id)
  #filter each subject
  temp <- filter(dist, subject_id1 == subject.id & subject_id2 == subject.id)
  table(temp$class1 == temp$class2)
  #create a comparasion group
  temp$compare <- paste(temp$class1, temp$class2, sep="_")
  #calculate mean and sd
  temp2 <- temp %>% group_by(compare) %>% mutate(dist.mean = mean(dist)) %>% mutate(dist.sd =sd(dist)) %>% select(subject_id1, compare, dist.mean,dist.sd) %>% unique()
  #store data
  distance_table <- rbind(distance_table, temp2)
}

#combine with metadata
dt.meta <- merge(distance_table, select(meta.data, SubjectID, FPG_class,IRIS,Gender,Ethnicity,Adj.age, BMI, Th17.Group,OGTT_Class,A1C_Class), by.x="subject_id1", by.y="SubjectID")


p.iris <- ggplot(filter(dt.meta, IRIS != "Unknown"), aes(x=IRIS, y=dist.mean, fill=IRIS)) + geom_point(size=0.5) + geom_boxplot(alpha=0.8)+ base_theme
p.iris <- p.iris + facet_wrap(.~compare, scales = "free") + scale_fill_manual(values = iris_color) + ggtitle("Within individual, between bodysite compare (IRIS)")
p.iris 
ggsave("./Analysis/Suppl.figure/supplitoFig3.IRIS.betweenbodysite.distance.pdf",p.iris, height = 4, width = 6, dpi=300)

wilcox.test(filter(dt.meta, IRIS== "IS"&compare=="nasal_oral")$dist.mean, filter(dt.meta, IRIS== "IR"&compare=="nasal_oral")$dist.mean,alternative ="greater")
wilcox.test(filter(dt.meta, IRIS== "IS"&compare=="skin_nasal")$dist.mean, filter(dt.meta, IRIS== "IR"&compare=="skin_nasal")$dist.mean,alternative ="greater")
wilcox.test(filter(dt.meta, IRIS== "IS"&compare=="skin_oral")$dist.mean, filter(dt.meta, IRIS== "IR"&compare=="skin_oral")$dist.mean,alternative ="greater")
wilcox.test(filter(dt.meta, IRIS== "IS"&compare=="stool_nasal")$dist.mean, filter(dt.meta, IRIS== "IR"&compare=="stool_nasal")$dist.mean,alternative ="greater")
wilcox.test(filter(dt.meta, IRIS== "IS"&compare=="stool_oral")$dist.mean, filter(dt.meta, IRIS== "IR"&compare=="stool_oral")$dist.mean,alternative ="greater")
wilcox.test(filter(dt.meta, IRIS== "IS"&compare=="stool_skin")$dist.mean, filter(dt.meta, IRIS== "IR"&compare=="stool_skin")$dist.mean,alternative ="greater")


p.bmi <- ggplot(dt.meta, aes(x=BMI, y=dist.mean)) + geom_point(size=0.5) + base_theme + geom_smooth(method="lm")
p.bmi <- p.bmi + facet_wrap(.~compare, scales = "free", ncol=3) + ggtitle("Within individual, between bodysite compare (BMI)")
p.bmi 
ggsave("./Analysis/Suppl.figure/supplitoFig3.BMI.betweenbodysite.distance.pdf",p.bmi, height = 4, width = 6, dpi=300)

summary(lm(dist.mean ~ BMI, data=filter(dt.meta,compare== "nasal_oral")))
summary(lm(dist.mean ~ BMI, data=filter(dt.meta,compare== "skin_nasal")))
summary(lm(dist.mean ~ BMI, data=filter(dt.meta,compare== "skin_oral")))
summary(lm(dist.mean ~ BMI, data=filter(dt.meta,compare== "stool_nasal")))
summary(lm(dist.mean ~ BMI, data=filter(dt.meta,compare== "stool_oral")))
summary(lm(dist.mean ~ BMI, data=filter(dt.meta,compare== "stool_skin")))

#measure global effect
# summary(lm(dist.mean ~ FPG_class+IRIS+Gender+Ethnicity+Adj.age+ BMI+ Th17.Group+OGTT_Class+A1C_Class, data=filter(dt.meta,compare== "nasal_oral")))
# summary(lm(dist.mean ~ FPG_class+IRIS+Gender+Ethnicity+Adj.age+ BMI+ Th17.Group+OGTT_Class+A1C_Class, data=filter(dt.meta,compare== "skin_nasal")))
# summary(lm(dist.mean ~ FPG_class+IRIS+Gender+Ethnicity+Adj.age+ BMI+ Th17.Group+OGTT_Class+A1C_Class, data=filter(dt.meta,compare== "skin_oral")))
# summary(lm(dist.mean ~ FPG_class+IRIS+Gender+Ethnicity+Adj.age+ BMI+ Th17.Group+OGTT_Class+A1C_Class, data=filter(dt.meta,compare== "stool_nasal")))
# summary(lm(dist.mean ~ FPG_class+IRIS+Gender+Ethnicity+Adj.age+ BMI+ Th17.Group+OGTT_Class+A1C_Class, data=filter(dt.meta,compare== "stool_oral")))
# summary(lm(dist.mean ~ FPG_class+IRIS+Gender+Ethnicity+Adj.age+ BMI+ Th17.Group+OGTT_Class+A1C_Class, data=filter(dt.meta,compare== "stool_skin")))

