#Distance Over Time
#This is to calculate BC distance over time
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

load("./Analysis/Robject/Ordination_Bray.RData")
#load("./Analysis/Robject/PhyloseqObject.RData")
load("./Analysis/Robject/DetailedPhyloseq.RData")
load("./Analysis/Robject/Revision_MultiOmes_0509.RData")

HMP.Date <- read.csv("./Analysis/ori_meta_table/DateInfo.csv", header = T, row.names = 1)
HMP.Date$Date <- as.Date(HMP.Date$Date, format = "%m/%d/%y")
UMAP.coordi <- read.csv(file = "./Analysis/ori_meta_table/UMAP.coordi.csv",header = T, row.names = 1)

# #this is distance by genus
# load("~/Library/CloudStorage/Box-Box/human_microbiome_project/data_analysis/combine_microbiome/distance/stool/stool_braydist_by_genus")
# stool_braydist_by_genus_within <- filter(stool_braydist_by_genus, type!="between")
# stool_braydist_by_genus_within
# rm(stool_braydist_by_genus)
# load("~/Library/CloudStorage/Box-Box/human_microbiome_project/data_analysis/combine_microbiome/distance/skin/skin_braydist_by_genus")
# skin_braydist_by_genus_within <- filter(skin_braydist_by_genus, type!="between")
# skin_braydist_by_genus_within
# rm(skin_braydist_by_genus)
# load("~/Library/CloudStorage/Box-Box/human_microbiome_project/data_analysis/combine_microbiome/distance/oral/oral_braydist_by_genus")
# oral_braydist_by_genus_within <- filter(oral_braydist_by_genus, type!="between")
# oral_braydist_by_genus_within
# rm(oral_braydist_by_genus)
# load("~/Library/CloudStorage/Box-Box/human_microbiome_project/data_analysis/combine_microbiome/distance/nasal/nasal_braydist_by_genus")
# nasal_braydist_by_genus_within <- filter(nasal_braydist_by_genus, type!="between")
# nasal_braydist_by_genus_within
# rm(nasal_braydist_by_genus)
# 
# table(stool_braydist_by_genus_within$family2)
# table(skin_braydist_by_genus_within$family2)
# table(oral_braydist_by_genus_within$family2)
# table(nasal_braydist_by_genus_within$family2)
# stool_braydist_by_sample

table(temp$subject_id1)

table(stool_braydist_by_genus_within$genus) %>% data.frame() %>% arrange(desc(Freq))

temp <- filter(stool_braydist_by_genus_within,genus == "Unclassified_Ruminococcaceae")
coefplot2(lmer(dist ~ diffdays + (1|subject_id1), temp))

# physeq_C_ord$vectors[1:5,1:5]
# dim(physeq_C_ord$vectors)
# plot(data.frame(physeq_C_ord$vectors)[,1],data.frame(physeq_C_ord$vectors)[,2])

min(HMP.Date$Date[!is.na(HMP.Date$Date)])
HMP.Date <- HMP.Date %>% filter(!is.na(HMP.Date$Date)) %>% mutate(daysfrombegining = difftime(Date , "2010-03-28", units="weeks"))
HMP.Date <- HMP.Date %>% group_by(Subject) %>% mutate(order.subject = rank(daysfrombegining))

ST.sample <- data.frame(sample_data(physeq_ST)) %>% select(RandomID,SampleID,SubjectID,SampleType) %>% mutate(SampleType = "Stool")
SK.sample <- data.frame(sample_data(physeq_SK)) %>% select(KitID,SampleID,SubjectID,SampleType)
OR.sample <- data.frame(sample_data(physeq_OR)) %>% select(KitID,SampleID,SubjectID,SampleType)
NS.sample <- data.frame(sample_data(physeq_NS)) %>% select(RandomID,SampleID,SubjectID,SampleType) %>% mutate(SampleType = "Nasal")

colnames(SK.sample)[1] <- "RandomID"
colnames(OR.sample)[1] <- "RandomID"

meta.sample <- rbind(ST.sample,SK.sample,OR.sample,NS.sample)
meta.sample$RandomID

dim(UMAP.coordi)
UMAP.dist <- dist(UMAP.coordi)
UMAP.dist.long <- melt(as.matrix(UMAP.dist)) 

dim(UMAP.dist.long)

colnames(UMAP.dist.long) <- c("sample.A", "sample.B", "E.distance")

UMAP.dist.long$bodysite.A <- "Stool"
UMAP.dist.long$bodysite.B <- "Stool"
UMAP.dist.long$bodysite.A[UMAP.dist.long$sample.A %in% SK.sample$RandomID] <- "Skin"
UMAP.dist.long$bodysite.A[UMAP.dist.long$sample.A %in% OR.sample$RandomID] <- "Oral"
UMAP.dist.long$bodysite.A[UMAP.dist.long$sample.A %in% NS.sample$RandomID] <- "Nasal"
UMAP.dist.long$bodysite.B[UMAP.dist.long$sample.B %in% SK.sample$RandomID] <- "Skin"
UMAP.dist.long$bodysite.B[UMAP.dist.long$sample.B %in% OR.sample$RandomID] <- "Oral"
UMAP.dist.long$bodysite.B[UMAP.dist.long$sample.B %in% NS.sample$RandomID] <- "Nasal"

UMAP.dist.long$SampleID.A <- meta.sample$SampleID[match(UMAP.dist.long$sample.A,meta.sample$RandomID)]
UMAP.dist.long$SampleID.B <- meta.sample$SampleID[match(UMAP.dist.long$sample.B,meta.sample$RandomID)]

UMAP.dist.long$Subject.A <- meta.sample$SubjectID[match(UMAP.dist.long$sample.A,meta.sample$RandomID)]
UMAP.dist.long$Subject.B <- meta.sample$SubjectID[match(UMAP.dist.long$sample.B,meta.sample$RandomID)]

UMAP.dist.long$Week.A <- HMP.Date$daysfrombegining[match(UMAP.dist.long$SampleID.A,HMP.Date$SampleID)]
UMAP.dist.long$Week.B <- HMP.Date$daysfrombegining[match(UMAP.dist.long$SampleID.B,HMP.Date$SampleID)]

UMAP.dist.long$week.diff <- abs(UMAP.dist.long$Week.A - UMAP.dist.long$Week.B) %>% as.numeric()
UMAP.dist.long

table(str_detect(UMAP.dist.long$bodysite.A, "Stool"))
table(str_detect(UMAP.dist.long$bodysite.A, "Skin"))
table(str_detect(UMAP.dist.long$bodysite.A, "Oral"))
table(str_detect(UMAP.dist.long$bodysite.A, "Nasal"))

#remove lower half of the distance matrix
df <- data.frame(V1 = UMAP.dist.long$sample.A, V2=UMAP.dist.long$sample.B)
#Warning: take 8 minutes to run next step
df <- data.frame(t(apply(df, 1, sort)))
#select duplicated would remove auto-pair
UMAP.dist.long.upper <- UMAP.dist.long[duplicated(df),]
#proved here:
table(UMAP.dist.long$sample.A == "GAOICDYI" & UMAP.dist.long$sample.B == "GAOICDYI")
table(UMAP.dist.long.upper$sample.A == "GAOICDYI" & UMAP.dist.long.upper$sample.B == "GAOICDYI")

table(UMAP.dist.long.upper$bodysite.A)
table(UMAP.dist.long.upper$bodysite.B)

UMAP.dist.long.upper$type <- "Between"
UMAP.dist.long.upper$type[UMAP.dist.long.upper$Subject.A == UMAP.dist.long.upper$Subject.B] <- "Within"

p.dist.wb <- UMAP.dist.long.upper[UMAP.dist.long.upper$bodysite.A == UMAP.dist.long.upper$bodysite.B,] %>% ggplot(aes(x=type, y= E.distance, group=type, color=type)) + geom_boxplot()
p.dist.wb <- p.dist.wb + facet_wrap(.~bodysite.A) + base_theme + ggtitle("Distance at Genus Level, Between/Within")
p.dist.wb 

#Get a list if need to remove infection time point
temp.list <- ls[!str_detect(ls$CL4, "Infection") &
#                  str_detect(ls$SubjectID, "69-063") &
                  !str_detect(ls$CL4, "Imz")&
                  !str_detect(ls$CL4, "Ant"),] %>% select(SampleID)

temp.list

p.dist.wt.data <- UMAP.dist.long.upper[UMAP.dist.long.upper$bodysite.A == UMAP.dist.long.upper$bodysite.B,] %>% 
   filter(type == "Within")       #%>% 
#  filter(bodysite.A == "Stool")  #%>% 
#  filter(Subject.A =="69-063")   #%>% 
#  filter(SampleID.A %in% temp.list$SampleID) %>% filter(SampleID.B %in% temp.list$SampleID)

#Insert Infection Annotation to the plot
p.dist.wt.data$CL4 <- "Healthy"
p.dist.wt.data$CL4[p.dist.wt.data$SampleID.A %in% ls[str_detect(ls$CL4, "Infection"),]$SampleID] <- "Infection"
p.dist.wt.data$CL4[p.dist.wt.data$SampleID.B %in% ls[str_detect(ls$CL4, "Infection"),]$SampleID] <- "Infection"
p.dist.wt.data$CL4[p.dist.wt.data$SampleID.A %in% ls[str_detect(ls$CL4, "Imz"),]$SampleID] <- "Immunization"
p.dist.wt.data$CL4[p.dist.wt.data$SampleID.B %in% ls[str_detect(ls$CL4, "Imz"),]$SampleID] <- "Immunization"
p.dist.wt.data$CL4[p.dist.wt.data$SampleID.A %in% ls[str_detect(ls$CL4, "Ant"),]$SampleID] <- "Antibiotics"
p.dist.wt.data$CL4[p.dist.wt.data$SampleID.B %in% ls[str_detect(ls$CL4, "Ant"),]$SampleID] <- "Antibiotics"
table(p.dist.wt.data$CL4)
p.dist.wt.data$CL4 <- factor(p.dist.wt.data$CL4, levels = c('Antibiotics',"Immunization","Infection","Healthy"))

p.dist.wt.data$IRIS <- "Unknown"
p.dist.wt.data$IRIS[p.dist.wt.data$Subject.A %in% filter(sc, IRIS == "IR")$SubjectID] <- "IR"
p.dist.wt.data$IRIS[p.dist.wt.data$Subject.A %in% filter(sc, IRIS == "IS")$SubjectID] <- "IS"
table(p.dist.wt.data$IRIS)


p.dist.wt <- ggplot(p.dist.wt.data, aes(x=week.diff, y= E.distance)) + geom_point()
p.dist.wt <- p.dist.wt + facet_wrap(. ~bodysite.A, scales = "free_x") + base_theme + ggtitle("Within Distance by time")
#p.dist.wt <- p.dist.wt + geom_smooth(method="loess", color="blue")
p.dist.wt <- p.dist.wt + geom_smooth(method="lm", color="red")
p.dist.wt

#compare time related difference in four body sites
lmer.dist.wt <- lmer(E.distance ~  week.diff:bodysite.A + (week.diff|Subject.A),REML=T,data = p.dist.wt.data)
#lmer.dist.wt <- lmer(E.distance ~ week.diff + (1|Subject.A),REML=F,data = p.dist.wt.data %>% filter(bodysite.A =="Stool"))
summ(lmer.dist.wt, digit=4)
coefplot2(lmer.dist.wt)

lmer.dist.wt.iris <- lmer(E.distance ~  week.diff:IRIS:bodysite.A + (week.diff|Subject.A),REML=T,data = p.dist.wt.data %>% filter(IRIS != "Unknown"))
#lmer.dist.wt <- lmer(E.distance ~ week.diff + (1|Subject.A),REML=F,data = p.dist.wt.data %>% filter(bodysite.A =="Stool"))
summ(lmer.dist.wt.iris, digit=4)
coefplot2(lmer.dist.wt.iris)

par(mfrow=c(1,2))
coefplot2(lmer.dist.wt)
coefplot2(lmer.dist.wt.iris)

lmer.dist.wt.if <- lmer(E.distance ~  week.diff:bodysite.A:CL4 + (CL4|Subject.A),REML=T,data = p.dist.wt.data)
#lmer.dist.wt <- lmer(E.distance ~ week.diff + (1|Subject.A),REML=F,data = p.dist.wt.data %>% filter(bodysite.A =="Stool"))
summ(lmer.dist.wt.if, digit=4)
coefplot2(lmer.dist.wt.if)

lmer.dist.wt.if2 <- lmer(E.distance ~  week.diff:CL4 + (1|Subject.A), REML=T,data = p.dist.wt.data %>% filter(bodysite.A =="Stool"))
summ(lmer.dist.wt.if2, digit=4)
coefplot2(lmer.dist.wt.if2)

dev.off()


par(mfrow=c(2,1))
coefplot2(lmer.dist.wt.h, main="healthy time point only", xlim=c(-0.005, 0.035))
coefplot2(lmer.dist.wt, main="all time point", xlim=c(-0.005, 0.035))

# #compare time related difference in infection status by bodysite
# par(mfrow=c(2,1))
# 
# lmer.dist.wt.st.h <- lmer(E.distance ~  week.diff + (week.diff|Subject.A),REML=T,data = p.dist.wt.data  %>% filter(bodysite.A =="Stool") %>% filter(CL4 == "Healthy"))
# summ(lmer.dist.wt.st.h, digit=4)
# coefplot2(lmer.dist.wt.st.h, main="Stool_Healthy", xlim=c(-0.02,0.02))
# 
# lmer.dist.wt.st <- lmer(E.distance ~  week.diff + (week.diff|Subject.A),REML=T,data = p.dist.wt.data  %>% filter(bodysite.A =="Stool"))
# summ(lmer.dist.wt.st, digit=4)
# coefplot2(lmer.dist.wt.st,main="Stool_All", xlim=c(-0.02,0.02))
# 
# par(mfrow=c(2,1))
# 
# lmer.dist.wt.sk.h <- lmer(E.distance ~  week.diff + (week.diff|Subject.A),REML=F,data = p.dist.wt.data  %>% filter(bodysite.A =="Skin") %>% filter(CL4 == "Healthy"))
# summ(lmer.dist.wt.sk.h, digit=4)
# coefplot2(lmer.dist.wt.sk.h,main="Skin_Healthy", xlim=c(-0.02,0.02))
# 
# lmer.dist.wt.sk <- lmer(E.distance ~  week.diff + (week.diff|Subject.A),REML=F,data = p.dist.wt.data  %>% filter(bodysite.A =="Skin"))
# summ(lmer.dist.wt.sk, digit=4)
# coefplot2(lmer.dist.wt.sk,main="Skin_All", xlim=c(-0.02,0.02))
# 
# par(mfrow=c(2,1))
# 
# lmer.dist.wt.or.h <- lmer(E.distance ~  week.diff + (week.diff|Subject.A),REML=T,data = p.dist.wt.data  %>% filter(bodysite.A =="Oral") %>% filter(CL4 == "Healthy"))
# summ(lmer.dist.wt.or, digit=4)
# coefplot2(lmer.dist.wt.or)
# 
# lmer.dist.wt.or <- lmer(E.distance ~  week.diff:CL4 + (week.diff|Subject.A),REML=T,data = p.dist.wt.data  %>% filter(bodysite.A =="Oral"))
# summ(lmer.dist.wt.or, digit=4)
# coefplot2(lmer.dist.wt.or)
# 
# lmer.dist.wt.ns <- lmer(E.distance ~  week.diff:CL4 + (week.diff|Subject.A),REML=T,data = p.dist.wt.data  %>% filter(bodysite.A =="Nasal"))
# summ(lmer.dist.wt.ns, digit=4)
# coefplot2(lmer.dist.wt.ns)

par(mfrow=c(2,2))
coefplot2(lmer.dist.wt.st,main="Stool")
coefplot2(lmer.dist.wt.sk,main="Skin")
coefplot2(lmer.dist.wt.or,main="Oral")
coefplot2(lmer.dist.wt.ns,main="Nasal")

summ(lm.dist.wt.st)

#Check the beta-coefficient increase with time difference 
lm.dist.wt.st <- lmer(E.distance ~ week.diff:bodysite.A + (1|Subject.A),REML=F, data = p.dist.wt.data)# %>% filter(bodysite.A =="Stool"))
lm.dist.wt.st
summ(lm.dist.wt.st, digits=4)
#export_summs(lm.dist.wt.st.3)
#effect_plot(lm.dist.wt.st.3)
#simple_slopes(lm.dist.wt.st.3)
#graph_model(lm.dist.wt.st.3, y=E.distance, x=week.diff,lines=CL4)
coefplot2(lm.dist.wt.st)

lm.dist.wt.st <- lmer(E.distance ~ week.diff + (1|Subject.A), data = p.dist.wt.data %>% filter(bodysite.A =="Stool"))

lm.dist.wt.st.cl4 <- lmer(E.distance ~  week.diff * CL4 + (1|Subject.A), data = p.dist.wt.data %>% filter(bodysite.A =="Stool"))
summary(lm.dist.wt.st.cl4)

anova(lm.dist.wt.st, lm.dist.wt.st.cl4)


lm.dist.wt <- lmer(E.distance ~ CL4 + bodysite.A:CL4 + (CL4|Subject.A), data = p.dist.wt.data) #%>% filter(bodysite.A =="Stool"))
lm.dist.wt
summary(lm.dist.wt)
coefplot2(lm.dist.wt)


p.dist.wt.CL4 <- ggplot(p.dist.wt.data, aes(x=CL4, y= E.distance)) + geom_boxplot()
p.dist.wt.CL4 <- p.dist.wt.CL4 + facet_wrap(. ~bodysite.A, scales = "free_x") + base_theme + ggtitle("Within Distance by Infection")
#p.dist.wt <- p.dist.wt + geom_smooth(method="loess", color="blue")
p.dist.wt.CL4 <- p.dist.wt.CL4 + geom_smooth(method="lm", color="red")
p.dist.wt.CL4

coefplot2(lm.dist.wt)


longdis <- UMAP.dist.long.upper %>% filter(E.distance > 15 & type == "Within")



UMAP.dist.long.upper[UMAP.dist.long.upper$bodysite.A == UMAP.dist.long.upper$bodysite.B,] %>% filter(type == "Within")

Dist.Test <- UMAP.dist.long.upper[UMAP.dist.long.upper$bodysite.A == UMAP.dist.long.upper$bodysite.B,] %>% 
  filter(type == "Within") %>%
  select(Subject.A,E.distance,week.diff) %>%
  nest(-Subject.A)%>% 
  mutate(cor=map(data,~cor.test(.x$E.distance, .x$week.diff, method = "pearson"))) %>%
  mutate(tidied = map(cor, tidy)) %>% 
  unnest(tidied, .drop = T) %>% mutate(p.FDR = p.adjust(p.value, method = "BH")) %>% 
  select(-data, -cor)

###################perform test body site wise########################
stool.list <- UMAP.dist.long.upper[UMAP.dist.long.upper$bodysite.A == UMAP.dist.long.upper$bodysite.B,] %>% 
  filter(type == "Within") %>%
  filter(bodysite.A == "Stool") %>% 
#  filter(SampleID.A %in% temp.list$SampleID) %>% filter(SampleID.B %in% temp.list$SampleID) %>% 
  select(Subject.A) %>%
  table() %>% 
  data.frame() %>% rename(Subject = ".") %>%
  filter(Freq > 5) %>% select(Subject)

skin.list <- UMAP.dist.long.upper[UMAP.dist.long.upper$bodysite.A == UMAP.dist.long.upper$bodysite.B,] %>% 
  filter(type == "Within") %>%
  filter(bodysite.A == "Skin") %>% 
#  filter(SampleID.A %in% temp.list$SampleID) %>% filter(SampleID.B %in% temp.list$SampleID) %>% 
  select(Subject.A) %>%
  table() %>% 
  data.frame() %>% rename(Subject = ".") %>%
  filter(Freq > 5) %>% select(Subject)

oral.list <- UMAP.dist.long.upper[UMAP.dist.long.upper$bodysite.A == UMAP.dist.long.upper$bodysite.B,] %>% 
  filter(type == "Within") %>%
  filter(bodysite.A == "Oral") %>% 
#  filter(SampleID.A %in% temp.list$SampleID) %>% filter(SampleID.B %in% temp.list$SampleID) %>% 
  select(Subject.A) %>%
  table() %>% 
  data.frame() %>% rename(Subject = ".") %>%
  filter(Freq > 5) %>% select(Subject)

nasal.list <- UMAP.dist.long.upper[UMAP.dist.long.upper$bodysite.A == UMAP.dist.long.upper$bodysite.B,] %>% 
  filter(type == "Within") %>%
  filter(bodysite.A == "Nasal") %>% 
#  filter(SampleID.A %in% temp.list$SampleID) %>% filter(SampleID.B %in% temp.list$SampleID) %>% 
  select(Subject.A) %>%
  table() %>% 
  data.frame() %>% rename(Subject = ".") %>%
  filter(Freq > 5) %>% select(Subject)


#looking for each bodysite, stool N=63
Dist.Test.ST <- UMAP.dist.long.upper[UMAP.dist.long.upper$bodysite.A == UMAP.dist.long.upper$bodysite.B,] %>% 
  filter(type == "Within") %>%
  filter(bodysite.A == "Stool") %>%
  filter(Subject.A %in% stool.list$Subject) %>%
#  filter(SampleID.A %in% temp.list$SampleID) %>% filter(SampleID.B %in% temp.list$SampleID) %>% 
  select(Subject.A,E.distance,week.diff) %>%
  nest(-Subject.A)%>% 
  mutate(cor=map(data,~cor.test(.x$E.distance, .x$week.diff, method = "pearson"))) %>%
  mutate(tidied = map(cor, tidy)) %>% 
  unnest(tidied, .drop = T) %>% 
  select(-data, -cor) %>% mutate(p.FDR = p.adjust(p.value, method = "BH")) %>% mutate(Bodysite= "Stool")

#N=76
Dist.Test.SK <- UMAP.dist.long.upper[UMAP.dist.long.upper$bodysite.A == UMAP.dist.long.upper$bodysite.B,] %>% 
  filter(type == "Within") %>%
  filter(bodysite.A == "Skin") %>%
  filter(Subject.A %in% skin.list$Subject) %>%
#  filter(SampleID.A %in% temp.list$SampleID) %>% filter(SampleID.B %in% temp.list$SampleID) %>% 
  select(Subject.A,E.distance,week.diff) %>%
  nest(-Subject.A)%>% 
  mutate(cor=map(data,~cor.test(.x$E.distance, .x$week.diff, method = "pearson"))) %>%
  mutate(tidied = map(cor, tidy)) %>% 
  unnest(tidied, .drop = T) %>% 
  select(-data, -cor) %>% mutate(p.FDR = p.adjust(p.value, method = "BH"))%>% mutate(Bodysite= "Skin")

#N=75
Dist.Test.OR <- UMAP.dist.long.upper[UMAP.dist.long.upper$bodysite.A == UMAP.dist.long.upper$bodysite.B,] %>% 
  filter(type == "Within") %>%
  filter(bodysite.A == "Oral") %>%
  filter(Subject.A %in% oral.list$Subject) %>%
#  filter(SampleID.A %in% temp.list$SampleID) %>% filter(SampleID.B %in% temp.list$SampleID) %>% 
  select(Subject.A,E.distance,week.diff) %>%
  nest(-Subject.A)%>% 
  mutate(cor=map(data,~cor.test(.x$E.distance, .x$week.diff, method = "pearson"))) %>%
  mutate(tidied = map(cor, tidy)) %>% 
  unnest(tidied, .drop = T) %>% 
  select(-data, -cor) %>% mutate(p.FDR = p.adjust(p.value, method = "BH"))%>% mutate(Bodysite= "Oral")

#N=65
Dist.Test.NS <- UMAP.dist.long.upper[UMAP.dist.long.upper$bodysite.A == UMAP.dist.long.upper$bodysite.B,] %>% 
  filter(type == "Within") %>%
  filter(bodysite.A == "Nasal") %>%
  filter(Subject.A %in% nasal.list$Subject) %>%
  #add or not add infection time point
#  filter(SampleID.A %in% temp.list$SampleID) %>% filter(SampleID.B %in% temp.list$SampleID) %>% 
  select(Subject.A,E.distance,week.diff) %>%
  nest(-Subject.A)%>% 
  mutate(cor=map(data,~cor.test(.x$E.distance, .x$week.diff, method = "pearson"))) %>%
  mutate(tidied = map(cor, tidy)) %>% 
  unnest(tidied, .drop = T) %>% 
  select(-data, -cor) %>% mutate(p.FDR = p.adjust(p.value, method = "BH"))%>% mutate(Bodysite= "Nasal")

Dist.Test.Bodysite <- rbind(Dist.Test.ST,Dist.Test.SK,Dist.Test.OR,Dist.Test.NS)
#write.csv(file = "./Analysis/Distance_By_Time/Distance.Test.csv",Dist.Test.Bodysite)

p.comb <- ggplot(Dist.Test, aes(p.FDR,estimate)) + geom_point() + geom_vline(xintercept = 0.05) + base_theme + ggtitle("correlation combined body sites (N=86)")
p.comb
#ggsave(filename = "./Analysis/Suppl.figure/3.BDbyTime.PBeta.pdf", p.comb, width = 4, height = 3, dpi = 300)

p.comb.body <- ggplot(Dist.Test.Bodysite, aes(p.FDR, estimate, color=Bodysite)) + geom_point(size=0.5) + geom_vline(xintercept = c(0.05, 0.2)) + base_theme + scale_color_manual(values=body_site_color)
p.comb.body
#ggsave(filename = "./Analysis/Suppl.figure/3.BDbyTime.PBeta.bybodtsite.noinfection.pdf", p.comb.body, width = 5, height = 3, dpi = 300)

Dist.Test.Bodysite %>% filter(Subject.A %in% c("69-063", "70-1014"))


Dist.Test.Bodysite.sig <- filter(Dist.Test.Bodysite,p.FDR < 0.2)

ggplot(Dist.Test.Bodysite, aes(x=Bodysite, y=estimate, fill=Bodysite)) + geom_boxplot() + 
  stat_compare_means() + base_theme + scale_fill_manual(values=body_site_color) + ggtitle("Beta estimate across bodysite")


ggplot(Dist.Test.Bodysite.sig, aes(x=Bodysite, y=estimate, fill=Bodysite)) + geom_boxplot() + 
  stat_compare_means() + base_theme + scale_fill_manual(values=body_site_color) + ggtitle("FDR0.2 Beta estimate across bodysite")


##########Here we will insert an infection per period concept on this map

##############################################
# Cross body-sites test
##############################################
colnames(Dist.Test.ST) <- paste(colnames(Dist.Test.ST), "ST", sep=".")
colnames(Dist.Test.SK) <- paste(colnames(Dist.Test.SK), "SK", sep=".")
colnames(Dist.Test.OR) <- paste(colnames(Dist.Test.OR), "OR", sep=".")
colnames(Dist.Test.NS) <- paste(colnames(Dist.Test.NS), "NS", sep=".")

pcutoff <- 0.2

ST.SK.Compare <- merge(Dist.Test.ST,Dist.Test.SK, by.x = "Subject.A.ST", by.y = "Subject.A.SK") %>% filter(p.FDR.ST < pcutoff & p.FDR.SK < pcutoff)
ST.OR.Compare <- merge(Dist.Test.ST,Dist.Test.OR, by.x = "Subject.A.ST", by.y = "Subject.A.OR") %>% filter(p.FDR.ST < pcutoff & p.FDR.OR < pcutoff)
ST.NS.Compare <- merge(Dist.Test.ST,Dist.Test.NS, by.x = "Subject.A.ST", by.y = "Subject.A.NS") %>% filter(p.FDR.ST < pcutoff & p.FDR.NS < pcutoff)
SK.OR.Compare <- merge(Dist.Test.SK, Dist.Test.OR, by.x = "Subject.A.SK", by.y = "Subject.A.OR") %>% filter(p.FDR.SK < pcutoff & p.FDR.OR < pcutoff)
SK.NS.Compare <- merge(Dist.Test.SK, Dist.Test.NS, by.x = "Subject.A.SK", by.y = "Subject.A.NS") %>% filter(p.FDR.SK < pcutoff & p.FDR.NS < pcutoff)
OR.NS.Compare <- merge(Dist.Test.OR, Dist.Test.NS, by.x = "Subject.A.OR", by.y = "Subject.A.NS") %>% filter(p.FDR.OR < pcutoff & p.FDR.NS < pcutoff)

temp <- cor.test(ST.SK.Compare$estimate.ST,ST.SK.Compare$estimate.SK,method = "pearson")
c(temp$p.value,temp$estimate)

cor.test(ST.OR.Compare$estimate.ST,ST.OR.Compare$estimate.OR,method = "pearson")
cor.test(ST.NS.Compare$estimate.ST,ST.NS.Compare$estimate.NS,method = "pearson")
cor.test(SK.OR.Compare$estimate.SK,SK.OR.Compare$estimate.OR,method = "pearson")
cor.test(SK.NS.Compare$estimate.SK,SK.NS.Compare$estimate.NS,method = "pearson")
cor.test(OR.NS.Compare$estimate.OR,OR.NS.Compare$estimate.NS,method = "pearson")

tableGrob(cor.test(ST.SK.Compare$estimate.ST,ST.SK.Compare$estimate.SK,method = "pearson"))

p1 <- ggscatter(ST.SK.Compare,x= "estimate.ST",y= "estimate.SK",add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), conf.int = TRUE) +  stat_cor(method = "pearson", label.x = 0, label.y = 1)
p2 <- ggscatter(ST.OR.Compare,x= "estimate.ST",y= "estimate.OR",add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), conf.int = TRUE) +  stat_cor(method = "pearson", label.x = 0, label.y = 1)
p3 <- ggscatter(ST.NS.Compare,x= "estimate.ST",y= "estimate.NS",add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), conf.int = TRUE) +  stat_cor(method = "pearson", label.x = 0, label.y = 2)
p4 <- ggscatter(SK.OR.Compare,x= "estimate.SK",y= "estimate.OR",add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), conf.int = TRUE) +  stat_cor(method = "pearson", label.x = 0.15, label.y = 1)
p5 <- ggscatter(SK.NS.Compare,x= "estimate.SK",y= "estimate.NS",add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), conf.int = TRUE) +  stat_cor(method = "pearson", label.x = 0.15, label.y = 1)
p6 <- ggscatter(OR.NS.Compare,x= "estimate.OR",y= "estimate.NS",add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), conf.int = TRUE) +  stat_cor(method = "pearson", label.x = 0.15, label.y = 2)

p1 + p2 + p3 + p4 + p5 + p6

#test the correlation for "within distance"
UMAP.dist.long.upper.within <- filter(UMAP.dist.long.upper, type=="Within")

unique(UMAP.dist.long.upper.within$Subject.A)

temp <- filter(UMAP.dist.long.upper.within, Subject.A == "69-001")
temp <- UMAP.dist.long.upper.within

temp.stool <- filter(temp, bodysite.A == "Stool" & bodysite.B=="Stool")
temp.skin <- filter(temp, bodysite.A == "Skin" & bodysite.B=="Skin")
temp.oral <- filter(temp, bodysite.A == "Oral" & bodysite.B=="Oral")
temp.nasal <- filter(temp, bodysite.A == "Nasal" & bodysite.B=="Nasal")

colnames(temp.stool) <- paste(colnames(temp.stool), "stool", sep=".")
colnames(temp.skin) <- paste(colnames(temp.skin), "skin", sep=".")
colnames(temp.oral) <- paste(colnames(temp.oral), "oral", sep=".")
colnames(temp.nasal) <- paste(colnames(temp.nasal), "nasal", sep=".")

temp.ST.SK <- merge(temp.stool,temp.skin, by.x = c("SampleID.A.stool", "SampleID.B.stool"), by.y=c("SampleID.A.skin","SampleID.B.skin"))
mixed.lm.sk.st <- lmer(E.distance.stool ~ E.distance.skin + (1|Subject.A.stool), data = temp.ST.SK)
summary(mixed.lm.sk.st)
anova(mixed.lm.sk.st)
plot(temp.ST.SK$E.distance.stool,temp.ST.SK$E.distance.skin)
cor.test(temp.ST.SK$E.distance.stool,temp.ST.SK$E.distance.skin)
identical(temp.ST.SK$week.diff.stool,temp.ST.SK$week.diff.skin)

temp.ST.OR <- merge(temp.stool,temp.oral, by.x = c("SampleID.A.stool", "SampleID.B.stool"), by.y=c("SampleID.A.oral","SampleID.B.oral"))
mixed.lm.st.or <- lmer(E.distance.stool ~ E.distance.oral + (1|Subject.A.stool), data = temp.ST.OR)
summary(mixed.lm.st.or)
anova(mixed.lm.st.or)
cor.test(temp.ST.OR$E.distance.stool,temp.ST.OR$E.distance.oral)
identical(temp.ST.SK$week.diff.stool,temp.ST.SK$week.diff.skin)

temp.ST.NS <- merge(temp.stool,temp.nasal, by.x = c("SampleID.A.stool", "SampleID.B.stool"), by.y=c("SampleID.A.nasal","SampleID.B.nasal"))
mixed.lm.st.ns <- lmer(E.distance.stool ~ E.distance.nasal + (1|Subject.A.stool), data = temp.ST.NS)
summary(mixed.lm.st.ns)
anova(mixed.lm.st.ns)
plot(temp.ST.NS$E.distance.stool,temp.ST.NS$E.distance.nasal)
cor.test(temp.ST.NS$E.distance.stool,temp.ST.NS$E.distance.nasal)
identical(temp.ST.SK$week.diff.stool,temp.ST.SK$week.diff.skin)

temp.SK.OR <- merge(temp.skin,temp.oral, by.x = c("SampleID.A.skin", "SampleID.B.skin"), by.y=c("SampleID.A.oral","SampleID.B.oral"))
mixed.lm.sk.or <- lmer(E.distance.skin ~  E.distance.oral + (1|Subject.A.skin), data = temp.SK.OR)
summary(mixed.lm.sk.or)
anova(mixed.lm.sk.or)
cor.test(temp.SK.OR$E.distance.skin,temp.SK.OR$E.distance.oral)
identical(temp.SK.OR$week.diff.skin,temp.SK.OR$week.diff.oral)

temp.SK.NS <- merge(temp.skin, temp.nasal, by.x=c("SampleID.A.skin","SampleID.B.skin"), by.y = c("SampleID.A.nasal", "SampleID.B.nasal"))
mixed.lm.sk.ns <- lmer(E.distance.skin ~ E.distance.nasal + (1|Subject.A.skin), data = temp.SK.NS)
summary(mixed.lm.sk.ns)
anova(mixed.lm.sk.ns)
cor.test(temp.SK.NS$E.distance.nasal,temp.SK.NS$E.distance.skin)
plot(temp.SK.NS$E.distance.nasal,temp.SK.NS$E.distance.skin)

temp.OR.NS <- merge(temp.oral, temp.nasal, by.x=c("SampleID.A.oral","SampleID.B.oral"), by.y = c("SampleID.A.nasal", "SampleID.B.nasal"))
mixed.lm.or.ns <- lmer(E.distance.oral ~ E.distance.nasal + (1|Subject.A.oral), data = temp.OR.NS)
summary(mixed.lm.or.ns)
anova(mixed.lm.or.ns)
cor.test(temp.OR.NS$E.distance.oral,temp.OR.NS$E.distance.nasal)
plot(temp.OR.NS$E.distance.oral,temp.OR.NS$E.distance.nasal)

anova(mixed.lm.sk.st)
anova(mixed.lm.st.or)
anova(mixed.lm.st.ns)
anova(mixed.lm.sk.or)
anova(mixed.lm.sk.ns)
anova(mixed.lm.or.ns)

#plot coefficient
coefplot2(mixed.lm.st.ns,xlim=c(-0.1, 0.2))
coefplot2(mixed.lm.sk.st,xlim=c(-0.1, 0.2), add=T, offset= 0.6)
coefplot2(mixed.lm.st.or,xlim=c(-0.1, 0.2), add=T, offset= 0.3)
coefplot2(mixed.lm.sk.or,xlim=c(-0.1, 0.2), add=T, offset= - 0.3)
coefplot2(mixed.lm.sk.ns,xlim=c(-0.1, 0.2), add=T, offset= - 0.6)
coefplot2(mixed.lm.or.ns,xlim=c(-0.1, 0.2), add=T, offset= - 0.9)






############this is a calculation based on BC distance, not ED from UMAP
load("~/Library/CloudStorage/Box-Box/human_microbiome_project/data_analysis/combine_microbiome/within_subject_distance_by_sample/stool/stool_braydist_by_sample")
load("~/Library/CloudStorage/Box-Box/human_microbiome_project/data_analysis/combine_microbiome/within_subject_distance_by_sample/skin/skin_braydist_by_sample")
load("~/Library/CloudStorage/Box-Box/human_microbiome_project/data_analysis/combine_microbiome/within_subject_distance_by_sample/oral/oral_braydist_by_sample")
load("~/Library/CloudStorage/Box-Box/human_microbiome_project/data_analysis/combine_microbiome/within_subject_distance_by_sample/nasal/nasal_braydist_by_sample")

braydist_by_sample <- rbind(stool_braydist_by_sample, skin_braydist_by_sample,oral_braydist_by_sample,nasal_braydist_by_sample)
braydist_by_sample

braydist_by_sample$CL4 <- "Healthy"
braydist_by_sample$CL4[braydist_by_sample$sample_id1 %in% ls[str_detect(ls$CL4, "Infection"),]$SampleID] <- "Infection"
braydist_by_sample$CL4[braydist_by_sample$sample_id2 %in% ls[str_detect(ls$CL4, "Infection"),]$SampleID] <- "Infection"
braydist_by_sample$CL4[braydist_by_sample$sample_id1 %in% ls[str_detect(ls$CL4, "Imz"),]$SampleID] <- "Immunization"
braydist_by_sample$CL4[braydist_by_sample$sample_id2 %in% ls[str_detect(ls$CL4, "Imz"),]$SampleID] <- "Immunization"
braydist_by_sample$CL4[braydist_by_sample$sample_id1 %in% ls[str_detect(ls$CL4, "Ant"),]$SampleID] <- "Antibiotics"
braydist_by_sample$CL4[braydist_by_sample$sample_id2 %in% ls[str_detect(ls$CL4, "Ant"),]$SampleID] <- "Antibiotics"

braydist_by_sample$purt <- "all"
braydist_by_sample$purt[(braydist_by_sample$sample_id1 %in% ls[str_detect(ls$CL4, "Infection"),]$SampleID)&
                          ((braydist_by_sample$sample_id2 %in% ls[str_detect(ls$CL4, "Infection"),]$SampleID))] <- "Infection_within"
braydist_by_sample$purt[(braydist_by_sample$sample_id1 %in% ls[str_detect(ls$CL4, "Imz"),]$SampleID)&
                          ((braydist_by_sample$sample_id2 %in% ls[str_detect(ls$CL4, "Imz"),]$SampleID))] <- "Imz_within"
braydist_by_sample$purt[(braydist_by_sample$sample_id1 %in% ls[str_detect(ls$CL4, "Ant"),]$SampleID)&
                          ((braydist_by_sample$sample_id2 %in% ls[str_detect(ls$CL4, "Ant"),]$SampleID))] <- "Antibiotics_within"
table(braydist_by_sample$purt)

braydist_by_sample$IRIS <- "Unknown"
braydist_by_sample$IRIS[braydist_by_sample$subject_id1 %in% filter(sc, IRIS == "IR")$SubjectID] <- "IR"
braydist_by_sample$IRIS[braydist_by_sample$subject_id1 %in% filter(sc, IRIS == "IS")$SubjectID] <- "IS"
table(braydist_by_sample$IRIS)

summary(lmer(diffdays ~ dist + (1|subject_id1), stool_braydist_by_sample))

stool.lmer <- lmer(dist ~ diffdays + (1|subject_id1), stool_braydist_by_sample)
skin.lmer <- lmer(dist ~ diffdays + (1|subject_id1), skin_braydist_by_sample)
oral.lmer <- lmer(dist ~ diffdays + (1|subject_id1), oral_braydist_by_sample)
nasal.lmer <- lmer(dist ~ diffdays + (1|subject_id1), nasal_braydist_by_sample)

summary(stool.lmer)
summary(skin.lmer)
summary(oral.lmer)
summary(nasal.lmer)

combined.lmer <- lmer(dist ~ dataset:diffdays + (1|subject_id1), braydist_by_sample)
summary(combined.lmer)
coefplot2(combined.lmer, col = 1)

combined.lmer.IRIS <- lmer(dist ~ IRIS:dataset:diffdays + (1|subject_id1), filter(braydist_by_sample, IRIS != "Unknown" & purt == "all"))
summary(combined.lmer.IRIS)
coefplot2(combined.lmer.IRIS)
anova(combined.lmer.IRIS)

table(braydist_by_sample$CL4)

braydist_by_sample

braydist_by_sample$CL4 <- factor(braydist_by_sample$CL4, levels = c("Healthy","Infection", "Immunization","Antibiotics"))
braydist_by_sample$dataset <- factor(braydist_by_sample$dataset, levels=c("stool", "skin", "oral", "nasal"))

combined.lmer.healthy <- lmer(dist ~ IRIS:dataset:diffdays + (1|subject_id1), REML = F,filter(braydist_by_sample, CL4 == "Healthy" & diffdays > 120 &IRIS != "Unknown" & purt == "all"))
coefplot2(combined.lmer.healthy, xlim = c(-0.0002, 0.0004))
summary(combined.lmer.healthy)

combined.lmer.infection <- lmer(dist ~ IRIS:dataset:diffdays + (1|subject_id1), filter(braydist_by_sample, CL4 == "Infection" & diffdays > 120 &IRIS != "Unknown"& purt == "all"))
coefplot2(combined.lmer.infection,add=T, col=2, offset = -0.1)
summary(combined.lmer.infection)

combined.lmer.IMZ <- lmer(dist ~ IRIS:dataset:diffdays + (1|subject_id1), filter(braydist_by_sample, CL4 == "Immunization" & diffdays > 120 &IRIS != "Unknown"& purt == "all"))
coefplot2(combined.lmer.IMZ,add=T, col=3, offset = -0.2)

combined.lmer.anti <- lmer(dist ~ IRIS:dataset:diffdays + (1|subject_id1), filter(braydist_by_sample, CL4 == "Antibiotics" & diffdays > 120 &IRIS != "Unknown"& purt == "all"))
coefplot2(combined.lmer.anti,add=T, col=4, offset = -0.3)

combined.lmer.healthy <- lmer(dist ~ dataset:diffdays + (1|subject_id1), REML = F,filter(braydist_by_sample, CL4 == "Healthy" & IRIS == "IS"))
coefplot2(combined.lmer.healthy, xlim = c(-0.0002, 0.0004))
summary(combined.lmer.healthy)

combined.lmer.infection <- lmer(dist ~ dataset:diffdays + (1|subject_id1), filter(braydist_by_sample, CL4 == "Infection" & diffdays > 120 & IRIS == "IS"))
coefplot2(combined.lmer.infection,add=T, col=2, offset = -0.1)
summary(combined.lmer.infection)

combined.lmer.IMZ <- lmer(dist ~ dataset:diffdays + (1|subject_id1), filter(braydist_by_sample, CL4 == "Immunization" & diffdays > 120  & IRIS == "IS"))
coefplot2(combined.lmer.IMZ,add=T, col=3, offset = -0.2)

combined.lmer.anti <- lmer(dist ~ dataset:diffdays + (1|subject_id1), filter(braydist_by_sample, CL4 == "Antibiotics" & diffdays > 120  & IRIS == "IS"))
coefplot2(combined.lmer.anti,add=T, col=4, offset = -0.3)


combined.lmer.healthy <- lmer(dist ~ dataset:diffdays + (diffdays|subject_id1), filter(braydist_by_sample, CL4 == "Healthy" & diffdays > 120 & IRIS == "IR"))
coefplot2(combined.lmer.healthy, xlim = c(-0.0002, 0.0004))
summary(combined.lmer.healthy)

combined.lmer.infection <- lmer(dist ~ dataset:diffdays + (1|subject_id1), filter(braydist_by_sample, CL4 == "Infection" & diffdays > 120 & IRIS == "IR"))
coefplot2(combined.lmer.infection,add=T, col=2, offset = -0.1)
summary(combined.lmer.infection)

combined.lmer.IMZ <- lmer(dist ~ dataset:diffdays + (1|subject_id1), filter(braydist_by_sample, CL4 == "Immunization" & diffdays > 120 & IRIS == "IR"))
coefplot2(combined.lmer.IMZ,add=T, col=3, offset = -0.2)

combined.lmer.anti <- lmer(dist ~ dataset:diffdays + (1|subject_id1), filter(braydist_by_sample, CL4 == "Antibiotics" & diffdays > 120 & IRIS == "IR"))
coefplot2(combined.lmer.anti,add=T, col=4, offset = -0.3)


coefplot2(lmer(dist ~ IRIS:dataset+(1|subject_id1), filter(braydist_by_sample, CL4 == "Antibiotics" & purt =="all" & IRIS != "Unknown")))
coefplot2(lmer(dist ~ IRIS:dataset:diffdays+(1|subject_id1), filter(braydist_by_sample, CL4 == "Immunization" & purt =="all" & IRIS != "Unknown")))
coefplot2(lmer(dist ~ IRIS:dataset:diffdays+(1|subject_id1), filter(braydist_by_sample, CL4 == "Infection" & purt =="all" & IRIS != "Unknown")))
coefplot2(lmer(dist ~ IRIS:dataset:diffdays+(1|subject_id1), filter(braydist_by_sample, CL4 == "Healthy" & purt =="all" & IRIS != "Unknown")))

summary(lmer(dist ~ dataset:diffdays + (1|subject_id1) + IRIS, filter(braydist_by_sample, CL4 == "Antibiotics" & purt =="all")))
summary(lmer(dist ~ dataset:diffdays + (1|subject_id1) + IRIS, filter(braydist_by_sample, CL4 == "Infection" & purt =="all")))

max(braydist_by_sample$diffdays)/120



filter(braydist_by_sample, purt!="all") %>% ggplot(aes(x=CL4, y=dist)) + geom_boxplot() + facet_wrap(.~dataset) + stat_compare_means(label = "p.signif",ref.group = "Healthy")    

braydist_by_sample$x <- paste(braydist_by_sample$IRIS, braydist_by_sample$CL4, sep=".")
braydist_by_sample$x <- factor(braydist_by_sample$x, levels=c("IS.Healthy","IR.Healthy","IS.Infection", "IR.Infection",
                                                              "IS.Immunization","IR.Immunization","IS.Antibiotics","IR.Antibiotics",
                                                              "Unknown.Antibiotics", "Unknown.Healthy", "Unknown.Immunization", "Unknown.Infection"))
table(braydist_by_sample$x)

stability_purt <- filter(braydist_by_sample, purt=="all"&IRIS != "Unknown" & diffdays < 120 ) %>% ggplot(aes(x=x, y=dist)) +
  geom_jitter() + geom_boxplot()+ facet_wrap(.~dataset) + stat_compare_means(label = "p.signif",ref.group = "Healthy") + base_theme   
stability_purt

stability_purt_IS <- filter(braydist_by_sample, purt=="all"&IRIS == "IS") %>% ggplot(aes(x=CL4, y=dist)) + 
  geom_jitter(aes(color=dataset)) + geom_boxplot() + facet_wrap(.~dataset) + stat_compare_means(label = "p.signif",ref.group = "Healthy") + base_theme       
stability_purt_IS

#########
stability_purt_all_90 <- filter(braydist_by_sample,  diffdays < 90 & IRIS != "Unknown"& purt=="all") %>% ggplot(aes(x=x, y=dist)) + 
  geom_jitter(aes(color=dataset),size=0.05) + geom_boxplot(alpha=0.0) + facet_wrap(CL4~dataset, scale="free")  + base_theme  + scale_color_manual(values=body_site_color2) #+ stat_compare_means(label = "p.signif",label.y = 0.8)
stability_purt_all_90
#ggsave(filename = "./Analysis/Suppl.figure/stability_purt_all_90.pdf",stability_purt_all_90, width = 9, height = 10, dpi=300)

#test  
wilcox.test(filter(braydist_by_sample, diffdays < 90 & purt=="all" &  x == "IS.Healthy" & dataset == "stool")$dist, filter(braydist_by_sample, diffdays < 90 & purt=="all" &  x == "IR.Healthy" & dataset == "stool")$dist)
wilcox.test(filter(braydist_by_sample, diffdays < 90 & purt=="all" &  x == "IS.Healthy" & dataset == "skin")$dist, filter(braydist_by_sample, diffdays < 90 & purt=="all" &  x == "IR.Healthy" & dataset == "skin")$dist)
wilcox.test(filter(braydist_by_sample, diffdays < 90 & purt=="all" &  x == "IS.Healthy" & dataset == "oral")$dist, filter(braydist_by_sample, diffdays < 90 & purt=="all" &  x == "IR.Healthy" & dataset == "oral")$dist)
wilcox.test(filter(braydist_by_sample, diffdays < 90 & purt=="all" &  x == "IS.Healthy" & dataset == "nasal")$dist, filter(braydist_by_sample, diffdays < 90 & purt=="all" &  x == "IR.Healthy" & dataset == "nasal")$dist)

wilcox.test(filter(braydist_by_sample, diffdays < 90 & purt=="all" &  x == "IS.Infection" & dataset == "stool")$dist, filter(braydist_by_sample, diffdays < 90 & purt=="all" &  x == "IR.Infection" & dataset == "stool")$dist)
wilcox.test(filter(braydist_by_sample, diffdays < 90 & purt=="all" &  x == "IS.Infection" & dataset == "skin")$dist, filter(braydist_by_sample, diffdays < 90 & purt=="all" &  x == "IR.Infection" & dataset == "skin")$dist)
wilcox.test(filter(braydist_by_sample, diffdays < 90 & purt=="all" &  x == "IS.Infection" & dataset == "oral")$dist, filter(braydist_by_sample, diffdays < 90 & purt=="all" &  x == "IR.Infection" & dataset == "oral")$dist)
wilcox.test(filter(braydist_by_sample, diffdays < 90 & purt=="all" &  x == "IS.Infection" & dataset == "nasal")$dist, filter(braydist_by_sample, diffdays < 90 & purt=="all" &  x == "IR.Infection" & dataset == "nasal")$dist)

wilcox.test(filter(braydist_by_sample, diffdays < 90 & purt=="all" &  x == "IS.Immunization" & dataset == "stool")$dist, filter(braydist_by_sample, diffdays < 90 & purt=="all" &  x == "IR.Immunization" & dataset == "stool")$dist)
wilcox.test(filter(braydist_by_sample, diffdays < 90 & purt=="all" &  x == "IS.Immunization" & dataset == "skin")$dist, filter(braydist_by_sample, diffdays < 90 & purt=="all" &  x == "IR.Immunization" & dataset == "skin")$dist)
wilcox.test(filter(braydist_by_sample, diffdays < 90 & purt=="all" &  x == "IS.Immunization" & dataset == "oral")$dist, filter(braydist_by_sample, diffdays < 90 & purt=="all" &  x == "IR.Immunization" & dataset == "oral")$dist)
wilcox.test(filter(braydist_by_sample, diffdays < 90 & purt=="all" &  x == "IS.Immunization" & dataset == "nasal")$dist, filter(braydist_by_sample, diffdays < 90 & purt=="all" &  x == "IR.Immunization" & dataset == "nasal")$dist)

wilcox.test(filter(braydist_by_sample, diffdays < 90 & purt=="all" &  x == "IS.Antibiotics" & dataset == "stool")$dist, filter(braydist_by_sample, diffdays < 90 & purt=="all" &  x == "IR.Antibiotics" & dataset == "stool")$dist)
wilcox.test(filter(braydist_by_sample, diffdays < 90 & purt=="all" &  x == "IS.Antibiotics" & dataset == "skin")$dist, filter(braydist_by_sample, diffdays < 90 & purt=="all" &  x == "IR.Antibiotics" & dataset == "skin")$dist)
wilcox.test(filter(braydist_by_sample, diffdays < 90 & purt=="all" &  x == "IS.Antibiotics" & dataset == "oral")$dist, filter(braydist_by_sample, diffdays < 90 & purt=="all" &  x == "IR.Antibiotics" & dataset == "oral")$dist)
wilcox.test(filter(braydist_by_sample, diffdays < 90 & purt=="all" &  x == "IS.Antibiotics" & dataset == "nasal")$dist, filter(braydist_by_sample, diffdays < 90 & purt=="all" &  x == "IR.Antibiotics" & dataset == "nasal")$dist)


stability_purt_all <- filter(braydist_by_sample,  Time_Interval == "Y_01" & IRIS != "Unknown") %>% ggplot(aes(x=x, y=dist)) + 
  geom_jitter(aes(color=dataset)) + geom_boxplot(alpha=0.1) + facet_wrap(CL4~dataset, scale="free") + stat_compare_means(label = "p.signif") + base_theme  + scale_color_manual(values=body_site_color2)
stability_purt_all

IS.density <- filter(braydist_by_sample, purt=="all"&IRIS == "IS") %>% ggplot(aes(x=dist, color=CL4)) + geom_density() + facet_wrap(.~dataset) + ggtitle("IS") + base_theme
IR.density <- filter(braydist_by_sample, purt=="all"&IRIS == "IR") %>% ggplot(aes(x=dist, color=CL4)) + geom_density(linetype = "dashed") + facet_wrap(.~dataset) + ggtitle("IR")+ base_theme
density.plot <- IS.density + IR.density
#ggsave(filename = "./Analysis/Suppl.figure/density_distance.pdf",density.plot, width = 12, height = 5, dpi=300)

########################
stability_purt_all <- filter(braydist_by_sample, diffdays > 720  & IRIS != "Unknown") %>% ggplot(aes(x=x, y=dist)) + 
  geom_jitter(aes(color=dataset), size=1) + geom_boxplot(alpha=0.1) + facet_wrap(CL4~dataset, scale="free") + stat_compare_means(label = "p.signif") + base_theme  + scale_color_manual(values=body_site_color2)
stability_purt_all

table(braydist_by_sample$subject_id1)

stability_purt_IR <- filter(braydist_by_sample, purt=="all"&IRIS == "IR") %>% ggplot(aes(x=CL4, y=diffdays)) + 
  geom_jitter() + facet_wrap(.~dataset) + stat_compare_means(label = "p.signif",ref.group = "Healthy") + base_theme   
stability_purt_IR

stability_purt_IS <- filter(braydist_by_sample, purt=="all"&IRIS == "IS") %>% ggplot(aes(x=CL4, y=diffdays)) + 
  geom_jitter() + facet_wrap(.~dataset) + stat_compare_means(label = "p.signif",ref.group = "Healthy") + base_theme       
stability_purt_IS

diffday_purt_all <- filter(braydist_by_sample, diffdays > 720 & IRIS != "Unknown") %>% ggplot(aes(x=x, y=diffdays)) + 
  geom_boxplot() + facet_wrap(CL4~dataset, scale="free") + stat_compare_means(label = "p.signif") + base_theme   
diffday_purt_all

filter(braydist_by_sample, purt=="all" & IRIS != "Unknown") %>% ggplot(aes(x=diffdays, y=dist, group=IRIS)) + 
  geom_point() + geom_smooth(method="lm",aes(color=IRIS)) + facet_wrap(.~dataset)


summary(combined.lmer.IRIS)

plot(stool_braydist_by_sample$diffdays, stool_braydist_by_sample$dist)
plot(skin_braydist_by_sample$diffdays, skin_braydist_by_sample$dist)
plot(oral_braydist_by_sample$diffdays, oral_braydist_by_sample$dist)
plot(nasal_braydist_by_sample$diffdays, nasal_braydist_by_sample$dist)

braydist_by_sample

braydist_by_sample$Time_Interval <-paste("HY",formatC((findInterval(braydist_by_sample$diffdays, 1: 10* 180) + 1), width=2, flag="0"), sep="_")

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE),
      sem= (sd(x[[col]], na.rm=TRUE))/sqrt(length(x[[col]])))
       }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

df2 <-filter(braydist_by_sample, IRIS=="IS") %>%  
  data_summary(varname="dist", groupnames=c("CL4", "dataset","Time_Interval"))
df2$IRIS <- "IS"

df1 <-filter(braydist_by_sample,IRIS=="IR") %>%  
  data_summary(varname="dist", groupnames=c("CL4", "dataset","Time_Interval"))
df1$IRIS <- "IR"
#df2 <- filter(df2, IRIS != "Unknown")

p <- ggplot(df2, aes(x=Time_Interval, y=dist, group=CL4, color=dataset, shape=CL4)) + 
  geom_line(aes(linetype=CL4)) + facet_wrap(.~dataset)+base_theme+
  geom_point() + scale_color_manual(values = body_site_color2) + scale_linetype_manual(values=c("twodash", "solid","dotted","dashed"))  #+
 # geom_errorbar(aes(ymin=dist-sd, ymax=dist+sd), width=.2, position=position_dodge(0.05))
p

df3 <- rbind(df1,df2)
ggplot(df3, aes(x=Time_Interval, y=dist,group=IRIS, color=dataset)) + geom_jitter() + geom_boxplot(alpha=-0.1) + facet_wrap(dataset~CL4, ncol = 2, scales = "free")
ggplot(df3, aes(x=IRIS, y=dist, color=dataset)) + geom_jitter() + geom_boxplot(alpha=-0.1) + facet_wrap(dataset~., ncol = 2, scales = "free")

df3$time <- str_remove(df3$Time_Interval, "HY_")
df3$month <- as.numeric(df3$time) * 6

df3$CL4 <- factor(df3$CL4, levels = c("Healthy","Infection", "Immunization","Antibiotics"))
df3$dataset <- factor(df3$dataset, levels=c("stool", "skin", "oral", "nasal"))
p3 <- ggplot(df3, aes(x=month, y=dist, group=IRIS, color=dataset, shape=IRIS)) + 
  geom_line(aes(linetype=IRIS), size=0.5) + facet_grid(CL4~dataset, scale="free") + base_theme +
  geom_point(position=position_dodge(0.7), size= 0.7) + scale_color_manual(values = body_site_color2) + scale_linetype_manual(values=c("dashed", "solid")) +
 geom_errorbar(aes(ymin= dist-sem, ymax= dist+sem), width=.2, position=position_dodge(0.7), alpha=0.5)
p3
#ggsave(filename = "../../../human_microbiome_project/Figures/Figure3/Purtabation.figure.pdf",p3, width = 6, height = 4, dpi=300)

write.csv(file = "~/Desktop/Figure3Purtabation.table.csv",braydist_by_sample)

test <- filter(df3, CL4=="Infection" & dataset=="nasal")
test
test.plot <- ggline(test, x = "month", y = "dist", color = "IRIS",palette = c("#00AFBB", "#E7B800"))
test.plot
#####perform Two Way Anova
library(rstatix)
my_anova <- aov(dist ~ month * IRIS, data = test)
Anova(my_anova, type = "III")

#TukeyHSD(my_anova, which = "IRIS")
plot(my_anova, 2) 
shapiro.test(x = residuals(object = my_anova))


braydist_by_sample$Month_Interval <-paste("M",formatC((findInterval(braydist_by_sample$diffdays, 1: 60* 30) + 1), width=2, flag="0"), sep="_")
braydist_by_sample$Month_Interval

df2.1 <-filter(braydist_by_sample, IRIS=="IS") %>%  
  data_summary(varname="dist", groupnames=c("CL4", "dataset","Month_Interval"))
df2.1$IRIS <- "IS"

df1.1 <-filter(braydist_by_sample,IRIS=="IR") %>%  
  data_summary(varname="dist", groupnames=c("CL4", "dataset","Month_Interval"))
df1.1$IRIS <- "IR"

df3.1 <- rbind(df1.1,df2.1)

df3.1$time <- str_remove(df3.1$Month_Interval, "M_")
df3.1$month <- as.numeric(df3.1$time) * 1

p3.1 <- filter(df3.1,month < 4) %>% ggplot(aes(x=month, y=dist, group=IRIS, color=dataset, shape=IRIS)) + 
  geom_line(aes(linetype=IRIS), size=0.5) + facet_grid(CL4~dataset, scale="free") + base_theme +
  geom_point(position=position_dodge(0.1), size= 0.7) + scale_color_manual(values = body_site_color2) + scale_linetype_manual(values=c("dashed", "solid")) +
  geom_errorbar(aes(ymin= dist-sem, ymax= dist+sem), width=.2, position=position_dodge(0.1), alpha=0.5)
p3.1
ggsave(filename = "../../../human_microbiome_project/Figures/Figure3/Purtabation.figure_3month.pdf", p3.1, width = 6, height = 4, dpi=300)


status <- "Antibiotics"
bodysite <- "stool"

test.1 <- filter(df3.1, CL4==status & dataset==bodysite)
test.1

test.1$distlog <- log(test.1$dist + 1)

my_anova1 <- aov(dist ~ month * IRIS , data = (test.1))
pvalue <- Anova(my_anova1, type = "III")
pvalue
pdf(file = paste0("~/Library/CloudStorage/Box-Box/human_microbiome_project/Figures/Figure3/Purtabation/QQplot/", status,"_",bodysite,".pdf"), width = 5, height = 5)
plot(my_anova1, 2) 
text(1.3, -1.5, paste0("shapiro.test:",round(shapiro.test(x = residuals(object = my_anova1))$p, 4))) 
text(1.3, -1, paste0("Anova_IRIS:",round(pvalue$`Pr(>F)`[3], 4))) 
text(1.3, -0.5, paste0("beta_IRIS:",round(my_anova1$coefficients[3], 4))) 
dev.off()



pvalue$`Pr(>F)`[3]

shapiro.test(x = residuals(object = my_anova1))


combined.lmer.IRIS <- lmer(dist ~ IRIS:dataset:diffdays + (1|subject_id1), filter(braydist_by_sample, IRIS != "Unknown"))

summary(combined.lmer.IRIS)
#pdf(file = "~/Library/CloudStorage/Box-Box/human_microbiome_project/Figures/Figure3/Purtabation/IRIS_beta_plot.pdf", width = 8, height = 5)
coefplot2(combined.lmer.IRIS)
#dev.off()
anova(combined.lmer.IRIS)

combined.lmer.IRIS2 <- lmer(dist ~ dataset:diffdays + (1|subject_id1), filter(braydist_by_sample, IRIS != "Unknown"))
summary(combined.lmer.IRIS2)
coefplot2(combined.lmer.IRIS2)
anova(combined.lmer.IRIS2)

