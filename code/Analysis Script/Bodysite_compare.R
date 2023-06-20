#This is to compare stability between different individuals
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

setwd("~/Library/CloudStorage/Box-Box/XinZhouFiles/Projects/ZXP2_HMP2_Core/")

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

load("~/Library/CloudStorage/Box-Box/human_microbiome_project/data_analysis/combine_microbiome/within_subject_distance_by_sample/stool/stool_braydist_by_sample")
load("~/Library/CloudStorage/Box-Box/human_microbiome_project/data_analysis/combine_microbiome/within_subject_distance_by_sample/skin/skin_braydist_by_sample")
load("~/Library/CloudStorage/Box-Box/human_microbiome_project/data_analysis/combine_microbiome/within_subject_distance_by_sample/oral/oral_braydist_by_sample")
load("~/Library/CloudStorage/Box-Box/human_microbiome_project/data_analysis/combine_microbiome/within_subject_distance_by_sample/nasal/nasal_braydist_by_sample")

source("./Analysis/Analysis Script/sxt.tools.R")
load("./Analysis/Robject/Revision_MultiOmes_0509.RData")

dim(stool_braydist_by_sample)

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

braydist_by_sample$Time_Interval <-paste("Y",formatC((findInterval(braydist_by_sample$diffdays, 1: 10* 180) + 1), width=2, flag="0"), sep="_")

table(braydist_by_sample$subject_id1) %>% sort()

temp <- filter(braydist_by_sample, subject_id1 == "69-001")

ggplot(temp, aes(x=Time_Interval, y=dist, group=dataset)) + geom_boxplot()  + facet_wrap(.~dataset)


STSK <- merge(stool_braydist_by_sample, skin_braydist_by_sample, by=c("sample_id1","sample_id2"))
STOR <- merge(stool_braydist_by_sample, oral_braydist_by_sample, by=c("sample_id1","sample_id2"))
STNS <- merge(stool_braydist_by_sample, nasal_braydist_by_sample, by=c("sample_id1","sample_id2"))
SKOR <- merge(skin_braydist_by_sample, oral_braydist_by_sample, by=c("sample_id1","sample_id2"))
SKNS <- merge(skin_braydist_by_sample, nasal_braydist_by_sample, by=c("sample_id1","sample_id2"))
ORNS <- merge(oral_braydist_by_sample, nasal_braydist_by_sample, by=c("sample_id1","sample_id2"))

dim(stool_braydist_by_sample)
dim(skin_braydist_by_sample)
dim(oral_braydist_by_sample)
dim(nasal_braydist_by_sample)

table(STSK$subject_id1.x) %>% sort()
table(STOR$subject_id1.x) %>% sort()
table(STNS$subject_id1.x) %>% sort()
table(SKOR$subject_id1.x) %>% sort()
table(SKNS$subject_id1.x) %>% sort()
table(ORNS$subject_id1.x) %>% sort()

ggplot(temp2, aes(x=dist.x, y=dist.y, group=dataset.x)) + geom_point() + geom_smooth()
STSK.lm <- lmer(dist.x~dist.y + (1|subject_id1.x), data = STSK)
STOR.lm <- lmer(dist.x~dist.y + (1|subject_id1.x), data = STOR)
STNS.lm <- lmer(dist.x~dist.y + (1|subject_id1.x), data = STNS)
SKOR.lm <- lmer(dist.x~dist.y + (1|subject_id1.x), data = SKOR)
SKNS.lm <- lmer(dist.x~dist.y + (1|subject_id1.x), data = SKNS)
ORNS.lm <- lmer(dist.x~dist.y + (1|subject_id1.x), data = ORNS)

summary(STSK.lm)
summary(STOR.lm)
summary(STNS.lm)
summary(SKOR.lm)
summary(SKNS.lm)
summary(ORNS.lm)

coefplot2(STSK.lm,xlim=c(-0.1, 0.15))
coefplot2(STOR.lm,xlim=c(-0.1, 0.15), add=T, offset= 0.6)
coefplot2(STNS.lm,xlim=c(-0.1, 0.15), add=T, offset= 0.3)
coefplot2(SKOR.lm,xlim=c(-0.1, 0.15), add=T, offset= - 0.3)
coefplot2(SKNS.lm,xlim=c(-0.1, 0.15), add=T, offset= - 0.6)
coefplot2(ORNS.lm,xlim=c(-0.1, 0.15), add=T, offset= - 0.9)

summary(lm(dist.x~dist.y, data = filter(STSK, subject_id1.x == "69-001")))
summary(lm(dist.x~dist.y, data = filter(STOR, subject_id1.x == "69-001")))
summary(lm(dist.x~dist.y, data = filter(STNS, subject_id1.x == "69-001")))
summary(lm(dist.x~dist.y, data = filter(SKOR, subject_id1.x == "69-001")))
summary(lm(dist.x~dist.y, data = filter(SKNS, subject_id1.x == "69-001")))
summary(lm(dist.x~dist.y, data = filter(ORNS, subject_id1.x == "69-001")))

summary(lm(dist.x~dist.y, data = filter(STSK, subject_id1.x == "69-053")))
summary(lm(dist.x~dist.y, data = filter(STOR, subject_id1.x == "69-053")))
summary(lm(dist.x~dist.y, data = filter(STNS, subject_id1.x == "69-053")))
summary(lm(dist.x~dist.y, data = filter(SKOR, subject_id1.x == "69-053")))
summary(lm(dist.x~dist.y, data = filter(SKNS, subject_id1.x == "69-053")))
summary(lm(dist.x~dist.y, data = filter(ORNS, subject_id1.x == "69-053")))

ggplot(SKNS, aes(x=dist.x, y=dist.y)) + geom_point(size=0.05) + geom_smooth(method = "lm") + xlab("Skin") + ylab("Nasal")
ggplot(STSK, aes(x=dist.x, y=dist.y)) + geom_point(size=0.05) + geom_smooth(method = "lm") + xlab("Stool") + ylab("Skin")


#find the most unstable person 
length(as.character((table(stool_braydist_by_sample$subject_id1) %>% data.frame())$Var1))
length(as.character((table(skin_braydist_by_sample$subject_id1) %>% data.frame())$Var1))
length(as.character((table(oral_braydist_by_sample$subject_id1) %>% data.frame())$Var1))
length(as.character((table(nasal_braydist_by_sample$subject_id1) %>% data.frame())$Var1))

temp <- filter(stool_braydist_by_sample, subject_id1== "69-023")
cor.test(temp$dist, temp$diffdays)

#N=55
stool.estimate <- stool_braydist_by_sample %>% 
  select(subject_id1, dist, diffdays) %>%
  nest(-subject_id1) %>% 
  mutate(cor=map(data,~cor.test(.x$dist, .x$diffdays, method = "pearson"))) %>%
  mutate(tidied = map(cor, tidy)) %>% 
  unnest(tidied, .drop = T) %>% 
  select(-data, -cor) %>% mutate(p.FDR = p.adjust(p.value, method = "BH")) %>% mutate(Bodysite= "Stool")

#N=69
skin.estimate <- skin_braydist_by_sample %>% 
  select(subject_id1, dist, diffdays) %>%
  nest(-subject_id1) %>% 
  mutate(cor=map(data,~cor.test(.x$dist, .x$diffdays, method = "pearson"))) %>%
  mutate(tidied = map(cor, tidy)) %>% 
  unnest(tidied, .drop = T) %>% 
  select(-data, -cor) %>% mutate(p.FDR = p.adjust(p.value, method = "BH")) %>% mutate(Bodysite= "Skin")

#N=67
oral.estimate <- oral_braydist_by_sample %>% 
  select(subject_id1, dist, diffdays) %>%
  nest(-subject_id1) %>% 
  mutate(cor=map(data,~cor.test(.x$dist, .x$diffdays, method = "pearson"))) %>%
  mutate(tidied = map(cor, tidy)) %>% 
  unnest(tidied, .drop = T) %>% 
  select(-data, -cor) %>% mutate(p.FDR = p.adjust(p.value, method = "BH")) %>% mutate(Bodysite= "Oral")

#N=58
nasal.estimate <- nasal_braydist_by_sample %>% 
  select(subject_id1, dist, diffdays) %>%
  nest(-subject_id1) %>% 
  mutate(cor=map(data,~ cor.test(.x$dist, .x$diffdays, method = "pearson"))) %>%
  mutate(tidied = map(cor, tidy)) %>% 
  unnest(tidied, .drop = T) %>% 
  select(-data, -cor) %>% mutate(p.FDR = p.adjust(p.value, method = "BH")) %>% mutate(Bodysite= "Nasal")


stool.estimate
skin.estimate

p.stool <- ggplot(stool.estimate, aes(x=p.FDR, y=estimate)) + geom_point() + base_theme
p.stool

p.skin <- ggplot(skin.estimate, aes(x=p.FDR, y=estimate)) + geom_point() + base_theme
p.skin

p.oral <- ggplot(oral.estimate, aes(x=p.FDR, y=estimate)) + geom_point() + base_theme
p.oral

p.nasal <- ggplot(nasal.estimate, aes(x=p.FDR, y=estimate)) + geom_point() + base_theme
p.nasal

a <- mutate(stool.estimate, stool = estimate) %>% select(subject_id1,stool)
b <- mutate(skin.estimate, skin = estimate) %>% select(subject_id1,skin)
c <- mutate(oral.estimate, oral = estimate) %>% select(subject_id1,oral)
d <- mutate(nasal.estimate, nasal = estimate) %>% select(subject_id1,nasal)

estimation.matrix <- full_join(full_join(full_join(a, b, by="subject_id1"), c, by="subject_id1"),d, by="subject_id1")
e <- select(sc, SubjectID,IRIS)

estimation.matrix <- inner_join(estimation.matrix, e, by=c("subject_id1"="SubjectID"))


write.csv(file = "./Analysis/tables/beta.estimation.fourbodysite.csv",estimation.matrix)
