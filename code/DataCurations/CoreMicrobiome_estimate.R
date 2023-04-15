#Core Microbiome estimate
library(dplyr)
library(tidyverse)
library(summarytools)
library(Hmisc)
library(reshape2)
library(patchwork)
library(corrplot)

setwd("~/Box/XinZhouFiles/Projects/ZXP2_HMP2_Core/")

load("./Analysis/Robject/DetailedPhyloseq.RData")
load("./Analysis/Robject/Revision_MultiOmes_0509.RData")

#updated metadata, use the new version containing OGTT
sampleFreq <- read.csv(file = "./Analysis/ori_meta_table/SampleFrequency.csv", header = T)
#metadata <- read.csv("./OLD HMP Files/Analysis/Analysis/datamatrix_0425/threegroup0425_3cytokine_SD_Scale.csv", header = T, row.names = 1)
#metadata <- merge(metadata, sampleFreq, by.x="Row.names", by.y="SubjectID")
#colnames(metadata)[1] <- "SubjectID"
metadata <- read.csv("./Analysis/ori_meta_table/metadata.subject.csv", header = T,row.names = 1) 
metadata$SSPG <- as.numeric(metadata$SSPG)
metadata$FPG <- as.numeric(metadata$FPG)

ST_ASV_meta <- data.frame(sample_data(physeq_ST_freq)) %>% select(SampleID,SubjectID) %>%  data.frame(otu_table(physeq_ST_freq)) 
SK_ASV_meta <- data.frame(sample_data(physeq_SK_freq)) %>% select(SampleID,SubjectID) %>%  data.frame(otu_table(physeq_SK_freq)) 
OR_ASV_meta <- data.frame(sample_data(physeq_OR_freq)) %>% select(SampleID,SubjectID) %>%  data.frame(otu_table(physeq_OR_freq)) 
NS_ASV_meta <- data.frame(sample_data(physeq_NS_freq)) %>% select(SampleID,SubjectID) %>%  data.frame(otu_table(physeq_NS_freq)) 

ST_ASV_meta[1:10,1:10]
SK_ASV_meta[1:10,1:10]

table(duplicated(ST_ASV_meta$SampleID))
table(duplicated(SK_ASV_meta$SampleID))
table(duplicated(OR_ASV_meta$SampleID))
table(duplicated(NS_ASV_meta$SampleID))

####################################################################################
# Understand the low prevalance sample and its implications
####################################################################################
#stool
pre.abd_ST <- data.frame(colSums(ST_ASV_meta[-(1:2)] != 0), colSums(ST_ASV_meta[-(1:2)])/colSums(!!ST_ASV_meta[-(1:2)]))
colnames(pre.abd_ST) <- c("prevalance", "abundance")

stool.rankabundant <- data.frame()
for (i in 1:9) {
  print(i)
  stool.rankabundant[1,1:2] <- table(colSums(ST_ASV_meta[-(1:2)]) != 0)
  stool.rankabundant[(i + 1),1:2] <- table(colSums(ST_ASV_meta[-(1:2)] !=0) > i)
}
stool.rankabundant$V3 <- factor(c("all", 1:9), levels =c("all", 1:9) )
stool.rankabundant$V4 <- "Stool"
p.st.pr <- ggplot(stool.rankabundant, aes(x=V3, y=V2)) + geom_point()
p.st.pr

skin.rankabundant <- data.frame()
for (i in 1:9) {
  print(i)
  skin.rankabundant[1,1:2] <- table(colSums(SK_ASV_meta[-(1:2)]) != 0)
  skin.rankabundant[(i + 1),1:2] <- table(colSums(SK_ASV_meta[-(1:2)] !=0) > i)
}
skin.rankabundant$V3 <- factor(c("all", 1:9), levels =c("all", 1:9))
skin.rankabundant$V4 <- "Skin"
p.sk.pr <- ggplot(skin.rankabundant, aes(x=V3, y=V2)) + geom_point()
p.sk.pr

oral.rankabundant <- data.frame()
for (i in 1:9) {
  print(i)
  oral.rankabundant[1,1:2] <- table(colSums(OR_ASV_meta[-(1:2)]) != 0)
  oral.rankabundant[(i + 1),1:2] <- table(colSums(OR_ASV_meta[-(1:2)] !=0) > i)
}
oral.rankabundant$V3 <- factor(c("all", 1:9), levels =c("all", 1:9))
oral.rankabundant$V4 <- "Oral"
p.or.pr <- ggplot(oral.rankabundant, aes(x=V3, y=V2)) + geom_point()
p.or.pr

nasal.rankabundant <- data.frame()
for (i in 1:9) {
  print(i)
  nasal.rankabundant[1,1:2] <- table(colSums(NS_ASV_meta[-(1:2)]) != 0)
  nasal.rankabundant[(i + 1),1:2] <- table(colSums(NS_ASV_meta[-(1:2)] !=0) > i)
}
nasal.rankabundant$V3 <- factor(c("all", 1:9), levels =c("all", 1:9))
nasal.rankabundant$V4 <- "Nasal"
p.ns.pr <- ggplot(nasal.rankabundant, aes(x=V3, y=V2)) + geom_point()
p.ns.pr

rank.abundance <- rbind(stool.rankabundant,skin.rankabundant,oral.rankabundant,nasal.rankabundant)

p.rank <- ggplot(rank.abundance, aes(x=V3, y=V2, color=V4, group=V4)) + geom_point() + geom_line()
p.rank

p.pa.ST.1 <- ggplot(filter(pre.abd_ST, prevalance<10 & prevalance !=0 | prevalance == 30), aes(x=as.character(prevalance), y=log(abundance))) + geom_violin()
p.pa.ST.1 <- p.pa.ST.1 + ggtitle("Stool")
p.pa.ST.1

p.pa.ST.2 <- ggplot(filter(pre.abd_ST, prevalance !=0), aes(x=(prevalance), y=log(abundance))) + geom_point()
p.pa.ST.2 <- p.pa.ST.2 + ggtitle("Stool Microbiome Prevalence Abundance Plot")
p.pa.ST.2

#skin
pre.abd_SK <- data.frame(colSums(SK_ASV_meta[-(1:2)] != 0), colSums(SK_ASV_meta[-(1:2)])/colSums(!!SK_ASV_meta[-(1:2)]))
colnames(pre.abd_SK) <- c("prevalance", "abundance")

p.pa.SK.1 <- ggplot(filter(pre.abd_SK, prevalance<10 & prevalance !=0 | prevalance == 20 ), aes(x=as.character(prevalance), y=log(abundance))) + geom_violin()
p.pa.SK.1 <- p.pa.SK.1 + ggtitle("Skin")
p.pa.SK.1

p.pa.SK.2 <- ggplot(filter(pre.abd_SK, prevalance !=0), aes(x=(prevalance), y=log(abundance))) + geom_point()
p.pa.SK.2 <- p.pa.SK.2 + ggtitle("Skin Microbiome Prevalence Abundance Plot")
p.pa.SK.2

#oral
pre.abd_OR<- data.frame(colSums(OR_ASV_meta[-(1:2)] != 0), colSums(OR_ASV_meta[-(1:2)])/colSums(!!OR_ASV_meta[-(1:2)]))
colnames(pre.abd_OR) <- c("prevalance", "abundance")

p.pa.OR.1 <- ggplot(filter(pre.abd_OR, prevalance<10 & prevalance !=0 | prevalance == 20), aes(x=as.character(prevalance), y=log(abundance))) + geom_violin()
p.pa.OR.1 <- p.pa.OR.1 + ggtitle("Oral")
p.pa.OR.1

p.pa.OR.2 <- ggplot(filter(pre.abd_OR, prevalance !=0), aes(x=(prevalance), y=log(abundance))) + geom_point()
p.pa.OR.2 <- p.pa.OR.2 + ggtitle("Oral Microbiome Prevalence Abundance Plot")
p.pa.OR.2

#nasal
pre.abd_NS<- data.frame(colSums(NS_ASV_meta[-(1:2)] != 0), colSums(NS_ASV_meta[-(1:2)])/colSums(!!NS_ASV_meta[-(1:2)]))
colnames(pre.abd_NS) <- c("prevalance", "abundance")

p.pa.NS.1 <- ggplot(filter(pre.abd_NS, prevalance<10 & prevalance !=0 | prevalance == 30), aes(x=as.character(prevalance), y=log(abundance))) + geom_violin()
p.pa.NS.1 <- p.pa.NS.1 + ggtitle("Nasal")
p.pa.NS.1

p.pa.NS.2 <- ggplot(filter(pre.abd_NS, prevalance !=0), aes(x=(prevalance), y=log(abundance))) + geom_point()
p.pa.NS.2 <- p.pa.NS.2 + ggtitle("Nasal Microbiome Prevalence Abundance Plot")
p.pa.NS.2

p.pa.ST.1 + p.pa.SK.1 + p.pa.OR.1 + p.pa.NS.1
p.pa.ST.2 + p.pa.SK.2 + p.pa.OR.2 + p.pa.NS.2

##################################################################################################
#Calculate core genus
##################################################################################################

#########################################################
#stool
#########################################################
ST_genus_meta <- data.frame(sample_data(physeqGenus_ST)) %>% select(SampleID,SubjectID) %>%  data.frame(otu_table(physeqGenus_ST)) 
colnames(ST_genus_meta) <- tax_table(physeqGenus_ST)[,6][match(colnames(ST_genus_meta), rownames(tax_table(physeqGenus_ST)))]
colnames(ST_genus_meta)[1:2] <- c("SampleID", "SubjectID")

table(duplicated(ST_genus_meta$SampleID))
#subset a genus list with equal or more than 5 samples
ST_genus_meta_5 <- subset(ST_genus_meta, SubjectID %in% subset(sampleFreq, ST > 4)$SubjectID)
list5 <- as.character(subset(sampleFreq, ST > 4)$SubjectID)
ST.Pr <- data.frame()
for (i in 1: length(list5)) {
  print (i)
  print(list5[i])
  oper <- subset(ST_genus_meta_5, ST_genus_meta_5$SubjectID==list5[i])
  oper2 <- (colSums(oper!=0))/(max(colSums(oper==0)))
  ST.Pr[i,1:1233] <- oper2
}

colnames(ST.Pr) <- names(oper2)
row.names(ST.Pr) <- list5
ST.Pr <- select(ST.Pr, -SubjectID, -SampleID)

hist(unlist(ST.Pr[ST.Pr!=0]), breaks = 17, xlab="prevalence", main = "Prevalence of Stool Bacteria")
hist(unlist(ST.Pr[ST.Pr!=0 & ST.Pr!=1]), breaks = 17, xlab="prevalence", main = "Prevalence of Stool Bacteria")

a <- select(sampleFreq,SubjectID,ST)
b <- data.frame(rowSums(ST.Pr>=0.8),rowSums(ST.Pr != 0))
gdata <- merge(a,b, by.x="SubjectID", by.y="row.names")
colnames(gdata) <- c("SubjectID", "SampleTimes", "stability", "Observed_Genera")

pstool <- ggplot(gdata, aes(x=SampleTimes, y=stability)) + geom_point()
pstool

gdatameta <- merge(gdata,metadata,by="SubjectID")
gdatameta$ratio <- gdatameta$stability/gdatameta$Observed_Genera

p.ST.sta <- ggplot(gdatameta, aes(x=OGTT, y=stability)) + geom_point()  # + geom_point(aes(size=SampleTimes, color=as.factor(group)))
p.ST.sta <- p.ST.sta + geom_smooth(method="lm",se = F) + cowplot::theme_cowplot(12) + ylim(0, 60)
p.ST.sta

ST.lm <- lm(formula = FPG ~  stability + Observed_Genera + Adj.age + Gender , data=gdatameta)
summary(ST.lm)

p.ST.sta2 <- ggplot(subset(gdatameta, SubjectID != "69-000"), aes(x=stability, y=Observed_Genera)) + geom_point()  # + geom_point(aes(size=SampleTimes, color=as.factor(group)))
p.ST.sta2 <- p.ST.sta2 + geom_smooth(method="lm",se = 0.95) + cowplot::theme_cowplot(12) #+ ylim(0, 60)
p.ST.sta2

ST.lm2 <- lm(formula = stability ~ Observed_Genera, data=gdatameta)
summary(ST.lm2)

#########################################################
#skin
#########################################################
SK_genus_meta <- data.frame(sample_data(physeqGenus_SK)) %>% select(SampleID,SubjectID) %>%  data.frame(otu_table(physeqGenus_SK)) 
colnames(SK_genus_meta) <- tax_table(physeqGenus_SK)[,6][match(colnames(SK_genus_meta), rownames(tax_table(physeqGenus_SK)))]
colnames(SK_genus_meta)[1:2] <- c("SampleID", "SubjectID")
SK_genus_meta[1:5, 1:10]

#found a duplicate sample
table(duplicated(SK_genus_meta$SampleID))

#subset a genus list with equal or more than 5 samples
SK_genus_meta_5 <- subset(SK_genus_meta, SubjectID %in% subset(sampleFreq, SK > 4)$SubjectID)
liSK5 <- as.character(subset(sampleFreq, SK > 4)$SubjectID)
SK.Pr <- data.frame()
for (i in 1: length(liSK5)) {
  print (i)
  print(liSK5[i])
  oper <- subset(SK_genus_meta_5, SK_genus_meta_5$SubjectID==liSK5[i])
  oper2 <- (colSums(oper!=0))/(max(colSums(oper==0)))
  SK.Pr[i,1:1130] <- oper2
}

colnames(SK.Pr) <- names(oper2)
row.names(SK.Pr) <- liSK5
SK.Pr <- select(SK.Pr, -SubjectID, -SampleID)

hist(unlist(SK.Pr[SK.Pr!=0]), breaks = 17, xlab="prevalence", main = "Prevalence of Skin Bacteria")
hist(unlist(SK.Pr[SK.Pr!=0 & SK.Pr!=1]), breaks = 17, xlab="prevalence", main = "Prevalence of Skinl Bacteria")

a <- select(sampleFreq,SubjectID,SK)
b <- data.frame(rowSums(SK.Pr>=0.8),rowSums(SK.Pr != 0))
gdata_SK <- merge(a,b, by.x="SubjectID", by.y="row.names")
colnames(gdata_SK) <- c("SubjectID", "SampleTimes", "stability", "Observed_Genera")

pSKool <- ggplot(gdata_SK, aes(x=SampleTimes, y=stability)) + geom_point()
pSKool

gdatameta_SK <- merge(gdata_SK,metadata,by="SubjectID")
gdatameta_SK$ratio <- gdatameta_SK$stability/gdatameta_SK$Observed_Genera

p.SK.sta <- ggplot(gdatameta_SK, aes(x=OGTT, y=stability)) + geom_point()  # + geom_point(aes(size=SampleTimes, color=as.factor(group)))
p.SK.sta <- p.SK.sta + geom_smooth(method="lm",se = 0.95) + cowplot::theme_cowplot(12) + ylim(0, 60)
p.SK.sta

SK.lm <- lm(formula = OGTT ~ stability + Observed_Genera + Adj.age + Gender, data=gdatameta_SK)
summary(SK.lm)

p.SK.sta2 <- ggplot(subset(gdatameta_SK, SubjectID != "69-000"), aes(x=stability, y=Observed_Genera)) + geom_point()  # + geom_point(aes(size=SampleTimes, color=as.factor(group)))
p.SK.sta2 <- p.SK.sta2 + geom_smooth(method="lm",se = 0.95) + cowplot::theme_cowplot(12) #+ ylim(0, 0.25)
p.SK.sta2

SK.lm2 <- lm(formula = stability ~ Observed_Genera, data=subset(gdatameta_SK, SubjectID != "69-000"))
summary(SK.lm2)

#########################################################
#oral
#########################################################
OR_genus_meta <- data.frame(sample_data(physeqGenus_OR)) %>% select(SampleID,SubjectID) %>%  data.frame(otu_table(physeqGenus_OR)) 
colnames(OR_genus_meta) <- tax_table(physeqGenus_OR)[,6][match(colnames(OR_genus_meta), rownames(tax_table(physeqGenus_OR)))]
colnames(OR_genus_meta)[1:2] <- c("SampleID", "SubjectID")
OR_genus_meta[1:5, 1:10]

#found duplicate sample
table(duplicated(OR_genus_meta$SampleID))

#subset a genus list with equal or more than 5 samples
OR_genus_meta_5 <- subset(OR_genus_meta, SubjectID %in% subset(sampleFreq, OR > 4)$SubjectID)
liOR5 <- as.character(subset(sampleFreq, OR > 4)$SubjectID)
OR.Pr <- data.frame()
for (i in 1: length(liOR5)) {
  print (i)
  print(liOR5[i])
  oper <- subset(OR_genus_meta_5, OR_genus_meta_5$SubjectID==liOR5[i])
  oper2 <- (colSums(oper!=0))/(max(colSums(oper==0)))
  OR.Pr[i,1:1130] <- oper2
}

colnames(OR.Pr) <- names(oper2)
row.names(OR.Pr) <- liOR5
OR.Pr <- select(OR.Pr, -SubjectID, -SampleID)

hist(unlist(OR.Pr[OR.Pr!=0]), breaks = 17, xlab="prevalence", main = "Prevalence of Oral Bacteria")
hist(unlist(OR.Pr[OR.Pr!=0 & OR.Pr!=1]), breaks = 17, xlab="prevalence", main = "Prevalence of Oral Bacteria")

a <- select(sampleFreq,SubjectID,OR)
b <- data.frame(rowSums(OR.Pr>=0.8),rowSums(OR.Pr != 0))
gdata_OR <- merge(a,b, by.x="SubjectID", by.y="row.names")
colnames(gdata_OR) <- c("SubjectID", "SampleTimes", "stability", "Observed_Genera")

pORool <- ggplot(gdata_OR, aes(x=SampleTimes, y=stability)) + geom_point()
pORool

gdatameta_OR <- merge(gdata_OR,metadata,by="SubjectID")
gdatameta_OR$ratio <- gdatameta_OR$stability/gdatameta_OR$Observed_Genera

p.OR.sta <- ggplot(gdatameta_OR, aes(x=OGTT, y=stability)) + geom_point()  # + geom_point(aes(size=SampleTimes, color=as.factor(group)))
p.OR.sta <- p.OR.sta + geom_smooth(method="lm",se = 0.95) + cowplot::theme_cowplot(12) + ylim(0, 60)
p.OR.sta

OR.lm <- lm(formula = OGTT ~ stability + Observed_Genera + Adj.age + Gender, data=gdatameta_OR)
summary(OR.lm)

p.OR.sta2 <- ggplot(subset(gdatameta_OR, SubjectID != "69-000"), aes(x=stability, y=Observed_Genera)) + geom_point()  # + geom_point(aes(size=SampleTimes, color=as.factor(group)))
p.OR.sta2 <- p.OR.sta2 + geom_smooth(method="lm",se = 0.95) + cowplot::theme_cowplot(12)# + ylim(0, 60)
p.OR.sta2

OR.lm2 <- lm(formula = stability ~  Observed_Genera, data=subset(gdatameta_OR, SubjectID != "69-000"))
summary(OR.lm2)

#########################################################
#nasal
#########################################################
NS_genus_meta <- data.frame(sample_data(physeqGenus_NS)) %>% select(SampleID,SubjectID) %>%  data.frame(otu_table(physeqGenus_NS)) 
colnames(NS_genus_meta) <- tax_table(physeqGenus_NS)[,6][match(colnames(NS_genus_meta), rownames(tax_table(physeqGenus_NS)))]
colnames(NS_genus_meta)[1:2] <- c("SampleID", "SubjectID")
NS_genus_meta[1:5, 1:10]

#found duplicate sample
table(duplicated(NS_genus_meta$SampleID))

#subset a genus list with equal NS mNSe than 5 samples
NS_genus_meta_5 <- subset(NS_genus_meta, SubjectID %in% subset(sampleFreq, NS > 4)$SubjectID)
liNS5 <- as.character(subset(sampleFreq, NS > 4)$SubjectID)
NS.Pr <- data.frame()
for (i in 1: length(liNS5)) {
  print (i)
  print(liNS5[i])
  oper <- subset(NS_genus_meta_5, NS_genus_meta_5$SubjectID==liNS5[i])
  oper2 <- (colSums(oper!=0))/(max(colSums(oper==0)))
  NS.Pr[i,1:1233] <- oper2
}

colnames(NS.Pr) <- names(oper2)
row.names(NS.Pr) <- liNS5
NS.Pr <- select(NS.Pr, -SubjectID, -SampleID)

hist(unlist(NS.Pr[NS.Pr!=0]), breaks = 17, xlab="prevalence", main = "Prevalence of Nasal Bacteria")
hist(unlist(NS.Pr[NS.Pr!=0 & NS.Pr!=1]), breaks = 17, xlab="prevalence", main = "Prevalence of Nasal Bacteria")

a <- select(sampleFreq,SubjectID,NS)
b <- data.frame(rowSums(NS.Pr>=0.8),rowSums(NS.Pr != 0))
gdata_NS <- merge(a,b, by.x="SubjectID", by.y="row.names")
colnames(gdata_NS) <- c("SubjectID", "SampleTimes", "stability", "Observed_Genera")

pNSool <- ggplot(gdata_NS, aes(x=SampleTimes, y=stability)) + geom_point()
pNSool

gdatameta_NS <- merge(gdata_NS,metadata,by="SubjectID")
gdatameta_NS$ratio <- gdatameta_NS$stability/gdatameta_NS$Observed_Genera

p.NS.sta <- ggplot(gdatameta_NS, aes(x=OGTT, y=stability)) + geom_point()  # + geom_point(aes(size=SampleTimes, color=as.factor(group)))
p.NS.sta <- p.NS.sta + geom_smooth(method="lm",se = 0.95) + cowplot::theme_cowplot(12) + ylim(0, 60)
p.NS.sta

NS.lm <- lm(formula = OGTT ~ stability + Observed_Genera + Adj.age + Gender, data=gdatameta_NS)
summary(NS.lm)

p.NS.sta2 <- ggplot(subset(gdatameta_NS, SubjectID != "69-000"), aes(x=stability, y=Observed_Genera)) + geom_point()  # + geom_point(aes(size=SampleTimes, color=as.factor(group)))
p.NS.sta2 <- p.NS.sta2 + geom_smooth(method="lm",se = 0.95) + cowplot::theme_cowplot(12) 
#p.NS.sta2 <- p.NS.sta2 + facet_wrap(.~IRIS)
p.NS.sta2

NS.lm2 <- lm(formula = stability ~ Observed_Genera, data=gdatameta_NS)
summary(NS.lm2)

p.ST.sta + p.SK.sta + p.OR.sta + p.NS.sta

p.ST.sta2 + p.SK.sta2 + p.OR.sta2 + p.NS.sta2


par(mfrow=c(4,1))
hist(unlist(ST.Pr[ST.Pr!=0 & ST.Pr!=1]), breaks = 17, xlab="prevalence", main = "Prevalence of Stool Bacteria")
hist(unlist(SK.Pr[SK.Pr!=0 & SK.Pr!=1]), breaks = 17, xlab="prevalence", main = "Prevalence of Skin Bacteria")
hist(unlist(OR.Pr[OR.Pr!=0 & OR.Pr!=1]), breaks = 17, xlab="prevalence", main = "Prevalence of Oral Bacteria")
hist(unlist(NS.Pr[NS.Pr!=0 & NS.Pr!=1]), breaks = 17, xlab="prevalence", main = "Prevalence of Nasal Bacteria")

summary(ST.lm)
summary(SK.lm)
summary(OR.lm)
summary(NS.lm)

a1 <- select(gdatameta,SubjectID,stability,Observed_Genera,ratio)
colnames(a1) <- c("SubjectID", "Stool_NCore", "Stool_NGenera", "Stool_Ratio")
b1 <- select(gdatameta_SK,SubjectID,stability,Observed_Genera, ratio)
colnames(b1) <- c("SubjectID", "Skin_NCore", "Skin_NGenera", "Skin_Ratio")
c1 <- select(gdatameta_OR,SubjectID,stability,Observed_Genera, ratio)
colnames(c1) <- c("SubjectID", "Oral_NCore", "Oral_NGenera", "Oral_Ratio")
d1 <- select(gdatameta_NS,SubjectID,stability,Observed_Genera, ratio)
colnames(d1) <- c("SubjectID", "Nasal_NCore", "Nasal_NGenera", "Nasal_Ratio")

NCore <- merge(merge(a1,b1, by="SubjectID"), merge(c1,d1, by="SubjectID"), by="SubjectID")
NCore

pcore <- ggplot(NCore, aes(x=Stool_Ratio, y=Skin_Ratio)) + geom_point() + geom_smooth(method="lm",color="red")
pcore

pcore2 <- ggplot(NCore, aes(x=Stool_Ratio, y=Oral_Ratio)) + geom_point() + geom_smooth(method="lm")
pcore2

pcore3 <- ggplot(NCore, aes(x=Stool_Ratio, y=Nasal_Ratio)) + geom_point() + geom_smooth(method="lm")
pcore3

pcore4 <- ggplot(NCore, aes(x=Skin_Ratio, y=Oral_Ratio)) + geom_point() + geom_smooth(method="lm")
pcore4

pcore5 <- ggplot(NCore, aes(x=Skin_Ratio, y=Nasal_Ratio)) + geom_point() + geom_smooth(method="lm")
pcore5

pcore6 <- ggplot(NCore, aes(x=Oral_Ratio, y=Nasal_Ratio)) + geom_point() + geom_smooth(method="lm")
pcore6

pcore + pcore2 + pcore3 + pcore4 + pcore5 + pcore6

pobser<- ggplot(NCore, aes(x=Stool_NGenera, y=Skin_NGenera)) + geom_point() + geom_smooth(method="lm",color="red")
pobser

pobser2 <- ggplot(NCore, aes(x=Stool_NGenera, y=Oral_NGenera)) + geom_point() + geom_smooth(method="lm")
pobser2

pobser3 <- ggplot(NCore, aes(x=Stool_NGenera, y=Nasal_NGenera)) + geom_point() + geom_smooth(method="lm")
pobser3

pobser4 <- ggplot(NCore, aes(x=Skin_NGenera, y=Oral_NGenera)) + geom_point() + geom_smooth(method="lm")
pobser4

pobser5 <- ggplot(NCore, aes(x=Skin_NGenera, y=Nasal_NGenera)) + geom_point() + geom_smooth(method="lm")
pobser5

pobser6 <- ggplot(NCore, aes(x=Oral_NGenera, y=Nasal_NGenera)) + geom_point() + geom_smooth(method="lm")
pobser6

pobser + pobser2 + pobser3 + pobser4 + pobser5 + pobser6

row.names(NCore) <- NCore$SubjectID
NCore.meta <- merge(NCore, select(metadata,SubjectID,Th17.Group,SSPG,FPG,FPG_Mean,OGTT,Adj.age,BMI), by="SubjectID")
row.names(NCore.meta) <- NCore.meta$SubjectID

p.BMI.SSPG <- ggplot(NCore.meta, aes(x=SSPG, y=BMI)) + geom_point() + theme_cowplot() + geom_smooth(method = "lm")
p.BMI.SSPG

cor.Ncore <- rcorr(as.matrix(NCore.meta[,-1]), type="spearman")

corrplot(cor.Ncore$r, type="upper", order="original", p.mat = cor.Ncore$P, sig.level = 0.01, insig = "blank")

############################################################################################################################################################################################################

############################################################################################################################################################################################################

table(colSums(ST.Pr) == 0)
table(colSums(SK.Pr) == 0)
table(colSums(OR.Pr) == 0)
table(colSums(NS.Pr) == 0)

ST.Pr.Char <- ST.Pr[,colSums(ST.Pr) != 0]
ST.Pr.Char[ST.Pr.Char >= 0.8] <- "core"
ST.Pr.Char[ST.Pr.Char < 0.8 & ST.Pr.Char> 0.2] <- "interm"
ST.Pr.Char[ ST.Pr.Char<= 0.2] <- "oppor"

dat.ST <- data.frame()
dat.ST[1:3,1] <- c("core", "interm", "oppor")
colnames(dat.ST)[1] <- "Var1"
dim(ST.Pr.Char)
for (i in 1:360){
  print(i)
  dat1 <- data.frame(table(sapply(ST.Pr.Char[,i], function(x) x)))
  colnames(dat1)[2] <- colnames(ST.Pr.Char)[i]
  dat.ST <- merge(dat.ST,dat1, by="Var1", all=T)
  }

dat.ST[is.na(dat.ST)] <- 0

dat.ST.melt <- melt(dat.ST,id.vars = "Var1")
listst <- as.character(dat.ST.melt$variable[dat.ST.melt$value==55&dat.ST.melt$Var1=="oppor"])
dat.ST.melt$cater <- "cat1"
dat.ST.melt$cater[which(dat.ST.melt$variable %in% listst)] <- "cat2"

order.ST <- data.frame(t(dat.ST))[2:361,]
order.ST$X1 <- as.numeric(order.ST$X1)
order.ST$X2 <- as.numeric(order.ST$X2)
order.ST$X3 <- as.numeric(order.ST$X3)
order.ST <- order.ST[with(order.ST, order(X1, X2, X3)),]
order.ST$genus <- row.names(order.ST)
order.ST$order <- 1:360

dat.ST.melt$order <- order.ST$order[match(dat.ST.melt$variable,order.ST$genus)]

pcore.percent1 <- ggplot(dat.ST.melt, mapping=aes(y=reorder(variable,order), x=value, fill=Var1)) + geom_bar(position = "stack",stat="identity")
pcore.percent1 <- pcore.percent1 + facet_wrap(.~cater, scales = "free_y")
pcore.percent1

#ggsave(filename = "./Analysis/Suppl.figure/Core.microbiome.Pre.ST.pdf",pcore.percent1, width=12, height=24, dpi=300)

#seperate IRIS
ST.Pr.Char.IS <- ST.Pr.Char[row.names(ST.Pr.Char) %in% subset(sc,IRIS=="IS")$SubjectID,]
ST.Pr.Char.IR <- ST.Pr.Char[row.names(ST.Pr.Char) %in% subset(sc,IRIS=="IR")$SubjectID,]

dat.ST.IS <- data.frame()
dat.ST.IS[1:3,1] <- c("core", "interm", "oppor")
colnames(dat.ST.IS)[1] <- "Var1"
dim(ST.Pr.Char.IS)
for (i in 1:360){
  print(i)
  dat1 <- data.frame(table(sapply(ST.Pr.Char.IS[,i], function(x) x)))
  colnames(dat1)[2] <- colnames(ST.Pr.Char.IS)[i]
  dat.ST.IS <- merge(dat.ST.IS,dat1, by="Var1", all=T)
}
dat.ST.IS[is.na(dat.ST.IS)] <- 0
dat.ST.IS$IRIS <- "IS"

dat.ST.IR <- data.frame()
dat.ST.IR[1:3,1] <- c("core", "interm", "oppor")
colnames(dat.ST.IR)[1] <- "Var1"
dim(ST.Pr.Char.IR)
for (i in 1:360){
  print(i)
  dat1 <- data.frame(table(sapply(ST.Pr.Char.IR[,i], function(x) x)))
  colnames(dat1)[2] <- colnames(ST.Pr.Char.IR)[i]
  dat.ST.IR <- merge(dat.ST.IR,dat1, by="Var1", all=T)
}
dat.ST.IR[is.na(dat.ST.IR)] <- 0
dat.ST.IR$IRIS <- "IR"

#combine two dataset
dat.ST.melt.IRIS <- melt(rbind(dat.ST.IS,dat.ST.IR),id.vars = c("Var1", "IRIS"))
dat.ST.melt.IRIS$cater <- "cat1"
dat.ST.melt.IRIS$cater[which(dat.ST.melt.IRIS$variable %in% listst)] <- "cat2"

#adjust the plot order
orderIS <- data.frame(t(dat.ST.IS))[2:361,]
orderIS$X1 <- as.numeric(orderIS$X1)
orderIS$X2 <- as.numeric(orderIS$X2)
orderIS$X3 <- as.numeric(orderIS$X3)
orderIS <- orderIS[with(orderIS, order(X1, X2, X3)),]
orderIS$genus <- row.names(orderIS)
orderIS$order <- 1:360

#dat.ST.melt.IRIS$order <- dat.ST.melt.IRIS$variable

dat.ST.melt.IRIS$order <- orderIS$order[match(dat.ST.melt.IRIS$variable,orderIS$genus)]
dat.ST.melt.IRIS$IRIS <- factor(dat.ST.melt.IRIS$IRIS, levels=c("IS", "IR"))

pcore.percent3 <- filter(dat.ST.melt.IRIS, cater=="cat1") %>% ggplot(mapping=aes(y=reorder(variable,order), x=value, fill=Var1)) + geom_bar(position = "stack",stat="identity")
pcore.percent3 <- pcore.percent3 + facet_wrap(.~IRIS)
pcore.percent3

#ggsave(filename = "./Analysis/Suppl.figure/Core.microbiome.Pre.ST.IRIS.pdf",pcore.percent3, width=12, height=24, dpi=300)

subset(dat.ST.melt.IRIS, variable == "Coprococcus")$value[1:3]

for (i in 1:360){
  print(i)
  taxa.i <- as.character(unique(dat.ST.melt.IRIS$variable))[i]
  print(taxa.i)
  
  taxa.chi <- chisq.test(cbind(subset(dat.ST.melt.IRIS, variable == taxa.i)$value[1:3],
                   subset(dat.ST.melt.IRIS, variable == taxa.i)$value[4:6]))
}
subset(dat.ST.melt.IRIS, variable == "Butyrivibrio")$value

taxa.chi <- chisq.test(cbind(subset(dat.ST.melt.IRIS, variable == "Barnesiella")$value[1:3],
                             subset(dat.ST.melt.IRIS, variable == "Barnesiella")$value[4:6]))
taxa.chi

taxa.fisher <- fisher.test(cbind(subset(dat.ST.melt.IRIS, variable == "Barnesiella")$value[1:3],
                  subset(dat.ST.melt.IRIS, variable == "Barnesiella")$value[4:6]))

#Coprococcus
fisher.test(cbind(c(20,1),c(13,7)))
chisq.test(cbind(c(20,1),c(13,7)))

#Parasutterella
fisher.test(cbind(c(14,7),c(6,14)))
chisq.test(cbind(c(14,7),c(6,14)))

#Intestinimonas
fisher.test(cbind(c(10,11),c(8,12)))
chisq.test(cbind(c(10,11),c(8,12)))

#Butyricimonas
fisher.test(cbind(c(10,11),c(5,15)))
chisq.test(cbind(c(10,11),c(5,15)))

#Butyricicoccus
fisher.test(cbind(c(9,12),c(2,18)))
chisq.test(cbind(c(9,12),c(2,18)))

#Butyrivibrio
fisher.test(cbind(c(4,17),c(0,20)))
chisq.test(cbind(c(4,17),c(0,20)))


taxa.chi$p.value

chisq.test(cbind(subset(dat.ST.melt.IRIS, variable == "Odoribacter")$value[1:3],
           subset(dat.ST.melt.IRIS, variable == "Odoribacter")$value[4:6]))


#subset(data.frame(tax_table(physeq_ST)),Genus=="Unclassified_Alphaproteobacteria")
#subset(data.frame(tax_table(physeq_ST)),Genus=="Unclassified_Lachnospiraceae")
#subset(data.frame(tax_table(physeq_ST)),Genus=="Unclassified_Ruminococcaceae")
#subset(data.frame(tax_table(physeq_ST)),Genus=="Unclassified_Clostridiales")

#save(ST.Pr, SK.Pr, OR.Pr, NS.Pr, file = "./Analysis/Robject/Prevalance.RData")
#write.csv(file = "./Analysis/ori_meta_table/CoreMicrobiome.csv",NCore)

#dat.ST.melt[dat.ST.melt$value==55&dat.ST.melt$Var1=="oppor",]
########################################################################################################
#Proveide data for ZJU
########################################################################################################

HMP_meta$Random.Subject.ID <- sub("\\-.*", "",HMP_meta$ori.ID)
colnames(HMP_meta)
table(HMP_meta$Random.Subject.ID)

colnames(ns.df)

HMP_meta_ZJU <- select(HMP_meta, RandomID,Random.Subject.ID,Date,batch,SampleType)
HMP_meta_ZJU 
OTU_ZJU <- data.frame(otu_table(physeq_HMP))
TAX_ZJU <- data.frame(tax_table(physeq_HMP))

st.df.zju <- st.df
st.df.zju$SubjectID <- sc$rand_subject_id[match(st.df.zju$SubjectID, sc$SubjectID)]
st.df.zju <- select(st.df.zju, -HostSampleID)

ns.df.zju <- ns.df
ns.df.zju$SubjectID <- sc$rand_subject_id[match(ns.df.zju$SubjectID, sc$SubjectID)]
ns.df.zju <- select(ns.df.zju, -HostSampleID)

#save(st.df.zju,ns.df.zju,file = "~/Desktop/HMP_data_ZJU.RData")
#rm(HMP_meta_ZJU, OTU_ZJU, TAX_ZJU)
