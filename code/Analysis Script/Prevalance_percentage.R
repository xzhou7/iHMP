
#Prevalance Analysis (figure 1/3)
#this analysis is intent to generate core microbiome percentage
#Author: Xin Zhou, Ph.D.
#Last updated: Nov.14.2021
library(x)
library(phyloseq)
library(dplyr)
library(stringr)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(corrplot)
library(Hmisc)
library(reshape2)
library(ggcorrplot)

base_theme = theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 13),
        panel.grid.minor = element_blank())
body_site_color = c(
  "Stool" = ggsci::pal_jama()(n=7)[2],
  "Skin" = ggsci::pal_jama()(n=7)[3],
  "Oral" = ggsci::pal_jama()(n=7)[4],
  "Nasal" = ggsci::pal_jama()(n=7)[5])

setwd("~/Library/CloudStorage/Box-Box/XinZhouFiles/Projects/ZXP2_HMP2_Core/")
load("./Analysis/Robject/Prevalance.RData")
load("./Analysis/Robject/Revision_MultiOmes_0509.RData")
source("../../../human_microbiome_project/code/tools.R")

corenumber <- read.csv("./Analysis/ori_meta_table/Core_microbiome/CoreMicrobiome_observed_Subject.csv")
Genus_ST <- read.csv("./Analysis/Genus Table/ST/Genus_ST.csv")
Genus_SK <- read.csv("./Analysis/Genus Table/SK/Genus_SK.csv")
Genus_OR <- read.csv("./Analysis/Genus Table/OR/Genus_OR.csv")
Genus_NS <- read.csv("./Analysis/Genus Table/NS/Genus_NS.csv")
sc <- read.csv("./Analysis/ori_meta_table/metadata.subject.csv")

ST.Pr.Char <- ST.Pr[,colSums(ST.Pr) != 0]
ST.Pr.Char[ST.Pr.Char >= 0.8] <- "core"
ST.Pr.Char[ST.Pr.Char < 0.8 & ST.Pr.Char> 0.2] <- "interm"
ST.Pr.Char[ST.Pr.Char<= 0.2] <- "oppor"

SK.Pr.Char <- SK.Pr[,colSums(SK.Pr) != 0]
SK.Pr.Char[SK.Pr.Char >= 0.8] <- "core"
SK.Pr.Char[SK.Pr.Char < 0.8 & SK.Pr.Char> 0.2] <- "interm"
SK.Pr.Char[SK.Pr.Char<= 0.2] <- "oppor"

OR.Pr.Char <- OR.Pr[,colSums(OR.Pr) != 0]
OR.Pr.Char[OR.Pr.Char >= 0.8] <- "core"
OR.Pr.Char[OR.Pr.Char < 0.8 & OR.Pr.Char> 0.2] <- "interm"
OR.Pr.Char[OR.Pr.Char<= 0.2] <- "oppor"

NS.Pr.Char <- NS.Pr[,colSums(NS.Pr) != 0]
NS.Pr.Char[NS.Pr.Char >= 0.8] <- "core"
NS.Pr.Char[NS.Pr.Char < 0.8 & NS.Pr.Char> 0.2] <- "interm"
NS.Pr.Char[NS.Pr.Char<= 0.2] <- "oppor"

######################################################################

##Stool Prevalance

######################################################################
Genus_ST5 <- subset(Genus_ST, Genus_ST$SubjectID %in% row.names(ST.Pr))
#detect 0 in subsetted data
Genus_ST5[-(1:9)][(colSums(Genus_ST5[-(1:9)])==0)]
#remove unwanted varibles and 0
Genus_ST5.1 <- select(Genus_ST5, -X, -RandomID,-Date,-JAXID,-SampleType,
       -ori.ID, -batch,-colnames(Genus_ST5[-(1:9)][(colSums(Genus_ST5[-(1:9)])==0)]))
#replace "/" with "."
colnames(ST.Pr.Char) <- gsub("/",".",colnames(ST.Pr.Char))
#check if two dataframe is identical 
identical(colnames(Genus_ST5.1)[-(1:2)],colnames(ST.Pr.Char))

ST.pre.list <- unique(Genus_ST5.1$SubjectID)
coreprev.st <- data.frame()
coreprev.st[1,1:3] <- "start"
colnames(coreprev.st) <- c("X1", "X2", "X3")

for (i in 1:length(ST.pre.list)){
  j <- ST.pre.list[i]
  print(j)
x <- subset(Genus_ST5.1, SubjectID == j)
rownames(x) <- x$SampleID
x <- select(x, -SubjectID,-SampleID)
tx <- as_tibble(t(x))
x0 <- subset(ST.Pr.Char,row.names(ST.Pr.Char)==j)
tx$core <- t(x0)
tx.aggre <- select(tx, -core) %>% 
  aggregate(by = list(tx$core), FUN = sum, na.rm=T)
x.aggre <- data.frame(t(tx.aggre))
x.aggre <- filter(x.aggre, X1 != "core")
coreprev.st <- rbind(coreprev.st,x.aggre)
}

colnames(coreprev.st) <- c("core", "interm", "oppor")
coreprev.st.meta <- merge(select(Genus_ST, SampleID, SubjectID, Date), coreprev.st, by.x= "SampleID", by.y="row.names")

#coreprev.st$SubjectID <- sub('^([^-]+-[^-]+).*', '\\1', row.names(coreprev.st))

coreprev.st.meta$Date <- as.Date(coreprev.st.meta$Date, format ="%Y-%m-%d")
coreprev.st.meta$core <- as.numeric(coreprev.st.meta$core)
coreprev.st.meta$interm <- as.numeric(coreprev.st.meta$interm)
coreprev.st.meta$oppor <- as.numeric(coreprev.st.meta$oppor)
sc

coreprev.st.meta$IRIS <- "Unknown"
coreprev.st.meta$IRIS[coreprev.st.meta$SubjectID %in% (filter(sc,IRIS=="IS")$SubjectID)] <- "IS"
coreprev.st.meta$IRIS[coreprev.st.meta$SubjectID %in% (filter(sc,IRIS=="IR")$SubjectID)] <- "IR"

p <- ggplot(coreprev.st.meta, aes(x=Date, y=core, group=SubjectID))  + geom_line(size=0.2) #+ geom_point(size=0.5)
p <- p + geom_line(aes(y=oppor), color="red", size=0.2)+ scale_y_continuous(name = "Core", sec.axis = sec_axis( trans=~.*1, name="Second Axis"))
p <- p + facet_wrap(.~IRIS, scale="free", ncol=1) 
p <- p + ggtitle("Stool_Core") + base_theme + scale_y_continuous(labels = scales::percent, limits=c(0,1))
p
#ggsave("./Analysis/Suppl.figure/Core.Microbiome/Stool.core.opp.pdf",p, height = 2, width = 4, dpi = 300)


p2 <- ggplot(coreprev.st.meta, aes(x=Date, y=interm)) + geom_point() + geom_line()
p2 <- p2 + facet_wrap(.~SubjectID, scale="free")
p2

p3 <- ggplot(coreprev.st.meta, aes(x=Date, y=oppor)) + geom_point() + geom_line()
p3 <- p3 + facet_wrap(.~SubjectID, scale="free_x")
p3

clinic.df1 <- subset(clinic.df, clinic.df$SubjectID %in% rownames(ST.Pr))
clinic.df1$Date <- as.Date(clinic.df1$CollectionDate, format = "%m/%d/%y")
p4 <- ggplot(clinic.df1, aes(x=Date, y=A1C)) + geom_point() + geom_line()
#p4 <- p4 + facet_wrap(.~SubjectID, scale="free")
p4

ck.df1 <- subset(ck.df, ck.df$SubjectID %in% rownames(ST.Pr))
ck.df1$Date <- as.Date(ck.df1$CollectionDate, format = "%m/%d/%y")
p5 <- ggplot(ck.df1, aes(x=Date, y=CHEX4)) + geom_point() + geom_line()
p5 <- p5 + facet_wrap(.~SubjectID, scale="free_x")
p5

# pCHEX1 <- ggplot(ck.df1, aes(x=Date, y=CHEX1)) + geom_point() 
# pCHEX1 <- pCHEX1 + geom_hline(yintercept = mean(ck.df1$CHEX1),color="red")
# pCHEX1 <- pCHEX1 + geom_hline(yintercept = (mean(ck.df1$CHEX1) + 5 * sd(ck.df1$CHEX1)), color="red", linetype="dashed")
# pCHEX1 <- pCHEX1 + geom_hline(yintercept = (mean(ck.df1$CHEX1) - 5 * sd(ck.df1$CHEX1)), color="red", linetype="dashed")
# pCHEX1
# 
# pCHEX2 <- ggplot(ck.df1, aes(x=Date, y=CHEX2)) + geom_point() 
# pCHEX2 <- pCHEX2 + geom_hline(yintercept = mean(ck.df1$CHEX2),color="red")
# pCHEX2 <- pCHEX2 + geom_hline(yintercept = (mean(ck.df1$CHEX2) + 5 * sd(ck.df1$CHEX2)), color="red", linetype="dashed")
# pCHEX2 <- pCHEX2 + geom_hline(yintercept = (mean(ck.df1$CHEX2) - 5 * sd(ck.df1$CHEX2)), color="red", linetype="dashed")
# pCHEX2
# 
# pCHEX3 <- ggplot(ck.df1, aes(x=Date, y=CHEX3)) + geom_point() 
# pCHEX3 <- pCHEX3 + geom_hline(yintercept = mean(ck.df1$CHEX3),color="red")
# pCHEX3 <- pCHEX3 + geom_hline(yintercept = (mean(ck.df1$CHEX3) + 5 * sd(ck.df1$CHEX3)), color="red", linetype="dashed")
# pCHEX3 <- pCHEX3 + geom_hline(yintercept = (mean(ck.df1$CHEX3) - 5 * sd(ck.df1$CHEX3)), color="red", linetype="dashed")
# pCHEX3
# 
# pCHEX4 <- ggplot(ck.df1, aes(x=Date, y=CHEX4)) + geom_point() 
# pCHEX4 <- pCHEX4 + geom_hline(yintercept = mean(ck.df1$CHEX4),color="red")
# pCHEX4 <- pCHEX4 + geom_hline(yintercept = (mean(ck.df1$CHEX4) + 5 * sd(ck.df1$CHEX4)), color="red", linetype="dashed")
# pCHEX4 <- pCHEX4 + geom_hline(yintercept = (mean(ck.df1$CHEX4) - 5 * sd(ck.df1$CHEX4)), color="red", linetype="dashed")
# pCHEX4
# 
# pCHEX1+pCHEX2+pCHEX3+pCHEX4

core.bysubject <- select(coreprev.st.meta,core,interm,oppor) %>%
  aggregate(by=list(coreprev.st.meta$SubjectID), mean)
                          
core.observed <- data.frame(rowSums(ST.Pr>=0.8))
colnames(core.observed) <- "observed"
core.bysubject.meta <- merge(core.bysubject,core.observed, by.x="Group.1", by.y="row.names")
core.bysubject.meta <- merge(core.bysubject.meta, sc, by.x="Group.1",by.y="SubjectID")
colnames(core.bysubject.meta)[1] <- "SubjectID"

core.cv.bysubject <- select(coreprev.st.meta,core,interm,oppor) %>%
  aggregate(by=list(coreprev.st.meta$SubjectID), function(x) mean(x)/sd(x))

colnames(core.cv.bysubject) <- c("SubjectID", "core_stability", "interm_stability", "oppor_stability")
core.bysubject.meta <- merge(core.bysubject.meta, core.cv.bysubject, by="SubjectID")

core.bysubject.meta[core.bysubject.meta$core_stability>100,]

p7 <- ggplot(filter(core.bysubject.meta, IRIS != "Unknown"), aes(y=core_stability, x=observed, color=IRIS)) + geom_point()
p7 <- p7  + geom_smooth(method="lm", se=F) #+ facet_wrap(.~IRIS) #, scale="free_x"
p7
# 
# cor.test(filter(core.bysubject.meta, IRIS != "Unknown")$observed, filter(core.bysubject.meta, IRIS != "Unknown")$SSPG, method = "spearman")
# cor.test(filter(core.bysubject.meta.sk, IRIS != "Unknown")$observed, filter(core.bysubject.meta.sk, IRIS != "Unknown")$SSPG, method = "spearman")
# cor.test(filter(core.bysubject.meta.or, IRIS != "Unknown")$observed, filter(core.bysubject.meta.or, IRIS != "Unknown")$SSPG, method = "spearman")
# cor.test(filter(core.bysubject.meta.ns, IRIS != "Unknown")$observed, filter(core.bysubject.meta.ns, IRIS != "Unknown")$SSPG, method = "spearman")
# 
# cor.test(filter(core.bysubject.meta, IRIS != "Unknown")$observed, filter(core.bysubject.meta, IRIS != "Unknown")$BMI, method = "spearman")
# cor.test(filter(core.bysubject.meta.sk, IRIS != "Unknown")$observed, filter(core.bysubject.meta.sk, IRIS != "Unknown")$BMI, method = "spearman")
# cor.test(filter(core.bysubject.meta.or, IRIS != "Unknown")$observed, filter(core.bysubject.meta.or, IRIS != "Unknown")$BMI, method = "spearman")
# cor.test(filter(core.bysubject.meta.ns, IRIS != "Unknown")$observed, filter(core.bysubject.meta.ns, IRIS != "Unknown")$BMI, method = "spearman")
# 
# p7 <- ggplot(filter(core.bysubject.meta), aes(y=core_stability, x=oppor_stability)) + geom_point()
# p7 <- p7  + geom_smooth(method="lm", se=F) #+ facet_wrap(.~IRIS) #, scale="free_x"
# p7

ST.lm.is <- lm(formula =  core_stability ~ observed + Adj.age + Gender + BMI, data=subset(core.bysubject.meta, IRIS=="IS"))
summary(ST.lm.is)
ST.lm.ir <- lm(formula =  core_stability ~ observed + Adj.age + Gender + BMI, data=subset(core.bysubject.meta, IRIS=="IR"))
summary(ST.lm.ir)

p7 <- ggplot(filter(core.bysubject.meta, IRIS != "Unknown"), aes(y=core, x=BMI, color=IRIS)) + geom_point()
p7 <- p7  + geom_smooth(method="lm", se=F) 
p7

p7 <- ggplot(filter(core.bysubject.meta, IRIS != "Unknown"), aes(y=core, x=IRIS)) + geom_point()
p7 <- p7  + stat_compare_means()
p7

colnames(core.bysubject.meta)

core.bysubject.meta$Th17.Group

core.bysubject.meta.nu <- select(core.bysubject.meta, core,interm,oppor,observed,SSPG,
                                 FPG, Adj.age, BMI,OGTT, core_stability,interm_stability,oppor_stability)

core.bysubject.meta.nu$SSPG <- as.numeric(core.bysubject.meta.nu$SSPG)
core.bysubject.meta.nu$FPG <- as.numeric(core.bysubject.meta.nu$FPG)

res <- rcorr(as.matrix(core.bysubject.meta.nu), type = "pearson")
res

stool.corrematrix <- do.call(cbind, res)
colnames(stool.corrematrix) <- paste(c(rep("B", 12), rep("N", 12), rep("P",12)), colnames(stool.corrematrix), sep="_")
stool.corrematrix

#write.csv(file = "./Analysis/ori_meta_table/Core_microbiome/Core_ST_Precentage_Bysample.csv",coreprev.st.meta)
#write.csv(file = "./Analysis/ori_meta_table/Core_microbiome/Core_ST_Precentage_Bysubject.csv",core.bysubject.meta)


######################################################################

##Saliva Prevalance

######################################################################
Genus_OR5 <- subset(Genus_OR, Genus_OR$SubjectID %in% row.names(OR.Pr))
#detect 0 in subsetted data
Genus_OR5[-(1:7)][(colSums(Genus_OR5[-(1:7)])==0)]
#remove unwanted varibles and 0, leave sample ID and subject ID only
Genus_OR5.1 <- select(Genus_OR5, -X, -KitID,-Date, -SampleType,
                      -batch,-colnames(Genus_OR5[-(1:7)][(colSums(Genus_OR5[-(1:7)])==0)]))
#replace "/" with "."
colnames(OR.Pr.Char) <- gsub("/",".",colnames(OR.Pr.Char))
#replace Unclassified_Actinobacteria.1 with Unclassified_Actinobacteria in OR.Pr.Char
colnames(OR.Pr.Char)[colnames(OR.Pr.Char) == "Unclassified_Actinobacteria.1"] <- "Unclassified_Actinobacteria"
#check if two dataframe is identical, must be ture to excute next step
identical(colnames(Genus_OR5.1)[-(1:2)],colnames(OR.Pr.Char))

OR.pre.list <- unique(Genus_OR5.1$SubjectID)
coreprev.or <- data.frame()
coreprev.or[1,1:3] <- "start"
colnames(coreprev.or) <- c("X1", "X2", "X3")

for (i in 1:length(OR.pre.list)){
  j <- OR.pre.list[i]
  print(j)
  x <- subset(Genus_OR5.1, SubjectID == j)
  rownames(x) <- x$SampleID
  x <- select(x, -SubjectID,-SampleID)
  tx <- as_tibble(t(x))
  x0 <- subset(OR.Pr.Char,row.names(OR.Pr.Char)==j)
  tx$core <- t(x0)
  tx.aggre <- select(tx, -core) %>% 
    aggregate(by = list(tx$core), FUN = sum, na.rm=T)
  x.aggre <- data.frame(t(tx.aggre))
  x.aggre <- filter(x.aggre, X1 != "core")
  coreprev.or <- rbind(coreprev.or,x.aggre)
}

colnames(coreprev.or) <- c("core", "interm", "oppor")
coreprev.or.meta <- merge(select(Genus_OR, SampleID, SubjectID, Date), coreprev.or, by.x= "SampleID", by.y="row.names")
coreprev.or.meta
#coreprev.st$SubjectID <- sub('^([^-]+-[^-]+).*', '\\1', row.names(coreprev.st))

coreprev.or.meta$Date <- as.Date(coreprev.or.meta$Date, format ="%Y-%m-%d")
coreprev.or.meta$core <- as.numeric(coreprev.or.meta$core)
coreprev.or.meta$interm <- as.numeric(coreprev.or.meta$interm)
coreprev.or.meta$oppor <- as.numeric(coreprev.or.meta$oppor)

coreprev.or.meta$IRIS <- "Unknown"
coreprev.or.meta$IRIS[coreprev.or.meta$SubjectID %in% (filter(sc,IRIS=="IS")$SubjectID)] <- "IS"
coreprev.or.meta$IRIS[coreprev.or.meta$SubjectID %in% (filter(sc,IRIS=="IR")$SubjectID)] <- "IR"

p <- ggplot(coreprev.or.meta, aes(x=Date, y=core, group=SubjectID)) + geom_point() + geom_line()
p <- p + geom_line(aes(y=oppor), color="red")+ scale_y_continuous(name = "Core", sec.axis = sec_axis( trans=~.*1, name="Second Axis"))
p <- p + facet_wrap(.~IRIS, scale="free", ncol=1) + ggtitle("Oral_Core")
p

p2 <- ggplot(coreprev.or.meta, aes(x=Date, y=interm)) + geom_point() + geom_line()
p2 <- p2 + facet_wrap(.~SubjectID, scale="free_x")
p2

p3 <- ggplot(coreprev.or.meta, aes(x=Date, y=oppor)) + geom_point() + geom_line()
p3 <- p3 + facet_wrap(.~SubjectID, scale="free_x")
p3

clinic.df1 <- subset(clinic.df, clinic.df$SubjectID %in% rownames(OR.Pr))
clinic.df1$Date <- as.Date(clinic.df1$CollectionDate, format = "%m/%d/%y")
p4 <- ggplot(clinic.df1, aes(x=Date, y=LYM)) + geom_point() + geom_line()
p4 <- p4 + facet_wrap(.~SubjectID, scale="free_x")
p4

ck.df1 <- subset(ck.df, ck.df$SubjectID %in% rownames(OR.Pr))
ck.df1$Date <- as.Date(ck.df1$CollectionDate, format = "%m/%d/%y")
p5 <- ggplot(ck.df1, aes(x=Date, y=CHEX4)) + geom_point() + geom_line()
p5 <- p5 + facet_wrap(.~SubjectID, scale="free_x")
p5

core.bysubject.or <- select(coreprev.or.meta,core,interm,oppor) %>%
  aggregate(by=list(coreprev.or.meta$SubjectID), mean)

core.observed.or <- data.frame(rowSums(OR.Pr>=0.8))
colnames(core.observed.or) <- "observed"
core.bysubject.meta.or <- merge(core.bysubject.or,core.observed.or, by.x="Group.1", by.y="row.names")
core.bysubject.meta.or <- merge(core.bysubject.meta.or, sc, by.x="Group.1",by.y="SubjectID")
colnames(core.bysubject.meta.or)[1] <- "SubjectID"

core.cv.bysubject.or <- select(coreprev.or.meta,core,interm,oppor) %>%
  aggregate(by=list(coreprev.or.meta$SubjectID), function(x) mean(x)/sd(x))

colnames(core.cv.bysubject.or) <- c("SubjectID", "core_stability", "interm_stability", "oppor_stability")
core.bysubject.meta.or <- merge(core.bysubject.meta.or, core.cv.bysubject.or, by="SubjectID")

OR.lm.is <- lm(formula =  core_stability ~ observed, data=subset(core.bysubject.meta.or, IRIS=="IS"))
summary(OR.lm.is)
OR.lm.ir <- lm(formula =  core_stability ~ observed + Adj.age + Gender + BMI , data=subset(core.bysubject.meta.or, IRIS=="IR"))
summary(OR.lm.ir)

p7 <- ggplot(filter(core.bysubject.meta.or, IRIS != "Unknown"), aes(y=observed, x=Adj.age, color=IRIS)) + geom_point()
p7 <- p7  + geom_smooth(method="lm", se=F) 
p7

p7 <- ggplot(filter(core.bysubject.meta.or, IRIS != "IA"), aes(y=core_stability, x=OGTT)) + geom_point()
p7 <- p7  + geom_smooth(method="loess", se=F) 
p7

p7 <- ggplot(filter(core.bysubject.meta.or, IRIS != "Unknown"), aes(y=observed, x=IRIS)) + geom_boxplot()
p7 <- p7  + stat_compare_means(method = "t.test")
p7

colnames(core.bysubject.meta.or)

core.bysubject.meta.or.nu <- select(core.bysubject.meta.or, core,interm,oppor,observed,SSPG,
                                 FPG, Adj.age, BMI,OGTT, core_stability,interm_stability,oppor_stability)

core.bysubject.meta.or.nu$SSPG <- as.numeric(core.bysubject.meta.or.nu$SSPG)
core.bysubject.meta.or.nu$FPG <- as.numeric(core.bysubject.meta.or.nu$FPG)

res <- rcorr(as.matrix(core.bysubject.meta.or.nu),type = "pearson")
corrplot(res$r, type="upper", order="alphabet", p.mat = res$P, sig.level = 0.01, insig = "blank")

#write.csv(file = "./Analysis/ori_meta_table/Core_microbiome/Core_OR_Precentage_Bysample.csv",coreprev.or.meta)
#write.csv(file = "./Analysis/ori_meta_table/Core_microbiome/Core_OR_Precentage_Bysubject.csv",core.bysubject.meta.or)

######################################################################

##Skin Prevalance

######################################################################
Genus_SK5 <- subset(Genus_SK, Genus_SK$SubjectID %in% row.names(SK.Pr))
#detect 0 in subsetted data
Genus_SK5[-(1:7)][(colSums(Genus_SK5[-(1:7)])==0)]
#remove unwanted varibles and 0, leave sample ID and subject ID only
Genus_SK5.1 <- select(Genus_SK5, -X, -KitID,-Date, -SampleType,
                      -batch,-colnames(Genus_SK5[-(1:7)][(colSums(Genus_SK5[-(1:7)])==0)]))
#replace "/" with "."
colnames(SK.Pr.Char) <- gsub("/",".",colnames(SK.Pr.Char))
#replace Unclassified_Actinobacteria.1 with Unclassified_Actinobacteria in SK.Pr.Char
colnames(SK.Pr.Char)[colnames(SK.Pr.Char) == "Unclassified_candidate_division_WPS-1"] <- "Unclassified_candidate_division_WPS.1"
colnames(SK.Pr.Char)[650] <- "Unclassified_Actinobacteria.1"

#check if two dataframe is identical, must be ture to excute next step
identical(colnames(Genus_SK5.1)[-(1:2)],colnames(SK.Pr.Char))

SK.pre.list <- unique(Genus_SK5.1$SubjectID)
coreprev.sk <- data.frame()
coreprev.sk[1,1:3] <- "start"
colnames(coreprev.sk) <- c("X1", "X2", "X3")

for (i in 1:length(SK.pre.list)){
  j <- SK.pre.list[i]
  print(j)
  x <- subset(Genus_SK5.1, SubjectID == j)
  rownames(x) <- x$SampleID
  x <- select(x, -SubjectID,-SampleID)
  tx <- as_tibble(t(x))
  x0 <- subset(SK.Pr.Char,row.names(SK.Pr.Char)==j)
  tx$core <- t(x0)
  tx.aggre <- select(tx, -core) %>% 
    aggregate(by = list(tx$core), FUN = sum, na.rm=T)
  x.aggre <- data.frame(t(tx.aggre))
  x.aggre <- filter(x.aggre, X1 != "core")
  coreprev.sk <- rbind(coreprev.sk,x.aggre)
}

colnames(coreprev.sk) <- c("core", "interm", "oppor")
coreprev.sk.meta <- merge(select(Genus_SK, SampleID, SubjectID, Date), coreprev.sk, by.x= "SampleID", by.y="row.names")

#coreprev.st$SubjectID <- sub('^([^-]+-[^-]+).*', '\\1', row.names(coreprev.st))

coreprev.sk.meta$Date <- as.Date(coreprev.sk.meta$Date, format ="%Y-%m-%d")
coreprev.sk.meta$core <- as.numeric(coreprev.sk.meta$core)
coreprev.sk.meta$interm <- as.numeric(coreprev.sk.meta$interm)
coreprev.sk.meta$oppor <- as.numeric(coreprev.sk.meta$oppor)
p <- ggplot(coreprev.sk.meta, aes(x=Date, y=core)) + geom_point() + geom_line()
p

clinic.df1 <- subset(clinic.df, clinic.df$SubjectID %in% rownames(SK.Pr))
clinic.df1$Date <- as.Date(clinic.df1$CollectionDate, format = "%m/%d/%y")
p4 <- ggplot(clinic.df1, aes(x=Date, y=LYM)) + geom_point() + geom_line()
p4 <- p4 + facet_wrap(.~SubjectID, scale="free_x")
p4

ck.df1 <- subset(ck.df, ck.df$SubjectID %in% rownames(SK.Pr))
ck.df1$Date <- as.Date(ck.df1$CollectionDate, format = "%m/%d/%y")
p5 <- ggplot(ck.df1, aes(x=Date, y=CHEX4)) + geom_point() + geom_line()
p5 <- p5 + facet_wrap(.~SubjectID, scale="free_x")
p5

core.bysubject.sk <- select(coreprev.sk.meta,core,interm,oppor) %>%
  aggregate(by=list(coreprev.sk.meta$SubjectID), mean)

core.observed.sk <- data.frame(rowSums(SK.Pr>=0.8))
colnames(core.observed.sk) <- "observed"
core.bysubject.meta.sk <- merge(core.bysubject.sk,core.observed.sk, by.x="Group.1", by.y="row.names")
core.bysubject.meta.sk <- merge(core.bysubject.meta.sk, sc, by.x="Group.1",by.y="SubjectID")
colnames(core.bysubject.meta.sk)[1] <- "SubjectID"

core.cv.bysubject.sk <- select(coreprev.sk.meta,core,interm,oppor) %>%
  aggregate(by=list(coreprev.sk.meta$SubjectID), function(x) mean(x)/sd(x))

colnames(core.cv.bysubject.sk) <- c("SubjectID", "core_stability", "interm_stability", "oppor_stability")
core.bysubject.meta.sk <- merge(core.bysubject.meta.sk, core.cv.bysubject.sk, by="SubjectID")

p7 <- ggplot(filter(core.bysubject.meta.sk, IRIS != "Unknown"), aes(y=core_stability, x=observed, color=IRIS)) + geom_point()
p7 <- p7  + geom_smooth(method="lm", se=F) + facet_wrap(.~IRIS) #, scale="free_x"
p7

SK.lm.is <- lm(formula =  core_stability ~ observed + Adj.age + Gender + BMI, data=subset(core.bysubject.meta.sk, IRIS=="IS"))
summary(SK.lm.is)
SK.lm.ir <- lm(formula =  core_stability ~ observed + Adj.age + Gender + BMI, data=subset(core.bysubject.meta.sk, IRIS=="IR"))
summary(SK.lm.ir)

p7 <- ggplot(filter(core.bysubject.meta.sk, IRIS != "Unknown"), aes(y=observed, x=Adj.age, color=IRIS)) + geom_point()
p7 <- p7  + geom_smooth(method="lm", se=F) 
p7

p7 <- ggplot(filter(core.bysubject.meta.sk, IRIS != "Unknown"), aes(y=observed, x=IRIS)) + geom_boxplot()
p7 <- p7  + stat_compare_means(method = "t.test")
p7

core.bysubject.meta.sk.nu <- select(core.bysubject.meta.sk, core,interm,oppor,observed,SSPG,
                                    FPG, Adj.age, BMI,OGTT, core_stability,interm_stability,oppor_stability)

core.bysubject.meta.sk.nu$SSPG <- as.numeric(core.bysubject.meta.sk.nu$SSPG)
core.bysubject.meta.sk.nu$FPG <- as.numeric(core.bysubject.meta.sk.nu$FPG)

res <- rcorr(as.matrix(core.bysubject.meta.sk.nu),type = "pearson")
corrplot(res$r, type="upper", order="alphabet", p.mat = res$P, sig.level = 0.01, insig = "blank")

#write.csv(file = "./Analysis/ori_meta_table/Core_microbiome/Core_SK_Precentage_Bysample.csv",coreprev.sk.meta)
#write.csv(file = "./Analysis/ori_meta_table/Core_microbiome/Core_SK_Precentage_Bysubject.csv",core.bysubject.meta.sk)

######################################################################

##Nasal Prevalance

######################################################################
Genus_NS5 <- subset(Genus_NS, Genus_NS$SubjectID %in% row.names(NS.Pr))
#detect 0 in subsetted data
Genus_NS5[-(1:9)][(colSums(Genus_NS5[-(1:9)])==0)]
#remove unwanted varibles and 0, leave sample ID and subject ID only
Genus_NS5.1 <- select(Genus_NS5, -X, -RandomID,-Date,-JAXID,-SampleType,
                      -ori.ID, -batch,-colnames(Genus_NS5[-(1:9)][(colSums(Genus_NS5[-(1:9)])==0)]))
#replace "/" with "."
colnames(NS.Pr.Char) <- gsub("/",".",colnames(NS.Pr.Char))
#replace Unclassified_Actinobacteria.1 with Unclassified_Actinobacteria in NS.Pr.Char
colnames(NS.Pr.Char)[colnames(NS.Pr.Char) == "Unclassified_candidate_division_WPS-1"] <- "Unclassified_candidate_division_WPS.1"
colnames(NS.Pr.Char)[650] <- "Unclassified_Actinobacteria.1"

#check if two dataframe is identical, must be ture to excute next step
identical(colnames(Genus_NS5.1)[-(1:2)],colnames(NS.Pr.Char))

NS.pre.list <- unique(Genus_NS5.1$SubjectID)
coreprev.ns <- data.frame()
coreprev.ns[1,1:3] <- "start"
colnames(coreprev.ns) <- c("X1", "X2", "X3")

for (i in 1:length(NS.pre.list)){
  j <- NS.pre.list[i]
  print(j)
  x <- subset(Genus_NS5.1, SubjectID == j)
  rownames(x) <- x$SampleID
  x <- select(x, -SubjectID,-SampleID)
  tx <- as_tibble(t(x))
  x0 <- subset(NS.Pr.Char,row.names(NS.Pr.Char)==j)
  tx$core <- t(x0)
  tx.aggre <- select(tx, -core) %>% 
    aggregate(by = list(tx$core), FUN = sum, na.rm=T)
  x.aggre <- data.frame(t(tx.aggre))
  x.aggre <- filter(x.aggre, X1 != "core")
  coreprev.ns <- rbind(coreprev.ns,x.aggre)
}

colnames(coreprev.ns) <- c("core", "interm", "oppor")
coreprev.ns.meta <- merge(select(Genus_NS, SampleID, SubjectID, Date), coreprev.ns, by.x= "SampleID", by.y="row.names")

coreprev.ns.meta$Date <- as.Date(coreprev.ns.meta$Date, format ="%Y-%m-%d")
coreprev.ns.meta$core <- as.numeric(coreprev.ns.meta$core)
coreprev.ns.meta$interm <- as.numeric(coreprev.ns.meta$interm)
coreprev.ns.meta$oppor <- as.numeric(coreprev.ns.meta$oppor)

clinic.df1 <- subset(clinic.df, clinic.df$SubjectID %in% rownames(NS.Pr))
clinic.df1$Date <- as.Date(clinic.df1$CollectionDate, format = "%m/%d/%y")

ck.df1 <- subset(ck.df, ck.df$SubjectID %in% rownames(NS.Pr))
ck.df1$Date <- as.Date(ck.df1$CollectionDate, format = "%m/%d/%y")

core.bysubject.ns <- select(coreprev.ns.meta,core,interm,oppor) %>%
  aggregate(by=list(coreprev.ns.meta$SubjectID), mean)

core.observed.ns <- data.frame(rowSums(NS.Pr>=0.8))
colnames(core.observed.ns) <- "observed"
core.bysubject.meta.ns <- merge(core.bysubject.ns,core.observed.ns, by.x="Group.1", by.y="row.names")
core.bysubject.meta.ns <- merge(core.bysubject.meta.ns, sc, by.x="Group.1",by.y="SubjectID")
colnames(core.bysubject.meta.ns)[1] <- "SubjectID"

core.cv.bysubject.ns <- select(coreprev.ns.meta,core,interm,oppor) %>%
  aggregate(by=list(coreprev.ns.meta$SubjectID), function(x) mean(x)/sd(x))

colnames(core.cv.bysubject.ns) <- c("SubjectID", "core_stability", "interm_stability", "oppor_stability")
core.bysubject.meta.ns <- merge(core.bysubject.meta.ns, core.cv.bysubject.ns, by="SubjectID")

NS.lm.is <- lm(formula =  core_stability ~ observed + Adj.age + Gender, data=subset(core.bysubject.meta.ns, IRIS=="IS"))
summary(NS.lm.is)
NS.lm.ir <- lm(formula =  core_stability ~ observed + Adj.age + Gender, data=subset(core.bysubject.meta.ns, IRIS=="IR"))
summary(NS.lm.ir)

core.bysubject.meta.ns.nu <- select(core.bysubject.meta.ns, core,interm,oppor,observed,SSPG,
                                    FPG, Adj.age, BMI,OGTT, core_stability,interm_stability,oppor_stability)

core.bysubject.meta.ns.nu$SSPG <- as.numeric(core.bysubject.meta.ns.nu$SSPG)
core.bysubject.meta.ns.nu$FPG <- as.numeric(core.bysubject.meta.ns.nu$FPG)

res <- rcorr(as.matrix(core.bysubject.meta.ns.nu),type = "pearson")
corrplot(res$r, type="upper", order="alphabet", p.mat = res$P, sig.level = 0.01, insig = "blank")

#write.csv(file = "./Analysis/ori_meta_table/Core_microbiome/Core_NS_Precentage_Bysample.csv",coreprev.ns.meta)
#write.csv(file = "./Analysis/ori_meta_table/Core_microbiome/Core_NS_Precentage_Bysubject.csv",core.bysubject.meta.ns)

################################################################################

#Overall Characteristics

################################################################################
dim(ST.Pr)
dim(SK.Pr)
dim(OR.Pr)
dim(NS.Pr)

ST.Pr = ST.Pr[colSums(ST.Pr) != 0]
SK.Pr = SK.Pr[colSums(SK.Pr) != 0]
OR.Pr = OR.Pr[colSums(OR.Pr) != 0]
NS.Pr = NS.Pr[colSums(NS.Pr) != 0]

median.ST.core0 <- lapply(ST.Pr, median)
median.ST.core <- as.data.frame(do.call(rbind, median.ST.core0)) %>% arrange(desc(V1))

median.SK.core0 <- lapply(SK.Pr, median)
median.SK.core <- as.data.frame(do.call(rbind, median.SK.core0)) %>% arrange(desc(V1))

median.OR.core0 <- lapply(OR.Pr, median)
median.OR.core <- as.data.frame(do.call(rbind, median.OR.core0)) %>% arrange(desc(V1))

median.NS.core0 <- lapply(NS.Pr, median)
median.NS.core <- as.data.frame(do.call(rbind, median.NS.core0)) %>% arrange(desc(V1))

sort(median.NS.core$V1, decreasing = T)

#rank by median value and take the first 100 for rank prevalance curve
ST.Pr.100 <- ST.Pr[,rownames(median.ST.core)[1:100]]
SK.Pr.100 <- SK.Pr[,rownames(median.SK.core)[1:100]]
OR.Pr.100 <- OR.Pr[,rownames(median.OR.core)[1:100]]
NS.Pr.100 <- NS.Pr[,rownames(median.NS.core)[1:100]]

smooth.data.pr <- cbind((median.ST.core$V1)[1:100], 
                        (median.SK.core$V1)[1:100], 
                        (median.OR.core$V1)[1:100],
                        (median.NS.core$V1)[1:100])
colnames(smooth.data.pr) <- c("Stool", "Skin", "Oral", "Nasal")
rownames(smooth.data.pr) <- 1:100

p.smooth <- melt(smooth.data.pr) %>% ggplot(aes(x=Var1, y=value, group=Var2)) + geom_smooth(aes(color=Var2))
p.smooth <- p.smooth + base_theme + scale_color_manual(values = body_site_color) + xlab("Rank") + ylab("Prevalence") + scale_y_continuous(labels = scales::percent, limits=c(0,1))
p.smooth <- p.smooth + ggtitle("Rank Prevalence Curve of the microbiome on four body sites")
p.smooth

# write.csv(file = "./Analysis/tables/Rank_Prevalen_Curve.csv",smooth.data.pr)
# ggsave2(filename = "./Analysis/Suppl.figure/RankPrevaCurve.pdf", p.smooth, width = 6, height = 2.5, dpi = 300)

#we can also sort by mean value instead of median value
# par(mfrow=c(4,1))
# ST.Pr.100 <- ST.Pr[,sort(colSums(ST.Pr),decreasing=T)[1:100] %>% names()]
# SK.Pr.100 <- SK.Pr[,sort(colSums(SK.Pr),decreasing=T)[1:100] %>% names()]
# OR.Pr.100 <- OR.Pr[,sort(colSums(OR.Pr),decreasing=T)[1:100] %>% names()]
# NS.Pr.100 <- NS.Pr[,sort(colSums(NS.Pr),decreasing=T)[1:100] %>% names()]

p.st.pre <- melt(ST.Pr.100) %>% mutate(variable=factor(variable, levels=rev(colnames(ST.Pr.100)))) %>% 
  ggplot(aes(y=variable, x=value)) + ggtitle("Stool") + geom_boxplot() + geom_vline(xintercept = 0.8, color="blue") +
  geom_vline(xintercept = 0.2, color="red") + scale_x_reverse() + base_theme   
p.sk.pre <- melt(SK.Pr.100) %>% mutate(variable=factor(variable, levels=rev(colnames(SK.Pr.100)))) %>% 
  ggplot(aes(y=variable, x=value)) + ggtitle("Skin") + geom_boxplot() + geom_vline(xintercept = 0.8, color="blue") +
  geom_vline(xintercept = 0.2, color="red") + scale_x_reverse() + base_theme 
p.or.pre <- melt(OR.Pr.100) %>% mutate(variable=factor(variable, levels=rev(colnames(OR.Pr.100)))) %>% 
  ggplot(aes(y=variable, x=value)) + ggtitle("Oral") + geom_boxplot() + geom_vline(xintercept = 0.8, color="blue") +
  geom_vline(xintercept = 0.2, color="red") + scale_x_reverse() + base_theme  
p.ns.pre <- melt(NS.Pr.100) %>% mutate(variable=factor(variable, levels=rev(colnames(NS.Pr.100)))) %>% 
  ggplot(aes(y=variable, x=value)) + ggtitle("Nasal") + geom_boxplot() + geom_vline(xintercept = 0.8, color="blue") +
  geom_vline(xintercept = 0.2, color="red") + scale_x_reverse() + base_theme 

table(colMeans(ST.Pr.100) > 0.8)
table(colMeans(SK.Pr.100) > 0.8)
table(colMeans(OR.Pr.100) > 0.8)
table(colMeans(NS.Pr.100) > 0.8)

table(colMeans(ST.Pr.100) < 0.8 & colMeans(ST.Pr.100) > 0.2)
table(colMeans(SK.Pr.100) < 0.8 & colMeans(SK.Pr.100) > 0.2)
table(colMeans(OR.Pr.100) < 0.8 & colMeans(OR.Pr.100) > 0.2)
table(colMeans(NS.Pr.100) < 0.8 & colMeans(NS.Pr.100) > 0.2)

# ggsave(filename = "./Analysis/Suppl.figure/Pre100.ST.pdf",p.st.pre, height = 13.5, width = 6, dpi=300)
# ggsave(filename = "./Analysis/Suppl.figure/Pre100.SK.pdf",p.sk.pre, height = 13.5, width = 6, dpi=300)
# ggsave(filename = "./Analysis/Suppl.figure/Pre100.OR.pdf",p.or.pre, height = 13.5, width = 6, dpi=300)
# ggsave(filename = "./Analysis/Suppl.figure/Pre100.NS.pdf",p.ns.pre, height = 13.5, width = 6, dpi=300)
# 
# write.csv(file = "./Analysis/tables/ST.pr.csv",ST.Pr)
# write.csv(file = "./Analysis/tables/SK.pr.csv",SK.Pr)
# write.csv(file = "./Analysis/tables/OR.pr.csv",OR.Pr)
# write.csv(file = "./Analysis/tables/NS.pr.csv",NS.Pr)

######################################################
#####separate the rank prev curve by iris#############
ST.Pr.IS <- filter(ST.Pr,rownames(ST.Pr) %in% filter(sc, IRIS=="IS")$SubjectID)
ST.Pr.IR <- filter(ST.Pr,rownames(ST.Pr) %in% filter(sc, IRIS=="IR")$SubjectID)

SK.Pr.IS <- filter(SK.Pr,rownames(SK.Pr) %in% filter(sc, IRIS=="IS")$SubjectID)
SK.Pr.IR <- filter(SK.Pr,rownames(SK.Pr) %in% filter(sc, IRIS=="IR")$SubjectID)

OR.Pr.IS <- filter(OR.Pr,rownames(OR.Pr) %in% filter(sc, IRIS=="IS")$SubjectID)
OR.Pr.IR <- filter(OR.Pr,rownames(OR.Pr) %in% filter(sc, IRIS=="IR")$SubjectID)

NS.Pr.IS <- filter(NS.Pr,rownames(NS.Pr) %in% filter(sc, IRIS=="IS")$SubjectID)
NS.Pr.IR <- filter(NS.Pr,rownames(NS.Pr) %in% filter(sc, IRIS=="IR")$SubjectID)

median.ST.IS.core0 <- lapply(ST.Pr.IS, median)
median.ST.IS.core <- as.data.frame(do.call(rbind, median.ST.IS.core0)) %>% arrange(desc(V1))
median.ST.IR.core0 <- lapply(ST.Pr.IR, median)
median.ST.IR.core <- as.data.frame(do.call(rbind, median.ST.IR.core0)) %>% arrange(desc(V1))

median.SK.IS.core0 <- lapply(SK.Pr.IS, median)
median.SK.IS.core <- as.data.frame(do.call(rbind, median.SK.IS.core0)) %>% arrange(desc(V1))
median.SK.IR.core0 <- lapply(SK.Pr.IR, median)
median.SK.IR.core <- as.data.frame(do.call(rbind, median.SK.IR.core0)) %>% arrange(desc(V1))

median.OR.IS.core0 <- lapply(OR.Pr.IS, median)
median.OR.IS.core <- as.data.frame(do.call(rbind, median.OR.IS.core0)) %>% arrange(desc(V1))
median.OR.IR.core0 <- lapply(OR.Pr.IR, median)
median.OR.IR.core <- as.data.frame(do.call(rbind, median.OR.IR.core0)) %>% arrange(desc(V1))

median.NS.IS.core0 <- lapply(NS.Pr.IS, median)
median.NS.IS.core <- as.data.frame(do.call(rbind, median.NS.IS.core0)) %>% arrange(desc(V1))
median.NS.IR.core0 <- lapply(NS.Pr.IR, median)
median.NS.IR.core <- as.data.frame(do.call(rbind, median.NS.IR.core0)) %>% arrange(desc(V1))

ST.Pr.IS.100 <- ST.Pr.IS[,rownames(median.ST.IS.core)[1:100]]
ST.Pr.IR.100 <- ST.Pr.IR[,rownames(median.ST.IR.core)[1:100]]

SK.Pr.IS.100 <- SK.Pr.IS[,rownames(median.SK.IS.core)[1:100]]
SK.Pr.IR.100 <- SK.Pr.IR[,rownames(median.SK.IR.core)[1:100]]

OR.Pr.IS.100 <- OR.Pr.IS[,rownames(median.OR.IS.core)[1:100]]
OR.Pr.IR.100 <- OR.Pr.IR[,rownames(median.OR.IR.core)[1:100]]

NS.Pr.IS.100 <- NS.Pr.IS[,rownames(median.NS.IS.core)[1:100]]
NS.Pr.IR.100 <- NS.Pr.IR[,rownames(median.NS.IR.core)[1:100]]


smooth.data.pr.IRIS <- cbind((median.ST.IS.core$V1)[1:100],
                             (median.ST.IR.core$V1)[1:100], 
                             (median.SK.IS.core$V1)[1:100], 
                             (median.SK.IR.core$V1)[1:100], 
                             (median.OR.IS.core$V1)[1:100],
                             (median.OR.IR.core$V1)[1:100],
                             (median.NS.IS.core$V1)[1:100],
                             (median.NS.IR.core$V1)[1:100])

colnames(smooth.data.pr.IRIS) <- c("Stool_IS","Stool_IR", "Skin_IS", "Skin_IR", "Oral_IS","Oral_IR", "Nasal_IS","Nasal_IR")
rownames(smooth.data.pr.IRIS) <- 1:100

body_site_color_IRIS <- c("#DF8F44FF", "#F1deb6","#00A1D5FF","#B6edf1", "#B24745FF","#F1c8b6", "#79AF97FF", "#B6f1b7")
names(body_site_color_IRIS) <- c("Stool_IS","Stool_IR", "Skin_IS", "Skin_IR", "Oral_IS","Oral_IR", "Nasal_IS","Nasal_IR")

p.smooth.iris <- melt(smooth.data.pr.IRIS) %>% ggplot(aes(x=Var1, y=value, group=Var2)) + geom_smooth(aes(color=Var2), se = 0)
p.smooth.iris <- p.smooth.iris + base_theme + scale_color_manual(values = body_site_color_IRIS) + xlab("Rank") + ylab("Prevalence") + scale_y_continuous(labels = scales::percent, limits=c(0,1))
p.smooth.iris <- p.smooth.iris + ggtitle("Rank Prevalence Curve of the microbiome on four body sites (IRIS seperated)")
p.smooth.iris

#ggsave(filename = "./Analysis/Suppl.figure/RankPrevaCurve_IRIS.pdf", p.smooth.iris, width = 6, height = 2.5, dpi = 300)


genus.mean.st <- aggregate(Genus_ST5.1, by= list(Genus_ST5.1$SubjectID), mean) %>% select(-SubjectID, -SampleID) %>% melt()
colnames(genus.mean.st) <- c("Var1", "Var2","Prev")
genus.mean.sk <- aggregate(Genus_SK5.1, by= list(Genus_SK5.1$SubjectID), mean) %>% select(-SubjectID, -SampleID) %>% melt()
colnames(genus.mean.sk) <- c("Var1", "Var2","Prev")
genus.mean.or <- aggregate(Genus_OR5.1, by= list(Genus_OR5.1$SubjectID), mean) %>% select(-SubjectID, -SampleID) %>% melt()
colnames(genus.mean.or) <- c("Var1", "Var2","Prev")
genus.mean.ns <- aggregate(Genus_NS5.1, by= list(Genus_NS5.1$SubjectID), mean) %>% select(-SubjectID, -SampleID) %>% melt()
colnames(genus.mean.ns) <- c("Var1", "Var2","Prev")

stool.pre.ab <- ST.Pr %>% as.matrix() %>% melt() %>% merge(genus.mean.st, by=c("Var1", "Var2"))
skin.pre.ab <- SK.Pr %>% as.matrix() %>% melt() %>% merge(genus.mean.sk, by=c("Var1", "Var2"))
oral.pre.ab <- OR.Pr %>% as.matrix() %>% melt() %>% merge(genus.mean.or, by=c("Var1", "Var2"))
nasal.pre.ab <- NS.Pr %>% as.matrix() %>% melt() %>% merge(genus.mean.ns, by=c("Var1", "Var2"))

p.st.pre.ab <- ggplot(stool.pre.ab, aes(x=value, y=Prev)) + geom_point(color="#DF8F44FF",size=0.6)  + base_theme + scale_y_log10() + ggtitle("Stool Microbiome Ecology")
p.st.pre.ab <- p.st.pre.ab + ylab("Relative Abundance (longitudinal mean)") + xlab("Longitudinal Prevalence") + geom_hline(yintercept=0.001, linetype=3)
p.st.pre.ab

p.sk.pre.ab <- ggplot(skin.pre.ab, aes(x=value, y=Prev)) + geom_point(color="#00A1D5FF",size=0.6)  + base_theme + scale_y_log10() + ggtitle("Skin Microbiome Ecology")
p.sk.pre.ab <- p.sk.pre.ab + ylab("Relative Abundance (longitudinal mean)") + xlab("Longitudinal Prevalence")+ geom_hline(yintercept=0.001, linetype=3)
p.sk.pre.ab

p.or.pre.ab <- ggplot(oral.pre.ab, aes(x=value, y=Prev)) + geom_point(color="#B24745FF",size=0.6)  + base_theme + scale_y_log10() + ggtitle("Oral Microbiome Ecology")
p.or.pre.ab <- p.or.pre.ab + ylab("Relative Abundance (longitudinal mean)") + xlab("Longitudinal Prevalence")+ geom_hline(yintercept=0.001, linetype=3)
p.or.pre.ab

p.ns.pre.ab <- ggplot(nasal.pre.ab, aes(x=value, y=Prev)) + geom_point(color="#79AF97FF",size=0.6)  + base_theme + scale_y_log10() + ggtitle("Nasal Microbiome Ecology")
p.ns.pre.ab <- p.ns.pre.ab + ylab("Relative Abundance (longitudinal mean)") + xlab("Longitudinal Prevalence")+ geom_hline(yintercept=0.001, linetype=3)
p.ns.pre.ab

p.pre.ab <- p.st.pre.ab + p.sk.pre.ab + p.or.pre.ab + p.ns.pre.ab
p.pre.ab
#ggsave(filename = "./Analysis/Suppl.figure/Pre_Abun.pdf", p.pre.ab, width = 12, height = 8, dpi = 300)

p.st.pre.ab.iris <- merge(stool.pre.ab, sc, by.x=, "Var1", by.y="SubjectID", all.x = T) %>% filter(IRIS != "Unknown") %>% 
  ggplot(aes(x=value, y=Prev, color=IRIS))  + base_theme + scale_y_log10() + ggtitle("Stool Microbiome Ecology") + scale_color_manual(values = iris_color)
p.st.pre.ab.iris <- p.st.pre.ab.iris + geom_smooth(method = "loess") + ylab("Relative Abundance (longitudinal mean)") + xlab("Longitudinal Prevalence") + geom_hline(yintercept=0.001, linetype=3)
p.st.pre.ab.iris

p.sk.pre.ab.iris <- merge(skin.pre.ab, sc, by.x=, "Var1", by.y="SubjectID", all.x = T) %>% filter(IRIS != "Unknown") %>% 
  ggplot(aes(x=value, y=Prev, color=IRIS))  + base_theme + scale_y_log10() + ggtitle("Skin Microbiome Ecology") + scale_color_manual(values = iris_color)
p.sk.pre.ab.iris <- p.sk.pre.ab.iris + geom_smooth(method = "loess") + ylab("Relative Abundance (longitudinal mean)") + xlab("Longitudinal Prevalence") + geom_hline(yintercept=0.001, linetype=3)
p.sk.pre.ab.iris

p.or.pre.ab.iris <- merge(oral.pre.ab, sc, by.x=, "Var1", by.y="SubjectID", all.x = T) %>% filter(IRIS != "Unknown") %>% 
  ggplot(aes(x=value, y=Prev, color=IRIS))  + base_theme + scale_y_log10() + ggtitle("Oral Microbiome Ecology") + scale_color_manual(values = iris_color)
p.or.pre.ab.iris <- p.or.pre.ab.iris + geom_smooth(method = "loess") + ylab("Relative Abundance (longitudinal mean)") + xlab("Longitudinal Prevalence") + geom_hline(yintercept=0.001, linetype=3)
p.or.pre.ab.iris

p.ns.pre.ab.iris <- merge(nasal.pre.ab, sc, by.x=, "Var1", by.y="SubjectID", all.x = T) %>% filter(IRIS != "Unknown") %>% 
  ggplot(aes(x=value, y=Prev, color=IRIS)) + base_theme + scale_y_log10() + ggtitle("Nasal Microbiome Ecology") + scale_color_manual(values = iris_color)
p.ns.pre.ab.iris <- p.ns.pre.ab.iris + geom_smooth(method = "loess") + ylab("Relative Abundance (longitudinal mean)") + xlab("Longitudinal Prevalence") + geom_hline(yintercept=0.001, linetype=3)
p.ns.pre.ab.iris

p.pre.ab.iris <- p.st.pre.ab.iris + p.sk.pre.ab.iris + p.or.pre.ab.iris + p.ns.pre.ab.iris
p.pre.ab.iris

ggsave(filename = "./Analysis/Suppl.figure/Pre_Abun_IRIS_nopoint.pdf", p.pre.ab.iris, width = 12, height = 8, dpi = 300)

stool.pre.ab %>%  ggplot(aes(x=value, y=Prev)) + geom_point()

stool.pre.ab.test <- stool.pre.ab %>% group_by(Var2) %>% mutate(meanvalue=mean(value)) %>% mutate(meanprev = mean(Prev)) %>% select(Var2,meanvalue,meanprev) %>% unique()
stool.pre.ab.test
write.csv(file="~/Desktop/stool.pre.ab.test.csv",stool.pre.ab.test)


