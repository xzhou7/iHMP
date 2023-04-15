#Cytokine Analysis
#Author: Xin Zhou

library(stringr)
library(ggplot2)
library(cowplot)
library(patchwork)
library(tidyverse)
library(dplyr)

setwd("~/Box/XinZhouFiles/Projects/ZXP2_HMP2_Core/")
Cytokine.df <- read.csv("./Analysis/ori_meta_table/Cytokinetable/00.HMP_Cytokine_value.csv", header = T)
standard <- read.csv("./Analysis/ori_meta_table/Cytokinetable/01.cytoinestandardcurve.csv",header = T)
load("./Analysis/Robject/Revision_MultiOmes_0509.RData")
infectionstate <- read.csv("./Analysis/ori_meta_table/hmp2_infection_metadata_0323.csv", header = T)
sampleFreq <- read.csv("./Analysis/ori_meta_table/SampleFrequency.csv",header = T)


#########################################################################################################
#calculate the percent of 0 in pg/ul data
#########################################################################################################
setdiff(ck.df$SampleID, infectionstate$SampleID)

table(str_detect(ck.df$SampleID, "69-028"))

sum(str_count(Cytokine.df$MIP1A, pattern = "<"))
sum(str_count(Cytokine.df$MIP1A, pattern = "<")==0)

percent.df <- data.frame()

#count how many measurement are in format of < #### value
for (i in 1:62){
  percent.df[i,1] <- sum(str_count(Cytokine.df[,i+2], pattern = "<"))
}

percent.df$cytokine <- colnames(Cytokine.df)[3:64]
percent.df$total <- 1083
percent.df$ULOD <-  percent.df$V1/percent.df$total

percent.df$cytokine <- factor(percent.df$cytokine, levels = percent.df$cytokine[order(percent.df[4])])
p <- ggplot(percent.df[order(percent.df[4]),], aes(y=cytokine, x=ULOD)) + geom_bar(stat="identity")
p <- p + theme(axis.text.x=element_text(angle=90,hjust=1)) + labs(x="Percent of Samples Under LOD") + theme_cowplot()
p

Data <- as.data.frame(apply(Cytokine.df[,(3:64)], 2, as.numeric))
Data[is.na(Data)] <- 0

var_nonzero <- function(x) var(x[!x == 0])
sd_nonzero <- function(x) sd(x[!x == 0])
mean_nonzero <- function(x) mean(x[!x == 0])

ck.variance <- as.data.frame(apply(Data, 2, var_nonzero))
colnames(ck.variance) <- "variance"

percent.df1 <- merge(percent.df, ck.variance, by.x="cytokine", by.y="row.names")

p <- ggplot(percent.df1, aes(x=log(variance), y=ULOD)) + geom_point()
p

#########################################################################################################
#Plot CHEX 
#########################################################################################################

standard_0 <- filter(standard, Type=="MFI" & Standard == "SB")

colMeans(standard_0[-(1:3)])

llod <- data.frame(sapply(standard_0[-(1:3)], mean, simplify = TRUE))
llod[2] <-  sapply(standard_0[-(1:3)], sd, simplify = TRUE)
colnames(llod) <- c("mean", "sd")
llod$llod.mfi <- llod$mean + 3 * llod$sd 

ck.df1 <- ck.df
ck.df1$Date <- as.Date(ck.df$CollectionDate, format = "%m/%d/%y")
#######################################################################################################################################
pCHEX1 <- ggplot(ck.df1, aes(x=Date, y=CHEX1, color=Plate)) + geom_point() 
pCHEX1 <- pCHEX1 + geom_hline(yintercept = mean(ck.df1$CHEX1),color="red")
pCHEX1 <- pCHEX1 + geom_hline(yintercept = (mean(ck.df1$CHEX1) + 5 * sd(ck.df1$CHEX1)), color="red", linetype="dotted")
pCHEX1 <- pCHEX1 + geom_hline(yintercept = (mean(ck.df1$CHEX1) - 5 * sd(ck.df1$CHEX1)), color="red", linetype="dotted")
pCHEX1 <- pCHEX1 + geom_hline(yintercept = 10071, color="blue", linetype="dashed") + ggtitle("instrument performance")
pCHEX1

pCHEX2 <- ggplot(ck.df1, aes(x=Date, y=CHEX2, color=Plate)) + geom_point() 
pCHEX2 <- pCHEX2 + geom_hline(yintercept = mean(ck.df1$CHEX2),color="red")
pCHEX2 <- pCHEX2 + geom_hline(yintercept = (mean(ck.df1$CHEX2) + 5 * sd(ck.df1$CHEX2)), color="red", linetype="dotted")
pCHEX2 <- pCHEX2 + geom_hline(yintercept = (mean(ck.df1$CHEX2) - 5 * sd(ck.df1$CHEX2)), color="red", linetype="dotted")
pCHEX2 <- pCHEX2 + geom_hline(yintercept = 9076, color="blue", linetype="dashed")+ ggtitle("detection antibody")
pCHEX2

pCHEX3 <- ggplot(ck.df1, aes(x=Date, y=CHEX3, color=Plate)) + geom_point() 
pCHEX3 <- pCHEX3 + geom_hline(yintercept = mean(ck.df1$CHEX3),color="red")
pCHEX3 <- pCHEX3 + geom_hline(yintercept = (mean(ck.df1$CHEX3) + 5 * sd(ck.df1$CHEX3)), color="red", linetype="dotted")
pCHEX3 <- pCHEX3 + geom_hline(yintercept = (mean(ck.df1$CHEX3) - 5 * sd(ck.df1$CHEX3)), color="red", linetype="dotted")
pCHEX3 <- pCHEX3 + geom_hline(yintercept = 2862, color="blue", linetype="dashed")+ ggtitle("fluorescent reporter")
pCHEX3

pCHEX4 <- ggplot(ck.df1, aes(x=Date, y=CHEX4, color=Plate)) + geom_point() 
pCHEX4 <- pCHEX4 + geom_hline(yintercept = mean(ck.df1$CHEX4),color="red")
pCHEX4 <- pCHEX4 + geom_hline(yintercept = (mean(ck.df1$CHEX4) + 5 * sd(ck.df1$CHEX4)), color="red", linetype="dotted")
pCHEX4 <- pCHEX4 + geom_hline(yintercept = (mean(ck.df1$CHEX4) - 5 * sd(ck.df1$CHEX4)), color="red", linetype="dotted")
pCHEX4 <- pCHEX4 + geom_hline(yintercept = 14, color="blue", linetype="dashed")+ ggtitle("nonspecific binding")
pCHEX4

pCHEX1 + pCHEX2 + pCHEX3 + pCHEX4 + plot_layout(guides = "collect")  & theme(legend.position = 'bottom')

#######################################################################################################################################

ck.df0 <- subset(ck.df, CHEX1 != 0 & CHEX4 < 50)
ck.df0$SubjectID[ck.df0$SubjectID == "?69-106"] <- "69-106"
#get cytokine geommean for 965 samples
op.df <- ck.df0[,4:70]
peoplelist <- unique(op.df[,"SubjectID"])
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
z1 <- data.frame()
for (j in 1:length(peoplelist)){
  i <- peoplelist[j]
  print(i)
  x <- subset(op.df, SubjectID==i)
  y <- x[,1:66]
  y1 <- round(y,digits=3)
  y2 <- lapply(y1, gm_mean)
  z <- as.data.frame(y2)
  z1 <- rbind(z,z1)
}
row.names(z1) <- rev(peoplelist)
rm(i,j,x,y,y1,y2,z)

z1

op.df[op.df$SubjectID == "69-005",]

xy <- aggregate(op.df, by = list(op.df$SubjectID),FUN=var)
z2 <- table(ck.df$SubjectID)
z2 <- as.data.frame(z2)
var0 <- merge(z2, xy, by.x = "Var1", by.y <- "Group.1")
row.names(var0) <- var0$Var1
rm(xy,z2)

colnames(var0)[2] <- "NumberofSample"
colnames(var0)[3:68] <- paste(colnames(var0)[3:68],"Var", sep="_")
#write.csv(file = "./Analysis/ori_meta_table/llod.cytokine.csv", llod)
#write.csv(file="./Analysis/ori_meta_table/Cytokinetable/Cytokine.Bysample.csv", ck.df0)
#write.csv(file = "./Analysis/ori_meta_table/Cytokinetable/Cytokine.BySubject.csv", z1)
#write.csv(file = "./Analysis/ori_meta_table/Cytokinetable/Cytokine.Var.csv", var0)

ck.df0 <- subset(ck.df0, ck.df0$SubjectID %in% sampleFreq$SubjectID)

op.df1 <- ck.df0[,c(1,4:71, 75)]
op.df1$CollectionDate <- as.Date(op.df1$CollectionDate, format = "%m/%d/%y")

#calculate the z score based on MFI value (none-log transformed)
dat = op.df1 %>% 
  gather(variable, value,-SampleID, -SubjectID,-CollectionDate,-CL4) %>%
  group_by(SubjectID, variable) %>% 
  mutate(z_score_group = (value - mean(value)) / sd(value)) %>% 
  mutate(z_score_group2 = scale(value)) %>%
  mutate(Z_mean = mean(z_score_group)) %>% 
  ungroup %>% 
  mutate(z_score_ungrouped = (value - mean(value)) / sd(value)) %>%
  mutate(Z_sd = sd(z_score_ungrouped)) 

p.zcytokine <- filter(dat, variable == "IP10") %>% ggplot(aes(x=CollectionDate, y= value, color=CL4, group=SubjectID)) +geom_line() + geom_point()
p.zcytokine <- p.zcytokine + facet_wrap(.~SubjectID, scale="free") #+ geom_hline(yintercept = 2)
p.zcytokine

#calculate the z score based on MFI value (log2 transformed)
op.df2 <- op.df1
op.df2[2:67] <- log(op.df2[2:67])

dat2 = op.df2 %>% 
  gather(variable, value,-SampleID, -SubjectID,-CollectionDate,-CL4) %>%
  group_by(SubjectID, variable) %>% 
  mutate(z_score_group = (value - mean(value)) / sd(value)) %>% 
  ungroup %>% 
  mutate(z_score_ungrouped = (value - mean(value)) / sd(value)) 

p.zcytokine_log <- filter(dat2, variable == "IP10") %>% ggplot(aes(x=CollectionDate, y= z_score_group,color=CL4, group=SubjectID)) +geom_line() + geom_point()
p.zcytokine_log <- p.zcytokine_log + facet_wrap(.~SubjectID,scales = "free") + geom_hline(yintercept = 1.5)
p.zcytokine_log

p.zcytokine_log <- filter(dat2, variable == "LEPTIN") %>% ggplot(aes(x=CollectionDate, y= value,color=CL4, group=SubjectID)) +geom_line() + geom_point()
p.zcytokine_log <- p.zcytokine_log + facet_wrap(.~SubjectID) + geom_hline(yintercept = 1.5)
p.zcytokine_log

p.cytokine_log <- filter(dat2, variable == "IP10") %>% ggplot(aes(x=CollectionDate, y= z_score_group,color=CL4, group=SubjectID)) +geom_line() + geom_point()
p.cytokine_log <- p.cytokine_log + facet_wrap(.~SubjectID, scale="free") + geom_hline(yintercept = 110)
p.cytokine_log
 
p.zcytokine_log <- ggplot(dat2, aes(x= z_score_ungrouped)) + geom_density()
p.zcytokine_log <- p.zcytokine_log + facet_wrap(.~variable, scales="free") 
p.zcytokine_log

dat2

#display limit of detection in the MFI data
vlinedata <- select(llod, -mean, - sd)
vlinedata$variable <- row.names(vlinedata)
p.zcytokine_llod <- ggplot(dat, aes(x= value)) + geom_density() + theme_cowplot(font_size=14)
p.zcytokine_llod <- p.zcytokine_llod + facet_wrap(.~variable, scales="free") + geom_vline(aes(xintercept=llod.mfi), vlinedata, color="red")
p.zcytokine_llod
ggsave(filename = "./Analysis/Suppl.figure/cytokine.llod.pdf", p.zcytokine_llod, width = 25, height=12, dpi=300)


pbar <- ggplot(dat3, aes(x=variable, y=z_score_group)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
pbar <- pbar +  geom_hline(yintercept = 1.5)
pbar

llod

circulating.ck.list <- c("BNDF", "EOTAXIN", "HGF", "ICAM1", "Il12P40","IL7","LEPTIN", "MCP1", "PAI1", "PDGFBB","RANTES","RESISTIN", "TNFA","VCAM")

count(dat,"SubjectID")

clinic.df$BASO + clinic.df$NEUT + clinic.df$EOS + clinic.df$LYM + clinic.df$MONO

colnames(clinic.df)[str_detect(colnames(clinic.df), "AB")]





