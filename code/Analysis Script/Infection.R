#Infection/Immunization 

#This is to carry out analysis based on time points
library(ggplot2)
library(dplyr)
library(tidyverse)
library(Mfuzz)
library(reshape2)

base_theme =
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 13),
        panel.grid.minor = element_blank())

setwd("~/Box/XinZhouFiles/Projects/ZXP2_HMP2_Core/Analysis/Infection/")

load("../Robject/DetailedPhyloseq.RData")
load("../Robject/Revision_MultiOmes_0509.RData")
load("../Robject/Prevalance.RData")

meta_Subject_HMP <- read.csv("../ori_meta_table/metadata.subject.csv", header=T, row.names = 1)
meta_Sample_HMP <- read.csv("../ori_meta_table/hmp2_infection_metadata_0323.csv", header=T)
Dateinfo <- read.csv("../ori_meta_table/DateInfo.csv", header = T) %>% select(-X)
Dateinfo$Date <- as.Date(Dateinfo$Date, format = "%m/%d/%y")
meta_Sample_HMP <- Dateinfo %>% merge(meta_Sample_HMP, by="SampleID",all = T) %>% select(-SubjectID, -CollectionDate)
meta_Sample_HMP <- filter(meta_Sample_HMP,!is.na(Date))

table(meta_Sample_HMP$CL4)
table(meta_Sample_HMP$InfectionState_jax)
meta_Sample_HMP <- meta_Sample_HMP %>% group_by(Subject) %>% mutate(rank = rank(Date))
meta_Sample_HMP <- meta_Sample_HMP[with(meta_Sample_HMP, order(Subject, rank)),]
meta_Sample_HMP$rownames <- 1:length(meta_Sample_HMP$SampleID)

#meta_Sample_HMP[meta_Sample_HMP$InfectionState_jax=="post",]

#INTP=Infection Time Point
#INTP <- meta_Sample_HMP %>% filter(InfectionState_jax %in% c("pre", "early","late", "recovery","post"))
INTP2 <- meta_Sample_HMP %>% filter(InfectionState_jax %in% c("early","late", "recovery"))

#IMTP=Infection Time Point
IMTP <- meta_Sample_HMP %>% filter(CL4 %in% c("Imz", "Imz_L"))

#Antibiotics Time Point
ATP <- meta_Sample_HMP %>% filter(CL4 %in% c("Ant", "Ant_L"))

INTP1 <- INTP2 %>% group_by(Subject) %>% group_by(seq_id = cumsum(c(1, diff(rownames)) != 1)) %>% nest() 
INTP1
unlist(lapply(lapply(INTP1$data,"[[", c("SampleID")), head, n = 1L))
unlist(lapply(lapply(INTP1$data,"[[", c("SampleID")), tail, n = 1L))

IMTP1 <- IMTP %>% group_by(Subject) %>% group_by(seq_id = cumsum(c(1, diff(rank)) != 1)) %>% nest() 
unlist(lapply(lapply(IMTP1$data,"[[", c("SampleID")), head, n = 1L))
unlist(lapply(lapply(IMTP1$data,"[[", c("SampleID")), tail, n = 1L))

ATP1 <- ATP %>% group_by(Subject) %>% group_by(seq_id = cumsum(c(1, diff(rank)) != 1)) %>% nest() 
unlist(lapply(lapply(ATP1$data,"[[", c("SampleID")), head, n = 1L))
unlist(lapply(lapply(ATP1$data,"[[", c("SampleID")), tail, n = 1L))

meta_Sample_HMP$IN_event <- "other"
meta_Sample_HMP$IM_event <- "other"
meta_Sample_HMP$AT_event <- "other"

for (i in 01:51){
  temp.sampleid <- unlist(lapply(INTP1$data,"[[", c("SampleID"))[i])
  print(temp.sampleid)
  meta_Sample_HMP$IN_event[meta_Sample_HMP$SampleID %in% temp.sampleid] <- formatC(format = "d", i, flag= "0", width=2)
  print(formatC(format = "d", i, flag= "0", width=2))
}

for (i in 01:33){
  temp.sampleid <- unlist(lapply(IMTP1$data,"[[", c("SampleID"))[i])
  print(temp.sampleid)
  meta_Sample_HMP$IM_event[meta_Sample_HMP$SampleID %in% temp.sampleid] <- formatC(format = "d", i, flag= "0", width=2)
  print(formatC(format = "d", i, flag= "0", width=2))
}

for (i in 01:07){
  temp.sampleid <- unlist(lapply(ATP1$data,"[[", c("SampleID"))[i])
  print(temp.sampleid)
  meta_Sample_HMP$AT_event[meta_Sample_HMP$SampleID %in% temp.sampleid] <- formatC(format = "d", i, flag= "0", width=2)
  print(formatC(format = "d", i, flag= "0", width=2))
}

meta_Sample_HMP$INTP_jax <- meta_Sample_HMP$InfectionState_jax
meta_Sample_HMP$IMTP_jax <- meta_Sample_HMP$CL4
meta_Sample_HMP$ATP_jax <- meta_Sample_HMP$CL4

meta_Sample_HMP$INTP_jax[! (meta_Sample_HMP$InfectionState_jax %in% c("early","late", "recovery")) ] <- "other"
meta_Sample_HMP$IMTP_jax[! (meta_Sample_HMP$CL4 %in% c("Imz","Imz_L")) ] <- "other"
meta_Sample_HMP$ATP_jax[! (meta_Sample_HMP$CL4 %in% c("Ant","Ant_L")) ] <- "other"

meta_Sample_HMP$INTP_jax[meta_Sample_HMP$rownames %in% (meta_Sample_HMP$rownames[meta_Sample_HMP$SampleID %in% unlist(lapply(lapply(INTP1$data,"[[", c("SampleID")), head, n = 1L))] - 1)] <- "pre"
meta_Sample_HMP$INTP_jax[meta_Sample_HMP$rownames %in% (meta_Sample_HMP$rownames[meta_Sample_HMP$SampleID %in% unlist(lapply(lapply(INTP1$data,"[[", c("SampleID")), tail, n = 1L))] + 1)] <- "post"
table(meta_Sample_HMP$INTP_jax)

meta_Sample_HMP$IMTP_jax[meta_Sample_HMP$rownames %in% (meta_Sample_HMP$rownames[meta_Sample_HMP$SampleID %in% unlist(lapply(lapply(IMTP1$data,"[[", c("SampleID")), head, n = 1L))] - 1)] <- "pre"
meta_Sample_HMP$IMTP_jax[meta_Sample_HMP$rownames %in% (meta_Sample_HMP$rownames[meta_Sample_HMP$SampleID %in% unlist(lapply(lapply(IMTP1$data,"[[", c("SampleID")), tail, n = 1L))] + 1)] <- "post"
table(meta_Sample_HMP$IMTP_jax)

meta_Sample_HMP$IMTP_jax[meta_Sample_HMP$SampleID == c("70-1011-01","69-027-01")] <- "other"

meta_Sample_HMP$ATP_jax[meta_Sample_HMP$rownames %in% (meta_Sample_HMP$rownames[meta_Sample_HMP$SampleID %in% unlist(lapply(lapply(ATP1$data,"[[", c("SampleID")), head, n = 1L))] - 1)] <- "pre"
meta_Sample_HMP$ATP_jax[meta_Sample_HMP$rownames %in% (meta_Sample_HMP$rownames[meta_Sample_HMP$SampleID %in% unlist(lapply(lapply(ATP1$data,"[[", c("SampleID")), tail, n = 1L))] + 1)] <- "post"
table(meta_Sample_HMP$ATP_jax)

INTP <- meta_Sample_HMP %>% filter(INTP_jax != "other")

#IMTP=Infection Time Point
IMTP <- meta_Sample_HMP %>% filter(IMTP_jax != "other")

#Antibiotics Time Point
ATP <- meta_Sample_HMP %>% filter(ATP_jax != "other")

table(meta_Sample_HMP$INTP_jax)
table(meta_Sample_HMP$IMTP_jax)
table(meta_Sample_HMP$ATP_jax)

ggplot(INTP, aes(x=Date, y=Subject, fill=INTP_jax)) + geom_tile()
ggplot(IMTP, aes(x=Date, y=Subject, fill=IMTP_jax)) + geom_tile()
ggplot(ATP, aes(x=Date, y=Subject, fill=ATP_jax)) + geom_tile()

write.csv(meta_Sample_HMP,"../ori_meta_table/hmp2_infection_metadata_1027.csv")

#prepare sample data
sample_stool_0 <- data.frame(sample_data(physeqGenus_ST))
sample_skin_0 <- data.frame(sample_data(physeqGenus_SK))
sample_oral_0 <- data.frame(sample_data(physeqGenus_OR))
sample_nasal_0 <- data.frame(sample_data(physeqGenus_NS))

sample_stool <- merge(meta_Sample_HMP, sample_stool_0, by="SampleID",all=T)
sample_skin <- merge(meta_Sample_HMP, sample_skin_0, by="SampleID", all=T)
sample_oral <- merge(meta_Sample_HMP, sample_oral_0, by="SampleID", all=T)
sample_nasal <- merge(meta_Sample_HMP, sample_nasal_0, by="SampleID", all=T)

table(sample_stool$InfectionState_jax)

Sample_metadata <- select(sample_stool_0, SampleID, RandomID, Date) %>% rename(Date_stool=Date,ST=RandomID) %>% merge(Dateinfo, by="SampleID", all.y=T)
Sample_metadata <- select(sample_skin_0, SampleID, KitID, Date) %>% rename(Date_skin=Date,SK=KitID)  %>% merge(Sample_metadata, by="SampleID", all.y=T)
Sample_metadata <- select(sample_oral_0, SampleID, KitID, Date) %>% rename(Date_oral=Date,OR=KitID)  %>% merge(Sample_metadata, by="SampleID", all.y=T)
Sample_metadata <- select(sample_nasal_0, SampleID, RandomID, Date) %>% rename(Date_nasal=Date,NS=RandomID)  %>% merge(Sample_metadata, by="SampleID", all.y=T)


head(Dateinfo)

#prepare stool genus table
stool.genus <- data.frame(otu_table(physeqGenus_ST))
stool.genus <- stool.genus[,colSums(stool.genus)!=0]
tax.stool <- data.frame(tax_table(physeqGenus_ST))
tax.stool$ASV <- rownames(tax.stool)
colnames(stool.genus) <- tax.stool$Genus[match(colnames(stool.genus),tax.stool$ASV)]

stool.long.data <- sample_stool %>% select(Subject, RandomID) %>% merge(stool.genus, by.x="RandomID", by.y ="row.names") %>% gather(variable, value, -Subject, -RandomID) %>% group_by(Subject,variable) %>% mutate(zscore = scale(value)) %>% ungroup
stool.long.data$timepoint <- "other"
stool.long.data$timepoint[stool.long.data$RandomID %in% sample_stool$RandomID[sample_stool$InfectionState_jax== "pre"]] <- "1_pre"
stool.long.data$timepoint[stool.long.data$RandomID %in% sample_stool$RandomID[sample_stool$InfectionState_jax== "early"]] <- "2_ely"
stool.long.data$timepoint[stool.long.data$RandomID %in% sample_stool$RandomID[sample_stool$InfectionState_jax== "late"]] <- "3_lat"
stool.long.data$timepoint[stool.long.data$RandomID %in% sample_stool$RandomID[sample_stool$InfectionState_jax== "recovery"]] <- "4_rcv"
stool.long.data$timepoint[stool.long.data$RandomID %in% sample_stool$RandomID[sample_stool$InfectionState_jax== "post"]] <- "5_pst"
stool.long.data <- filter(stool.long.data, timepoint != "other")
stool.long.data$genus <- paste(stool.long.data$variable,"stool", sep="_")
stool.long.data$bodysite <- "stool"
stool.long.data <- filter(stool.long.data, stool.long.data$zscore != "NaN")
stool.long.data <- stool.long.data %>% group_by(genus, timepoint) %>% mutate(abundance = mean(zscore))

stool.anti.long.data <- sample_stool %>% select(Subject, RandomID) %>% merge(stool.genus, by.x="RandomID", by.y ="row.names") %>% gather(variable, value, -Subject, -RandomID) %>% group_by(Subject,variable) %>% mutate(zscore = scale(value)) %>% ungroup
stool.anti.long.data$timepoint <- "other"
stool.anti.long.data$timepoint[stool.anti.long.data$RandomID %in% sample_stool$RandomID[sample_stool$ATP_jax== "pre"]] <- "1_pre"
stool.anti.long.data$timepoint[stool.anti.long.data$RandomID %in% sample_stool$RandomID[sample_stool$ATP_jax== "Ant"]] <- "2_ely"
stool.anti.long.data$timepoint[stool.anti.long.data$RandomID %in% sample_stool$RandomID[sample_stool$ATP_jax== "Ant_L"]] <- "3_lat"
stool.anti.long.data$timepoint[stool.anti.long.data$RandomID %in% sample_stool$RandomID[sample_stool$ATP_jax== "post"]] <- "4_pst"
stool.anti.long.data <- filter(stool.anti.long.data, timepoint != "other")
stool.anti.long.data$genus <- paste(stool.anti.long.data$variable,"stool", sep="_")
stool.anti.long.data$bodysite <- "stool"
stool.anti.long.data <- filter(stool.anti.long.data, stool.anti.long.data$zscore != "NaN")
stool.anti.long.data <- stool.anti.long.data %>% group_by(genus, timepoint) %>% mutate(abundance = mean(zscore))

stool.imz.long.data <- sample_stool %>% select(Subject, RandomID) %>% merge(stool.genus, by.x="RandomID", by.y ="row.names") %>% gather(variable, value, -Subject, -RandomID) %>% group_by(Subject,variable) %>% mutate(zscore = scale(value)) %>% ungroup
stool.imz.long.data$timepoint <- "other"
stool.imz.long.data$timepoint[stool.imz.long.data$RandomID %in% sample_stool$RandomID[sample_stool$IMTP_jax== "pre"]] <- "1_pre"
stool.imz.long.data$timepoint[stool.imz.long.data$RandomID %in% sample_stool$RandomID[sample_stool$IMTP_jax== "Imz"]] <- "2_ely"
stool.imz.long.data$timepoint[stool.imz.long.data$RandomID %in% sample_stool$RandomID[sample_stool$IMTP_jax== "Imz_L"]] <- "3_lat"
stool.imz.long.data$timepoint[stool.imz.long.data$RandomID %in% sample_stool$RandomID[sample_stool$IMTP_jax== "post"]] <- "4_pst"
stool.imz.long.data <- filter(stool.imz.long.data, timepoint != "other")
stool.imz.long.data$genus <- paste(stool.imz.long.data$variable,"stool", sep="_")
stool.imz.long.data$bodysite <- "stool"
stool.imz.long.data <- filter(stool.imz.long.data, stool.imz.long.data$zscore != "NaN")
stool.imz.long.data <- stool.imz.long.data %>% group_by(genus, timepoint) %>% mutate(abundance = mean(zscore))

#Skin data
skin.genus <- data.frame(otu_table(physeqGenus_SK))
skin.genus <- skin.genus[,colSums(skin.genus)!=0]
tax.skin <- data.frame(tax_table(physeqGenus_SK))
tax.skin$ASV <- rownames(tax.skin)
colnames(skin.genus) <- tax.skin$Genus[match(colnames(skin.genus),tax.skin$ASV)]

skin.long.data <- sample_skin %>% select(Subject, KitID) %>% merge(skin.genus, by.x="KitID", by.y ="row.names") %>% gather(variable, value, -Subject, -KitID) %>% group_by(Subject,variable) %>% mutate(zscore = scale(value)) %>% ungroup
skin.long.data$timepoint <- "other"
skin.long.data$timepoint[skin.long.data$KitID %in% sample_skin$KitID[sample_skin$InfectionState_jax== "pre"]] <- "1_pre"
skin.long.data$timepoint[skin.long.data$KitID %in% sample_skin$KitID[sample_skin$InfectionState_jax== "early"]] <- "2_ely"
skin.long.data$timepoint[skin.long.data$KitID %in% sample_skin$KitID[sample_skin$InfectionState_jax== "late"]] <- "3_lat"
skin.long.data$timepoint[skin.long.data$KitID %in% sample_skin$KitID[sample_skin$InfectionState_jax== "recovery"]] <- "4_rcv"
skin.long.data$timepoint[skin.long.data$KitID %in% sample_skin$KitID[sample_skin$InfectionState_jax== "post"]] <- "5_pst"
skin.long.data <- filter(skin.long.data, timepoint != "other")
skin.long.data$genus <- paste(skin.long.data$variable,"skin", sep="_")
skin.long.data$bodysite <- "skin"
skin.long.data <- filter(skin.long.data, skin.long.data$zscore != "NaN")
skin.long.data <- skin.long.data %>% group_by(genus, timepoint) %>% mutate(abundance = mean(zscore))

skin.anti.long.data <- sample_skin %>% select(Subject, KitID) %>% merge(skin.genus, by.x="KitID", by.y ="row.names") %>% gather(variable, value, -Subject, -KitID) %>% group_by(Subject,variable) %>% mutate(zscore = scale(value)) %>% ungroup
skin.anti.long.data$timepoint <- "other"
skin.anti.long.data$timepoint[skin.anti.long.data$KitID %in% sample_skin$KitID[sample_skin$ATP_jax== "pre"]] <- "1_pre"
skin.anti.long.data$timepoint[skin.anti.long.data$KitID %in% sample_skin$KitID[sample_skin$ATP_jax== "Ant"]] <- "2_ely"
skin.anti.long.data$timepoint[skin.anti.long.data$KitID %in% sample_skin$KitID[sample_skin$ATP_jax== "Ant_L"]] <- "3_lat"
skin.anti.long.data$timepoint[skin.anti.long.data$KitID %in% sample_skin$KitID[sample_skin$ATP_jax== "post"]] <- "4_pst"
skin.anti.long.data <- filter(skin.anti.long.data, timepoint != "other")
skin.anti.long.data$genus <- paste(skin.anti.long.data$variable,"skin", sep="_")
skin.anti.long.data$bodysite <- "skin"
skin.anti.long.data <- filter(skin.anti.long.data, skin.anti.long.data$zscore != "NaN")
skin.anti.long.data <- skin.anti.long.data %>% group_by(genus, timepoint) %>% mutate(abundance = mean(zscore))

skin.imz.long.data <- sample_skin %>% select(Subject, KitID) %>% merge(skin.genus, by.x="KitID", by.y ="row.names") %>% gather(variable, value, -Subject, -KitID) %>% group_by(Subject,variable) %>% mutate(zscore = scale(value)) %>% ungroup
skin.imz.long.data$timepoint <- "other"
skin.imz.long.data$timepoint[skin.imz.long.data$KitID %in% sample_skin$KitID[sample_skin$IMTP_jax== "pre"]] <- "1_pre"
skin.imz.long.data$timepoint[skin.imz.long.data$KitID %in% sample_skin$KitID[sample_skin$IMTP_jax== "Imz"]] <- "2_ely"
skin.imz.long.data$timepoint[skin.imz.long.data$KitID %in% sample_skin$KitID[sample_skin$IMTP_jax== "Imz_L"]] <- "3_lat"
skin.imz.long.data$timepoint[skin.imz.long.data$KitID %in% sample_skin$KitID[sample_skin$IMTP_jax== "post"]] <- "4_pst"
skin.imz.long.data <- filter(skin.imz.long.data, timepoint != "other")
skin.imz.long.data$genus <- paste(skin.imz.long.data$variable,"skin", sep="_")
skin.imz.long.data$bodysite <- "skin"
skin.imz.long.data <- filter(skin.imz.long.data, skin.imz.long.data$zscore != "NaN")
skin.imz.long.data <- skin.imz.long.data %>% group_by(genus, timepoint) %>% mutate(abundance = mean(zscore))

#oral data
oral.genus <- data.frame(otu_table(physeqGenus_OR))
oral.genus <- oral.genus[,colSums(oral.genus)!=0]
tax.oral <- data.frame(tax_table(physeqGenus_OR))
tax.oral$ASV <- rownames(tax.oral)
colnames(oral.genus) <- tax.oral$Genus[match(colnames(oral.genus),tax.oral$ASV)]

oral.long.data <- sample_oral %>% select(Subject, KitID) %>% merge(oral.genus, by.x="KitID", by.y ="row.names") %>% gather(variable, value, -Subject, -KitID) %>% group_by(Subject,variable) %>% mutate(zscore = scale(value)) %>% ungroup
oral.long.data$timepoint <- "other"
oral.long.data$timepoint[oral.long.data$KitID %in% sample_oral$KitID[sample_oral$InfectionState_jax== "pre"]] <- "1_pre"
oral.long.data$timepoint[oral.long.data$KitID %in% sample_oral$KitID[sample_oral$InfectionState_jax== "early"]] <- "2_ely"
oral.long.data$timepoint[oral.long.data$KitID %in% sample_oral$KitID[sample_oral$InfectionState_jax== "late"]] <- "3_lat"
oral.long.data$timepoint[oral.long.data$KitID %in% sample_oral$KitID[sample_oral$InfectionState_jax== "recovery"]] <- "4_rcv"
oral.long.data$timepoint[oral.long.data$KitID %in% sample_oral$KitID[sample_oral$InfectionState_jax== "post"]] <- "5_pst"
oral.long.data <- filter(oral.long.data, timepoint != "other")
oral.long.data$genus <- paste(oral.long.data$variable,"oral", sep="_")
oral.long.data$bodysite <- "oral"
oral.long.data <- filter(oral.long.data, oral.long.data$zscore != "NaN")
oral.long.data <- oral.long.data %>% group_by(genus, timepoint) %>% mutate(abundance = mean(zscore))

oral.anti.long.data <- sample_oral %>% select(Subject, KitID) %>% merge(oral.genus, by.x="KitID", by.y ="row.names") %>% gather(variable, value, -Subject, -KitID) %>% group_by(Subject,variable) %>% mutate(zscore = scale(value)) %>% ungroup
oral.anti.long.data$timepoint <- "other"
oral.anti.long.data$timepoint[oral.anti.long.data$KitID %in% sample_oral$KitID[sample_oral$ATP_jax== "pre"]] <- "1_pre"
oral.anti.long.data$timepoint[oral.anti.long.data$KitID %in% sample_oral$KitID[sample_oral$ATP_jax== "Ant"]] <- "2_ely"
oral.anti.long.data$timepoint[oral.anti.long.data$KitID %in% sample_oral$KitID[sample_oral$ATP_jax== "Ant_L"]] <- "3_lat"
oral.anti.long.data$timepoint[oral.anti.long.data$KitID %in% sample_oral$KitID[sample_oral$ATP_jax== "post"]] <- "4_pst"
oral.anti.long.data <- filter(oral.anti.long.data, timepoint != "other")
oral.anti.long.data$genus <- paste(oral.anti.long.data$variable,"oral", sep="_")
oral.anti.long.data$bodysite <- "oral"
oral.anti.long.data <- filter(oral.anti.long.data, oral.anti.long.data$zscore != "NaN")
oral.anti.long.data <- oral.anti.long.data %>% group_by(genus, timepoint) %>% mutate(abundance = mean(zscore))

oral.imz.long.data <- sample_oral %>% select(Subject, KitID) %>% merge(oral.genus, by.x="KitID", by.y ="row.names") %>% gather(variable, value, -Subject, -KitID) %>% group_by(Subject,variable) %>% mutate(zscore = scale(value)) %>% ungroup
oral.imz.long.data$timepoint <- "other"
oral.imz.long.data$timepoint[oral.imz.long.data$KitID %in% sample_oral$KitID[sample_oral$IMTP_jax== "pre"]] <- "1_pre"
oral.imz.long.data$timepoint[oral.imz.long.data$KitID %in% sample_oral$KitID[sample_oral$IMTP_jax== "Imz"]] <- "2_ely"
oral.imz.long.data$timepoint[oral.imz.long.data$KitID %in% sample_oral$KitID[sample_oral$IMTP_jax== "Imz_L"]] <- "3_lat"
oral.imz.long.data$timepoint[oral.imz.long.data$KitID %in% sample_oral$KitID[sample_oral$IMTP_jax== "post"]] <- "4_pst"
oral.imz.long.data <- filter(oral.imz.long.data, timepoint != "other")
oral.imz.long.data$genus <- paste(oral.imz.long.data$variable,"oral", sep="_")
oral.imz.long.data$bodysite <- "oral"
oral.imz.long.data <- filter(oral.imz.long.data, oral.imz.long.data$zscore != "NaN")
oral.imz.long.data <- oral.imz.long.data %>% group_by(genus, timepoint) %>% mutate(abundance = mean(zscore))

#nasal
nasal.genus <- data.frame(otu_table(physeqGenus_NS))
nasal.genus <- nasal.genus[,colSums(nasal.genus)!=0]
tax.nasal <- data.frame(tax_table(physeqGenus_NS))
tax.nasal$ASV <- rownames(tax.nasal)
colnames(nasal.genus) <- tax.nasal$Genus[match(colnames(nasal.genus),tax.nasal$ASV)]

nasal.long.data <- sample_nasal %>% select(Subject, RandomID) %>% merge(nasal.genus, by.x="RandomID", by.y ="row.names") %>% gather(variable, value, -Subject, -RandomID) %>% group_by(Subject,variable) %>% mutate(zscore = scale(value)) %>% ungroup
nasal.long.data$timepoint <- "other"
nasal.long.data$timepoint[nasal.long.data$RandomID %in% sample_nasal$RandomID[sample_nasal$InfectionState_jax== "pre"]] <- "1_pre"
nasal.long.data$timepoint[nasal.long.data$RandomID %in% sample_nasal$RandomID[sample_nasal$InfectionState_jax== "early"]] <- "2_ely"
nasal.long.data$timepoint[nasal.long.data$RandomID %in% sample_nasal$RandomID[sample_nasal$InfectionState_jax== "late"]] <- "3_lat"
nasal.long.data$timepoint[nasal.long.data$RandomID %in% sample_nasal$RandomID[sample_nasal$InfectionState_jax== "recovery"]] <- "4_rcv"
nasal.long.data$timepoint[nasal.long.data$RandomID %in% sample_nasal$RandomID[sample_nasal$InfectionState_jax== "post"]] <- "5_pst"
nasal.long.data <- filter(nasal.long.data, timepoint != "other")
nasal.long.data$genus <- paste(nasal.long.data$variable,"nasal", sep="_")
nasal.long.data$bodysite <- "nasal"
nasal.long.data <- filter(nasal.long.data, nasal.long.data$zscore != "NaN")
nasal.long.data <- nasal.long.data %>% group_by(genus, timepoint) %>% mutate(abundance = mean(zscore))

nasal.anti.long.data <- sample_nasal %>% select(Subject, RandomID) %>% merge(nasal.genus, by.x="RandomID", by.y ="row.names") %>% gather(variable, value, -Subject, -RandomID) %>% group_by(Subject,variable) %>% mutate(zscore = scale(value)) %>% ungroup
nasal.anti.long.data$timepoint <- "other"
nasal.anti.long.data$timepoint[nasal.anti.long.data$RandomID %in% sample_nasal$RandomID[sample_nasal$ATP_jax== "pre"]] <- "1_pre"
nasal.anti.long.data$timepoint[nasal.anti.long.data$RandomID %in% sample_nasal$RandomID[sample_nasal$ATP_jax== "Ant"]] <- "2_ely"
nasal.anti.long.data$timepoint[nasal.anti.long.data$RandomID %in% sample_nasal$RandomID[sample_nasal$ATP_jax== "Ant_L"]] <- "3_lat"
nasal.anti.long.data$timepoint[nasal.anti.long.data$RandomID %in% sample_nasal$RandomID[sample_nasal$ATP_jax== "post"]] <- "4_pst"
nasal.anti.long.data <- filter(nasal.anti.long.data, timepoint != "other")
nasal.anti.long.data$genus <- paste(nasal.anti.long.data$variable,"nasal", sep="_")
nasal.anti.long.data$bodysite <- "nasal"
nasal.anti.long.data <- filter(nasal.anti.long.data, nasal.anti.long.data$zscore != "NaN")
nasal.anti.long.data <- nasal.anti.long.data %>% group_by(genus, timepoint) %>% mutate(abundance = mean(zscore))

nasal.imz.long.data <- sample_nasal %>% select(Subject, RandomID) %>% merge(nasal.genus, by.x="RandomID", by.y ="row.names") %>% gather(variable, value, -Subject, -RandomID) %>% group_by(Subject,variable) %>% mutate(zscore = scale(value)) %>% ungroup
nasal.imz.long.data$timepoint <- "other"
nasal.imz.long.data$timepoint[nasal.imz.long.data$RandomID %in% sample_nasal$RandomID[sample_nasal$IMTP_jax== "pre"]] <- "1_pre"
nasal.imz.long.data$timepoint[nasal.imz.long.data$RandomID %in% sample_nasal$RandomID[sample_nasal$IMTP_jax== "Imz"]] <- "2_ely"
nasal.imz.long.data$timepoint[nasal.imz.long.data$RandomID %in% sample_nasal$RandomID[sample_nasal$IMTP_jax== "Imz_L"]] <- "3_lat"
nasal.imz.long.data$timepoint[nasal.imz.long.data$RandomID %in% sample_nasal$RandomID[sample_nasal$IMTP_jax== "post"]] <- "4_pst"
nasal.imz.long.data <- filter(nasal.imz.long.data, timepoint != "other")
nasal.imz.long.data$genus <- paste(nasal.imz.long.data$variable,"nasal", sep="_")
nasal.imz.long.data$bodysite <- "nasal"
nasal.imz.long.data <- filter(nasal.imz.long.data, nasal.imz.long.data$zscore != "NaN")
nasal.imz.long.data <- nasal.imz.long.data %>% group_by(genus, timepoint) %>% mutate(abundance = mean(zscore))

stool.wide.data <- stool.long.data %>% select(genus, timepoint, abundance) %>% dcast(genus ~ timepoint, value.var="abundance",median)
skin.wide.data <- skin.long.data %>% select(genus, timepoint, abundance) %>% dcast(genus ~ timepoint, value.var="abundance",median)
oral.wide.data <- oral.long.data %>% select(genus, timepoint, abundance) %>% dcast(genus ~ timepoint, value.var="abundance",median)
nasal.wide.data <- nasal.long.data %>% select(genus, timepoint, abundance) %>% dcast(genus ~ timepoint, value.var="abundance",median)

stool.anti.wide.data <- stool.anti.long.data %>% select(genus, timepoint, abundance) %>% dcast(genus ~ timepoint, value.var="abundance",median)
skin.anti.wide.data <- skin.anti.long.data %>% select(genus, timepoint, abundance) %>% dcast(genus ~ timepoint, value.var="abundance",median)
oral.anti.wide.data <- oral.anti.long.data %>% select(genus, timepoint, abundance) %>% dcast(genus ~ timepoint, value.var="abundance",median)
nasal.anti.wide.data <- nasal.anti.long.data %>% select(genus, timepoint, abundance) %>% dcast(genus ~ timepoint, value.var="abundance",median)

stool.imz.wide.data <- stool.imz.long.data %>% select(genus, timepoint, abundance) %>% dcast(genus ~ timepoint, value.var="abundance",median)
skin.imz.wide.data <- skin.imz.long.data %>% select(genus, timepoint, abundance) %>% dcast(genus ~ timepoint, value.var="abundance",median)
oral.imz.wide.data <- oral.imz.long.data %>% select(genus, timepoint, abundance) %>% dcast(genus ~ timepoint, value.var="abundance",median)
nasal.imz.wide.data <- nasal.imz.long.data %>% select(genus, timepoint, abundance) %>% dcast(genus ~ timepoint, value.var="abundance",median)

all.wide.data <- rbind(stool.wide.data,skin.wide.data,oral.wide.data,nasal.wide.data)
all.wide.data.anti <- rbind(stool.anti.wide.data, skin.anti.wide.data, oral.anti.wide.data, nasal.anti.wide.data)
all.wide.data.imz <- rbind(stool.imz.wide.data,skin.imz.wide.data,oral.imz.wide.data,nasal.imz.wide.data)

wide.data.m <- select(all.wide.data, -genus) %>% data.matrix()
rownames(wide.data.m) <- all.wide.data$genus

wide.data.anti.m <- select(all.wide.data.anti, -genus) %>% data.matrix()
rownames(wide.data.anti.m) <- all.wide.data.anti$genus

wide.data.imz.m <- select(all.wide.data.imz, -genus) %>% data.matrix()
rownames(wide.data.imz.m) <- all.wide.data.imz$genus


#############################
#MFUZZ model
############################
expression.data <- new("ExpressionSet", exprs=wide.data.m)

Infection.fuzz <- filter.NA(expression.data, thres=0.50)
Infection.fuzz <- fill.NA(Infection.fuzz,mode="mean")
Infection.fuzz <- filter.std(Infection.fuzz,min.std=0)
Infection.fuzz <- standardise(Infection.fuzz)

rownames(Infection.fuzz)

#c=number of clusters
m <- mestimate(Infection.fuzz)
Dmin(Infection.fuzz, m, crange=seq(2,10,1), repeats = 3,visu = T)

#Infection:  0.9376453 2.0304539 2.5071461 2.6505353 1.4824889 1.0796205 0.8667652 0.5954440 0.4999062
#Antibiotics: 1.3809674 2.0518644 2.3356519 1.3251860 1.3304577 1.1510846 1.0815840 0.7426700 0.5194115
#Immunization:  1.6084656 2.2588112 2.4557611 1.2065085 1.0263797 0.8512610 0.6694625 0.6101928 0.5798568

Infection.cl <- mfuzz(Infection.fuzz, 4, m)

pdf(file = "./Infection.pdf")
mfuzz.plot2(Infection.fuzz,cl=Infection.cl,min.mem=0.7,mfrow=c(3,3),time.labels=c("pre","early","late","rec","post"),x11 = FALSE)
dev.off()

infection.cl <- Infection.cl
#antibiotics.cl <- Infection.cl
#immuni.cl <- Infection.cl

save(infection.cl,antibiotics.cl,immuni.cl, file = "./Fuzz.result.RData")

Infection.cl$membership
clus_1 <- names(Infection.cl$cluster[Infection.cl$cluster == 1])
clus_2 <- names(Infection.cl$cluster[Infection.cl$cluster == 2])
clus_3 <- names(Infection.cl$cluster[Infection.cl$cluster == 3])
clus_4 <- names(Infection.cl$cluster[Infection.cl$cluster == 4])
clus_5 <- names(Infection.cl$cluster[Infection.cl$cluster == 5])
clus_6 <- names(Infection.cl$cluster[Infection.cl$cluster == 6])

table(str_detect(all.wide.data.anti$genus, "nasal")) 

infection <- c(293, 803, 277, 884)
antibiotics <- c(230, 483, 182, 648)
immunization <- c(288, 640, 251,845)

factor <- infection

stool_1 <- table(str_detect(clus_1, "stool"))["TRUE"]/factor[1]
skin_1 <- table(str_detect(clus_1, "skin"))["TRUE"]/factor[2]
oral_1 <- table(str_detect(clus_1, "oral"))["TRUE"]/factor[3]
nasal_1 <- table(str_detect(clus_1, "nasal"))["TRUE"]/factor[4]
clust1 <- c(stool_1,skin_1,oral_1,nasal_1)
names(clust1) <- c("Stool", "Skin", "Oral", "Nasal")
plot_clust <- melt(clust1)
plot_clust$name <- rownames(plot_clust)
p <- ggplot(plot_clust, aes(x=name,y=value)) + geom_bar(stat = "identity", width=0.5)
p

stool_2 <- table(str_detect(clus_2, "stool"))["TRUE"]/factor[1]
skin_2 <- table(str_detect(clus_2, "skin"))["TRUE"]/factor[2]
oral_2 <- table(str_detect(clus_2, "oral"))["TRUE"]/factor[3]
nasal_2 <- table(str_detect(clus_2, "nasal"))["TRUE"]/factor[4]
clust2 <- c(stool_2,skin_2,oral_2,nasal_2)
names(clust2) <- c("Stool", "Skin", "Oral", "Nasal")
plot_clust <- melt(clust2)
plot_clust$name <- rownames(plot_clust)
p <- ggplot(plot_clust, aes(x=name,y=value)) + geom_bar(stat = "identity", width=0.5)
p

stool_3 <- table(str_detect(clus_3, "stool"))["TRUE"]/factor[1]
skin_3 <- table(str_detect(clus_3, "skin"))["TRUE"]/factor[2]
oral_3 <- table(str_detect(clus_3, "oral"))["TRUE"]/factor[3]
nasal_3 <- table(str_detect(clus_3, "nasal"))["TRUE"]/factor[4]
clust3 <- c(stool_3,skin_3,oral_3,nasal_3)
names(clust3) <- c("Stool", "Skin", "Oral", "Nasal")
plot_clust <- melt(clust3)
plot_clust$name <- rownames(plot_clust)
p <- ggplot(plot_clust, aes(x=name,y=value)) + geom_bar(stat = "identity", width=0.5)
p

stool_4 <- table(str_detect(clus_4, "stool"))["TRUE"]/factor[1]
skin_4 <- table(str_detect(clus_4, "skin"))["TRUE"]/factor[2]
oral_4 <- table(str_detect(clus_4, "oral"))["TRUE"]/factor[3]
nasal_4 <- table(str_detect(clus_4, "nasal"))["TRUE"]/factor[4]
clust4 <- c(stool_4,skin_4,oral_1,nasal_4)
names(clust4) <- c("Stool", "Skin", "Oral", "Nasal")
plot_clust <- melt(clust4)
plot_clust$name <- rownames(plot_clust)
p <- ggplot(plot_clust, aes(x=name,y=value)) + geom_bar(stat = "identity", width=0.5)
p

stool_5 <- table(str_detect(clus_5, "stool"))["TRUE"]/factor[1]
skin_5 <- table(str_detect(clus_5, "skin"))["TRUE"]/factor[2]
oral_5 <- table(str_detect(clus_5, "oral"))["TRUE"]/factor[3]
nasal_5 <- table(str_detect(clus_5, "nasal"))["TRUE"]/factor[4]
clust5 <- c(stool_5,skin_5,oral_5,nasal_5)
names(clust5) <- c("Stool", "Skin", "Oral", "Nasal")
plot_clust <- melt(clust5)
plot_clust$name <- rownames(plot_clust)
p <- ggplot(plot_clust, aes(x=name,y=value)) + geom_bar(stat = "identity", width=0.5)
p

stool_6 <- table(str_detect(clus_6, "stool"))["TRUE"]/factor[1]
skin_6 <- table(str_detect(clus_6, "skin"))["TRUE"]/factor[2]
oral_6 <- table(str_detect(clus_6, "oral"))["TRUE"]/factor[3]
nasal_6 <- table(str_detect(clus_6, "nasal"))["TRUE"]/factor[4]
clust6 <- c(stool_6,skin_6,oral_6,nasal_6)
names(clust6) <- c("Stool", "Skin", "Oral", "Nasal")
plot_clust <- melt(clust6)
plot_clust$name <- rownames(plot_clust)
p <- ggplot(plot_clust, aes(x=name,y=value)) + geom_bar(stat = "identity", width=0.5)
p

Infection.cl$membership[rownames(Infection.cl$membership)=="Parasutterella_stool",]

infection.mem <- as.data.frame(Infection.cl$membership)
data.frame(rowMax(infection.mem))
infection.mem$max <- rowMax(infection.mem)

#Diversity Analysis
load("../Robject/Diversity_Datatable.RData")
identical(row.names(stool.nasal.ASV.diversity),row.names(stool.nasal.Genus.diversity))
colnames(stool.nasal.Genus.diversity) = paste(colnames(stool.nasal.Genus.diversity),"G", sep = "_")
colnames(skin.oral.Genus.diversity) = paste(colnames(skin.oral.Genus.diversity),"G", sep = "_")

SN_global_diversity <- cbind(stool.nasal.ASV.diversity,stool.nasal.Genus.diversity)
SO_global_diversity <- cbind(skin.oral.ASV.diversity,skin.oral.Genus.diversity)

colnames(SN_global_diversity)
colnames(SO_global_diversity)

ggplot(SN_global_diversity, aes(x=Observed,y=Observed_G, color=SampleType)) + geom_point()

#stool.infection
stool.diversity <- SN_global_diversity %>% filter(SampleType=="ST") %>% merge(sample_stool, by="RandomID") %>% filter(INTP_jax %in% c("pre", "early","late", "recovery", "post"))
table(stool.diversity$INTP_jax)
stool.diversity$INTP_jax <- factor(stool.diversity$INTP_jax, levels=c("pre", "early","late", "recovery", "post"))#"healthy_other",

stool.diversity <- stool.diversity %>% group_by(IN_event) %>% mutate(Observed_Z = scale(Observed))
colnames(stool.diversity)
p.stool.di <- ggplot(stool.diversity, aes(x=INTP_jax,y=Shannon)) + geom_point() + geom_boxplot() + stat_compare_means(label = "p.signif",ref.group ="pre", method = "t.test", method.args = list(alternative = "less"))
p.stool.di <- p.stool.di + base_theme + theme(axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Shannon Diverisity Change During Infection")
p.stool.di

p.stool.di.2 <- ggplot(stool.diversity, aes(x=INTP_jax,y=pielou)) + geom_point() + geom_boxplot() + stat_compare_means(label = "p.signif",ref.group ="pre", method = "t.test", method.args = list(alternative = "less"))
p.stool.di.2 <- p.stool.di.2 + base_theme + theme(axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Pielou Evenness Change During Infection")
p.stool.di.2

#skin.infection
skin.diversity <- SO_global_diversity %>% filter(SampleType=="Skin") %>% merge(sample_skin, by="KitID") %>% filter(INTP_jax %in% c("pre", "early","late", "recovery", "post"))
table(skin.diversity$INTP_jax)
skin.diversity$INTP_jax <- factor(skin.diversity$INTP_jax, levels=c("pre", "early","late", "recovery", "post"))#"healthy_other",

skin.diversity <- skin.diversity %>% group_by(IN_event) %>% mutate(Observed_Z = scale(Observed))

p.skin.di <- ggplot(skin.diversity, aes(x=INTP_jax,y=Shannon)) + geom_point() + geom_boxplot() + stat_compare_means(label = "p.signif",ref.group ="pre", method = "t.test", method.args = list(alternative = "greater"))
p.skin.di <- p.skin.di + base_theme + theme(axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Skin Shannon Diverisity Change During Infection")
p.skin.di

p.skin.di.2 <- ggplot(skin.diversity, aes(x=INTP_jax,y=pielou)) + geom_point() + geom_boxplot() + stat_compare_means(label = "p.signif",ref.group ="pre", method = "t.test", method.args = list(alternative = "greater"))
p.skin.di.2 <- p.skin.di.2 + base_theme + theme(axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Skin Pielou Evenness Change During Infection")
p.skin.di.2

#oral.infection
oral.diversity <- SO_global_diversity %>% filter(SampleType=="Oral") %>% merge(sample_oral, by="KitID") %>% filter(INTP_jax %in% c("pre", "early","late", "recovery", "post"))
table(oral.diversity$INTP_jax)
oral.diversity$INTP_jax <- factor(oral.diversity$INTP_jax, levels=c("pre", "early","late", "recovery", "post"))#"healthy_other",

oral.diversity <- oral.diversity %>% group_by(IN_event) %>% mutate(Observed_Z = scale(Observed))

p.oral.di <- ggplot(oral.diversity, aes(x=INTP_jax,y=Shannon)) + geom_point() + geom_boxplot() + stat_compare_means(label = "p.signif",ref.group ="pre", method = "t.test", method.args = list(alternative = "greater"))
p.oral.di <- p.oral.di + base_theme + theme(axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Oral Shannon Diverisity Change During Infection")
p.oral.di

p.oral.di.2 <- ggplot(oral.diversity, aes(x=INTP_jax,y=pielou)) + geom_point() + geom_boxplot() + stat_compare_means(label = "p.signif",ref.group ="pre", method = "t.test", method.args = list(alternative = "greater"))
p.oral.di.2 <- p.oral.di.2 + base_theme + theme(axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Oral Pielou Evenness Change During Infection")
p.oral.di.2

#nasal.infection
nasal.diversity <- SN_global_diversity %>% filter(SampleType=="NS") %>% merge(sample_nasal, by="RandomID") %>% filter(INTP_jax %in% c("pre", "early","late", "recovery", "post", "other"))
table(nasal.diversity$INTP_jax)
nasal.diversity$INTP_jax <- factor(nasal.diversity$INTP_jax, levels=c("other","pre", "early","late", "recovery", "post")) #"healthy_other",

p.nasal.di <- ggplot(nasal.diversity, aes(x=INTP_jax,y=Shannon)) + geom_point() + geom_boxplot() + stat_compare_means(label = "p.format", ref.group ="other", method.args = list(alternative = "greater"))
p.nasal.di <- p.nasal.di + base_theme + theme(axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Nasal Shannon Diverisity Change During Infection")
p.nasal.di

p.nasal.di.2 <- ggplot(nasal.diversity, aes(x=INTP_jax,y=pielou)) + geom_point() + geom_boxplot() + stat_compare_means(label = "p.format", ref.group ="other", method.args = list(alternative = "greater"))
p.nasal.di.2 <- p.nasal.di.2 + base_theme + theme(axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Nasal Pielou Evenness Change During Infection")
p.nasal.di.2

#stool.immunization
stool.diversity.IM <- SN_global_diversity %>% filter(SampleType=="ST") %>% merge(sample_stool, by="RandomID") %>% filter(IMTP_jax %in% c("pre", "Imz","Imz_L", "post"))
table(stool.diversity.IM$IMTP_jax)
stool.diversity.IM$IMTP_jax <- factor(stool.diversity.IM$IMTP_jax, levels=c("pre", "Imz","Imz_L", "post"))#"healthy_other",

colnames(stool.diversity.IM)
p.stool.di.im <- ggplot(stool.diversity.IM, aes(x=IMTP_jax,y=Shannon)) + geom_point() + geom_boxplot() + stat_compare_means(label = "p.signif",ref.group ="pre", method = "t.test", method.args = list(alternative = "less"))
p.stool.di.im <- p.stool.di.im + base_theme + theme(axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Stool Shannon Diverisity Change During Immunization")
p.stool.di.im

p.stool.di.im.2 <- ggplot(stool.diversity.IM, aes(x=IMTP_jax,y=pielou)) + geom_point() + geom_boxplot() + stat_compare_means(label = "p.signif",ref.group ="pre", method = "t.test", method.args = list(alternative = "less"))
p.stool.di.im.2 <- p.stool.di.im.2 + base_theme + theme(axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Stool Pielou Evenness Change During Immunization")
p.stool.di.im.2

p.stool.di.im + p.stool.di.im.2

#skin.immunization
skin.diversity.IM <- SO_global_diversity %>% filter(SampleType=="Skin") %>% merge(sample_skin, by="KitID") %>% filter(IMTP_jax %in% c("pre", "Imz","Imz_L", "post"))
table(skin.diversity.IM$IMTP_jax)
skin.diversity.IM$IMTP_jax <- factor(skin.diversity.IM$IMTP_jax, levels=c("pre", "Imz","Imz_L", "post"))#"healthy_other",

colnames(skin.diversity.IM)
p.skin.di.im <- ggplot(skin.diversity.IM, aes(x=IMTP_jax,y=Shannon)) + geom_point() + geom_boxplot() + stat_compare_means(label = "p.signif",ref.group ="pre", method = "t.test", method.args = list(alternative = "less"))
p.skin.di.im <- p.skin.di.im + base_theme + theme(axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Skin Shannon Diverisity Change During Immunization")
p.skin.di.im

p.skin.di.im.2 <- ggplot(stool.diversity.IM, aes(x=IMTP_jax,y=pielou)) + geom_point() + geom_boxplot() + stat_compare_means(label = "p.signif",ref.group ="pre", method = "t.test", method.args = list(alternative = "less"))
p.skin.di.im.2 <- p.skin.di.im.2 + base_theme + theme(axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Skin Pielou Evenness Change During Immunization")
p.skin.di.im.2

p.skin.di.im + p.skin.di.im.2


#oral.immunization
oral.diversity.IM <- SO_global_diversity %>% filter(SampleType=="Oral") %>% merge(sample_oral, by="KitID") %>% filter(IMTP_jax %in% c("pre", "Imz","Imz_L", "post"))
table(oral.diversity.IM$IMTP_jax)
oral.diversity.IM$IMTP_jax <- factor(oral.diversity.IM$IMTP_jax, levels=c("pre", "Imz","Imz_L", "post"))#"healthy_other",

colnames(oral.diversity.IM)
p.oral.di.im <- ggplot(oral.diversity.IM, aes(x=IMTP_jax,y=Shannon)) + geom_point() + geom_boxplot() + stat_compare_means(label = "p.signif",ref.group ="pre", method = "t.test", method.args = list(alternative = "less"))
p.oral.di.im <- p.oral.di.im + base_theme + theme(axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Oral Shannon Diverisity Change During Immunization")
p.oral.di.im

p.oral.di.im.2 <- ggplot(stool.diversity.IM, aes(x=IMTP_jax,y=pielou)) + geom_point() + geom_boxplot() + stat_compare_means(label = "p.signif",ref.group ="pre", method = "t.test", method.args = list(alternative = "less"))
p.oral.di.im.2 <- p.oral.di.im.2 + base_theme + theme(axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Oral Pielou Evenness Change During Immunization")
p.oral.di.im.2

p.oral.di.im + p.oral.di.im.2

#nasal.immunization
nasal.diversity.IM <- SN_global_diversity %>% filter(SampleType=="NS") %>% merge(sample_nasal, by="RandomID") %>% filter(IMTP_jax %in% c("pre", "Imz","Imz_L", "post"))
table(nasal.diversity.IM$IMTP_jax)
nasal.diversity.IM$IMTP_jax <- factor(nasal.diversity.IM$IMTP_jax, levels=c("pre", "Imz","Imz_L", "post"))#"healthy_other",

colnames(nasal.diversity.IM)
p.nasal.di.im <- ggplot(nasal.diversity.IM, aes(x=IMTP_jax,y=Shannon)) + geom_point() + geom_boxplot() + stat_compare_means(label = "p.signif",ref.group ="pre", method = "t.test", method.args = list(alternative = "less"))
p.nasal.di.im <- p.nasal.di.im + base_theme + theme(axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Nasal Shannon Diverisity Change During Immunization")
p.nasal.di.im

p.nasal.di.im.2 <- ggplot(nasal.diversity.IM, aes(x=IMTP_jax,y=pielou)) + geom_point() + geom_boxplot() + stat_compare_means(label = "p.signif",ref.group ="pre", method = "t.test", method.args = list(alternative = "less"))
p.nasal.di.im.2 <- p.nasal.di.im.2 + base_theme + theme(axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Nasal Pielou Evenness Change During Immunization")
p.nasal.di.im.2

p.nasal.di.im + p.nasal.di.im.2


#stool.antibiotics
stool.diversity.AT <- SN_global_diversity %>% filter(SampleType=="ST") %>% merge(sample_stool, by="RandomID") %>% filter(ATP_jax %in% c("pre", "Ant","Ant_L", "post"))
table(stool.diversity.AT$ATP_jax)
stool.diversity.AT$ATP_jax <- factor(stool.diversity.AT$ATP_jax, levels=c("pre", "Ant","Ant_L", "post"))#"healthy_other",

colnames(stool.diversity.IM)
p.stool.di.at <- ggplot(stool.diversity.AT, aes(x=ATP_jax,y=Shannon)) + geom_point() + geom_boxplot() + stat_compare_means(label = "p.format",ref.group ="pre", method.args = list(alternative = "less"))
p.stool.di.at <- p.stool.di.at + base_theme + theme(axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Stool Shannon Diverisity Change During Antibiotics")
p.stool.di.at

p.stool.di.at.2 <- ggplot(stool.diversity.AT, aes(x=ATP_jax,y=pielou)) + geom_point() + geom_boxplot() + stat_compare_means(label = "p.signif",ref.group ="pre", method = "t.test", method.args = list(alternative = "less"))
p.stool.di.at.2 <- p.stool.di.at.2 + base_theme + theme(axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Stool Pielou Evenness Change During Antibiotics")
p.stool.di.at.2

p.stool.di.at + p.stool.di.at

#skin.antibiotics
skin.diversity.AT <- SO_global_diversity %>% filter(SampleType=="Skin") %>% merge(sample_skin, by="KitID") %>% filter(ATP_jax %in% c("pre", "Ant","Ant_L", "post"))
table(skin.diversity.AT$ATP_jax)
skin.diversity.AT$ATP_jax <- factor(skin.diversity.AT$ATP_jax, levels=c("pre", "Ant","Ant_L", "post"))#"healthy_other",

colnames(skin.diversity.AT)
p.skin.di.at <- ggplot(skin.diversity.AT, aes(x=ATP_jax,y=Shannon)) + geom_point() + geom_boxplot() + stat_compare_means(label = "p.format",ref.group ="pre",method.args = list(alternative = "greater"))
p.skin.di.at <- p.skin.di.at + base_theme + theme(axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Skin Shannon Diverisity Change During Antibiotics")
p.skin.di.at

p.skin.di.at.2 <- ggplot(skin.diversity.AT, aes(x=ATP_jax,y=pielou)) + geom_point() + geom_boxplot() + stat_compare_means(label = "p.signif",ref.group ="pre", method.args = list(alternative = "greater"))
p.skin.di.at.2 <- p.skin.di.at.2 + base_theme + theme(axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Skin Pielou Evenness Change During Antibiotics")
p.skin.di.at.2

p.skin.di.at + p.skin.di.at.2


#Oral.antibiotics
oral.diversity.AT <- SO_global_diversity %>% filter(SampleType=="Oral") %>% merge(sample_oral, by="KitID") %>% filter(ATP_jax %in% c("pre", "Ant","Ant_L", "post"))
table(oral.diversity.AT$ATP_jax)
oral.diversity.AT$ATP_jax <- factor(oral.diversity.AT$ATP_jax, levels=c("pre", "Ant","Ant_L", "post"))#"healthy_other",

colnames(oral.diversity.AT)
p.oral.di.at <- ggplot(oral.diversity.AT, aes(x=ATP_jax,y=Shannon)) + geom_point() + geom_boxplot() + stat_compare_means(label = "p.format",ref.group ="pre",method.args = list(alternative = "greater"))
p.oral.di.at <- p.oral.di.at + base_theme + theme(axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Oral Shannon Diverisity Change During Antibiotics")
p.oral.di.at

p.oral.di.at.2 <- ggplot(oral.diversity.AT, aes(x=ATP_jax,y=pielou)) + geom_point() + geom_boxplot() + stat_compare_means(label = "p.format",ref.group ="pre", method.args = list(alternative = "greater"))
p.oral.di.at.2 <- p.oral.di.at.2 + base_theme + theme(axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("oral Pielou Evenness Change During Antibiotics")
p.oral.di.at.2

p.oral.di.at + p.oral.di.at.2

#nasal.antibiotics
nasal.diversity.AT <- SN_global_diversity %>% filter(SampleType=="NS") %>% merge(sample_nasal, by="RandomID") %>% filter(ATP_jax %in% c("pre", "Ant","Ant_L", "post"))
table(nasal.diversity.AT$ATP_jax)
nasal.diversity.AT$ATP_jax <- factor(nasal.diversity.AT$ATP_jax, levels=c("pre", "Ant","Ant_L", "post"))#"healthy_other",

colnames(nasal.diversity.AT)
p.nasal.di.at <- ggplot(nasal.diversity.AT, aes(x=ATP_jax,y=Shannon)) + geom_point() + geom_boxplot() + stat_compare_means(label = "p.format",ref.group ="pre", method.args = list(alternative = "less"))
p.nasal.di.at <- p.nasal.di.at + base_theme + theme(axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("nasal Shannon Diverisity Change During Antibiotics")
p.nasal.di.at

p.nasal.di.at.2 <- ggplot(stool.diversity.AT, aes(x=ATP_jax,y=pielou)) + geom_point() + geom_boxplot() + stat_compare_means(label = "p.signif",ref.group ="pre", method = "t.test", method.args = list(alternative = "less"))
p.nasal.di.at.2 <- p.nasal.di.at.2 + base_theme + theme(axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("nasal Pielou Evenness Change During Antibiotics")
p.nasal.di.at.2

p.nasal.di.at + p.nasal.di.at.2

save(stool.diversity, stool.diversity.IM, stool.diversity.AT, file = "./Diversity_By_timepoint.RData")


