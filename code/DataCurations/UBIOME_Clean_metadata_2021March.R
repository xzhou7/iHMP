#Ubiome metadata cleaning
#Xin Zhou
#Aug. 07, 2019
#Updated 03-18-2021
library(tidyverse)
library(taxize)

setwd("~/Box/XinZhouFiles/Projects/ZXP2_HMP2_Core/Analysis/")
########################################################################################################################################################
MIX <- read.csv("~/Box/XinZhouFiles/Projects/ZXP2_HMP2_Core/History/UBIOME Data/06uBiomeMeta/XZ_2016_Mix_WZ.csv", header = T)
Skin2017 <- read.csv("~/Box/XinZhouFiles/Projects/ZXP2_HMP2_Core/History/UBIOME Data/06uBiomeMeta/XZ_2017_skin_WZ.csv", header = T)
Skin2018 <- read.csv("~/Box/XinZhouFiles/Projects/ZXP2_HMP2_Core/History/UBIOME Data/06uBiomeMeta/XZ_2018_Skin_WZ.csv", header = T)
Tongue2017  <- read.csv("~/Box/XinZhouFiles/Projects/ZXP2_HMP2_Core/History/UBIOME Data/06uBiomeMeta/XZ_2017_tongue_WZ.csv", header = T)
Tongue2018  <- read.csv("~/Box/XinZhouFiles/Projects/ZXP2_HMP2_Core/History/UBIOME Data/06uBiomeMeta/XZ_2018_Tongue_WZ.csv", header = T)

Metafile_01 <- rbind(MIX, Skin2017, Skin2018, Tongue2017, Tongue2018)
Metafile_02 <- read.csv("~/Box/XinZhouFiles/Projects/ZXP2_HMP2_Core/History/UBIOME Data/06uBiomeMeta/XZ_KitID_201807.csv", header = T)
batchinglist <- read.csv("~/Desktop/QC result/batchinglist.csv", header = T)
########################################################################################################################################################

length(unique(Metafile_01$KitID)) 

length(intersect(Metafile_02$Kit.id, Metafile_01$KitID)) 

setdiff(Metafile_01$KitID, Metafile_02$Kit.id)

#compare UBM_ID with exsiting sample
load("/Users/xzhou7/Box/XinZhouFiles/Projects/ZXP2_HMP2_Core/Analysis/Robject/OTUandTAX.RData")

UBM_clean[1:10,1:10]

#root.table$KitID <- paste("UBMkid",root.table$barcode, sep="")
#root.table$KitID

#collect those sample that has fastq file, but not sample metadata
missingmeta <- setdiff(UBM_clean$X,Metafile_01$KitID)

Metafile_01[duplicated(Metafile_01$KitID),]

#see a lot duplicated in this file, so call the unquie value of this matrix
finalmetafile <- unique(Metafile_01[which(Metafile_01$KitID %in% UBM_clean$X),])
finalmetafile_bodysite <- merge(finalmetafile, Metafile_02, by.x = "KitID", by.y= "Kit.id", all.x = T)

#two sample type unknown, need to put them into PCA
subset(finalmetafile_bodysite, KitID == "UBMkid502104873")
subset(finalmetafile_bodysite, KitID == "UBMkid535104878")
finalmetafile_bodysite$SampleType[finalmetafile_bodysite$KitID %in% c("UBMkid535104878", "UBMkid502104873")] <- "UNK"
       
setdiff(finalmetafile$KitID,finalmetafile_bodysite$KitID)

write.csv(finalmetafile_bodysite, "./UbiomSAMPLEmeta.CSV")








