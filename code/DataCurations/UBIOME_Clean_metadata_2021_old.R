#Ubiome metadata cleaning
#Xin Zhou
#Aug. 07, 2019
#Updated 03-17-2021
library(tidyverse)
library(taxize)

MIX <- read.csv("~/Box/XinZhouFiles/Projects/ZXP2_HMP2_Core/History/UBIOME Data/06uBiomeMeta/XZ_2016_Mix_WZ.csv", header = T)
Skin2017 <- read.csv("~/Box/XinZhouFiles/Projects/ZXP2_HMP2_Core/History/UBIOME Data/06uBiomeMeta/XZ_2017_skin_WZ.csv", header = T)
Skin2018 <- read.csv("~/Box/XinZhouFiles/Projects/ZXP2_HMP2_Core/History/UBIOME Data/06uBiomeMeta/XZ_2018_Skin_WZ.csv", header = T)
Tongue2017  <- read.csv("~/Box/XinZhouFiles/Projects/ZXP2_HMP2_Core/History/UBIOME Data/06uBiomeMeta/XZ_2017_tongue_WZ.csv", header = T)
Tongue2018  <- read.csv("~/Box/XinZhouFiles/Projects/ZXP2_HMP2_Core/History/UBIOME Data/06uBiomeMeta/XZ_2018_Tongue_WZ.csv", header = T)

Metafile_01 <- rbind(MIX, Skin2017, Skin2018, Tongue2017, Tongue2018)
Metafile_02 <- read.csv("~/Box/XinZhouFiles/Projects/ZXP2_HMP2_Core/History/UBIOME Data/06uBiomeMeta/XZ_KitID_201807.csv", header = T)

length(unique(Metafile_01$KitID)) 

length(intersect(Metafile_02$Kit.id, Metafile_01$KitID)) 

setdiff(Metafile_01$KitID, Metafile_02$Kit.id)

#compare UBM_ID with exsiting sample
ASV.table<- read.csv("~/Desktop/Ubiome result/Tax Table/Ubiome_genus_count.csv", header = T, row.names = 1)


#root.table$KitID <- paste("UBMkid",root.table$barcode, sep="")
#root.table$KitID

#collect those sample that has fastq file, but not sample metadata
missingmeta <- setdiff(row.names(genus.table),Metafile_01$KitID)

finalmetafile <- Metafile_01[which(Metafile_01$KitID %in% row.names(genus.table)),]

#https://cran.r-project.org/web/packages/taxize/vignettes/taxize_vignette.html

tax.list <- colnames(genus.table)
q <- data.frame()
for (i in 1:993){
print(tax.list[i])
q[i,1:6] <- tax_name(tax.list[i], get = c("phylum", "class", "order", "family"), db="itis")
print(tax.list[i])
  }
q

#write.csv(q, "~/Desktop/UbiomTAX.TEMP.CSV")
OTUTABLE <- read.csv("~/Desktop/Ubiome result/Tax Table/Ubiome_genus_count.csv", header=T, row.names = 1)
metafile <- Metafile_01[which(Metafile_01$KitID %in% rownames(OTUTABLE)),]
metafile <- unique(metafile)
metafile2 <- merge(metafile, Metafile_02, by.x="KitID", by.y= "Kit.id")
SAMPLE <- unique(metafile2)
row.names(SAMPLE) <- SAMPLE$KitID
batchinglist <- read.csv("~/Desktop/QC result/batchinglist.csv", header = T)

setdiff(SAMPLE$KitID,batchinglist$KitID)
SAMPLE_V02 <- merge(SAMPLE, batchinglist, by = "KitID", all.y = F) 
#write.csv(SAMPLE, "~/Desktop/UbiomSAMPLE_V01.CSV")
#write.csv(SAMPLE_V02, "~/Desktop/UbiomSAMPLE_V02.csv")







