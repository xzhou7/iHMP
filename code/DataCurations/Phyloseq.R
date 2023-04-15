
library(phyloseq)
library(stringr)
library(ggplot2)
library(cowplot)
library(ShortRead)
library(patchwork)
library(tidyverse)

setwd("~/Box/XinZhouFiles/Projects/ZXP2_HMP2_Core/Analysis/")
load("./Robject/Meta_MultiOmes0413.RData")
load("./Robject/OTUandTAX.RData")
DateInfo <- read.csv("./ori_meta_table/DateInfo.csv", header = T)
load("./Robject/Revision_MultiOmes_0509.RData")

DateInfo_clean <- DateInfo[!is.na(DateInfo$Date),]
DateInfo_clean <- DateInfo_clean[!str_detect(DateInfo_clean$SampleID, "s"),]
DateInfo_missing <-  DateInfo[is.na(DateInfo$Date),]
DateInfo_clean$Date <- as.Date(DateInfo_clean$Date, format = "%m/%d/%y")
DateInfo_clean[DateInfo_clean$Subject=="70-1014",] 

#69-028 is 70-1014
HMP_meta_final$SubjectID[str_detect(HMP_meta_final$V7, "ZMWEIX1")] <- "70-1014"
HMP_meta_final[str_detect(HMP_meta_final$V7, "ZMWEIX1"),]
HMP_meta_final$SampleID[str_detect(HMP_meta_final$V7, "ZMWEIX1")] <- paste("69-028", HMP_meta_final$V2.1[str_detect(HMP_meta_final$V7, "ZMWEIX1")],sep="-")
HMP_meta_final[str_detect(HMP_meta_final$V7, "ZMWEIX1"),]

HMP_meta_final[str_detect(HMP_meta_final$SubjectID, "69-023"),]
HMP_meta_final[is.na(HMP_meta_final$SubjectID),]

pdate <- ggplot(DateInfo_clean,aes(x=Date, y=Subject)) + geom_point()
pdate

#ggsave(filename = "./Suppl.figure/subject.before.trim.pdf", height = 13, width = 6, dpi=300)

UBM_clean[1:10, 1:10]
HMP_final[1:10, 1:10]

row.names(UBM_clean) <- UBM_clean$X
UBM_clean <- UBM_clean[,-1]

UBM_otu_table <- UBM_clean[which(row.names(UBM_clean) %in% UBM_cleanmeta$KitID),]
UBM_otu_table <- UBM_otu_table[, which(colSums(UBM_otu_table) != 0)]
dim(UBM_otu_table)

rarespecieslist <- colnames(UBM_clean[,colSums(UBM_clean) < 10])
rarespeciestaxa <- subset(UBM_Taxa_clean, X %in% rarespecieslist)

rarespecieslist2 <- colnames(HMP_final[, colSums(HMP_final) < 10])
rarespeciestaxa2 <- subset(HMP_Taxa_clean, X %in% rarespecieslist2)

p <- data.frame()
for (i in 1:7){
  p[i,1] <- sum(str_count(rarespeciestaxa[,i], pattern = "Unclassified"))/4396
}

for (i in 1:7){
  p[i+7,1] <- sum(str_count(rarespeciestaxa2[,i], pattern = "Unclassified"))/28969
}

p$V2 <- rep(colnames(HMP_Taxa_clean)[1:7],2)
p$V2 <- factor(p$V2, levels = colnames(HMP_Taxa_clean)[1:7])
p$V3 <- c(rep("Stool/Nasal", 7), rep("Skin/Oral", 7))

p.ident <- ggplot(p, aes(x=V2, y=(1-V1), color=V3, group=V3)) + geom_point() + geom_line() + theme_cowplot()
p.ident <- p.ident + ggtitle("Identification rate of rare_taxa (ASV < 10)  in each phylogenetic level")+scale_color_manual(values=c("red","blue"))
p.ident

#ggsave(p.ident, filename = "./Suppl.figure/identification_rate_rare.pdf", width = 7, height = 3, dpi = 300)

row.names(UBM_Taxa_clean) <- UBM_Taxa_clean$X
UBM_Taxa_clean <- UBM_Taxa_clean[,-1]

#sequencing depth more than 1000, and remove empty ASV
UBM_OTU <- UBM_otu_table[rowSums(UBM_otu_table) > 1000,]
UBM_OTU <- UBM_OTU[,which(colSums(UBM_OTU) != 0)]

HMP_OTU <- HMP_final[rowSums(HMP_final) > 1000,]
HMP_OTU <- HMP_OTU[,which(colSums(HMP_OTU) != 0)]

rownames(HMP_OTU)
rownames(UBM_OTU)

dim(HMP_OTU)
dim(UBM_OTU)

HMP_Sample <- HMP_meta_final[which(HMP_meta_final$X %in% rownames(HMP_OTU)),]
UBM_Sample <- UBM_cleanmeta[which(UBM_cleanmeta$KitID %in% rownames(UBM_OTU)),]
row.names(UBM_Sample) <- UBM_Sample$KitID

UBM_Sample$SubjectID[str_detect(UBM_Sample$SampleID, "69-028")] <- "70-1014"
UBM_Sample[str_detect(UBM_Sample$SampleID, "69-001-6023"),]

HMP_Sample[str_detect(HMP_Sample$SampleID, "69-001-20"),]

#Change HMP sample ID
UBM_Sample$SampleID[UBM_Sample$KitID == "UBMkid114128237"] <- "69-032-1053"
UBM_Sample$SubjectID[UBM_Sample$KitID == "UBMkid114128237"] <- "69-032"
UBM_Sample$SampleID[UBM_Sample$KitID=="UBMkid460128188"] <- "69-012-09"
UBM_Sample$SampleID[UBM_Sample$KitID=="UBMkid407107574"] <- "69-001-1012"
UBM_Sample <- UBM_Sample[UBM_Sample$KitID!= "UBMkid877819249",]

#remove sequencing depth outlier based on LOF analysis， UBMkid577819081 is a duplicate sample, so added here
removeUBMlist <-  c("UBMkid381128669", "UBMkid415128106", "UBMkid424128727", "UBMkid477128903",
                    "UBMkid569127769", "UBMkid656129034", "UBMkid709127971", "UBMkid766128994", 
                    "UBMkid777128286", "UBMkid377129270", "UBMkid469129164", "UBMkid577819081")

UBM_Sample <- UBM_Sample[! UBM_Sample$KitID %in% removeUBMlist,]

HMP_Sample$SampleID[HMP_Sample$V7=="ZVTCAK9-01X"] <- "69-061-01"
HMP_Sample$SampleID[HMP_Sample$V7=="ZKFV71L-05X"] <-"70-1010-05"
HMP_Sample$SampleID[HMP_Sample$V7=="ZOZOW1T-44ST"] <-"69-001-44"
HMP_Sample$SampleID[HMP_Sample$V7=="ZOZOW1T-44N"] <-"69-001-44"
HMP_Sample$SampleID[HMP_Sample$V7=="ZOZOW1T-1015R"] <-"69-001-1015"
HMP_Sample$SampleID[HMP_Sample$V7=="ZOZOW1T-48bST"] <-"69-001-48"
HMP_Sample$SampleID[HMP_Sample$V7=="ZOZOW1T-48bN"] <-"69-001-48"
HMP_Sample$SampleID[HMP_Sample$V7=="ZOZOW1T-49ST"] <-"69-001-49"
HMP_Sample$SampleID[HMP_Sample$V7=="ZOZOW1T-49N"] <-"69-001-49"
HMP_Sample$SampleID[HMP_Sample$V7=="ZOZOW1T-53ST"] <-"69-001-53"
HMP_Sample$SampleID[HMP_Sample$V7=="ZOZOW1T-54bN"] <-"69-001-54"
HMP_Sample$SampleID[HMP_Sample$V7=="ZOZOW1T-54bST"] <-"69-001-54"
HMP_Sample$SampleID[HMP_Sample$V7=="ZOZOW1T-55ST"] <-"69-001-55"
HMP_Sample$SampleID[HMP_Sample$V7=="ZOZOW1T-55N"] <-"69-001-55"
HMP_Sample$SampleID[HMP_Sample$V7=="ZOZOW1T-56ST"] <-"69-001-56"
HMP_Sample$SampleID[HMP_Sample$V7=="ZOZOW1T-56N"] <-"69-001-56"
HMP_Sample$SampleID[HMP_Sample$V7=="ZOZOW1T-57ST"] <-"69-001-57"
HMP_Sample$SampleID[HMP_Sample$V7=="ZOZOW1T-57N"] <-"69-001-57"
HMP_Sample$SampleID[HMP_Sample$V7=="ZOZOW1T-58ST"] <-"69-001-58"
HMP_Sample$SampleID[HMP_Sample$V7=="ZOZOW1T-58N"] <-"69-001-58"
HMP_Sample$SampleID[HMP_Sample$V7=="ZOZOW1T-59ST"] <-"69-001-59"
HMP_Sample$SampleID[HMP_Sample$V7=="ZOZOW1T-59N"] <-"69-001-59"
HMP_Sample$SampleID[HMP_Sample$V7=="ZOZOW1T-6013R"] <-"69-001-6013"
HMP_Sample$SampleID[HMP_Sample$V7=="ZOZOW1T-6013A"] <-"69-001-6013_TR2"
HMP_Sample$SampleID[HMP_Sample$V7=="ZOZOW1T-6015R"] <-"69-001-6015"
HMP_Sample$SampleID[HMP_Sample$V7=="ZOZOW1T-6016R"] <-"69-001-6016"
HMP_Sample$SampleID[HMP_Sample$V7=="ZOZOW1T-6021R"] <-"69-001-6021"
HMP_Sample$SampleID[HMP_Sample$V7=="ZOZOW1T-6023bR"] <-"69-001-6023"
HMP_Sample$SampleID[HMP_Sample$V7=="ZOZOW1T-61R"] <-"69-001-61"
HMP_Sample$SampleID[HMP_Sample$V7=="ZOZOW1T-63R"] <-"69-001-63"
HMP_Sample$SampleID[HMP_Sample$V7=="ZOZOW1T-65R"] <-"69-001-65"
HMP_Sample$SampleID[HMP_Sample$V7=="ZOZOW1T-66R"] <-"69-001-66"
HMP_Sample$SampleID[HMP_Sample$V7=="ZOZOW1T-68R"] <-"69-001-68"
HMP_Sample$SampleID[HMP_Sample$V7=="ZOZOW1T-7013A"] <-"69-001-7013"
HMP_Sample$SampleID[HMP_Sample$V7=="ZOZOW1T-71R"] <-"69-001-71"
HMP_Sample$SampleID[HMP_Sample$V7=="ZOZOW1T-80R"] <-"69-001-80"
HMP_Sample$SampleID[HMP_Sample$V7=="ZOZOW1T-83R"] <-"69-001-83"
HMP_Sample$SampleID[HMP_Sample$V7=="ZOZOW1T-84R"] <-"69-001-84"

HMP_Sample$batch[HMP_Sample$X=="BCIXQPED"] <-"B86HB"

#two sample have incorrect collection date
HMP_Sample <- subset(HMP_Sample, X != "RHAOCUPN")
HMP_Sample <- subset(HMP_Sample, X != "PRLMUBZT")

#one sample has too deep sequencing depth as outlier
HMP_Sample <- subset(HMP_Sample, X != "KTCVOQCL")

HMP_Sample[str_detect(HMP_Sample$SubjectID, "70-1014"),]
UBM_Sample[str_detect(UBM_Sample$SubjectID, "70-1014"),]
UBM_Sample[str_detect(UBM_Sample$KitID, "UBMkid942104956"),]

#remove one duplicated sample



UBM_Sample <- merge(UBM_Sample, DateInfo, by="SampleID")
UBM_Sample <- select(UBM_Sample, KitID,SampleID,SampleType,batch,Date.y, Subject)
row.names(UBM_Sample) <- UBM_Sample$KitID
colnames(UBM_Sample)[6] <- "SubjectID"
colnames(UBM_Sample)[5] <- "Date"
UBM_Sample <- subset(UBM_Sample, SampleType != "UNK")
UBM_Sample$SampleType[UBM_Sample$SampleType == "Mouth"] <- "Oral"
UBM_Sample$Date <- as.Date(UBM_Sample$Date, format = "%m/%d/%y")
UBM_Sample <- UBM_Sample[!str_detect(UBM_Sample$SubjectID, "s"),]

HMP_Sample <- merge(HMP_Sample, DateInfo, by="SampleID")
HMP_Sample <- HMP_Sample[!is.na(HMP_Sample$Date),]
row.names(HMP_Sample) <- HMP_Sample$X.x
HMP_Sample <- subset(HMP_Sample, duplicated==F)
HMP_Sample <- select(HMP_Sample, X.x, SampleID, Subject, Date,V1, V3,V7,V9)
colnames(HMP_Sample) <- c("RandomID", "SampleID", "SubjectID", "Date", "JAXID", "SampleType", "ori.ID", "batch")
HMP_Sample$Date <- as.Date(HMP_Sample$Date, format = "%m/%d/%y")

HMP_Sample[str_detect(HMP_Sample$SampleID, "69-028"),]

#Overview of the sample
HMP_meta <- HMP_Sample
UBM_meta <- UBM_Sample

HMP_meta_ST <- subset(HMP_meta, SampleType == "ST")
HMP_meta_NS <- subset(HMP_meta, SampleType == "NS")

UBM_meta_OR <- subset(UBM_meta, SampleType == "Oral")
UBM_meta_SK <- subset(UBM_meta, SampleType == "Skin")

HMP_meta2 <- select(HMP_meta, RandomID,SampleID,SampleType,batch,Date,SubjectID)
colnames(HMP_meta2)[1] <- "KitID"
Combined_meta <- rbind(HMP_meta2, UBM_meta)

p <- ggplot(Combined_meta, aes(x=Date, y=SubjectID, color=SampleType)) + geom_point() #+ geom_jitter()
p <- p + theme_cowplot() + facet_wrap(.~SampleType, ncol=4)
p

sampleFreq <- merge(merge(data.frame(table(HMP_meta_ST$SubjectID)),data.frame(table(HMP_meta_NS$SubjectID)), by="Var1", all=T), 
                    merge(data.frame(table(UBM_meta_OR$SubjectID)),data.frame(table(UBM_meta_SK$SubjectID)), by="Var1",all=T), by= "Var1",all=T) 

colnames(sampleFreq) <- c("SubjectID", "ST", "NS", "OR", "SK")
sampleFreq[is.na(sampleFreq)] <- 0
listA <- as.character(sampleFreq$SubjectID[! sampleFreq$SubjectID %in% sc$SubjectID])
listB <- as.character(sampleFreq$SubjectID[rowSums(sampleFreq[,2:5] < 3) == 4])
removelist <- unique(c(listA,listB))
removelist <- c(removelist, "69-004")

HMP_Sample <- HMP_Sample[!HMP_Sample$SubjectID %in% removelist,]
UBM_Sample <- UBM_Sample[!UBM_Sample$SubjectID %in% removelist,]

keeplist <- as.character(sampleFreq$SubjectID)[!as.character(sampleFreq$SubjectID) %in% removelist]

datadura <- data.frame()
for (i in 1:length(keeplist)){
  print(i)
  data <- subset(Combined_meta, SubjectID == keeplist[i])
  datadura[i,1] <- keeplist[i]
  datadura[i,2] <- min(data$Date)
  datadura[i,3] <- max(data$Date)
  datadura[i,4] <- (as.numeric(max(data$Date) - min(data$Date)))
}

colnames(datadura) <- c("SubjectID", "Start_Date", "End_Date", "Duration_Day")
#write.csv(file = "./ori_meta_table/SubjectDuration.csv",datadura)

HMP_TAXA <- HMP_Taxa_clean[which(HMP_Taxa_clean$X %in% colnames(HMP_OTU)),1:7]
UBM_TAXA <- UBM_Taxa_clean[which(rownames(UBM_Taxa_clean) %in% colnames(UBM_OTU)),]

HMP_OTU <- HMP_OTU[which(rownames(HMP_OTU) %in% HMP_Sample$RandomID),]
HMP_OTU <- HMP_OTU[,which(colSums(HMP_OTU) != 0)]

OTU = otu_table(as.matrix(HMP_OTU),taxa_are_rows = F)
TAX = tax_table(as.matrix(HMP_TAXA))
map1 = sample_data(HMP_Sample)
physeq_HMP = phyloseq(OTU,TAX,map1)

UBM_OTU <- UBM_OTU[which(rownames(UBM_OTU) %in% UBM_Sample$KitID),]
UBM_OTU <- UBM_OTU[,which(colSums(UBM_OTU) != 0)]

OTU1 = otu_table(as.matrix(UBM_OTU),taxa_are_rows = F)
TAX1 = tax_table(as.matrix(UBM_TAXA))
map2 = sample_data(UBM_Sample)
physeq_UBM = phyloseq(OTU1,TAX1,map2)

save(physeq_HMP, physeq_UBM, file = "./Robject/PhyloseqObject.RData")





