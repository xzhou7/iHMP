#integret the HMP and UMB metadata
#Xin Zhou 
#03-18-2021

library(dplyr)
library(ggplot2)
library(stringr)
library(factoextra)

setwd("~/Box/XinZhouFiles/Projects/ZXP2_HMP2_Core/Analysis/")
HMP_meta <- read.csv("./ori_meta_table/Metadata.HMP.csv", header = T)
UBM_meta <- read.csv("./ori_meta_table/UbiomSAMPLEmeta.CSV", header = T)
batchlist <- read.csv("./ori_meta_table/batchinglist.csv", header = T)
lipidomicsmeta <- read.csv("./ori_meta_table/Lipidomics_clinical.csv", header = T)
JAXmeta <- read.csv("./ori_meta_table/hmp2_infection_metadata_0323.csv", header = T)

load("../OLD HMP Files/Analysis/Analysis/Revision_MultiOmes_0509.RData")

lipidomicsmeta$CollectionDate <- as.Date(lipidomicsmeta$CollectionDate, format = "%m/%d/%y")

lipidomicsmeta[lipidomicsmeta$CollectionDate!="2010-01-12",]

lipidomicsmeta$HSCRP[lipidomicsmeta$HSCRP == "<0.2"] <- "0.19"
lipidomicsmeta$HSCRP <- as.numeric(lipidomicsmeta$HSCRP)

lipidomicsmeta$AG[lipidomicsmeta$AG == "<2"] <- "1.5"
lipidomicsmeta$AG <- as.numeric(lipidomicsmeta$AG)

lipidomicsmeta$ALCRU[lipidomicsmeta$ALCRU == "<10.00"] <- "9.9"
lipidomicsmeta$ALCRU <- as.numeric(lipidomicsmeta$ALCRU)

lipidomicsmeta$ALT[lipidomicsmeta$ALT == "<10"] <- "9"
lipidomicsmeta$ALT[lipidomicsmeta$ALT == "<15"] <- "9"
lipidomicsmeta$ALT <- as.numeric(lipidomicsmeta$ALT)

lipidomicsmeta$AST[lipidomicsmeta$AST == "<10"] <- "9"
lipidomicsmeta$AST <- as.numeric(lipidomicsmeta$AST)

lipidomicsmeta$IGM[lipidomicsmeta$IGM == "<6"] <- "5"
lipidomicsmeta$IGM <- as.numeric(lipidomicsmeta$IGM)

lipidomicsmeta$TBIL[lipidomicsmeta$TBIL == "<0.5"] <- "0.4"
lipidomicsmeta$TBIL <- as.numeric(lipidomicsmeta$TBIL)

lipidomicsmeta$UALB[lipidomicsmeta$UALB == "<5"] <- "4"
lipidomicsmeta$UALB <- as.numeric(lipidomicsmeta$UALB)

lipidomicsmeta$EGFR[lipidomicsmeta$EGFR == ">60"] <- "61"
lipidomicsmeta$EGFR <- as.numeric(lipidomicsmeta$EGFR)

lipidomicsmeta <- select(lipidomicsmeta, -UALBCR,-Collection.Date )

for (i in 1:54) {
  print(colnames(lipidomicsmeta)[i])
  print(is.numeric(lipidomicsmeta[,i]))
}

p.a1c <- ggplot(lipidomicsmeta, aes(x=CollectionDate, y=EGFR, group=Subject)) + geom_line(na.rm= TRUE)
p.a1c

p.crp <- ggplot(lipidomicsmeta, aes(x=CollectionDate, y=HSCRP, group=Subject)) + geom_line(na.rm= TRUE)
p.crp

#3.5 to 5 as cutoffs
p.ldlhdl <- ggplot(lipidomicsmeta, aes(x=CollectionDate, y=LDLHDL, group=Subject)) + geom_line(na.rm= TRUE)
p.ldlhdl

batchlist <- batchlist[,-1]

#get subject ID by split sample ID varible
HMP_meta[,12:15] <- as.data.frame(str_split(HMP_meta$V7, "-", simplify = T))

#replace random ID with  subject ID
HMP_meta$SubjectID <- sc$SubjectID[match(unlist(HMP_meta$V1.1), sc$rand_subject_id)]  
HMP_meta[is.na(HMP_meta$SubjectID),]

#list unmapped subject
HMP_unmapped <- HMP_meta$V1.1[is.na(HMP_meta$SubjectID)]
HMP_unmapped

sum(table(HMP_meta$V1.1))
sum(table(HMP_meta$SubjectID)) 

#get sample ID after replacement
HMP_meta$SampleID <- paste(HMP_meta$SubjectID, HMP_meta$V2.1, sep = "-")

#seperate Stool sample with Nasal sample
HMP_ST_meta <- subset(HMP_meta, V3== "ST")
HMP_NS_meta <- subset(HMP_meta, V3== "NS")
table(duplicated(HMP_ST_meta$SampleID))
table(duplicated(HMP_NS_meta$SampleID))

HMP_ST_meta[duplicated(HMP_ST_meta$SampleID),]
#check duplicated sample (techniqual control)
subset(HMP_ST_meta, SampleID=="69-001-1022")

table(HMP_meta$V4.1)
subset(HMP_meta, V3.1==1 | V4.1==1)

table(HMP_meta$SampleID)[table(HMP_meta$SampleID)>2]

duplicatedST <- names(table(HMP_ST_meta$SampleID)[table(HMP_ST_meta$SampleID)>1])
duplicatedNS <- names(table(HMP_NS_meta$SampleID)[table(HMP_NS_meta$SampleID)>1])

subset(HMP_ST_meta, SampleID %in% duplicatedST)

HMP_ST_meta$duplicated <- duplicated(HMP_ST_meta$SampleID)
HMP_NS_meta$duplicated <- duplicated(HMP_NS_meta$SampleID)

HMP_ST_meta$metadata <- HMP_ST_meta$SampleID %in% ls$SampleID
HMP_NS_meta$metadata <- HMP_NS_meta$SampleID %in% ls$SampleID

#this is to check if the metadata would cover all samples, seems yes
#ls[,ls$SampleID=="69-064-6036"]
#ls[str_detect(ls$SampleID,"69-091"),]
#HMP_ST_meta[str_detect(HMP_ST_meta$SampleID, "69-091"),]
#setdiff(ls$SampleID, HMP_NS_meta$SampleID)
#setdiff(HMP_NS_meta$SampleID,ls$SampleID)

UBM_OL_meta <- subset(UBM_meta, SampleType=="Mouth")
UBM_SK_meta <- subset(UBM_meta, SampleType=="Skin")

length(intersect(UBM_OL_meta$SampleID,ls$SampleID))
length(intersect(UBM_SK_meta$SampleID,ls$SampleID))

length(intersect(HMP_NS_meta$SampleID,ns.df$HostSampleID))
setdiff(HMP_ST_meta$SampleID,st.df$HostSampleID)
setdiff(st.df$HostSampleID,HMP_ST_meta$SampleID)

HMP_meta_final <- rbind(HMP_ST_meta, HMP_NS_meta)
row.names(HMP_meta_final) <- HMP_meta_final$X
table(duplicated(HMP_meta_final$V1))

#UBM_meta[UBM_meta$KitID == "UBMkid908819636", ]
UBM_meta[UBM_meta$SampleID ==  "70-1003-13", ]

ls[str_detect(ls$SampleID, "69-105"), ]
ls[str_detect(ls$SubjectID, "69-106"), ]

UBM_meta[str_detect(UBM_meta$SampleID, "69-028"),]
HMP_meta_final[str_detect(HMP_meta_final$SampleID, "69-028"),]

HMP_unmapped
HMP_meta_final[str_detect(HMP_meta_final$V1.1, "ZMWEIX1"),]

#write.csv(HMP_meta_final, file = "~/Desktop/HMP_meta.csv")

######
#check 24 unmapped sample
#table(HMP_unmapped)
#for (i in 1:105){
#  print(i)
#  print(do.call(setdiff, strsplit(c("ZMWEIX1", sc$rand_subject_id[i]), split = "")))
#}
#####

table(UBM_meta$SampleID)[table(UBM_meta$SampleID) >2]
#For UB sample, adjust metadata, TR means it is a technical replicate
UBM_meta$SampleID[UBM_meta$KitID == "UBMkid987128845"] <- "69-001-7024_TR"
UBM_meta$SampleID[UBM_meta$KitID == "UBMkid432128562"] <- "69-028-6014_2"
UBM_meta$SampleID[UBM_meta$KitID == "UBMkid265104610"] <- "69-028-6015_2"
UBM_meta$SampleID[UBM_meta$KitID == "UBMkid742104607"] <- "69-028-6034_2"
UBM_meta$SampleID[UBM_meta$KitID == "UBMkid815128098"] <- "69-028-6034_2"
UBM_meta$SampleID[UBM_meta$KitID == "UBMkid951128874"] <- "69-123-3014_TR"
UBM_meta$SampleID[UBM_meta$KitID == "UBMkid992129229"] <- "69-033-08_2"
UBM_meta$SampleID[UBM_meta$KitID == "UBMkid798127504"] <- "69-033-08_2"
UBM_meta$SampleID[UBM_meta$KitID == "UBMkid406129149"] <- "69-036-04_2"
UBM_meta$SampleID[UBM_meta$KitID == "UBMkid413107467"] <- "69-036-04_2"
UBM_meta$SampleID[UBM_meta$KitID == "UBMkid535128181"] <- "69-053-6031_2"
UBM_meta$SampleID[UBM_meta$KitID == "UBMkid841105105"] <- "69-053-6031_2"
UBM_meta$SampleID[UBM_meta$KitID == "UBMkid606127816"] <- "69-064-6014_2"
UBM_meta$SampleID[UBM_meta$KitID == "UBMkid576105155"] <- "69-064-6014_2"
UBM_meta$SampleID[UBM_meta$KitID == "UBMkid656129404"] <- "69-064-6024_2"
UBM_meta$SampleID[UBM_meta$KitID == "UBMkid949105681"] <- "69-064-6024_2"
UBM_meta$SampleID[UBM_meta$KitID == "UBMkid626127817"] <- "69-066-05_2"
UBM_meta$SampleID[UBM_meta$KitID == "UBMkid808129075"] <- "69-069-6032_2"
UBM_meta$SampleID[UBM_meta$KitID == "UBMkid850105490"] <- "69-069-6032_2"
UBM_meta$SampleID[UBM_meta$KitID == "UBMkid510128587"] <- "69-069-6033_2"
UBM_meta$SampleID[UBM_meta$KitID == "UBMkid797105461"] <- "69-069-6033_2"
UBM_meta$SampleID[UBM_meta$KitID == "UBMkid712128617"] <- "69-069-6034_2"
UBM_meta$SampleID[UBM_meta$KitID == "UBMkid726105611"] <- "69-069-6034_2"
UBM_meta$SampleID[UBM_meta$KitID == "UBMkid281105195"] <- "69-069-6035_2"
UBM_meta$SampleID[UBM_meta$KitID == "UBMkid509128684"] <- "69-069-6035_2"
UBM_meta$SampleID[UBM_meta$KitID == "UBMkid673128605"] <- "69-069-6036_2"
UBM_meta$SampleID[UBM_meta$KitID == "UBMkid634105596"] <- "69-069-6036_2"
UBM_meta$SampleID[UBM_meta$KitID == "UBMkid520128856"] <- "69-074-6033_2"
UBM_meta$SampleID[UBM_meta$KitID == "UBMkid265105003"] <- "69-074-6033_2"
UBM_meta$SampleID[UBM_meta$KitID == "UBMkid724129219"] <- "69-074-6014_2"
UBM_meta$SampleID[UBM_meta$KitID == "UBMkid286105109"] <- "69-074-6014_2"
UBM_meta$SampleID[UBM_meta$KitID == "UBMkid962129382"] <- "69-090-1031_TR"
UBM_meta$SampleID[UBM_meta$KitID == "UBMkid550128730"] <- "69-095-05_2"
UBM_meta$SampleID[UBM_meta$KitID == "UBMkid716105406"] <- "69-095-05_2"
UBM_meta$SampleID[UBM_meta$KitID == "UBMkid983129431"] <- "69-099-02_2"
UBM_meta$SampleID[UBM_meta$KitID == "UBMkid983105280"] <- "69-099-02_2"
UBM_meta$SampleID[UBM_meta$KitID == "UBMkid772127985"] <- "69-111-6015_2"
UBM_meta$SampleID[UBM_meta$KitID == "UBMkid197107357"] <- "69-111-6015_2"
UBM_meta$SampleID[UBM_meta$KitID == "UBMkid153129358"] <- "69-113-03_2"
UBM_meta$SampleID[UBM_meta$KitID == "UBMkid681105562"] <- "69-113-03_2"
UBM_meta$SampleID[UBM_meta$KitID == "UBMkid812128646"] <- "69-114-03_2"
UBM_meta$SampleID[UBM_meta$KitID == "UBMkid665105530"] <- "69-116-03_2"
UBM_meta$SampleID[UBM_meta$KitID == "UBMkid215129073"] <- "69-116-03_2"
UBM_meta$SampleID[UBM_meta$KitID == "UBMkid249127949"] <- "70-1003-07_TR"
UBM_meta$SampleID[UBM_meta$KitID == "UBMkid535104878"] <- "70-1003-13_TR"
UBM_meta$SampleID[UBM_meta$KitID == "UBMkid201128600"] <- "70-1006-09_2"
UBM_meta$SampleID[UBM_meta$KitID == "UBMkid794105575"] <- "70-1006-09_2"
UBM_meta$SampleID[UBM_meta$KitID == "UBMkid746128962"] <- "70-1006-2012_2"
UBM_meta$SampleID[UBM_meta$KitID == "UBMkid165105021"] <- "70-1006-2012_2"
UBM_meta$SampleID[UBM_meta$KitID == "UBMkid502104873"] <- "70-1008-2024_UNK"
UBM_meta$SampleID[UBM_meta$KitID == "UBMkid138127920"] <- "70-1010-04_2"
UBM_meta$SampleID[UBM_meta$KitID == "UBMkid181107507"] <- "70-1010-04_2"
UBM_meta$SampleID[UBM_meta$KitID == "UBMkid234129400"] <- "69-001-1012_2"
UBM_meta$SampleID[UBM_meta$KitID == "UBMkid407107574"] <- "69-001-1012_2"
UBM_meta$SampleID[UBM_meta$KitID == "UBMkid908819636"] <- "69-001-6013_TR"
UBM_meta$SampleID[UBM_meta$KitID == "UBMkid570105026"] <- "70-1004-05"
UBM_meta$SampleID[UBM_meta$KitID == "UBMkid983128793"] <- "69-046-10"

#match 69-001 with date
table(UBM_meta$Date[UBM_meta$SampleID ==  "69-001"] %in% ls$CollectionDate[str_detect(ls$SampleID, "69-001")])
UBM_meta$Date[UBM_meta$SampleID ==  "69-001"][!UBM_meta$Date[UBM_meta$SampleID ==  "69-001"] %in% ls$CollectionDate[str_detect(ls$SampleID, "69-001")]]
UBM_no69001 <- subset(UBM_meta, SampleID !=  "69-001")
UBM_no69001[str_detect(UBM_no69001$SampleID, "69-001"),]
UBM69001 <- subset(UBM_meta, SampleID == "69-001")
ls_69001 <- subset(ls, SubjectID== "69-001")
UBM69001$SampleID <- ls_69001$SampleID[match(UBM69001$Date, ls_69001$CollectionDate)]
UBM69001$Date <- as.Date(UBM69001$Date, format = "%m/%d/%y")
UBM_no69001$Date <- as.Date(UBM_no69001$Date,format = "%m/%d/%y")
ls_69001$CollectionDate <- as.Date(ls_69001$CollectionDate,format = "%m/%d/%y")

#match sample ID whose date are close
#20130605/20130606
UBM69001$SampleID[UBM69001$KitID == "UBMkid930819174"] <- "69-001-61"
#20140911/20140912
UBM69001$SampleID[UBM69001$KitID == "UBMkid877819249"] <- "69-001-4011"
#20150125/20150126
UBM69001$SampleID[UBM69001$KitID == "UBMkid577819081"] <-"69-001-6023"

#these IDs do not have metadata in ls
UBM69001_clean <- UBM69001[! UBM69001$KitID %in% c("UBMkid751819609", "UBMkid108819837","UBMkid187819213","UBMkid448819207","UBMkid467819525","UBMkid548819477"),]
UBM_cleanmeta <- rbind(UBM69001_clean,UBM_no69001)
table(UBM_cleanmeta$SampleType)
UBM_cleanmeta <- merge(UBM_cleanmeta, batchlist, by = "KitID")
table(gsub("^([^-]*-[^-]*)-.*$", "\\1", UBM_cleanmeta$SampleID))
UBM_cleanmeta$SubjectID <- gsub("^([^-]*-[^-]*)-.*$", "\\1", UBM_cleanmeta$SampleID)
table(UBM_cleanmeta$SubjectID)


setdiff(JAXmeta$SampleID, ls$SampleID)
setdiff(ls$SampleID, JAXmeta$SampleID)


save(UBM_cleanmeta, HMP_meta_final, JAXmeta, lipidomicsmeta, ls, metb.curated, metbcr.df, rnaseq.log.df, sc, file = "./Robject/Meta_MultiOmes0413.RData")


