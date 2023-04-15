# make metadata for HMP and uBIOME data 
library(stringr)
library(ggplot2)
library(cowplot)
library(dplyr)
library(gplots)
library(stringi)

setwd("~/Box/XinZhouFiles/Projects/ZXP2_HMP2_Core/")

#set input file
#################################################################################################################################
HMPotu <- read.csv("~/Desktop/HMP_DADA2_result/otu.table_20200207.csv", row.names = 1, header = T)
HMPtaxa <- read.csv("~/Desktop/HMP_DADA2_result/tax.table_20200207.csv", row.names = 1, header = T)

UBMotu <- read.csv("~/Desktop/UBIOME_DADA2_result_V18/SkinOralotu.table_F_20210217.csv", header = T)
UBMtaxa <- read.csv("~/Desktop/UBIOME_DADA2_result_V18/SkinOraltax.table_F_20210217.csv", header = T)

OldHMPiden <- read.csv("./OLD HMP Files/Analysis/Analysis/lefseTAX.csv", header = T)
#################################################################################################################################

p <- data.frame()
for (i in 1:7){
  p[i,1] <- sum(table(HMPtaxa[i]))/61610
}
for (i in 1:7){
  p[i+7,1] <- sum(table(UBMtaxa[i]))/14513
}


p$V2 <- rep(colnames(HMPtaxa),2)
p$V2 <- factor(p$V2, levels = colnames(HMPtaxa))
p$V3 <- c(rep("Stool/Nasal", 7), rep("Skin/Oral", 7))

for (i in 1:6){
  p [i+14, 1] <- 1 -  (sum(str_count(OldHMPiden[,i+1],"unclassified"))/1953)
}

p$V2[15:20] <- p$V2[8:14]
p$V3[15:20] <- "Stool_2019"

p.ident <- ggplot(p, aes(x=V2, y=V1, color=V3, group=V3)) + geom_point() + geom_line() + theme_cowplot()
p.ident <- p.ident + ggtitle("Identification rate of data in each phylogenetic level")
p.ident

ggsave(p.ident, filename = "./Analysis/Suppl.figure/identification_rate.pdf",width = 7, height = 3, dpi = 300)

##############################################################################
# remove duplicated sample
##############################################################################
UBM_dup <- UBMotu$X[duplicated(UBMotu$X)]
UBM_duplicated <- UBMotu[which(UBMotu$X %in% UBM_dup),]
UBM_clean <- UBMotu[which(! UBMotu$X %in% UBM_dup),]

#check how many ASV are lost due to subset samples
table(colSums(UBMotu[,-1])==0)
table(colSums(UBM_clean[,-1])==0)
UBM_final <- UBM_clean[,colSums(UBM_clean[,-1])!=0]
218/14513

#same for HMP data
HMPmeta <- as.data.frame(str_split(row.names(HMPotu), pattern = "_", simplify = T))
HMPmeta$V9 <- gsub("\\..*","",HMPmeta$V8)
DuplicatedJAXID <- HMPmeta$V1[duplicated(HMPmeta$V1)]
#disect metadata sample name
DuplicatedHMPmeta <- HMPmeta[which(HMPmeta$V1 %in% DuplicatedJAXID), ]

table(DuplicatedHMPmeta$V1, DuplicatedHMPmeta$V3)  
DuplicatedHMPmeta$V10 <- paste(DuplicatedHMPmeta$V1,DuplicatedHMPmeta$V3, sep=".")

#emove sample with duplicated JAXID (V10) in HMPotu file, should be 2086 samples left
HMP_clean <- HMPotu[-as.numeric(row.names(DuplicatedHMPmeta[which(duplicated(DuplicatedHMPmeta$V10)),])),]

HMPmeta2 <- HMPmeta
row.names(HMPmeta2) <- row.names(HMPotu)
HMPmeta_clean <- HMPmeta2[-as.numeric(row.names(DuplicatedHMPmeta[which(duplicated(DuplicatedHMPmeta$V10)),])),]
table(duplicated(HMPmeta_clean))
#check how many ASV are lost dur to subset samples
table(colSums(HMPotu)==0)
table(colSums(HMP_clean)==0)
HMP_final <- HMP_clean[,colSums(HMP_clean)!=0]

#create random ID for HMP data to deidentify patient
HMPmeta_clean$ID <- stri_rand_strings(nrow(HMP_clean), 8, pattern="[A-Z]")

rownames(HMP_clean) <- HMPmeta_clean$ID
rownames(HMPmeta_clean) <- HMPmeta_clean$ID
dim(HMPmeta_clean)

#replace HMP otu ID
row.names(HMP_final) <- HMPmeta_clean$ID

##############################################################################
# add unclassified to TAX table
##############################################################################
HMP_Taxa <- mutate_all(HMPtaxa, as.character)

HMP_Taxa$Species[is.na(HMP_Taxa$Species)&!is.na(HMP_Taxa$Genus)] <- paste("Unclassified", HMP_Taxa$Genus[is.na(HMP_Taxa$Species)&!is.na(HMP_Taxa$Genus)], sep="_")

HMP_Taxa$Species[is.na(HMP_Taxa$Species)&is.na(HMP_Taxa$Genus)&!is.na(HMP_Taxa$Family)] <- paste("Unclassified", HMP_Taxa$Family[is.na(HMP_Taxa$Species)&is.na(HMP_Taxa$Genus)&!is.na(HMP_Taxa$Family)], sep="_")
HMP_Taxa$Genus[is.na(HMP_Taxa$Genus)&!is.na(HMP_Taxa$Family)] <- paste("Unclassified", HMP_Taxa$Family[is.na(HMP_Taxa$Genus)&!is.na(HMP_Taxa$Family)], sep="_")

HMP_Taxa$Species[is.na(HMP_Taxa$Species)&is.na(HMP_Taxa$Genus)&is.na(HMP_Taxa$Family)&!is.na(HMP_Taxa$Order)] <- paste("Unclassified", HMP_Taxa$Order[is.na(HMP_Taxa$Species)&is.na(HMP_Taxa$Genus)&is.na(HMP_Taxa$Family)&!is.na(HMP_Taxa$Order)], sep="_")
HMP_Taxa$Genus[is.na(HMP_Taxa$Genus)&is.na(HMP_Taxa$Family)&!is.na(HMP_Taxa$Order)] <- paste("Unclassified", HMP_Taxa$Order[is.na(HMP_Taxa$Genus)&is.na(HMP_Taxa$Family)&!is.na(HMP_Taxa$Order)], sep="_")
HMP_Taxa$Family[is.na(HMP_Taxa$Family)&!is.na(HMP_Taxa$Order)] <- paste("Unclassified", HMP_Taxa$Order[is.na(HMP_Taxa$Family)&!is.na(HMP_Taxa$Order)], sep="_")

HMP_Taxa$Species[is.na(HMP_Taxa$Species)&is.na(HMP_Taxa$Genus)&is.na(HMP_Taxa$Family)&is.na(HMP_Taxa$Order)&!is.na(HMP_Taxa$Class)] <- paste("Unclassified", HMP_Taxa$Class[is.na(HMP_Taxa$Species)&is.na(HMP_Taxa$Genus)&is.na(HMP_Taxa$Family)&is.na(HMP_Taxa$Order)&!is.na(HMP_Taxa$Class)], sep="_")
HMP_Taxa$Genus[is.na(HMP_Taxa$Genus)&is.na(HMP_Taxa$Family)&is.na(HMP_Taxa$Order)&!is.na(HMP_Taxa$Class)] <- paste("Unclassified", HMP_Taxa$Class[is.na(HMP_Taxa$Genus)&is.na(HMP_Taxa$Family)&is.na(HMP_Taxa$Order)&!is.na(HMP_Taxa$Class)], sep="_")
HMP_Taxa$Family[is.na(HMP_Taxa$Family)&is.na(HMP_Taxa$Order)&!is.na(HMP_Taxa$Class)] <- paste("Unclassified", HMP_Taxa$Class[is.na(HMP_Taxa$Family)&is.na(HMP_Taxa$Order)&!is.na(HMP_Taxa$Class)], sep="_")
HMP_Taxa$Order[is.na(HMP_Taxa$Order)&!is.na(HMP_Taxa$Class)] <- paste("Unclassified", HMP_Taxa$Class[is.na(HMP_Taxa$Order)&!is.na(HMP_Taxa$Class)], sep="_")

HMP_Taxa$Species[is.na(HMP_Taxa$Species)&is.na(HMP_Taxa$Genus)&is.na(HMP_Taxa$Family)&is.na(HMP_Taxa$Order)&is.na(HMP_Taxa$Class)&!is.na(HMP_Taxa$Phylum)] <- paste("Unclassified", HMP_Taxa$Phylum[is.na(HMP_Taxa$Species)&is.na(HMP_Taxa$Genus)&is.na(HMP_Taxa$Family)&is.na(HMP_Taxa$Order)&is.na(HMP_Taxa$Class)&!is.na(HMP_Taxa$Phylum)], sep="_")
HMP_Taxa$Genus[is.na(HMP_Taxa$Genus)&is.na(HMP_Taxa$Family)&is.na(HMP_Taxa$Order)&is.na(HMP_Taxa$Class)&!is.na(HMP_Taxa$Phylum)] <- paste("Unclassified", HMP_Taxa$Phylum[is.na(HMP_Taxa$Genus)&is.na(HMP_Taxa$Family)&is.na(HMP_Taxa$Order)&is.na(HMP_Taxa$Class)&!is.na(HMP_Taxa$Phylum)], sep="_")
HMP_Taxa$Family[is.na(HMP_Taxa$Family)&is.na(HMP_Taxa$Order)&is.na(HMP_Taxa$Class)&!is.na(HMP_Taxa$Phylum)] <- paste("Unclassified", HMP_Taxa$Phylum[is.na(HMP_Taxa$Family)&is.na(HMP_Taxa$Order)&is.na(HMP_Taxa$Class)&!is.na(HMP_Taxa$Phylum)], sep="_")
HMP_Taxa$Order[is.na(HMP_Taxa$Order)&is.na(HMP_Taxa$Class)&!is.na(HMP_Taxa$Phylum)] <- paste("Unclassified", HMP_Taxa$Phylum[is.na(HMP_Taxa$Order)&is.na(HMP_Taxa$Class)&!is.na(HMP_Taxa$Phylum)], sep="_")
HMP_Taxa$Class[is.na(HMP_Taxa$Class)&!is.na(HMP_Taxa$Phylum)]<- paste("Unclassified", HMP_Taxa$Phylum[is.na(HMP_Taxa$Class)&!is.na(HMP_Taxa$Phylum)], sep="_")

HMP_Taxa$Species[is.na(HMP_Taxa$Species)&is.na(HMP_Taxa$Genus)&is.na(HMP_Taxa$Family)&is.na(HMP_Taxa$Order)&is.na(HMP_Taxa$Class)&is.na(HMP_Taxa$Phylum)&!is.na(HMP_Taxa$Kingdom)] <- paste("Unclassified", HMP_Taxa$Kingdom[is.na(HMP_Taxa$Species)&is.na(HMP_Taxa$Genus)&is.na(HMP_Taxa$Family)&is.na(HMP_Taxa$Order)&is.na(HMP_Taxa$Class)&is.na(HMP_Taxa$Phylum)&!is.na(HMP_Taxa$Kingdom)], sep="_")
HMP_Taxa$Genus[is.na(HMP_Taxa$Genus)&is.na(HMP_Taxa$Family)&is.na(HMP_Taxa$Order)&is.na(HMP_Taxa$Class)&is.na(HMP_Taxa$Phylum)&!is.na(HMP_Taxa$Kingdom)] <- paste("Unclassified", HMP_Taxa$Kingdom[is.na(HMP_Taxa$Genus)&is.na(HMP_Taxa$Family)&is.na(HMP_Taxa$Order)&is.na(HMP_Taxa$Class)&is.na(HMP_Taxa$Phylum)&!is.na(HMP_Taxa$Kingdom)], sep="_")
HMP_Taxa$Family[is.na(HMP_Taxa$Family)&is.na(HMP_Taxa$Order)&is.na(HMP_Taxa$Class)&is.na(HMP_Taxa$Phylum)&!is.na(HMP_Taxa$Kingdom)] <- paste("Unclassified", HMP_Taxa$Kingdom[is.na(HMP_Taxa$Family)&is.na(HMP_Taxa$Order)&is.na(HMP_Taxa$Class)&is.na(HMP_Taxa$Phylum)&!is.na(HMP_Taxa$Kingdom)], sep="_")
HMP_Taxa$Order[is.na(HMP_Taxa$Order)&is.na(HMP_Taxa$Class)&is.na(HMP_Taxa$Phylum)&!is.na(HMP_Taxa$Kingdom)] <- paste("Unclassified", HMP_Taxa$Kingdom[is.na(HMP_Taxa$Order)&is.na(HMP_Taxa$Class)&is.na(HMP_Taxa$Phylum)&!is.na(HMP_Taxa$Kingdom)], sep="_")
HMP_Taxa$Class[is.na(HMP_Taxa$Class)&is.na(HMP_Taxa$Phylum)&!is.na(HMP_Taxa$Kingdom)] <- paste("Unclassified", HMP_Taxa$Kingdom[is.na(HMP_Taxa$Class)&is.na(HMP_Taxa$Phylum)&!is.na(HMP_Taxa$Kingdom)], sep="_")
HMP_Taxa$Phylum[is.na(HMP_Taxa$Phylum)&!is.na(HMP_Taxa$Kingdom)] <- paste("Unclassified", HMP_Taxa$Kingdom[is.na(HMP_Taxa$Phylum)&!is.na(HMP_Taxa$Kingdom)], sep="_")

HMP_Taxa[is.na(HMP_Taxa$Species),]
HMP_Taxa[10990:10999,]
HMP_Taxa[is.na(HMP_Taxa$Species),] <- "Unidentified_ASV"
HMP_Taxa[10990:10999,]
HMP_Taxa$X <- row.names(HMPtaxa)

#after this step, there should be 58409 ASV left in the data
HMP_Taxa_clean <- subset(HMP_Taxa, X %in% colnames(HMP_clean)[colSums(HMP_clean)!=0])

#same for Ubiome sample
UBM_Taxa <- mutate_all(UBMtaxa, as.character)

UBM_Taxa$Species[is.na(UBM_Taxa$Species)&!is.na(UBM_Taxa$Genus)] <- paste("Unclassified", UBM_Taxa$Genus[is.na(UBM_Taxa$Species)&!is.na(UBM_Taxa$Genus)], sep="_")

UBM_Taxa$Species[is.na(UBM_Taxa$Species)&is.na(UBM_Taxa$Genus)&!is.na(UBM_Taxa$Family)] <- paste("Unclassified", UBM_Taxa$Family[is.na(UBM_Taxa$Species)&is.na(UBM_Taxa$Genus)&!is.na(UBM_Taxa$Family)], sep="_")
UBM_Taxa$Genus[is.na(UBM_Taxa$Genus)&!is.na(UBM_Taxa$Family)] <- paste("Unclassified", UBM_Taxa$Family[is.na(UBM_Taxa$Genus)&!is.na(UBM_Taxa$Family)], sep="_")

UBM_Taxa$Species[is.na(UBM_Taxa$Species)&is.na(UBM_Taxa$Genus)&is.na(UBM_Taxa$Family)&!is.na(UBM_Taxa$Order)] <- paste("Unclassified", UBM_Taxa$Order[is.na(UBM_Taxa$Species)&is.na(UBM_Taxa$Genus)&is.na(UBM_Taxa$Family)&!is.na(UBM_Taxa$Order)], sep="_")
UBM_Taxa$Genus[is.na(UBM_Taxa$Genus)&is.na(UBM_Taxa$Family)&!is.na(UBM_Taxa$Order)] <- paste("Unclassified", UBM_Taxa$Order[is.na(UBM_Taxa$Genus)&is.na(UBM_Taxa$Family)&!is.na(UBM_Taxa$Order)], sep="_")
UBM_Taxa$Family[is.na(UBM_Taxa$Family)&!is.na(UBM_Taxa$Order)] <- paste("Unclassified", UBM_Taxa$Order[is.na(UBM_Taxa$Family)&!is.na(UBM_Taxa$Order)], sep="_")

UBM_Taxa$Species[is.na(UBM_Taxa$Species)&is.na(UBM_Taxa$Genus)&is.na(UBM_Taxa$Family)&is.na(UBM_Taxa$Order)&!is.na(UBM_Taxa$Class)] <- paste("Unclassified", UBM_Taxa$Class[is.na(UBM_Taxa$Species)&is.na(UBM_Taxa$Genus)&is.na(UBM_Taxa$Family)&is.na(UBM_Taxa$Order)&!is.na(UBM_Taxa$Class)], sep="_")
UBM_Taxa$Genus[is.na(UBM_Taxa$Genus)&is.na(UBM_Taxa$Family)&is.na(UBM_Taxa$Order)&!is.na(UBM_Taxa$Class)] <- paste("Unclassified", UBM_Taxa$Class[is.na(UBM_Taxa$Genus)&is.na(UBM_Taxa$Family)&is.na(UBM_Taxa$Order)&!is.na(UBM_Taxa$Class)], sep="_")
UBM_Taxa$Family[is.na(UBM_Taxa$Family)&is.na(UBM_Taxa$Order)&!is.na(UBM_Taxa$Class)] <- paste("Unclassified", UBM_Taxa$Class[is.na(UBM_Taxa$Family)&is.na(UBM_Taxa$Order)&!is.na(UBM_Taxa$Class)], sep="_")
UBM_Taxa$Order[is.na(UBM_Taxa$Order)&!is.na(UBM_Taxa$Class)] <- paste("Unclassified", UBM_Taxa$Class[is.na(UBM_Taxa$Order)&!is.na(UBM_Taxa$Class)], sep="_")

UBM_Taxa$Species[is.na(UBM_Taxa$Species)&is.na(UBM_Taxa$Genus)&is.na(UBM_Taxa$Family)&is.na(UBM_Taxa$Order)&is.na(UBM_Taxa$Class)&!is.na(UBM_Taxa$Phylum)] <- paste("Unclassified", UBM_Taxa$Phylum[is.na(UBM_Taxa$Species)&is.na(UBM_Taxa$Genus)&is.na(UBM_Taxa$Family)&is.na(UBM_Taxa$Order)&is.na(UBM_Taxa$Class)&!is.na(UBM_Taxa$Phylum)], sep="_")
UBM_Taxa$Genus[is.na(UBM_Taxa$Genus)&is.na(UBM_Taxa$Family)&is.na(UBM_Taxa$Order)&is.na(UBM_Taxa$Class)&!is.na(UBM_Taxa$Phylum)] <- paste("Unclassified", UBM_Taxa$Phylum[is.na(UBM_Taxa$Genus)&is.na(UBM_Taxa$Family)&is.na(UBM_Taxa$Order)&is.na(UBM_Taxa$Class)&!is.na(UBM_Taxa$Phylum)], sep="_")
UBM_Taxa$Family[is.na(UBM_Taxa$Family)&is.na(UBM_Taxa$Order)&is.na(UBM_Taxa$Class)&!is.na(UBM_Taxa$Phylum)] <- paste("Unclassified", UBM_Taxa$Phylum[is.na(UBM_Taxa$Family)&is.na(UBM_Taxa$Order)&is.na(UBM_Taxa$Class)&!is.na(UBM_Taxa$Phylum)], sep="_")
UBM_Taxa$Order[is.na(UBM_Taxa$Order)&is.na(UBM_Taxa$Class)&!is.na(UBM_Taxa$Phylum)] <- paste("Unclassified", UBM_Taxa$Phylum[is.na(UBM_Taxa$Order)&is.na(UBM_Taxa$Class)&!is.na(UBM_Taxa$Phylum)], sep="_")
UBM_Taxa$Class[is.na(UBM_Taxa$Class)&!is.na(UBM_Taxa$Phylum)]<- paste("Unclassified", UBM_Taxa$Phylum[is.na(UBM_Taxa$Class)&!is.na(UBM_Taxa$Phylum)], sep="_")

UBM_Taxa$Species[is.na(UBM_Taxa$Species)&is.na(UBM_Taxa$Genus)&is.na(UBM_Taxa$Family)&is.na(UBM_Taxa$Order)&is.na(UBM_Taxa$Class)&is.na(UBM_Taxa$Phylum)&!is.na(UBM_Taxa$Kingdom)] <- paste("Unclassified", UBM_Taxa$Kingdom[is.na(UBM_Taxa$Species)&is.na(UBM_Taxa$Genus)&is.na(UBM_Taxa$Family)&is.na(UBM_Taxa$Order)&is.na(UBM_Taxa$Class)&is.na(UBM_Taxa$Phylum)&!is.na(UBM_Taxa$Kingdom)], sep="_")
UBM_Taxa$Genus[is.na(UBM_Taxa$Genus)&is.na(UBM_Taxa$Family)&is.na(UBM_Taxa$Order)&is.na(UBM_Taxa$Class)&is.na(UBM_Taxa$Phylum)&!is.na(UBM_Taxa$Kingdom)] <- paste("Unclassified", UBM_Taxa$Kingdom[is.na(UBM_Taxa$Genus)&is.na(UBM_Taxa$Family)&is.na(UBM_Taxa$Order)&is.na(UBM_Taxa$Class)&is.na(UBM_Taxa$Phylum)&!is.na(UBM_Taxa$Kingdom)], sep="_")
UBM_Taxa$Family[is.na(UBM_Taxa$Family)&is.na(UBM_Taxa$Order)&is.na(UBM_Taxa$Class)&is.na(UBM_Taxa$Phylum)&!is.na(UBM_Taxa$Kingdom)] <- paste("Unclassified", UBM_Taxa$Kingdom[is.na(UBM_Taxa$Family)&is.na(UBM_Taxa$Order)&is.na(UBM_Taxa$Class)&is.na(UBM_Taxa$Phylum)&!is.na(UBM_Taxa$Kingdom)], sep="_")
UBM_Taxa$Order[is.na(UBM_Taxa$Order)&is.na(UBM_Taxa$Class)&is.na(UBM_Taxa$Phylum)&!is.na(UBM_Taxa$Kingdom)] <- paste("Unclassified", UBM_Taxa$Kingdom[is.na(UBM_Taxa$Order)&is.na(UBM_Taxa$Class)&is.na(UBM_Taxa$Phylum)&!is.na(UBM_Taxa$Kingdom)], sep="_")
UBM_Taxa$Class[is.na(UBM_Taxa$Class)&is.na(UBM_Taxa$Phylum)&!is.na(UBM_Taxa$Kingdom)] <- paste("Unclassified", UBM_Taxa$Kingdom[is.na(UBM_Taxa$Class)&is.na(UBM_Taxa$Phylum)&!is.na(UBM_Taxa$Kingdom)], sep="_")
UBM_Taxa$Phylum[is.na(UBM_Taxa$Phylum)&!is.na(UBM_Taxa$Kingdom)] <- paste("Unclassified", UBM_Taxa$Kingdom[is.na(UBM_Taxa$Phylum)&!is.na(UBM_Taxa$Kingdom)], sep="_")

UBM_Taxa[is.na(UBM_Taxa$Species),]
UBM_Taxa[5270:5279,]
UBM_Taxa[is.na(UBM_Taxa$Species),] <- "Unidentified_ASV"
UBM_Taxa[5270:5279,]

UBM_Taxa_clean <- subset(UBM_Taxa, X %in% colnames(UBM_clean)[colSums(UBM_clean[,-1])!=0]) 

length(table(UBM_Taxa_clean$Genus))

length(table(HMP_Taxa_clean$Genus))
#################################################################################################################################
save(HMP_final, HMP_Taxa_clean, UBM_clean, UBM_Taxa_clean, file = "./Analysis/Robject/OTUandTAX.RData")

write.csv(HMPmeta_clean, "./Analysis/ori_meta_table/Metadata.HMP.csv")


