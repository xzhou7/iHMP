
#Load Library

library(tidyr)
library(dplyr)
library(tidyverse)

#------------------------------------------------------------------------#
#pool uBiome File together
#xin zhou
#Aug. 02. 2019
#Updated Sep. 04. 2019
#------------------------------------------------------------------------#

tax01 <- read.csv("~/Desktop/01_uBiome_DataSharing_Sept-2017 /01_taxonomy_annotation.tsv", sep="\t")
tax02 <- read.csv("~/Desktop/02_uBiome_DataSharing-4-Jan-2018/02_taxonomy_annotation.tsv", sep="\t")
tax03 <- read.csv("~/Desktop/03_uBiome_DataSharing-27-Mar-2018/03_taxonomy_annotation.tsv", sep="\t")
tax04 <- read.csv("~/Desktop/04_uBiome-DataSharing-18-07-2018/04_taxonomy_annotation.tsv", sep="\t")
tax05 <- read.csv("~/Desktop/05_uBiome-DataSharing-29-Nov-2018/05_taxonomy_annotation.tsv", sep="\t")

#Add barcode to Tax01
tax01$barcodeexp <- gsub("(_[^_]+)$","", tax01$Sample)
tax01$barcode <- gsub("(_[^_]_[^_]+)$","", tax01$Sample)
tax01$exp <- gsub("([^_]+_)","", tax01$barcodeexp)
tax01$site <- gsub("([^_]+_)","", tax01$Sample)
table(tax01$site)


#remove two duplicated samples
#354107326_NA0008329011_4_mouth in tax03&tax01 (only keep it in Tax04)
#277094135_NA0009416083_4_mouth in tax02 (Only keep it in Tax05)
#remove any sample that is started with 756351017 (duplicates without metadata)
tax01.1 <- subset(tax01, Sample != "354107326_4_mouth")
tax02.1 <- subset(tax02, Sample != "277094135_NA0009416083_4_mouth")
tax02.2 <- subset(tax02.1, !grepl("756351017*", tax02.1$Sample))
tax03.1 <- subset(tax03, Sample != "354107326_NA0008329011_4_mouth")



#add batch and barcode to Tax02~Tax05 
#tax01.1$batch <- "Batch01"
#tax02.2$batch <- "Batch02"
#tax03.1$batch <- "Batch03"
#tax04$batch <- "Batch04"
#tax05$batch <- "Batch05"

tax025 <- rbind(tax02.2,tax03.1,tax04,tax05)
tax025$Sample <- gsub("_f$", "", tax025$Sample)
tax025$barcodeexp <- gsub("(_[^_]+)$","", tax025$Sample)
tax025$barcode <- gsub("(_[^_]+_[^_]+_[^_]+)$","", tax025$Sample)
tax025$exp <- gsub("([^_]+_)","", tax025$barcodeexp)
tax025$site <- gsub("([^_]+_)","", tax025$Sample)

taxfile <- rbind(tax01.1, tax025)
#rm(tax01,tax02,tax025, tax03, tax04, tax05)
taxfile$KitID <- paste("UBMkid", taxfile$barcode, sep = "")

table(taxfile$site) 

root.table <- subset(taxfile, Rank=="root")

#write a batching list
#batchingfile <- unique(select(taxfile, KitID, batch))
#write.csv(batchingfile, "~/Desktop/QC result/batchinglist.csv")

#------------------------------------------------------------------------#
#Ubiome metadata cleaning
#Xin Zhou
#Aug. 07, 2019
#------------------------------------------------------------------------#

MIX <- read.csv("~/Desktop/06uBiomeMeta/XZ_2016_Mix_WZ.csv", header = T)
Skin2017 <- read.csv("~/Desktop/06uBiomeMeta/XZ_2017_skin_WZ.csv", header = T)
Skin2018 <- read.csv("~/Desktop/06uBiomeMeta/XZ_2018_Skin_WZ.csv", header = T)
Tongue2017  <- read.csv("~/Desktop/06uBiomeMeta/XZ_2017_tongue_WZ.csv", header = T)
Tongue2018  <- read.csv("~/Desktop/06uBiomeMeta/XZ_2018_Tongue_WZ.csv", header = T)

Metafile_01 <- rbind(MIX, Skin2017, Skin2018, Tongue2017, Tongue2018)
Metafile_02 <- read.csv("~/Desktop/06uBiomeMeta/XZ_KitID_201807.csv", header = T)

#merge two metafiles
metafile <- merge(unique(Metafile_01), unique(Metafile_02), by.x = "KitID", by.y = "Kit.id")

#Find Sequencing File that are not part of iPOP study
nometa <- setdiff(root.table$KitID, metafile$KitID)

nometaroot <- root.table[which(root.table$KitID %in% nometa),]
nometaskin <- nometaroot[nometaroot$site=='skin',]

#get Kit ID for sample counts less than 1000
kitless1000 <- root.table[root.table$Counts<1000,]$KitID
taxfile1000 <- taxfile[-which(taxfile$KitID %in% kitless1000),]

#------------------------------------------------------------------------#
#get tax table by rank
#------------------------------------------------------------------------#
root.table <- subset(taxfile1000, Taxon=="root")
table(root.table$site)
#identify duplicates (should return <0 rows>)
root.table[duplicated(gsub("_[^_]{12}_","_",root.table$barcodeexp)),]
#identify duplicated barcode (which is not an issue, since these are not ipop data)
table(root.table$barcode)[table(root.table$barcode)>1]
root.table.ipop <- root.table[-which(root.table$KitID %in% nometa),]

#get phylum table
phylumall<- subset(taxfile1000, Rank==("phylum"))
phylum.with.meta <- phylumall[-which(phylumall$KitID %in% nometa),]
phylum.ipop <- phylum.with.meta[,c("KitID", "Taxon", "Counts")] %>% spread(key="Taxon", value="Counts")
phylum.ipop.freq <- phylum.with.meta[,c("KitID", "Taxon", "Percent")] %>% spread(key="Taxon", value="Percent")
#replace all NAs to 0
phylum.ipop[is.na(phylum.ipop)] <- 0
phylum.ipop.freq[is.na(phylum.ipop.freq)] <- 0
table(is.na(phylum.ipop))
table(is.na(phylum.ipop.freq))
rowSums(phylum.ipop.freq[,-1])
phylum.ipop.freq$Unclassified_phylum <- 100 -  (rowSums(phylum.ipop.freq[,-1]))
#check if all rows are sum to 100%
summary((rowSums(phylum.ipop.freq[,-1]))==100)

#find duplication in phylum data (no need, because no duplicated barcode in ipop)
#phylum.ipop[duplicated(gsub("_[^_]{12}_","_",phylum.ipop$barcodeexp)),]
#Did not find metadata for 756351017, thus remove this file
#get duplicate data for 354107326*, check if two samples are similar
#subset(phylumcast, grepl("354107326*", phylumcast$barcodeexp))

phylum.finaltable <- phylum.ipop.freq
phylum.finaltable$Unclassified_phylum <- 100 - rowSums(phylum.ipop.freq[,-1])
phylum.finalcount <- phylum.finaltable[,-1] * (rowSums(phylum.ipop[,-1]))
row.names(phylum.finalcount) <- phylum.finaltable[,1]

write.csv(phylum.finaltable, file="~/Desktop/Ubiome result/Tax Table/Ubiome_phylum_freq.csv")
write.csv(phylum.finalcount, file="~/Desktop/Ubiome result/Tax Table/Ubiome_phylum_count.csv")

#class
classall <- subset(taxfile1000, Rank==("class"))
class.with.meta <- classall[-which(classall$KitID %in% nometa),]
class.ipop <- class.with.meta[,c("KitID", "Taxon", "Counts")] %>% spread(key="Taxon", value="Counts")
class.ipop.freq <- class.with.meta[,c("KitID", "Taxon", "Percent")] %>% spread(key="Taxon", value="Percent")
class.ipop[is.na(class.ipop)] <- 0
class.ipop.freq[is.na(class.ipop.freq)] <- 0
rowSums(class.ipop.freq[,-1]) %>% 
  density() %>% 
  plot(main="Class Level", xlab = "Percent of Taxonomy Identified in 2118 Samples") -> classplot
classplot

class.finaltable <- class.ipop.freq
class.finaltable$Unclassified_class <- 100 - rowSums(class.ipop.freq[,-1])
class.finalcount <- class.finaltable[,-1] * (rowSums(phylum.ipop[,-1]))
row.names(class.finalcount) <- class.finaltable[,1]

write.csv(class.finaltable, file="~/Desktop/Ubiome result/Tax Table/Ubiome_class_freq.csv")
write.csv(class.finalcount, file="~/Desktop/Ubiome result/Tax Table/Ubiome_class_count.csv")


#order
orderall <- subset(taxfile1000, Rank==("order"))
order.with.meta <- orderall[-which(orderall$KitID %in% nometa),]
order.ipop <- order.with.meta[,c("KitID", "Taxon", "Counts")] %>% spread(key="Taxon", value="Counts")
order.ipop.freq <- order.with.meta[,c("KitID", "Taxon", "Percent")] %>% spread(key="Taxon", value="Percent")
order.ipop[is.na(order.ipop)] <- 0
order.ipop.freq[is.na(order.ipop.freq)] <- 0
rowSums(order.ipop.freq[,-1]) %>% density() -> orderplot
plot(orderplot, main="Order Level", xlab = "Percent of Taxonomy Identified in 2118 Samples")

order.finaltable <- order.ipop.freq
order.finaltable$Unclassified_order <- 100 - rowSums(order.ipop.freq[,-1])
order.finalcount <- order.finaltable[,-1] * (rowSums(phylum.ipop[,-1]))
row.names(order.finalcount) <- order.finaltable[,1]

write.csv(order.finaltable, file="~/Desktop/Ubiome result/Tax Table/Ubiome_order_freq.csv")
write.csv(order.finalcount, file="~/Desktop/Ubiome result/Tax Table/Ubiome_order_count.csv")


#family
familyall <- subset(taxfile1000, Rank==("family"))
family.with.meta <- familyall[-which(familyall$KitID %in% nometa),]
family.ipop <- family.with.meta[,c("KitID", "Taxon", "Counts")] %>% spread(key="Taxon", value="Counts")
family.ipop.freq <- family.with.meta[,c("KitID", "Taxon", "Percent")] %>% spread(key="Taxon", value="Percent")
family.ipop[is.na(family.ipop)] <- 0
family.ipop.freq[is.na(family.ipop.freq)] <- 0
rowSums(family.ipop.freq[,-1]) %>% density() -> familyplot
plot(familyplot, main="Family Level", xlab = "Percent of Taxonomy Identified in 2118 Samples")

family.finaltable <- family.ipop.freq
family.finaltable$Unclassified_family <- 100 - rowSums(family.ipop.freq[,-1])
family.finalcount <- family.finaltable[,-1] * (rowSums(phylum.ipop[,-1]))/100
row.names(family.finalcount) <- family.finaltable[,1]
family.finalcount$Unclassified_family[family.finalcount$Unclassified_family < 1] <- 0

write.csv(family.finaltable, file="~/Desktop/Ubiome result/Tax Table/Ubiome_family_freq.csv")
write.csv(family.finalcount, file="~/Desktop/Ubiome result/Tax Table/Ubiome_family_count.csv")

#genus
genusall <- subset(taxfile1000, Rank==("genus"))
genus.with.meta <- genusall[-which(genusall$KitID %in% nometa),]
genus.ipop <- genus.with.meta[,c("KitID", "Taxon", "Counts")] %>% spread(key="Taxon", value="Counts")
genus.ipop.freq <- genus.with.meta[,c("KitID", "Taxon", "Percent")] %>% spread(key="Taxon", value="Percent")
genus.ipop[is.na(genus.ipop)] <- 0
genus.ipop.freq[is.na(genus.ipop.freq)] <- 0
rowSums(genus.ipop.freq[,-1]) %>% density() -> genusplot
plot(genusplot, main="Genus Level", xlab = "Percent of Taxonomy Identified in 2118 Samples")
table(colSums(genus.ipop.freq[,-1])==0)

genus.finaltable <- genus.ipop.freq
genus.finaltable$Unclassified_genus <- 100 - (rowSums(genus.ipop.freq[,-1]))
genus.finalcount <- genus.finaltable[,-1] * (rowSums(phylum.ipop[,-1]))/100
row.names(genus.finalcount) <- genus.finaltable[,1]
genus.finalcount$Unclassified_genus[genus.finalcount$Unclassified_genus < 1] <- 0

rowSums(genus.ipop[,-1])
rowSums(genus.finalcount)

write.csv(genus.finaltable, file="~/Desktop/Ubiome result/Tax Table/Ubiome_genus_freq.csv")
write.csv(genus.finalcount, file="~/Desktop/Ubiome result/Tax Table/Ubiome_genus_count.csv")

#species
speciesall <- subset(taxfile1000, Rank==("species"))
species.with.meta <- speciesall[-which(speciesall$KitID %in% nometa),]
species.ipop <- species.with.meta[,c("KitID", "Taxon", "Counts")] %>% spread(key="Taxon", value="Counts")
species.ipop.freq <- species.with.meta[,c("KitID", "Taxon", "Percent")] %>% spread(key="Taxon", value="Percent")
species.ipop[is.na(species.ipop)] <- 0
species.ipop.freq[is.na(species.ipop.freq)] <- 0

rowSums(species.ipop.freq[,-1]) %>% density() -> speciesplot
plot(speciesplot, main="Species Level", xlab = "Percent of Taxonomy Identified in 2118 Samples")

species.finaltable <- species.ipop.freq
species.finaltable$Unclassified_Species <- 100 - rowSums(species.finaltable[,-1])
species.finalcount <- species.finaltable[,-1] * (rowSums(phylum.ipop[,-1]))
row.names(species.finalcount) <- species.finaltable[,1]

write.csv(species.finaltable, file="~/Desktop/Ubiome result/Tax Table/Ubiome_Species_freq.csv")
write.csv(species.finalcount, file="~/Desktop/Ubiome result/Tax Table/Ubiome_Species_count.csv")
#see if genes table and phylum table contains same number of subjects
setdiff(phylum.with.meta$KitID,genus.ipop$KitID)

#prevalence abundace ratio
plot(log(((colSums(genus.finalcount!=0)))/2118),log(apply(genus.finaltable[,-1],2,function(x) mean(x[x!=0]))), xlab = "log_prevalance", ylab="log_mean_abundance")


