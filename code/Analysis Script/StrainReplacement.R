#This analysis is intended to look for strain replacement rate in each individuals
library(x)
library(phyloseq)
library(dplyr)
library(tidyverse)
library(lubridate)

setwd("~/Library/CloudStorage/Box-Box/XinZhouFiles/Projects/ZXP2_HMP2_Core/Analysis/")

load("./Robject/DetailedPhyloseq.RData")

######################################################################################################
#Stool replacement
######################################################################################################

#prepare data for prevalence calculation 
Genus_ST <- data.frame(otu_table(physeqGenus_ST))
Taxa_ST <- data.frame(tax_table(physeqGenus_ST))
Sample_ST <- data.frame(sample_data(physeqGenus_ST))
Taxa_ST$ASV <- rownames(Taxa_ST)
subject_list <- unique(Sample_ST$SubjectID)

#replace ASV relative abundance with presents absents measurement
#1 means present, 0 means absent
Genus_ST[Genus_ST>0] <- 1

#replace the ASV name with actual genus name
colnames(Genus_ST) <- Taxa_ST$Genus[match(colnames(Genus_ST),Taxa_ST$ASV)]

#remove genus that only show up once in entire dataset
Genus_ST <- Genus_ST[colSums(Genus_ST) > 1] 
dim(Genus_ST)

#Create list of taxa for loop
st.taxa.list <- colnames(Genus_ST)

Genus_ST_meta <- merge(select(Sample_ST, RandomID, SampleID, SubjectID, Date, IRIS), Genus_ST, by.x="RandomID", by.y="row.names")
Genus_ST_meta[1:15]

#start loop from this line
# Stool_Strain_replace <-data.frame()
# for (taxa_i in st.taxa.list){
#   print(taxa_i)
#   for (subject_j in subject_list){
#     print(subject_j)
#     temp <- filter(Genus_ST_meta, SubjectID == subject_j) %>% select(Date,SampleID, taxa_i)
#     temp <- temp[order(as.Date(temp$Date, format="%Y/%m/%d")),]
#     temp$order <- rank(as.Date(temp$Date, format="%Y/%m/%d"))
#     rle_x = rle(temp[,3])
#     end = cumsum(rle_x$lengths)
#     start = c(1, lag(end)[-1] + 1)
#     result.df <- data.frame(start, end,rle_x$values)
#     
#     result.df$Start_Date <- temp$Date[match(result.df$start, temp$order)]
#     result.df$Start_Date <- temp$Date[match(result.df$start, temp$order)]
#     result.df$End_Date <- temp$Date[match(result.df$end, temp$order)]
#     
#     result.df$Start_Date_before <- temp$Date[match((result.df$start-1), temp$order)]
#     result.df$End_Date_after <- temp$Date[match((result.df$end+1), temp$order)]
#     
#     result.df$Start <- temp$SampleID[match((result.df$start-1), temp$order)]
#     result.df$End <- temp$SampleID[match((result.df$end+1), temp$order)]
#     
#     #result.df_absent <- filter(result.df, rle_x.values == 0)
#     result.df$taxa <- taxa_i
#     result.df$SubjectID <- subject_j
#     result.df <- mutate(result.df, timediff = End_Date_after-Start_Date_before)
#     Stool_Strain_replace <- rbind(Stool_Strain_replace, result.df)
#   }
# }
# 
# Stool_Strain_replace_absent <- filter(Stool_Strain_replace, rle_x.values == 0 & timediff != "NA days")
# Stool_replacement_table <- select(Stool_Strain_replace_absent, Start, End, taxa, SubjectID, timediff)
# Stool_replacement_table$Bodysite <- "Stool"
# 
# write.csv(file = "./tables/strainreplacement/stool.replacement.csv",Stool_replacement_table)
Stool_replacement_table <- read.csv("./tables/strainreplacement/stool.replacement.csv", header = T, row.names = 1)
Stool_replacement_table



######################################################################################################
#Skin replacement
######################################################################################################

#prepare data for prevalence calculation 
Genus_SK <- data.frame(otu_table(physeqGenus_SK))
Taxa_SK <- data.frame(tax_table(physeqGenus_SK))
Sample_SK <- data.frame(sample_data(physeqGenus_SK))
Taxa_SK$ASV <- rownames(Taxa_SK)
subject_list <- unique(Sample_SK$SubjectID)

#replace ASV relative abundance with presents absents measurement
#1 means present, 0 means absent
Genus_SK[Genus_SK>0] <- 1

#replace the ASV name with actual genus name
colnames(Genus_SK) <- Taxa_SK$Genus[match(colnames(Genus_SK),Taxa_SK$ASV)]

#remove genus that only show up once in entire dataset
Genus_SK <- Genus_SK[colSums(Genus_SK) > 1] 
dim(Genus_SK)

#Create list of taxa for loop
sk.taxa.list <- colnames(Genus_SK)

Genus_SK_meta <- merge(select(Sample_SK, KitID, SampleID, SubjectID, Date), Genus_SK, by.x="KitID", by.y="row.names")
Genus_SK_meta[1:15]
# 
# #start loop from this line
# Skin_Strain_replace <-data.frame()
# for (taxa_i in sk.taxa.list){
#   print(taxa_i)
#   for (subject_j in subject_list){
#     print(subject_j)
#     temp <- filter(Genus_SK_meta, SubjectID == subject_j) %>% select(Date,SampleID, taxa_i)
#     temp <- temp[order(as.Date(temp$Date, format="%Y/%m/%d")),]
#     temp$order <- rank(as.Date(temp$Date, format="%Y/%m/%d"))
#     rle_x = rle(temp[,3])
#     end = cumsum(rle_x$lengths)
#     start = c(1, lag(end)[-1] + 1)
#     result.df <- data.frame(start, end,rle_x$values)
#     
#     result.df$Start_Date <- temp$Date[match(result.df$start, temp$order)]
#     result.df$Start_Date <- temp$Date[match(result.df$start, temp$order)]
#     result.df$End_Date <- temp$Date[match(result.df$end, temp$order)]
#     
#     result.df$Start_Date_before <- temp$Date[match((result.df$start-1), temp$order)]
#     result.df$End_Date_after <- temp$Date[match((result.df$end+1), temp$order)]
#     
#     result.df$Start <- temp$SampleID[match((result.df$start-1), temp$order)]
#     result.df$End <- temp$SampleID[match((result.df$end+1), temp$order)]
#     
#     #result.df_absent <- filter(result.df, rle_x.values == 0)
#     result.df$taxa <- taxa_i
#     result.df$SubjectID <- subject_j
#     result.df <- mutate(result.df, timediff = End_Date_after-Start_Date_before)
#     Skin_Strain_replace <- rbind(Skin_Strain_replace, result.df)
#   }
# }
# 
# Skin_Strain_replace_absent <- filter(Skin_Strain_replace, rle_x.values == 0 & timediff != "NA days")
# Skin_replacement_table <- select(Skin_Strain_replace_absent, Start, End, taxa, SubjectID, timediff)
# Skin_replacement_table$Bodysite <- "Skin"

#write.csv(file = "./tables/strainreplacement/skin.replacement.csv",Skin_replacement_table)
Skin_replacement_table <- read.csv("./tables/strainreplacement/skin.replacement.csv", header = T,row.names = 1)

######################################################################################################
#Oral replacement
######################################################################################################
#prepare data for prevalence calculation 
Genus_OR <- data.frame(otu_table(physeqGenus_OR))
Taxa_OR <- data.frame(tax_table(physeqGenus_OR))
Sample_OR <- data.frame(sample_data(physeqGenus_OR))
Taxa_OR$ASV <- rownames(Taxa_OR)
subject_list <- unique(Sample_OR$SubjectID)

#replace ASV relative abundance with presents absents measurement
#1 means present, 0 means absent
Genus_OR[Genus_OR>0] <- 1

#replace the ASV name with actual genus name
colnames(Genus_OR) <- Taxa_OR$Genus[match(colnames(Genus_OR),Taxa_OR$ASV)]

#remove genus that only show up once in entire dataset
Genus_OR <- Genus_OR[colSums(Genus_OR) > 1] 
dim(Genus_OR)

#Create list of taxa for loop
or.taxa.list <- colnames(Genus_OR)

Genus_OR_meta <- merge(select(Sample_OR, KitID, SampleID, SubjectID, Date), Genus_OR, by.x="KitID", by.y="row.names")
Genus_OR_meta[1:15]

# #start loop from this line
# Oral_Strain_replace <-data.frame()
# for (taxa_i in or.taxa.list){
#   print(taxa_i)
#   for (subject_j in subject_list){
#     print(subject_j)
#     temp <- filter(Genus_OR_meta, SubjectID == subject_j) %>% select(Date,SampleID, taxa_i)
#     temp <- temp[order(as.Date(temp$Date, format="%Y/%m/%d")),]
#     temp$order <- rank(as.Date(temp$Date, format="%Y/%m/%d"))
#     rle_x = rle(temp[,3])
#     end = cumsum(rle_x$lengths)
#     start = c(1, lag(end)[-1] + 1)
#     result.df <- data.frame(start, end,rle_x$values)
#     
#     result.df$Start_Date <- temp$Date[match(result.df$start, temp$order)]
#     result.df$Start_Date <- temp$Date[match(result.df$start, temp$order)]
#     result.df$End_Date <- temp$Date[match(result.df$end, temp$order)]
#     
#     result.df$Start_Date_before <- temp$Date[match((result.df$start-1), temp$order)]
#     result.df$End_Date_after <- temp$Date[match((result.df$end+1), temp$order)]
#     
#     result.df$Start <- temp$SampleID[match((result.df$start-1), temp$order)]
#     result.df$End <- temp$SampleID[match((result.df$end+1), temp$order)]
#     
#     #result.df_absent <- filter(result.df, rle_x.values == 0)
#     result.df$taxa <- taxa_i
#     result.df$SubjectID <- subject_j
#     result.df <- mutate(result.df, timediff = End_Date_after-Start_Date_before)
#     Oral_Strain_replace <- rbind(Oral_Strain_replace, result.df)
#   }
# }
# 
# Oral_Strain_replace_absent <- filter(Oral_Strain_replace, rle_x.values == 0 & timediff != "NA days")
# Oral_replacement_table <- select(Oral_Strain_replace_absent, Start, End, taxa, SubjectID, timediff)
# Oral_replacement_table$Bodysite <- "Oral"
# Oral_replacement_table
# 
# write.csv(file = "./tables/strainreplacement/oral.replacement.csv",Oral_replacement_table)
Oral_replacement_table <- read.csv("./tables/strainreplacement/oral.replacement.csv", header = T,row.names = 1)

######################################################################################################
#Nasal replacement
######################################################################################################
#prepare data for prevalence calculation 
Genus_NS <- data.frame(otu_table(physeqGenus_NS))
Taxa_NS <- data.frame(tax_table(physeqGenus_NS))
Sample_NS <- data.frame(sample_data(physeqGenus_NS))
Taxa_NS$ASV <- rownames(Taxa_NS)
subject_list <- unique(Sample_NS$SubjectID)

#replace ASV relative abundance with presents absents measurement
#1 means present, 0 means absent
Genus_NS[Genus_NS>0] <- 1

#replace the ASV name with actual genus name
colnames(Genus_NS) <- Taxa_NS$Genus[match(colnames(Genus_NS),Taxa_NS$ASV)]

#remove genus that only show up once in entire dataset
Genus_NS <- Genus_NS[colSums(Genus_NS) > 1] 
dim(Genus_NS)

#Create list of taxa for loop
ns.taxa.list <- colnames(Genus_NS)

Genus_NS_meta <- merge(select(Sample_NS, RandomID, SampleID, SubjectID, Date), Genus_NS, by.x="RandomID", by.y="row.names")
Genus_NS_meta[1:15]

# #start loop from this line
# Nasal_Strain_replace <-data.frame()
# for (taxa_i in ns.taxa.list){
#   print(taxa_i)
#   for (subject_j in subject_list){
#     print(subject_j)
#     temp <- filter(Genus_NS_meta, SubjectID == subject_j) %>% select(Date,SampleID, taxa_i)
#     temp <- temp[order(as.Date(temp$Date, format="%Y/%m/%d")),]
#     temp$order <- rank(as.Date(temp$Date, format="%Y/%m/%d"))
#     rle_x = rle(temp[,3])
#     end = cumsum(rle_x$lengths)
#     start = c(1, lag(end)[-1] + 1)
#     result.df <- data.frame(start, end,rle_x$values)
#     
#     result.df$Start_Date <- temp$Date[match(result.df$start, temp$order)]
#     result.df$Start_Date <- temp$Date[match(result.df$start, temp$order)]
#     result.df$End_Date <- temp$Date[match(result.df$end, temp$order)]
#     
#     result.df$Start_Date_before <- temp$Date[match((result.df$start-1), temp$order)]
#     result.df$End_Date_after <- temp$Date[match((result.df$end+1), temp$order)]
#     
#     result.df$Start <- temp$SampleID[match((result.df$start-1), temp$order)]
#     result.df$End <- temp$SampleID[match((result.df$end+1), temp$order)]
#     
#     #result.df_absent <- filter(result.df, rle_x.values == 0)
#     result.df$taxa <- taxa_i
#     result.df$SubjectID <- subject_j
#     result.df <- mutate(result.df, timediff = End_Date_after-Start_Date_before)
#     Nasal_Strain_replace <- rbind(Nasal_Strain_replace, result.df)
#   }
# }
# 
# Nasal_Strain_replace_absent <- filter(Nasal_Strain_replace, rle_x.values == 0 & timediff != "NA days")
# Nasal_replacement_table <- select(Nasal_Strain_replace_absent, Start, End, taxa, SubjectID, timediff)
# Nasal_replacement_table$Bodysite <- "Nasal"

#write.csv(file = "./tables/strainreplacement/nasal.replacement.csv",Nasal_replacement_table)
Nasal_replacement_table <- read.csv("./tables/strainreplacement/nasal.replacement.csv", header = T,row.names = 1)
Nasal_replacement_table


#############appendix: try it for a single bacteria and single individual
taxa_i <- "Staphylococcus"
subject_j <- "69-001"

temp <- filter(Genus_ST_meta, SubjectID == subject_j) %>% select(Date,SampleID, taxa_i)
temp <- temp[order(as.Date(temp$Date, format="%Y/%m/%d")),]
temp$order <- rank(as.Date(temp$Date, format="%Y/%m/%d"))
rle_x = rle(temp[,3])
end = cumsum(rle_x$lengths)
start = c(1, lag(end)[-1] + 1)
result.df <- data.frame(start, end,rle_x$values)

result.df$Start_Date <- temp$Date[match(result.df$start, temp$order)]
result.df$Start_Date <- temp$Date[match(result.df$start, temp$order)]
result.df$End_Date <- temp$Date[match(result.df$end, temp$order)]

result.df$Start_Date_before <- temp$Date[match((result.df$start-1), temp$order)]
result.df$End_Date_after <- temp$Date[match((result.df$end+1), temp$order)]

result.df$Start <- temp$SampleID[match((result.df$start-1), temp$order)]
result.df$End <- temp$SampleID[match((result.df$end+1), temp$order)]

#result.df_absent <- filter(result.df, rle_x.values == 0)
result.df$taxa <- taxa_i
result.df$SubjectID <- subject_j
result.df <- mutate(result.df, timediff = End_Date_after-Start_Date_before)
Stool_Strain_replace <- rbind(Stool_Strain_replace, result.df)
Stool_Strain_replace
