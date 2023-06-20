# 00_prepare-dataset
# 
# Dan Spakowicz (& modified by JJ)
# 
############################################################################################################
# JJ Create a single dataset containing either IL17A, IL17F, or IL22, plus metadata 


library(tidyverse)


## Read in the data --------------------------------------------------------

df.tax <- read.table("../data/stool_16s_data.tsv",
                     header = TRUE, sep = "\t")

df.cyt <- read.table('../data/stool_cytokine_data.tsv',
                     header=TRUE, sep='\t')

df.meta <- read.table('../data/stool_meta_data.tsv',
                      header=TRUE, sep='\t', stringsAsFactors=FALSE)


df.clus <- read.csv("data/threegroup0425_3cytokine_SD_Scale.csv") %>%
  dplyr::select(Row.names, V1, V2, V3) %>%
  mutate(cluster = apply(., 1, function(x) which.max(x) - 1)) %>%
  dplyr::select(Row.names, cluster) %>%
  rename("SubjectID" = "Row.names") %>%
  mutate(cluster = as.factor(cluster))


###############################################################################
# Fix the metadata so that there is a healthy vs. sick column
unique(df.meta$InfectionState_jax)
df.meta$Infection <- df.meta$InfectionState_jax
df.meta$Infection[df.meta$InfectionState_jax %in% c("early", "late", "recovery", "unknown_infection")] <- 'sick'
df.meta$Infection[df.meta$InfectionState_jax %in% c("healthy_other", "pre", "other", "post")] <- 'healthy'
unique(df.meta$Infection)

###############################################################################
# Check that the data are in the same order in each of the three tables
identical(sort(row.names(df.cyt)), sort(row.names(df.meta)))
identical(sort(row.names(df.cyt)), sort(row.names(df.tax)))

df.meta <- df.meta[row.names(df.cyt),]
df.tax <- df.tax[row.names(df.cyt),]

identical(row.names(df.meta), row.names(df.tax))

###############################################################################
# Save the data.frame of taxa 
###############################################################################
save(df.tax, file = "microbes.Rdata")

##########################################################################################
# Fetch the relevant data from the relevant tables for IL17F
df.tmp <- cbind(df.meta[, c("SampleID", 'SubjectID', 'Season', 'Infection', "CollectionDate")],
                df.cyt[, c('IL17F')])
names(df.tmp) <- c("SampleID", 'SubjectID', 'Season', 'Infection', "CollectionDate", 'IL17F')

# Set Factors
df.tmp$Season <- as.factor(df.tmp$Season)
df.tmp$Infection <- as.factor(df.tmp$Infection)
head(df.tmp)

# Format dates
df.tmp$CollectionDate <- as.Date(df.tmp$CollectionDate, format = "%F")

# Create normalized days variables
model.df <- 
  df.tmp %>%
  mutate(CollectionDate = as.Date(CollectionDate, format = "%F")) %>%
  group_by(SubjectID) %>%
  mutate(days = CollectionDate - min(CollectionDate)) %>%
  mutate(days2 = as.numeric(days)^2) %>%
  mutate(days3 = as.numeric(days)^3) %>%
  left_join(., df.clus) %>%
  ungroup

model.df$SubjectID <- as.factor(model.df$SubjectID)

###
### Save output for modeling
###

save(model.df, file = "IL17F.Rdata")

##########################################################################################
# Fetch the relevant data from the relevant tables for IL17A
df.tmp <- cbind(df.meta[, c("SampleID", 'SubjectID', 'Season', 'Infection', "CollectionDate")],
                  df.cyt[, c('IL17A')])
names(df.tmp) <- c("SampleID", 'SubjectID', 'Season', 'Infection', "CollectionDate", 'IL17A')

# Set Factors
df.tmp$Season <- as.factor(df.tmp$Season)
df.tmp$Infection <- as.factor(df.tmp$Infection)
head(df.tmp)

# Format dates
df.tmp$CollectionDate <- as.Date(df.tmp$CollectionDate, format = "%F")

# Create normalized days variables
model.df <- 
  df.tmp %>%
  mutate(CollectionDate = as.Date(CollectionDate, format = "%F")) %>%
  group_by(SubjectID) %>%
  mutate(days = CollectionDate - min(CollectionDate)) %>%
  mutate(days2 = as.numeric(days)^2) %>%
  mutate(days3 = as.numeric(days)^3) %>%
  left_join(., df.clus) %>%
  ungroup

model.df$SubjectID <- as.factor(model.df$SubjectID)

###
### Save output for modeling
###

save(model.df, file = "IL17A.Rdata")

##########################################################################################
# Fetch the relevant data from the relevant tables for IL22
df.tmp <- cbind(df.meta[, c("SampleID", 'SubjectID', 'Season', 'Infection', "CollectionDate")],
                  df.cyt[, c('IL22')])
names(df.tmp) <- c("SampleID", 'SubjectID', 'Season', 'Infection', "CollectionDate", 'IL22')

# Set Factors
df.tmp$Season <- as.factor(df.tmp$Season)
df.tmp$Infection <- as.factor(df.tmp$Infection)
head(df.tmp)

# Format dates
df.tmp$CollectionDate <- as.Date(df.tmp$CollectionDate, format = "%F")

# Create normalized days variables
model.df <- 
  df.tmp %>%
  mutate(CollectionDate = as.Date(CollectionDate, format = "%F")) %>%
  group_by(SubjectID) %>%
  mutate(days = CollectionDate - min(CollectionDate)) %>%
  mutate(days2 = as.numeric(days)^2) %>%
  mutate(days3 = as.numeric(days)^3) %>%
  left_join(., df.clus) %>%
  ungroup

model.df$SubjectID <- as.factor(model.df$SubjectID)

###
### Save output for modeling
###

save(model.df, file = "IL22.Rdata")

