#### Remove all variables from the workspace ####
rm(list = ls())

#### Set directory ####
setwd('/Users/hmallick/Dropbox (Personal)/Repos/hmp2_analysis')

# Load required utility functions
source("./differential_abundance/src/core_DA_functions.R")

library(plyr)
library(dplyr)
library(nlme)
library(tibble)

#########
# HMP 2 #
#########

source("./common/disease_colors.r")
source("./common/load_bugs.r")
source("./common/nonredundant_samples.r")

# Remove tech reps in HMP2
hmp2_bugs.pcl.unfiltered<-nonredundant_samples(bugs.pcl)
# bugs.pcl$ns - hmp2_bugs.pcl.unfiltered$ns  # 42 removed
# hmp2_bugs.pcl.unfiltered$ns; hmp2_bugs.pcl.unfiltered$nf # 1578 samples, 1479 features

# Remove Pediatrics Samples
hmp2_bugs.pcl0 <- hmp2_bugs.pcl.unfiltered %>% pcl.filter.s(keep=which(hmp2_bugs.pcl.unfiltered$meta$consent_age >= 18))
# hmp2_bugs.pcl0$ns; hmp2_bugs.pcl0$nf # 823 samples, 1479 features

# Remove Intial Samples
hmp2_bugs.pcl <- hmp2_bugs.pcl0 %>% pcl.filter.s(keep=which(hmp2_bugs.pcl0$meta$week_num > 20))
# hmp2_bugs.pcl$ns; hmp2_bugs.pcl$nf # 410 samples, 1479 features


############
# HMP 1-II #
############

hmp12_bugs_unfiltered.pcl <- pcl.read(file.path(HMP2_data, "hmp1-ii/hmp1-II_metaphlan2-mtd-qcd.tsv"), metadata.rows=8) %>%
  pcl.filter.s(STSite == "Stool")

# Remove tech reps
dt<-hmp12_bugs_unfiltered.pcl$meta
dt_notechrep<-dt[!duplicated(dt[1:2]),]
# nrow(dt) - nrow(dt_notechrep) # About 82 removed

# Subset and match
hmp12_bugs.pcl<-list()
hmp12_bugs.pcl$x<-hmp12_bugs_unfiltered.pcl$x[rownames(dt_notechrep),]
hmp12_bugs.pcl$meta<-hmp12_bugs_unfiltered.pcl$meta[rownames(dt_notechrep),]
all(rownames(hmp12_bugs.pcl$x)==rownames(hmp12_bugs.pcl$meta))
hmp12_bugs.pcl$ns<-nrow(hmp12_bugs.pcl$x)
hmp12_bugs.pcl$nf<-ncol(hmp12_bugs.pcl$x)

# Merge the datasets
hmp12_bugs.pcl$meta$diagnosis <- "HMP1-II"
hmp12_bugs.pcl$meta$VISNO[hmp12_bugs.pcl$meta$VISNO=="02S"] <- "2"
hmp12_bugs.pcl$meta$VISNO[hmp12_bugs.pcl$meta$VISNO=="03S"] <- "3"
hmp12_bugs.pcl$meta$VISNO <- as.numeric(as.character(hmp12_bugs.pcl$meta$VISNO))
merged.pcl <- pcl.merge(hmp2_bugs.pcl, hmp12_bugs.pcl)

merged.pcl$meta$merged_subj <- as.character(merged.pcl$meta$subject)
merged.pcl$meta$merged_subj[is.na(merged.pcl$meta$merged_subj)] <-
  as.character(merged.pcl$meta$RANDSID)[is.na(merged.pcl$meta$merged_subj)]
merged.pcl$meta$merged_subj <- factor(merged.pcl$meta$merged_subj)

# Restrict to nonIBD Only
merged.pcl.healthy <- merged.pcl %>%
    pcl.filter.s(diagnosis == "nonIBD" || diagnosis == "HMP1-II") %>%
    pcl.only(rank="s") %>% pcl.nicenames %>%
    pcl.normalize

################################################################
# Species-level DA to distinguish HMP1-II from non-IBD in HMP2 #
################################################################

# Filter out features with no variance or with >90% zeros
species.pcl <- merged.pcl.healthy %>% pcl.filter.f(any(is.finite(x)) && var(x, na.rm=T) > 0 && mean(!is.na(x) & (x>0)) >= 0.1)
species.pcl$x<- species.pcl$x[ , colSums(is.na(species.pcl$x)) == 0] # Delete columns with NA's
species.pcl$nf<-ncol(species.pcl$x)
species.pcl$ns<-nrow(species.pcl$x)

# Filter out features with little variance
sds <- pcl.apply.f(species.pcl, sd(x, na.rm=T))
species.pcl <- species.pcl %>% pcl.filter.f(keep=sds > median(sds) / 2)
species.pcl <- pcl.map.fnames(species.pcl, gsub("^.*?([^\\|]+)$", "\\1", Name))
colnames(species.pcl$x)<-gsub(".*s__", "", colnames(species.pcl$x))

# Run Analysis
run_DA_HMP1_II(species.pcl, transformation = 'CLR_perfeature')
run_DA_HMP1_II(species.pcl, transformation = 'CLR_persample')
run_DA_HMP1_II(species.pcl, transformation = 'CLR_overall')
run_DA_HMP1_II(species.pcl, transformation = 'Vanilla')


