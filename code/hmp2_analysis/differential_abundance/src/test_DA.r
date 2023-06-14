#### Remove all variables from the workspace ####
rm(list = ls())

#### Set directory ####
setwd('/Users/hmallick/Dropbox (Personal)/Repos/hmp2_analysis')

###### Required package list ######
library(plyr); library(dplyr); library(vegan); library(car);
library(nlme); library(lme4); library(data.table);
library(tibble); library(stringr); library(magrittr)
# library(devtools); #install_github("zdk123/SpiecEasi")
library(SpiecEasi)

# Load specific data loaders (takes a while)
source("./common/load_bugs.r")
source("./common/load_metabolites.r")
source("./common/load_proteins.r")
source("./common/load_pathways.r")
source("./common/load_ecs.r")
source("./common/load_kos.r")
source("./common/merge_metadata.r")
source("./common/match_datasets.r")
source("./common/disease_activity.r") # By default merges disease activity info in the bugs table.
# To merge the disease activity info to other files, use the function merge_disease_activity(pcl).

# Load required utility functions
source("./differential_abundance/src/core_DA_functions.R")

################
# Prepare data #
################

# Subset to species-level taxa only
species.pcl <- pcl.only(bugs.pcl, rank="s")
species.pcl.vst <- pcl.map.fnames(species.pcl, gsub("^.*?([^\\|]+)$", "\\1", Name))
colnames(species.pcl.vst$x)<-gsub(".*s__", "", colnames(species.pcl.vst$x))

######################################################################
# Dump CLR-transformed Species Profiles for Ordination and PERMANOVA #
######################################################################

# Filter Low-Prevalence Bugs
species_CLR.pcl <- species.pcl.vst %>% pcl.filter.f(any(is.finite(x)) && var(x, na.rm=T) > 0 && mean(!is.na(x) & (x>0)) >= 0.1)

# CLR Transformation
species_CLR_perfeature.pcl<-  species_CLR.pcl %>% CLR_perfeature
species_CLR_persample.pcl<-  species_CLR.pcl %>% CLR_persample
species_CLR_overall.pcl<- species_CLR.pcl %>% CLR_overall

# Dump
save(species_CLR_perfeature.pcl, file = paste(HMP2_root, "/data/mgx/species_CLR_perfeature.RData", sep=''))
save(species_CLR_persample.pcl, file = paste(HMP2_root, "/data/mgx/species_CLR_persample.RData", sep=''))
save(species_CLR_overall.pcl, file = paste(HMP2_root, "/data/mgx/species_CLR_overall.RData", sep=''))

# Binarize and filter out all-zero or near-zero viral features
viruses.pcl$x<-ifelse(viruses.pcl$x > 0, 1, 0)
viruses.pcl.vst <- viruses.pcl %>% pcl.filter.f(keep=which(colMeans(viruses.pcl$x) > 0.05))
viruses.pcl.vst <- pcl.map.fnames(viruses.pcl.vst, gsub("^.*?([^\\|]+)$", "\\1", Name))
colnames(viruses.pcl.vst$x)<-gsub(".*s__", "", colnames(viruses.pcl.vst$x))

# Unique metabolites with names
metabolites.pcl.vst<-metabolites.named.pcl

# Full metabolites for Sena
metabolites_Sena.pcl.vst<-metabolites.pcl


# Filter KO DNAs and RNAs (flag KOs highly correlated with a single bug abundance, threshold 0.6)
ko.dna.unstrat.pcl.filtered<-filterKO(ko.dna.unstrat.pcl, species.pcl)
ko.rna.unstrat.pcl.filtered<-filterKO(ko.rna.unstrat.pcl, species.pcl)

# Create PCL for RNA/DNA ratios from KO bug roll ups
bugkoRNADNAs.pcl<-createPCL(bugs.fromko.rna.pcl, bugs.fromko.dna.pcl)

#############################################
# Define all omic names and transformations #
#############################################
transMethods<-c('CLR_perfeature', 'CLR_persample', 'CLR_overall', 'Vanilla')

# Rename rest of the pcl files similarly for consistency
proteins.pcl.vst<-proteins.pcl
proteins.ecs.pcl.vst<-proteins.ecs.pcl
proteins.kos.pcl.vst<-proteins.kos.pcl
pwy.dna.unstrat.pcl.vst<-pwy.dna.unstrat.pcl
pwy.rna.unstrat.pcl.vst<-pwy.rna.unstrat.pcl
pwy.rnadna.unstrat.pcl.vst<-pwy.rnadna.unstrat.pcl
pwy.rnadna.unstrat.pcl.vst$x[!is.finite(pwy.rnadna.unstrat.pcl.vst$x)] <- NA
ec.dna.unstrat.pcl.vst<-ec.dna.unstrat.pcl
ec.rna.unstrat.pcl.vst<-ec.rna.unstrat.pcl
ko.dna.unstrat.pcl.filtered.vst<-ko.dna.unstrat.pcl.filtered
ko.rna.unstrat.pcl.filtered.vst<-ko.rna.unstrat.pcl.filtered
bugkoRNADNAs.pcl.vst<-bugkoRNADNAs.pcl
bugkoRNADNAs.pcl.vst$x[!is.finite(bugkoRNADNAs.pcl.vst$x)] <- NA

################################
# Analysis 1: Perform DA tests #
################################

# Run for every transformation and omics feature combination
for (i in 1:length(transMethods)){
    run_DA(species.pcl.vst, "bugs", transformation = transMethods[i])
    run_DA(bugkoRNADNAs.pcl.vst, "bugkoRNADNAs", transMethods[i])
    run_DA(proteins.pcl.vst, "proteins", transMethods[i])
    run_DA(proteins.ecs.pcl.vst, "proteinECs", transMethods[i])
    run_DA(proteins.kos.pcl.vst, "proteinKOs", transMethods[i])
    run_DA(pwy.dna.unstrat.pcl.vst, "pwyDNAs", transMethods[i])
    run_DA(pwy.rna.unstrat.pcl.vst, "pwyRNAs", transMethods[i])
    run_DA(pwy.rnadna.unstrat.pcl.vst, "pwyRNADNAs", transMethods[i])
    run_DA(ec.dna.unstrat.pcl.vst, "ecDNAs", transMethods[i])
    run_DA(ec.rna.unstrat.pcl.vst, "ecRNAs", transMethods[i])
    run_DA(ko.dna.unstrat.pcl.filtered.vst, "koDNAs", transMethods[i])
    run_DA(ko.rna.unstrat.pcl.filtered.vst, "koRNAs", transMethods[i])
    run_DA(metabolites.pcl.vst, "metabolites", transMethods[i])
    library(lme4); run_DA_binomial(viruses.pcl.vst, "viruses",  transMethods[i])
    }

#################################################
# Analysis 2: Perform DA tests (Activity Index) #
#################################################

# Just run for Vanilla
run_DA_activity_index(species.pcl.vst, "bugs",  transformation = 'Vanilla')
run_DA_activity_index(bugkoRNADNAs.pcl.vst, "bugkoRNADNAs",  transformation = 'Vanilla')
run_DA_activity_index(proteins.pcl.vst, "proteins",  transformation = 'Vanilla')
run_DA_activity_index(proteins.ecs.pcl.vst, "proteinECs",  transformation = 'Vanilla')
run_DA_activity_index(proteins.kos.pcl.vst, "proteinKOs",  transformation = 'Vanilla')
run_DA_activity_index(pwy.dna.unstrat.pcl.vst, "pwyDNAs",  transformation = 'Vanilla')
run_DA_activity_index(pwy.rna.unstrat.pcl.vst, "pwyRNAs",  transformation = 'Vanilla')
run_DA_activity_index(pwy.rnadna.unstrat.pcl.vst, "pwyRNADNAs",  transformation = 'Vanilla')
run_DA_activity_index(ec.dna.unstrat.pcl.vst, "ecDNAs",  transformation = 'Vanilla')
run_DA_activity_index(ec.rna.unstrat.pcl.vst, "ecRNAs",  transformation = 'Vanilla')
run_DA_activity_index(ko.dna.unstrat.pcl.filtered.vst, "koDNAs",  transformation = 'Vanilla')
run_DA_activity_index(ko.rna.unstrat.pcl.filtered.vst, "koRNAs",  transformation = 'Vanilla')
run_DA_activity_index(metabolites.pcl.vst, "metabolites",  transformation = 'Vanilla')
library(lme4); run_DA_activity_index_binomial(viruses.pcl.vst, "viruses",  transformation = 'Vanilla')

###########################################
# Analysis 3: Perform DA tests (baseline) #
###########################################

# Just run for Vanilla

run_baseline_DA(species.pcl.vst, "bugs", transformation = 'Vanilla')
run_baseline_DA(bugkoRNADNAs.pcl.vst, "bugkoRNADNAs", transformation = 'Vanilla')
run_baseline_DA(proteins.pcl.vst, "proteins", transformation = 'Vanilla')
run_baseline_DA(proteins.ecs.pcl.vst, "proteinECs", transformation = 'Vanilla')
run_baseline_DA(proteins.kos.pcl.vst, "proteinKOs", transformation = 'Vanilla')
run_baseline_DA(pwy.dna.unstrat.pcl.vst, "pwyDNAs", transformation = 'Vanilla')
run_baseline_DA(pwy.rna.unstrat.pcl.vst, "pwyRNAs", transformation = 'Vanilla')
run_baseline_DA(pwy.rnadna.unstrat.pcl.vst, "pwyRNADNAs", transformation = 'Vanilla')
run_baseline_DA(ec.dna.unstrat.pcl.vst, "ecDNAs", transformation = 'Vanilla')
run_baseline_DA(ec.rna.unstrat.pcl.vst, "ecRNAs", transformation = 'Vanilla')
run_baseline_DA(ko.dna.unstrat.pcl.filtered.vst, "koDNAs", transformation = 'Vanilla')
run_baseline_DA(ko.rna.unstrat.pcl.filtered.vst, "koRNAs", transformation = 'Vanilla')
run_baseline_DA(metabolites.pcl.vst, "metabolites", transformation = 'Vanilla')
library(lme4); run_baseline_DA_binomial(viruses.pcl.vst, "viruses", transformation = 'Vanilla')


############################################
# Analysis 4: Perform DA tests (PreActive) #
############################################

# Load Prective Indices
PreActive_mgx<-data.table::fread(paste(HMP2_root,'/analysis/active_disease_definition/active_samples_mgx.tsv', sep=''), header=FALSE)
PreActive_mgx<-as.data.frame(PreActive_mgx[PreActive_mgx$V2!='Active',])
PreActive_mgx<-column_to_rownames(PreActive_mgx, 'V1')
PreActive_mtx<-data.table::fread(paste(HMP2_root,'/analysis/active_disease_definition/active_samples_mtx.tsv', sep=''), header=FALSE)
PreActive_mtx<-as.data.frame(PreActive_mtx[PreActive_mtx$V2!='Active',])
PreActive_mtx<-column_to_rownames(PreActive_mtx, 'V1')
PreActive_mpx<-data.table::fread(paste(HMP2_root,'/analysis/active_disease_definition/active_samples_mpx.tsv', sep=''), header=FALSE)
PreActive_mpx<-as.data.frame(PreActive_mpx[PreActive_mpx$V2!='Active',])
PreActive_mpx<-column_to_rownames(PreActive_mpx, 'V1')
PreActive_mbx<-data.table::fread(paste(HMP2_root,'/analysis/active_disease_definition/active_samples_mbx.tsv', sep=''), header=FALSE)
PreActive_mbx<-as.data.frame(PreActive_mbx[PreActive_mbx$V2!='Active',])
PreActive_mbx<-column_to_rownames(PreActive_mbx, 'V1')
PreActive_ratio<-subset(PreActive_mgx, rownames(PreActive_mgx) %in% intersect(rownames(PreActive_mgx), rownames(PreActive_mtx)))


# Perform PreActive VS InActive DA analysis

# Just run for Vanilla

run_DA_PreActive(species.pcl.vst, "bugs", PreActive_mgx, transformation = 'Vanilla')
run_DA_PreActive(bugkoRNADNAs.pcl.vst, "bugkoRNADNAs", PreActive_ratio, transformation = 'Vanilla')
# run_DA_PreActive(proteins.pcl.vst, "proteins", PreActive_mpx, transformation = 'Vanilla') # Currently does not work due to mismatch in row names
run_DA_PreActive(proteins.ecs.pcl.vst, "proteinECs", PreActive_mpx, transformation = 'Vanilla')
run_DA_PreActive(proteins.kos.pcl.vst, "proteinKOs", PreActive_mpx, transformation = 'Vanilla')
run_DA_PreActive(pwy.dna.unstrat.pcl.vst, "pwyDNAs", PreActive_mgx, transformation = 'Vanilla')
run_DA_PreActive(pwy.rna.unstrat.pcl.vst, "pwyRNAs", PreActive_mtx, transformation = 'Vanilla')
run_DA_PreActive(pwy.rnadna.unstrat.pcl.vst, "pwyRNADNAs", PreActive_ratio, transformation = 'Vanilla')
run_DA_PreActive(ec.dna.unstrat.pcl.vst, "ecDNAs", PreActive_mgx, transformation = 'Vanilla')
run_DA_PreActive(ec.rna.unstrat.pcl.vst, "ecRNAs", PreActive_mtx, transformation = 'Vanilla')
run_DA_PreActive(ko.dna.unstrat.pcl.filtered.vst, "koDNAs", PreActive_mgx, transformation = 'Vanilla')
run_DA_PreActive(ko.rna.unstrat.pcl.filtered.vst, "koRNAs", PreActive_mtx, transformation = 'Vanilla')
run_DA_PreActive(metabolites.pcl.vst, "metabolites", PreActive_mbx, transformation = 'Vanilla')



###################################
# Some Host Analyses for Response #
###################################

#######
# HBI #
#######

# Just run for Vanilla

run_DA_HBI(species.pcl.vst, "bugs", transformation = 'Vanilla')
run_DA_HBI(bugkoRNADNAs.pcl.vst, "bugkoRNADNAs", transformation = 'Vanilla')
run_DA_HBI(proteins.pcl.vst, "proteins", transformation = 'Vanilla')
run_DA_HBI(proteins.ecs.pcl.vst, "proteinECs", transformation = 'Vanilla')
run_DA_HBI(proteins.kos.pcl.vst, "proteinKOs", transformation = 'Vanilla')
run_DA_HBI(pwy.dna.unstrat.pcl.vst, "pwyDNAs", transformation = 'Vanilla')
run_DA_HBI(pwy.rna.unstrat.pcl.vst, "pwyRNAs", transformation = 'Vanilla')
run_DA_HBI(pwy.rnadna.unstrat.pcl.vst, "pwyRNADNAs", transformation = 'Vanilla')
run_DA_HBI(ec.dna.unstrat.pcl.vst, "ecDNAs", transformation = 'Vanilla')
run_DA_HBI(ec.rna.unstrat.pcl.vst, "ecRNAs", transformation = 'Vanilla')
run_DA_HBI(ko.dna.unstrat.pcl.filtered.vst, "koDNAs", transformation = 'Vanilla')
run_DA_HBI(ko.rna.unstrat.pcl.filtered.vst, "koRNAs", transformation = 'Vanilla')
run_DA_HBI(metabolites.pcl.vst, "metabolites", transformation = 'Vanilla')
library(lme4); run_DA_HBI_binomial(viruses.pcl.vst, "viruses",  transformation = 'Vanilla')

#################
# HBI Objective #
#################

# Just run for Vanilla

run_DA_HBIobjective(species.pcl.vst, "bugs", transformation = 'Vanilla')
run_DA_HBIobjective(bugkoRNADNAs.pcl.vst, "bugkoRNADNAs", transformation = 'Vanilla')
run_DA_HBIobjective(proteins.pcl.vst, "proteins", transformation = 'Vanilla')
run_DA_HBIobjective(proteins.ecs.pcl.vst, "proteinECs", transformation = 'Vanilla')
run_DA_HBIobjective(proteins.kos.pcl.vst, "proteinKOs", transformation = 'Vanilla')
run_DA_HBIobjective(pwy.dna.unstrat.pcl.vst, "pwyDNAs", transformation = 'Vanilla')
run_DA_HBIobjective(pwy.rna.unstrat.pcl.vst, "pwyRNAs", transformation = 'Vanilla')
run_DA_HBIobjective(pwy.rnadna.unstrat.pcl.vst, "pwyRNADNAs", transformation = 'Vanilla')
run_DA_HBIobjective(ec.dna.unstrat.pcl.vst, "ecDNAs", transformation = 'Vanilla')
run_DA_HBIobjective(ec.rna.unstrat.pcl.vst, "ecRNAs", transformation = 'Vanilla')
run_DA_HBIobjective(ko.dna.unstrat.pcl.filtered.vst, "koDNAs", transformation = 'Vanilla')
run_DA_HBIobjective(ko.rna.unstrat.pcl.filtered.vst, "koRNAs", transformation = 'Vanilla')
run_DA_HBIobjective(metabolites.pcl.vst, "metabolites", transformation = 'Vanilla')
library(lme4); run_DA_HBIobjective_binomial(viruses.pcl.vst, "viruses",  transformation = 'Vanilla')


#########
# SCCAI #
#########

# Just run for Vanilla

run_DA_SCCAI(species.pcl.vst, "bugs", transformation = 'Vanilla')
run_DA_SCCAI(bugkoRNADNAs.pcl.vst, "bugkoRNADNAs", transformation = 'Vanilla')
run_DA_SCCAI(proteins.pcl.vst, "proteins", transformation = 'Vanilla')
run_DA_SCCAI(proteins.ecs.pcl.vst, "proteinECs", transformation = 'Vanilla')
run_DA_SCCAI(proteins.kos.pcl.vst, "proteinKOs", transformation = 'Vanilla')
run_DA_SCCAI(pwy.dna.unstrat.pcl.vst, "pwyDNAs", transformation = 'Vanilla')
run_DA_SCCAI(pwy.rna.unstrat.pcl.vst, "pwyRNAs", transformation = 'Vanilla')
run_DA_SCCAI(pwy.rnadna.unstrat.pcl.vst, "pwyRNADNAs", transformation = 'Vanilla')
run_DA_SCCAI(ec.dna.unstrat.pcl.vst, "ecDNAs", transformation = 'Vanilla')
run_DA_SCCAI(ec.rna.unstrat.pcl.vst, "ecRNAs", transformation = 'Vanilla')
run_DA_SCCAI(ko.dna.unstrat.pcl.filtered.vst, "koDNAs", transformation = 'Vanilla')
run_DA_SCCAI(ko.rna.unstrat.pcl.filtered.vst, "koRNAs", transformation = 'Vanilla')
run_DA_SCCAI(metabolites.pcl.vst, "metabolites", transformation = 'Vanilla')
library(lme4); run_DA_SCCAI_binomial(viruses.pcl.vst, "viruses",  transformation = 'Vanilla')


##############################################################
# Separate Analysis for Sena (Don't run - Takes a Long Time) #
##############################################################
# Just run for Vanilla

# Transformation
# run_DA(metabolites_Sena.pcl.vst, "metabolites_Sena", 'Vanilla')

###################################################################
# CROSSED VS NESTED RANDOM EFFECTS COMPARISON - EXAMPLE WITH BUGS #
###################################################################

# # The following produce identical results (except the last one)
# pcl<-species.pcl.vst; transformation = 'Vanilla'
#
# # Add antibiotics
# pcl <- merge_metadata(pcl, fields=c("Antibiotics"))
#
# # Transformation
# pcl<-pcl.transform(pcl, name, transformation)
#
# # Extract relevant metadata
# df <- pcl$meta
# df$active<-as.factor(df$active*1) # Convert to binary
# df$active <- relevel(df$active, ref = "0") # Set 'inactive' as reference
# i<-1
# df$y <- pcl$x[rownames(df),i]
#
# lme(y ~ diagnosis + diagnosis %in% active + Antibiotics + consent_age, random= ~ 1 | site_name/subject, data=df, na.action = na.omit)
# lme(y ~ diagnosis + diagnosis %in% active + Antibiotics + consent_age, random= list(~ 1 | site_name, ~1 | subject), data=df, na.action = na.omit)
# lmer(y ~ diagnosis + diagnosis %in% active + Antibiotics + consent_age + (1 | site_name/subject), data=df, na.action = na.omit)
# lmer(y ~ diagnosis + diagnosis %in% active + Antibiotics + consent_age + (1 | site_name) + (1|subject), data=df, na.action = na.omit)
# lme(formula, random=pdBlocked(list(pdIdent(~ 0 + site_name),pdIdent(~ 0 + subject))), data = df, na.action = na.omit) # DOES NOT WORK!!

####################
# CLR Sanity Check #
####################

# D<-matrix(c(0.1, 0.2, 0.7, 0.4, 1, 0.1, 0.1, 2.6, 0.3, 2, 5, 6), nrow = 3, byrow=TRUE)
# Hotelling::clr(D)
# t(SpiecEasi::clr(D, mar = 1, base = exp(1)))
# chemometrics::clr(D)
# pseudocount_col <-c(1, 2, 3, 4)
# pseudocount_row <-c(1, 2, 3)
# D + rep(0.5*pseudocount_col, each = nrow(D))
# D + rep(0.5*pseudocount_row, ncol(D))
# D

