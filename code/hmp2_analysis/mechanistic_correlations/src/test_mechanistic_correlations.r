#### Remove all variables from the workspace ####
rm(list = ls())

#### Set directory ####
setwd('/Users/hmallick/Dropbox (Personal)/Repos/hmp2_analysis')

###### Required package list ######
library(plyr); library(dplyr); library(vegan); library(car);
library(nlme); library(lme4); library(data.table);
library(tibble); library(stringr); library(magrittr)

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
source("./mechanistic_correlations/src/core_mechanistic_functions.R")

# Subset to species-level taxa only
species.pcl <- pcl.only(bugs.pcl, rank="s")

# Binarize and filter out all-zero or near-zero viral features
viruses.pcl$x<-ifelse(viruses.pcl$x > 0, 1, 0)
viruses.pcl.vst <- viruses.pcl %>% pcl.filter.f(keep=which(colMeans(viruses.pcl$x) > 0.05))
viruses.pcl.vst <- pcl.map.fnames(viruses.pcl.vst, gsub("^.*?([^\\|]+)$", "\\1", Name))
colnames(viruses.pcl.vst$x)<-gsub(".*s__", "", colnames(viruses.pcl.vst$x))

# Unique metabolites with names
metabolites.pcl.vst<-metabolites.named.pcl

# Filter KO DNAs and RNAs (flag KOs highly correlated with a single bug abundance, threshold 0.6)
ko.dna.unstrat.pcl.filtered<-filterKO(ko.dna.unstrat.pcl, species.pcl)
ko.rna.unstrat.pcl.filtered<-filterKO(ko.rna.unstrat.pcl, species.pcl)

# Create PCL for RNA/DNA ratios from KO bug roll ups
bugkoRNADNAs.pcl<-createPCL(bugs.fromko.rna.pcl, bugs.fromko.dna.pcl)

#######################################
# Transformation (Ad hoc VST) Summary #
#######################################
# Counts -> log (base 2) with pseudocount 1
# Proportions (between 0 and 1) -> Arc Sine Square Root
# Ratios -> log (base 2)

# Bugs
species.pcl.vst<-species.pcl
species.pcl.vst$x <- asin(sqrt(species.pcl.vst$x)) # proportions
species.pcl.vst <- pcl.map.fnames(species.pcl.vst, gsub("^.*?([^\\|]+)$", "\\1", Name))
colnames(species.pcl.vst$x)<-gsub(".*s__", "", colnames(species.pcl.vst$x))


# Other omics
metabolites.pcl.vst$x <- log2(metabolites.pcl.vst$x + 1) # counts
proteins.pcl.vst <- proteins.pcl
proteins.pcl.vst$x <- log2(proteins.pcl.vst$x + 1) # counts
proteins.ecs.pcl.vst <- proteins.ecs.pcl
proteins.ecs.pcl.vst$x <- log2(proteins.ecs.pcl.vst$x+1) # counts
proteins.kos.pcl.vst <- proteins.kos.pcl
proteins.kos.pcl.vst$x <- log2(proteins.kos.pcl.vst$x+1)  # counts
pwy.dna.unstrat.pcl.vst <- pwy.dna.unstrat.pcl
pwy.dna.unstrat.pcl.vst$x <- asin(sqrt(pwy.dna.unstrat.pcl.vst$x)) # proportions
pwy.rna.unstrat.pcl.vst <- pwy.rna.unstrat.pcl
pwy.rna.unstrat.pcl.vst$x <- asin(sqrt(pwy.rna.unstrat.pcl.vst$x)) # proportions
pwy.rnadna.unstrat.pcl.vst <- pwy.rnadna.unstrat.pcl
pwy.rnadna.unstrat.pcl.vst$x <- log2(pwy.rnadna.unstrat.pcl.vst$x) # ratios
pwy.rnadna.unstrat.pcl.vst$x[!is.finite(pwy.rnadna.unstrat.pcl.vst$x)] <- NA
ec.dna.unstrat.pcl.vst <- ec.dna.unstrat.pcl
ec.dna.unstrat.pcl.vst$x <- asin(sqrt(ec.dna.unstrat.pcl.vst$x)) # proportions
ec.rna.unstrat.pcl.vst <- ec.rna.unstrat.pcl
ec.rna.unstrat.pcl.vst$x <- asin(sqrt(ec.rna.unstrat.pcl.vst$x)) # proportions
ko.dna.unstrat.pcl.filtered.vst <- ko.dna.unstrat.pcl.filtered
ko.dna.unstrat.pcl.filtered.vst$x <- asin(sqrt(ko.dna.unstrat.pcl.filtered.vst$x)) # proportions
ko.rna.unstrat.pcl.filtered.vst <- ko.rna.unstrat.pcl.filtered
ko.rna.unstrat.pcl.filtered.vst$x <- asin(sqrt(ko.rna.unstrat.pcl.filtered.vst$x)) # proportions
bugkoRNADNAs.pcl.vst <- bugkoRNADNAs.pcl
bugkoRNADNAs.pcl.vst$x <- log2(bugkoRNADNAs.pcl.vst$x) # ratios
bugkoRNADNAs.pcl.vst$x[!is.finite(bugkoRNADNAs.pcl.vst$x)] <- NA

# Extract subject random effects and residuals for each data type (both full and empty models)
library(lme4)
bugs.gl<-extract_residuals(species.pcl.vst, "bugs")
bugkoRNADNAs.gl<-extract_residuals(bugkoRNADNAs.pcl.vst, "bugkoRNADNAs")
viruses.gl<-extract_residuals(viruses.pcl.vst, "viruses")
proteins.gl<-extract_residuals(proteins.pcl.vst, "proteins")
proteinECs.gl<-extract_residuals(proteins.ecs.pcl.vst, "proteinECs")
proteinKOs.gl<-extract_residuals(proteins.kos.pcl.vst, "proteinKOs")
pwyDNAs.gl<-extract_residuals(pwy.dna.unstrat.pcl.vst, "pwyDNAs")
pwyRNAs.gl<-extract_residuals(pwy.rna.unstrat.pcl.vst, "pwyRNAs")
pwyRNADNAs.gl<-extract_residuals(pwy.rnadna.unstrat.pcl.vst, "pwyRNADNAs")
ecDNAs.gl<-extract_residuals(ec.dna.unstrat.pcl.vst, "ecDNAs")
ecRNAs.gl<-extract_residuals(ec.rna.unstrat.pcl.vst, "ecRNAs")
koDNAs.gl<-extract_residuals(ko.dna.unstrat.pcl.filtered.vst, "koDNAs")
koRNAs.gl<-extract_residuals(ko.rna.unstrat.pcl.filtered.vst, "koRNAs")
metabolites.gl<-extract_residuals(metabolites.pcl.vst, "metabolites")

# Match residuals for HAllA and AllA pairwise association testing
match_residuals(bugs.gl, metabolites.gl)
match_residuals(bugkoRNADNAs.gl, metabolites.gl)
match_residuals(viruses.gl, metabolites.gl)
match_residuals(proteins.gl, metabolites.gl)
match_residuals(proteinECs.gl, metabolites.gl)
match_residuals(proteinKOs.gl, metabolites.gl)
match_residuals(pwyDNAs.gl, metabolites.gl)
match_residuals(pwyRNAs.gl, metabolites.gl)
match_residuals(pwyRNADNAs.gl, metabolites.gl)
match_residuals(ecDNAs.gl, metabolites.gl)
match_residuals(ecRNAs.gl, metabolites.gl)
match_residuals(koDNAs.gl, metabolites.gl)
match_residuals(koRNAs.gl, metabolites.gl)

# Added on January 4, 2018
match_residuals(proteinECs.gl, bugs.gl)
match_residuals(ecDNAs.gl, bugs.gl)
match_residuals(ecRNAs.gl, bugs.gl)

# Added on January 17, 2018
fecalcal.gl<-extract_residuals_fecalcal()
match_residuals(bugkoRNADNAs.gl, fecalcal.gl)
match_residuals(ecDNAs.gl, fecalcal.gl)
match_residuals(ecRNAs.gl, fecalcal.gl)
match_residuals(proteinECs.gl, fecalcal.gl)
match_residuals(bugs.gl, fecalcal.gl)
match_residuals(metabolites.gl, fecalcal.gl)






