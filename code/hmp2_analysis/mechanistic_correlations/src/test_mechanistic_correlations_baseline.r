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

# Set of unique metabolites for which we have names
metabolites.pcl.vst<-metabolites.named.pcl

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
proteins.ecs.pcl.vst <- proteins.ecs.pcl
proteins.ecs.pcl.vst$x <- log2(proteins.ecs.pcl.vst$x+1) # counts
ec.dna.unstrat.pcl.vst <- ec.dna.unstrat.pcl
ec.dna.unstrat.pcl.vst$x <- asin(sqrt(ec.dna.unstrat.pcl.vst$x)) # proportions
ec.rna.unstrat.pcl.vst <- ec.rna.unstrat.pcl
ec.rna.unstrat.pcl.vst$x <- asin(sqrt(ec.rna.unstrat.pcl.vst$x)) # proportions
bugkoRNADNAs.pcl.vst <- bugkoRNADNAs.pcl
bugkoRNADNAs.pcl.vst$x <- log2(bugkoRNADNAs.pcl.vst$x) # ratios
bugkoRNADNAs.pcl.vst$x[!is.finite(bugkoRNADNAs.pcl.vst$x)] <- NA

# Extract residuals for each microbial data type (both full and empty models)
bugs.gl<-extract_residuals_baseline(species.pcl.vst, "bugs")
bugkoRNADNAs.gl<-extract_residuals_baseline(bugkoRNADNAs.pcl.vst, "bugkoRNADNAs")
proteinECs.gl<-extract_residuals_baseline(proteins.ecs.pcl.vst, "proteinECs")
ecDNAs.gl<-extract_residuals_baseline(ec.dna.unstrat.pcl.vst, "ecDNAs")
ecRNAs.gl<-extract_residuals_baseline(ec.rna.unstrat.pcl.vst, "ecRNAs")
metabolites.gl<-extract_residuals_baseline(metabolites.pcl.vst, "metabolites")


# Host data types to correlate: serology, fecalcal, HTX-rectum and HTX-ileum

# Serology
source("./common/load_serology.r")

# HTX

# Filter by DE features
library(readxl)
read_excel_allsheets <- function(filename) {
  sheets <- readxl::excel_sheets(filename)
  x <-    lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  names(x) <- sheets
  x
}

DElist<-read_excel_allsheets("/Users/hmallick/Dropbox (Huttenhower Lab)/HMP2/analysis/host_transcriptomics/HTX_Diff_Expression.xlsx")
Ileum_features<-unique(unlist(lapply(DElist[1:2], "[", 1), use.names=FALSE))
Rectum_features<-unique(unlist(lapply(DElist[3:5], "[", 1), use.names=FALSE))


source("./common/load_biopsies.r")

# HTX-ileum

HTXileum.pcl = biopsy_htx.counts.pcl %>% pcl.filter.s(biopsy_location=='Ileum' & !is.na(biopsy_location))
HTXileum.pcl$x<-HTXileum.pcl$x[, Ileum_features]
HTXileum.pcl$nf<-ncol(HTXileum.pcl$x)

# HTX-rectum
HTXrectum.pcl = biopsy_htx.counts.pcl %>% pcl.filter.s(biopsy_location=='Rectum' & !is.na(biopsy_location))
HTXrectum.pcl$x<-HTXrectum.pcl$x[, Rectum_features]
HTXrectum.pcl$nf<-ncol(HTXrectum.pcl$x)


# VST
serology.pcl.vst<-serology.eu.pcl
serology.pcl.vst$x <- log2(serology.pcl.vst$x + 1)
HTXrectum.pcl.vst<-HTXrectum.pcl
HTXrectum.pcl.vst$x <- log2(HTXrectum.pcl.vst$x + 1)
HTXileum.pcl.vst<-HTXileum.pcl
HTXileum.pcl.vst$x <- log2(HTXileum.pcl.vst$x + 1)

# Extract residuals
serology.gl<-extract_residuals_baseline(serology.pcl.vst, "serology")
HTXrectum.gl<-extract_residuals_baseline(HTXrectum.pcl.vst, "HTXrectum")
HTXileum.gl<-extract_residuals_baseline(HTXileum.pcl, "HTXileum")

# Match residuals
match_residuals_baseline(bugs.gl, serology.gl)
match_residuals_baseline(bugs.gl, HTXrectum.gl)
match_residuals_baseline(bugs.gl, HTXileum.gl)
match_residuals_baseline(metabolites.gl, serology.gl)
match_residuals_baseline(metabolites.gl, HTXrectum.gl)
match_residuals_baseline(metabolites.gl, HTXileum.gl)
match_residuals_baseline(ecRNAs.gl, serology.gl)
match_residuals_baseline(ecRNAs.gl, HTXrectum.gl)
match_residuals_baseline(ecRNAs.gl, HTXileum.gl)
match_residuals_baseline(ecDNAs.gl, serology.gl)
match_residuals_baseline(ecDNAs.gl, HTXrectum.gl)
match_residuals_baseline(ecDNAs.gl, HTXileum.gl)
match_residuals_baseline(proteinECs.gl, serology.gl)
match_residuals_baseline(proteinECs.gl, HTXrectum.gl)
match_residuals_baseline(proteinECs.gl, HTXileum.gl)
match_residuals_baseline(bugkoRNADNAs.gl, serology.gl)
match_residuals_baseline(bugkoRNADNAs.gl, HTXrectum.gl)
match_residuals_baseline(bugkoRNADNAs.gl, HTXileum.gl)

