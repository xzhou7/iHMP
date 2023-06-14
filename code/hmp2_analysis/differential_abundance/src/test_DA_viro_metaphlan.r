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

# Load required utility functions
source("./differential_abundance/src/core_DA_functions.R")
source("./common/merge_metadata.r")
source("./common/match_datasets.r")
source("./common/disease_activity.r") # By default merges disease activity info in the bugs table.

# Binarize and filter out all-zero or near-zero viral features
viruses.pcl$x<-ifelse(viruses.pcl$x > 0, 1, 0)
viruses.pcl.vst <- viruses.pcl %>% pcl.filter.f(keep=which(colMeans(viruses.pcl$x) > 0.05))
viruses.pcl.vst <- pcl.map.fnames(viruses.pcl.vst, gsub("^.*?([^\\|]+)$", "\\1", Name))
colnames(viruses.pcl.vst$x)<-gsub(".*s__", "", colnames(viruses.pcl.vst$x))

################################
# Analysis 1: Perform DA tests #
################################
library(lme4); run_DA_binomial(viruses.pcl.vst, "viruses_metaphlan", transformation = "Vanilla")

#################################################
# Analysis 2: Perform DA tests (Activity Index) #
#################################################
library(lme4); run_DA_activity_index_binomial(viruses.pcl.vst, "viruses_metaphlan", transformation = "Vanilla")

###########################################
# Analysis 3: Perform DA tests (baseline) #
###########################################
library(lme4); run_baseline_DA_binomial(viruses.pcl.vst, "viruses_metaphlan", transformation = "Vanilla")




