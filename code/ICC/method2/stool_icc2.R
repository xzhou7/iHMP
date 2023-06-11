#' ---
#' title: "stool microbiome metabolome correlation"
#' author: 
#'   - name: "Xiaotao Shen" 
#'     url: https://www.shenxt.info/
#'     affiliation: Stanford School of Medicine
#' date: "`r Sys.Date()`"
#' site: distill::distill_website
#' output: 
#'   distill::distill_article:
#'     code_folding: false
#' ---

#+ r setup, echo=TRUE, eval = TRUE, include = TRUE

# no_function()
# set work directory

setwd(masstools::get_project_wd())
library(tidyverse)
rm(list = ls())

####function from wenyu nature paper
##Variance decomposition
#First build a custom function with a dataframe containing variables, A1C, SubjectID, CollectionDate and CL4
#Providing the function with ds and n (1:m for variables)
myVD <- function(ds, trans){
  library(lme4)
  htb = {} 
  for (subject in unique(ds$SubjectID)) {
    hset = ds[(ds$SubjectID == subject & ds$CL4 == "Healthy"), ]
    hset$CollectionDate <- as.Date(hset$CollectionDate, "%m/%d/%y")
    hset = hset[order(hset$CollectionDate),]
    hset$Days = difftime(hset$CollectionDate, hset$CollectionDate[1], units = "days")
    htb = rbind(htb, hset)
  }
  # htb[, 1:m] <- apply(htb[, 1:m], 2, function(x) as.numeric(x))
  htb$Days <- as.numeric(htb$Days)
  htb = as.data.frame(htb)
  vd = {}
  for (i in 8:ncol(htb)) {
    if (trans == "Original") {model <- tryCatch(lmer(scale(htb[, i]) ~ 1 + Days + A1C + SSPG + FPG + (1|SubjectID), data = htb, REML = FALSE), error=function(err) NA)}
    if (trans == "Log10") {model <- tryCatch(lmer(scale(log10(htb[, i])) ~ 1 + Days + A1C + SSPG + FPG + (1|SubjectID), data = htb, REML = FALSE), error=function(err) NA)}
    if (trans == "Arcsin") {model <- tryCatch(lmer(scale(asin(sqrt(htb[, i]))) ~ 1 + Days + A1C + SSPG + FPG + (1|SubjectID), data = htb, REML = FALSE), error=function(err) NA)}
    if (!is.na(model)) {
      v.var = as.data.frame(VarCorr(model), comp=c("Variance"))[, 4] #Extracting random effect variances
      v.var <- c(v.var, as.data.frame(anova(model))[, 2]) #Extracting Sum Sq for fixed effects
    } else {
      #for variable that has too many zeros or for any other error situations
      v.var = rep(NA, 6)
    }
    vd <- rbind(vd, v.var)
  }
  rownames(vd) <- colnames(htb)[-c(1:7)]
  vd <- as.data.frame(vd)
  colnames(vd) <- c("random_Subject", "random_Residual", "fix_Days", "fix_A1C", "fix_SSPG", "fix_FPG")
  return(vd)
} #End of function


####----------------------------------------------------------------------------
####asv level
#######
####load data
load(here::here("data/from_xin/Revision_MultiOmes_0509.RData"))
clinic.sc = merge(clinic.df, sc, by = "SubjectID")
clinic.sc[, 62:63] <- apply(clinic.sc[, 62:63], 2, function(x) as.numeric(x))
a1c.df = clinic.sc %>% 
  dplyr::select(c(SubjectID,SampleID,A1C, SSPG, FPG, CL4, CollectionDate))

stool_microbiome_expression_data = 
  data.table::fread(here::here("data/from_xin/Genus Table/ST/ASV_ST.csv")) %>% 
  dplyr::select(SampleID, everything()) %>% 
  dplyr::select(-c(V1:batch)) %>% 
  as.data.frame()

idx = 
  stool_microbiome_expression_data %>% 
  dplyr::select(-SampleID) %>% 
  apply(2, function(x){
    sum(x != 0)/nrow(stool_microbiome_expression_data)
  }) %>% 
  `>=`(0.1) %>% 
  which() %>% 
  `+`(1)

length(idx)

stool_microbiome_expression_data = 
  stool_microbiome_expression_data[,c(1, idx)]

#16S_stool
ds =
  stool_microbiome_expression_data %>% 
  dplyr::left_join(a1c.df, by = "SampleID") %>% 
  dplyr::filter(!is.na(SSPG) & 
                  !is.na(A1C) & 
                  !is.na(FPG)) %>% 
  dplyr::select(SubjectID, SampleID,A1C, SSPG, FPG, CL4, CollectionDate, everything())

vd_stool_asv2 = myVD(ds = ds, trans = "Arcsin")
dir.create("data_analysis/stool_microbiome/ICC/")
save(vd_stool_asv2, file = "data_analysis/stool_microbiome/ICC/vd_stool_asv2")

###load variable_info

load(here::here("data_analysis/stool_microbiome/data_preparation/variable_info"))

# variable_info =
# variable_info %>% 
#   dplyr::filter(variable_id %in% rownames(vd_stool_asv2))

####----------------------------------------------------------------------------
####Phylum level
#######
####load data
stool_microbiome_expression_data = 
  data.table::fread(here::here("data/from_xin/Genus Table/ST/Phylum_ST.csv")) %>% 
  dplyr::select(SampleID, everything()) %>% 
  dplyr::select(-c(V1:batch)) %>% 
  as.data.frame()

idx = 
which(colnames(stool_microbiome_expression_data) %in% variable_info$Phylum)
length(idx)

stool_microbiome_expression_data = 
  stool_microbiome_expression_data[,c(1, idx)]

#16S_stool
ds =
stool_microbiome_expression_data %>% 
  dplyr::left_join(a1c.df, by = "SampleID") %>% 
  dplyr::filter(!is.na(SSPG) & 
                  !is.na(A1C) & 
                  !is.na(FPG)) %>% 
  dplyr::select(SubjectID, SampleID,A1C, SSPG, FPG, CL4, CollectionDate, everything())

vd_stool_phylum2 = myVD(ds = ds, trans = "Arcsin")
dir.create("data_analysis/stool_microbiome/ICC/")
save(vd_stool_phylum2, file = "data_analysis/stool_microbiome/ICC/vd_stool_phylum2")


####----------------------------------------------------------------------------
####order level
#######
####load data
load(here::here("data/from_xin/Revision_MultiOmes_0509.RData"))
clinic.sc = merge(clinic.df, sc, by = "SubjectID")
clinic.sc[, 62:63] <- apply(clinic.sc[, 62:63], 2, function(x) as.numeric(x))
a1c.df = clinic.sc %>% 
  dplyr::select(c(SubjectID,SampleID,A1C, SSPG, FPG, CL4, CollectionDate))

stool_microbiome_expression_data = 
  data.table::fread(here::here("data/from_xin/Genus Table/ST/Order_ST.csv")) %>% 
  dplyr::select(SampleID, everything()) %>% 
  dplyr::select(-c(V1:batch)) %>% 
  as.data.frame()

idx = 
  which(colnames(stool_microbiome_expression_data) %in% variable_info$Order) 
  
stool_microbiome_expression_data = 
  stool_microbiome_expression_data[,c(1, idx)]

#16S_stool
ds =
  stool_microbiome_expression_data %>% 
  dplyr::left_join(a1c.df, by = "SampleID") %>% 
  dplyr::filter(!is.na(SSPG) & 
                  !is.na(A1C) & 
                  !is.na(FPG)) %>% 
  dplyr::select(SubjectID, SampleID,A1C, SSPG, FPG, CL4, CollectionDate, everything())

vd_stool_order2 = myVD(ds = ds, trans = "Arcsin")
dir.create("data_analysis/stool_microbiome/ICC/")
save(vd_stool_order2, file = "data_analysis/stool_microbiome/ICC/vd_stool_order2")

####----------------------------------------------------------------------------
####genus level
#######
####load data
load(here::here("data/from_xin/Revision_MultiOmes_0509.RData"))
clinic.sc = merge(clinic.df, sc, by = "SubjectID")
clinic.sc[, 62:63] <- apply(clinic.sc[, 62:63], 2, function(x) as.numeric(x))
a1c.df = clinic.sc %>% 
  dplyr::select(c(SubjectID,SampleID,A1C, SSPG, FPG, CL4, CollectionDate))

stool_microbiome_expression_data = 
  data.table::fread(here::here("data/from_xin/Genus Table/ST/Genus_ST.csv")) %>% 
  dplyr::select(SampleID, everything()) %>% 
  dplyr::select(-c(V1:batch)) %>% 
  as.data.frame()

idx = 
  which(colnames(stool_microbiome_expression_data) %in% variable_info$Genus)
  
length(idx)
stool_microbiome_expression_data = 
  stool_microbiome_expression_data[,c(1, idx)]

#16S_stool
ds =
  stool_microbiome_expression_data %>% 
  dplyr::left_join(a1c.df, by = "SampleID") %>% 
  dplyr::filter(!is.na(SSPG) & 
                  !is.na(A1C) & 
                  !is.na(FPG)) %>% 
  dplyr::select(SubjectID, SampleID,A1C, SSPG, FPG, CL4, CollectionDate, everything())

vd_stool_genus2 = myVD(ds = ds, trans = "Arcsin")
dir.create("data_analysis/stool_microbiome/ICC/")
save(vd_stool_genus2, file = "data_analysis/stool_microbiome/ICC/vd_stool_genus2")


####----------------------------------------------------------------------------
####family level
#######
####load data
load(here::here("data/from_xin/Revision_MultiOmes_0509.RData"))
clinic.sc = merge(clinic.df, sc, by = "SubjectID")
clinic.sc[, 62:63] <- apply(clinic.sc[, 62:63], 2, function(x) as.numeric(x))
a1c.df = clinic.sc %>% 
  dplyr::select(c(SubjectID,SampleID,A1C, SSPG, FPG, CL4, CollectionDate))

stool_microbiome_expression_data = 
  data.table::fread(here::here("data/from_xin/Genus Table/ST/Family_ST.csv")) %>% 
  dplyr::select(SampleID, everything()) %>% 
  dplyr::select(-c(V1:batch)) %>% 
  as.data.frame()

idx = 
  which(colnames(stool_microbiome_expression_data) %in% variable_info$Family)
  
length(idx)
stool_microbiome_expression_data = 
  stool_microbiome_expression_data[,c(1, idx)]

#16S_stool
ds =
  stool_microbiome_expression_data %>% 
  dplyr::left_join(a1c.df, by = "SampleID") %>% 
  dplyr::filter(!is.na(SSPG) & 
                  !is.na(A1C) & 
                  !is.na(FPG)) %>% 
  dplyr::select(SubjectID, SampleID,A1C, SSPG, FPG, CL4, CollectionDate, everything())

vd_stool_family2 = myVD(ds = ds, trans = "Arcsin")
dir.create("data_analysis/stool_microbiome/ICC/")
save(vd_stool_family2, file = "data_analysis/stool_microbiome/ICC/vd_stool_family2")


####----------------------------------------------------------------------------
####class level
#######
####load data
load(here::here("data/from_xin/Revision_MultiOmes_0509.RData"))
clinic.sc = merge(clinic.df, sc, by = "SubjectID")
clinic.sc[, 62:63] <- apply(clinic.sc[, 62:63], 2, function(x) as.numeric(x))
a1c.df = clinic.sc %>% 
  dplyr::select(c(SubjectID,SampleID,A1C, SSPG, FPG, CL4, CollectionDate))

stool_microbiome_expression_data = 
  data.table::fread(here::here("data/from_xin/Genus Table/ST/Class_ST.csv")) %>% 
  dplyr::select(SampleID, everything()) %>% 
  dplyr::select(-c(V1:batch)) %>% 
  as.data.frame()

idx = 
  which(colnames(stool_microbiome_expression_data) %in% variable_info$Class)
  
stool_microbiome_expression_data = 
  stool_microbiome_expression_data[,c(1, idx)]
length(idx)
#16S_stool
ds =
  stool_microbiome_expression_data %>% 
  dplyr::left_join(a1c.df, by = "SampleID") %>% 
  dplyr::filter(!is.na(SSPG) & 
                  !is.na(A1C) & 
                  !is.na(FPG)) %>% 
  dplyr::select(SubjectID, SampleID,A1C, SSPG, FPG, CL4, CollectionDate, everything())

vd_stool_class2 = myVD(ds = ds, trans = "Arcsin")
dir.create("data_analysis/stool_microbiome/ICC/")
save(vd_stool_class2, file = "data_analysis/stool_microbiome/ICC/vd_stool_class2")






