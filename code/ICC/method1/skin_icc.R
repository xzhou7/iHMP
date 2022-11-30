#' ---
#' title: "nasal microbiome metabolome correlation"
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

no_function()
# set work directory

masstools::setwd_project()
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
####Phylum level
#######
####load data
load(here::here("data/from_xin/Revision_MultiOmes_0509.RData"))
clinic.sc = merge(clinic.df, sc, by = "SubjectID")
clinic.sc[, 62:63] <- apply(clinic.sc[, 62:63], 2, function(x) as.numeric(x))
a1c.df = clinic.sc %>% 
  dplyr::select(c(SubjectID,SampleID,A1C, SSPG, FPG, CL4, CollectionDate))

skin_microbiome_expression_data = 
  data.table::fread(here::here("data/from_xin/Genus Table/SK/Phylum_SK.csv")) %>% 
  dplyr::select(SampleID, everything()) %>% 
  dplyr::select(-c(V1:SubjectID)) %>% 
  as.data.frame()

idx = 
  skin_microbiome_expression_data %>% 
  dplyr::select(-SampleID) %>% 
  apply(2, function(x){
    sum(x != 0)/(nrow(skin_microbiome_expression_data))
  }) %>% 
  `>=`(0.2) %>% 
  which() %>% 
  `+`(1)
length(idx)
skin_microbiome_expression_data = 
  skin_microbiome_expression_data[,c(1, idx)]

#16S_nasal
ds =
  skin_microbiome_expression_data %>% 
  dplyr::left_join(a1c.df, by = "SampleID") %>% 
  dplyr::filter(!is.na(SSPG) & 
                  !is.na(A1C) & 
                  !is.na(FPG)) %>% 
  dplyr::select(SubjectID, SampleID,A1C, SSPG, FPG, CL4, CollectionDate, everything())

vd_skin_phylum = myVD(ds = ds, trans = "Arcsin")
dir.create("data_analysis/skin_microbiome/ICC/")
save(vd_skin_phylum, file = "data_analysis/skin_microbiome/ICC/vd_skin_phylum")


####----------------------------------------------------------------------------
####order level
#######
####load data
load(here::here("data/from_xin/Revision_MultiOmes_0509.RData"))
clinic.sc = merge(clinic.df, sc, by = "SubjectID")
clinic.sc[, 62:63] <- apply(clinic.sc[, 62:63], 2, function(x) as.numeric(x))
a1c.df = clinic.sc %>% 
  dplyr::select(c(SubjectID,SampleID,A1C, SSPG, FPG, CL4, CollectionDate))

skin_microbiome_expression_data = 
  data.table::fread(here::here("data/from_xin/Genus Table/SK/Order_SK.csv")) %>% 
  dplyr::select(SampleID, everything()) %>% 
  dplyr::select(-c(V1:SubjectID)) %>% 
  as.data.frame()

idx = 
  skin_microbiome_expression_data %>% 
  dplyr::select(-SampleID) %>% 
  apply(2, function(x){
    sum(x != 0)/(nrow(skin_microbiome_expression_data))
  }) %>% 
  `>=`(0.2) %>% 
  which() %>% 
  `+`(1)

skin_microbiome_expression_data = 
  skin_microbiome_expression_data[,c(1, idx)]

#16S_nasal
ds =
  skin_microbiome_expression_data %>% 
  dplyr::left_join(a1c.df, by = "SampleID") %>% 
  dplyr::filter(!is.na(SSPG) & 
                  !is.na(A1C) & 
                  !is.na(FPG)) %>% 
  dplyr::select(SubjectID, SampleID,A1C, SSPG, FPG, CL4, CollectionDate, everything())

vd_skin_order = myVD(ds = ds, trans = "Arcsin")
dir.create("data_analysis/skin_microbiome/ICC/")
save(vd_skin_order, file = "data_analysis/skin_microbiome/ICC/vd_skin_order")


####----------------------------------------------------------------------------
####genus level
#######
####load data
load(here::here("data/from_xin/Revision_MultiOmes_0509.RData"))
clinic.sc = merge(clinic.df, sc, by = "SubjectID")
clinic.sc[, 62:63] <- apply(clinic.sc[, 62:63], 2, function(x) as.numeric(x))
a1c.df = clinic.sc %>% 
  dplyr::select(c(SubjectID,SampleID,A1C, SSPG, FPG, CL4, CollectionDate))

skin_microbiome_expression_data = 
  data.table::fread(here::here("data/from_xin/Genus Table/SK/Genus_SK.csv")) %>% 
  dplyr::select(SampleID, everything()) %>% 
  dplyr::select(-c(V1:SubjectID)) %>% 
  as.data.frame()

idx = 
  skin_microbiome_expression_data %>% 
  dplyr::select(-SampleID) %>% 
  apply(2, function(x){
    sum(x != 0)/(nrow(skin_microbiome_expression_data))
  }) %>% 
  `>=`(0.2) %>% 
  which() %>% 
  `+`(1)

length(idx)
skin_microbiome_expression_data = 
  skin_microbiome_expression_data[,c(1, idx)]

#16S_nasal
ds =
  skin_microbiome_expression_data %>% 
  dplyr::left_join(a1c.df, by = "SampleID") %>% 
  dplyr::filter(!is.na(SSPG) & 
                  !is.na(A1C) & 
                  !is.na(FPG)) %>% 
  dplyr::select(SubjectID, SampleID,A1C, SSPG, FPG, CL4, CollectionDate, everything())

vd_skin_genus = myVD(ds = ds, trans = "Arcsin")
dir.create("data_analysis/skin_microbiome/ICC/")
save(vd_skin_genus, file = "data_analysis/skin_microbiome/ICC/vd_skin_genus")



####----------------------------------------------------------------------------
####family level
#######
####load data
load(here::here("data/from_xin/Revision_MultiOmes_0509.RData"))
clinic.sc = merge(clinic.df, sc, by = "SubjectID")
clinic.sc[, 62:63] <- apply(clinic.sc[, 62:63], 2, function(x) as.numeric(x))
a1c.df = clinic.sc %>% 
  dplyr::select(c(SubjectID,SampleID,A1C, SSPG, FPG, CL4, CollectionDate))

skin_microbiome_expression_data = 
  data.table::fread(here::here("data/from_xin/Genus Table/SK/Family_SK.csv")) %>% 
  dplyr::select(SampleID, everything()) %>% 
  dplyr::select(-c(V1:SubjectID)) %>% 
  as.data.frame()

idx = 
  skin_microbiome_expression_data %>% 
  dplyr::select(-SampleID) %>% 
  apply(2, function(x){
    sum(x != 0)/(nrow(skin_microbiome_expression_data))
  }) %>% 
  `>=`(0.2) %>% 
  which() %>% 
  `+`(1)
length(idx)
skin_microbiome_expression_data = 
  skin_microbiome_expression_data[,c(1, idx)]

#16S_nasal
ds =
  skin_microbiome_expression_data %>% 
  dplyr::left_join(a1c.df, by = "SampleID") %>% 
  dplyr::filter(!is.na(SSPG) & 
                  !is.na(A1C) & 
                  !is.na(FPG)) %>% 
  dplyr::select(SubjectID, SampleID,A1C, SSPG, FPG, CL4, CollectionDate, everything())

vd_skin_family = myVD(ds = ds, trans = "Arcsin")
dir.create("data_analysis/skin_microbiome/ICC/")
save(vd_skin_family, file = "data_analysis/skin_microbiome/ICC/vd_skin_family")


####----------------------------------------------------------------------------
####class level
#######
####load data
load(here::here("data/from_xin/Revision_MultiOmes_0509.RData"))
clinic.sc = merge(clinic.df, sc, by = "SubjectID")
clinic.sc[, 62:63] <- apply(clinic.sc[, 62:63], 2, function(x) as.numeric(x))
a1c.df = clinic.sc %>% 
  dplyr::select(c(SubjectID,SampleID,A1C, SSPG, FPG, CL4, CollectionDate))

skin_microbiome_expression_data = 
  data.table::fread(here::here("data/from_xin/Genus Table/SK/Class_SK.csv")) %>% 
  dplyr::select(SampleID, everything()) %>% 
  dplyr::select(-c(V1:SubjectID)) %>% 
  as.data.frame()

idx = 
  skin_microbiome_expression_data %>% 
  dplyr::select(-SampleID) %>% 
  apply(2, function(x){
    sum(x != 0)/(nrow(skin_microbiome_expression_data))
  }) %>% 
  `>=`(0.2) %>% 
  which() %>% 
  `+`(1)
length(idx)
skin_microbiome_expression_data = 
  skin_microbiome_expression_data[,c(1, idx)]
length(idx)
#16S_nasal
ds =
  skin_microbiome_expression_data %>% 
  dplyr::left_join(a1c.df, by = "SampleID") %>% 
  dplyr::filter(!is.na(SSPG) & 
                  !is.na(A1C) & 
                  !is.na(FPG)) %>% 
  dplyr::select(SubjectID, SampleID,A1C, SSPG, FPG, CL4, CollectionDate, everything())

vd_skin_class = myVD(ds = ds, trans = "Arcsin")
dir.create("data_analysis/skin_microbiome/ICC/")
save(vd_skin_class, file = "data_analysis/skin_microbiome/ICC/vd_skin_class")


####----------------------------------------------------------------------------
####asv level
#######
####load data
load(here::here("data/from_xin/Revision_MultiOmes_0509.RData"))
clinic.sc = merge(clinic.df, sc, by = "SubjectID")
clinic.sc[, 62:63] <- apply(clinic.sc[, 62:63], 2, function(x) as.numeric(x))
a1c.df = clinic.sc %>% 
  dplyr::select(c(SubjectID,SampleID,A1C, SSPG, FPG, CL4, CollectionDate))

skin_microbiome_expression_data = 
  data.table::fread(here::here("data/from_xin/Genus Table/SK/ASV_SK.csv")) %>% 
  dplyr::select(SampleID, everything()) %>% 
  dplyr::select(-c(V1:SubjectID)) %>% 
  as.data.frame()

idx = 
  skin_microbiome_expression_data %>% 
  dplyr::select(-SampleID) %>% 
  apply(2, function(x){
    sum(x != 0)/nrow(skin_microbiome_expression_data)
  }) %>% 
  `>=`(0.2) %>% 
  which() %>% 
  `+`(1)
length(idx)
skin_microbiome_expression_data = 
  skin_microbiome_expression_data[,c(1, idx)]

#16S_nasal
ds =
  skin_microbiome_expression_data %>% 
  dplyr::left_join(a1c.df, by = "SampleID") %>% 
  dplyr::filter(!is.na(SSPG) & 
                  !is.na(A1C) & 
                  !is.na(FPG)) %>% 
  dplyr::select(SubjectID, SampleID,A1C, SSPG, FPG, CL4, CollectionDate, everything())

vd_skin_asv = myVD(ds = ds, trans = "Arcsin")
dir.create("data_analysis/skin_microbiome/ICC/")
save(vd_skin_asv, file = "data_analysis/skin_microbiome/ICC/vd_skin_asv")


