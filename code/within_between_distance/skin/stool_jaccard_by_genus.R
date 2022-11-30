###no source
masstools::setwd_project()
rm(list = ls())
library(tidyverse)

source("code/tools.R")
load("data/from_xin/physeq_clean.rda")
load("data_analysis/stool_microbiome/data_preparation/sample_info")
family_info = readxl::read_xlsx("data/from_xin/Related People In iPOP_from.sophia.xlsx", 
                                col_names = FALSE) 

data = readr::read_csv("Figures/Figure3/strainreplacement/stool.replacement.csv")
setwd("data_analysis/combine_microbiome/distance/stool/")

jaccard_genus <- function(a, b) {
  a[a>0] = 1
  b[b>0] = 1
  temp_data = 
    data.frame(a, b) %>% 
    dplyr::filter(a != 0 | b != 0)
  if(nrow(temp_data) == 0){
    return(NA)
  }
  sum(rowSums(temp_data) == 2)/nrow(temp_data)
}

colnames(family_info) = c("family", "subject_id", "family_role")

family_info =
    family_info %>%
    dplyr::mutate(subject_id = stringr::str_replace(subject_id, "1636\\-", ""))

subject_info =
    sample_info %>%
    dplyr::select(subject_id,
                  IRIS,
                  SSPG,
                  FPG,
                  SSPG.Date,
                  Class,
                  Gender,
                  Ethnicity,
                  Adj.age)

library(phyloseq)
library(tidyverse)
library(cowplot)
library(parallel)

ls()

physeq_stool = physeq_ST

variable_info =
    phyloseq::tax_table(physeq_stool) %>%
    as.data.frame()

sample_info =
    suppressWarnings(as_tibble(sample_data(physeq_stool))) %>%
    dplyr::select(RandomID:Date, batch) %>% 
    dplyr::left_join(family_info, by = c("SubjectID" = "subject_id"))

allGenus = as.vector(tax_table(physeq_stool)[, "Genus"])

uniqueGenus = unique(allGenus)

message(length(uniqueGenus), " genus to calculate")

####remain the distance == 1 which are not 0 in both two samples
library(future)
library(tictoc)

stool_jaccard_by_genus = 
  furrr::future_map(
  as.data.frame(t(data)),
  .f = function(i) {
    cat(i[2], "_", i[3]," ")
    temp_genus = i[4]
    
    subdat = tryCatch(prune_taxa(allGenus == temp_genus, physeq_stool), 
                      error = function(e) return(NA))
    subdat = tryCatch(prune_taxa(taxa_sums(subdat) > 0, subdat), 
                      error = function(e) return(NA))
    subdat = tryCatch(prune_samples(sample_sums(subdat) > 0, subdat),
                      error = function(e) return(NA))
    subdat = tryCatch(prune_samples(sample_names(subdat) %in% sample_info$RandomID[match(i[2:3], sample_info$SampleID)], subdat),
                      error = function(e) return(NA))
    
    if(is.na(subdat)){
      jaccard = return(NA)
      asv_number = NA
    }else{
      temp_data = as.data.frame(subdat@otu_table)
      jaccard = 
        jaccard_genus(a = as.numeric(temp_data[1,]),
                      b = as.numeric(temp_data[2,]))
      temp_data2 = 
        data.frame(a = as.numeric(temp_data[1,]),
                   b = as.numeric(temp_data[2,])) %>% 
        dplyr::filter(a !=0 | b != 0)  
      asv_number = nrow(temp_data2)
    }
    data.frame(Start = i[2],
               End = i[3],
               taxa = i[4],
               SubjectID = i[5], 
               jaccard = jaccard,
               asv_number = asv_number)
  }, .progress = TRUE
) %>% 
  do.call(rbind, .) %>% 
  as.data.frame()

save(stool_jaccard_by_genus,
     file = "stool_jaccard_by_genus",
     compress = "xz")
