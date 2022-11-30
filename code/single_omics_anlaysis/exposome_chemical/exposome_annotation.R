masstools::setwd_project()
setwd("data_20200511/exposome/")
library(tidyverse)
rm(list=ls())

peak_table <-
  readxl::read_xlsx("SmallMoleculeProfiling_20161129Summary_neg2.xlsx", sheet = 2)

peak_table <- peak_table[-1,]

peak_table <- peak_table[,c(1,3,4)]

colnames(peak_table)[1] <- "name"
peak_table$rt <- peak_table$rt * 60

peak_table_pos <- 
  peak_table %>% 
  dplyr::filter(stringr::str_detect(name, "PM"))

peak_table_neg <- 
  peak_table %>% 
  dplyr::filter(stringr::str_detect(name, "NM"))

dir.create("annotation")
setwd("annotation/")

# write.csv(peak_table_pos, "peak_table_pos.csv", row.names = FALSE)
# 
# write.csv(peak_table_neg, "peak_table_neg.csv", row.names = FALSE)
library(metID)
result_pos <- identify_metabolites(
  ms1.data = "peak_table_pos.csv",
  polarity = "positive",
  ce = 'all',
  column = "rp",
  total.score.tol = 0.5,
  candidate.num = 3,
  threads = 3,
  database = "list1_ms1_database"
)

result_neg <- identify_metabolites(
  ms1.data = "peak_table_neg.csv",
  polarity = "negative",
  ce = 'all',
  column = "rp",
  total.score.tol = 0.5,
  candidate.num = 3,
  threads = 3,
  database = "list1_ms1_database"
)

annotation_table_pos <-
  get_identification_table(result_pos,
                           type = "old",
                           candidate.num = 3)

annotation_table_neg <-
  get_identification_table(result_neg,
                           type = "old",
                           candidate.num = 3)

annotation_table_pos$Identification <-
  annotation_table_pos$Identification %>% 
  lapply(function(x){
    if(is.na(x)){
      return(NA)
    }else{
      x <- stringr::str_split(x, "\\{\\}")[[1]] 
      x <- grep("\\(M\\+H\\)", x, value = TRUE)
      if(length(x) == 0){
        return(NA)
      }else{
        paste(x, collapse = "{}")
      }
    }
  }) %>% 
  unlist()

annotation_table_neg$Identification <-
  annotation_table_neg$Identification %>% 
  lapply(function(x){
    if(is.na(x)){
      return(NA)
    }else{
      x <- stringr::str_split(x, "\\{\\}")[[1]] 
      x <- grep("\\(M\\-H\\)", x, value = TRUE)
      if(length(x) == 0){
        return(NA)
      }else{
        paste(x, collapse = "{}")
      }
    }
  }) %>% 
  unlist()


annotation_table <- rbind(annotation_table_pos, annotation_table_neg)

annotation_table <-
  annotation_table %>% 
  dplyr::filter(!is.na(Identification))

write.csv(annotation_table, "annotation_table.csv", row.names = FALSE)
