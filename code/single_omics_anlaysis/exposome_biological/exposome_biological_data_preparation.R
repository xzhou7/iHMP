##avoid source
no_function()

##CPM
masstools::setwd_project()
library(tidyverse)
setwd("data_20200511/microbiome/")
rm(list = ls())

dna_family <- read.table("Exposome_select_peng_DNA-FAMILY.txt", sep = ",", header = TRUE)
dna_genus <- read.table("Exposome_select_peng_DNA-GENUS.txt", sep = ",", header = TRUE)
dna_phylum <- read.table("Exposome_select_peng_DNA-PHYLUM.txt", sep = ",", header = TRUE)
dna_species <- read.table("Exposome_select_peng_DNA-SPECIES.txt", sep = ",", header = TRUE)
dna_sample_info <- readr::read_csv("Exposome_select_peng_metadata-DNA.csv")

rna_phylum <- read.table("Exposome_select_peng_RNA-PHYLUM.txt", sep = ",", header = TRUE)
rna_family <- read.table("Exposome_select_peng_RNA-FAMILY.txt", sep = ",", header = TRUE)
rna_genus <- read.table("Exposome_select_peng_RNA-GENUS.txt", sep = ",", header = TRUE)
rna_species <- read.table("Exposome_select_peng_RNA-SPECIES.txt", sep = ",", header = TRUE)
rna_sample_info <- readr::read_csv("Exposome_select_peng_metadata-RNA.csv")


##DNA
colnames(dna_family)
colnames(dna_genus)
colnames(dna_phylum)
colnames(dna_species)

dna_sample_info$samplenames

dna_family <- 
dna_family %>% 
  tibble::column_to_rownames(var = "X") %>% 
  as.data.frame()

dna_genus <- 
  dna_genus %>% 
  tibble::column_to_rownames(var = "X") %>% 
  as.data.frame()

dna_phylum <- 
  dna_phylum %>% 
  tibble::column_to_rownames(var = "X") %>% 
  as.data.frame()

dna_species <- 
  dna_species %>% 
  tibble::column_to_rownames(var = "X") %>% 
  as.data.frame()

# ##change to percentage
# dna_family <- 
# dna_family %>% 
#   apply(2, function(x){
#     x*100/sum(x)
#   }) %>% 
#   as.data.frame()
# 
# dna_genus <- 
#   dna_genus %>% 
#   apply(2, function(x){
#     x*100/sum(x)
#   }) %>% 
#   as.data.frame()
# 
# dna_phylum <- 
#   dna_phylum %>% 
#   apply(2, function(x){
#     x*100/sum(x)
#   }) %>% 
#   as.data.frame()
# 
# dna_species <- 
#   dna_species %>% 
#   apply(2, function(x){
#     x*100/sum(x)
#   }) %>% 
#   as.data.frame()

 dna_variable_info <- 
  rbind(
    data.frame(id = rownames(dna_family), class = 'family', stringsAsFactors = FALSE),
    data.frame(id = rownames(dna_genus), class = 'genus', stringsAsFactors = FALSE),
    data.frame(id = rownames(dna_phylum), class = 'phylum', stringsAsFactors = FALSE),
    data.frame(id = rownames(dna_species), class = 'species', stringsAsFactors = FALSE)
  )

sum(
  dna_variable_info$id == c(rownames(dna_family),
                            rownames(dna_genus),
                            rownames(dna_phylum),
                            rownames(dna_species))  
)

dna_variable_info

dna_expression_data <-
  rbind(dna_family, dna_genus, dna_phylum, dna_species) %>%
  as.data.frame()

dna_sample_info <- 
  dna_sample_info %>% 
  dplyr::select(-X1) %>% 
  dplyr::mutate(sample_id = samplenames) %>% 
  dplyr::select(-samplenames) %>% 
  dplyr::select(sample_id, everything())


dna_sample_info$sample_id == colnames(dna_expression_data)

dna_variable_info <-
  dna_variable_info %>% 
  dplyr::mutate(variable_id = paste(class, id, sep = "_")) %>% 
    dplyr::select(variable_id, short_name = id, level = class)

dna_variable_info$short_name == rownames(dna_expression_data)

rownames(dna_expression_data) <- dna_variable_info$variable_id

colnames(dna_expression_data) == dna_sample_info$sample_id




library(openxlsx)
wb = createWorkbook()
modifyBaseFont(wb, fontSize = 12, fontName = "Arial Narrow")
addWorksheet(wb, sheetName = "Sample information", gridLines = TRUE)
addWorksheet(wb, sheetName = "Variable information", gridLines = TRUE)
addWorksheet(wb, sheetName = "Expression data", gridLines = TRUE)
freezePane(wb, sheet = 1, firstRow = TRUE, firstCol = TRUE) 
freezePane(wb, sheet = 2, firstRow = TRUE, firstCol = TRUE) 
freezePane(wb, sheet = 3, firstRow = TRUE, firstCol = FALSE) 
writeDataTable(wb, sheet = 1, x = dna_sample_info,
               colNames = TRUE, rowNames = FALSE)
writeDataTable(wb, sheet = 2, x = dna_variable_info,
               colNames = TRUE, rowNames = FALSE)
writeDataTable(wb, sheet = 3, x = dna_expression_data,
               colNames = TRUE, rowNames = FALSE)
saveWorkbook(wb, "exposome_biological_data.xlsx", overwrite = TRUE)



save(dna_expression_data, file = "dna_expression_data")
save(dna_sample_info, file = "dna_sample_info")
save(dna_variable_info, file = "dna_variable_info")



##rna
colnames(rna_family)
colnames(rna_genus)
colnames(rna_phylum)
colnames(rna_species)

rna_sample_info$samplenames

rna_family <- 
  rna_family %>% 
  tibble::column_to_rownames(var = "X") %>% 
  as.data.frame()

rna_genus <- 
  rna_genus %>% 
  tibble::column_to_rownames(var = "X") %>% 
  as.data.frame()

rna_phylum <- 
  rna_phylum %>% 
  tibble::column_to_rownames(var = "X") %>% 
  as.data.frame()

rna_species <- 
  rna_species %>% 
  tibble::column_to_rownames(var = "X") %>% 
  as.data.frame()

##change to percentage
# rna_family <- 
#   rna_family %>% 
#   apply(2, function(x){
#     x*100/sum(x)
#   }) %>% 
#   as.data.frame()
# 
# rna_genus <- 
#   rna_genus %>% 
#   apply(2, function(x){
#     x*100/sum(x)
#   }) %>% 
#   as.data.frame()
# 
# rna_phylum <- 
#   rna_phylum %>% 
#   apply(2, function(x){
#     x*100/sum(x)
#   }) %>% 
#   as.data.frame()
# 
# rna_species <- 
#   rna_species %>% 
#   apply(2, function(x){
#     x*100/sum(x)
#   }) %>% 
#   as.data.frame()

rna_variable_info <- 
  rbind(
    data.frame(id = rownames(rna_family), class = 'family', stringsAsFactors = FALSE),
    data.frame(id = rownames(rna_genus), class = 'genus', stringsAsFactors = FALSE),
    data.frame(id = rownames(rna_phylum), class = 'phylum', stringsAsFactors = FALSE),
    data.frame(id = rownames(rna_species), class = 'species', stringsAsFactors = FALSE)
  )

sum(
  rna_variable_info$id == c(rownames(rna_family),
                            rownames(rna_genus),
                            rownames(rna_phylum),
                            rownames(rna_species))  
)

rna_variable_info

rna_expression_data <-
  rbind(rna_family, rna_genus, rna_phylum, rna_species) %>%
  as.data.frame()

rna_sample_info <- 
  rna_sample_info %>% 
  dplyr::select(-X1) %>% 
  dplyr::mutate(sample_id = samplenames) %>% 
  dplyr::select(-samplenames) %>% 
  dplyr::select(sample_id, everything())

rna_sample_info$sample_id == colnames(rna_expression_data)

rna_variable_info <-
  rna_variable_info %>% 
  dplyr::mutate(variable_id = paste(class, id, sep = "_")) %>% 
  dplyr::select(variable_id, short_name = id, level = class)

rna_variable_info$short_name == rownames(rna_expression_data)

rownames(rna_expression_data) <- rna_variable_info$variable_id

colnames(rna_expression_data) == rna_sample_info$sample_id

save(rna_expression_data, file = "rna_expression_data")
save(rna_sample_info, file = "rna_sample_info")
save(rna_variable_info, file = "rna_variable_info")

colnames(dna_expression_data)
colnames(rna_expression_data)

dna_sample_info$ownership
rna_sample_info$ownership

dna_sample_info$date.start
rna_sample_info$date.start
















