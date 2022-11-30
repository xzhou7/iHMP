#' ---
#' title: "Skin microbiome lipidome correlation"
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

source("code/tools.R")

######work directory
masstools::setwd_project()
setwd("data_analysis/correlation_network/whole_data_set/inter_omics_correlation_network/summary")

####load data
###skin microbiome
{
  load(here::here("data_analysis/skin_microbiome/data_preparation/expression_data"))
  load(here::here("data_analysis/skin_microbiome/data_preparation/sample_info"))
  load(here::here("data_analysis/skin_microbiome/data_preparation/variable_info"))
  
  skin_microbiome_expression_data = expression_data
  skin_microbiome_sample_info = sample_info
  skin_microbiome_variable_info = variable_info
  
  ###read genus table
  expression_data =
    data.table::fread(here::here("data/from_xin/Genus Table/SK/Genus_SK.csv")) %>% 
    as.data.frame() %>% 
    tibble::column_to_rownames(var = "SampleID") %>%
    dplyr::select(-c(V1:SubjectID)) %>%
    t() %>%
    as.data.frame()
  
  skin_microbiome_variable_info =
    skin_microbiome_variable_info[match(rownames(expression_data), skin_microbiome_variable_info$Genus),]
  
  skin_microbiome_variable_info$Genus == rownames(expression_data)
  
  ###remove the variables which Genus are NA
  remove_idx = which(is.na(skin_microbiome_variable_info$Genus))
  remove_idx
  if(length(remove_idx) > 0){
    skin_microbiome_variable_info = skin_microbiome_variable_info[-remove_idx,]
    expression_data = expression_data[-remove_idx,]
  }
  
  rownames(expression_data) = skin_microbiome_variable_info$variable_id
  
  skin_microbiome_expression_data =
    expression_data
  
  skin_microbiome_variable_info =
    skin_microbiome_variable_info %>%
    dplyr::filter(!stringr::str_detect(Genus, "Unclassified_Bacteria"))
  
  skin_microbiome_expression_data = 
    skin_microbiome_expression_data[skin_microbiome_variable_info$variable_id,]
  
  dim(skin_microbiome_sample_info)
  dim(skin_microbiome_variable_info)
  dim(skin_microbiome_expression_data)
  
  rownames(skin_microbiome_expression_data) == skin_microbiome_variable_info$variable_id
  colnames(skin_microbiome_expression_data) == skin_microbiome_sample_info$sample_id
  
  remian_subject_id =
    skin_microbiome_sample_info %>%
    dplyr::group_by(subject_id) %>%
    dplyr::summarise(n = n()) %>%
    dplyr::ungroup() %>%
    dplyr::filter(n > 5) %>%
    dplyr::pull(subject_id)

  skin_microbiome_sample_info =
    skin_microbiome_sample_info %>%
    dplyr::filter(subject_id %in% remian_subject_id)

  skin_microbiome_expression_data =
    skin_microbiome_expression_data[,skin_microbiome_sample_info$sample_id]
  
}

###stool microbiome
{
  load(here::here("data_analysis/stool_microbiome/data_preparation/expression_data"))
  load(here::here("data_analysis/stool_microbiome/data_preparation/sample_info"))
  load(here::here("data_analysis/stool_microbiome/data_preparation/variable_info"))
  
  stool_microbiome_expression_data = expression_data
  stool_microbiome_sample_info = sample_info
  stool_microbiome_variable_info = variable_info
  
  ###read genus table
  expression_data =
    data.table::fread(here::here("data/from_xin/Genus Table/ST/Genus_ST.csv")) %>% 
    as.data.frame() %>% 
    tibble::column_to_rownames(var = "SampleID") %>%
    dplyr::select(-c(V1:batch)) %>%
    t() %>%
    as.data.frame()
  
  stool_microbiome_variable_info =
    stool_microbiome_variable_info[match(rownames(expression_data), stool_microbiome_variable_info$Genus),]
  
  stool_microbiome_variable_info$Genus == rownames(expression_data)
  
  rownames(expression_data) = stool_microbiome_variable_info$variable_id
  
  stool_microbiome_expression_data =
    expression_data
  
  stool_microbiome_variable_info = 
    stool_microbiome_variable_info %>% 
    dplyr::filter(!stringr::str_detect(Genus, "Unclassified_Bacteria"))
  
  stool_microbiome_expression_data = 
    stool_microbiome_expression_data[stool_microbiome_variable_info$variable_id,]
  
  dim(stool_microbiome_sample_info)
  dim(stool_microbiome_variable_info)
  
  rownames(stool_microbiome_expression_data) == stool_microbiome_variable_info$variable_id
  colnames(stool_microbiome_expression_data) == stool_microbiome_sample_info$sample_id
  
  
  remian_subject_id =
    stool_microbiome_sample_info %>%
    dplyr::group_by(subject_id) %>%
    dplyr::summarise(n = n()) %>%
    dplyr::ungroup() %>%
    dplyr::filter(n > 5) %>%
    dplyr::pull(subject_id)
  
  stool_microbiome_sample_info =
    stool_microbiome_sample_info %>%
    dplyr::filter(subject_id %in% remian_subject_id)
  
  stool_microbiome_expression_data =
    stool_microbiome_expression_data[,stool_microbiome_sample_info$sample_id]
  
}

###nasal microbiome
{
  load(here::here("data_analysis/nasal_microbiome/data_preparation/expression_data"))
  load(here::here("data_analysis/nasal_microbiome/data_preparation/sample_info"))
  load(here::here("data_analysis/nasal_microbiome/data_preparation/variable_info"))
  
  
  nasal_microbiome_expression_data = expression_data
  nasal_microbiome_sample_info = sample_info
  nasal_microbiome_variable_info = variable_info
  
  ###read genus table
  expression_data =
    data.table::fread(here::here("data/from_xin/Genus Table/NS/Genus_NS.csv")) %>% 
    as.data.frame() %>% 
    tibble::column_to_rownames(var = "SampleID") %>%
    dplyr::select(-c(V1:batch)) %>%
    t() %>%
    as.data.frame()
  
  nasal_microbiome_variable_info =
    nasal_microbiome_variable_info[match(rownames(expression_data), nasal_microbiome_variable_info$Genus),]
  
  nasal_microbiome_variable_info$Genus == rownames(expression_data)
  
  ###remove the variables which Genus are NA
  remove_idx = which(is.na(nasal_microbiome_variable_info$Genus))
  remove_idx
  if(length(remove_idx) > 0){
    nasal_microbiome_variable_info = nasal_microbiome_variable_info[-remove_idx,]
    expression_data = expression_data[-remove_idx,]
  }
  
  rownames(expression_data) = nasal_microbiome_variable_info$variable_id
  
  nasal_microbiome_expression_data =
    expression_data
  
  nasal_microbiome_variable_info = 
    nasal_microbiome_variable_info %>% 
    dplyr::filter(!stringr::str_detect(Genus, "Unclassified_Bacteria"))
  
  nasal_microbiome_expression_data = 
    nasal_microbiome_expression_data[nasal_microbiome_variable_info$variable_id,]
  
  dim(nasal_microbiome_sample_info)
  dim(nasal_microbiome_variable_info)
  
  rownames(nasal_microbiome_expression_data) == nasal_microbiome_variable_info$variable_id
  colnames(nasal_microbiome_expression_data) == nasal_microbiome_sample_info$sample_id
  
  remian_subject_id =
    nasal_microbiome_sample_info %>%
    dplyr::group_by(subject_id) %>%
    dplyr::summarise(n = n()) %>%
    dplyr::ungroup() %>%
    dplyr::filter(n > 5) %>%
    dplyr::pull(subject_id)
  
  nasal_microbiome_sample_info =
    nasal_microbiome_sample_info %>%
    dplyr::filter(subject_id %in% remian_subject_id)
  
  nasal_microbiome_expression_data =
    nasal_microbiome_expression_data[,nasal_microbiome_sample_info$sample_id]
  
}

###stool microbiome
{
  load(here::here("data_analysis/oral_microbiome/data_preparation/expression_data"))
  load(here::here("data_analysis/oral_microbiome/data_preparation/sample_info"))
  load(here::here("data_analysis/oral_microbiome/data_preparation/variable_info"))
  
  oral_microbiome_expression_data = expression_data
  oral_microbiome_sample_info = sample_info
  oral_microbiome_variable_info = variable_info
  
  ###read genus table
  expression_data =
    data.table::fread(here::here("data/from_xin/Genus Table/OR/Genus_OR.csv")) %>% 
    as.data.frame() %>% 
    tibble::column_to_rownames(var = "SampleID") %>%
    dplyr::select(-c(V1:SubjectID)) %>%
    t() %>%
    as.data.frame()
  
  oral_microbiome_variable_info =
    oral_microbiome_variable_info[match(rownames(expression_data), oral_microbiome_variable_info$Genus),]
  
  oral_microbiome_variable_info$Genus == rownames(expression_data)
  
  ###remove the variables which Genus are NA
  remove_idx = which(is.na(oral_microbiome_variable_info$Genus))
  remove_idx
  if(length(remove_idx) > 0){
    oral_microbiome_variable_info = oral_microbiome_variable_info[-remove_idx,]
    expression_data = expression_data[-remove_idx,]
  }
  
  rownames(expression_data) = oral_microbiome_variable_info$variable_id
  
  oral_microbiome_expression_data =
    expression_data
  
  oral_microbiome_variable_info = 
    oral_microbiome_variable_info %>% 
    dplyr::filter(!stringr::str_detect(Genus, "Unclassified_Bacteria"))
  
  oral_microbiome_expression_data = 
    oral_microbiome_expression_data[oral_microbiome_variable_info$variable_id,]
  
  dim(oral_microbiome_sample_info)
  dim(oral_microbiome_variable_info)
  
  rownames(oral_microbiome_expression_data) == oral_microbiome_variable_info$variable_id
  colnames(oral_microbiome_expression_data) == oral_microbiome_sample_info$sample_id
  
  remian_subject_id =
    oral_microbiome_sample_info %>%
    dplyr::group_by(subject_id) %>%
    dplyr::summarise(n = n()) %>%
    dplyr::ungroup() %>%
    dplyr::filter(n > 5) %>%
    dplyr::pull(subject_id)
  
  oral_microbiome_sample_info =
    oral_microbiome_sample_info %>%
    dplyr::filter(subject_id %in% remian_subject_id)
  
  oral_microbiome_expression_data =
    oral_microbiome_expression_data[,oral_microbiome_sample_info$sample_id]
  
}

###plasma lipidome
{
  load(here::here("data_analysis/lipidome/data_preparation/new_expression_data"))
  load(here::here("data_analysis/lipidome/data_preparation/new_sample_info"))
  load(here::here("data_analysis/lipidome/data_preparation/new_variable_info"))  
  
  lipidome_expression_data = new_expression_data
  lipidome_sample_info = new_sample_info
  lipidome_variable_info = new_variable_info
  
  ##remove QC samples
  lipidome_sample_info =
    lipidome_sample_info %>%
    dplyr::filter(subject_id != "QC")
  
  lipidome_expression_data = 
    lipidome_expression_data[,lipidome_sample_info$sample_id]
  
  dim(lipidome_expression_data)
  length(unique(lipidome_sample_info$subject_id))
  
  ###match samples
  dim(skin_microbiome_sample_info)
  dim(lipidome_sample_info)
  
  length(skin_microbiome_sample_info$subject_id)
  length(unique(skin_microbiome_sample_info$subject_id))
  
  ###just matched samples according to sample id, only 1 missed
  intersect_sample_id =
    intersect(skin_microbiome_sample_info$sample_id,
              lipidome_sample_info$sample_id)
  
  length(intersect_sample_id)
  
  skin_microbiome_expression_data =
    skin_microbiome_expression_data[,intersect_sample_id]
  
  lipidome_expression_data =
    lipidome_expression_data[,intersect_sample_id]
  
  skin_microbiome_sample_info =
    skin_microbiome_sample_info[match(intersect_sample_id, skin_microbiome_sample_info$sample_id),]
  
  lipidome_sample_info =
    lipidome_sample_info[match(intersect_sample_id, lipidome_sample_info$sample_id),]
  
  length(unique(lipidome_sample_info$subject_id))
  
  sum(lipidome_sample_info$sample_id == skin_microbiome_sample_info$sample_id)
  sum(lipidome_sample_info$subject_id == skin_microbiome_sample_info$subject_id)
  
  
  
  remian_subject_id =
    lipidome_sample_info %>%
    dplyr::group_by(subject_id) %>%
    dplyr::summarise(n = n()) %>%
    dplyr::ungroup() %>%
    dplyr::filter(n > 5) %>%
    dplyr::pull(subject_id)
  
  lipidome_sample_info =
    lipidome_sample_info %>%
    dplyr::filter(subject_id %in% remian_subject_id)
  
  lipidome_expression_data =
    lipidome_expression_data[,lipidome_sample_info$sample_id]
  
}

###plasma metabolome
{
  load(here::here(
    "data_analysis/metabolome/data_preparation/expression_data"
  ))
  load(here::here("data_analysis/metabolome/data_preparation/sample_info"))
  load(here::here(
    "data_analysis/metabolome/data_preparation/variable_info"
  ))
  
  metabolome_expression_data = expression_data
  metabolome_sample_info = sample_info
  metabolome_variable_info = variable_info
  
  metabolome_sample_info$CollectionDate =
    as.Date(metabolome_sample_info$CollectionDate, "%m/%d/%y")
  
  remian_subject_id =
    metabolome_sample_info %>%
    dplyr::group_by(subject_id) %>%
    dplyr::summarise(n = n()) %>%
    dplyr::ungroup() %>%
    dplyr::filter(n > 5) %>%
    dplyr::pull(subject_id)
  
  metabolome_sample_info =
    metabolome_sample_info %>%
    dplyr::filter(subject_id %in% remian_subject_id)
  
  metabolome_expression_data =
    metabolome_expression_data[,metabolome_sample_info$sample_id]
  
}


###plasma proteome
{
  load(here::here(
    "data_analysis/proteome/data_preparation/expression_data"
  ))
  load(here::here("data_analysis/proteome/data_preparation/sample_info"))
  load(here::here(
    "data_analysis/proteome/data_preparation/variable_info"
  ))
  
  proteome_expression_data = expression_data
  proteome_sample_info = sample_info
  proteome_variable_info = variable_info
  
  proteome_sample_info$CollectionDate =
    as.Date(proteome_sample_info$CollectionDate, "%m/%d/%y")
  
  dim(proteome_expression_data)
  length(unique(proteome_sample_info$subject_id))
  
  
  remian_subject_id =
    proteome_sample_info %>%
    dplyr::group_by(subject_id) %>%
    dplyr::summarise(n = n()) %>%
    dplyr::ungroup() %>%
    dplyr::filter(n > 5) %>%
    dplyr::pull(subject_id)
  
  proteome_sample_info =
    proteome_sample_info %>%
    dplyr::filter(subject_id %in% remian_subject_id)
  
  proteome_expression_data =
    proteome_expression_data[,proteome_sample_info$sample_id]
  
}

###get the correlation data
###stool microbiome
{
  stool_microbiome_metabolome_edge =
    readxl::read_xlsx(
      here::here(
        "data_analysis/correlation_network/whole_data_set/stool_microbiome_vs_metabolome/stool_microbiome_vs_metabolome.xlsx"
      ),
      sheet = 2
    )
  stool_microbiome_metabolome_edge$from = paste("stool", stool_microbiome_metabolome_edge$from, sep = "_")
  
  stool_microbiome_metabolome_node =
    readxl::read_xlsx(
      here::here(
        "data_analysis/correlation_network/whole_data_set/stool_microbiome_vs_metabolome/stool_microbiome_vs_metabolome.xlsx"
      ),
      sheet = 1
    )
  
  stool_microbiome_metabolome_node$node[grep("microbiome", stool_microbiome_metabolome_node$class)] =
    paste("stool", stool_microbiome_metabolome_node$node[grep("microbiome", stool_microbiome_metabolome_node$class)],
          sep = "_")
  
  
  stool_microbiome_lipidome_edge =
    readxl::read_xlsx(
      here::here(
        "data_analysis/correlation_network/whole_data_set/stool_microbiome_vs_lipidome/stool_microbiome_vs_lipidome.xlsx"
      ),
      sheet = 2
    )
  stool_microbiome_lipidome_edge$from = paste("stool", stool_microbiome_lipidome_edge$from, sep = "_")
  
  stool_microbiome_lipidome_node =
    readxl::read_xlsx(
      here::here(
        "data_analysis/correlation_network/whole_data_set/stool_microbiome_vs_lipidome/stool_microbiome_vs_lipidome.xlsx"
      ),
      sheet = 1
    )
  
  stool_microbiome_lipidome_node$node[grep("microbiome", stool_microbiome_lipidome_node$class)] =
    paste("stool", stool_microbiome_lipidome_node$node[grep("microbiome", stool_microbiome_lipidome_node$class)],
          sep = "_")
  
  
  stool_microbiome_proteome_edge =
    readxl::read_xlsx(
      here::here(
        "data_analysis/correlation_network/whole_data_set/stool_microbiome_vs_proteome/stool_microbiome_vs_proteome.xlsx"
      ),
      sheet = 2
    )
  stool_microbiome_proteome_edge$from = paste("stool", stool_microbiome_proteome_edge$from, sep = "_")
  
  stool_microbiome_proteome_node =
    readxl::read_xlsx(
      here::here(
        "data_analysis/correlation_network/whole_data_set/stool_microbiome_vs_proteome/stool_microbiome_vs_proteome.xlsx"
      ),
      sheet = 1
    )
  
  stool_microbiome_proteome_node$node[grep("microbiome", stool_microbiome_proteome_node$class)] =
    paste("stool", stool_microbiome_proteome_node$node[grep("microbiome", stool_microbiome_proteome_node$class)],
          sep = "_")
  
}


###skin microbiome
{
  skin_microbiome_metabolome_edge =
    readxl::read_xlsx(
      here::here(
        "data_analysis/correlation_network/whole_data_set/skin_microbiome_vs_metabolome/skin_microbiome_vs_metabolome.xlsx"
      ),
      sheet = 2
    )
  skin_microbiome_metabolome_edge$from = paste("skin", skin_microbiome_metabolome_edge$from, sep = "_")
  
  skin_microbiome_metabolome_node =
    readxl::read_xlsx(
      here::here(
        "data_analysis/correlation_network/whole_data_set/skin_microbiome_vs_metabolome/skin_microbiome_vs_metabolome.xlsx"
      ),
      sheet = 1
    )
  
  skin_microbiome_metabolome_node$node[grep("microbiome", skin_microbiome_metabolome_node$class)] =
    paste("skin", skin_microbiome_metabolome_node$node[grep("microbiome", skin_microbiome_metabolome_node$class)],
          sep = "_")
  
  
  skin_microbiome_lipidome_edge =
    readxl::read_xlsx(
      here::here(
        "data_analysis/correlation_network/whole_data_set/skin_microbiome_vs_lipidome/skin_microbiome_vs_lipidome.xlsx"
      ),
      sheet = 2
    )
  skin_microbiome_lipidome_edge$from = paste("skin", skin_microbiome_lipidome_edge$from, sep = "_")
  
  skin_microbiome_lipidome_node =
    readxl::read_xlsx(
      here::here(
        "data_analysis/correlation_network/whole_data_set/skin_microbiome_vs_lipidome/skin_microbiome_vs_lipidome.xlsx"
      ),
      sheet = 1
    )
  
  skin_microbiome_lipidome_node$node[grep("microbiome", skin_microbiome_lipidome_node$class)] =
    paste("skin", skin_microbiome_lipidome_node$node[grep("microbiome", skin_microbiome_lipidome_node$class)],
          sep = "_")
  
  skin_microbiome_proteome_edge =
    readxl::read_xlsx(
      here::here(
        "data_analysis/correlation_network/whole_data_set/skin_microbiome_vs_proteome/skin_microbiome_vs_proteome.xlsx"
      ),
      sheet = 2
    )
  skin_microbiome_proteome_edge$from = paste("skin", skin_microbiome_proteome_edge$from, sep = "_")
  
  skin_microbiome_proteome_node =
    readxl::read_xlsx(
      here::here(
        "data_analysis/correlation_network/whole_data_set/skin_microbiome_vs_proteome/skin_microbiome_vs_proteome.xlsx"
      ),
      sheet = 1
    )
  
  skin_microbiome_proteome_node$node[grep("microbiome", skin_microbiome_proteome_node$class)] =
    paste("skin", skin_microbiome_proteome_node$node[grep("microbiome", skin_microbiome_proteome_node$class)],
          sep = "_")
  
  # skin_microbiome_proteome_node$class = "Skin microbiome"
}


###nasal microbiome
{
  nasal_microbiome_metabolome_edge =
    readxl::read_xlsx(
      here::here(
        "data_analysis/correlation_network/whole_data_set/nasal_microbiome_vs_metabolome/nasal_microbiome_vs_metabolome.xlsx"
      ),
      sheet = 2
    )
  nasal_microbiome_metabolome_edge$from = paste("nasal", nasal_microbiome_metabolome_edge$from, sep = "_")
  
  nasal_microbiome_metabolome_node =
    readxl::read_xlsx(
      here::here(
        "data_analysis/correlation_network/whole_data_set/nasal_microbiome_vs_metabolome/nasal_microbiome_vs_metabolome.xlsx"
      ),
      sheet = 1
    )
  
  nasal_microbiome_metabolome_node$node[grep("microbiome", nasal_microbiome_metabolome_node$class)] =
    paste("nasal", nasal_microbiome_metabolome_node$node[grep("microbiome", nasal_microbiome_metabolome_node$class)],
          sep = "_")
  
  
  nasal_microbiome_lipidome_edge =
    readxl::read_xlsx(
      here::here(
        "data_analysis/correlation_network/whole_data_set/nasal_microbiome_vs_lipidome/nasal_microbiome_vs_lipidome.xlsx"
      ),
      sheet = 2
    )
  nasal_microbiome_lipidome_edge$from = paste("nasal", nasal_microbiome_lipidome_edge$from, sep = "_")
  
  nasal_microbiome_lipidome_node =
    readxl::read_xlsx(
      here::here(
        "data_analysis/correlation_network/whole_data_set/nasal_microbiome_vs_lipidome/nasal_microbiome_vs_lipidome.xlsx"
      ),
      sheet = 1
    )
  
  nasal_microbiome_lipidome_node$node[grep("microbiome", nasal_microbiome_lipidome_node$class)] =
    paste("nasal", nasal_microbiome_lipidome_node$node[grep("microbiome", nasal_microbiome_lipidome_node$class)],
          sep = "_")
  
  
  nasal_microbiome_proteome_edge =
    readxl::read_xlsx(
      here::here(
        "data_analysis/correlation_network/whole_data_set/nasal_microbiome_vs_proteome/nasal_microbiome_vs_proteome.xlsx"
      ),
      sheet = 2
    )
  nasal_microbiome_proteome_edge$from = paste("nasal", nasal_microbiome_proteome_edge$from, sep = "_")
  
  nasal_microbiome_proteome_node =
    readxl::read_xlsx(
      here::here(
        "data_analysis/correlation_network/whole_data_set/nasal_microbiome_vs_proteome/nasal_microbiome_vs_proteome.xlsx"
      ),
      sheet = 1
    )
  
  nasal_microbiome_proteome_node$node[grep("microbiome", nasal_microbiome_proteome_node$class)] =
    paste("nasal", nasal_microbiome_proteome_node$node[grep("microbiome", nasal_microbiome_proteome_node$class)],
          sep = "_")
  
}

###oral microbiome
{
  oral_microbiome_metabolome_edge =
    readxl::read_xlsx(
      here::here(
        "data_analysis/correlation_network/whole_data_set/oral_microbiome_vs_metabolome/oral_microbiome_vs_metabolome.xlsx"
      ),
      sheet = 2
    )
  oral_microbiome_metabolome_edge$from = paste("oral", oral_microbiome_metabolome_edge$from, sep = "_")
  
  oral_microbiome_metabolome_node =
    readxl::read_xlsx(
      here::here(
        "data_analysis/correlation_network/whole_data_set/oral_microbiome_vs_metabolome/oral_microbiome_vs_metabolome.xlsx"
      ),
      sheet = 1
    )
  
  oral_microbiome_metabolome_node$node[grep("microbiome", oral_microbiome_metabolome_node$class)] =
    paste("oral", oral_microbiome_metabolome_node$node[grep("microbiome", oral_microbiome_metabolome_node$class)],
          sep = "_")
  
  
  oral_microbiome_lipidome_edge =
    readxl::read_xlsx(
      here::here(
        "data_analysis/correlation_network/whole_data_set/oral_microbiome_vs_lipidome/oral_microbiome_vs_lipidome.xlsx"
      ),
      sheet = 2
    )
  oral_microbiome_lipidome_edge$from = paste("oral", oral_microbiome_lipidome_edge$from, sep = "_")
  
  oral_microbiome_lipidome_node =
    readxl::read_xlsx(
      here::here(
        "data_analysis/correlation_network/whole_data_set/oral_microbiome_vs_lipidome/oral_microbiome_vs_lipidome.xlsx"
      ),
      sheet = 1
    )
  
  oral_microbiome_lipidome_node$node[grep("microbiome", oral_microbiome_lipidome_node$class)] =
    paste("oral", oral_microbiome_lipidome_node$node[grep("microbiome", oral_microbiome_lipidome_node$class)],
          sep = "_")
  
  
  
  
  oral_microbiome_proteome_edge =
    readxl::read_xlsx(
      here::here(
        "data_analysis/correlation_network/whole_data_set/oral_microbiome_vs_proteome/oral_microbiome_vs_proteome.xlsx"
      ),
      sheet = 2
    )
  oral_microbiome_proteome_edge$from = paste("oral", oral_microbiome_proteome_edge$from, sep = "_")
  
  oral_microbiome_proteome_node =
    readxl::read_xlsx(
      here::here(
        "data_analysis/correlation_network/whole_data_set/oral_microbiome_vs_proteome/oral_microbiome_vs_proteome.xlsx"
      ),
      sheet = 1
    )
  
  oral_microbiome_proteome_node$node[grep("microbiome", oral_microbiome_proteome_node$class)] =
    paste("oral", oral_microbiome_proteome_node$node[grep("microbiome", oral_microbiome_proteome_node$class)],
          sep = "_")
  
}

edge_data =
  rbind(
    stool_microbiome_metabolome_edge,
    stool_microbiome_lipidome_edge,
    stool_microbiome_proteome_edge,
    skin_microbiome_metabolome_edge,
    skin_microbiome_lipidome_edge,
    skin_microbiome_proteome_edge,
    nasal_microbiome_metabolome_edge,
    nasal_microbiome_lipidome_edge,
    nasal_microbiome_proteome_edge,
    oral_microbiome_metabolome_edge,
    oral_microbiome_lipidome_edge,
    oral_microbiome_proteome_edge
  ) %>% 
  dplyr::mutate(direction = 
                  case_when(cor > 0 ~ "positive",
                            cor < 0 ~ "negative")) %>% 
  dplyr::mutate(pos_neg = direction)

node_data = 
  rbind(
    stool_microbiome_metabolome_node,
    stool_microbiome_lipidome_node,
    stool_microbiome_proteome_node,
    skin_microbiome_metabolome_node,
    skin_microbiome_lipidome_node,
    skin_microbiome_proteome_node,
    nasal_microbiome_metabolome_node,
    nasal_microbiome_lipidome_node,
    nasal_microbiome_proteome_node,
    oral_microbiome_metabolome_node,
    oral_microbiome_lipidome_node,
    oral_microbiome_proteome_node
  ) %>% 
  dplyr::distinct(node, .keep_all = TRUE)

##remove Archaea
remove_name = 
  c("Zea",
    stool_microbiome_variable_info %>% 
      dplyr::filter(Kingdom == "Archaea") %>% 
      dplyr::pull(Genus),
    skin_microbiome_variable_info %>% 
      dplyr::filter(Kingdom == "Archaea") %>% 
      dplyr::pull(Genus),
    oral_microbiome_variable_info %>% 
      dplyr::filter(Kingdom == "Archaea") %>% 
      dplyr::pull(Genus),
    nasal_microbiome_variable_info %>% 
      dplyr::filter(Kingdom == "Archaea") %>% 
      dplyr::pull(Genus)    
  ) %>% 
  unique()

if(length(remove_name) > 0) {
  node_data = 
    node_data %>% 
    dplyr::filter(!true_name %in% remove_name)
}

node_data$node

library(ggraph)
library(igraph)
library(tidygraph)

edge_data =
  edge_data %>%
  dplyr::filter(from %in% node_data$node & to %in% node_data$node)

node_data =
  node_data %>%
  dplyr::filter(node %in% edge_data$from | node %in% edge_data$to)

dim(edge_data)
dim(node_data)

edge_data$edge_name = paste(edge_data$from, edge_data$to, sep = ";")

graph <- 
  tidygraph::tbl_graph(nodes = node_data, 
                       edges = edge_data,
                       directed = FALSE) %>% 
  dplyr::mutate(Degree = centrality_degree(mode = 'all'))

####plot some examples
dim(edge_data)
edge_data$from_true_name

idx = which(edge_data$to_true_name == "Alpha-N-Phenylacetyl-L-glutamine")
edge_data$from_true_name[idx]

edge_data[idx,]

for(i in idx){
  from = edge_data[i,]$from 
  to = edge_data[i,]$to
  
  from_type = node_data$class[node_data$node == from]
  to_type = node_data$class[node_data$node == to]
  
  from_data =
    switch(
      EXPR = from_type,
      "Stool microbiome" = stool_microbiome_expression_data,
      "Skin microbiome" = skin_microbiome_expression_data,
      "Oral microbiome" = oral_microbiome_expression_data,
      "Nasal microbiome" = nasal_microbiome_expression_data,
      "Metabolite" = metabolome_expression_data,
      "Lipidome" = lipidome_expression_data,
      "Proteome" = proteome_expression_data
    )
  
  from_sample_info =
    switch(
      EXPR = from_type,
      "Stool microbiome" = stool_microbiome_sample_info,
      "Skin microbiome" = skin_microbiome_sample_info,
      "Oral microbiome" = oral_microbiome_sample_info,
      "Nasal microbiome" = nasal_microbiome_sample_info,
      "Metabolite" = metabolome_sample_info,
      "Lipidome" = lipidome_sample_info,
      "Proteome" = proteome_sample_info
    )
  
  to_data =
    switch(
      EXPR = to_type,
      "Stool microbiome" = stool_microbiome_expression_data,
      "Skin microbiome" = skin_microbiome_expression_data,
      "Oral microbiome" = oral_microbiome_expression_data,
      "Nasal microbiome" = nasal_microbiome_expression_data,
      "Metabolite" = metabolome_expression_data,
      "Lipidome" = lipidome_expression_data,
      "Proteome" = proteome_expression_data
    )
  
  to_sample_info =
    switch(
      EXPR = to_type,
      "Stool microbiome" = stool_microbiome_sample_info,
      "Skin microbiome" = skin_microbiome_sample_info,
      "Oral microbiome" = oral_microbiome_sample_info,
      "Nasal microbiome" = nasal_microbiome_sample_info,
      "Metabolite" = metabolome_sample_info,
      "Lipidome" = lipidome_sample_info,
      "Proteome" = proteome_sample_info
    )
  
  ###match samples
  ###just matched samples according to sample id, only 1 missed
  intersect_sample_id =
    intersect(from_sample_info$sample_id,
              to_sample_info$sample_id)

  from_data =
    from_data[,intersect_sample_id]

  to_data =
    to_data[,intersect_sample_id]

  from_sample_info =
    from_sample_info[match(intersect_sample_id, from_sample_info$sample_id),]

  to_sample_info =
    to_sample_info[match(intersect_sample_id, to_sample_info$sample_id),]

  sum(from_sample_info$sample_id == to_sample_info$sample_id)
  sum(from_sample_info$subject_id == to_sample_info$subject_id)
  

  dim(from_data)
  dim(to_data)

  from =
    from %>%
    stringr::str_replace_all(., "stool_", "") %>%
    stringr::str_replace_all(., "skin_", "") %>%
    stringr::str_replace_all(., "oral_", "") %>%
    stringr::str_replace_all(., "nasal_", "")

  plot(as.numeric(from_data[from,]),
       as.numeric(to_data[to,]))
}




