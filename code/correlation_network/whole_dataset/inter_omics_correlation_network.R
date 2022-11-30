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
setwd("data_analysis/correlation_network/whole_data_set/inter_omics_correlation_network")

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


###output plot
#####network to show the correlation between microbiome and metabolite
dim(stool_microbiome_metabolome_edge)
dim(stool_microbiome_lipidome_edge)

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

###example network
library(plyr)

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

###output result
microbiome_info =
  edge_data %>%
  plyr::dlply(.variables = .(from)) %>%
  purrr::map(function(z) {
    x =
      z %>%
      dplyr::select(to) %>%
      dplyr::left_join(node_data, by = c("to" = "node")) %>%
      dplyr::pull(class)
    
    data.frame(
      microbiome_id = z$from[1],
      microbiome_name = z$from_true_name[1],
      microbiome_class = node_data$class[match(z$from[1],node_data$node)],
      total_number = length(x),
      Proteome = sum(x == "Proteome") * 100 / length(x),
      Metabolite = sum(x == "Metabolite") * 100 / length(x),
      Lipidome = sum(x == "Lipidome") * 100 / length(x)
    )
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame() %>% 
  dplyr::arrange(microbiome_class, desc(total_number)) %>% 
  dplyr::mutate(direction = case_when(
    Proteome >= 50 ~ "Proteome",
    Metabolite >= 50 ~ "Metabolite",
    Lipidome >= 50 ~ "Lipidome",
    TRUE ~ "Mixed"
  ))

microbiome_info %>%
  plyr::dlply(.variables = .(microbiome_class)) %>%
  purrr::map(function(x) {
    table(x$direction) * 100 / nrow(x)
  })

molecular_info =
  edge_data %>%
  plyr::dlply(.variables = .(to)) %>%
  purrr::map(function(z) {
    x =
      z %>%
      dplyr::select(from) %>%
      dplyr::left_join(node_data, by = c("from" = "node")) %>%
      dplyr::pull(class)
    
    data.frame(
      molecular_id = z$to[1],
      molecular_name = z$to_true_name[1],
      molecular_class = node_data$class[match(z$to[1],node_data$node)],
      total_number = length(x),
      Stool = sum(x == "Stool microbiome") * 100 / length(x),
      Skin = sum(x == "Skin microbiome") * 100 / length(x),
      Nasal = sum(x == "Nasal microbiome") * 100 / length(x),
      Oral = sum(x == "Oral microbiome") * 100 / length(x)
    )
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame() %>% 
  dplyr::arrange(molecular_class, desc(total_number)) %>% 
  dplyr::mutate(direction = case_when(
    Stool >= 50 ~ "Stool",
    Skin >= 50 ~ "Skin",
    Nasal >= 50 ~ "Nasal",
    Oral >= 50 ~ "Oral",
    TRUE ~ "Mixed"
  ))

molecular_info %>%
  plyr::dlply(.variables = .(molecular_class)) %>%
  purrr::map(function(x) {
    table(x$direction) * 100 / nrow(x)
  })

microbiome_info

####subnetwork according to microbiome
#####for each microbiome, give it a class according to connection to which omics
node = igraph::vertex_attr(graph)$node
module = 
data.frame(node = node) %>% 
  dplyr::left_join(microbiome_info[,c("microbiome_id", "total_number","direction")],
                   by = c('node' = "microbiome_id"))

module$direction[is.na(module$direction)] = "Other"
module$direction[module$direction == "Mixed"] = "Other"
module$direction[which(module$total_number < 6)] = "Other"

module = module$direction

graph =
  graph %>%
  tidygraph::activate(what = "nodes") %>%
  dplyr::mutate(module = module)

# subnetworks <-
#   igraph::cluster_fast_greedy(graph = graph,
#                                    weights = abs(edge_attr(graph,
#                                                            "cor")))
# 
# plot =
#   modularity_plot(subnetworks = subnetworks)
# 
# plot
# 
# table(igraph::membership(communities = subnetworks))
# 
# ggsave(
#   plot,
#   filename = "all_modularity.pdf",
#   width = 9,
#   height = 7
# )
# 
# membership =
#   igraph::cut_at(communities = subnetworks, no = 5)
# 
# table(membership)
# 
# module = paste("Module", membership, sep = "_")
#     
# remove_module_name =
#     data.frame(module) %>%
#       dplyr::group_by(module) %>%
#       dplyr::summarise(n = n()) %>%
#       dplyr::filter(n < 10) %>%
#       dplyr::pull(module)
# 
# module[module %in% remove_module_name] = "Other"
# table(module)
# 
# graph =
#   graph %>%
#   tidygraph::activate(what = "nodes") %>%
#   dplyr::mutate(module = module)

angle <- 360 * (c(1:nrow(node_data)) - 0.5)/nrow(node_data)
hjust <- ifelse(angle > 180, 1, 0)
angle <- ifelse(angle > 180, 90 - angle + 180, 90 - angle)

node_data =
  igraph::vertex_attr(graph) %>% 
  dplyr::bind_rows()

label_node =
  node_data %>%
  plyr::dlply(.variables = .(class)) %>%
  purrr::map(function(x) {
    x %>%
      dplyr::arrange(desc(Degree)) %>%
      # dplyr::filter(Degree > 10) %>%
      head(5) %>%
      pull(node)
  }) %>%
  unlist() %>%
  unname()

label_node_stool = 
node_data %>%
  dplyr::filter(class == "Stool microbiome") %>%
  dplyr::filter(
    true_name %in% c(
      "Bacteroides",
      "Prevotella",
      "Phocaeicola",
      "Unclassified_Ruminococcaceae"
    )
  ) %>% 
  dplyr::pull(node)

label_node_skin = 
  node_data %>%
  dplyr::filter(class == "Skin microbiome") %>%
  dplyr::filter(
    true_name %in% c(
      "Corynebacterium",
      "Staphylococcus",
      "Cutibacterium",
      "Anaerococcus"
    )
  ) %>% 
  dplyr::pull(node)

label_node_nasal = 
  node_data %>%
  dplyr::filter(class == "Nasal microbiome") %>%
  dplyr::filter(
    true_name %in% c(
      "Corynebacterium",
      "Staphylococcus",
      "Cutibacterium",
      "Pertoniphilus"
    )
  ) %>% 
  dplyr::pull(node)

label_node_oral = 
  node_data %>%
  dplyr::filter(class == "Oral microbiome") %>%
  dplyr::filter(
    true_name %in% c(
      "Prevotella",
      "Neisseria",
      "Haemophilus",
      "Streptococcus"
    )
  ) %>% 
  dplyr::pull(node)

label_node =
  c(label_node,
    label_node_skin,
    label_node_stool,
    label_node_oral,
    label_node_nasal) %>% 
  unique()
    
g = graph

coords <-
  create_layout(g, layout = "fr") %>%
  dplyr::select(node, x, y)

my_graph <-
  create_layout(
    graph = g,
    layout = "manual",
    x = coords$x,
    y = coords$y
  )

plot <-
  ggraph(graph,
         layout = "manual", 
         x = coords$x,
         y = coords$y,
         circular = FALSE) +
  ggforce::geom_mark_hull(
    aes(
      x = x,
      y = y,
      group = module,
      fill = module,
      filter = module != 'Other'
    ),
    concavity = 200,
    show.legend = FALSE,
    alpha = 0.1
  ) +
  geom_edge_link(
    aes(
      color = pos_neg,
      width = -log(p_adjust, 10)
      # linetype = significance
    ),
    alpha = 1,
    show.legend = TRUE
  ) +
  scale_edge_linetype_manual(values = c("p.adj<0.05" = 1,
                                        "0.05<p.adj<0.2" = 2)) +
  ggstar::geom_star(aes(
    x = x,
    y = y,
    size = Degree,
    fill = class,
    starshape = class
  )) +
  ggstar::scale_starshape_manual(values = omics_shape) +
  # geom_node_text(aes(x = x ,
  #                    y = y ,
  #                    color = "black",
  #                    label = true_name),
  #                repel = TRUE,
  #                size = 2) +
  shadowtext::geom_shadowtext(aes(x = x, y = y,
                                  label = ifelse(node %in% label_node, true_name, NA),
                                  color = class), 
                              size = 5,
                              check_overlap = TRUE,
                              bg.colour = "white",
                              show.legend = FALSE) +
  guides(edge_color = guide_edge_colorbar(title = "Correlation"),
         color = guide_legend(title = "Class")) +
  scale_edge_width_continuous(range = c(0.1, 1)) +
  scale_size_continuous(range = c(2, 8)) +
  scale_fill_manual(values = omics_color) +
  scale_color_manual(values = omics_color) +
  scale_edge_color_manual(values = c("positive" = viridis::inferno(n = 10)[5],
                                     "negative" = viridis::inferno(n = 10)[2])) +
  ggraph::theme_graph() +
  theme(plot.background = element_rect(fill = "transparent", color = NA), 
        panel.background = element_rect(fill = "transparent", color = NA),
        legend.position = "right",
        legend.background = element_rect(fill = "transparent", color = NA))

plot

# # extrafont::loadfonts()
# ggsave(plot,
#        filename = "microbiome_molecular_cor_plot.pdf",
#        width = 17,
#        height = 14)

###only highlight stool microbiome
idx = 
edge_data %>% 
  dplyr::left_join(node_data[,c("node", "class")],
                   by = c("from" = "node")) %>% 
  dplyr::rename(from_class = class) %>% 
  dplyr::pull(from_class) %>% 
  `==`("Stool microbiome") %>% 
  which()

hightlight_node =
  c(unique(edge_data$from[idx]),
    unique(edge_data$to[idx]))

hightlight_node = 
node_data[,c("node")] %>% 
  dplyr::mutate(highlight = 
                  case_when(node %in% hightlight_node ~ "yes",
                            TRUE ~ 'no'))

temp_graph  =
  graph %>% 
  tidygraph::activate(what = "nodes") %>% 
  tidygraph::left_join(hightlight_node, by  ="node") %>% 
  tidygraph::mutate(new_class = case_when(
    highlight == "yes" ~ class,
    highlight != "yes" ~ "no",
  )) 

temp_graph =
  temp_graph %>%
  tidygraph::activate(what = "edges") %>%
  tidygraph::mutate(from_id = stringr::str_split(edge_name, ";") %>% lapply(function(x)x[1]) %>% unlist(),
                    to_id = stringr::str_split(edge_name, ";") %>% lapply(function(x)x[2]) %>% unlist()) %>% 
  tidygraph::mutate(
    new_edge_class =
      case_when(
        from_id %in% hightlight_node$node[hightlight_node$highlight == "no"] ~ "no",
        to_id %in% hightlight_node$node[hightlight_node$highlight == "no"] ~ "no",
        TRUE ~ "yes"
      )
  )

new_class_alpha = c(
  "Metabolite" = 1,
  "Lipidome" = 1,
  "Proteome" = 1,
  "Stool microbiome" = 1,
  "Skin microbiome" = 1,
  "Oral microbiome" = 1,
  "Nasal microbiome" = 1,
  "no" = 0.3
)



plot <-
  ggraph(temp_graph,
         layout = "manual", 
         x = coords$x,
         y = coords$y,
         circular = FALSE) +
  geom_edge_link(
    aes(alpha = new_edge_class),
    color = "black",
    show.legend = TRUE
  ) +
  scale_edge_alpha_manual(values = c("yes" = 1, "no" = 0.05)) +
  ggstar::geom_star(aes(
    x = x,
    y = y,
    size = Degree,
    fill = new_class,
    alpha = new_class,
    starshape = class
  ), 
  color = NA,
  show.legend = TRUE) +
  ggstar::scale_starshape_manual(values = omics_shape) +
  guides(edge_color = guide_edge_colorbar(title = "Correlation"),
         color = guide_legend(title = "Class")) +
  scale_size_continuous(range = c(2, 8)) +
  scale_fill_manual(values = c(omics_color, no = "grey")) +
  scale_alpha_manual(values = new_class_alpha) +
  ggraph::theme_graph() +
  theme(plot.background = element_rect(fill = "transparent", color = NA), 
        panel.background = element_rect(fill = "transparent", color = NA),
        legend.position = "right",
        legend.background = element_rect(fill = "transparent", color = NA))

plot

# extrafont::loadfonts()
ggsave(plot,
       filename = "highlight_stool_microbiome_molecular_cor_plot.pdf",
       width = 17,
       height = 14)







###only highlight skin microbiome
idx = 
  edge_data %>% 
  dplyr::left_join(node_data[,c("node", "class")],
                   by = c("from" = "node")) %>% 
  dplyr::rename(from_class = class) %>% 
  dplyr::pull(from_class) %>% 
  `==`("Skin microbiome") %>% 
  which()

hightlight_node =
  c(unique(edge_data$from[idx]),
    unique(edge_data$to[idx]))

hightlight_node = 
  node_data[,c("node")] %>% 
  dplyr::mutate(highlight = 
                  case_when(node %in% hightlight_node ~ "yes",
                            TRUE ~ 'no'))

temp_graph  =
  graph %>% 
  tidygraph::activate(what = "nodes") %>% 
  tidygraph::left_join(hightlight_node, by  ="node") %>% 
  tidygraph::mutate(new_class = case_when(
    highlight == "yes" ~ class,
    highlight != "yes" ~ "no",
  )) 

temp_graph =
  temp_graph %>%
  tidygraph::activate(what = "edges") %>%
  tidygraph::mutate(from_id = stringr::str_split(edge_name, ";") %>% lapply(function(x)x[1]) %>% unlist(),
                    to_id = stringr::str_split(edge_name, ";") %>% lapply(function(x)x[2]) %>% unlist()) %>% 
  tidygraph::mutate(
    new_edge_class =
      case_when(
        from_id %in% hightlight_node$node[hightlight_node$highlight == "no"] ~ "no",
        to_id %in% hightlight_node$node[hightlight_node$highlight == "no"] ~ "no",
        TRUE ~ "yes"
      )
  )

plot <-
  ggraph(temp_graph,
         layout = "manual", 
         x = coords$x,
         y = coords$y,
         circular = FALSE) +
  geom_edge_link(
    aes(alpha = new_edge_class),
    color = "black",
    show.legend = TRUE
  ) +
  scale_edge_alpha_manual(values = c("yes" = 1, "no" = 0.05)) +
  ggstar::geom_star(aes(
    x = x,
    y = y,
    size = Degree,
    fill = new_class,
    alpha = new_class,
    starshape = class
  ), 
  color = NA,
  show.legend = TRUE) +
  ggstar::scale_starshape_manual(values = omics_shape) +
  guides(edge_color = guide_edge_colorbar(title = "Correlation"),
         color = guide_legend(title = "Class")) +
  scale_size_continuous(range = c(2, 8)) +
  scale_fill_manual(values = c(omics_color, no = "grey")) +
  scale_alpha_manual(values = new_class_alpha) +
  ggraph::theme_graph() +
  theme(plot.background = element_rect(fill = "transparent", color = NA), 
        panel.background = element_rect(fill = "transparent", color = NA),
        legend.position = "right",
        legend.background = element_rect(fill = "transparent", color = NA))

plot

# extrafont::loadfonts()
ggsave(plot,
       filename = "highlight_skin_microbiome_molecular_cor_plot.pdf",
       width = 17,
       height = 14)

###only highlight oral microbiome
idx = 
  edge_data %>% 
  dplyr::left_join(node_data[,c("node", "class")],
                   by = c("from" = "node")) %>% 
  dplyr::rename(from_class = class) %>% 
  dplyr::pull(from_class) %>% 
  `==`("Oral microbiome") %>% 
  which()

hightlight_node =
  c(unique(edge_data$from[idx]),
    unique(edge_data$to[idx]))

hightlight_node = 
  node_data[,c("node")] %>% 
  dplyr::mutate(highlight = 
                  case_when(node %in% hightlight_node ~ "yes",
                            TRUE ~ 'no'))

temp_graph  =
  graph %>% 
  tidygraph::activate(what = "nodes") %>% 
  tidygraph::left_join(hightlight_node, by  ="node") %>% 
  tidygraph::mutate(new_class = case_when(
    highlight == "yes" ~ class,
    highlight != "yes" ~ "no",
  )) 

temp_graph =
  temp_graph %>%
  tidygraph::activate(what = "edges") %>%
  tidygraph::mutate(from_id = stringr::str_split(edge_name, ";") %>% lapply(function(x)x[1]) %>% unlist(),
                    to_id = stringr::str_split(edge_name, ";") %>% lapply(function(x)x[2]) %>% unlist()) %>% 
  tidygraph::mutate(
    new_edge_class =
      case_when(
        from_id %in% hightlight_node$node[hightlight_node$highlight == "no"] ~ "no",
        to_id %in% hightlight_node$node[hightlight_node$highlight == "no"] ~ "no",
        TRUE ~ "yes"
      )
  )

plot <-
  ggraph(temp_graph,
         layout = "manual", 
         x = coords$x,
         y = coords$y,
         circular = FALSE) +
  geom_edge_link(
    aes(alpha = new_edge_class),
    color = "black",
    show.legend = TRUE
  ) +
  scale_edge_alpha_manual(values = c("yes" = 1, "no" = 0.05)) +
  ggstar::geom_star(aes(
    x = x,
    y = y,
    size = Degree,
    fill = new_class,
    alpha = new_class,
    starshape = class
  ), 
  color = NA,
  show.legend = TRUE) +
  ggstar::scale_starshape_manual(values = omics_shape) +
  guides(edge_color = guide_edge_colorbar(title = "Correlation"),
         color = guide_legend(title = "Class")) +
  scale_size_continuous(range = c(2, 8)) +
  scale_fill_manual(values = c(omics_color, no = "grey")) +
  scale_alpha_manual(values = new_class_alpha) +
  ggraph::theme_graph() +
  theme(plot.background = element_rect(fill = "transparent", color = NA), 
        panel.background = element_rect(fill = "transparent", color = NA),
        legend.position = "right",
        legend.background = element_rect(fill = "transparent", color = NA))

plot

# extrafont::loadfonts()
ggsave(plot,
       filename = "highlight_oral_microbiome_molecular_cor_plot.pdf",
       width = 17,
       height = 14)






###only highlight nasal microbiome
idx = 
  edge_data %>% 
  dplyr::left_join(node_data[,c("node", "class")],
                   by = c("from" = "node")) %>% 
  dplyr::rename(from_class = class) %>% 
  dplyr::pull(from_class) %>% 
  `==`("Nasal microbiome") %>% 
  which()

hightlight_node =
  c(unique(edge_data$from[idx]),
    unique(edge_data$to[idx]))

hightlight_node = 
  node_data[,c("node")] %>% 
  dplyr::mutate(highlight = 
                  case_when(node %in% hightlight_node ~ "yes",
                            TRUE ~ 'no'))

temp_graph  =
  graph %>% 
  tidygraph::activate(what = "nodes") %>% 
  tidygraph::left_join(hightlight_node, by  ="node") %>% 
  tidygraph::mutate(new_class = case_when(
    highlight == "yes" ~ class,
    highlight != "yes" ~ "no",
  )) 

temp_graph =
  temp_graph %>%
  tidygraph::activate(what = "edges") %>%
  tidygraph::mutate(from_id = stringr::str_split(edge_name, ";") %>% lapply(function(x)x[1]) %>% unlist(),
                    to_id = stringr::str_split(edge_name, ";") %>% lapply(function(x)x[2]) %>% unlist()) %>% 
  tidygraph::mutate(
    new_edge_class =
      case_when(
        from_id %in% hightlight_node$node[hightlight_node$highlight == "no"] ~ "no",
        to_id %in% hightlight_node$node[hightlight_node$highlight == "no"] ~ "no",
        TRUE ~ "yes"
      )
  )

plot <-
  ggraph(temp_graph,
         layout = "manual", 
         x = coords$x,
         y = coords$y,
         circular = FALSE) +
  geom_edge_link(
    aes(alpha = new_edge_class),
    color = "black",
    show.legend = TRUE
  ) +
  scale_edge_alpha_manual(values = c("yes" = 1, "no" = 0.05)) +
  ggstar::geom_star(aes(
    x = x,
    y = y,
    size = Degree,
    fill = new_class,
    alpha = new_class,
    starshape = class
  ), 
  color = NA,
  show.legend = TRUE) +
  ggstar::scale_starshape_manual(values = omics_shape) +
  guides(edge_color = guide_edge_colorbar(title = "Correlation"),
         color = guide_legend(title = "Class")) +
  scale_size_continuous(range = c(2, 8)) +
  scale_fill_manual(values = c(omics_color, no = "grey")) +
  scale_alpha_manual(values = new_class_alpha) +
  ggraph::theme_graph() +
  theme(plot.background = element_rect(fill = "transparent", color = NA), 
        panel.background = element_rect(fill = "transparent", color = NA),
        legend.position = "right",
        legend.background = element_rect(fill = "transparent", color = NA))

plot

# extrafont::loadfonts()
ggsave(plot,
       filename = "highlight_nasal_microbiome_molecular_cor_plot.pdf",
       width = 17,
       height = 14)


######----------------------------------------------------------
molecular_info

dim(node_data)
dim(edge_data)

table(node_data$class)

edge_data %>%
  dplyr::select(from, to) %>%
  dplyr::left_join(node_data[, c("node", "class")], by = c("from" = "node")) %>%
  dplyr::rename(from_class = class) %>% 
  dplyr::left_join(node_data[, c("node", "class")], by = c("to" = "node")) %>%
  dplyr::rename(to_class = class) %>% 
  dplyr::group_by(from_class, to_class) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::arrange(from_class, desc(n))

microbiome_info %>% 
  dplyr::group_by(microbiome_class) %>% 
  dplyr::summarise(n = sum(total_number))

# edge distributation
library(ggsankey)
library(dplyr)
temp_data = 
  edge_data %>% 
  dplyr::left_join(node_data[,c("node", "class")], by = c("from" = "node")) %>% 
  dplyr::rename(from_class = class) %>% 
  dplyr::left_join(node_data[,c("node", "class")], by = c("to" = "node")) %>% 
  dplyr::rename(to_class = class) %>% 
  dplyr::select(Microbiome = from_class, Molecule = to_class) %>% 
  ggsankey::make_long(Molecule, Microbiome)

color = c(omics_color)

plot = 
  ggplot(
    temp_data,
    aes(
      x = x,
      next_x = next_x,
      node = node,
      next_node = next_node,
      fill = factor(node),
      label = node
    )
  ) +
  geom_sankey(flow.alpha = 0.6,
              node.color = "black", 
              width = 0.3) +
  geom_sankey_label(size = 4, 
                    aes(color = factor(node)),
                    # color = "white", 
                    fill = "white") +
  scale_fill_manual(values = color) +
  scale_color_manual(values = color) +
  theme_sankey(base_size = 10) +
  labs(x = NULL) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = .5),
        plot.background = element_blank(),
        panel.background = element_blank())

plot

# ggsave(
#   plot,
#   filename = "edge_distributation.pdf",
#   width = 12,
#   height = 7
# )

###node distributation
##mosaic plot
temp_data = 
node_data %>%
  dplyr::select(node, class) %>%
  dplyr::group_by(class) %>%
  dplyr::summarise(n = n())  %>%
  dplyr::ungroup() %>% 
  dplyr::mutate(class = factor(
    class,
    levels = c(
      "Stool microbiome",
      "Skin microbiome",
      "Oral microbiome",
      "Nasal microbiome",
      "Proteome",
      "Metabolite",
      "Lipidome"
    )
  ))

rownames(temp_data) = as.character(temp_data$class)

temp_data = 
temp_data[rev(levels(temp_data$class)),]

temp_data$y_position = 
purrr::map(1:nrow(temp_data), 
           function(idx){
             if(idx == 1){
               return(temp_data$n[idx]/2)}
             sum(temp_data$n[1:(idx-1)]) + temp_data$n[idx]/2
           }) %>% 
  unlist()

temp_data$label = paste(temp_data$n, "(", round(temp_data$n*100/sum(temp_data$n), 2), '%)',
                        sep = "")

plot = 
temp_data %>% 
  ggplot(aes(x = 1, y = n)) +
  geom_bar(aes(fill = class), stat = "identity",
           color = "white",
           show.legend = FALSE) +
  scale_fill_manual(values = omics_color) +
  geom_text(aes(x = 1, y = y_position, label = label),
            color = "white") +
  scale_x_continuous(limits = c(0, 1.5)) +
  coord_polar(theta = "y") +
  theme_void()

plot

# ggsave(plot, filename = "node_class.pdf", width = 9, height = 7)


#####explore the distributation 
edge_data %>%
  dplyr::left_join(node_data[, c("node", "class")], 
                   by = c("from" = "node")) %>%
  dplyr::rename(from_class = class) %>%
  dplyr::left_join(node_data[, c("node", "class")], 
                   by = c("to" = "node")) %>%
  dplyr::rename(to_class = class) %>%
  dplyr::filter(to_class == "Lipidome") %>%
  ggplot() +
  geom_density(mapping = aes(cor, 
                             color = from_class)) +
  scale_color_manual(values = omics_color[c("Stool microbiome",
                                            "Skin microbiome",
                                            "Oral microbiome",
                                            "Nasal microbiome")]) +
  facet_grid(cols = vars(direction), scales = "free_x") +
  base_theme +
  labs(x = "Correlation", y = "Density")

median_int = 
skin_microbiome_expression_data[match(node_data$node[grep("skin", node_data$node)],
                                      paste("skin", rownames(skin_microbiome_expression_data), sep = "_")),] %>% 
  apply(1, median)

temp_data = 
data.frame(
  node = node_data$node[grep("skin", node_data$node)],
  median_int = unname(median_int))

median_cor = 
temp_data$node %>%
  purrr::map(function(x) {
    edge_data %>%
      dplyr::filter(from == x) %>%
      dplyr::filter(to %in% lipidome_variable_info$variable_id) %>% 
      dplyr::pull(cor) %>% 
      median()
  }) %>% 
  unlist()

plot(median_int, median_cor)

#######lipidome
library(dplyr)
library(ggplot2)
library(ggmosaic)
library(ggmosaic)

df =
  edge_data %>%
  dplyr::left_join(node_data[, c("node", "class")],
                   by = c("from" = "node")) %>%
  dplyr::rename(from_class = class) %>%
  dplyr::left_join(node_data[, c("node", "class")],
                   by = c("to" = "node")) %>%
  dplyr::rename(to_class = class) %>%
  dplyr::filter(to_class == "Lipidome") %>%
  dplyr::mutate(from_class = factor(
    from_class,
    levels = c(
      "Skin microbiome",
      "Oral microbiome",
      "Stool microbiome",
      "Nasal microbiome"
    )
  ))

plot =
  ggplot(df) +
  geom_mosaic(aes(x = product(direction, from_class), fill = direction),
              show.legend = FALSE) +
  scale_fill_manual(values = c(
    "positive" = ggsci::pal_aaas()(n = 10)[2],
    "negative" = ggsci::pal_aaas()(n = 10)[1]
  )) +
  labs(x = "", y = "") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = "transparent"),
    panel.background = element_rect(fill = "transparent")
  )

plot

temp = 
df %>% 
  dplyr::group_by(from_class, direction) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(prop = round(n*100/sum(n), 2))

text_data =
  ggplot_build(plot)$data %>%
  as.data.frame %>%
  transmute(
    x__from_class,
    x__fill__direction,
    x.position = (xmax + xmin) / 2,
    y.position = (ymax + ymin) / 2
  ) %>%
  dplyr::rename(from_class = x__from_class,
                direction = x__fill__direction) %>%
  dplyr::left_join(temp, by = c("from_class", "direction"))

plot = 
plot +
  geom_text(aes(x.position, y.position, label = prop),
            data = text_data,
            color = "white")
plot
# ggsave(
#   plot = plot,
#   filename = "proteome_direction_distributation.pdf",
#   width = 9,
#   height = 7
# )



# #####output result
# library(openxlsx)
# wb <- createWorkbook()
# modifyBaseFont(wb, fontSize = 12, fontName = "Time New Roma")
# addWorksheet(wb, sheetName = "Node data", gridLines = TRUE)
# addWorksheet(wb, sheetName = "Edge data", gridLines = TRUE)
# addWorksheet(wb, sheetName = "Microbiome information", gridLines = TRUE)
# addWorksheet(wb, sheetName = "Molecular information", gridLines = TRUE)
# freezePane(wb, sheet = 1, firstRow = TRUE, firstCol = TRUE)
# freezePane(wb, sheet = 2, firstRow = TRUE, firstCol = TRUE)
# freezePane(wb, sheet = 3, firstRow = TRUE, firstCol = TRUE)
# freezePane(wb, sheet = 4, firstRow = TRUE, firstCol = TRUE)
# writeDataTable(wb, sheet = 1, x = as.data.frame(node_data) %>% dplyr::arrange(class, desc(Degree)),
#                colNames = TRUE, rowNames = FALSE)
# writeDataTable(wb, sheet = 2, x = as.data.frame(edge_data),
#                colNames = TRUE, rowNames = FALSE)
# writeDataTable(wb, sheet = 3, x = microbiome_info,
#                colNames = TRUE, rowNames = FALSE)
# writeDataTable(wb, sheet = 4, x = molecular_info,
#                colNames = TRUE, rowNames = FALSE)
# saveWorkbook(wb = wb, file = "inter_omics_cor_network.xlsx", overwrite = TRUE)


######pathway enrichment for proteins in each microbiome
stool_proteome =
  edge_data %>%
  dplyr::left_join(node_data[, c("node", "class")],
                   by = c("from" = "node")) %>%
  dplyr::rename(from_class = class) %>%
  dplyr::left_join(node_data[, c("node", "class")],
                   by = c("to" = "node")) %>%
  dplyr::rename(to_class = class) %>%
  dplyr::filter(from_class == "Stool microbiome" &
                  to_class == "Proteome") %>%
  # dplyr::filter(direction == "positive") %>%
  dplyr::select(to) %>%
  dplyr::distinct(.keep_all = TRUE) %>% 
  dplyr::left_join(proteome_variable_info, by = c("to" = "variable_id")) %>% 
  dplyr::filter(!is.na(UNIPROT) & !is.na(ENTREZID))


skin_proteome =
  edge_data %>%
  dplyr::left_join(node_data[, c("node", "class")],
                   by = c("from" = "node")) %>%
  dplyr::rename(from_class = class) %>%
  dplyr::left_join(node_data[, c("node", "class")],
                   by = c("to" = "node")) %>%
  dplyr::rename(to_class = class) %>%
  dplyr::filter(from_class == "Skin microbiome" &
                  to_class == "Proteome") %>%
  # dplyr::filter(direction == "positive") %>%
  dplyr::select(to) %>%
  dplyr::distinct(.keep_all = TRUE) %>% 
  dplyr::left_join(proteome_variable_info, by = c("to" = "variable_id")) %>% 
  dplyr::filter(!is.na(UNIPROT) & !is.na(ENTREZID))

oral_proteome =
  edge_data %>%
  dplyr::left_join(node_data[, c("node", "class")],
                   by = c("from" = "node")) %>%
  dplyr::rename(from_class = class) %>%
  dplyr::left_join(node_data[, c("node", "class")],
                   by = c("to" = "node")) %>%
  dplyr::rename(to_class = class) %>%
  dplyr::filter(from_class == "Oral microbiome" &
                  to_class == "Proteome") %>%
  # dplyr::filter(direction == "positive") %>%
  dplyr::select(to) %>%
  dplyr::distinct(.keep_all = TRUE) %>% 
  dplyr::left_join(proteome_variable_info, by = c("to" = "variable_id")) %>% 
  dplyr::filter(!is.na(UNIPROT) & !is.na(ENTREZID))

nasal_proteome =
  edge_data %>%
  dplyr::left_join(node_data[, c("node", "class")],
                   by = c("from" = "node")) %>%
  dplyr::rename(from_class = class) %>%
  dplyr::left_join(node_data[, c("node", "class")],
                   by = c("to" = "node")) %>%
  dplyr::rename(to_class = class) %>%
  dplyr::filter(from_class == "Nasal microbiome" &
                  to_class == "Proteome") %>%
  # dplyr::filter(direction == "positive") %>%
  dplyr::select(to) %>%
  dplyr::distinct(.keep_all = TRUE) %>% 
  dplyr::left_join(proteome_variable_info, by = c("to" = "variable_id")) %>% 
  dplyr::filter(!is.na(UNIPROT) & !is.na(ENTREZID))

library(clusterProfiler)

# stool_proteome_go =
#   clusterProfiler::enrichGO(
#     gene = stool_proteome$ENTREZID,
#     keyType = "ENTREZID",
#     ont = "ALL",
#     pvalueCutoff = 0.05,
#     pAdjustMethod = "BH",
#     qvalueCutoff = 0.05,
#     OrgDb = org.Hs.eg.db
#   )
# 
# save(stool_proteome_go, file = "pathway_enrichment/stool_proteome_go")
# 
# 
# stool_proteome_kegg =
#   clusterProfiler::enrichKEGG(
#     gene = stool_proteome$UNIPROT,
#     organism = "hsa",
#     keyType = "uniprot",
#     pvalueCutoff = 0.05,
#     pAdjustMethod = "BH",
#     qvalueCutoff = 0.05
#   )
# 
# save(stool_proteome_kegg, file = "pathway_enrichment/stool_proteome_kegg")
# 
# stool_proteome_kegg@result$p.adjust
# 
# 
# 
# skin_proteome_go =
#   clusterProfiler::enrichGO(
#     gene = skin_proteome$ENTREZID,
#     keyType = "ENTREZID",
#     ont = "ALL",
#     pvalueCutoff = 0.05,
#     pAdjustMethod = "BH",
#     qvalueCutoff = 0.05,
#     OrgDb = org.Hs.eg.db
#   )
# 
# save(skin_proteome_go, file = "pathway_enrichment/skin_proteome_go")
# 
# skin_proteome_kegg =
#   clusterProfiler::enrichKEGG(
#     gene = skin_proteome$UNIPROT,
#     organism = "hsa",
#     keyType = "uniprot",
#     pvalueCutoff = 0.05,
#     pAdjustMethod = "BH",
#     qvalueCutoff = 0.05
#   )
# 
# save(skin_proteome_kegg, file = "pathway_enrichment/skin_proteome_kegg")
# 
# skin_proteome_kegg@result$p.adjust
# 
# 
# oral_proteome_go =
#   clusterProfiler::enrichGO(
#     gene = oral_proteome$ENTREZID,
#     keyType = "ENTREZID",
#     ont = "ALL",
#     pvalueCutoff = 0.05,
#     pAdjustMethod = "BH",
#     qvalueCutoff = 0.05,
#     OrgDb = org.Hs.eg.db
#   )
# 
# save(oral_proteome_go, file = "pathway_enrichment/oral_proteome_go")
# 
# oral_proteome_kegg =
#   clusterProfiler::enrichKEGG(
#     gene = oral_proteome$UNIPROT,
#     organism = "hsa",
#     keyType = "uniprot",
#     pvalueCutoff = 0.05,
#     pAdjustMethod = "BH",
#     qvalueCutoff = 0.05
#   )
# 
# save(oral_proteome_kegg, file = "pathway_enrichment/oral_proteome_kegg")
# 
# oral_proteome_kegg@result$p.adjust
# 
# 
# nasal_proteome_go =
#   clusterProfiler::enrichGO(
#     gene = nasal_proteome$ENTREZID,
#     keyType = "ENTREZID",
#     ont = "ALL",
#     pvalueCutoff = 0.05,
#     pAdjustMethod = "BH",
#     qvalueCutoff = 0.05,
#     OrgDb = org.Hs.eg.db
#   )
# 
# save(nasal_proteome_go, file = "pathway_enrichment/nasal_proteome_go")
# 
# nasal_proteome_kegg =
#   clusterProfiler::enrichKEGG(
#     gene = nasal_proteome$UNIPROT,
#     organism = "hsa",
#     keyType = "uniprot",
#     pvalueCutoff = 0.05,
#     pAdjustMethod = "BH",
#     qvalueCutoff = 0.05
#   )
# 
# save(nasal_proteome_kegg, file = "pathway_enrichment/nasal_proteome_kegg")
# 
# nasal_proteome_kegg@result$p.adjust

load("pathway_enrichment/proteome/nasal/nasal_proteome_go")
load("pathway_enrichment/proteome/nasal/nasal_proteome_kegg")

load("pathway_enrichment/proteome/oral/oral_proteome_go")
load("pathway_enrichment/proteome/oral/oral_proteome_kegg")

load("pathway_enrichment/proteome/skin/skin_proteome_go")
load("pathway_enrichment/proteome/skin/skin_proteome_kegg")

load("pathway_enrichment/proteome/stool/stool_proteome_go")
load("pathway_enrichment/proteome/stool/stool_proteome_kegg")

nasal_proteome_go@result$geneName =
  nasal_proteome_go@result$geneID %>%
  purrr::map(function(x) {
    x = stringr::str_split(x, "/")[[1]]
    paste(proteome_variable_info$variable_id[match(x, proteome_variable_info$ENTREZID)], collapse = "/")
  }) %>%
  unlist()

oral_proteome_go@result$geneName =
  oral_proteome_go@result$geneID %>%
  purrr::map(function(x) {
    x = stringr::str_split(x, "/")[[1]]
    paste(proteome_variable_info$variable_id[match(x, proteome_variable_info$ENTREZID)], collapse = "/")
  }) %>%
  unlist()

stool_proteome_go@result$geneName =
  stool_proteome_go@result$geneID %>%
  purrr::map(function(x) {
    x = stringr::str_split(x, "/")[[1]]
    paste(proteome_variable_info$variable_id[match(x, proteome_variable_info$ENTREZID)], collapse = "/")
  }) %>%
  unlist()

skin_proteome_go@result$geneName =
  skin_proteome_go@result$geneID %>%
  purrr::map(function(x) {
    x = stringr::str_split(x, "/")[[1]]
    paste(proteome_variable_info$variable_id[match(x, proteome_variable_info$ENTREZID)], collapse = "/")
  }) %>%
  unlist()

nasal_proteome_kegg@result$geneName =
  nasal_proteome_kegg@result$geneID %>%
  purrr::map(function(x) {
    x = stringr::str_split(x, "/")[[1]]
    paste(proteome_variable_info$variable_id[match(x, proteome_variable_info$UNIPROT)], collapse = "/")
  }) %>%
  unlist()

oral_proteome_kegg@result$geneName =
  oral_proteome_kegg@result$geneID %>%
  purrr::map(function(x) {
    x = stringr::str_split(x, "/")[[1]]
    paste(proteome_variable_info$variable_id[match(x, proteome_variable_info$UNIPROT)], collapse = "/")
  }) %>%
  unlist()

stool_proteome_kegg@result$geneName =
  stool_proteome_kegg@result$geneID %>%
  purrr::map(function(x) {
    x = stringr::str_split(x, "/")[[1]]
    paste(proteome_variable_info$variable_id[match(x, proteome_variable_info$UNIPROT)], collapse = "/")
  }) %>%
  unlist()

skin_proteome_kegg@result$geneName =
  skin_proteome_kegg@result$geneID %>%
  purrr::map(function(x) {
    x = stringr::str_split(x, "/")[[1]]
    paste(proteome_variable_info$variable_id[match(x, proteome_variable_info$UNIPROT)], collapse = "/")
  }) %>%
  unlist()

######pathway enrichment for metabolites in each microbiome
stool_metabolome =
  edge_data %>%
  dplyr::left_join(node_data[, c("node", "class")],
                   by = c("from" = "node")) %>%
  dplyr::rename(from_class = class) %>%
  dplyr::left_join(node_data[, c("node", "class")],
                   by = c("to" = "node")) %>%
  dplyr::rename(to_class = class) %>%
  dplyr::filter(from_class == "Stool microbiome" &
                  to_class == "Metabolite") %>%
  # dplyr::filter(direction == "positive") %>%
  dplyr::select(to) %>%
  dplyr::distinct(.keep_all = TRUE) %>%
  dplyr::left_join(metabolome_variable_info[, c("variable_id", "HMDB", "KEGG")],
                   by = c("to" = "variable_id")) %>%
  dplyr::filter(!is.na(HMDB) | !is.na(KEGG)) %>%
  dplyr::filter(HMDB != "" | KEGG != "")


skin_metabolome =
  edge_data %>%
  dplyr::left_join(node_data[, c("node", "class")],
                   by = c("from" = "node")) %>%
  dplyr::rename(from_class = class) %>%
  dplyr::left_join(node_data[, c("node", "class")],
                   by = c("to" = "node")) %>%
  dplyr::rename(to_class = class) %>%
  dplyr::filter(from_class == "Skin microbiome" &
                  to_class == "Metabolite") %>%
  # dplyr::filter(direction == "positive") %>%
  dplyr::select(to) %>%
  dplyr::distinct(.keep_all = TRUE) %>% 
  dplyr::left_join(metabolome_variable_info[,c("variable_id", "HMDB", "KEGG")], 
                   by = c("to" = "variable_id")) %>% 
  dplyr::filter(!is.na(HMDB) | !is.na(KEGG)) %>% 
  dplyr::filter(HMDB != "" | KEGG != "")

oral_metabolome =
  edge_data %>%
  dplyr::left_join(node_data[, c("node", "class")],
                   by = c("from" = "node")) %>%
  dplyr::rename(from_class = class) %>%
  dplyr::left_join(node_data[, c("node", "class")],
                   by = c("to" = "node")) %>%
  dplyr::rename(to_class = class) %>%
  dplyr::filter(from_class == "Oral microbiome" &
                  to_class == "Metabolite") %>%
  # dplyr::filter(direction == "positive") %>%
  dplyr::select(to) %>%
  dplyr::distinct(.keep_all = TRUE) %>% 
  dplyr::left_join(metabolome_variable_info[,c("variable_id", "HMDB", "KEGG")], 
                   by = c("to" = "variable_id")) %>% 
  dplyr::filter(!is.na(HMDB) | !is.na(KEGG)) %>% 
  dplyr::filter(HMDB != "" | KEGG != "")

nasal_metabolome =
  edge_data %>%
  dplyr::left_join(node_data[, c("node", "class")],
                   by = c("from" = "node")) %>%
  dplyr::rename(from_class = class) %>%
  dplyr::left_join(node_data[, c("node", "class")],
                   by = c("to" = "node")) %>%
  dplyr::rename(to_class = class) %>%
  dplyr::filter(from_class == "Nasal microbiome" &
                  to_class == "Metabolite") %>%
  # dplyr::filter(direction == "positive") %>%
  dplyr::select(to) %>%
  dplyr::distinct(.keep_all = TRUE) %>% 
  dplyr::left_join(metabolome_variable_info[,c("variable_id", "HMDB", "KEGG")], 
                   by = c("to" = "variable_id")) %>% 
  dplyr::filter(!is.na(HMDB) | !is.na(KEGG)) %>% 
  dplyr::filter(HMDB != "" | KEGG != "")

library(clusterProfiler)
library(metPath)

stool_metabolome

data("hmdb_pathway", package = "metPath")
hmdb_pathway

get_pathway_class(hmdb_pathway)

remain_idx =
  hmdb_pathway@pathway_class %>%
  unlist() %>%
  stringr::str_detect("Metabolic;primary_pathway") %>%
  which()

hmdb_pathway =
  filter_pathway(object = hmdb_pathway, remain_idx = remain_idx)

data("kegg_hsa_pathway", package = "metPath")

remain_idx =
  kegg_hsa_pathway@pathway_class %>%
  unlist() %>%
  stringr::str_detect("Disease") %>%
  `!`() %>%
  which()

pathway_database =
  filter_pathway(object = kegg_hsa_pathway, remain_idx = remain_idx)

# stool_metabolome_hmdb =
#   enrich_hmdb(
#     query_id = unique(unlist(stringr::str_split(stool_metabolome$HMDB[!is.na(stool_metabolome$HMDB)], "\\|"))),
#     query_type = "compound",
#     id_type = "HMDB",
#     pathway_database = hmdb_pathway,
#     only_primary_pathway = TRUE,
#     p_cutoff = 0.05,
#     p_adjust_method = "BH",
#     threads = 3
#   )
#
# save(stool_metabolome_hmdb, file = "pathway_enrichment/stool_metabolome_hmdb")
#
#
# stool_metabolome_kegg =
#   enrich_kegg(query_id = unique(unlist(stringr::str_split(stool_metabolome$KEGG[!is.na(stool_metabolome$KEGG)], "\\|"))),
#               query_type = "compound",
#               id_type = "KEGG",
#               pathway_database = pathway_database,
#               p_cutoff = 0.05,
#               p_adjust_method = "BH",
#               threads = 3)
#
# save(stool_metabolome_kegg, file = "pathway_enrichment/stool_metabolome_kegg")
#
#
#
#
# skin_metabolome_hmdb =
#   enrich_hmdb(
#     query_id = unique(unlist(stringr::str_split(skin_metabolome$HMDB[!is.na(skin_metabolome$HMDB)], "\\|"))),
#     query_type = "compound",
#     id_type = "HMDB",
#     pathway_database = hmdb_pathway,
#     only_primary_pathway = TRUE,
#     p_cutoff = 0.05,
#     p_adjust_method = "BH",
#     threads = 3
#   )
#
# save(skin_metabolome_hmdb, file = "pathway_enrichment/skin_metabolome_hmdb")
#
#
# skin_metabolome_kegg =
#   enrich_kegg(query_id = unique(unlist(stringr::str_split(skin_metabolome$KEGG[!is.na(skin_metabolome$KEGG)], "\\|"))),
#               query_type = "compound",
#               id_type = "KEGG",
#               pathway_database = pathway_database,
#               p_cutoff = 0.05,
#               p_adjust_method = "BH",
#               threads = 3)
#
# save(skin_metabolome_kegg, file = "pathway_enrichment/skin_metabolome_kegg")
#
#
# oral_metabolome_hmdb =
#   enrich_hmdb(
#     query_id = unique(unlist(stringr::str_split(oral_metabolome$HMDB[!is.na(oral_metabolome$HMDB)], "\\|"))),
#     query_type = "compound",
#     id_type = "HMDB",
#     pathway_database = hmdb_pathway,
#     only_primary_pathway = TRUE,
#     p_cutoff = 0.05,
#     p_adjust_method = "BH",
#     threads = 3
#   )
#
# save(oral_metabolome_hmdb, file = "pathway_enrichment/oral_metabolome_hmdb")
#
#
# oral_metabolome_kegg =
#   enrich_kegg(query_id = unique(unlist(stringr::str_split(oral_metabolome$KEGG[!is.na(oral_metabolome$KEGG)], "\\|"))),
#               query_type = "compound",
#               id_type = "KEGG",
#               pathway_database = pathway_database,
#               p_cutoff = 0.05,
#               p_adjust_method = "BH",
#               threads = 3)
#
# save(oral_metabolome_kegg, file = "pathway_enrichment/oral_metabolome_kegg")
#
#
# nasal_metabolome_hmdb =
#   enrich_hmdb(
#     query_id = unique(unlist(stringr::str_split(nasal_metabolome$HMDB[!is.na(nasal_metabolome$HMDB)], "\\|"))),
#     query_type = "compound",
#     id_type = "HMDB",
#     pathway_database = hmdb_pathway,
#     only_primary_pathway = TRUE,
#     p_cutoff = 0.05,
#     p_adjust_method = "BH",
#     threads = 3
#   )
#
# save(nasal_metabolome_hmdb, file = "pathway_enrichment/nasal_metabolome_hmdb")
#
#
# nasal_metabolome_kegg =
#   enrich_kegg(query_id = unique(unlist(stringr::str_split(nasal_metabolome$KEGG[!is.na(nasal_metabolome$KEGG)], "\\|"))),
#               query_type = "compound",
#               id_type = "KEGG",
#               pathway_database = pathway_database,
#               p_cutoff = 0.05,
#               p_adjust_method = "BH",
#               threads = 3)
#
# save(nasal_metabolome_kegg, file = "pathway_enrichment/nasal_metabolome_kegg")

load("pathway_enrichment/metabolome/nasal/nasal_metabolome_hmdb")
load("pathway_enrichment/metabolome/nasal/nasal_metabolome_kegg")

load("pathway_enrichment/metabolome/oral/oral_metabolome_hmdb")
load("pathway_enrichment/metabolome/oral/oral_metabolome_kegg")

load("pathway_enrichment/metabolome/skin/skin_metabolome_hmdb")
load("pathway_enrichment/metabolome/skin/skin_metabolome_kegg")

load("pathway_enrichment/metabolome/stool/stool_metabolome_hmdb")
load("pathway_enrichment/metabolome/stool/stool_metabolome_kegg")

nasal_metabolome_hmdb@result$mapped_metabolite =
  nasal_metabolome_hmdb@result$mapped_id %>%
  purrr::map(function(x) {
    if (all(x == "")) {
      return(x)
    }
    x = stringr::str_split(x, ";")[[1]]
    paste(metabolome_variable_info$Metabolite[match(x, metabolome_variable_info$HMDB)], collapse = ";")
  }) %>%
  unlist()

oral_metabolome_hmdb@result$mapped_metabolite =
  oral_metabolome_hmdb@result$mapped_id %>%
  purrr::map(function(x) {
    if (all(x == "")) {
      return(x)
    }
    x = stringr::str_split(x, ";")[[1]]
    paste(metabolome_variable_info$Metabolite[match(x, metabolome_variable_info$HMDB)], collapse = ";")
  }) %>%
  unlist()

stool_metabolome_hmdb@result$mapped_metabolite =
  stool_metabolome_hmdb@result$mapped_id %>%
  purrr::map(function(x) {
    if (all(x == "")) {
      return(x)
    }
    x = stringr::str_split(x, ";")[[1]]
    paste(metabolome_variable_info$Metabolite[match(x, metabolome_variable_info$HMDB)], collapse = ";")
  }) %>%
  unlist()

skin_metabolome_hmdb@result$mapped_metabolite =
  skin_metabolome_hmdb@result$mapped_id %>%
  purrr::map(function(x) {
    if (all(x == "")) {
      return(x)
    }
    x = stringr::str_split(x, ";")[[1]]
    paste(metabolome_variable_info$Metabolite[match(x, metabolome_variable_info$HMDB)], collapse = ";")
  }) %>%
  unlist()

nasal_metabolome_kegg@result$mapped_metabolite =
  nasal_metabolome_kegg@result$mapped_id %>%
  purrr::map(function(x) {
    if (all(x == "")) {
      return(x)
    }
    x = stringr::str_split(x, ";")[[1]]
    paste(metabolome_variable_info$Metabolite[match(x, metabolome_variable_info$KEGG)], collapse = ";")
  }) %>%
  unlist()

oral_metabolome_kegg@result$mapped_metabolite =
  oral_metabolome_kegg@result$mapped_id %>%
  purrr::map(function(x) {
    if (all(x == "")) {
      return(x)
    }
    x = stringr::str_split(x, ";")[[1]]
    paste(metabolome_variable_info$Metabolite[match(x, metabolome_variable_info$KEGG)], collapse = ";")
  }) %>%
  unlist()

stool_metabolome_kegg@result$mapped_metabolite =
  stool_metabolome_kegg@result$mapped_id %>%
  purrr::map(function(x) {
    if (all(x == "")) {
      return(x)
    }
    x = stringr::str_split(x, ";")[[1]]
    paste(metabolome_variable_info$Metabolite[match(x, metabolome_variable_info$KEGG)], collapse = ";")
  }) %>%
  unlist()

skin_metabolome_kegg@result$mapped_metabolite =
  skin_metabolome_kegg@result$mapped_id %>%
  purrr::map(function(x) {
    if (all(x == "")) {
      return(x)
    }
    x = stringr::str_split(x, ";")[[1]]
    paste(metabolome_variable_info$Metabolite[match(x, metabolome_variable_info$KEGG)], collapse = ";")
  }) %>%
  unlist()

# ####output result (proteome)
# library(openxlsx)
# wb <- createWorkbook()
# modifyBaseFont(wb, fontSize = 12, fontName = "Times New Roma")
# addWorksheet(wb, sheetName = "GO Stool microbiome", gridLines = TRUE)
# addWorksheet(wb, sheetName = "GO Skin microbiome", gridLines = TRUE)
# addWorksheet(wb, sheetName = "GO Oral microbiome", gridLines = TRUE)
# addWorksheet(wb, sheetName = "GO Nasal microbiome", gridLines = TRUE)
# addWorksheet(wb, sheetName = "KEGG Stool microbiome", gridLines = TRUE)
# addWorksheet(wb, sheetName = "KEGG Skin microbiome", gridLines = TRUE)
# addWorksheet(wb, sheetName = "KEGG Oral microbiome", gridLines = TRUE)
# addWorksheet(wb, sheetName = "KEGG Nasal microbiome", gridLines = TRUE)
# freezePane(wb, sheet = 1, firstRow = TRUE, firstCol = TRUE)
# freezePane(wb, sheet = 2, firstRow = TRUE, firstCol = TRUE)
# freezePane(wb, sheet = 3, firstRow = TRUE, firstCol = TRUE)
# freezePane(wb, sheet = 4, firstRow = TRUE, firstCol = TRUE)
# freezePane(wb, sheet = 5, firstRow = TRUE, firstCol = TRUE)
# freezePane(wb, sheet = 6, firstRow = TRUE, firstCol = TRUE)
# freezePane(wb, sheet = 7, firstRow = TRUE, firstCol = TRUE)
# freezePane(wb, sheet = 8, firstRow = TRUE, firstCol = TRUE)
# 
# writeDataTable(
#   wb,
#   sheet = 1,
#   x = stool_proteome_go@result %>% dplyr::arrange(p.adjust) %>% dplyr::filter(Count > 0),
#   colNames = TRUE,
#   rowNames = FALSE
# )
# 
# writeDataTable(
#   wb,
#   sheet = 2,
#   x = skin_proteome_go@result %>% dplyr::arrange(p.adjust) %>% dplyr::filter(Count > 0),
#   colNames = TRUE,
#   rowNames = FALSE
# )
# 
# writeDataTable(
#   wb,
#   sheet = 3,
#   x = oral_proteome_go@result %>% dplyr::arrange(p.adjust) %>% dplyr::filter(Count > 0),
#   colNames = TRUE,
#   rowNames = FALSE
# )
# 
# writeDataTable(
#   wb,
#   sheet = 4,
#   x = nasal_proteome_go@result %>% dplyr::arrange(p.adjust) %>% dplyr::filter(Count > 0),
#   colNames = TRUE,
#   rowNames = FALSE
# )
# 
# writeDataTable(
#   wb,
#   sheet = 5,
#   x = stool_proteome_kegg@result %>% dplyr::arrange(p.adjust) %>% dplyr::filter(Count > 0),
#   colNames = TRUE,
#   rowNames = FALSE
# )
# 
# writeDataTable(
#   wb,
#   sheet = 6,
#   x = skin_proteome_kegg@result %>% dplyr::arrange(p.adjust) %>% dplyr::filter(Count > 0),
#   colNames = TRUE,
#   rowNames = FALSE
# )
# 
# writeDataTable(
#   wb,
#   sheet = 7,
#   x = oral_proteome_kegg@result %>% dplyr::arrange(p.adjust) %>% dplyr::filter(Count > 0),
#   colNames = TRUE,
#   rowNames = FALSE
# )
# 
# writeDataTable(
#   wb,
#   sheet = 8,
#   x = nasal_proteome_kegg@result %>% dplyr::arrange(p.adjust) %>% dplyr::filter(Count > 0),
#   colNames = TRUE,
#   rowNames = FALSE
# )
# 
# saveWorkbook(wb, "pathway_enrichment/proteome_pathway_enrichment_result.xlsx", overwrite = TRUE)


# ###output result (metabolome)
# library(openxlsx)
# wb <- createWorkbook()
# modifyBaseFont(wb, fontSize = 12, fontName = "Times New Roma")
# addWorksheet(wb, sheetName = "Stool microbiome HMDB", gridLines = TRUE)
# addWorksheet(wb, sheetName = "Skin microbiome HMDB", gridLines = TRUE)
# addWorksheet(wb, sheetName = "Oral microbiome HMDB", gridLines = TRUE)
# addWorksheet(wb, sheetName = "Nasal microbiome HMDB", gridLines = TRUE)
# addWorksheet(wb, sheetName = "Stool microbiome KEGG", gridLines = TRUE)
# addWorksheet(wb, sheetName = "Skin microbiome KEGG", gridLines = TRUE)
# addWorksheet(wb, sheetName = "Oral microbiome KEGG", gridLines = TRUE)
# addWorksheet(wb, sheetName = "Nasal microbiome KEGG", gridLines = TRUE)
# freezePane(wb, sheet = 1, firstRow = TRUE, firstCol = TRUE)
# freezePane(wb, sheet = 2, firstRow = TRUE, firstCol = TRUE)
# freezePane(wb, sheet = 3, firstRow = TRUE, firstCol = TRUE)
# freezePane(wb, sheet = 4, firstRow = TRUE, firstCol = TRUE)
# freezePane(wb, sheet = 5, firstRow = TRUE, firstCol = TRUE)
# freezePane(wb, sheet = 6, firstRow = TRUE, firstCol = TRUE)
# freezePane(wb, sheet = 7, firstRow = TRUE, firstCol = TRUE)
# freezePane(wb, sheet = 8, firstRow = TRUE, firstCol = TRUE)
# 
# writeDataTable(
#   wb,
#   sheet = 1,
#   x = stool_metabolome_hmdb@result %>% dplyr::arrange(p_value_adjust) %>% dplyr::filter(mapped_number > 1),
#   colNames = TRUE,
#   rowNames = FALSE
# )
# 
# writeDataTable(
#   wb,
#   sheet = 2,
#   x = skin_metabolome_hmdb@result %>% dplyr::arrange(p_value_adjust) %>% dplyr::filter(mapped_number > 1),
#   colNames = TRUE,
#   rowNames = FALSE
# )
# 
# writeDataTable(
#   wb,
#   sheet = 3,
#   x = oral_metabolome_hmdb@result %>% dplyr::arrange(p_value_adjust) %>% dplyr::filter(mapped_number > 1),
#   colNames = TRUE,
#   rowNames = FALSE
# )
# 
# writeDataTable(
#   wb,
#   sheet = 4,
#   x = nasal_metabolome_hmdb@result %>% dplyr::arrange(p_value_adjust) %>% dplyr::filter(mapped_number > 1),
#   colNames = TRUE,
#   rowNames = FALSE
# )
# 
# writeDataTable(
#   wb,
#   sheet = 5,
#   x = stool_metabolome_kegg@result %>% dplyr::arrange(p_value_adjust) %>% dplyr::filter(mapped_number > 1),
#   colNames = TRUE,
#   rowNames = FALSE
# )
# 
# writeDataTable(
#   wb,
#   sheet = 6,
#   x = skin_metabolome_kegg@result %>% dplyr::arrange(p_value_adjust) %>% dplyr::filter(mapped_number > 1),
#   colNames = TRUE,
#   rowNames = FALSE
# )
# 
# writeDataTable(
#   wb,
#   sheet = 7,
#   x = oral_metabolome_kegg@result %>% dplyr::arrange(p_value_adjust) %>% dplyr::filter(mapped_number > 1),
#   colNames = TRUE,
#   rowNames = FALSE
# )
# 
# writeDataTable(
#   wb,
#   sheet = 8,
#   x = nasal_metabolome_kegg@result %>% dplyr::arrange(p_value_adjust) %>% dplyr::filter(mapped_number > 1),
#   colNames = TRUE,
#   rowNames = FALSE
# )
# 
# saveWorkbook(wb, "pathway_enrichment/metabolome_pathway_enrichment_result.xlsx", overwrite = TRUE)



######for proteome, we need to get the similarity between GO terms
########for stool microbiome
stool_go_result = readxl::read_xlsx("pathway_enrichment/proteome/proteome_pathway_enrichment_result_manual.xlsx",
                                    sheet = 1)
skin_go_result = readxl::read_xlsx("pathway_enrichment/proteome/proteome_pathway_enrichment_result_manual.xlsx",
                                   sheet = 2)
oral_go_result = readxl::read_xlsx("pathway_enrichment/proteome/proteome_pathway_enrichment_result_manual.xlsx",
                                   sheet = 3)
nasal_go_result = readxl::read_xlsx("pathway_enrichment/proteome/proteome_pathway_enrichment_result_manual.xlsx",
                                    sheet = 4)
# stool_sim_matrix = get_go_sim(result = stool_go_result, sim_cutoff = 0.7)
# skin_sim_matrix = get_go_sim(result = skin_go_result, sim_cutoff = 0.7)
# oral_sim_matrix = get_go_sim(result = oral_go_result, sim_cutoff = 0.7)
# nasal_sim_matrix = get_go_sim(result = nasal_go_result, sim_cutoff = 0.7)
# 
# save(stool_sim_matrix, file = "pathway_enrichment/proteome/stool/stool_sim_matrix")
# save(skin_sim_matrix, file = "pathway_enrichment/proteome/skin/skin_sim_matrix")
# save(oral_sim_matrix, file = "pathway_enrichment/proteome/oral/oral_sim_matrix")
# save(nasal_sim_matrix, file = "pathway_enrichment/proteome/nasal/nasal_sim_matrix")

load("pathway_enrichment/proteome/stool/stool_sim_matrix")
load("pathway_enrichment/proteome/skin/skin_sim_matrix")
load("pathway_enrichment/proteome/oral/oral_sim_matrix")
load("pathway_enrichment/proteome/nasal/nasal_sim_matrix")

######stool microbiome to find the cluster of GO terms
# get_go_cluster(result = stool_go_result,
#                sim_matrix = stool_sim_matrix,
#                output_path = "pathway_enrichment/proteome/stool/")
# 
# get_go_cluster(result = skin_go_result,
#                sim_matrix = skin_sim_matrix,
#                output_path = "pathway_enrichment/proteome/skin/")
# 
# get_go_cluster(result = oral_go_result,
#                sim_matrix = oral_sim_matrix,
#                output_path = "pathway_enrichment/proteome/oral/")
# 
# get_go_cluster(result = nasal_go_result,
#                sim_matrix = nasal_sim_matrix,
#                output_path = "pathway_enrichment/proteome/nasal/")

load("pathway_enrichment/proteome/stool/result")
load("pathway_enrichment/proteome/stool/result_cluster")
stool_result = result
stool_result_cluster = result_cluster

load("pathway_enrichment/proteome/skin/result")
load("pathway_enrichment/proteome/skin/result_cluster")
skin_result = result
skin_result_cluster = result_cluster

load("pathway_enrichment/proteome/oral/result")
load("pathway_enrichment/proteome/oral/result_cluster")
oral_result = result
oral_result_cluster = result_cluster

load("pathway_enrichment/proteome/nasal/result")
load("pathway_enrichment/proteome/nasal/result_cluster")
nasal_result = result
nasal_result_cluster = result_cluster

###output the cluster annotation for each cluster
dir.create("pathway_enrichment/proteome/stool/GO_graph")
dir.create("pathway_enrichment/proteome/skin/GO_graph")
dir.create("pathway_enrichment/proteome/oral/GO_graph")
dir.create("pathway_enrichment/proteome/nasal/GO_graph")

#####output GO result
# output_go_structure(result_cluster = stool_result_cluster,
#                     output_path = "pathway_enrichment/proteome/stool/GO_graph/")
# 
# output_go_structure(result_cluster = skin_result_cluster,
#                     output_path = "pathway_enrichment/proteome/skin/GO_graph/")
# 
# output_go_structure(result_cluster = oral_result_cluster,
#                     output_path = "pathway_enrichment/proteome/oral/GO_graph/")
# 
# output_go_structure(result_cluster = nasal_result_cluster,
#                     output_path = "pathway_enrichment/proteome/nasal/GO_graph/")

dim(stool_result_cluster)

stool_result_cluster =
  readxl::read_xlsx("pathway_enrichment/proteome/stool/result_cluster.xlsx")
skin_result_cluster =
  readxl::read_xlsx("pathway_enrichment/proteome/skin/result_cluster.xlsx")
oral_result_cluster =
  readxl::read_xlsx("pathway_enrichment/proteome/oral/result_cluster.xlsx")
nasal_result_cluster =
  readxl::read_xlsx("pathway_enrichment/proteome/nasal/result_cluster.xlsx")

###only show the top 10 enriched pathways?
head(stool_result_cluster$cluster_annotation, 10)
head(skin_result_cluster$cluster_annotation, 10)
head(oral_result_cluster$cluster_annotation, 10)
head(nasal_result_cluster$cluster_annotation, 10)

#####output them to a xlsx
# library(openxlsx)
# wb <- createWorkbook()
# modifyBaseFont(wb, fontSize = 12, fontName = "Times New Roma")
# addWorksheet(wb, sheetName = "Stool microbiome", gridLines = TRUE)
# addWorksheet(wb, sheetName = "Skin microbiome", gridLines = TRUE)
# addWorksheet(wb, sheetName = "Oral microbiome", gridLines = TRUE)
# addWorksheet(wb, sheetName = "Nasal microbiome", gridLines = TRUE)
# freezePane(wb, sheet = 1, firstRow = TRUE, firstCol = TRUE)
# freezePane(wb, sheet = 2, firstRow = TRUE, firstCol = TRUE)
# freezePane(wb, sheet = 3, firstRow = TRUE, firstCol = TRUE)
# freezePane(wb, sheet = 4, firstRow = TRUE, firstCol = TRUE)
# 
# writeDataTable(
#   wb,
#   sheet = 1,
#   x = stool_result_cluster,
#   colNames = TRUE,
#   rowNames = FALSE
# )
# 
# writeDataTable(
#   wb,
#   sheet = 2,
#   x = skin_result_cluster,
#   colNames = TRUE,
#   rowNames = FALSE
# )
# 
# writeDataTable(
#   wb,
#   sheet = 3,
#   x = oral_result_cluster,
#   colNames = TRUE,
#   rowNames = FALSE
# )
# 
# writeDataTable(
#   wb,
#   sheet = 4,
#   x = nasal_result_cluster,
#   colNames = TRUE,
#   rowNames = FALSE
# )
# 
# saveWorkbook(wb, "pathway_enrichment/proteome/proteome_pathway_cluster.xlsx", overwrite = TRUE)

#######network to show the enrichment results protein
###only the BP and MF
temp_stool_result = 
stool_result_cluster %>% 
  dplyr::filter(ONTOLOGY != "CC") %>% 
  head(10) %>%
  dplyr::select(node) %>% 
  tidyr::separate_rows(node, sep = ";") %>% 
  dplyr::left_join(stool_result[,c("node", "p.adjust")], by = "node") %>% 
  dplyr::mutate(Stool = -log(p.adjust, 10)) %>% 
  dplyr::select(-p.adjust)

temp_skin_result = 
  skin_result_cluster %>% 
  dplyr::filter(ONTOLOGY != "CC") %>% 
  head(10) %>%
  dplyr::select(node) %>% 
  tidyr::separate_rows(node, sep = ";") %>% 
  dplyr::left_join(skin_result[,c("node", "p.adjust")], by = "node") %>% 
  dplyr::mutate(Skin = -log(p.adjust, 10)) %>% 
  dplyr::select(-p.adjust)

temp_oral_result = 
  oral_result_cluster %>% 
  dplyr::filter(ONTOLOGY != "CC") %>% 
  head(10) %>%
  dplyr::select(node) %>% 
  tidyr::separate_rows(node, sep = ";") %>% 
  dplyr::left_join(oral_result[,c("node", "p.adjust")], by = "node") %>% 
  dplyr::mutate(Oral = -log(p.adjust, 10)) %>% 
  dplyr::select(-p.adjust)

temp_nasal_result = 
  nasal_result_cluster %>% 
  dplyr::filter(ONTOLOGY != "CC") %>% 
  head(10) %>%
  dplyr::select(node) %>% 
  tidyr::separate_rows(node, sep = ";") %>% 
  dplyr::left_join(nasal_result[,c("node", "p.adjust")], by = "node") %>% 
  dplyr::mutate(Nasal = -log(p.adjust, 10)) %>% 
  dplyr::select(-p.adjust)

node_data =
  temp_stool_result %>%
  dplyr::full_join(temp_skin_result, by = "node") %>%
  dplyr::full_join(temp_oral_result, by = "node") %>%
  dplyr::full_join(temp_nasal_result, by = "node")

node_data$Stool[is.na(node_data$Stool)] = 0
node_data$Skin[is.na(node_data$Skin)] = 0
node_data$Oral[is.na(node_data$Oral)] = 0
node_data$Nasal[is.na(node_data$Nasal)] = 0

node_data = 
node_data %>%
  dplyr::left_join(
    rbind(stool_result, skin_result, nasal_result, oral_result) %>%
      dplyr::select(node, ONTOLOGY, Description) %>%
      dplyr::distinct(node, .keep_all = TRUE),
    by = "node"
  )

node_data_bp = 
  node_data %>% 
  dplyr::filter(ONTOLOGY == "BP")

node_data_mf = 
  node_data %>% 
  dplyr::filter(ONTOLOGY == "MF")

# calculate go similarity
# bp_sim =
#   simplifyEnrichment::GO_similarity(go_id = node_data_bp$node,
#                                     ont = "BP",
#                                     measure = "Wang") %>%
#   as.data.frame() %>%
#   tibble::rownames_to_column(var = "name1") %>%
#   tidyr::pivot_longer(cols = -name1,
#                       names_to = "name2",
#                       values_to = "sim") %>%
#   dplyr::filter(name1 != name2)
# 
# 
# mf_sim =
#   simplifyEnrichment::GO_similarity(go_id = node_data_mf$node,
#                                     ont = "MF",
#                                     measure = "Wang") %>%
#   as.data.frame() %>%
#   tibble::rownames_to_column(var = "name1") %>%
#   tidyr::pivot_longer(cols = -name1,
#                       names_to = "name2",
#                       values_to = "sim") %>%
#   dplyr::filter(name1 != name2)
# save(bp_sim, file = "pathway_enrichment/proteome/bp_sim")
# save(mf_sim, file = "pathway_enrichment/proteome/mf_sim")

load("pathway_enrichment/proteome/bp_sim")
load("pathway_enrichment/proteome/mf_sim")

bp_sim =
  bp_sim %>% 
  dplyr::filter(sim > 0.5)

mf_sim =
  mf_sim %>% 
  dplyr::filter(sim > 0.5)

###calculate the jaccard index between bp and mf
library(org.Hs.eg.db)
library(GO.db)

###get the gene id for each term
bp_term_gene = 
  purrr::map(.x = node_data_bp$node, function(node1) {
    unname(get(node1, org.Hs.egGO2ALLEGS))
  })

names(bp_term_gene) = node_data_bp$node

mf_term_gene = 
  purrr::map(.x = node_data_mf$node, function(node1) {
    unname(get(node1, org.Hs.egGO2ALLEGS))
  })

names(mf_term_gene) = node_data_mf$node

gene_number_bp =
  bp_term_gene %>%
  lapply(length) %>%
  unlist() %>%
  data.frame(gene_number = .) %>%
  tibble::rownames_to_column(var = "node")

gene_number_mf =
  mf_term_gene %>%
  lapply(length) %>%
  unlist() %>%
  data.frame(gene_number = .) %>%
  tibble::rownames_to_column(var = "node")

node_data =
  node_data %>%
  dplyr::left_join(rbind(gene_number_bp,
                         gene_number_mf), by = "node")


bp_mf_sim =
  purrr::map(.x = node_data_bp$node, function(node1) {
    purrr::map(.x = node_data_mf$node, function(node2) {
      gene_set1 = bp_term_gene[[which(names(bp_term_gene) == node1)]]
      gene_set2 = mf_term_gene[[which(names(mf_term_gene) == node2)]]
      
      data.frame(
        name1 = node1,
        name2 = node2,
        sim = length(intersect(gene_set1, gene_set2)) / length(union(gene_set1, gene_set2))
      )
    }) %>%
      do.call(rbind, .) %>%
      as.data.frame()
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()

sum(bp_mf_sim$sim > 0.5)

bp_mf_sim = 
  bp_mf_sim %>% 
  dplyr::filter(sim > 0.5)

edge_data = 
  rbind(data.frame(bp_sim, edge_type = "BP"),
        data.frame(mf_sim, edge_type = "MF"),
        data.frame(bp_mf_sim, edge_type = "BP_MF"))

dim(node_data)
dim(edge_data)  

library(ggraph)
library(tidygraph)

graph <- 
  tidygraph::tbl_graph(nodes = node_data, 
                       edges = edge_data,
                       directed = FALSE) %>% 
  dplyr::mutate(Degree = centrality_degree(mode = 'all'))

library(igraph)
library(graphlayouts)
xy = 
ggraph::create_layout(graph, layout = "kk")
V(graph)$x <- xy$x
V(graph)$y <- xy$y
library(scatterpie)
plot <-
ggraph(
  graph,
  layout = 'manual',
  x = V(graph)$x,
  y = V(graph)$y,
  circular = FALSE
) +
  geom_edge_link(aes(width = sim,
                   color = edge_type),
               alpha = 1,
               show.legend = TRUE) +
  geom_scatterpie(
    data = as_data_frame(graph, "vertices"),
    cols = c("Stool", "Skin", "Oral", "Nasal")
  ) +
  coord_equal() +
  geom_node_text(aes(label = Description),
                 repel = TRUE,
                 size = 2) +
  scale_fill_manual(values = body_site_color) +
  # shadowtext::geom_shadowtext(aes(x = x, y = y,
  #                                 label = ifelse(node %in% label_node, true_name, NA),
  #                                 color = class), 
  #                             size = 5,
  #                             check_overlap = TRUE,
  #                             bg.colour = "white",
  #                             show.legend = FALSE) +
  guides(edge_color = guide_edge_colorbar(title = "Edge type"),
         color = guide_legend(title = "Microbiome")) +
  scale_edge_width_continuous(range = c(0.1, 1)) +
  scale_size_continuous(range = c(2, 8)) +
  scale_edge_color_manual(values = c("BP" = "black",
                                     "MF" = "black",
                                     "BP_MF" = "grey")) +
  theme(plot.background = element_rect(fill = "transparent", color = NA), 
        panel.background = element_rect(fill = "transparent", color = NA),
        legend.position = "right",
        legend.background = element_rect(fill = "transparent", color = NA))

plot

# ggsave(plot, filename = "pathway_enrichment/proteome/enriched_pathways.pdf", width = 9, height = 7)



#######the CC
#######network to show the enrichment results protein
###only the CC
temp_stool_result = 
  stool_result_cluster %>% 
  dplyr::filter(ONTOLOGY == "CC") %>% 
  head(10) %>%
  dplyr::select(node) %>% 
  tidyr::separate_rows(node, sep = ";") %>% 
  dplyr::left_join(stool_result[,c("node", "p.adjust")], by = "node") %>% 
  dplyr::mutate(Stool = -log(p.adjust, 10)) %>% 
  dplyr::select(-p.adjust)

temp_skin_result = 
  skin_result_cluster %>% 
  dplyr::filter(ONTOLOGY == "CC") %>% 
  head(10) %>%
  dplyr::select(node) %>% 
  tidyr::separate_rows(node, sep = ";") %>% 
  dplyr::left_join(skin_result[,c("node", "p.adjust")], by = "node") %>% 
  dplyr::mutate(Skin = -log(p.adjust, 10)) %>% 
  dplyr::select(-p.adjust)

temp_oral_result = 
  oral_result_cluster %>% 
  dplyr::filter(ONTOLOGY == "CC") %>% 
  head(10) %>%
  dplyr::select(node) %>% 
  tidyr::separate_rows(node, sep = ";") %>% 
  dplyr::left_join(oral_result[,c("node", "p.adjust")], by = "node") %>% 
  dplyr::mutate(Oral = -log(p.adjust, 10)) %>% 
  dplyr::select(-p.adjust)

temp_nasal_result = 
  nasal_result_cluster %>% 
  dplyr::filter(ONTOLOGY == "CC") %>% 
  head(10) %>%
  dplyr::select(node) %>% 
  tidyr::separate_rows(node, sep = ";") %>% 
  dplyr::left_join(nasal_result[,c("node", "p.adjust")], by = "node") %>% 
  dplyr::mutate(Nasal = -log(p.adjust, 10)) %>% 
  dplyr::select(-p.adjust)

node_data =
  temp_stool_result %>%
  dplyr::full_join(temp_skin_result, by = "node") %>%
  dplyr::full_join(temp_oral_result, by = "node") %>%
  dplyr::full_join(temp_nasal_result, by = "node")

node_data$Stool[is.na(node_data$Stool)] = 0
node_data$Skin[is.na(node_data$Skin)] = 0
node_data$Oral[is.na(node_data$Oral)] = 0
node_data$Nasal[is.na(node_data$Nasal)] = 0

node_data = 
  node_data %>%
  dplyr::left_join(
    rbind(stool_result, skin_result, nasal_result, oral_result) %>%
      dplyr::select(node, ONTOLOGY, Description) %>%
      dplyr::distinct(node, .keep_all = TRUE),
    by = "node"
  )


#calculate similarity
cc_sim =
  simplifyEnrichment::GO_similarity(go_id = node_data$node,
                                    ont = "CC",
                                    measure = "Wang") %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "name1") %>%
  tidyr::pivot_longer(cols = -name1,
                      names_to = "name2",
                      values_to = "sim") %>%
  dplyr::filter(name1 != name2) %>%
  dplyr::filter(sim > 0.5)


###calculate the jaccard index between bp and mf
library(org.Hs.eg.db)
library(GO.db)

###get the gene id for each term
cc_term_gene = 
  purrr::map(.x = node_data$node, function(node1) {
    unname(get(node1, org.Hs.egGO2ALLEGS))
  })

names(cc_term_gene) = node_data$node

gene_number_cc =
  cc_term_gene %>%
  lapply(length) %>%
  unlist() %>%
  data.frame(gene_number = .) %>%
  tibble::rownames_to_column(var = "node")

node_data =
  node_data %>%
  dplyr::left_join(rbind(gene_number_cc), by = "node")


edge_data = 
  rbind(data.frame(cc_sim, edge_type = "CC"))

dim(node_data)
dim(edge_data)  

library(ggraph)
library(tidygraph)

graph <-
  tidygraph::tbl_graph(nodes = node_data,
                       edges = edge_data,
                       directed = FALSE) %>%
  dplyr::mutate(Degree = centrality_degree(mode = 'all'))

library(igraph)
library(graphlayouts)
xy = 
  ggraph::create_layout(graph, layout = "kk")
V(graph)$x <- xy$x
V(graph)$y <- xy$y

plot <-
  ggraph(
    graph,
    layout = 'manual',
    x = V(graph)$x,
    y = V(graph)$y,
    circular = FALSE
  ) +
  geom_edge_link(aes(width = sim,
                     color = edge_type),
                 alpha = 1,
                 show.legend = TRUE) +
  geom_scatterpie(
    data = as_data_frame(graph, "vertices"),
    cols = c("Stool", "Skin", "Oral", "Nasal")
  ) +
  coord_equal() +
  geom_node_text(aes(label = Description),
                 repel = TRUE,
                 size = 2) +
  scale_fill_manual(values = body_site_color) +
  # shadowtext::geom_shadowtext(aes(x = x, y = y,
  #                                 label = ifelse(node %in% label_node, true_name, NA),
  #                                 color = class), 
  #                             size = 5,
  #                             check_overlap = TRUE,
  #                             bg.colour = "white",
  #                             show.legend = FALSE) +
  guides(edge_color = guide_edge_colorbar(title = "Edge type"),
         color = guide_legend(title = "Microbiome")) +
  scale_edge_width_continuous(range = c(0.1, 1)) +
  scale_size_continuous(range = c(2, 8)) +
  scale_edge_color_manual(values = c("CC" = "black")) +
  theme(plot.background = element_rect(fill = "transparent", color = NA), 
        panel.background = element_rect(fill = "transparent", color = NA),
        legend.position = "right",
        legend.background = element_rect(fill = "transparent", color = NA))

plot

# ggsave(plot, filename = "pathway_enrichment/proteome/enriched_pathways_cc.pdf", width = 9, height = 7)



############metabolome
###because no significant pathways have been enriched, so here we try to use the class of metabolites
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

library(plyr)

edge_data =
  edge_data %>%
  dplyr::filter(from %in% node_data$node & to %in% node_data$node)

node_data =
  node_data %>%
  dplyr::filter(node %in% edge_data$from | node %in% edge_data$to)

node_data$node

dim(edge_data)
dim(node_data)

########metabolite
stool_metabolome =
  edge_data %>%
  dplyr::left_join(node_data[, c("node", "class")],
                   by = c("from" = "node")) %>%
  dplyr::rename(from_class = class) %>%
  dplyr::left_join(node_data[, c("node", "class")],
                   by = c("to" = "node")) %>%
  dplyr::rename(to_class = class) %>%
  dplyr::filter(from_class == "Stool microbiome" &
                  to_class == "Metabolite") %>%
  dplyr::select(to) %>%
  dplyr::distinct(.keep_all = TRUE) %>%
  dplyr::left_join(metabolome_variable_info,
                   by = c("to" = "variable_id"))

skin_metabolome =
  edge_data %>%
  dplyr::left_join(node_data[, c("node", "class")],
                   by = c("from" = "node")) %>%
  dplyr::rename(from_class = class) %>%
  dplyr::left_join(node_data[, c("node", "class")],
                   by = c("to" = "node")) %>%
  dplyr::rename(to_class = class) %>%
  dplyr::filter(from_class == "Skin microbiome" &
                  to_class == "Metabolite") %>%
  dplyr::select(to) %>%
  dplyr::distinct(.keep_all = TRUE) %>%
  dplyr::left_join(metabolome_variable_info,
                   by = c("to" = "variable_id"))

oral_metabolome =
  edge_data %>%
  dplyr::left_join(node_data[, c("node", "class")],
                   by = c("from" = "node")) %>%
  dplyr::rename(from_class = class) %>%
  dplyr::left_join(node_data[, c("node", "class")],
                   by = c("to" = "node")) %>%
  dplyr::rename(to_class = class) %>%
  dplyr::filter(from_class == "Oral microbiome" &
                  to_class == "Metabolite") %>%
  dplyr::select(to) %>%
  dplyr::distinct(.keep_all = TRUE) %>%
  dplyr::left_join(metabolome_variable_info,
                   by = c("to" = "variable_id"))

nasal_metabolome =
  edge_data %>%
  dplyr::left_join(node_data[, c("node", "class")],
                   by = c("from" = "node")) %>%
  dplyr::rename(from_class = class) %>%
  dplyr::left_join(node_data[, c("node", "class")],
                   by = c("to" = "node")) %>%
  dplyr::rename(to_class = class) %>%
  dplyr::filter(from_class == "Nasal microbiome" &
                  to_class == "Metabolite") %>%
  dplyr::select(to) %>%
  dplyr::distinct(.keep_all = TRUE) %>%
  dplyr::left_join(metabolome_variable_info,
                   by = c("to" = "variable_id"))

table(stool_metabolome$Super.pathway)
table(skin_metabolome$Super.pathway)
table(oral_metabolome$Super.pathway)
table(nasal_metabolome$Super.pathway)

#####mosaic plot to show the distributation
df = 
  rbind(
    data.frame(stool_metabolome[,c("to", "Super.pathway")], type = "Stool") %>% dplyr::filter(!is.na(Super.pathway)),
    data.frame(skin_metabolome[,c("to", "Super.pathway")], type = "Skin") %>% dplyr::filter(!is.na(Super.pathway)),
    data.frame(oral_metabolome[,c("to", "Super.pathway")], type = "Oral") %>% dplyr::filter(!is.na(Super.pathway)),
    data.frame(nasal_metabolome[,c("to", "Super.pathway")], type = "Nasal") %>% dplyr::filter(!is.na(Super.pathway))
  ) %>% 
  dplyr::mutate(type = factor(type, levels = c("Stool", "Skin", "Oral", "Nasal")))
  
plot =
  ggplot(df) +
  geom_mosaic(aes(x = product(Super.pathway, type), fill = Super.pathway),
              show.legend = TRUE) +
  scale_fill_manual(values = metabolite_class_color) +
  labs(x = "", y = "") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = "transparent"),
    panel.background = element_rect(fill = "transparent")
  )

plot

temp = 
  df %>% 
  dplyr::group_by(type, Super.pathway) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::mutate(prop = round(n*100/sum(n), 2))

text_data =
  ggplot_build(plot)$data %>%
  as.data.frame %>%
  transmute(
    x__type,
    x__fill__Super.pathway,
    x.position = (xmax + xmin) / 2,
    y.position = (ymax + ymin) / 2
  ) %>%
  dplyr::rename(type = x__type,
                Super.pathway = x__fill__Super.pathway) %>%
  dplyr::left_join(temp, by = c("type", "Super.pathway"))

plot = 
  plot +
  geom_text(aes(x.position, y.position, label = prop),
            data = text_data,
            color = "white")
plot

# ggsave(plot,
#        filename = "pathway_enrichment/metabolome/metabolite_distributation.pdf",
#        width = 11,
#        height = 7)


#####lipid pathway enrichment

########Lipidome
stool_lipidome =
  edge_data %>%
  dplyr::left_join(node_data[, c("node", "class")],
                   by = c("from" = "node")) %>%
  dplyr::rename(from_class = class) %>%
  dplyr::left_join(node_data[, c("node", "class")],
                   by = c("to" = "node")) %>%
  dplyr::rename(to_class = class) %>%
  dplyr::filter(from_class == "Stool microbiome" &
                  to_class == "Lipidome") %>%
  dplyr::select(to, direction) %>%
  dplyr::distinct(.keep_all = TRUE) %>%
  dplyr::left_join(lipidome_variable_info,
                   by = c("to" = "variable_id"))

skin_lipidome =
  edge_data %>%
  dplyr::left_join(node_data[, c("node", "class")],
                   by = c("from" = "node")) %>%
  dplyr::rename(from_class = class) %>%
  dplyr::left_join(node_data[, c("node", "class")],
                   by = c("to" = "node")) %>%
  dplyr::rename(to_class = class) %>%
  dplyr::filter(from_class == "Skin microbiome" &
                  to_class == "Lipidome") %>%
  dplyr::select(to) %>%
  dplyr::distinct(.keep_all = TRUE) %>%
  dplyr::left_join(lipidome_variable_info,
                   by = c("to" = "variable_id"))

oral_lipidome =
  edge_data %>%
  dplyr::left_join(node_data[, c("node", "class")],
                   by = c("from" = "node")) %>%
  dplyr::rename(from_class = class) %>%
  dplyr::left_join(node_data[, c("node", "class")],
                   by = c("to" = "node")) %>%
  dplyr::rename(to_class = class) %>%
  dplyr::filter(from_class == "Oral microbiome" &
                  to_class == "Lipidome") %>%
  dplyr::select(to) %>%
  dplyr::distinct(.keep_all = TRUE) %>%
  dplyr::left_join(lipidome_variable_info,
                   by = c("to" = "variable_id"))

nasal_lipidome =
  edge_data %>%
  dplyr::left_join(node_data[, c("node", "class")],
                   by = c("from" = "node")) %>%
  dplyr::rename(from_class = class) %>%
  dplyr::left_join(node_data[, c("node", "class")],
                   by = c("to" = "node")) %>%
  dplyr::rename(to_class = class) %>%
  dplyr::filter(from_class == "Nasal microbiome" &
                  to_class == "Lipidome") %>%
  dplyr::select(to) %>%
  dplyr::distinct(.keep_all = TRUE) %>%
  dplyr::left_join(lipidome_variable_info,
                   by = c("to" = "variable_id"))

table(unlist(stringr::str_split(stool_lipidome$LipidClass, ';')))
table(unlist(stringr::str_split(skin_lipidome$LipidClass, ';')))
table(unlist(stringr::str_split(oral_lipidome$LipidClass, ';')))
table(unlist(stringr::str_split(nasal_lipidome$LipidClass, ';')))

#####mosaic plot to show the distributation
df =
  rbind(
    data.frame(stool_lipidome[, c("to", "LipidClass")], type = "Stool") %>% dplyr::filter(!is.na(LipidClass)),
    data.frame(skin_lipidome[, c("to", "LipidClass")], type = "Skin") %>% dplyr::filter(!is.na(LipidClass)),
    data.frame(oral_lipidome[, c("to", "LipidClass")], type = "Oral") %>% dplyr::filter(!is.na(LipidClass)),
    data.frame(nasal_lipidome[, c("to", "LipidClass")], type = "Nasal") %>% dplyr::filter(!is.na(LipidClass))
  ) %>% 
  tidyr::separate_rows(LipidClass, sep = ";") %>% 
  dplyr::mutate(type = factor(type, levels = c("Stool", "Skin", "Oral", "Nasal")))

plot =
  ggplot(df) +
  geom_mosaic(aes(x = product(LipidClass, type), fill = LipidClass),
              show.legend = TRUE) +
  scale_fill_manual(values = lipid_class_color) +
  labs(x = "", y = "") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = "transparent"),
    panel.background = element_rect(fill = "transparent")
  )

plot

temp = 
  df %>% 
  dplyr::group_by(type, LipidClass) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::mutate(prop = round(n*100/sum(n), 2))

text_data =
  ggplot_build(plot)$data %>%
  as.data.frame %>%
  transmute(
    x__type,
    x__fill__LipidClass,
    x.position = (xmax + xmin) / 2,
    y.position = (ymax + ymin) / 2
  ) %>%
  dplyr::rename(type = x__type,
                LipidClass = x__fill__LipidClass) %>%
  dplyr::left_join(temp, by = c("type", "LipidClass"))

plot = 
  plot +
  geom_text(aes(x.position, y.position, label = prop),
            data = text_data,
            color = "white")
plot

# ggsave(plot,
#        filename = "pathway_enrichment/lipidome/lipid_distributation.pdf",
#        width = 11,
#        height = 7)


edge_data %>%
  dplyr::left_join(node_data[, c("node", "class")],
                   by = c("from" = "node")) %>%
  dplyr::rename(from_class = class) %>%
  dplyr::left_join(node_data[, c("node", "class")],
                   by = c("to" = "node")) %>%
  dplyr::rename(to_class = class) %>%
  dplyr::filter(from_class == "Skin microbiome" &
                  to_class == "Lipidome") %>%
  dplyr::select(to, direction) %>%
  dplyr::left_join(lipidome_variable_info[,c("variable_id", "LipidClass")],
                   by = c("to" = "variable_id")) %>% 
  tidyr::separate_rows(LipidClass, sep = ";") %>% 
  dplyr::group_by(LipidClass, direction) %>% 
  dplyr::summarise(n = n())


######lipid minion analysis
dir.create("pathway_enrichment/lipidome/lipidminion")
###all lipids
all_lipid =
  lipidome_variable_info %>% 
  dplyr::select(Lipid_Name) %>% 
  tidyr::separate_rows(Lipid_Name, sep = ";") %>% 
  dplyr::mutate(
    Lipid_Name =  stringr::str_replace_all(Lipid_Name, "\\_", "\\/") %>%
      stringr::str_replace_all("TAG", "TG") %>%
      stringr::str_replace_all("DAG", "DG")
  )

all_lipid$Lipid_Name[grep("TG", all_lipid$Lipid_Name)] =
  all_lipid$Lipid_Name[grep("TG", all_lipid$Lipid_Name)] %>%
  purrr::map(function(x) {
    stringr::str_extract(x, "TG[0-9]{1,3}\\.[0-9]{1,2}") %>%
      stringr::str_replace("TG", "") %>%
      stringr::str_split("\\.") %>%
      `[[`(1) %>%
      paste(collapse = ":") %>%
      paste("TG(", ., ")", sep = "")
  }) %>%
  unlist()

all_lipid =
  all_lipid %>%
  dplyr::distinct() %>% 
  dplyr::rename(lipid = Lipid_Name)
  
###stool 
stool_lipid =
  stool_lipidome %>%
  dplyr::select(Lipid_Name) %>%
  tidyr::separate_rows(Lipid_Name, sep = ";") %>%
  dplyr::mutate(
    Lipid_Name =  stringr::str_replace_all(Lipid_Name, "\\_", "\\/") %>%
      stringr::str_replace_all("TAG", "TG") %>%
      stringr::str_replace_all("DAG", "DG")
  )

stool_lipid$Lipid_Name[grep("TG", stool_lipid$Lipid_Name)] =
  stool_lipid$Lipid_Name[grep("TG", stool_lipid$Lipid_Name)] %>%
  purrr::map(function(x) {
    stringr::str_extract(x, "TG[0-9]{1,3}\\.[0-9]{1,2}") %>%
      stringr::str_replace("TG", "") %>%
      stringr::str_split("\\.") %>%
      `[[`(1) %>%
      paste(collapse = ":") %>%
      paste("TG(", ., ")", sep = "")
  }) %>%
  unlist()

stool_lipid =
  stool_lipid %>%
  dplyr::distinct() %>% 
  dplyr::rename(lipid = Lipid_Name)

# write.table(
#   all_lipid,
#   file = "pathway_enrichment/lipidome/stool/all_lipid.txt",
#   row.names = FALSE,
#   col.names = TRUE,
#   quote = FALSE
# )
# 
# write.table(
#   stool_lipid,
#   file = "pathway_enrichment/lipidome/stool/stool_lipid.txt",
#   row.names = FALSE,
#   col.names = TRUE,
#   quote = FALSE
# )







###skin 
skin_lipid =
  skin_lipidome %>%
  dplyr::select(Lipid_Name) %>%
  tidyr::separate_rows(Lipid_Name, sep = ";") %>%
  dplyr::mutate(
    Lipid_Name =  stringr::str_replace_all(Lipid_Name, "\\_", "\\/") %>%
      stringr::str_replace_all("TAG", "TG") %>%
      stringr::str_replace_all("DAG", "DG")
  )

skin_lipid$Lipid_Name[grep("TG", skin_lipid$Lipid_Name)] =
  skin_lipid$Lipid_Name[grep("TG", skin_lipid$Lipid_Name)] %>%
  purrr::map(function(x) {
    stringr::str_extract(x, "TG[0-9]{1,3}\\.[0-9]{1,2}") %>%
      stringr::str_replace("TG", "") %>%
      stringr::str_split("\\.") %>%
      `[[`(1) %>%
      paste(collapse = ":") %>%
      paste("TG(", ., ")", sep = "")
  }) %>%
  unlist()

skin_lipid =
  skin_lipid %>%
  dplyr::distinct() %>% 
  dplyr::rename(lipid = Lipid_Name)

# write.table(
#   all_lipid,
#   file = "pathway_enrichment/lipidome/skin/all_lipid.txt",
#   row.names = FALSE,
#   col.names = TRUE,
#   quote = FALSE
# )
# 
# write.table(
#   skin_lipid,
#   file = "pathway_enrichment/lipidome/skin/skin_lipid.txt",
#   row.names = FALSE,
#   col.names = TRUE,
#   quote = FALSE
# )


###then open the lipidminon website https://omicstools.pnnl.gov/shiny/lipid-mini-on/


###no important classifications are found


