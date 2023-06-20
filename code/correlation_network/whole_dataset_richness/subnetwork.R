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

setwd(masstools::get_project_wd())
library(tidyverse)
rm(list = ls())

source("code/tools.R")

######work directory
setwd(masstools::get_project_wd())
setwd("data_analysis/correlation_network/whole_data_set_enrichness/inter_omics_correlation_network")

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
        "data_analysis/correlation_network/whole_data_set_enrichness/stool_microbiome_vs_metabolome/stool_microbiome_vs_metabolome.xlsx"
      ),
      sheet = 2
    )
  stool_microbiome_metabolome_edge$from = paste("stool", stool_microbiome_metabolome_edge$from, sep = "_")
  
  stool_microbiome_metabolome_node =
    readxl::read_xlsx(
      here::here(
        "data_analysis/correlation_network/whole_data_set_enrichness/stool_microbiome_vs_metabolome/stool_microbiome_vs_metabolome.xlsx"
      ),
      sheet = 1
    )
  
  stool_microbiome_metabolome_node$node[grep("microbiome", stool_microbiome_metabolome_node$class)] =
    paste("stool", stool_microbiome_metabolome_node$node[grep("microbiome", stool_microbiome_metabolome_node$class)],
          sep = "_")
  
  
  stool_microbiome_lipidome_edge =
    readxl::read_xlsx(
      here::here(
        "data_analysis/correlation_network/whole_data_set_enrichness/stool_microbiome_vs_lipidome/stool_microbiome_vs_lipidome.xlsx"
      ),
      sheet = 2
    )
  stool_microbiome_lipidome_edge$from = paste("stool", stool_microbiome_lipidome_edge$from, sep = "_")
  
  stool_microbiome_lipidome_node =
    readxl::read_xlsx(
      here::here(
        "data_analysis/correlation_network/whole_data_set_enrichness/stool_microbiome_vs_lipidome/stool_microbiome_vs_lipidome.xlsx"
      ),
      sheet = 1
    )
  
  stool_microbiome_lipidome_node$node[grep("microbiome", stool_microbiome_lipidome_node$class)] =
    paste("stool", stool_microbiome_lipidome_node$node[grep("microbiome", stool_microbiome_lipidome_node$class)],
          sep = "_")
  
  
  stool_microbiome_proteome_edge =
    readxl::read_xlsx(
      here::here(
        "data_analysis/correlation_network/whole_data_set_enrichness/stool_microbiome_vs_proteome/stool_microbiome_vs_proteome.xlsx"
      ),
      sheet = 2
    )
  stool_microbiome_proteome_edge$from = paste("stool", stool_microbiome_proteome_edge$from, sep = "_")
  
  stool_microbiome_proteome_node =
    readxl::read_xlsx(
      here::here(
        "data_analysis/correlation_network/whole_data_set_enrichness/stool_microbiome_vs_proteome/stool_microbiome_vs_proteome.xlsx"
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
        "data_analysis/correlation_network/whole_data_set_enrichness/skin_microbiome_vs_metabolome/skin_microbiome_vs_metabolome.xlsx"
      ),
      sheet = 2
    )
  skin_microbiome_metabolome_edge$from = paste("skin", skin_microbiome_metabolome_edge$from, sep = "_")
  
  skin_microbiome_metabolome_node =
    readxl::read_xlsx(
      here::here(
        "data_analysis/correlation_network/whole_data_set_enrichness/skin_microbiome_vs_metabolome/skin_microbiome_vs_metabolome.xlsx"
      ),
      sheet = 1
    )
  
  skin_microbiome_metabolome_node$node[grep("microbiome", skin_microbiome_metabolome_node$class)] =
    paste("skin", skin_microbiome_metabolome_node$node[grep("microbiome", skin_microbiome_metabolome_node$class)],
          sep = "_")
  
  
  skin_microbiome_lipidome_edge =
    readxl::read_xlsx(
      here::here(
        "data_analysis/correlation_network/whole_data_set_enrichness/skin_microbiome_vs_lipidome/skin_microbiome_vs_lipidome.xlsx"
      ),
      sheet = 2
    )
  skin_microbiome_lipidome_edge$from = paste("skin", skin_microbiome_lipidome_edge$from, sep = "_")
  
  skin_microbiome_lipidome_node =
    readxl::read_xlsx(
      here::here(
        "data_analysis/correlation_network/whole_data_set_enrichness/skin_microbiome_vs_lipidome/skin_microbiome_vs_lipidome.xlsx"
      ),
      sheet = 1
    )
  
  skin_microbiome_lipidome_node$node[grep("microbiome", skin_microbiome_lipidome_node$class)] =
    paste("skin", skin_microbiome_lipidome_node$node[grep("microbiome", skin_microbiome_lipidome_node$class)],
          sep = "_")
  
  skin_microbiome_proteome_edge =
    readxl::read_xlsx(
      here::here(
        "data_analysis/correlation_network/whole_data_set_enrichness/skin_microbiome_vs_proteome/skin_microbiome_vs_proteome.xlsx"
      ),
      sheet = 2
    )
  skin_microbiome_proteome_edge$from = paste("skin", skin_microbiome_proteome_edge$from, sep = "_")
  
  skin_microbiome_proteome_node =
    readxl::read_xlsx(
      here::here(
        "data_analysis/correlation_network/whole_data_set_enrichness/skin_microbiome_vs_proteome/skin_microbiome_vs_proteome.xlsx"
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
        "data_analysis/correlation_network/whole_data_set_enrichness/nasal_microbiome_vs_metabolome/nasal_microbiome_vs_metabolome.xlsx"
      ),
      sheet = 2
    )
  nasal_microbiome_metabolome_edge$from = paste("nasal", nasal_microbiome_metabolome_edge$from, sep = "_")
  
  nasal_microbiome_metabolome_node =
    readxl::read_xlsx(
      here::here(
        "data_analysis/correlation_network/whole_data_set_enrichness/nasal_microbiome_vs_metabolome/nasal_microbiome_vs_metabolome.xlsx"
      ),
      sheet = 1
    )
  
  nasal_microbiome_metabolome_node$node[grep("microbiome", nasal_microbiome_metabolome_node$class)] =
    paste("nasal", nasal_microbiome_metabolome_node$node[grep("microbiome", nasal_microbiome_metabolome_node$class)],
          sep = "_")
  
  
  nasal_microbiome_lipidome_edge =
    readxl::read_xlsx(
      here::here(
        "data_analysis/correlation_network/whole_data_set_enrichness/nasal_microbiome_vs_lipidome/nasal_microbiome_vs_lipidome.xlsx"
      ),
      sheet = 2
    )
  nasal_microbiome_lipidome_edge$from = paste("nasal", nasal_microbiome_lipidome_edge$from, sep = "_")
  
  nasal_microbiome_lipidome_node =
    readxl::read_xlsx(
      here::here(
        "data_analysis/correlation_network/whole_data_set_enrichness/nasal_microbiome_vs_lipidome/nasal_microbiome_vs_lipidome.xlsx"
      ),
      sheet = 1
    )
  
  nasal_microbiome_lipidome_node$node[grep("microbiome", nasal_microbiome_lipidome_node$class)] =
    paste("nasal", nasal_microbiome_lipidome_node$node[grep("microbiome", nasal_microbiome_lipidome_node$class)],
          sep = "_")
  
  
  nasal_microbiome_proteome_edge =
    readxl::read_xlsx(
      here::here(
        "data_analysis/correlation_network/whole_data_set_enrichness/nasal_microbiome_vs_proteome/nasal_microbiome_vs_proteome.xlsx"
      ),
      sheet = 2
    )
  nasal_microbiome_proteome_edge$from = paste("nasal", nasal_microbiome_proteome_edge$from, sep = "_")
  
  nasal_microbiome_proteome_node =
    readxl::read_xlsx(
      here::here(
        "data_analysis/correlation_network/whole_data_set_enrichness/nasal_microbiome_vs_proteome/nasal_microbiome_vs_proteome.xlsx"
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
        "data_analysis/correlation_network/whole_data_set_enrichness/oral_microbiome_vs_metabolome/oral_microbiome_vs_metabolome.xlsx"
      ),
      sheet = 2
    )
  oral_microbiome_metabolome_edge$from = paste("oral", oral_microbiome_metabolome_edge$from, sep = "_")
  
  oral_microbiome_metabolome_node =
    readxl::read_xlsx(
      here::here(
        "data_analysis/correlation_network/whole_data_set_enrichness/oral_microbiome_vs_metabolome/oral_microbiome_vs_metabolome.xlsx"
      ),
      sheet = 1
    )
  
  oral_microbiome_metabolome_node$node[grep("microbiome", oral_microbiome_metabolome_node$class)] =
    paste("oral", oral_microbiome_metabolome_node$node[grep("microbiome", oral_microbiome_metabolome_node$class)],
          sep = "_")
  
  
  oral_microbiome_lipidome_edge =
    readxl::read_xlsx(
      here::here(
        "data_analysis/correlation_network/whole_data_set_enrichness/oral_microbiome_vs_lipidome/oral_microbiome_vs_lipidome.xlsx"
      ),
      sheet = 2
    )
  oral_microbiome_lipidome_edge$from = paste("oral", oral_microbiome_lipidome_edge$from, sep = "_")
  
  oral_microbiome_lipidome_node =
    readxl::read_xlsx(
      here::here(
        "data_analysis/correlation_network/whole_data_set_enrichness/oral_microbiome_vs_lipidome/oral_microbiome_vs_lipidome.xlsx"
      ),
      sheet = 1
    )
  
  oral_microbiome_lipidome_node$node[grep("microbiome", oral_microbiome_lipidome_node$class)] =
    paste("oral", oral_microbiome_lipidome_node$node[grep("microbiome", oral_microbiome_lipidome_node$class)],
          sep = "_")
  
  
  
  
  oral_microbiome_proteome_edge =
    readxl::read_xlsx(
      here::here(
        "data_analysis/correlation_network/whole_data_set_enrichness/oral_microbiome_vs_proteome/oral_microbiome_vs_proteome.xlsx"
      ),
      sheet = 2
    )
  oral_microbiome_proteome_edge$from = paste("oral", oral_microbiome_proteome_edge$from, sep = "_")
  
  oral_microbiome_proteome_node =
    readxl::read_xlsx(
      here::here(
        "data_analysis/correlation_network/whole_data_set_enrichness/oral_microbiome_vs_proteome/oral_microbiome_vs_proteome.xlsx"
      ),
      sheet = 1
    )
  
  oral_microbiome_proteome_node$node[grep("microbiome", oral_microbiome_proteome_node$class)] =
    paste("oral", oral_microbiome_proteome_node$node[grep("microbiome", oral_microbiome_proteome_node$class)],
          sep = "_")
  
}


###output plot
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

edge_data =
  edge_data %>%
  dplyr::filter(from %in% node_data$node & to %in% node_data$node)

node_data =
  node_data %>%
  dplyr::filter(node %in% edge_data$from | node %in% edge_data$to)

dim(edge_data)
dim(node_data)

all_graph <- 
  tidygraph::tbl_graph(nodes = node_data, 
                       edges = edge_data,
                       directed = FALSE) %>% 
  dplyr::mutate(Degree = centrality_degree(mode = 'all'))

######some subnetwork
###ethyl glucuronide
temp_edge_data = 
  edge_data %>% 
  dplyr::filter(to_true_name == "ethyl glucuronide")

temp_node_data = 
  node_data %>% 
  dplyr::filter(node %in% c(temp_edge_data$from, temp_edge_data$to))

all_graph <-
  tidygraph::tbl_graph(nodes = node_data,
                       edges = edge_data,
                       directed = FALSE) %>%
  dplyr::mutate(Degree = centrality_degree(mode = 'all'))

# igraph::subgraph(graph = all_graph, v = temp_node_data$node)

temp_graph =
  all_graph %>%
  tidygraph::activate(what = "nodes") %>%
  tidygraph::filter(node %in% temp_node_data$node)

g <- temp_graph
library(igraph)
library(ggraph)

V(g)$type <- bipartite_mapping(g)$type

coords <-
  create_layout(g, layout = "bipartite") %>%
  dplyr::select(node, class, true_name, x, y)

coords =
  coords %>%
  dplyr::left_join(
    edge_data %>% dplyr::filter(to_true_name == "ethyl glucuronide") %>% dplyr::select(from, cor, p_adjust),
    by = c("node" = "from")
  )

coords$cor[is.na(coords$cor)] = 0

coords <-
  coords %>%
  dplyr::mutate(y1 = x) %>%
  dplyr::mutate(x1 = case_when(
    y == 1 ~ 0,
    y == 0 ~ 1
  )) %>%
  dplyr::select(-c(x, y)) %>%
  dplyr::select(x = x1, y = y1, everything())

coords$x[coords$cor < 0] = -1
coords$x[coords$cor > 0] = 1
coords$x[coords$class == "Metabolite"] = 0

my_graph <-
  create_layout(
    graph = g,
    layout = "manual",
    x = coords$x,
    y = coords$y
    # node.position = coords
  )

plot = 
  ggraph(temp_graph, 
         layout = "manual", 
         x = coords$x,
         y = coords$y) +
  geom_edge_diagonal(aes(
    color = direction,
    width = -log(p_adjust, 10)
  ),
  alpha = 1,
  show.legend = TRUE) +
  ggstar::geom_star(aes(
    x = x,
    y = y,
    size = Degree,
    fill = class,
    starshape = class
  )) +
  ggstar::scale_starshape_manual(values = omics_shape) +
  shadowtext::geom_shadowtext(
    aes(
      x = x,
      y = y,
      label = true_name,
      color = class
    ),
    bg.color = "white",
    size = 5,
    show.legend = FALSE
  ) +
  ggraph::theme_graph() +
  guides(edge_color = guide_edge_colorbar(title = "Correlation"),
         color = guide_legend(title = "Class")) +
  scale_edge_width_continuous(range = c(0.3, 1)) +
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

# ggsave(plot,
#        filename = "ethyl glucuronide_subnetwork.pdf",
#        width = 7,
#        height = 7)



####
###Akkermansia subnetwork
temp_edge_data = 
  edge_data %>% 
  dplyr::filter(from_true_name == "Akkermansia")

temp_node_data = 
  node_data %>% 
  dplyr::filter(node %in% c(temp_edge_data$from, temp_edge_data$to))

all_graph <-
  tidygraph::tbl_graph(nodes = node_data,
                       edges = edge_data,
                       directed = FALSE) %>%
  dplyr::mutate(Degree = centrality_degree(mode = 'all'))

temp_graph =
  all_graph %>%
  tidygraph::activate(what = "nodes") %>%
  tidygraph::filter(node %in% temp_node_data$node)

g <- temp_graph
library(igraph)
library(ggraph)

V(g)$type <- bipartite_mapping(g)$type

coords <-
  create_layout(g, layout = "bipartite") %>%
  dplyr::select(node, class, true_name, x, y)

coords =
  coords %>%
  dplyr::left_join(
    edge_data %>% dplyr::filter(from_true_name == "Akkermansia") %>% dplyr::select(to, cor, p_adjust),
    by = c("node" = "to")
  )

coords$cor[is.na(coords$cor)] = 0

coords <-
  coords %>%
  dplyr::mutate(y1 = x) %>%
  dplyr::mutate(x1 = case_when(
    y == 1 ~ 0,
    y == 0 ~ 1
  )) %>%
  dplyr::select(-c(x, y)) %>%
  dplyr::select(x = x1, y = y1, everything())

coords$x[coords$cor < 0] = -1
coords$x[coords$cor > 0] = 1
coords$x[coords$class == "Stool microbiome"] = 0

my_graph <-
  create_layout(
    graph = g,
    layout = "manual",
    x = coords$x,
    y = coords$y
    # node.position = coords
  )

plot = 
  ggraph(temp_graph, 
         layout = "manual", 
         x = coords$x,
         y = coords$y) +
  geom_edge_diagonal(aes(
    color = direction,
    width = -log(p_adjust, 10)
  ),
  alpha = 1,
  show.legend = TRUE) +
  ggstar::geom_star(aes(
    x = x,
    y = y,
    size = Degree,
    fill = class,
    starshape = class
  )) +
  ggstar::scale_starshape_manual(values = omics_shape) +
  shadowtext::geom_shadowtext(
    aes(
      x = x,
      y = y,
      label = true_name,
      color = class
    ),
    bg.color = "white",
    size = 5,
    show.legend = FALSE
  ) +
  ggraph::theme_graph() +
  guides(edge_color = guide_edge_colorbar(title = "Correlation"),
         color = guide_legend(title = "Class")) +
  scale_edge_width_continuous(range = c(0.3, 1)) +
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

# ggsave(plot,
#        filename = "Akkermansia_subnetwork.pdf",
#        width = 7,
#        height = 7)


###gamma-glutamylthreonine
temp_edge_data = 
  edge_data %>% 
  dplyr::filter(to_true_name == "gamma-glutamylthreonine")

temp_node_data = 
  node_data %>% 
  dplyr::filter(node %in% c(temp_edge_data$from, temp_edge_data$to))

temp_graph =
  all_graph %>%
  tidygraph::activate(what = "nodes") %>%
  tidygraph::filter(node %in% temp_node_data$node)

g <- temp_graph
library(igraph)
library(ggraph)

V(g)$type <- bipartite_mapping(g)$type

coords <-
  create_layout(g, layout = "bipartite") %>%
  dplyr::select(node, class, true_name, x, y)

coords =
  coords %>%
  dplyr::left_join(
    edge_data %>% dplyr::filter(to_true_name == "gamma-glutamylthreonine") %>% dplyr::select(from, cor, p_adjust),
    by = c("node" = "from")
  )

coords$cor[is.na(coords$cor)] = 0

coords <-
  coords %>%
  dplyr::mutate(y1 = x) %>%
  dplyr::mutate(x1 = case_when(
    y == 1 ~ 0,
    y == 0 ~ 1
  )) %>%
  dplyr::select(-c(x, y)) %>%
  dplyr::select(x = x1, y = y1, everything())

coords$x[coords$cor < 0] = -1
coords$x[coords$cor > 0] = 1
coords$x[coords$class == "Metabolite"] = 0

my_graph <-
  create_layout(
    graph = g,
    layout = "manual",
    x = coords$x,
    y = coords$y
    # node.position = coords
  )

plot = 
  ggraph(temp_graph, 
         layout = "manual", 
         x = coords$x,
         y = coords$y) +
  geom_edge_diagonal(aes(
    color = direction,
    width = -log(p_adjust, 10)
  ),
  alpha = 1,
  show.legend = TRUE) +
  ggstar::geom_star(aes(
    x = x,
    y = y,
    size = Degree,
    fill = class,
    starshape = class
  )) +
  ggstar::scale_starshape_manual(values = omics_shape) +
  shadowtext::geom_shadowtext(
    aes(
      x = x,
      y = y,
      label = true_name,
      color = class
    ),
    bg.color = "white",
    size = 5,
    show.legend = FALSE
  ) +
  ggraph::theme_graph() +
  guides(edge_color = guide_edge_colorbar(title = "Correlation"),
         color = guide_legend(title = "Class")) +
  scale_edge_width_continuous(range = c(0.3, 1)) +
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

# ggsave(plot,
#        filename = "gamma-glutamylthreonine_subnetwork.pdf",
#        width = 7,
#        height = 7)



###peptide subnetwork
peptide_name =
  metabolome_variable_info$Metabolite[which(metabolome_variable_info$Super.pathway == "Peptide")]

temp_edge_data = 
  edge_data %>% 
  dplyr::filter(to_true_name %in% peptide_name)

temp_node_data = 
  node_data %>% 
  dplyr::filter(node %in% c(temp_edge_data$from, temp_edge_data$to))

temp_graph =
  all_graph %>%
  tidygraph::activate(what = "nodes") %>%
  tidygraph::filter(node %in% temp_node_data$node)

plot = 
  ggraph(temp_graph, 
         layout = "kk") +
  geom_edge_link(aes(
    color = pos_neg,
    width = -log(p_adjust, 10)
  ),
  alpha = 1,
  show.legend = TRUE) +
  ggstar::geom_star(aes(
    x = x,
    y = y,
    size = Degree,
    fill = class,
    starshape = class
  )) +
  ggstar::scale_starshape_manual(values = omics_shape) +
  shadowtext::geom_shadowtext(
    aes(
      x = x,
      y = y,
      label = ifelse(class == "Metabolite", true_name, NA),
      color = class
    ),
    bg.color = "white",
    size = 3,
    show.legend = FALSE
  ) +
  ggraph::theme_graph() +
  guides(edge_color = guide_edge_colorbar(title = "Correlation"),
         color = guide_legend(title = "Class")) +
  scale_edge_width_continuous(range = c(0.3, 1)) +
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

# ggsave(plot,
#        filename = "peptide_subnetwork.pdf",
#        width = 7,
#        height = 7)




###p-Cresol sulfate and p-Cresol glucuronide
temp_edge_data = 
  edge_data %>% 
  dplyr::filter(to_true_name %in% c("p-Cresol sulfate", "p-Cresol glucuronide"))

temp_node_data = 
  node_data %>% 
  dplyr::filter(node %in% c(temp_edge_data$from, temp_edge_data$to))

temp_graph =
  all_graph %>%
  tidygraph::activate(what = "nodes") %>%
  tidygraph::filter(node %in% temp_node_data$node)

plot = 
  ggraph(temp_graph, 
         layout = "kk") +
  geom_edge_link(aes(
    color = pos_neg,
    width = -log(p_adjust, 10)
  ),
  alpha = 1,
  show.legend = TRUE) +
  ggstar::geom_star(aes(
    x = x,
    y = y,
    size = Degree,
    fill = class,
    starshape = class
  )) +
  ggstar::scale_starshape_manual(values = omics_shape) +
  shadowtext::geom_shadowtext(
    aes(
      x = x,
      y = y,
      label = true_name,
      color = class
    ),
    bg.color = "white",
    size = 3,
    show.legend = FALSE
  ) +
  ggraph::theme_graph() +
  guides(edge_color = guide_edge_colorbar(title = "Correlation"),
         color = guide_legend(title = "Class")) +
  scale_edge_width_continuous(range = c(0.3, 1)) +
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

# ggsave(plot,
#        filename = "p-Cresol_subnetwork.pdf",
#        width = 7,
#        height = 7)


