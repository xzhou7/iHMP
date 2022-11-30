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
setwd(
  "data_analysis/correlation_network/whole_data_set/inter_omics_correlation_network"
)

####load data
###skin microbiome
{
  load(here::here(
    "data_analysis/skin_microbiome/data_preparation/expression_data"
  ))
  load(here::here(
    "data_analysis/skin_microbiome/data_preparation/sample_info"
  ))
  load(here::here(
    "data_analysis/skin_microbiome/data_preparation/variable_info"
  ))
  
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
    skin_microbiome_variable_info[match(rownames(expression_data),
                                        skin_microbiome_variable_info$Genus), ]
  
  skin_microbiome_variable_info$Genus == rownames(expression_data)
  
  ###remove the variables which Genus are NA
  remove_idx = which(is.na(skin_microbiome_variable_info$Genus))
  remove_idx
  if (length(remove_idx) > 0) {
    skin_microbiome_variable_info = skin_microbiome_variable_info[-remove_idx, ]
    expression_data = expression_data[-remove_idx, ]
  }
  
  rownames(expression_data) = skin_microbiome_variable_info$variable_id
  
  skin_microbiome_expression_data =
    expression_data
  
  skin_microbiome_variable_info =
    skin_microbiome_variable_info %>%
    dplyr::filter(!stringr::str_detect(Genus, "Unclassified_Bacteria"))
  
  skin_microbiome_expression_data =
    skin_microbiome_expression_data[skin_microbiome_variable_info$variable_id, ]
  
  dim(skin_microbiome_sample_info)
  dim(skin_microbiome_variable_info)
  dim(skin_microbiome_expression_data)
  
  rownames(skin_microbiome_expression_data) == skin_microbiome_variable_info$variable_id
  colnames(skin_microbiome_expression_data) == skin_microbiome_sample_info$sample_id
}


###stool microbiome
{
  load(
    here::here(
      "data_analysis/stool_microbiome/data_preparation/expression_data"
    )
  )
  load(here::here(
    "data_analysis/stool_microbiome/data_preparation/sample_info"
  ))
  load(here::here(
    "data_analysis/stool_microbiome/data_preparation/variable_info"
  ))
  
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
    stool_microbiome_variable_info[match(rownames(expression_data),
                                         stool_microbiome_variable_info$Genus), ]
  
  stool_microbiome_variable_info$Genus == rownames(expression_data)
  
  rownames(expression_data) = stool_microbiome_variable_info$variable_id
  
  stool_microbiome_expression_data =
    expression_data
  
  stool_microbiome_variable_info =
    stool_microbiome_variable_info %>%
    dplyr::filter(!stringr::str_detect(Genus, "Unclassified_Bacteria"))
  
  stool_microbiome_expression_data =
    stool_microbiome_expression_data[stool_microbiome_variable_info$variable_id, ]
  
  dim(stool_microbiome_sample_info)
  dim(stool_microbiome_variable_info)
  
  rownames(stool_microbiome_expression_data) == stool_microbiome_variable_info$variable_id
  colnames(stool_microbiome_expression_data) == stool_microbiome_sample_info$sample_id
}

###nasal microbiome
{
  load(
    here::here(
      "data_analysis/nasal_microbiome/data_preparation/expression_data"
    )
  )
  load(here::here(
    "data_analysis/nasal_microbiome/data_preparation/sample_info"
  ))
  load(here::here(
    "data_analysis/nasal_microbiome/data_preparation/variable_info"
  ))
  
  
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
    nasal_microbiome_variable_info[match(rownames(expression_data),
                                         nasal_microbiome_variable_info$Genus), ]
  
  nasal_microbiome_variable_info$Genus == rownames(expression_data)
  
  ###remove the variables which Genus are NA
  remove_idx = which(is.na(nasal_microbiome_variable_info$Genus))
  remove_idx
  if (length(remove_idx) > 0) {
    nasal_microbiome_variable_info = nasal_microbiome_variable_info[-remove_idx, ]
    expression_data = expression_data[-remove_idx, ]
  }
  
  rownames(expression_data) = nasal_microbiome_variable_info$variable_id
  
  nasal_microbiome_expression_data =
    expression_data
  
  nasal_microbiome_variable_info =
    nasal_microbiome_variable_info %>%
    dplyr::filter(!stringr::str_detect(Genus, "Unclassified_Bacteria"))
  
  nasal_microbiome_expression_data =
    nasal_microbiome_expression_data[nasal_microbiome_variable_info$variable_id, ]
  
  dim(nasal_microbiome_sample_info)
  dim(nasal_microbiome_variable_info)
  
  rownames(nasal_microbiome_expression_data) == nasal_microbiome_variable_info$variable_id
  colnames(nasal_microbiome_expression_data) == nasal_microbiome_sample_info$sample_id
}

###stool microbiome
{
  load(here::here(
    "data_analysis/oral_microbiome/data_preparation/expression_data"
  ))
  load(here::here(
    "data_analysis/oral_microbiome/data_preparation/sample_info"
  ))
  load(here::here(
    "data_analysis/oral_microbiome/data_preparation/variable_info"
  ))
  
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
    oral_microbiome_variable_info[match(rownames(expression_data),
                                        oral_microbiome_variable_info$Genus), ]
  
  oral_microbiome_variable_info$Genus == rownames(expression_data)
  
  ###remove the variables which Genus are NA
  remove_idx = which(is.na(oral_microbiome_variable_info$Genus))
  remove_idx
  if (length(remove_idx) > 0) {
    oral_microbiome_variable_info = oral_microbiome_variable_info[-remove_idx, ]
    expression_data = expression_data[-remove_idx, ]
  }
  
  rownames(expression_data) = oral_microbiome_variable_info$variable_id
  
  oral_microbiome_expression_data =
    expression_data
  
  oral_microbiome_variable_info =
    oral_microbiome_variable_info %>%
    dplyr::filter(!stringr::str_detect(Genus, "Unclassified_Bacteria"))
  
  oral_microbiome_expression_data =
    oral_microbiome_expression_data[oral_microbiome_variable_info$variable_id, ]
  
  dim(oral_microbiome_sample_info)
  dim(oral_microbiome_variable_info)
  
  rownames(oral_microbiome_expression_data) == oral_microbiome_variable_info$variable_id
  colnames(oral_microbiome_expression_data) == oral_microbiome_sample_info$sample_id
}

###plasma lipidome
{
  load(here::here(
    "data_analysis/lipidome/data_preparation/new_expression_data"
  ))
  load(here::here(
    "data_analysis/lipidome/data_preparation/new_sample_info"
  ))
  load(here::here(
    "data_analysis/lipidome/data_preparation/new_variable_info"
  ))
  
  lipidome_expression_data = new_expression_data
  lipidome_sample_info = new_sample_info
  lipidome_variable_info = new_variable_info
  
  ##remove QC samples
  lipidome_sample_info =
    lipidome_sample_info %>%
    dplyr::filter(subject_id != "QC")
  
  lipidome_expression_data =
    lipidome_expression_data[, lipidome_sample_info$sample_id]
  
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
    skin_microbiome_expression_data[, intersect_sample_id]
  
  lipidome_expression_data =
    lipidome_expression_data[, intersect_sample_id]
  
  skin_microbiome_sample_info =
    skin_microbiome_sample_info[match(intersect_sample_id,
                                      skin_microbiome_sample_info$sample_id), ]
  
  lipidome_sample_info =
    lipidome_sample_info[match(intersect_sample_id, lipidome_sample_info$sample_id), ]
  
  length(unique(lipidome_sample_info$subject_id))
  
  sum(lipidome_sample_info$sample_id == skin_microbiome_sample_info$sample_id)
  sum(lipidome_sample_info$subject_id == skin_microbiome_sample_info$subject_id)
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
  )

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

edge_data =
  edge_data %>%
  dplyr::mutate(direction =
                  case_when(cor > 0 ~ "positive",
                            cor < 0 ~ "negative")) %>%
  dplyr::mutate(pos_neg = direction)

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

temp_data <-
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
      microbiome_class = node_data$class[match(z$from[1], node_data$node)],
      total_number = length(x),
      Proteome = sum(x == "Proteome") * 100 / length(x),
      Metabolite = sum(x == "Metabolite") * 100 / length(x),
      Lipidome = sum(x == "Lipidome") * 100 / length(x)
    )
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame() %>%
  dplyr::arrange(microbiome_class, desc(total_number)) %>%
  dplyr::mutate(
    direction = case_when(
      Proteome >= 50 ~ "Proteome",
      Metabolite >= 50 ~ "Metabolite",
      Lipidome >= 50 ~ "Lipidome",
      TRUE ~ "Mixed"
    )
  )

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
      molecular_class = node_data$class[match(z$to[1], node_data$node)],
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
  dplyr::mutate(
    direction = case_when(
      Stool >= 50 ~ "Stool",
      Skin >= 50 ~ "Skin",
      Nasal >= 50 ~ "Nasal",
      Oral >= 50 ~ "Oral",
      TRUE ~ "Mixed"
    )
  )

molecular_info %>%
  plyr::dlply(.variables = .(molecular_class)) %>%
  purrr::map(function(x) {
    table(x$direction) * 100 / nrow(x)
  })


microbiome_info

#####output result
library(openxlsx)
wb <- createWorkbook()
modifyBaseFont(wb, fontSize = 12, fontName = "Time New Roma")
addWorksheet(wb, sheetName = "Node data", gridLines = TRUE)
addWorksheet(wb, sheetName = "Edge data", gridLines = TRUE)
addWorksheet(wb, sheetName = "Microbiome information", gridLines = TRUE)
addWorksheet(wb, sheetName = "Molecular information", gridLines = TRUE)
freezePane(wb,
           sheet = 1,
           firstRow = TRUE,
           firstCol = TRUE)
freezePane(wb,
           sheet = 2,
           firstRow = TRUE,
           firstCol = TRUE)
freezePane(wb,
           sheet = 3,
           firstRow = TRUE,
           firstCol = TRUE)
freezePane(wb,
           sheet = 4,
           firstRow = TRUE,
           firstCol = TRUE)
writeDataTable(
  wb,
  sheet = 1,
  x = as.data.frame(node_data),
  colNames = TRUE,
  rowNames = FALSE
)
writeDataTable(
  wb,
  sheet = 2,
  x = as.data.frame(edge_data),
  colNames = TRUE,
  rowNames = FALSE
)
writeDataTable(
  wb,
  sheet = 3,
  x = microbiome_info,
  colNames = TRUE,
  rowNames = FALSE
)
writeDataTable(
  wb,
  sheet = 4,
  x = molecular_info,
  colNames = TRUE,
  rowNames = FALSE
)
saveWorkbook(wb = wb,
             file = "inter_omics_cor_network.xlsx",
             overwrite = TRUE)


####subnetwork
#####for each microbiome, give it a class according to connection to which omics
node = igraph::vertex_attr(temp_data)$node
module =
  data.frame(node = node) %>%
  dplyr::left_join(microbiome_info[, c("microbiome_id", "total_number", "direction")],
                   by = c('node' = "microbiome_id"))

module$direction[is.na(module$direction)] = "Other"
module$direction[module$direction == "Mixed"] = "Other"
module$direction[which(module$total_number < 6)] = "Other"

module = module$direction

temp_data =
  temp_data %>%
  tidygraph::activate(what = "nodes") %>%
  dplyr::mutate(module = module)

# subnetworks <-
#   igraph::cluster_fast_greedy(graph = temp_data,
#                                    weights = abs(edge_attr(temp_data,
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
# temp_data =
#   temp_data %>%
#   tidygraph::activate(what = "nodes") %>%
#   dplyr::mutate(module = module)

angle <- 360 * (c(1:nrow(node_data)) - 0.5) / nrow(node_data)
hjust <- ifelse(angle > 180, 1, 0)
angle <- ifelse(angle > 180, 90 - angle + 180, 90 - angle)

node_data =
  igraph::vertex_attr(temp_data) %>%
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

plot <-
  ggraph(temp_data,
         layout = 'fr',
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
  geom_edge_link(aes(color = pos_neg,
                     width = -log(p_adjust, 10)
                     # linetype = significance),
                     alpha = 1,
                     show.legend = TRUE) +
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
                   shadowtext::geom_shadowtext(
                     aes(
                       x = x,
                       y = y,
                       label = ifelse(node %in% label_node, true_name, NA),
                       color = class
                     ),
                     size = 5,
                     check_overlap = TRUE,
                     bg.colour = "white",
                     show.legend = FALSE
                   ) +
                   guides(edge_color = guide_edge_colorbar(title = "Correlation"),
                          color = guide_legend(title = "Class")) +
                   scale_edge_width_continuous(range = c(0.1, 1)) +
                   scale_size_continuous(range = c(2, 8)) +
                   scale_fill_manual(values = omics_color) +
                   scale_color_manual(values = omics_color) +
                   scale_edge_color_manual(values = c(
                     "positive" = viridis::inferno(n = 10)[2],
                     "negative" = viridis::inferno(n = 10)[5]
                   )) +
                   ggraph::theme_graph() +
                   theme(
                     plot.background = element_rect(fill = "transparent", color = NA),
                     panel.background = element_rect(fill = "transparent", color = NA),
                     legend.position = "right",
                     legend.background = element_rect(fill = "transparent", color = NA)
                   )
                 
                 plot
                 
                 # extrafont::loadfonts()
                 ggsave(plot,
                        filename = "microbiome_molecular_cor_plot.pdf",
                        width = 17,
                        height = 14)
                 