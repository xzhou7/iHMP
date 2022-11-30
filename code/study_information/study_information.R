
###
no_function()

masstools::setwd_project()
library(tidyverse)
library(phyloseq)
rm(list = ls())

####load raw data
source(here::here("code/tools.R"))

#####clinical information
load("data/from_xin/ls")

ls$CL1

###proteome
load("data_analysis/proteome/data_preparation/sample_info")
proteome_sample_info <-
  sample_info

load("data_analysis/proteome/data_preparation/variable_info")
proteome_variable_info <-
  variable_info

###metabolome
load("data_analysis/metabolome/data_preparation/sample_info")
metabolome_sample_info <-
  sample_info

load("data_analysis/metabolome/data_preparation/variable_info")
metabolome_variable_info <-
  variable_info

###lipidome
load("data_analysis/lipidome/data_preparation/sample_info")
lipidome_sample_info <-
  sample_info

load("data_analysis/lipidome/data_preparation/variable_info")
lipidome_variable_info <-
  variable_info

###cytokine
load("data_analysis/cytokine/data_preparation/sample_info")
cytokine_sample_info <-
  sample_info

load("data_analysis/cytokine/data_preparation/variable_info")
cytokine_variable_info <-
  variable_info

###stool_microbiome
load("data_analysis/stool_microbiome/data_preparation/sample_info")
stool_microbiome_sample_info <-
  sample_info

load("data_analysis/stool_microbiome/data_preparation/variable_info")
stool_microbiome_variable_info <-
  variable_info

###skin_microbiome
load("data_analysis/skin_microbiome/data_preparation/sample_info")
skin_microbiome_sample_info <-
  sample_info

load("data_analysis/skin_microbiome/data_preparation/variable_info")
skin_microbiome_variable_info <-
  variable_info

###oral_microbiome
load("data_analysis/oral_microbiome/data_preparation/sample_info")
oral_microbiome_sample_info <-
  sample_info

load("data_analysis/oral_microbiome/data_preparation/variable_info")
oral_microbiome_variable_info <-
  variable_info

###nasal_microbiome
load("data_analysis/nasal_microbiome/data_preparation/sample_info")
nasal_microbiome_sample_info <-
  sample_info

load("data_analysis/nasal_microbiome/data_preparation/variable_info")
nasal_microbiome_variable_info <-
  variable_info

###clinical lab test
load("data/from_xin/Revision_MultiOmes_0509.RData")
clinic.df$SampleID

lipidome_sample_info2 <- 
  lipidome_sample_info %>% 
  dplyr::filter(!stringr::str_detect(sample_id, "QC"))

##plasma sample number
c(proteome_sample_info$sample_id,
  metabolome_sample_info$sample_id,
  lipidome_sample_info2$sample_id,
  cytokine_sample_info$SampleID,
  clinic.df$SampleID) %>% 
  unique() %>% 
  length()

length(stool_microbiome_sample_info$sample_id)
length(skin_microbiome_sample_info$sample_id)
length(oral_microbiome_sample_info$sample_id)
length(nasal_microbiome_sample_info$sample_id)


##visit number
c(paste(proteome_sample_info$subject_id,as.Date(proteome_sample_info$CollectionDate, format = "%m/%d/%y")),
  paste(metabolome_sample_info$subject_id, as.Date(metabolome_sample_info$CollectionDate, format = "%m/%d/%y")),
  paste(lipidome_sample_info2$subject_id, lipidome_sample_info2$date),
  paste(cytokine_sample_info$subject_id, as.Date(cytokine_sample_info$CollectionDate, format = "%m/%d/%y")),
  paste(clinic.df$subject_id, as.Date(clinic.df$CollectionDate, format = "%m/%d/%y")),
  paste(stool_microbiome_sample_info$subject_id, stool_microbiome_sample_info$Date),
  paste(skin_microbiome_sample_info$subject_id, skin_microbiome_sample_info$Date),
  paste(oral_microbiome_sample_info$subject_id, oral_microbiome_sample_info$Date),
  paste(nasal_microbiome_sample_info$subject_id, nasal_microbiome_sample_info$Date)
) %>% 
  unique() %>% 
  length()

nrow(proteome_variable_info) * nrow(proteome_sample_info) +
724 * nrow(metabolome_sample_info) +
nrow(lipidome_variable_info) * nrow(lipidome_sample_info2) +
62 * nrow(cytokine_sample_info) +
nrow(clinic.df) * (ncol(clinic.df) - 7) +
nrow(stool_microbiome_variable_info) * nrow(stool_microbiome_sample_info) +
nrow(skin_microbiome_variable_info) * nrow(skin_microbiome_sample_info) +
nrow(oral_microbiome_variable_info) * nrow(oral_microbiome_sample_info) +
nrow(nasal_microbiome_variable_info) * nrow(nasal_microbiome_sample_info)



dir.create("data_analysis/study_information")
setwd("data_analysis/study_information")

library(ComplexUpset)




temp_data =
  rbind(
    data.frame(sample_id = proteome_sample_info[, c("sample_id")], class = "Proteome"),
    data.frame(sample_id = metabolome_sample_info[, c("sample_id")], class = "Metabolite"),
    data.frame(sample_id = lipidome_sample_info[, c("sample_id")], class = "Lipidome"),
    data.frame(
      sample_id = cytokine_sample_info[, c("SampleID"), drop = FALSE] %>% dplyr::rename(sample_id = "SampleID"),
      class = "Cytokine"
    ),
    data.frame(sample_id = stool_microbiome_sample_info[, c("sample_id")], class = "Stool microbiome"),
    data.frame(sample_id = skin_microbiome_sample_info[, c("sample_id")], class = "Skin microbiome"),
    data.frame(sample_id = nasal_microbiome_sample_info[, c("sample_id")], class = "Nasal microbiome"),
    data.frame(sample_id = oral_microbiome_sample_info[, c("sample_id")], class = "Oral microbiome")
  ) %>% 
  dplyr::mutate(value = TRUE)

temp_data = 
temp_data %>%
  tidyr::pivot_wider(
    names_from = "class",
    values_from = "value"
  ) %>% 
  tibble::column_to_rownames(var = "sample_id")

temp_data[is.na(temp_data)] = FALSE

total_plot = 
upset(
  data = temp_data,
  intersect = colnames(temp_data),
  name = 'genre',
  width_ratio = 0.1,
  # min_size=20,
  set_sizes=(
    upset_set_size()
    + geom_text(aes(label=..count..), hjust=1.1, stat='count')
    # you can also add annotations on top of bars:
    + annotate(geom='text', label='@', x='Drama', y=850, color='white', size=3)
    + expand_limits(y=1100)
    + theme(axis.text.x=element_text(angle=90))
  )
)

total_plot

# ggsave(total_plot, filename = "total_plot.pdf", width = 21, height = 7)

plot = 
  upset(
    data = temp_data,
    intersect = colnames(temp_data),
    name = 'genre',
    width_ratio = 0.1,
    min_size=20,
    set_sizes=(
      upset_set_size()
      + geom_text(aes(label=..count..), hjust=1.1, stat='count')
      # you can also add annotations on top of bars:
      + annotate(geom='text', label='@', x='Drama', y=850, color='white', size=3)
      + expand_limits(y=1100)
      + theme(axis.text.x=element_text(angle=90))
    )
  )

plot

# ggsave(plot, filename = "plot_upset.pdf", width = 8, height = 6)


list <-
  c(
    "69-001",
    "69-003",
    "69-010",
    "69-012",
    "69-013",
    "69-016",
    "69-021",
    "69-022",
    "69-023",
    "69-026",
    "69-027",
    "70-1014",
    "69-031" ,
    "69-032" ,
    "69-033" ,
    "69-034",
    "69-035" ,
    "69-036" ,
    "69-037" ,
    "69-039" ,
    "69-040" ,
    "69-041" ,
    "69-043" ,
    "69-045" ,
    "69-046" ,
    "69-047",
    "69-048" ,
    "69-052",
    "69-053",
    "69-055" ,
    "69-056"  ,
    "69-057"  ,
    "69-058"  ,
    "69-060",
    "69-063" ,
    "69-064",
    "69-065",
    "69-066" ,
    "69-068",
    "69-069",
    "69-070",
    "69-071" ,
    "69-073",
    "69-074"  ,
    "69-076" ,
    "69-077" ,
    "69-078" ,
    "69-079" ,
    "69-080",
    "69-081",
    "69-085",
    "69-086" ,
    "69-087" ,
    "69-090" ,
    "69-091" ,
    "69-095",
    "69-097" ,
    "69-099" ,
    "69-100" ,
    "69-103"  ,
    "69-104" ,
    "69-107" ,
    "69-109" ,
    "69-110" ,
    "69-111",
    "69-112",
    "69-113",
    "69-114",
    "69-116",
    "69-118" ,
    "69-120",
    "69-121" ,
    "69-122",
    "69-123" ,
    "69-124" ,
    "70-1001" ,
    "70-1002",
    "70-1003" ,
    "70-1004",
    "70-1005",
    "70-1006",
    "70-1008",
    "70-1010",
    "70-1011",
    "70-1012",
    "70-1015"
  )


temp_data = 
metabolome_sample_info %>% 
  dplyr::select(subject_id,IRIS) %>% 
  dplyr::filter(subject_id %in% list) %>% 
  dplyr::distinct(subject_id, .keep_all = TRUE) %>% 
  dplyr::count(IRIS) 

temp_data = 
temp_data %>% 
  dplyr::filter(!is.na(IRIS))

library(ggiraphExtra)

plot =
  ggiraphExtra::ggPie(
    temp_data,
    aes(pies = IRIS, count = n)
    # colour = c(iris_color, Unknown = "grey"),
    # fill = c(iris_color, Unknown = "grey"),
  )
plot

ggsave(plot, filename = "IRIS_pie.pdf", width = 7, height = 7)

library(plyr)

ls %>%
  plyr::dlply(.variables = .(SubjectID)) %>% 
  purrr::map(function(x){
    sum(stringr::str_detect("Infection",x$CL4))*100/nrow(x)
  }) %>% 
  unlist() %>% 
  plot()

ls_new <-
  ls %>% 
    dplyr::filter(SubjectID %in% list)

temp_data = 
  data.frame(n = c(sum(stringr::str_detect("Infection",ls_new$CL4), na.rm = TRUE),
                        nrow(ls_new) - sum(stringr::str_detect("Infection",ls_new$CL4), na.rm = TRUE)),
             status = c("Infection", "Healthy")) 

plot <-
  ggiraphExtra::ggPie(
    temp_data,
    aes(pies = status, count = n)
  )

plot

ggsave(plot, filename = "infection_pie.pdf", width = 7, height = 7)




  
