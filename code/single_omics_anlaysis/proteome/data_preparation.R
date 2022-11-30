###
no_function()

masstools::setwd_project()
library(tidyverse)
rm(list = ls())
####load raw data
load("data/from_xin/Revision_MultiOmes_0509.RData")
setwd("data_analysis/proteome/data_preparation")
expression_data = 
  swathprot.PCR.df

colnames(expression_data)

sample_info = 
  expression_data %>% 
  dplyr::select(SampleID, SubjectID:CL4)

expression_data =
  expression_data %>% 
  dplyr::select(-c(SampleID, SubjectID:CL4)) %>% 
  t() %>%
  as.data.frame()

colnames(expression_data) = sample_info$SampleID

variable_info =
  data.frame(variable_id = rownames(expression_data))

dim(expression_data)
dim(sample_info)
dim(variable_info)

colnames(expression_data) == sample_info$SampleID
rownames(expression_data)  == variable_info$variable_id

sample_info = 
  sample_info %>% 
  dplyr::rename(sample_id = SampleID,
                subject_id = SubjectID)

###add information to variable_info
variable_info$variable_id

library(clusterProfiler)
library(org.Hs.eg.db)
library(plyr)
library(UniprotR)

UNIPROT =
  clusterProfiler::bitr(
    geneID = variable_info$variable_id,
    fromType = "SYMBOL",
    toType = "UNIPROT",
    OrgDb = org.Hs.eg.db,
    drop = FALSE
  )

write.csv(UNIPROT, "UNIPROT.csv", row.names = FALSE)

UNIPROT =
  UNIPROT %>%
  plyr::dlply(.variables = .(SYMBOL)) %>%
  purrr::map(function(x) {
    if (nrow(x) == 1) {
      return(x)
    }
    
    temp2 = UniprotR::GetMiscellaneous(ProteinAccList = x$UNIPROT)
    
    x$status = temp2$Status
    
    x = 
      x %>% 
      dplyr::filter(status == "reviewed")
    
    if (nrow(x) == 1) {
      return(x)
    }
    x
  })


UNIPROT %>% lapply(nrow) %>% unlist() %>% plot

UNIPROT %>% lapply(nrow) %>% unlist() %>% `>`(1) %>% which()

UNIPROT[[48]] = 
  UNIPROT[[48]][2,]

UNIPROT[[262]] = 
  UNIPROT[[262]][2,]

UNIPROT = 
UNIPROT %>% 
  purrr::map(function(x){
    x[,1:2]
  }) %>% 
  do.call(rbind, .) %>% 
  as.data.frame()

which(is.na(UNIPROT$UNIPROT))

write.csv(UNIPROT, "UNIPROT.csv", row.names = FALSE)

UNIPROT = readr::read_csv("UNIPROT_manual.csv")

variable_info = 
variable_info %>% 
  dplyr::left_join(UNIPROT, by = c("variable_id" = "SYMBOL"))

variable_info$variable_id == rownames(expression_data)

colnames(variable_info)

###add entrizid
ENTREZID =
  clusterProfiler::bitr(
    geneID = variable_info$UNIPROT,
    fromType = "UNIPROT",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db,
    drop = FALSE
  ) %>% 
  dplyr::distinct(UNIPROT,.keep_all = TRUE)

variable_info = 
  variable_info %>% 
  dplyr::left_join(ENTREZID, by = c("UNIPROT"))

which(sample_info$sample_id == "69-066-05")
sample_info[538,]$subject_id = "69-066"

colnames(sample_info)
colnames(sc)

sample_info <-
  sample_info %>% 
  dplyr::left_join(sc, by = c("subject_id" = "SubjectID"))





# cor_data <-
#   purrr::map(1:(nrow(expression_data) - 1), function(idx1) {
#     cat(idx1, " ")
#     purrr::map((idx1 + 1):nrow(expression_data), function(idx2) {
#       x = as.numeric(expression_data[idx1, ])
#       y = as.numeric(expression_data[idx2, ])
#       temp =
#         cor.test(x, y, method = "spearman")
#       data.frame(
#         name1 = variable_info$variable_id[idx1],
#         name2 = variable_info$variable_id[idx2],
#         cor = unname(temp$estimate),
#         p = unname(temp$p.value)
#       )
#     }) %>%
#       do.call(rbind, .) %>%
#       as.data.frame()
#   }) %>%
#   do.call(rbind, .) %>%
#   as.data.frame()
# save(cor_data, file = "cor_data")

load("cor_data")

cor_data$p_adjust = p.adjust(cor_data$p, method = "BH")

plot <- 
  cor_data %>% 
  ggplot(aes(cor)) +
  geom_histogram(color = "black", binwidth = 0.05) +
  theme_bw() +
  labs(x = "Spearman correlation",
       y = "Count")
plot
# ggsave(plot, filename = "correlation distributation.pdf", width = 9, height = 7)


save(sample_info, file = "sample_info")
save(expression_data, file = "expression_data")
save(variable_info, file = "variable_info")

