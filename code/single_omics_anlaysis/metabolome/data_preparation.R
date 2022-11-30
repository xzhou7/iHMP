###
no_function()

masstools::setwd_project()
library(tidyverse)
rm(list = ls())
####load raw data
load("data/from_xin/Revision_MultiOmes_0509.RData")

variable_info = metb.curated %>%
  dplyr::select(Compounds_ID, dplyr::everything()) %>%
  dplyr::rename(variable_id = Compounds_ID)

expression_data = metbcr.df

sample_info =
  expression_data %>%
  dplyr::select(SampleID, SubjectID:CL4) %>%
  dplyr::rename(sample_id = SampleID,
                subject_id = SubjectID) %>%
  dplyr::left_join(sc, by = c("subject_id" = "SubjectID"))

expression_data =
  expression_data %>%
  dplyr::select(-c(SampleID, SubjectID:CL4)) %>%
  t() %>%
  as.data.frame()

dim(variable_info)
dim(expression_data)
dim(sample_info)

colnames(expression_data) = sample_info$sample_id

setdiff(rownames(expression_data),
        variable_info$variable_id)

setdiff(variable_info$variable_id,
        rownames(expression_data))

length(variable_info$variable_id)
length(unique(variable_info$variable_id))

which(duplicated(variable_info$variable_id))
variable_info %>%
  dplyr::filter(variable_id %in% variable_info$variable_id[c(383, 435)])

variable_info =
  variable_info %>%
  dplyr::distinct(variable_id, .keep_all = TRUE)

dim(variable_info)

variable_info$variable_id == rownames(expression_data)
sort(variable_info$variable_id) == sort(rownames(expression_data))
expression_data =
  expression_data[variable_info$variable_id, ]

rownames(expression_data) == variable_info$variable_id
colnames(expression_data) == sample_info$sample_id

masstools::setwd_project()
setwd("data_analysis/metabolome/data_preparation/")

variable_info =
  variable_info %>%
  dplyr::filter(!stringr::str_detect(Metabolite, "C[0-9]{1,2}H[0-9]{1,2}")) %>%
  dplyr::filter(!stringr::str_detect(Metabolite, "Unknown")) %>%
  dplyr::mutate(Metabolite = stringr::str_replace_all(Metabolite, '\"', "")) %>%
  dplyr::mutate(true_name =
                  stringr::str_replace_all(Metabolite, "\\([0-9]{1}\\)", ""))

mean_int =
  expression_data[variable_info$variable_id, ] %>%
  apply(1, mean)

variable_info$mean_int = mean_int

variable_info =
  variable_info %>%
  plyr::dlply(.variables = .(true_name)) %>%
  purrr::map(function(x) {
    x %>%
      dplyr::filter(MSMS == "YES") %>%
      dplyr::filter(mean_int == max(mean_int))
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()

variable_info$Metabolite = variable_info$true_name

expression_data = expression_data[variable_info$variable_id, ]

# write.csv(variable_info, "variable_info.csv", row.names = FALSE)

variable_info = readr::read_csv("variable_info_manual.csv")

variable_info$variable_id == rownames(expression_data)

library(masstools)

idx = which(!is.na(variable_info$HMDB) & is.na(variable_info$KEGG))

hmdb_kegg =
  purrr::map(
    variable_info$HMDB[idx],
    .f = function(x) {
      x = stringr::str_split(x, pattern = "\\|")[[1]][1]
      result = masstools::trans_ID(
        query = x,
        from = "Human Metabolome Database",
        to = "KEGG",
        top = 1
      )
      result
    }
  ) %>%
  do.call(rbind, .) %>%
  as.data.frame()

nrow(hmdb_kegg)
length(idx)

variable_info$KEGG[idx] = hmdb_kegg$KEGG[idx]


variable_info$HMDB =
  variable_info$HMDB %>%
  purrr::map(function(x) {
    if (is.na(x)) {
      return(x)
    }
    
    x = stringr::str_split(x, "\\|")[[1]]
    
    x =
      purrr::map(x, function(y) {
        if (nchar(y) == 9) {
          y = stringr::str_replace(y, "HMDB", "HMDB00")
        }
        y
      }) %>%
      unlist() %>%
      paste(., collapse = "|")
    x
  }) %>%
  unlist()



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
