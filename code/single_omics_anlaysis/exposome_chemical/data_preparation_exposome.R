###not source
not_exist()

masstools::setwd_project()
setwd("data_20200511/exposome/")
library(tidyverse)
rm(list=ls())

metabolite_table <-
  readxl::read_xlsx("SmallMoleculeProfiling_20161129Summary_neg2.xlsx", sheet = 1)

peak_table <-
  readxl::read_xlsx("SmallMoleculeProfiling_20161129Summary_neg2.xlsx", sheet = 2)

sample_info <-
  readxl::read_xlsx("SmallMoleculeProfiling_20161129Summary_neg2.xlsx", sheet = 3)

dim(metabolite_table)

colnames(metabolite_table)[1] <- "peak_ID"

metabolite_table2 <- 
  metabolite_table %>% 
  dplyr::filter(!is.na(peak_ID)) %>% 
  dplyr::filter(stringr::str_detect(peak_ID, "[PN]{1}M[0-9]{1,5}"))


variable_info <- 
  metabolite_table2 %>% 
  select(peak_ID:Fold_Enrich)


expression_data <- 
  metabolite_table2 %>% 
  select(-c(peak_ID:Fold_Enrich))

colnames(expression_data)

expression_data <- 
  expression_data %>% 
  dplyr::select(contains("B0"))

class(expression_data)
expression_data <- as.data.frame(expression_data)

# idx1 <- match(variable_info$peak_ID, peak_table$MetID)
# idx2 <- match(colnames(expression_data), colnames(peak_table))
# 
# expression_data <- peak_table[idx1, idx2]

colnames(expression_data)

rownames(expression_data) <- variable_info$peak_ID

sum(is.na(expression_data))

save(variable_info, file = "variable_info")
save(expression_data, file = "expression_data")


sample_info

##add the latitude and longitude for different locations
sample_info
colnames(sample_info)[1] <- "sample_id"
colnames(sample_info)[2] <- "start_date"
colnames(sample_info)[3] <- "end_date"
colnames(sample_info)[4] <- "comments"


sample_info <- 
  sample_info %>% 
  dplyr::mutate(location = comments)

latitude <- c(37.408899,
              37.408899,
              37.408899,
              37.408899,
              37.408899,
              38.106940,
              # 38.106940,
              37.408899,
              37.408899,
              37.408899,
              46.601860,
              29.764424,
              37.408899,
              37.408899,
              37.408899,
              42.362238,
              37.408899,
              41.876097,
              38.552765,
              38.984273,
              42.362238)

longitude <- c(-122.149814, 
               -122.149814, 
               -122.149814, 
               -122.149814, 
               -122.149814,
               -122.564569,
               # -122.564569,
               -122.149814,
               -122.149814,
               -122.149814,
               -112.038093,
               -95.368194,
               -122.149814, 
               -122.149814, 
               -122.149814, 
               -71.056961,
               -122.149814, 
               -87.627601,
               -121.451996,
               -77.094888,
               -71.056961
               )

  
sample_info <-
  data.frame(sample_info, 
             longitude, 
             latitude, 
             stringsAsFactors = FALSE)

sample_info$location[grep("Campus", sample_info$location)] <- "Campus"
sample_info$location[grep("Porter office", sample_info$location)] <- "Campus"


save(sample_info, file = "sample_info")


# ###read variable information from peng
# variable_info <- readxl::read_excel("exposome_variable_info_peng_revised_200512.xlsx")
# 
# 
# save(variable_info, file = "variable_info")


library(openxlsx)
wb = createWorkbook()
modifyBaseFont(wb, fontSize = 12, fontName = "Arial Narrow")
addWorksheet(wb, sheetName = "Sample information", gridLines = TRUE)
addWorksheet(wb, sheetName = "Variable information", gridLines = TRUE)
addWorksheet(wb, sheetName = "Expression data", gridLines = TRUE)
freezePane(wb, sheet = 1, firstRow = TRUE, firstCol = TRUE) 
freezePane(wb, sheet = 2, firstRow = TRUE, firstCol = TRUE) 
freezePane(wb, sheet = 3, firstRow = TRUE, firstCol = FALSE) 
writeDataTable(wb, sheet = 1, x = sample_info,
               colNames = TRUE, rowNames = FALSE)
writeDataTable(wb, sheet = 2, x = variable_info,
               colNames = TRUE, rowNames = FALSE)
writeDataTable(wb, sheet = 3, x = expression_data,
               colNames = TRUE, rowNames = FALSE)
saveWorkbook(wb, "exposome_chemical_data.xlsx", overwrite = TRUE)


