



library(dplyr)

pvalue <- read.csv("~/Desktop/FigureS7/p_values_comparision_1000_vs-10000.csv", header = T)
table.DMI <- read.csv("~/Desktop/FigureS7/New_DMI_TABLE.csv", header = T, row.names = 1)
table.FS <- read.csv("~/Desktop/FigureS7/Extended Data Table FS.csv", header = T)

table.New <- full_join(table.DMI, table.FS, by = c("genus" = "genus", "bodysite" = "class"))

table(table.New$bodysite, table.New$trimmed)

table.withP <- inner_join(pvalue, table.New, by = c("genus" = "genus", "data_class" = "bodysite"))

colnames(table.S7)
max(table.New$fc1_p_adjust_10000)

top_20_values <- table.withP %>%
  arrange(desc(fc1_p_adjust_10000)) %>%  
  slice(1:20)  
top_20_values

filter(table.New, data_class == "Nasal")
filter(table.S7, bodysite == "Nasal")

filter(table.withP, trimmed == "no")$fc1_p_adjust_10000
table(table.New$confidance)
table(table.New$data_class)

filter(table.withP , genus== "Campylobacter")
filter(table.withP , genus== "Leptotrichia")
filter(table.withP , genus== "Cutibacterium")
filter(table.withP , genus== "Cutibacterium")

write.csv(file = "~/Desktop/FigureS7/tableS7_2024.csv", table.withP)

table.withP.sig <- filter(table.withP, trimmed == "no")

write.csv(file = "~/Desktop/FigureS7/tableS7_2024.sig.csv", table.withP.sig)

filter(pvalue , genus== "Campylobacter")
filter(pvalue , genus== "Leptotrichia")
filter(pvalue , genus== "Cutibacterium")

filter(table.S7 , genus== "Campylobacter")
filter(table.S7 , genus== "Leptotrichia")
filter(table.S7 , genus== "Cutibacterium")


