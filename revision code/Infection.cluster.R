library(readxl)
library(tidyverse)
library(stringr)

setwd("~/Library/CloudStorage/Box-Box/human_microbiome_project/")
source("./code/tools.R")

IR_Cluster1 <- read_xlsx("./data_analysis/combine_omics/clustering_IR/cluster_1/cluster1_with_p_value.xlsx")
IR_Cluster2 <- read_xlsx("./data_analysis/combine_omics/clustering_IR/cluster_2/cluster2_with_p_value.xlsx")
IR_Cluster3 <- read_xlsx("./data_analysis/combine_omics/clustering_IR/cluster_3/cluster3_with_p_value.xlsx")

IS_Cluster1 <- read_xlsx("./data_analysis/combine_omics/clustering_IS/cluster_1/cluster1_with_p_value.xlsx")
IS_Cluster2 <- read_xlsx("./data_analysis/combine_omics/clustering_IS/cluster_2/cluster2_with_p_value.xlsx")
IS_Cluster3 <- read_xlsx("./data_analysis/combine_omics/clustering_IS/cluster_3/cluster3_with_p_value.xlsx")

IR_Cluster1$bodysite <- str_extract(IR_Cluster1$variable_id,"skin|oral|nasal|stool")
IR_Cluster1$bodysite[is.na(IR_Cluster1$bodysite)] <- "other"
IR_Cluster1$Cluster <- "1"

IR_Cluster2$bodysite <- str_extract(IR_Cluster2$variable_id, "skin|oral|nasal|stool")
IR_Cluster2$bodysite[is.na(IR_Cluster2$bodysite)] <- "other"
IR_Cluster2$Cluster <- "2"

IR_Cluster3$bodysite <- str_extract(IR_Cluster3$variable_id, "skin|oral|nasal|stool")
IR_Cluster3$bodysite[is.na(IR_Cluster3$bodysite)] <- "other"
IR_Cluster3$Cluster <- "3"

IS_Cluster1$bodysite <- str_extract(IS_Cluster1$variable_id, "skin|oral|nasal|stool")
IS_Cluster1$bodysite[is.na(IS_Cluster1$bodysite)] <- "other"
IS_Cluster1$Cluster <- "1"

IS_Cluster2$bodysite <- str_extract(IS_Cluster2$variable_id, "skin|oral|nasal|stool")
IS_Cluster2$bodysite[is.na(IS_Cluster2$bodysite)] <- "other"
IS_Cluster2$Cluster <- "2"

IS_Cluster3$bodysite <- str_extract(IS_Cluster3$variable_id, "skin|oral|nasal|stool")
IS_Cluster3$bodysite[is.na(IS_Cluster3$bodysite)] <- "other"
IS_Cluster3$Cluster <- "3"

IR.table <- rbind(filter(IR_Cluster1, bodysite != "other"), filter(IR_Cluster2, bodysite != "other"), filter(IR_Cluster3, bodysite != "other")) %>% select(variable_id:Species,Cluster, p_value, bodysite)
IR.table$IRIS <- "IR"
IR.table <- mutate(IR.table, classification = case_when(
  Cluster == 1 ~ "B",
  Cluster == 2 ~ "C",
  Cluster == 3 ~ "A"))

IS.table <- rbind(filter(IS_Cluster1, bodysite != "other"), filter(IS_Cluster2, bodysite != "other"), filter(IS_Cluster3, bodysite != "other")) %>% select(variable_id:Species,Cluster, p_value, bodysite)
IS.table$IRIS <- "IS"

IS.table <- mutate(IS.table, classification = case_when(
  Cluster == 1 ~ "A",
  Cluster == 2 ~ "B",
  Cluster == 3 ~ "C"))

figure3E_taxa <- rbind(IR.table, IS.table)
table(figure3E_taxa$IRIS, figure3E_taxa$bodysite)

table(IR.table$classification, IR.table$bodysite)
table(IS.table$classification, IS.table$bodysite)

#output an supplimentary table 
Supplimentary_table <- rbind(IS.table, IR.table)
write.csv(file = "~/Desktop/infection.cluster.csv", Supplimentary_table)

hist(Supplimentary_table$p_value)

# Create separate data frames for IR and IS counts by bodysite and classification
IR_counts <- IR.table %>% count(bodysite, classification)
IS_counts <- IS.table %>% count(bodysite, classification)

# Add a new column to indicate the source
IR_counts$source <- "IR"
IS_counts$source <- "IS"

# Combine the IR and IS counts into a single data frame
combined_counts <- bind_rows(IR_counts, IS_counts)
combined_counts$bodysite <- factor(combined_counts$bodysite, levels = c("stool", "skin", "oral","nasal"))
combined_counts$source <- factor(combined_counts$source, levels = c("IS","IR"))
# Create the bar graph
p.compare <- combined_counts %>% filter(classification != "C") %>% group_by(source, bodysite) %>%
  mutate(percentage = n / sum(n) * 100) %>% 
  ggplot(aes(x = source, y = percentage, fill = classification)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(~ bodysite) +
  labs(title = "IR vs IS cluster association during infection",
       x = "Insulin Status",
       y = "percentage") +
  theme_minimal()
p.compare


#ggsave(filename = "~/Desktop/Suppli.to3E.pdf", p.compare, width = 4, height = 4, dpi = 300)

cat("Stool\n")
chisq_stool <- chisq.test(rbind(c(29, 17), c(21, 21)))
print(chisq_stool)

cat("Skin\n")
chisq_skin <- chisq.test(rbind(c(37, 12), c(15, 11)))
print(chisq_skin)

cat("Oral\n")
chisq_oral <- chisq.test(rbind(c(20, 7), c(9, 27)))
print(chisq_oral)

cat("/Nasal")
chisq_nasal <- chisq.test(rbind(c(8,23), c(28,12)))
print(chisq_nasal)

############################################
#compare IR IS 
############################################
percentage_table2 <- combined_table %>%
  group_by(IRIS, classification, bodysite) %>%
  summarize(count = n()) %>%
  mutate(percentage = count / sum(count) * 100)

percentage_table2$IRIS <- factor(percentage_table2$IRIS, levels = c("IS", "IR"))

# Plot the bar graph
ggplot(percentage_table2, aes(x = IRIS, y = percentage, fill = bodysite)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(~ classification) +
  labs(x = "IRIS", y = "Percentage") +
  theme_minimal() +
  scale_fill_manual(values = c("stool" = "#DF8F44FF", "skin" = "#00A1D5FF", "oral" = "#B24745FF", "nasal" = "#79AF97FF")) +
  theme(legend.title = element_blank())

# 
# data.S8 <- read.csv("~/Desktop/TableS8/tableS8.csv", header = T)
# # Assuming your data is in a data frame called 'data'
# wide_table <- data.S8 %>%
#   spread(key = "Correlation.Body.Sites.and.Directions", value = "Counts") %>%
#   mutate(RowSum = rowSums(select(., Nasal_Negative:Stool_Positive), na.rm = TRUE))
# 
# # Calculate column sums
# col_sums <- colSums(select(wide_table,Nasal_Negative:Stool_Positive), na.rm = TRUE)
# 
# # Add the column sums as a new row to the wide table
# wide_table_with_sums <- wide_table %>%
#   add_row(Cytokine = "ColSum", 
#           Nasal_Negative = col_sums["Nasal_Negative"], Nasal_Positive = col_sums["Nasal_Positive"],
#           Oral_Negative = col_sums["Oral_Negative"], Oral_Positive = col_sums["Oral_Positive"],
#           Skin_Negative = col_sums["Skin_Negative"], Skin_Positive = col_sums["Skin_Positive"],
#           Stool_Negative = col_sums["Stool_Negative"], Stool_Positive = col_sums["Stool_Positive"],
#           RowSum = sum(col_sums))
# 
# write.csv(file = "~/Desktop/TableS8/new.table.S8.csv",wide_table_with_sums)

data.s11 <- rbind(read_xlsx("~/Desktop/NewtableS11/stool_cytokine_result_estpval-fdr.xlsx"),
                  read_xlsx("~/Desktop/NewtableS11/skin_cytokine_result_estpval-fdr.xlsx"),
                  read_xlsx("~/Desktop/NewtableS11/oral_cytokine_result_estpval-fdr.xlsx"),
                  read_xlsx("~/Desktop/NewtableS11/nasal_cytokine_result_estpval-fdr.xlsx"))

my_table.s11 <- data.s11 %>%
  mutate(significance = ifelse(low * high > 0, "yes", "no"))

table(my_table.s11$site, my_table.s11$significance)
#perform test between body sites
contingency_table <- table(my_table.s11$site, my_table.s11$significance)
site_pairs <- combn(unique(my_table.s11$site), 2, simplify = FALSE)
pairwise_chi_squared <- lapply(site_pairs, function(pair) {
  sub_table <- contingency_table[pair,]
  chisq <- chisq.test(sub_table)
  list(pair = paste(pair, collapse = " vs. "), p_value = chisq$p.value)
})
results_df <- bind_rows(pairwise_chi_squared)
results_df

write.csv(file = "~/Desktop/NewtableS11/new_s11.csv", my_table.s11)


table(IR_Cluster1$p_value < 0.05)
table(IR_Cluster2$p_value < 0.05)
table(IR_Cluster3$p_value < 0.05)

table(IS_Cluster1$p_value < 0.05)
table(IS_Cluster2$p_value < 0.05)
table(IS_Cluster3$p_value < 0.05)

#adjust for p values
p_values <- c(0.01545, 6.77E-08, 6.75E-08, 6.78E-08, 6.76E-08, 6.78E-08, 6.77E-08,
              6.74E-08, 6.66E-08, 6.77E-08, 6.75E-08, 6.67E-08, 7.88E-08, 6.76E-08,
              0.2732, 6.64E-08, 9.11E-08, 6.73E-08, 0.4567, 6.65E-08, 2.21E-07)

# Adjust p-values using the Benjamini-Hochberg method
adjusted_p_values <- p.adjust(p_values, method = "BH")

# Print the adjusted p-values
print(adjusted_p_values)


