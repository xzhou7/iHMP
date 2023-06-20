#Generate Table For S9 and S10 

sample_intra <- read.csv("/Users/xzhou7/Library/CloudStorage/Box-Box/human_microbiome_project/data_analysis/correlation_network/whole_data_set/intra_microbiome_correlation/intra.csv", header = T)
sample_inter <- read.csv("/Users/xzhou7/Library/CloudStorage/Box-Box/human_microbiome_project/data_analysis/correlation_network/whole_data_set/intra_microbiome_correlation/inter.csv", header = T)

subject_intra <- read.csv("/Users/xzhou7/Library/CloudStorage/Box-Box/human_microbiome_project/data_analysis/correlation_network/whole_data_set/intra_microbiome_correlation/subject_intra.csv", header = T)
subject_inter <- read.csv("/Users/xzhou7/Library/CloudStorage/Box-Box/human_microbiome_project/data_analysis/correlation_network/whole_data_set/intra_microbiome_correlation/subject_inter.csv", header = T)

sample_intra$compare <- "within_bodysite"
sample_inter$compare <- "between_bodysite"
subject_intra$compare <- "within_bodysite"
subject_inter$compare <- "between_bodysite"


sample_wise_core <- rbind(sample_intra %>% select(from_Genus, to_Genus, from_class,  to_class, cor, p, p_adjust, compare),
                          sample_inter %>% select(from_Genus, to_Genus, from_class,  to_class, cor, p, p_adjust, compare))

subject_wise_core <- rbind(subject_intra %>% select(from_Genus, to_Genus, from_class,  to_class, cor, p_value, p_value_adjust, compare),
                           subject_inter %>% select(from_Genus, to_Genus, from_class,  to_class, cor, p_value, p_value_adjust, compare))

table(sample_intra$from_Genus == sample_intra$to_Genus)
table(subject_intra$from_Genus == subject_intra$to_Genus)

table(sample_inter$from_Genus == sample_inter$to_Genus)
table(subject_inter$from_Genus == subject_inter$to_Genus)

subject_intra[subject_intra$from_Genus == subject_intra$to_Genus,]


abs(sample_inter$p) %>% density() %>% plot()


abs(sample_wise_core$p_adjust) %>% max()

write.csv(file = "~/Desktop/1st.revision.CELL-S-22-04878/Supplementary_data/TableS9_samplewisecore.csv",sample_wise_core)
write.csv(file = "~/Desktop/1st.revision.CELL-S-22-04878/Supplementary_data/TableS10_subjectwisecore.csv",subject_wise_core)


sample_wise_core$cor

sort(sample_wise_core$cor, decreasing = TRUE)[1:10]

sample_wise_core[order(-sample_wise_core$cor), ][1:10, ]

sample_inter[order(-sample_inter$cor), ][1:10, ]

sample_wise_core[sample_wise_core$from_Genus == "Bacteroides",]

sample_wise_core[sample_wise_core$from_Genus == "Phocaeicola",]
sample_wise_core[sample_wise_core$to_Genus == "Phocaeicola",]

library(dplyr)

sample_wise_core %>% filter(from_Genus == "Bacteroides" | to_Genus == "Bacteroides")

print(filtered_rows)


sample_wise_core from

Phocaeicola
