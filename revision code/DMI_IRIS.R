#DMI and IRIS

library(phyloseq)

setwd("~/Library/CloudStorage/Box-Box/")
getwd()
source("./human_microbiome_project/code/tools.R")

load(file = "./XinZhouFiles/Projects/ZXP2_HMP2_Core/Analysis/Robject/DetailedPhyloseq.RData")

DMI_table <- read.csv("~/Desktop/IRIS_DMI/single_table_DMI_FS.csv", header = T)
participant.info <- read.csv(file = './XinZhouFiles/Projects/ZXP2_HMP2_Core/Analysis/ori_meta_table/metadata.subject.csv', header = T)

DMI_table_Stool_keep <- filter(DMI_table, bodysite=="Stool") 
DMI_table_Skin_keep <- filter(DMI_table, bodysite=="Skin")
DMI_table_Oral_keep <- filter(DMI_table, bodysite=="Oral")
DMI_table_Nasal_keep <- filter(DMI_table, bodysite=="Nasal")

#Stool data
stool.otu <- otu_table(physeqGenus_ST) %>% data.frame()
stool.otu[1:5,1:5]

stool.sample <- sample_data(physeqGenus_ST) %>% data.frame()
stool.sample

stool.taxa <- tax_table(physeqGenus_ST) %>% data.frame()
stool.taxa$ASV <- rownames(stool.taxa)
stool.sample

stool.otu.bysubject <- aggregate(stool.otu, by = list(stool.sample$SubjectID), FUN = mean)
stool.otu.bysubject[1:5, 1:20]

colnames(stool.otu.bysubject) <- stool.taxa$Genus[match(colnames(stool.otu.bysubject), stool.taxa$ASV)]
colnames(stool.otu.bysubject)[1] <- "SubjectID"

stool.DMI.genus <- select(stool.otu.bysubject,SubjectID, DMI_table_Stool_keep$genus)

DMI_table_Stool_keep

stool.long.DMI <- stool.DMI.genus %>%
  pivot_longer(-SubjectID, names_to = "genus", values_to = "abundance") %>% 
  inner_join(DMI_table_Stool_keep, by = "genus") %>% 
  mutate(abundance_dmi = abundance * dmi) %>%
  group_by(SubjectID) %>%
  summarise(total_abundance_dmi = sum(abundance_dmi)) %>% 
  left_join(select(participant.info, SubjectID, IRIS), by="SubjectID")

p.dmi.iris.stool <- ggplot(stool.long.DMI %>% filter(IRIS != "Unknown"), aes(x=IRIS, y=total_abundance_dmi)) + geom_boxplot() + stat_compare_means() + base_theme + ggtitle("Stool")
p.dmi.iris.stool
ggsave(filename = "~/Desktop/IRIS_DMI/p.dmi.iris.stool.pdf", p.dmi.iris.stool, width = 3, height = 4, dpi = 300)

wilcox.test(total_abundance_dmi ~ IRIS, data = stool.long.DMI %>% filter(IRIS != "Unknown"))

#Skin data
skin.otu <- otu_table(physeqGenus_SK) %>% data.frame()
skin.otu[1:5,1:5]

skin.sample <- sample_data(physeqGenus_SK) %>% data.frame()
skin.sample

skin.taxa <- tax_table(physeqGenus_SK) %>% data.frame()
skin.taxa$ASV <- rownames(skin.taxa)
skin.sample

skin.otu.bysubject <- aggregate(skin.otu, by = list(skin.sample$SubjectID), FUN = mean)
skin.otu.bysubject[1:5, 1:20]

colnames(skin.otu.bysubject) <- skin.taxa$Genus[match(colnames(skin.otu.bysubject), skin.taxa$ASV)]
colnames(skin.otu.bysubject)[1] <- "SubjectID"

skin.DMI.genus <- select(skin.otu.bysubject,SubjectID, DMI_table_Skin_keep$genus)

DMI_table_Skin_keep

skin.long.DMI <- skin.DMI.genus %>%
  pivot_longer(-SubjectID, names_to = "genus", values_to = "abundance") %>% 
  inner_join(DMI_table_Skin_keep, by = "genus") %>% 
  mutate(abundance_dmi = abundance * dmi) %>%
  group_by(SubjectID) %>%
  summarise(total_abundance_dmi = sum(abundance_dmi)) %>% 
  left_join(select(participant.info, SubjectID, IRIS), by="SubjectID")

p.dmi.iris.skin <- ggplot(skin.long.DMI %>% filter(IRIS != "Unknown"), aes(x=IRIS, y=total_abundance_dmi)) + geom_boxplot() + stat_compare_means() + base_theme + ggtitle("Skin")
p.dmi.iris.skin
ggsave(filename = "~/Desktop/IRIS_DMI/p.dmi.iris.skin.pdf", p.dmi.iris.skin, width = 3, height = 4, dpi = 300)

wilcox.test(total_abundance_dmi ~ IRIS, data = skin.long.DMI %>% filter(IRIS != "Unknown"))

#Oral data
oral.otu <- otu_table(physeqGenus_OR) %>% data.frame()
oral.otu[1:5,1:5]

oral.sample <- sample_data(physeqGenus_OR) %>% data.frame()
oral.sample

oral.taxa <- tax_table(physeqGenus_OR) %>% data.frame()
oral.taxa$ASV <- rownames(oral.taxa)
oral.sample

oral.otu.bysubject <- aggregate(oral.otu, by = list(oral.sample$SubjectID), FUN = mean)
oral.otu.bysubject[1:5, 1:20]

colnames(oral.otu.bysubject) <- oral.taxa$Genus[match(colnames(oral.otu.bysubject), oral.taxa$ASV)]
colnames(oral.otu.bysubject)[1] <- "SubjectID"

oral.DMI.genus <- select(oral.otu.bysubject,SubjectID, DMI_table_Oral_keep$genus)

DMI_table_Oral_keep

oral.long.DMI <- oral.DMI.genus %>%
  pivot_longer(-SubjectID, names_to = "genus", values_to = "abundance") %>% 
  inner_join(DMI_table_Oral_keep, by = "genus") %>% 
  mutate(abundance_dmi = abundance * dmi) %>%
  group_by(SubjectID) %>%
  summarise(total_abundance_dmi = sum(abundance_dmi)) %>% 
  left_join(select(participant.info, SubjectID, IRIS), by="SubjectID")

p.dmi.iris.oral <- ggplot(oral.long.DMI %>% filter(IRIS != "Unknown"), aes(x=IRIS, y=total_abundance_dmi)) + geom_boxplot() + stat_compare_means() + base_theme + ggtitle("Oral")
p.dmi.iris.oral
ggsave(filename = "~/Desktop/IRIS_DMI/p.dmi.iris.oral.pdf", p.dmi.iris.oral, width = 3, height = 4, dpi = 300)

wilcox.test(total_abundance_dmi ~ IRIS, data = oral.long.DMI %>% filter(IRIS != "Unknown"))

#Nasal data
nasal.otu <- otu_table(physeqGenus_NS) %>% data.frame()
nasal.otu[1:5,1:5]

nasal.sample <- sample_data(physeqGenus_NS) %>% data.frame()
nasal.sample

nasal.taxa <- tax_table(physeqGenus_NS) %>% data.frame()
nasal.taxa$ASV <- rownames(nasal.taxa)
nasal.sample

nasal.otu.bysubject <- aggregate(nasal.otu, by = list(nasal.sample$SubjectID), FUN = mean)
nasal.otu.bysubject[1:5, 1:20]

colnames(nasal.otu.bysubject) <- nasal.taxa$Genus[match(colnames(nasal.otu.bysubject), nasal.taxa$ASV)]
colnames(nasal.otu.bysubject)[1] <- "SubjectID"

nasal.DMI.genus <- select(nasal.otu.bysubject,SubjectID, DMI_table_Nasal_keep$genus)

DMI_table_Nasal_keep

nasal.long.DMI <- nasal.DMI.genus %>%
  pivot_longer(-SubjectID, names_to = "genus", values_to = "abundance") %>% 
  inner_join(DMI_table_Nasal_keep, by = "genus") %>% 
  mutate(abundance_dmi = abundance * dmi) %>%
  group_by(SubjectID) %>%
  summarise(total_abundance_dmi = sum(abundance_dmi)) %>% 
  left_join(select(participant.info, SubjectID, IRIS), by="SubjectID")

p.dmi.iris.nasal <- ggplot(nasal.long.DMI %>% filter(IRIS != "Unknown"), aes(x=IRIS, y=total_abundance_dmi)) + geom_boxplot() + stat_compare_means() + base_theme + ggtitle("Nasal")
p.dmi.iris.nasal
ggsave(filename = "~/Desktop/IRIS_DMI/p.dmi.iris.nasal.pdf", p.dmi.iris.nasal, width = 3, height = 4, dpi = 300)

wilcox.test(total_abundance_dmi ~ IRIS, data = nasal.long.DMI %>% filter(IRIS != "Unknown"))



# Stool data
stool.family_score.genus <- select(stool.otu.bysubject,SubjectID, DMI_table_Stool_keep$genus)

stool.long.family_score <- stool.family_score.genus %>%
  pivot_longer(-SubjectID, names_to = "genus", values_to = "abundance") %>% 
  inner_join(DMI_table_Stool_keep, by = "genus") %>% 
  mutate(abundance_family_score = abundance * family_score) %>%
  group_by(SubjectID) %>%
  summarise(total_abundance_family_score = sum(abundance_family_score)) %>% 
  left_join(select(participant.info, SubjectID, IRIS), by="SubjectID")

ggplot(stool.long.family_score %>% filter(IRIS != "Unknown"), aes(x=IRIS, y=total_abundance_family_score)) + geom_boxplot() + stat_compare_means() + base_theme

wilcox.test(total_abundance_family_score ~ IRIS, data = stool.long.family_score %>% filter(IRIS != "Unknown"))

# Skin data
skin.family_score.genus <- select(skin.otu.bysubject,SubjectID, DMI_table_Skin_keep$genus)

skin.long.family_score <- skin.family_score.genus %>%
  pivot_longer(-SubjectID, names_to = "genus", values_to = "abundance") %>% 
  inner_join(DMI_table_Skin_keep, by = "genus") %>% 
  mutate(abundance_family_score = abundance * family_score) %>%
  group_by(SubjectID) %>%
  summarise(total_abundance_family_score = sum(abundance_family_score)) %>% 
  left_join(select(participant.info, SubjectID, IRIS), by="SubjectID")

ggplot(skin.long.family_score %>% filter(IRIS != "Unknown"), aes(x=IRIS, y=total_abundance_family_score)) + geom_boxplot() + stat_compare_means() + base_theme

wilcox.test(total_abundance_family_score ~ IRIS, data = skin.long.family_score %>% filter(IRIS != "Unknown"))

# Oral data
oral.family_score.genus <- select(oral.otu.bysubject,SubjectID, DMI_table_Oral_keep$genus)

oral.long.family_score <- oral.family_score.genus %>%
  pivot_longer(-SubjectID, names_to = "genus", values_to = "abundance") %>% 
  inner_join(DMI_table_Oral_keep, by = "genus") %>% 
  mutate(abundance_family_score = abundance * family_score) %>%
  group_by(SubjectID) %>%
  summarise(total_abundance_family_score = sum(abundance_family_score)) %>% 
  left_join(select(participant.info, SubjectID, IRIS), by="SubjectID")

ggplot(oral.long.family_score %>% filter(IRIS != "Unknown"), aes(x=IRIS, y=total_abundance_family_score)) + geom_boxplot() + stat_compare_means() + base_theme

wilcox.test(total_abundance_family_score ~ IRIS, data = oral.long.family_score %>% filter(IRIS != "Unknown"))

#Nasal data
nasal.family_score.genus <- select(nasal.otu.bysubject,SubjectID, DMI_table_Nasal_keep$genus)

nasal.long.family_score <- nasal.family_score.genus %>%
  pivot_longer(-SubjectID, names_to = "genus", values_to = "abundance") %>% 
  inner_join(DMI_table_Nasal_keep, by = "genus") %>% 
  mutate(abundance_family_score = abundance * family_score) %>%
  group_by(SubjectID) %>%
  summarise(total_abundance_family_score = sum(abundance_family_score)) %>% 
  left_join(select(participant.info, SubjectID, IRIS), by="SubjectID")

ggplot(nasal.long.family_score %>% filter(IRIS != "Unknown"), aes(x=IRIS, y=total_abundance_family_score)) + geom_boxplot() + stat_compare_means() + base_theme

wilcox.test(total_abundance_family_score ~ IRIS, data = nasal.long.family_score %>% filter(IRIS != "Unknown"))




