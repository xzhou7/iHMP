#Check A1C with Akkermansia
library(phyloseq)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggstatsplot)
library(ggpubr)
library(patchwork)
library(cowplot)
library(ggrepel)


setwd("~/Library/CloudStorage/Box-Box/XinZhouFiles/Projects/ZXP2_HMP2_Core/Analysis/")

load("./Robject/DetailedPhyloseq.RData")
load("./Robject/Revision_MultiOmes_0509.RData")

otu_table <- otu_table(physeqGenus_ST) %>% data.frame()
meta_sample <- sample_data(physeqGenus_ST) %>% data.frame()
tax_table <- tax_table(physeqGenus_ST) %>% data.frame()
tax_table$ASV <- row.names(tax_table)

meta_sample_SampleID <- left_join(clinic.df, meta_sample, by=c("SampleID"= "SampleID"))
colnames(otu_table) <- tax_table$Genus[match(colnames(otu_table), tax_table$ASV)]
otu_table <- otu_table %>% select(-Unclassified_Actinobacteria)
otu_table$RandomID <- row.names(otu_table)
stool.OTU_meta <- left_join(meta_sample_SampleID, otu_table, by = "RandomID")
stool.OTU_meta[1:5, 1:5]

colnames(stool.OTU_meta)[str_detect(colnames(stool.OTU_meta),"Akk")]

stool.OTU_meta$CollectionDate <- as.Date(stool.OTU_meta$CollectionDate, format = "%m/%d/%y")

ggplot(stool.OTU_meta, aes(x=A1C, y=Akkermansia, group=SubjectID.x)) + geom_point()

# #########################
# individual.list <- unique(stool.OTU_meta$SubjectID.x)
# 
# for (i in 1:104){
# individual.label <- individual.list[i]
# print(individual.label)
# 
# title <- paste(individual.label, "A1C and Akkermansia relationship", sep=" ")
# 
# p.a1c <- stool.OTU_meta %>% filter(SubjectID.x %in% individual.label) %>%
#   ggplot(aes(x=CollectionDate)) +
#   geom_line(aes(y=A1C), size=2, color="red") +
#   theme_cowplot() + theme(axis.title.x = element_blank()) + ggtitle(title)
# p.a1c
# 
# ## p.akk <- stool.OTU_meta %>% filter(SubjectID.x %in% individual.label) %>% select(CollectionDate,SubjectID.x,A1C,Akkermansia) %>% #na.omit() %>%
# ##   ggplot(aes(x=CollectionDate,y=SubjectID.x)) +
# ##   geom_tile(aes(fill=Akkermansia), size=2) + scale_fill_continuous(low="#F39B7FFF", high="red",na.value="white") + theme_cowplot() + 
# ##   theme(axis.title.y = element_blank())
# ## p.akk
# 
# a1c.akk <- stool.OTU_meta %>% filter(SubjectID.x %in% individual.label) %>% select(CollectionDate,SubjectID.x,A1C,Akkermansia) 
# a1c.akk$Akkermansia[is.na(a1c.akk$Akkermansia)] <- 0
# a1c.akk$Akkermansia <- a1c.akk$Akkermansia +1
# 
# p.akk <-  a1c.akk %>%
#   ggplot(aes(x=CollectionDate, y=Akkermansia, fill="red"))+
#   geom_bar(stat = "identity") + theme_cowplot() +
#   theme(axis.title.y = element_blank()) + scale_y_continuous(trans='log10')
# p.akk
# 
# p.temp <- p.a1c/p.akk +  plot_layout(heights  = c(2, 1))
# 
# path <- paste0("./Suppl.figure/A1C_Akker/",individual.label, ".pdf")
# ggsave(path, p.temp, width = 7, height = 5, dpi = 300)
# }

#######check mean value
stool.OTU_meta %>% select(SubjectID.x,A1C,Akkermansia) %>%
  ggplot(aes(x=A1C,y=Akkermansia)) + geom_point() + geom_smooth(method = "lm") + geom_vline(xintercept = 5.7) + geom_vline(xintercept = 6.4 )

stool.OTU_meta %>% select(SubjectID.x,A1C,Akkermansia) %>% 
  group_by(SubjectID.x) %>% mutate(mean_A1C = mean(A1C)) %>% na.omit() %>%
  mutate(mean_Akk = mean(Akkermansia)) %>% 
  ggplot(aes(x=mean_A1C, y=mean_Akk)) + geom_point() + geom_smooth(method="lm")


stool.AKK_A1C.by.subject <- stool.OTU_meta %>% select(SubjectID.x, Akkermansia, A1C, IRIS) %>%
  group_by(SubjectID.x) %>% mutate(mean_A1C = mean(A1C)) %>% na.omit() %>%
  mutate(mean_Akk = mean(Akkermansia)) %>% select(SubjectID.x, mean_Akk, mean_A1C, IRIS) %>% unique()

stool.AKK_A1C.by.subject
p.individual <- ggplot(stool.AKK_A1C.by.subject, aes(x=mean_A1C, y=mean_Akk, color=IRIS, label=SubjectID.x)) + geom_point() + geom_smooth(method = "lm") + theme_cowplot() +
  geom_label_repel(data = filter(stool.AKK_A1C.by.subject, mean_Akk > 0.02), aes(x=mean_A1C, y=mean_Akk, color=IRIS, label=SubjectID.x))
p.individual

ggsave(filename = "./Suppl.figure/A1C_Akker/01.BY.INDIVIDUAl.pdf",p.individual, width = 7, height = 4, dpi=300)

table.to.save <- stool.AKK_A1C.by.subject %>% arrange(desc(mean_Akk))
write.csv(file = "./Suppl.figure/A1C_Akker/data.individual.csv", table.to.save)

########################################################################
individual.label <- "69-113"
title <- paste(individual.label, "A1C and Akkermansia relationship", sep=" ")
p.a1c <- stool.OTU_meta %>% filter(SubjectID.x %in% individual.label) %>%
  ggplot(aes(x=CollectionDate)) +
 geom_line(aes(y=A1C), size=2, color="red") +
 theme_cowplot() + theme(axis.title.x = element_blank()) + ggtitle(title)
p.a1c

a1c.akk <- stool.OTU_meta %>% filter(SubjectID.x %in% individual.label) %>% select(CollectionDate,SubjectID.x,A1C,Akkermansia) 
a1c.akk$Akkermansia[is.na(a1c.akk$Akkermansia)] <- 0
a1c.akk$Akkermansia <- a1c.akk$Akkermansia +1
a1c.akk

p.akk <-  a1c.akk %>%
ggplot(aes(x=CollectionDate, y=Akkermansia, fill="red"))+
geom_bar(stat = "identity", width = 1.5) + theme_cowplot() +
theme(axis.title.y = element_blank()) + scale_y_continuous(trans='log10')
p.akk

p.temp <- p.a1c/p.akk +  plot_layout(heights  = c(2, 1))
p.temp

########################################################################
p.akk.iris <-  stool.AKK_A1C.by.subject %>% filter(IRIS != "Unknown") %>% ggbetweenstats(IRIS, mean_Akk)
p.akk.iris

ggsave(filename = "./Suppl.figure/A1C_Akker/IRIS.pdf",p.akk.iris, width = 7, height = 5, dpi = 300 )


