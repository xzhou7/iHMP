setwd("~/Library/CloudStorage/Box-Box/XinZhouFiles/Projects/ZXP2_HMP2_Core/")
library(phyloseq)
library(ggplot2)
library(dplyr)
library(cowplot)
library(patchwork)

body_site_color = c(
  "Stool" = ggsci::pal_jama()(n=7)[2],
  "Skin" = ggsci::pal_jama()(n=7)[3],
  "Oral" = ggsci::pal_jama()(n=7)[4],
  "Nasal" = ggsci::pal_jama()(n=7)[5])

load("./Analysis/Robject/DetailedPhyloseq.RData")

physeq.data.st <- data.frame(otu_table(physeq_ST))
physeq.sample.st <- data.frame(sample_data(physeq_ST))

physeq.data.sk <- data.frame(otu_table(physeq_SK))
physeq.sample.sk <- data.frame(sample_data(physeq_SK))

physeq.data.or <- data.frame(otu_table(physeq_OR))
physeq.sample.or <- data.frame(sample_data(physeq_OR))

physeq.data.ns <- data.frame(otu_table(physeq_NS))
physeq.sample.ns <- data.frame(sample_data(physeq_NS))

#check if OTU table and sample table are in the same order
identical(rownames(physeq.data),physeq.sample$RandomID)

ST.physeq <- data.frame(rowSums(physeq.data.st))
colnames(ST.physeq) <- "depth"
ST.physeq$bodysite <- "Stool"
SK.physeq <- data.frame(rowSums(physeq.data.sk))
colnames(SK.physeq) <- "depth"
SK.physeq$bodysite <- "Skin"
OR.physeq <- data.frame(rowSums(physeq.data.or))
colnames(OR.physeq) <- "depth"
OR.physeq$bodysite <- "Oral"
NS.physeq <- data.frame(rowSums(physeq.data.ns))
colnames(NS.physeq) <- "depth"
NS.physeq$bodysite <- "Nasal"

depth <- rbind(ST.physeq,SK.physeq,OR.physeq,NS.physeq)
#write.csv(file = "./Analysis/ori_meta_table/SequencingDepth.csv", depth)
p1 <- filter(depth, bodysite %in% c("Stool", "Nasal")) %>% ggplot(aes(x=log10(depth), color=bodysite)) + geom_density() + geom_vline(xintercept = log10(7000))
p1 <- p1 +  geom_vline(xintercept = log10(70000)) + theme_cowplot() + scale_color_manual(values = body_site_color)
mean(ST.physeq$depth)
sd(ST.physeq$depth)
p1 <- p1 + annotate("text", x=3, y=1.3,hjust = 0, label= "ST: mean_23554 \n       sd_13548")
mean(NS.physeq$depth)
sd(NS.physeq$depth)
p1 <- p1 + annotate("text", x=3, y=0.9, hjust = 0, label= "NS: mean_24899 \n       sd_12456")
p1

p2 <- filter(depth, bodysite %in% c("Skin", "Oral")) %>% ggplot(aes(x=log10(depth), color=bodysite)) + geom_density() + geom_vline(xintercept = log10(7000))
p2 <- p2 +  geom_vline(xintercept = log10(500000)) + theme_cowplot() + scale_color_manual(values = body_site_color)
mean(SK.physeq$depth)
sd(SK.physeq$depth)
p2 <- p2 + annotate("text", x=3, y=0.9,hjust = 0, label= "SK: mean_74517 \n       sd_119445")
mean(OR.physeq$depth)
sd(OR.physeq$depth)
p2 <- p2 + annotate("text", x=3, y=0.6,hjust = 0, label= "OR: mean_132912 \n       sd_108622")
p2

sequencing_depth <- p1 / p2
sequencing_depth
#ggsave(filename = "./Analysis/Suppl.figure/Depth.pdf",sequencing_depth, height = 6, width = 8, dpi = 300)

