#Overview of sample collected
library(phyloseq)
library(stringr)
library(cowplot)
library(ggfortify)
library(vegan)
library(tidyverse)
library(microbiome)
library(patchwork)

setwd("/Users/xzhou7/Library/CloudStorage/Box-Box/XinZhouFiles/Projects/ZXP2_HMP2_Core/Analysis/")

load("./Robject/PhyloseqObject.RData")
load("./Robject/Revision_MultiOmes_0509.RData")
sc <- read.csv("./ori_meta_table/metadata.subject.csv", header = T, row.names = 1)

HMP_meta <- data.frame(sample_data(physeq_HMP))
UBM_meta <- data.frame(sample_data(physeq_UBM))

HMP_meta_ST <- subset(HMP_meta, SampleType == "ST")
HMP_meta_NS <- subset(HMP_meta, SampleType == "NS")

UBM_meta_OR <- subset(UBM_meta, SampleType == "Oral")
UBM_meta_SK <- subset(UBM_meta, SampleType == "Skin")

HMP_meta2 <- select(HMP_meta, RandomID,SampleID,SampleType,batch,Date,SubjectID)
colnames(HMP_meta2)[1] <- "KitID"
Combined_meta <- rbind(HMP_meta2, UBM_meta)

Combined_meta$RandomSubjectID <- sc$rand_subject_id[match(Combined_meta$SubjectID, sc$SubjectID)]
p <- ggplot(Combined_meta, aes(x=Date, y=RandomSubjectID, color=SampleType)) + geom_point() + geom_jitter()
p <- p +theme_cowplot() + facet_wrap(.~SampleType, ncol=4)
p

#ggsave(filename = "./Suppl.figure/subject.before.trim.pdf", p, height = 14, width = 8, dpi=300)

sampleFreq <- merge(merge(data.frame(table(HMP_meta_ST$SubjectID)),data.frame(table(HMP_meta_NS$SubjectID)), by="Var1", all=T), 
                     merge(data.frame(table(UBM_meta_OR$SubjectID)),data.frame(table(UBM_meta_SK$SubjectID)), by="Var1",all=T), by= "Var1",all=T) 
  
colnames(sampleFreq) <- c("SubjectID", "ST", "NS", "OR", "SK")
sampleFreq[is.na(sampleFreq)] <- 0

#################################################################################
#sample_data(physeq_HMP)$batch[sample_data(physeq_HMP)$RandomID == "BCIXQPED"] <- "B86HB"

physeq_ST <- subset_samples(physeq_HMP, SampleType == "ST")
physeq_ST_freq <- transform_sample_counts(physeq_ST, function(x) x / sum(x))

physeq_ord <- ordinate(physeq_ST_freq, "PCoA", dist = "bray")
physeq69001 <- subset_samples(physeq_ST_freq, SubjectID %in% c("69-001", "69-031"))

sample_data(physeq_ST_freq)$ASV31 <- as.numeric(otu_table(physeq_ST_freq)[,31])

sample_data(physeq_ST_freq)$IRIS <- sc$IRIS[match(sample_data(physeq_ST_freq)$Subject,sc$SubjectID)]

p1.1 <- plot_ordination(physeq_ST_freq, physeq_ord, type = "samples",color = "batch",  axes = c(1,2))  
p1.1 <- p1.1 + geom_point(size = 0.5) + theme_cowplot()
p1.1

#find out who is the outlier
sample <- rownames(physeq_ord$vectors)[physeq_ord$vectors[,2] < -0.2]
HMP_meta_ST[HMP_meta_ST$RandomID %in% sample,]

physeqGenus_ST <- tax_glom(physeq_ST_freq, taxrank = "Genus")
physeq_ord_G <- ordinate(physeqGenus_ST, "PCoA", dist = "bray")

sample_data(physeqGenus_ST)$Bacteroides <- as.numeric(otu_table(physeqGenus_ST)[,5])
sample_data(physeqGenus_ST)$Prevotella <- as.numeric(otu_table(physeqGenus_ST)[,9])
sample_data(physeqGenus_ST)$Agathobacter <- as.numeric(otu_table(physeqGenus_ST)[,8])
sample_data(physeqGenus_ST)$Ruminococcus <- as.numeric(otu_table(physeqGenus_ST)[,13])

p1.2 <- plot_ordination(physeq_ST_freq, physeq_ord_G, type = "samples", axes = c(1,2),color = "SubjectID") 
p1.2 <- p1.2 + geom_point(size = 0.5) + theme_cowplot()
p1.2

p1.3 <- plot_ordination(physeqGenus_ST, physeq_ord_G, type = "samples", axes = c(1,2),color = "Bacteroides") 
p1.3 <- p1.3 + geom_point(size = 0.01) + theme_cowplot() + scale_colour_gradient2(low = "grey", mid="orange", high="red",midpoint=0.5)
p1.3

p1.4 <- plot_ordination(physeqGenus_ST, physeq_ord_G, type = "samples", axes = c(1,2),color = "Prevotella") 
p1.4 <- p1.4 + geom_point(size = 0.01) + theme_cowplot() + scale_colour_gradient2(low = "grey", mid="orange", high="red",midpoint=0.5)
p1.4

p1.5 <- plot_ordination(physeqGenus_ST, physeq_ord_G, type = "samples", axes = c(1,2),color = "Agathobacter") 
p1.5 <- p1.5 + geom_point(size = 0.01) + theme_cowplot() + scale_colour_gradient2(low = "grey", mid="orange", high="red",midpoint=0.15)
p1.5

p1.6 <- plot_ordination(physeqGenus_ST, physeq_ord_G, type = "samples", axes = c(1,2),color = "Ruminococcus") 
p1.6 <- p1.6 + geom_point(size = 0.01) + theme_cowplot() + scale_colour_gradient2(low = "grey", mid="orange", high="red",midpoint=0.3)
p1.6

(p1.3 + p1.4) / (p1.5 + p1.6)

tax_table(physeqGenus_ST)[1:20]
otu_table(physeqGenus_ST)[1:10, 1:14]
colSums(otu_table(physeqGenus_ST))[1:20]
otu_table(physeq_ST_freq)[1:10, 30:35]
otu_table(physeqGenus_ST)[,6]

p1.d <- ggplot(otu_table(physeqGenus_ST), aes(x=ASV12)) + stat_density(geom="line", position = "identity") +xlim(0.00000000001,1)# coord_cartesian(xlim = c(0.1, 0.3))
p1.d <- p1.d  + stat_density(data=otu_table(physeqGenus_ST), mapping=aes(x=ASV5), geom="line", position = "identity") +xlim(0.00000000001,1)
p1.d <- p1.d  + stat_density(data=otu_table(physeqGenus_ST), mapping=aes(x=ASV26), geom="line", position = "identity") +xlim(0.00000000001,1)
p1.d

tax_table(physeqGenus_ST)

##########################################################################################################################
physeq_SK <- subset_samples(physeq_UBM, SampleType == "Skin")
physeq_SK_freq <- transform_sample_counts(physeq_SK,function(x) x / sum(x))

physeq_ord2 <- ordinate(physeq_SK_freq, "PCoA", dist = "bray")
#physeq_ord2_jsd <- ordinate(physeq_SK_freq, "PCoA", dist = "jsd")

p2.1 <- plot_ordination(physeq_SK_freq, physeq_ord2, type = "samples",color="batch", axes = c(1,2)) 
p2.1 <- p2.1  + geom_point(size = 0.5) + theme_cowplot()
p2.1

physeqGenus_SK <- tax_glom(physeq_SK_freq, taxrank = "Genus")
physeq_ord2_G <- ordinate(physeqGenus_SK, "PCoA", dist = "bray")

p2.2 <- plot_ordination(physeqGenus_SK, physeq_ord2_G, type = "samples",color="SubjectID", axes = c(1,2)) 
p2.2 <- p2.2  + geom_point(size = 0.5) + theme_cowplot()
p2.2

otu_table(physeqGenus_SK)[1:5,30:40]
tax_table(physeqGenus_SK)

sample_data(physeqGenus_SK)$Cutibacterium <- as.numeric(otu_table(physeqGenus_SK)[,1])
sample_data(physeqGenus_SK)$Neisseria <- as.numeric(otu_table(physeqGenus_SK)[,2])
sample_data(physeqGenus_SK)$Staphylococcus <- as.numeric(otu_table(physeqGenus_SK)[,3])
sample_data(physeqGenus_SK)$Corynebacterium <- as.numeric(otu_table(physeqGenus_SK)[,13])
sample_data(physeqGenus_SK)$Moraxella <- as.numeric(otu_table(physeqGenus_SK)[,40])

p2.3 <- plot_ordination(physeqGenus_SK, physeq_ord2_G, type = "samples",color="Cutibacterium", axes = c(1,2)) 
p2.3 <- p2.3 + geom_point(size = 0.01) + theme_cowplot() + scale_colour_gradient2(low = "grey", mid="orange", high="red",midpoint=0.6)
p2.3

p2.4 <- plot_ordination(physeqGenus_SK, physeq_ord2_G, type = "samples",color="Corynebacterium", axes = c(1,2)) 
p2.4 <- p2.4 + geom_point(size = 0.01) + theme_cowplot() + scale_colour_gradient2(low = "grey", mid="orange", high="red",midpoint=0.4)
p2.4

p2.5 <- plot_ordination(physeqGenus_SK, physeq_ord2_G, type = "samples",color="Neisseria", axes = c(1,2)) 
p2.5 <- p2.5 + geom_point(size = 0.01) + theme_cowplot() + scale_colour_gradient2(low = "grey", mid="orange", high="red",midpoint=0.4)
p2.5

p2.6 <- plot_ordination(physeqGenus_SK, physeq_ord2_G, type = "samples",color="Staphylococcus", axes = c(1,2)) 
p2.6 <- p2.6 + geom_point(size = 0.01) + theme_cowplot() + scale_colour_gradient2(low = "grey", mid="orange", high="red",midpoint=0.4)
p2.6

p2.7 <- plot_ordination(physeqGenus_SK, physeq_ord2_G, type = "samples",color="Moraxella", axes = c(1,2)) 
p2.7 <- p2.7 + geom_point(size = 0.01) + theme_cowplot() + scale_colour_gradient2(low = "grey", mid="orange", high="red",midpoint=0.015)
p2.7

p2.3 + p2.4 + p2.5 + p2.6

##########################################################################################################################
physeq_OR <- subset_samples(physeq_UBM, SampleType == "Oral")
physeq_OR_freq <- transform_sample_counts(physeq_OR,function(x) x / sum(x))
physeq_ord3 <- ordinate(physeq_OR_freq, "PCoA", dist = "bray")

p3.1 <- plot_ordination(physeq_OR_freq, physeq_ord3, type = "samples",color="SubjectID", axes = c(1,2)) 
p3.1 <- p3.1 + geom_point(size = 0.5) + theme_cowplot()
p3.1

physeqGenus_OR <- tax_glom(physeq_OR_freq, taxrank = "Genus")
physeq_ord3_G <- ordinate(physeqGenus_OR, "PCoA", dist = "bray")

otu_table(physeqGenus_OR)[1:10, 1:10]
colSums(otu_table(physeqGenus_OR))[1:10]
tax_table(physeq_OR_freq)[1:10]

sample_data(physeqGenus_OR)$Cutibacterium <- as.numeric(otu_table(physeqGenus_OR)[,2])
sample_data(physeqGenus_OR)$Prevotella <- as.numeric(otu_table(physeqGenus_OR)[,4])
sample_data(physeqGenus_OR)$Veillonella <- as.numeric(otu_table(physeqGenus_OR)[,5])
sample_data(physeqGenus_OR)$Haemophilus <- as.numeric(otu_table(physeqGenus_OR)[,6])

p3.2 <- plot_ordination(physeqGenus_OR, physeq_ord3_G, type = "samples",color="SubjectID", axes = c(1,2)) 
p3.2 <- p3.2 + geom_point(size = 0.5) + theme_cowplot()
p3.2

p3.3 <- plot_ordination(physeqGenus_OR, physeq_ord3_G, type = "samples",color="Cutibacterium", axes = c(1,2)) 
p3.3 <- p3.3 + geom_point(size = 0.01) + theme_cowplot() + scale_colour_gradient2(low = "grey", mid="orange", high="red",midpoint=0.4)
p3.3

p3.4 <- plot_ordination(physeqGenus_OR, physeq_ord3_G, type = "samples",color="Prevotella", axes = c(1,2)) 
p3.4 <- p3.4 + geom_point(size = 0.01) + theme_cowplot() + scale_colour_gradient2(low = "grey", mid="orange", high="red",midpoint=0.3)
p3.4

p3.5 <- plot_ordination(physeqGenus_OR, physeq_ord3_G, type = "samples",color="Veillonella", axes = c(1,2)) 
p3.5 <- p3.5 + geom_point(size = 0.01) + theme_cowplot() + scale_colour_gradient2(low = "grey", mid="orange", high="red",midpoint=0.3)
p3.5

p3.6 <- plot_ordination(physeqGenus_OR, physeq_ord3_G, type = "samples",color="Haemophilus", axes = c(1,2)) 
p3.6 <- p3.6 + geom_point(size = 0.01) + theme_cowplot() + scale_colour_gradient2(low = "grey", mid="orange", high="red",midpoint=0.4)
p3.6

p3.3 + p3.4 + p3.5 + p3.6

##########################################################################################################################
physeq_NS <- subset_samples(physeq_HMP, SampleType=="NS")
physeq_NS_freq <- transform_sample_counts(physeq_NS, function(x) x / sum(x))
physeq_ord4 <- ordinate(physeq_NS_freq, "PCoA", dist = "bray")

p4.1 <- plot_ordination(physeq_NS_freq, physeq_ord4, type = "sample",color="batch", axes = c(1,2)) 
p4.1 <- p4.1 + geom_point(size = 0.5) + theme_cowplot()
p4.1

physeqGenus_NS <- tax_glom(physeq_NS_freq, taxrank = "Genus")
physeq_ord4_G <- ordinate(physeqGenus_NS, "PCoA", dist = "bray")

p4.2 <- plot_ordination(physeqGenus_NS, physeq_ord4_G, type = "samples",color="SubjectID", axes = c(1,2)) 
p4.2 <- p4.2 + geom_point(size = 0.5) + theme_cowplot()
p4.2

tax_table(physeqGenus_NS)
otu_table(physeqGenus_NS)[1:10, 1:10]

colSums(otu_table(physeqGenus_NS))[1:10]

sample_data(physeqGenus_NS)$Cutibacterium <- as.numeric(otu_table(physeqGenus_NS)[,1])
sample_data(physeqGenus_NS)$Corynebacterium <- as.numeric(otu_table(physeqGenus_NS)[,2])
sample_data(physeqGenus_NS)$Staphylococcus <- as.numeric(otu_table(physeqGenus_NS)[,3])
sample_data(physeqGenus_NS)$Moraxella <- as.numeric(otu_table(physeqGenus_NS)[,5])

p4.3 <- plot_ordination(physeqGenus_NS, physeq_ord4_G, type = "samples", color="Cutibacterium", axes = c(1,2)) 
p4.3 <- p4.3 + geom_point(size = 0.01) + theme_cowplot() + scale_colour_gradient2(low = "grey", mid="orange", high="red",midpoint=0.4)
p4.3

p4.4 <- plot_ordination(physeqGenus_NS, physeq_ord4_G, type = "samples",color="Corynebacterium", axes = c(1,2)) 
p4.4 <- p4.4 + geom_point(size = 0.01) + theme_cowplot() + scale_colour_gradient2(low = "grey", mid="orange", high="red",midpoint=0.4)
p4.4

p4.5 <- plot_ordination(physeqGenus_NS, physeq_ord4_G, type = "samples",color="Staphylococcus", axes = c(1,2)) 
p4.5 <- p4.5 + geom_point(size = 0.01) + theme_cowplot() + scale_colour_gradient2(low = "grey", mid="orange", high="red",midpoint=0.4)
p4.5

p4.6 <- plot_ordination(physeqGenus_NS, physeq_ord4_G, type = "samples",color="Moraxella", axes = c(1,2)) 
p4.6 <- p4.6 + geom_point(size = 0.01) + theme_cowplot() + scale_colour_gradient2(low = "grey", mid="orange", high="red",midpoint=0.4)
p4.6

p4.3 + p4.4 + p4.5 + p4.6

save(physeq_NS, physeq_NS_freq, physeqGenus_SK, physeqGenus_NS,physeqGenus_ST, physeqGenus_OR,
     physeq_OR, physeq_OR_freq, physeq_SK, physeq_SK_freq, physeq_ST, physeq_ST_freq, file = "./Robject/DetailedPhyloseq.RData")

############################################################################################################
#permanova
df = as(sample_data(physeq_ST_freq), "data.frame")
d = phyloseq::distance(physeq_ST_freq,"bray", "sample")
physeq.ST.adonis = adonis(d ~ SubjectID + batch, df, permutations= 999)
physeq.ST.adonis

df2 = as(sample_data(physeq_SK_freq), "data.frame")
d2 = phyloseq::distance(physeq_SK_freq, "bray", "sample")
physeq.SK.adonis = adonis(d2 ~ SubjectID + batch, df2, permutations= 999)
physeq.SK.adonis

df3 = as(sample_data(physeq_OR_freq), "data.frame")
d3 = phyloseq::distance(physeq_OR_freq, "bray", "sample")
physeq.OR.adonis = adonis(d3 ~ SubjectID + batch, df3, permutations= 999)
physeq.OR.adonis

df4 = as(sample_data(physeq_NS_freq), "data.frame")
d4 = phyloseq::distance(physeq_NS_freq, "bray", "sample")
physeq.NS.adonis = adonis(d4 ~ SubjectID + batch, df4, permutations= 999)
physeq.NS.adonis

#######################################variation plot################################
physeqGenus_NS_sample <- as(sample_data(physeqGenus_NS), "data.frame")
colnames(physeqGenus_NS_sample)

#plot Cutibacterium
mean_B <- aggregate(physeqGenus_NS_sample$Moraxella, list(physeqGenus_NS_sample$SubjectID), mean)
colnames(mean_B) <- c("SubjectID","mean_B")
NS_sample <- merge(physeqGenus_NS_sample, merge(mean_B, sc, by="SubjectID"), by="SubjectID")  %>% 
  group_by(SubjectID) %>% filter(n() >= 3)
colnames(NS_sample)

p.ns.r1 <- ggplot(NS_sample , aes(x=log10(Moraxella), y=fct_reorder(SubjectID,mean_B), color=IRIS, group=SubjectID)) + geom_point(size=0.1) + geom_line()
p.ns.r1 <- p.ns.r1 + theme_cowplot() + scale_color_manual(values=c("red", "grey", "grey")) + theme(panel.grid = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.x = element_text(face="italic"))
p.ns.r1 <- p.ns.r1 + ylab("Subject") 
p.ns.r1

p.ns.r2 <- ggplot(NS_sample , aes(x=log(Moraxella), y=fct_reorder(SubjectID,mean_B), color=Class, group=SubjectID)) + geom_point(size=0.1) + geom_line()
p.ns.r2 <- p.ns.r2 + theme_cowplot() + scale_color_manual(values=c("grey", "blue", "red","orange")) + theme(panel.grid = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.x = element_text(face="italic"))
p.ns.r2 <- p.ns.r2 + ylab("Subject")
p.ns.r2

p.ns.r3 <- ggplot(NS_sample , aes(x=log10(Moraxella), y=fct_reorder(SubjectID, as.numeric(BMI)), color=IRIS, group=SubjectID)) + geom_point(size=0.1) + geom_line()
p.ns.r3 <- p.ns.r3 + theme_cowplot() +  scale_color_manual(values=c("red", "grey", "grey","orange")) + theme(panel.grid = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.x = element_text(face="italic"))
p.ns.r3 <- p.ns.r3 + ylab("Subject") + geom_hline(yintercept="69-036", color="red", linetype=3) + geom_hline(yintercept="69-123", color="red", linetype=3)
p.ns.r3

physeqGenus_NS_sample <- as(sample_data(physeqGenus_NS), "data.frame")
colnames(physeqGenus_NS_sample)

##############################################################################################################
#plot stool
sample_data(physeqGenus_ST)$Akkermansia <- as.numeric(otu_table(physeqGenus_ST)[,18])

physeqGenus_ST_sample <- as(sample_data(physeqGenus_ST), "data.frame")
colnames(physeqGenus_ST_sample)

mean_A <- aggregate(physeqGenus_ST_sample$Akkermansia, list(physeqGenus_ST_sample$SubjectID), mean)
colnames(mean_A) <- c("SubjectID","mean_A")
ST_sample <- merge(physeqGenus_ST_sample, merge(mean_A, sc, by="SubjectID"), by="SubjectID")  %>% 
  group_by(SubjectID) %>% filter(n() >= 3)
colnames(ST_sample)

p.ST.r1 <- ggplot(ST_sample , aes(x=log10(Akkermansia), y=fct_reorder(SubjectID,mean_A), color=IRIS.x, group=SubjectID)) + geom_point(size=0.1) + geom_line()
p.ST.r1 <- p.ST.r1 + theme_cowplot() + scale_color_manual(values=c("red", "grey", "grey")) + theme(panel.grid = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.x = element_text(face="italic"))
p.ST.r1 <- p.ST.r1 + ylab("Subject") 
p.ST.r1

p.ST.r2 <- ggplot(ST_sample , aes(x=log(Akkermansia), y=fct_reorder(SubjectID,mean_A), color=Class, group=SubjectID)) + geom_point(size=0.1) + geom_line()
p.ST.r2 <- p.ST.r2 + theme_cowplot() + scale_color_manual(values=c("grey", "blue", "red","orange")) + theme(panel.grid = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.x = element_text(face="italic"))
p.ST.r2 <- p.ST.r2 + ylab("Subject")
p.ST.r2

p.ST.r3 <- ggplot(ST_sample , aes(x=log10(Akkermansia), y=fct_reorder(SubjectID, as.numeric(BMI)), color=IRIS.x, group=SubjectID)) + geom_point(size=0.1) + geom_line()
p.ST.r3 <- p.ST.r3 + theme_cowplot() +  scale_color_manual(values=c("red", "grey", "grey","orange")) + theme(panel.grid = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.x = element_text(face="italic"))
p.ST.r3 <- p.ST.r3 + ylab("Subject") + geom_hline(yintercept="69-036", color="red", linetype=3) + geom_hline(yintercept="69-123", color="red", linetype=3)
p.ST.r3




####################estimate diversity###########################################
physeq_HMP_genus <- tax_glom(physeq_HMP, taxrank = "Genus")
physeq_UBM_genus <- tax_glom(physeq_UBM, taxrank = "Genus")

HMP_richness <- estimate_richness(physeq_HMP)
UBM_richness <- estimate_richness(physeq_UBM)

d1 <-  data.frame(HMP_richness, sample_data(physeq_HMP))
d2 <- data.frame(UBM_richness, sample_data(physeq_UBM))
d3 <- data.frame(estimate_richness(physeq_HMP_genus), sample_data(physeq_HMP))
d4 <- data.frame(estimate_richness(physeq_UBM_genus), sample_data(physeq_UBM))

d1$depth <- rowSums(otu_table(physeq_HMP))
d2$depth <- rowSums(otu_table(physeq_UBM))

d3$depth <- rowSums(otu_table(physeq_HMP))
d4$depth <- rowSums(otu_table(physeq_UBM))

d1$LOF <- select(d1, depth, Observed) %>% lof(minPts = 21)
plot(1:20, arrange(d1,desc(LOF))$LOF[1:20])

d2$LOF <- select(d2, depth, Observed) %>% lof(minPts = 21)
plot(1:20, arrange(d2,desc(LOF))$LOF[1:20])

pdi.1 <- ggplot(d1, aes(x=log10(depth),y=Observed, color=LOF)) + geom_point()
pdi.1 <- pdi.1 + facet_wrap(.~SampleType) + geom_smooth(method = "loess")
pdi.1

pdi.2 <- ggplot(d1, aes(x=log(depth), y=Shannon, color=LOF)) + geom_point()
pdi.2 <- pdi.2 + facet_wrap(.~SampleType) + geom_smooth(method = "loess")
pdi.2

pdi.3 <- ggplot(d2, aes(x=log10(depth),y=Observed, color=LOF)) + geom_point()
pdi.3 <- pdi.3 + facet_wrap(.~SampleType) + geom_smooth(method = "loess")
pdi.3

pdi.4 <- ggplot(d2, aes(x=log(depth), y=Shannon, color=LOF)) + geom_point()
pdi.4 <- pdi.4 + facet_wrap(.~SampleType)  + geom_smooth(method = "loess")
pdi.4

(pdi.1 + pdi.3) / (pdi.2 + pdi.4) + plot_layout(guides = 'collect')

stool.nasal.ASV.diversity <- data.frame(evenness(physeq_HMP),d1)
skin.oral.ASV.diversity <- data.frame(evenness(physeq_UBM),d2)

stool.nasal.Genus.diversity <- data.frame(evenness(physeq_HMP_genus),d3)
skin.oral.Genus.diversity <- data.frame(evenness(physeq_UBM_genus),d4)

save(stool.nasal.ASV.diversity,stool.nasal.Genus.diversity,skin.oral.ASV.diversity,skin.oral.Genus.diversity, file = "./Robject/Diversity_Datatable.RData")


