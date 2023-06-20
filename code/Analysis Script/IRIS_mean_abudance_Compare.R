#Find difference of IRIS in individuals mean Microbiome data
#https://www.bioconductor.org/packages/devel/bioc/vignettes/microbiomeMarker/inst/doc/microbiomeMarker-vignette.html#installation
library(phyloseq)
library(ggplot2)
library(dplyr)
library(lefser)
library(mia)
library(microbial)
library(microbiomeMarker)
library(stringr)
library(ggstatsplot)

setwd("~/Library/CloudStorage/Box-Box/XinZhouFiles/Projects/ZXP2_HMP2_Core/")

source("./Analysis/Analysis Script/sxt.tools.R")

load("./Analysis/Robject/Revision_MultiOmes_0509.RData")
load("./Analysis/Robject/DetailedPhyloseq.RData")
sample.meta <- read.csv("./Analysis/ori_meta_table/metadata.subject.csv", header = T)
rownames(sample.meta) <- sample.meta$SubjectID

subjectlist.sk <- sample_data(physeq_SK_freq) %>% data.frame()
list.sk <- subjectlist.sk %>% select(SubjectID) %>% unique()
list.sk

sc86 <- sc[sc$SubjectID %in% list.sk$SubjectID,]
sd(sc86$Adj.age)

physeq_ST <- prune_taxa(taxa_sums(physeq_ST) > 0, physeq_ST) 
physeq_SK <- prune_taxa(taxa_sums(physeq_SK) > 0, physeq_SK) 
physeq_OR <- prune_taxa(taxa_sums(physeq_OR) > 0, physeq_OR) 
physeq_NS <- prune_taxa(taxa_sums(physeq_NS) > 0, physeq_NS) 

rowSums(otu_table(physeq_ST)) %>% mean()
rowSums(otu_table(physeq_ST)) %>% sd()
rowSums(otu_table(physeq_SK)) %>% mean()
rowSums(otu_table(physeq_SK)) %>% sd()
rowSums(otu_table(physeq_OR)) %>% mean()
rowSums(otu_table(physeq_OR)) %>% sd()
rowSums(otu_table(physeq_NS)) %>% mean()
rowSums(otu_table(physeq_NS)) %>% sd()

st.count <- phyloseq::otu_table(physeq_ST) %>% data.frame()
sk.count <- phyloseq::otu_table(physeq_SK) %>% data.frame()
or.count <- phyloseq::otu_table(physeq_OR) %>% data.frame()
ns.count <- phyloseq::otu_table(physeq_NS) %>% data.frame()

st.tax <- tax_table(physeq_ST) %>% data.frame()
sk.tax <- tax_table(physeq_SK) %>% data.frame()
or.tax <- tax_table(physeq_OR) %>% data.frame()
ns.tax <- tax_table(physeq_NS) %>% data.frame()

st.sample <- sample_data(physeq_ST) %>% data.frame()
sk.sample <- sample_data(physeq_SK) %>% data.frame()
or.sample <- sample_data(physeq_OR) %>% data.frame()
ns.sample <- sample_data(physeq_NS) %>% data.frame()

st.count.mean <- aggregate(st.count, by= list(st.sample$SubjectID), mean)
sk.count.mean <- aggregate(sk.count, by= list(sk.sample$SubjectID), mean)
or.count.mean <- aggregate(or.count, by= list(or.sample$SubjectID), mean)
ns.count.mean <- aggregate(ns.count, by= list(ns.sample$SubjectID), mean)

row.names(st.count.mean) <- st.count.mean$Group.1
row.names(sk.count.mean) <- sk.count.mean$Group.1
row.names(or.count.mean) <- or.count.mean$Group.1
row.names(ns.count.mean) <- ns.count.mean$Group.1

st.count.mean <- select(st.count.mean, -Group.1)
sk.count.mean <- select(sk.count.mean, -Group.1)
or.count.mean <- select(or.count.mean, -Group.1)
ns.count.mean <- select(ns.count.mean, -Group.1)

sample.meta[sample.meta$SubjectID %in% rownames(st.count.mean),-1]


physeq.st.subject <- phyloseq(otu_table(as.matrix(st.count.mean),taxa_are_rows=F), 
                              tax_table(as.matrix(st.tax)), 
                              sample_data(sample.meta[sample.meta$SubjectID %in% rownames(st.count.mean),-1])) %>% subset_taxa(Kingdom == "Bacteria")

physeq.sk.subject <- phyloseq(otu_table(as.matrix(sk.count.mean),taxa_are_rows=F),
                              tax_table(as.matrix(sk.tax)),
                              sample_data(sample.meta[sample.meta$SubjectID %in% rownames(sk.count.mean),-1])) %>% subset_taxa(Kingdom == "Bacteria")

physeq.or.subject <- phyloseq(otu_table(as.matrix(or.count.mean),taxa_are_rows=F), 
                              tax_table(as.matrix(or.tax)),
                              sample_data(sample.meta[sample.meta$SubjectID %in% rownames(or.count.mean),-1])) %>% subset_taxa(Kingdom == "Bacteria")

physeq.ns.subject <- phyloseq(otu_table(as.matrix(ns.count.mean),taxa_are_rows=F),
                              tax_table(as.matrix(ns.tax)),
                              sample_data(sample.meta[sample.meta$SubjectID %in% rownames(ns.count.mean),-1])) %>% subset_taxa(Kingdom == "Bacteria")



st.lefse <- prune_samples(physeq.st.subject@sam_data$IRIS != "Unknown", physeq.st.subject) %>%
  run_lefse(
    wilcoxon_cutoff = 0.05,
    norm = "TSS",
    group = "IRIS",
    kw_cutoff = 0.05,
    lda_cutoff = 0)
marker_table(st.lefse)
plot_abundance(st.lefse, group = "IRIS")

sk.lefse <- prune_samples(physeq.sk.subject@sam_data$IRIS != "Unknown", physeq.sk.subject) %>%
  run_lefse(
    wilcoxon_cutoff = 0.05,
    norm = "TSS",
    group = "IRIS",
    kw_cutoff = 0.05,
    lda_cutoff = 0)

plot_abundance(sk.lefse, group = "IRIS")

or.lefse <- prune_samples(physeq.or.subject@sam_data$IRIS != "Unknown", physeq.or.subject) %>%
  run_lefse(
    wilcoxon_cutoff = 0.05,
    norm = "TSS",
    group = "IRIS",
    kw_cutoff = 0.05,
    lda_cutoff = 0)

plot_abundance(or.lefse, group = "IRIS")

ns.lefse <- prune_samples(physeq.ns.subject@sam_data$IRIS != "Unknown", physeq.ns.subject) %>%
  run_lefse(
    wilcoxon_cutoff = 0.05,
    norm = "TSS",
    group = "IRIS",
    kw_cutoff = 0.05,
    lda_cutoff = 0)

plot_abundance(ns.lefse, group = "IRIS")

plot_cladogram(st.lefse, color = c(IS = "#BC3C29FF", IR = "#0072B5FF")) 

st.marker <- marker_table(st.lefse)%>% data.frame()
st.marker$bodysite <- "Stool"

sk.marker <- marker_table(sk.lefse)%>% data.frame()
sk.marker$bodysite <- "Skin"

or.marker <- marker_table(or.lefse)%>% data.frame()
or.marker$bodysite <- "Oral"

ns.marker <- marker_table(ns.lefse)%>% data.frame()
ns.marker$bodysite <- "Nasal"

gghistostats(sk.marker, x=ef_lda,normal.curve= T)

markers.lefse <- rbind(st.marker, sk.marker, or.marker, ns.marker)
#write.csv(file = "../../../human_microbiome_project/Supplementary_data/LEFSE_IRIS.result.csv",markers.lefse)

markers.lefse$bodysite <- factor(markers.lefse$bodysite, levels = c("Stool", "Skin", "Oral","Nasal"))

p.ldacompare <- ggbetweenstats(filter(markers.lefse),
                               x=bodysite, 
                               y=ef_lda,
                               type = "parametric",
                               p.adjust.method = "BH",
                               title = "Effect Size Comparasion")
p.ldacompare <- p.ldacompare + scale_color_manual(values = body_site_color)
p.ldacompare
#ggsave(filename = "../../../human_microbiome_project/Figures/FigureS1/EffectsizeCompare.pdf", width = 6, height = 5, dpi = 300)

physeq.st.subject_TSS <- transform_sample_counts(physeq.st.subject, function(x) x / sum(x))
physeq.sk.subject_TSS <- transform_sample_counts(physeq.sk.subject, function(x) x / sum(x)) 
physeq.or.subject_TSS <- transform_sample_counts(physeq.or.subject, function(x) x / sum(x)) 
physeq.ns.subject_TSS <- transform_sample_counts(physeq.ns.subject, function(x) x / sum(x)) 

rowSums(otu_table(physeq.st.subject_TSS))

physeq.st_Genus <- tax_glom(physeq.st.subject_TSS, taxrank = "Genus")
physeq.sk_Genus <- tax_glom(physeq.sk.subject_TSS, taxrank = "Genus")
physeq.or_Genus <- tax_glom(physeq.or.subject_TSS, taxrank = "Genus")
physeq.ns_Genus <- tax_glom(physeq.ns.subject_TSS, taxrank = "Genus")

genus.by.subject.st <- otu_table(physeq.st_Genus) %>% data.frame()
genus.by.subject.sk <- otu_table(physeq.sk_Genus) %>% data.frame()
genus.by.subject.or <- otu_table(physeq.or_Genus) %>% data.frame()
genus.by.subject.ns <- otu_table(physeq.ns_Genus) %>% data.frame()

st.tax$ASV <- rownames(st.tax)
sk.tax$ASV <- rownames(sk.tax)
or.tax$ASV <- rownames(or.tax)
ns.tax$ASV <- rownames(ns.tax)

colnames(genus.by.subject.st) <-  st.tax$Genus[match(colnames(genus.by.subject.st),st.tax$ASV)]
colnames(genus.by.subject.sk) <-  sk.tax$Genus[match(colnames(genus.by.subject.sk),sk.tax$ASV)]
colnames(genus.by.subject.or) <-  or.tax$Genus[match(colnames(genus.by.subject.or),or.tax$ASV)]
colnames(genus.by.subject.ns) <-  ns.tax$Genus[match(colnames(genus.by.subject.ns),ns.tax$ASV)]

stool.bygenus.bysubject <- merge(sample.meta, genus.by.subject.st, by="row.names")
skin.bygenus.bysubject <- merge(sample.meta, genus.by.subject.sk, by="row.names")
oral.bygenus.bysubject <- merge(sample.meta, genus.by.subject.or, by="row.names")
nasal.bygenus.bysubject <- merge(sample.meta, genus.by.subject.ns, by="row.names")

save(file = "./Analysis/Genus Table/by.subject.RData",stool.bygenus.bysubject, skin.bygenus.bysubject,oral.bygenus.bysubject,nasal.bygenus.bysubject)

colnames(stool.bygenus.bysubject)[str_detect(colnames(stool.bygenus.bysubject),"Akkermansia")]

load(file = "./Analysis/Genus Table/by.subject.RData")

p.st <- filter(stool.bygenus.bysubject, IRIS != "Unknown") %>%
  ggbetweenstats(x=IRIS,
               y=Akkermansia, 
               type = "nonparametric") 
p.st

colnames(skin.bygenus.bysubject)[str_detect(colnames(skin.bygenus.bysubject),"Cory")]

p.sk <- filter(skin.bygenus.bysubject, IRIS != "Unknown") %>%
  ggbetweenstats(x=IRIS,
                 y=Cutibacterium, 
                 type = "nonparametric") 
p.sk

colnames(oral.bygenus.bysubject)[str_detect(colnames(oral.bygenus.bysubject),"Pre")]
p.or <- filter(oral.bygenus.bysubject, IRIS != "Unknown") %>%
  ggbetweenstats(x=IRIS,
                 y=Prevotella, 
                 type = "nonparametric") 
p.or

colnames(nasal.bygenus.bysubject)[str_detect(colnames(nasal.bygenus.bysubject),"Cory")]
p.ns <- filter(nasal.bygenus.bysubject, IRIS != "Unknown") %>%
  ggbetweenstats(x=IRIS,
                 y=Corynebacterium, 
                 type = "nonparametric") 
p.ns


data.frame.meancyto <- select(ck.df, IL17F:ENA78, SubjectID) %>% group_by(SubjectID) %>% mutate_at(vars(-SubjectID),mean) %>% unique()
skin_sample_data <- sample_data(sample.meta[sample.meta$SubjectID %in% rownames(sk.count.mean),-1])

data.frame.meancyto.meta <- left_join(data.frame(skin_sample_data), as.data.frame(data.frame.meancyto), by="SubjectID")
dim(data.frame.meancyto.meta)

shapiro.test(data.frame.meancyto.meta$LEPTIN)
shapiro.test(data.frame.meancyto.meta$GMCSF)

pBMI_GMCSF <- ggscatterstats(data.frame.meancyto.meta, BMI, GMCSF, type = "nonparametric")
pBMI_LEPTIN <- ggscatterstats(data.frame.meancyto.meta, BMI, LEPTIN,type = "nonparametric")
p.corre <- pBMI_GMCSF + pBMI_LEPTIN
ggsave(filename = "./Analysis/Suppl.figure/LEPTIN.GMCSF.BMI.pdf",p.corre, width = 14, height = 7, dpi=300)




