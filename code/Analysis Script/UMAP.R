#Seurat UMAP of Microbiome

library(Seurat)
library(phyloseq)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(ggsci)
library(picante)
library(ggbiplot)
library(vegan)
library(digest)
library(stringr)

body_site_color = c(
  "Stool" = ggsci::pal_jama()(n=7)[2],
  "Skin" = ggsci::pal_jama()(n=7)[3],
  "Oral" = ggsci::pal_jama()(n=7)[4],
  "Nasal" = ggsci::pal_jama()(n=7)[5])

body_site_color2 = c(
  "ST" = ggsci::pal_jama()(n=7)[2],
  "Skin" = ggsci::pal_jama()(n=7)[3],
  "Oral" = ggsci::pal_jama()(n=7)[4],
  "NS" = ggsci::pal_jama()(n=7)[5])

base_theme = theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 13),
        panel.grid.minor = element_blank())

#https://hihayk.github.io/scale/#2/7/27/47/57/52/18/38/FFFF00/68/57/46/white
umap_color_scale <- scale_color_gradientn(trans="log10",colors = c("#C0C0C0","#FFEDCF","#FFE411","orange","#FF0000","#BA0000","#8F0000"),
                                          values = scales::rescale(log10(c(0,10,
                                                                           10.000001, 100,
                                                                           100.000001, 1000,
                                                                           1000.0000001, 1e+04,
                                                                           10000.000001,1e+05,
                                                                           100000.000001, 5e+05,
                                                                           500000.000001,1e+06))))

setwd("/Users/xzhou7/Library/CloudStorage/Box-Box/XinZhouFiles/Projects/ZXP2_HMP2_Core/Analysis/SeuratUMAP/")

load("../Robject/DetailedPhyloseq.RData")
meta_Subject_HMP <- read.csv("../ori_meta_table/metadata.subject.csv", header=T)
meta_Sample_HMP <- read.csv("../ori_meta_table/hmp2_infection_metadata_0323.csv", header=T)

load("../Robject/Revision_MultiOmes_0509.RData")

physeq_ST_Gen <- tax_glom(physeq_ST, taxrank = "Genus")
physeq_SK_Gen <- tax_glom(physeq_SK, taxrank = "Genus")
physeq_OR_Gen <- tax_glom(physeq_OR, taxrank = "Genus")
physeq_NS_Gen <- tax_glom(physeq_NS, taxrank = "Genus")

OTU_ST <- data.frame(otu_table(physeq_ST_Gen))
OTU_SK <- data.frame(otu_table(physeq_SK_Gen))
OTU_OR <- data.frame(otu_table(physeq_OR_Gen))
OTU_NS <- data.frame(otu_table(physeq_NS_Gen))

OTU_ST <- OTU_ST[colSums(OTU_ST)!=0]
OTU_SK <- OTU_SK[colSums(OTU_SK)!=0]
OTU_OR <- OTU_OR[colSums(OTU_OR)!=0]
OTU_NS <- OTU_NS[colSums(OTU_NS)!=0]

TAX1 <- data.frame(tax_table(physeq_ST))
TAX1$ASV <- row.names(TAX1)
TAX2 <- data.frame(tax_table(physeq_SK))
TAX2$ASV <- row.names(TAX2)

tOTU_ST <- data.frame(t(OTU_ST))
tOTU_ST$Genus <- TAX1$Genus[match(colnames(OTU_ST), TAX1$ASV)]
tOTU_ST <- tOTU_ST %>% select(Genus, everything())

tOTU_NS <- data.frame(t(OTU_NS))
tOTU_NS$Genus <- TAX1$Genus[match(colnames(OTU_NS), TAX1$ASV)]
tOTU_NS <- tOTU_NS %>% select(Genus, everything())

tOTU_NS[row.names(tOTU_NS) == "ASV8491",1] <- "Unclassified_C_Actinobacteria"

tOTU_SK <- data.frame(t(OTU_SK))
tOTU_SK$Genus <- TAX2$Genus[match(colnames(OTU_SK), TAX2$ASV)]
tOTU_SK <- tOTU_SK %>% select(Genus, everything())
tOTU_SK$Genus[duplicated(tOTU_SK$Genus)]
filter(tOTU_SK, Genus== "Unclassified_Actinobacteria")[1:2, 1:10]
tOTU_SK[row.names(tOTU_SK) == "OTU_1156",1] <- "Unclassified_C_Actinobacteria"

tOTU_OR <- data.frame(t(OTU_OR))
tOTU_OR$Genus <- TAX2$Genus[match(colnames(OTU_OR), TAX2$ASV)]
tOTU_OR <- tOTU_OR %>% select(Genus, everything())

HMP2.1 <- merge(tOTU_ST, tOTU_NS, by="Genus", all=T)
HMP2.2 <- merge(tOTU_SK, tOTU_OR, by="Genus", all=T)
HMP2 <- merge(HMP2.1,HMP2.2, by="Genus",all=T)

HMP2[is.na(HMP2)] <- 0
plot(density(colSums(HMP2[-1])))

thmp <- data.frame(t(HMP2))
colnames(thmp) <- thmp[1,]
thmp <- thmp[-1,]

HMP_count <- sapply(thmp, as.numeric, simplify=T)
row.names(HMP_count) <- row.names(thmp)
rowSums(HMP_count)

#make sure nothing is messed up with the transpose, this step should return T
identical(colSums(HMP2[-1]),rowSums(HMP_count))

st.sample<- data.frame(sample_data(physeq_ST))
sk.sample<- data.frame(sample_data(physeq_SK))
or.sample<- data.frame(sample_data(physeq_OR))
ns.sample<- data.frame(sample_data(physeq_NS))

dim(st.sample)
dim(ns.sample)
hmp.sample2.1 <- rbind(st.sample, ns.sample)
colnames(hmp.sample2.1)
hmp.sample2.1[1:3, 1:8]
hmp.sample2.1_short <- select(hmp.sample2.1, RandomID, SampleID, SubjectID, SampleType,Date,batch)

dim(or.sample)
dim(sk.sample)
hmp.sample2.2 <- rbind(or.sample, sk.sample)
colnames(hmp.sample2.2)
hmp.sample2.2[1:3,1:6]
hmp.sample2.2_short <- select(hmp.sample2.2, KitID, SampleID, SubjectID, SampleType, Date, batch)
colnames(hmp.sample2.2_short)[1] <- "RandomID"

hmp.sample <- rbind(hmp.sample2.1_short,hmp.sample2.2_short)

HMP_count_seuratobject <- data.frame(t(HMP_count))
##############################################################################################################################
# create seurat object    
##############################################################################################################################
HMP_seurat <- CreateSeuratObject(HMP_count_seuratobject, project = "iHMP", min.cells = 1, min.features = 1)
HMP_seurat <- AddMetaData(HMP_seurat, hmp.sample)
HMP_seurat[["percent.unclassified"]] <- PercentageFeatureSet(HMP_seurat, pattern = "^Unclassified")

VlnPlot(HMP_seurat,features = c("nFeature_RNA", "nCount_RNA", "percent.unclassified"), group.by = "SampleType", pt.size = 0)

FeatureScatter(HMP_seurat,feature1 = "nFeature_RNA", feature2 =  "percent.unclassified", group.by = "SampleType", pt.size = 0.5)
FeatureScatter(HMP_seurat,feature1 = "nCount_RNA", feature2 =  "percent.unclassified", group.by = "SampleType", pt.size = 0.5)

#transform to relative count per million
HMP_seurat <- NormalizeData(HMP_seurat, normalization.method = "RC", scale.factor = 1000000)

HMP_seurat <- FindVariableFeatures(HMP_seurat, selection.method = "vst", nfeatures = 2000)
HMP_seurat
VariableFeaturePlot(object = HMP_seurat)
head(VariableFeatures(HMP_seurat), 10)

all.genes <- rownames(HMP_seurat)
HMP_seurat <- ScaleData(HMP_seurat, features = all.genes)

HMP_seurat <- RunPCA(HMP_seurat, features = VariableFeatures(object = HMP_seurat))
DimPlot(HMP_seurat, reduction = "pca", group.by = "SampleType")
ElbowPlot(HMP_seurat)
HMP_seurat <- RunUMAP(HMP_seurat, dims = 1:9)
DimPlot(HMP_seurat, reduction = "umap", group.by = "SampleType",label = T)
DimPlot(HMP_seurat, reduction = "umap", group.by = "SubjectID",label = F)

HMP_seurat <- FindNeighbors(HMP_seurat, dims = 1:10)
HMP_seurat <- FindClusters(HMP_seurat, resolution = 0.2)
p.umap.ra <- DimPlot(HMP_seurat, group.by = "SampleType")
p.umap.ra <- p.umap.ra +coord_equal() + scale_color_manual(values = body_site_color)
p.umap.ra

HMP_seurat$SampleType[HMP_seurat$SampleType == "NS"] <- "Nasal"
HMP_seurat$SampleType[HMP_seurat$SampleType == "ST"] <- "Stool"


FeaturePlot(HMP_seurat, features="Prevotella", pt.size = 0.1)

#allMarkers <- FindAllMarkers(HMP_seurat,test.use = "bimod", min.pct = 0.25, only.pos = T)
#allMarkers %>% filter(cluster==3)
#allMarkers

#write.csv(file = "../SeuratUMAP/All.Markers.csv",allMarkers)

VlnPlot(HMP_seurat, features = c("Lawsonella", "percent.unclassified"))
VlnPlot(HMP_seurat, features = c("Lautropia", "percent.unclassified"))

#HMP_seurat <- RunTSNE(object = HMP_seurat)
#DimPlot(HMP_seurat, reduction = "tsne", group.by = "SampleType")
#DimPlot(HMP_seurat, reduction = "tsne", group.by = "SubjectID")
#DimPlot(HMP_seurat, reduction = "tsne", group.by = "batch", split.by = "SampleType")

######################################################################################################
# Create a phyloseq object
######################################################################################################
TAX <- rbind(TAX1,TAX2)[1:6]
TAX <- unique(TAX)
TAX[row.names(TAX) == "ASV8023",6] <- "Unclassified_C_Actinobacteria"
rownames(TAX) <- TAX$Genus
dim(TAX)
colnames(TAX)
dim(HMP_count_seuratobject)
HMP_count_seuratobject[1:10,1:10]
dim(hmp.sample)

hmp.sample$IRIS <- "unknown"
hmp.sample$IRIS[which(hmp.sample$SubjectID %in% filter(meta_Subject_HMP, IRIS=="IR")$SubjectID)] <- "IR"
hmp.sample$IRIS[which(hmp.sample$SubjectID %in% filter(meta_Subject_HMP, IRIS=="IS")$SubjectID)] <- "IS"

hmp.sample$FPG_Class <- "Normal"
hmp.sample$FPG_Class[which(hmp.sample$SubjectID %in% filter(meta_Subject_HMP, FPG_class=="Diabetes")$SubjectID)] <- "Diabetes"
hmp.sample$FPG_Class[which(hmp.sample$SubjectID %in% filter(meta_Subject_HMP, FPG_class=="Prediabetes")$SubjectID)] <- "Prediabetes"

hmp.sample$OTGG_Class <- "Unknown"
hmp.sample$OTGG_Class[which(hmp.sample$SubjectID %in% filter(meta_Subject_HMP, OGTT_Class=="Normal")$SubjectID)] <- "Normal"
hmp.sample$OTGG_Class[which(hmp.sample$SubjectID %in% filter(meta_Subject_HMP, OGTT_Class=="Prediabetes")$SubjectID)] <- "Prediabetes"
hmp.sample$OTGG_Class[which(hmp.sample$SubjectID %in% filter(meta_Subject_HMP, OGTT_Class=="Diabetes")$SubjectID)] <- "Diabetes"

hmp.sample$A1C_Class <- "1.Normal"
hmp.sample$A1C_Class[which(hmp.sample$SubjectID %in% filter(meta_Subject_HMP, A1C_Class=="2.NP")$SubjectID)] <- "2.NP"
hmp.sample$A1C_Class[which(hmp.sample$SubjectID %in% filter(meta_Subject_HMP, A1C_Class=="3.PN")$SubjectID)] <- "3.PN"
hmp.sample$A1C_Class[which(hmp.sample$SubjectID %in% filter(meta_Subject_HMP, A1C_Class=="4.VNP")$SubjectID)] <- "4.VNP"
hmp.sample$A1C_Class[which(hmp.sample$SubjectID %in% filter(meta_Subject_HMP, A1C_Class=="5.P")$SubjectID)] <- "5.P"
hmp.sample$A1C_Class[which(hmp.sample$SubjectID %in% filter(meta_Subject_HMP, A1C_Class=="6.VDP")$SubjectID)] <- "6.VDP"

hmp.sample$A1C_Class_S <- "1.Normal"
hmp.sample$A1C_Class_S[which(hmp.sample$SubjectID %in% filter(meta_Subject_HMP, A1C_Class=="2.NP")$SubjectID)] <- "2.NP"
hmp.sample$A1C_Class_S[which(hmp.sample$SubjectID %in% filter(meta_Subject_HMP, A1C_Class=="3.PN")$SubjectID)] <- "3.PN"
hmp.sample$A1C_Class_S[which(hmp.sample$SubjectID %in% filter(meta_Subject_HMP, A1C_Class=="4.VNP")$SubjectID)] <- "5.P"
hmp.sample$A1C_Class_S[which(hmp.sample$SubjectID %in% filter(meta_Subject_HMP, A1C_Class=="5.P")$SubjectID)] <- "5.P"
hmp.sample$A1C_Class_S[which(hmp.sample$SubjectID %in% filter(meta_Subject_HMP, A1C_Class=="6.VDP")$SubjectID)] <- "5.P"

hmp.sample$Th17_Class <- "1" 
hmp.sample$Th17_Class[which(hmp.sample$SubjectID %in% filter(meta_Subject_HMP, Th17.Group=="2")$SubjectID)] <- "2"
hmp.sample$Th17_Class[which(hmp.sample$SubjectID %in% filter(meta_Subject_HMP, Th17.Group=="3")$SubjectID)] <- "3"

OTU <- otu_table(HMP_count_seuratobject, taxa_are_rows = T)
TAX0 <- tax_table(as.matrix(TAX))
physeq.Combined = phyloseq(OTU, TAX0,sample_data(hmp.sample))
physeq.Combined_freq <- transform_sample_counts(physeq.Combined, function(x) x / sum(x))

#plot_bar(physeq.Combined_freq, fill = "Phylum")

physeq_C_ord <- ordinate(physeq.Combined_freq, "PCoA", dist = "bray")
physeq_C_ord$vectors[3155:3160, 1:5]

dim(physeq_C_ord$vectors)

sample_data(physeq.Combined_freq)$Prevotella <- as.numeric(otu_table(physeq.Combined_freq)["Prevotella"])

p1.1 <- plot_ordination(physeq.Combined_freq, physeq_C_ord, type = "samples",color = "Prevotella",  axes = c(1,2))  
p1.1 <- p1.1 + geom_point(size = 0.5) + theme_cowplot()
p1.1

#save(physeq_C_ord, file = "../Robject/Ordination_Bray.RData")

########################################################################################################################################
#replace PCA in seurat using BC distance matrix
########################################################################################################################################
head(Embeddings(HMP_seurat, reduction = "pca")[, 1:5])
physeq_C_ord$vectors[1:5, 1090:1094]

dim(physeq_C_ord$vectors)
hmp_bray <- physeq_C_ord$vectors

HMP_seurat[["bray"]] <- CreateDimReducObject(embeddings = hmp_bray, key = "Axis.", assay = DefaultAssay(HMP_seurat))
HMP_seurat[["bray"]]


DimPlot(HMP_seurat, reduction = "bray", pt.size = 0.5)
DimPlot(HMP_seurat)
VariableFeaturePlot(object = HMP_seurat)

HMP_seurat


HMP_seurat <- RunUMAP(HMP_seurat, reduction="bray", dims = 1:10)
pumap.combined <- DimPlot(HMP_seurat, pt.size = 0.05, group.by ="SampleType", label=T) + coord_equal() + scale_color_manual(values = body_site_color)
pumap.combined <- pumap.combined 
pumap.combined
#ggsave(filename ="../Suppl.figure/UMAP_Comb.pdf", pumap.combined, height = 3, width = 4, dpi=300)

UMAP.coordi <- HMP_seurat[["umap"]]@cell.embeddings
plot(UMAP.coordi[,1],UMAP.coordi[,2])

rownames(HMP_seurat)[str_detect(rownames(HMP_seurat), "Rum")]

featurelist <- c("Corynebacterium","Staphylococcus","Cutibacterium","Bacteroides",
                 "Unclassified-Ruminococcaceae","Prevotella","Streptococcus","Veillonella",
                 "Haemophilus","Neisseria", "Moraxella","Unclassified-Clostridia")

featurelist2 <- c("Neisseria","Unclassified-Neisseriales")#,"Lawsonella","Bacillus")

featurelist_all <- row.names(HMP_seurat)

RNAseurat <- GetAssayData(HMP_seurat) %>% data.frame()
featurelist_2ormore <- featurelist_all[(rowSums(RNAseurat!=0) > 1)]

# for (i in 1:length(featurelist_2ormore)){
#   print(featurelist_2ormore[i])
#   taxa <- featurelist_2ormore[i]
#   p <- FeaturePlot(HMP_seurat, features=taxa, pt.size = 0.1, max.cutoff = 1e+06, order = T) & base_theme & theme(legend.position = "none") & umap_color_scale
#   substr(taxa, 1,1)
#   taxa.clean <- str_replace(taxa,"/","_")
#   path <- paste0("./graph_result/", substr(taxa, 1,1),"/",taxa.clean,".pdf")
#   ggsave2(filename = path, p, width = 4.2, height = 4, dpi=300)
# }

FeaturePlot(HMP_seurat, features=featurelist, pt.size = 0.1, max.cutoff = 1e+06) & base_theme & theme(legend.position = "none") & umap_color_scale & coord_equal()

FeaturePlot(HMP_seurat, features="Moraxella", pt.size = 0.1, max.cutoff = 1e+05) & base_theme &  umap_color_scale & coord_equal()

########################################################check individual
table(HMP_seurat$SubjectID)

#Deidentification
HMP_seurat$individual <- "Backgound"
HMP_seurat$individual[HMP_seurat$SubjectID %in% "69-090"] <- "Participant #4"

p.individual <- DimPlot(HMP_seurat, group.by = "individual", order = T) & coord_equal() & scale_color_manual(values = c("light grey", "red")) 
p.individual <- p.individual + ggtitle("Participant #4")
p.individual

########################################################
head(AverageExpression(HMP_seurat, group.by = "SampleType"))

RidgePlot(HMP_seurat,"Acinetobacter", group.by = "orig.ident", log=F)
VlnPlot(HMP_seurat,"Akkermansia",group.by = "SampleType", log=T)

#save(HMP_seurat, hmp.sample, file = "./Seurat.Object.RData")
#write.csv(file = "../ori_meta_table/UMAP.coordi.csv", UMAP.coordi)


HMP_seurat$IRIS <- "unknown"
HMP_seurat$IRIS[which(HMP_seurat$SubjectID %in% filter(meta_Subject_HMP, IRIS=="IS")$SubjectID)] <- "IS"
HMP_seurat$IRIS[which(HMP_seurat$SubjectID %in% filter(meta_Subject_HMP, IRIS=="IR")$SubjectID)] <- "IR"

iris_color = 
  c("IR" = ggsci::pal_nejm()(n=8)[2],
    "IS" = ggsci::pal_nejm()(n=8)[1],
    "unknown"="grey") 

UMAP.IRIS <- DimPlot(HMP_seurat, group.by = "IRIS",pt.size = 0.1) & scale_color_manual(values = iris_color) & base_theme
UMAP.IRIS

#ggsave(filename = "../Suppl.figure/UMAP_IRIS.pdf",UMAP.IRIS, width = 4.2, height = 3, dpi=300)

########################################################################################################################################
#RDA analysis
########################################################################################################################################
physeq_C_ST_freq <- subset_samples(physeq.Combined_freq, SampleType=="ST")
physeq_C_ST_freq_ord <- ordinate(physeq_C_ST_freq, "PCoA", dist="bray")

p1.2 <- plot_ordination(physeq_C_ST_freq, physeq_C_ST_freq_ord, type = "samples",color = "Prevotella",  axes = c(1,2))  
p1.2 <- p1.2 + geom_point(size = 0.5) + theme_cowplot()
p1.2

physeq_C_NS_freq <- subset_samples(physeq.Combined_freq, SampleType=="NS")
physeq_C_NS_freq_ord <- ordinate(physeq_C_NS_freq, "PCoA", dist="bray")

p1.3 <- plot_ordination(physeq_C_NS_freq, physeq_C_NS_freq_ord, type = "samples", axes = c(1,2))
p1.3 <- p1.3 + geom_point(size = 0.5) + theme_cowplot()
p1.3

physeq_C_OR_freq <- subset_samples(physeq.Combined_freq, SampleType=="Oral")
physeq_C_OR_freq_ord <- ordinate(physeq_C_OR_freq, "PCoA", dist="bray")

p1.4 <- plot_ordination(physeq_C_OR_freq, physeq_C_OR_freq_ord, type = "samples", axes = c(1,2))
p1.4 <- p1.4 + geom_point(size = 0.5) + theme_cowplot()
p1.4

physeq_C_SK_freq <- subset_samples(physeq.Combined_freq, SampleType=="Skin")
physeq_C_SK_freq_ord <- ordinate(physeq_C_SK_freq, "PCoA", dist="bray")

p1.5 <- plot_ordination(physeq_C_SK_freq, physeq_C_SK_freq_ord, type = "samples", axes = c(1,2))
p1.5 <- p1.5 + geom_point(size = 0.5) + theme_cowplot()
p1.5


###########################<Stool Microbiome>#############################################################################################################
RDA_ST = ordinate(physeq_C_ST_freq, "PCoA", "bray")
summary(RDA_ST)
# p2 = plot_ordination(physeq_C_ST_freq, RDA_ST, 
#                      type="samples", 
#                      color="IRIS", 
#                      #label="Well",
#                      # type="taxa",color="Phylum", 
#                      title="RDA based on weighted BC distance")+
#   geom_point(size=7, alpha=0.75)+ 
#   theme_bw()
# p2

#https://r-from-a-learners-perspective.readthedocs.io/en/latest/part5/
vegan_otu <- function(physeq) {
  OTU <- otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

physeq_C_ST_freq_filter <- filter_taxa(physeq_C_ST_freq, function(x) mean(x) > 0, TRUE)
otutable_ST <- vegan_otu(physeq_C_ST_freq_filter)

decorana(otutable_ST)

hmp.sample.st <- subset(hmp.sample, hmp.sample$RandomID %in% row.names(otutable_ST))

#CCA
ord_ST_all <- cca(otutable_ST ~ IRIS + FPG_Class + OTGG_Class + A1C_Class_S + Th17_Class, data=hmp.sample.st) 
ord_ST_all

CCA_ord_ST<- cca(otutable_ST ~ ., hmp.sample.st[,c("IRIS", "FPG_Class", "OTGG_Class","A1C_Class_S","Th17_Class")]) 

CCA_ord_ST

vif.cca(ord_ST_all)

fit1 <- envfit(CCA_ord_ST,hmp.sample.st[,c("IRIS", "OTGG_Class","A1C_Class","Th17_Class")], perm = 999, display = "lc", scaling = "sites")
plot(CCA_ord_ST)
plot(CCA_ord_ST, type="n")

plot(fit1, col = "red", cex=1.2, axis=TRUE, p.max = 0.05)

summary(CCA_ord_ST)
anova(CCA_ord_ST)

scrs1 <- as.data.frame(scores(CCA_ord_ST, display = "sites"))
scrs1

Bray_Ord_ST <- ordinate(physeq_C_ST_freq_filter, "PCoA", dist = "bray")
Bray_Ord_ST[1:5]
Bray_Ord_ST$vectors
##########################################
##########################################
#https://zhuanlan.zhihu.com/p/493961702
# 
# enviroPCA <- rda(bray_ST,  scale = TRUE, na.action = "na.omit")
# 
# eig <- enviroPCA$CA$eig
# percent_var <- eig * 100 / sum(eig)
# PC1_varEx <- round(percent_var[1], 1)
# PC2_varEx <- round(percent_var[2], 1)
# 
# contrib <- round(100*scores(enviroPCA, display = "sp", scaling = 0)[,2]^2, 3)
# 
# library(tidyverse)
# enviroPCA_vars  <- scores(enviroPCA, display = 'species', choices = c(1, 2), scaling = 2) %>% 
#   as.data.frame() %>% 
#   rownames_to_column(., var = 'var') %>% 
#   mutate(contrib = contrib)
# 
# enviroPCA_vars
# 
# enviroPCA_vars_1 <- filter(enviroPCA_vars, contrib > 0.1)
# 
# filter(enviroPCA_vars, var == "Prevotella")
# 
# library(wesanderson)
# pal <- wes_palette("Darjeeling1", 3, type = "continuous")
# 
# library(ggplot2)
# ng1 <- theme(aspect.ratio=0.7,panel.background = element_blank(),
#              panel.grid.major = element_blank(),
#              panel.grid.minor = element_blank(),
#              panel.border=element_blank(),
#              axis.line.x = element_line(color="black",size=1),
#              axis.line.y = element_line(color="black",size=1),
#              axis.ticks=element_line(size = 1, color="black"),
#              axis.ticks.length=unit(0.25, 'cm'),
#              axis.text=element_text(color="black",size=15),
#              axis.title=element_text(color="black",size=1),
#              axis.title.y=element_text(vjust=2,size=17),
#              axis.title.x=element_text(vjust=0.1,size=17),
#              axis.text.x=element_text(size=15),
#              axis.text.y=element_text(size=15),
#              strip.text.x = element_text(size = 10, colour = "black",face = "bold"),
#              strip.background = element_rect(colour="black"),
#              legend.position = "top", legend.direction="vertical",
#              legend.text=element_text(size=17), legend.key = element_rect(fill = "white"),
#              legend.title = element_text(size=17),legend.key.size = unit(1.0, "cm"))
# 
# enviroPCA_variableContrib <- ggplot() +
#   geom_hline(yintercept = 0, linetype = "dotted") +
#   geom_vline(xintercept = 0, linetype = "dotted") +
#   geom_segment(data = enviroPCA_vars_1, aes(x = 0, xend = PC1, y=0, yend = PC2, color = contrib), 
#                size = 2, arrow = arrow(length = unit(0.02, "npc")), alpha = 1) +
#   geom_text(data = enviroPCA_vars_1,
#             aes(x = PC1, y = PC2, label = var,
#                 hjust = "inward", vjust =  0.5 * (1 - sign(PC1))),
#             color = "black", size = 3.5) + 
#   xlab(sprintf("PC1 (%.1f%%)", PC1_varEx)) + ylab(sprintf("PC2 (%.1f%%)", PC2_varEx)) +
#   scale_colour_gradientn(colours = rev(pal), breaks = seq(from = 5, to = 25, by = 5)) +
#   # scale_x_continuous(breaks = seq(from = -1, to = 1, by = 0.25)) +
#   # scale_y_continuous(breaks = seq(from = -1, to = 1, by = 0.25)) +
#   ng1 + theme(legend.position = "top",
#               legend.direction="horizontal",
#               # legend.title = element_blank(),
#               legend.key.size = unit(0.5, "cm"),
#               legend.spacing.x = unit(0.1, "cm"),
#               legend.text = element_text(size=10)) +
#   guides(color = guide_colourbar(barwidth = 10, barheight = 0.5))
# 
# enviroPCA_variableContrib

##########################################
#bdRDA
##########################################
bray_ST <- vegdist(otutable_ST, method = "bray")
#RDA_ST_all <- capscale(bray_ST ~ IRIS + FPG_Class + OTGG_Class + A1C_Class + Th17_Class, hmp.sample.st, add = TRUE)
RDA_ST_all <- capscale(bray_ST ~ IRIS, hmp.sample.st, add = T)

plot(RDA_ST_all, type = "p", scaling = "sites") 
RsquareAdj(RDA_ST_all)

anova(RDA_ST_all, by="axis", perm.max=500)
anova(RDA_ST_all, by="terms", permu=200)

temp <- vif.cca(RDA_ST_all)
temp
select_para <- names(temp[temp < 10])
select_para

fit <- envfit(RDA_ST_all, hmp.sample.st[,c("IRIS","A1C_Class", "Th17_Class")], perm = 999, display = "lc", scaling = "sites")
fit$factors$pvals

plot(RDA_ST_all)
plot(fit, col = "red", cex=1.2, axis=TRUE, p.max = 0.05)

summary(RDA_ST_all)
anova(RDA_ST_all)

scrs <- as.data.frame(scores(RDA_ST_all, display = "sites"))

spp.scrs.lib <- data.frame(fit$factors$centroids)
#spp.scrs <- filter(spp.scrs.lib, ! rownames(spp.scrs.lib) %in% c("OTGG_ClassUnknown","OTGG_ClassNormal", "OTGG_ClassPrediabetes","OTGG_ClassDiabetes", "IRISunknown"))
spp.scrs <- filter(spp.scrs.lib, rownames(spp.scrs.lib) %in% c("A1C_Class1.Normal", "A1C_Class2.NP","A1C_Class3.PN","A1C_Class4.VNP","A1C_Class5.P","A1C_Class6.VDP"))
#spp.scrs <- filter(spp.scrs.lib, rownames(spp.scrs.lib) %in% c("Th17_Class1", "Th17_Class2","Th17_Class3"))
#spp.scrs <- filter(spp.scrs.lib, rownames(spp.scrs.lib) %in% c("IRISIR", "IRISIS"))
#spp.scrs <- filter(spp.scrs.lib, rownames(spp.scrs.lib) %in% c("OTGG_ClassNormal", "OTGG_ClassPrediabetes","OTGG_ClassDiabetes"))

#spp.scrs <- filter(spp.scrs.lib, rownames(spp.scrs.lib) %in% c("A1C_Class_S1.Normal", "A1C_Class_S2.NP","A1C_Class_S3.PN","A1C_Class_S5.P"))

p.st <- ggplot(scrs) + geom_point(mapping = aes(x = CAP1, y = CAP2), alpha = 0.4, size = 1, color="#DF8F44FF") +
  geom_segment(data = spp.scrs, aes(x = 0, xend = CAP1*25, y = 0, yend = CAP2*25),arrow = arrow(length = unit(0.25, "cm")), colour = "blue", size=1) +
  geom_label_repel(data = spp.scrs, aes(x = CAP1*25, y = CAP2*25, label = row.names(spp.scrs)), size = 5, color="black") + base_theme +
  theme(axis.title=element_text(size=14), axis.text=element_text(size=14), legend.text=element_text(size=14)) +
  scale_color_manual(values =c("red", "blue", "black"), guide = guide_legend(ncol=1)) + theme(legend.position="right") 
p.st  <- p.st + base_theme
p.st

#ggsave(filename = "../Suppl.figure/RDA_A1C_ST.pdf",p.st, width = 5, height = 5, dpi=300)


spp.scrs <- filter(spp.scrs.lib, rownames(spp.scrs.lib) %in% c("Th17_Class1", "Th17_Class2","Th17_Class3"))

p.st <- ggplot(scrs) + geom_point(mapping = aes(x = CAP1, y = CAP2), alpha = 0.4, size = 1, color="#DF8F44FF") +
  geom_segment(data = spp.scrs, aes(x = 0, xend = CAP1*25, y = 0, yend = CAP2*25),arrow = arrow(length = unit(0.25, "cm")), colour = "blue", size=1) +
  geom_label_repel(data = spp.scrs, aes(x = CAP1*25, y = CAP2*25, label = row.names(spp.scrs)), size = 5, color="black") + base_theme +
  theme(axis.title=element_text(size=14), axis.text=element_text(size=14), legend.text=element_text(size=14)) +
  scale_color_manual(values =c("red", "blue", "black"), guide = guide_legend(ncol=1)) + theme(legend.position="right") 
p.st  <- p.st + base_theme
p.st

ggsave(filename = "../Suppl.figure/RDA_Th17_ST.pdf",p.st, width = 5, height = 5, dpi=300)

#identify the most abundant genera in Stool sample
sort(colSums(otutable_ST),decreasing = T)[1:10]

p.Bacteroides.st <- cbind(scrs,otutable_ST) %>% ggplot(aes(x=CAP1, y=CAP2, color=Bacteroides)) + geom_point(size=0.3) + coord_equal() +
  theme_cowplot() + scale_colour_gradientn(values=c(1, .3, .2, .1, 0), colours=c("#770000","red", "orange", "grey", "grey")) + base_theme + theme(legend.position="top")
p.Bacteroides.st

p.Prevotella.st <- cbind(scrs,otutable_ST) %>% ggplot(aes(x=CAP1, y=CAP2, color=Prevotella)) + geom_point(size=0.3) + coord_equal() +
  theme_cowplot() + scale_colour_gradientn(values=c(1, .3, .2, .1, 0), colours=c("#770000","red", "orange", "grey", "grey")) + base_theme+ theme(legend.position="top")
p.Prevotella.st

p.Phocaeicola.st <- cbind(scrs,otutable_ST) %>%ggplot(aes(x=CAP1, y=CAP2, color=Phocaeicola)) + geom_point(size=0.3) + coord_equal() +
  theme_cowplot() + scale_colour_gradientn(values=c(1, .3, .2, .1, 0), colours=c("#770000","red", "orange", "grey", "grey")) + base_theme+ theme(legend.position="top")
p.Phocaeicola.st

p.Faecalibacterium.st <- cbind(scrs,otutable_ST) %>%ggplot(aes(x=CAP1, y=CAP2, color=Faecalibacterium)) + geom_point(size=0.3) + coord_equal() +
  theme_cowplot() + scale_colour_gradientn(values=c(1, .3, .2, .1, 0), colours=c("#770000","red", "orange", "grey", "grey")) + base_theme+ theme(legend.position="top")
p.Faecalibacterium.st

p.Unclassified_Lachnospiraceae.st <- cbind(scrs,otutable_ST) %>%ggplot(aes(x=CAP1, y=CAP2, color=Unclassified_Lachnospiraceae)) + geom_point(size=0.3)  + coord_equal() +
  theme_cowplot() + scale_colour_gradientn(values=c(1, .3, .2, .1, 0), colours=c("#770000","red", "orange", "grey", "grey")) + base_theme + theme(legend.position="top")
p.Unclassified_Lachnospiraceae.st

p.Unclassified_Ruminococcaceae.st <- cbind(scrs,otutable_ST) %>%ggplot(aes(x=CAP1, y=CAP2, color=Unclassified_Ruminococcaceae)) + geom_point(size=0.3) + coord_equal() +
  theme_cowplot() + scale_colour_gradientn(values=c(1, .3, .2, .1, 0), colours=c("#770000","red", "orange", "grey", "grey")) + base_theme + theme(legend.position="top")
p.Unclassified_Ruminococcaceae.st

p.Alistipes.st <- cbind(scrs,otutable_ST) %>% ggplot(aes(x=CAP1, y=CAP2, color=Alistipes)) + geom_point(size=0.3) + coord_equal() +
  theme_cowplot() + scale_colour_gradientn(values=c(1, .3, .2, .1, 0), colours=c("#770000","red", "orange", "grey", "grey")) + base_theme + theme(legend.position="top")
p.Alistipes.st

p.Blautia.st <- cbind(scrs,otutable_ST) %>% ggplot(aes(x=CAP1, y=CAP2, color=Blautia)) + geom_point(size=0.3) + coord_equal() +
  theme_cowplot() + scale_colour_gradientn(values=c(1, .3, .2, .1, 0), colours=c("#770000","red", "orange", "grey", "grey")) + base_theme + theme(legend.position="top")
p.Blautia.st

p.Oscillibacter.st <- cbind(scrs,otutable_ST) %>% ggplot(aes(x=CAP1, y=CAP2, color=Oscillibacter)) + geom_point(size=0.3) + coord_equal() +
  theme_cowplot() + scale_colour_gradientn(values=c(1, .3, .2, .1, 0), colours=c("#770000","red", "orange", "grey", "grey")) + base_theme + theme(legend.position="top")
p.Oscillibacter.st

p.Agathobacter.st <- cbind(scrs,otutable_ST) %>% ggplot(aes(x=CAP1, y=CAP2, color=Agathobacter)) + geom_point(size=0.3) + coord_equal() +
  theme_cowplot() + scale_colour_gradientn(values=c(1, .3, .2, .1, 0), colours=c("#770000","red", "orange", "grey", "grey")) + base_theme + theme(legend.position="top")
p.Agathobacter.st

p.Bacteroides.st/p.Prevotella.st/p.Phocaeicola.st

p.st.suppli.braycurtis <- p.Bacteroides.st + p.Prevotella.st + p.Phocaeicola.st + p.Unclassified_Ruminococcaceae.st
p.st.suppli.braycurtis
#ggsave2(filename = "../Suppl.figure/BC.rep.genus.Stool.pdf", height = 5, width = 9, dpi = 300)

ggsave(filename = "../Suppl.figure/PCOA/Bacteroides.st.pdf",p.Bacteroides.st, height = 5, width = 4, dpi = 300)
ggsave(filename = "../Suppl.figure/PCOA/Prevotella.st.st.pdf",p.Prevotella.st, height = 5, width = 4, dpi = 300)
ggsave(filename = "../Suppl.figure/PCOA/Phocaeicola.st.pdf",p.Phocaeicola.st, height = 5, width = 4, dpi = 300)
ggsave(filename = "../Suppl.figure/PCOA/Faecalibacterium.st.pdf",p.Faecalibacterium.st, height = 5, width = 4, dpi = 300)
ggsave(filename = "../Suppl.figure/PCOA/Unclassified_Lachnospiraceae.st.pdf",p.Unclassified_Lachnospiraceae.st, height = 5, width = 4, dpi = 300)
ggsave(filename = "../Suppl.figure/PCOA/Unclassified_Ruminococcaceae.st.pdf",p.Unclassified_Ruminococcaceae.st, height = 5, width = 4, dpi = 300)
ggsave(filename = "../Suppl.figure/PCOA/Alistipes.st.pdf",p.Alistipes.st, height = 5, width = 4, dpi = 300)
ggsave(filename = "../Suppl.figure/PCOA/Blautia.st.pdf",p.Blautia.st, height = 5, width = 4, dpi = 300)
ggsave(filename = "../Suppl.figure/PCOA/Oscillibacter.st.pdf",p.Oscillibacter.st, height = 5, width = 4, dpi = 300)
ggsave(filename = "../Suppl.figure/PCOA/Agathobacter.st.pdf",p.Agathobacter.st, height = 5, width = 4, dpi = 300)

##########################<Skin Microbiome>########################################################################
physeq_C_SK_freq_filter <- filter_taxa(physeq_C_SK_freq, function(x) mean(x) > 0, TRUE)
otutable_SK <- vegan_otu(physeq_C_SK_freq_filter)
table(colSums(otutable_SK) > 0)

decorana(otutable_SK)

#remove UBMkid222819843 because it equals 0 after trim
#sum(subset(otutable_SK, rownames(otutable_SK) %in% "UBMkid222819843"))

#otutable_SK <- subset(otutable_SK, rownames(otutable_SK) != "UBMkid222819843")

#get sample metadata
hmp.sample.SK <- subset(hmp.sample, hmp.sample$RandomID %in% row.names(otutable_SK))

bray_SK <- vegdist(otutable_SK, method = "bray")
#RDA_SK_all <- capscale(bray_SK ~ IRIS + FPG_Class + OTGG_Class + A1C_Class + Th17_Class, hmp.sample.SK)
RDA_SK_all <- capscale(bray_SK ~ IRIS,hmp.sample.SK)
plot(RDA_SK_all, type = "p", scaling = "sites") 

temp2 <- vif.cca(RDA_SK_all)
temp2
select_para2 <- names(temp2[temp2 < 10])
select_para2

fit.sk <- envfit(RDA_SK_all, hmp.sample.SK[,c("IRIS", "FPG_Class","OTGG_Class","A1C_Class", "Th17_Class")], perm = 999, display = "lc", scaling = "sites")
fit.sk$factors$pvals

plot(fit.sk, col = "red", cex=1.2, axis=TRUE, p.max = 0.05)

summary(RDA_SK_all)
anova(RDA_SK_all)

scrs.sk <- as.data.frame(scores(RDA_SK_all, display = "sites"))

scores(RDA_SK_all, display = "cn")
scrs.sk

spp.scrs.lib.sk <- data.frame(fit.sk$factors$centroids)
spp.scrs.lib.sk
# spp.scrs.sk <- filter(spp.scrs.lib.sk, ! rownames(spp.scrs.lib.sk) %in% c("OTGG_ClassUnknown", "IRISunknown"))
spp.scrs.sk <- filter(spp.scrs.lib.sk, rownames(spp.scrs.lib.sk) %in% c("A1C_Class1.Normal", "A1C_Class2.NP","A1C_Class3.PN","A1C_Class4.VNP","A1C_Class5.P","A1C_Class6.VDP"))
# spp.scrs.sk <- filter(spp.scrs.lib.sk, rownames(spp.scrs.lib.sk) %in% c("A1C_Class_S1.Normal", "A1C_Class_S2.NP","A1C_Class_S3.PN","A1C_Class_S5.P"))
# spp.scrs.sk <- filter(spp.scrs.lib.sk, rownames(spp.scrs.lib.sk) %in% c("Th17_Class1", "Th17_Class2","Th17_Class3"))
# spp.scrs.sk <- filter(spp.scrs.lib.sk, rownames(spp.scrs.lib.sk) %in% c("IRISIR", "IRISIS","IRISunknown"))
# spp.scrs.sk <- filter(spp.scrs.lib.sk, rownames(spp.scrs.lib.sk) %in% c("OTGG_ClassNormal", "OTGG_ClassPrediabetes","OTGG_ClassDiabetes"))
# spp.scrs.sk <- filter(spp.scrs.lib.sk, rownames(spp.scrs.lib.sk) %in% c("FPG_ClassNormal", "FPG_ClassPrediabetes","FPG_ClassDiabetes"))

spp.scrs.sk

p.sk <- ggplot(scrs.sk) + geom_point(mapping = aes(x = CAP1, y = CAP2), alpha = 0.4, size = 1, color="#00A1D5FF") +
  geom_segment(data = spp.scrs.sk, aes(x = 0, xend = CAP1*10, y = 0, yend = CAP2*10),arrow = arrow(length = unit(0.25, "cm")), colour = "blue", size=1) +
  geom_label_repel(data = spp.scrs.sk, aes(x = CAP1*10, y = CAP2*10, label = row.names(spp.scrs.sk)), size = 5, color="black") + theme_bw() +
  theme(axis.title=element_text(size=14), axis.text=element_text(size=14), legend.text=element_text(size=14)) +
  scale_color_manual(values =c("red", "blue", "black"), guide = guide_legend(ncol=1)) + theme(legend.position="right") 

p.sk + base_theme
p.sk

ggsave(filename = "../Suppl.figure/RDA_A1C_SK.pdf",p.sk, width = 5, height = 5, dpi=300)

#Th17
spp.scrs.sk <- filter(spp.scrs.lib.sk, rownames(spp.scrs.lib.sk) %in% c("Th17_Class1", "Th17_Class2","Th17_Class3"))
spp.scrs.sk

p.sk <- ggplot(scrs.sk) + geom_point(mapping = aes(x = CAP1, y = CAP2), alpha = 0.4, size = 1, color="#00A1D5FF") +
  geom_segment(data = spp.scrs.sk, aes(x = 0, xend = CAP1*10, y = 0, yend = CAP2*10),arrow = arrow(length = unit(0.25, "cm")), colour = "blue", size=1) +
  geom_label_repel(data = spp.scrs.sk, aes(x = CAP1*10, y = CAP2*10, label = row.names(spp.scrs.sk)), size = 5, color="black") + theme_bw() +
  theme(axis.title=element_text(size=14), axis.text=element_text(size=14), legend.text=element_text(size=14)) +
  scale_color_manual(values =c("red", "blue", "black"), guide = guide_legend(ncol=1)) + theme(legend.position="right") 

p.sk + base_theme
p.sk

ggsave(filename = "../Suppl.figure/RDA_Th17_SK.pdf",p.sk, width = 5, height = 5, dpi=300)

#FPG
spp.scrs.sk <- filter(spp.scrs.lib.sk, rownames(spp.scrs.lib.sk) %in% c("FPG_ClassNormal", "FPG_ClassPrediabetes","FPG_ClassDiabetes"))
spp.scrs.sk

p.sk <- ggplot(scrs.sk) + geom_point(mapping = aes(x = CAP1, y = CAP2), alpha = 0.4, size = 1, color="#00A1D5FF") +
  geom_segment(data = spp.scrs.sk, aes(x = 0, xend = CAP1*10, y = 0, yend = CAP2*10),arrow = arrow(length = unit(0.25, "cm")), colour = "blue", size=1) +
  geom_label_repel(data = spp.scrs.sk, aes(x = CAP1*10, y = CAP2*10, label = row.names(spp.scrs.sk)), size = 5, color="black") + theme_bw() +
  theme(axis.title=element_text(size=14), axis.text=element_text(size=14), legend.text=element_text(size=14)) +
  scale_color_manual(values =c("red", "blue", "black"), guide = guide_legend(ncol=1)) + theme(legend.position="right") 

p.sk + base_theme
p.sk

ggsave(filename = "../Suppl.figure/RDA_FPG_SK.pdf",p.sk, width = 5, height = 5, dpi=300)
rownames(otutable_SK)
sort(colSums(otutable_SK),decreasing = T)[1:10]
sort(colSums(otutable_SK!=0),decreasing = T)[1:20]

colSums(otutable_SK!=0)

#Plot 10 most abundant bacteria
p.Corynebacterium.sk <- cbind(scrs.sk,otutable_SK) %>%ggplot(aes(x=CAP1, y=CAP2, color=Corynebacterium)) + geom_point(size=0.3)  + coord_equal() +
  theme_cowplot() + scale_colour_gradientn(values=c(1, .3, .2, .1, 0), colours=c("#770000","red", "orange", "grey", "grey")) + base_theme + theme(legend.position="top")
p.Corynebacterium.sk  

p.Staphylococcus.sk <- cbind(scrs.sk,otutable_SK) %>%ggplot(aes(x=CAP1, y=CAP2, color=Staphylococcus)) + geom_point(size=0.3)  + coord_equal() +
  theme_cowplot() + scale_colour_gradientn(values=c(1, .3, .2, .1, 0), colours=c("#770000","red", "orange", "grey", "grey")) + base_theme + theme(legend.position="top")
p.Staphylococcus.sk  

p.Cutibacterium.sk <- cbind(scrs.sk,otutable_SK) %>%ggplot(aes(x=CAP1, y=CAP2, color=Cutibacterium)) + geom_point(size=0.3)  + coord_equal() +
  theme_cowplot() + scale_colour_gradientn(values=c(1, .5, .4, .2, 0), colours=c("#770000","red", "orange", "grey", "grey")) + base_theme + theme(legend.position="top")
p.Cutibacterium.sk  

p.Streptococcus.sk <- cbind(scrs.sk,otutable_SK) %>%ggplot(aes(x=CAP1, y=CAP2, color=Streptococcus)) + geom_point(size=0.3) + coord_equal() +
  theme_cowplot() + scale_colour_gradientn(values=c(1, .3, .2, .1, 0), colours=c("#770000","red", "orange", "grey", "grey")) + base_theme + theme(legend.position="top")
p.Streptococcus.sk

p.Anaerococcus.sk <- cbind(scrs.sk,otutable_SK) %>%ggplot(aes(x=CAP1, y=CAP2, color=Anaerococcus)) + geom_point(size=0.3)  + coord_equal() +
  theme_cowplot() + scale_colour_gradientn(values=c(1, .3, .2, .1, 0), colours=c("#770000","red", "orange", "grey", "grey")) + base_theme + theme(legend.position="top")
p.Anaerococcus.sk 

p.Delftia.sk <- cbind(scrs.sk,otutable_SK) %>%ggplot(aes(x=CAP1, y=CAP2, color=Delftia)) + geom_point(size=0.3) + coord_equal() +
  theme_cowplot() + scale_colour_gradientn(values=c(1, .3, .2, .1, 0), colours=c("#770000","red", "orange", "grey", "grey")) + base_theme + theme(legend.position="top")
p.Delftia.sk 

p.Lysinibacillus.sk <- cbind(scrs.sk,otutable_SK) %>%ggplot(aes(x=CAP1, y=CAP2, color=Lysinibacillus)) + geom_point(size=0.3)  + coord_equal() +
  theme_cowplot() + scale_colour_gradientn(values=c(1, .3, .2, .1, 0), colours=c("#770000","red", "orange", "grey", "grey")) + base_theme + theme(legend.position="top")
p.Lysinibacillus.sk 

p.Peptoniphilus.sk <- cbind(scrs.sk,otutable_SK) %>%ggplot(aes(x=CAP1, y=CAP2, color=Peptoniphilus)) + geom_point(size=0.3)  + coord_equal() +
  theme_cowplot() + scale_colour_gradientn(values=c(1, .3, .2, .1, 0), colours=c("#770000","red", "orange", "grey", "grey")) + base_theme + theme(legend.position="top")
p.Peptoniphilus.sk 

p.Prevotella.sk <- cbind(scrs.sk,otutable_SK) %>%ggplot(aes(x=CAP1, y=CAP2, color=Prevotella)) + geom_point(size=0.3)  + coord_equal() +
  theme_cowplot() + scale_colour_gradientn(values=c(1, .3, .2, .1, 0), colours=c("#770000","red", "orange", "grey", "grey")) + base_theme + theme(legend.position="top")
p.Prevotella.sk 

p.Brucella.sk <- cbind(scrs.sk,otutable_SK) %>%ggplot(aes(x=CAP1, y=CAP2, color=Brucella)) + geom_point(size=0.3)  + coord_equal() +
  theme_cowplot() + scale_colour_gradientn(values=c(1, .3, .2, .1, 0), colours=c("#770000","red", "orange", "grey", "grey")) + base_theme + theme(legend.position="top")
p.Brucella.sk 

p.Moraxella.sk <- cbind(scrs.sk,otutable_SK) %>%ggplot(aes(x=CAP1, y=CAP2, color=Moraxella)) + geom_point(size=0.3)  + coord_equal() +
  theme_cowplot() + scale_colour_gradientn(values=c(1, .3, .2, .1, 0), colours=c("#770000","red", "orange", "grey", "grey")) + base_theme + theme(legend.position="top")
p.Moraxella.sk 

#save as supplementary figure
p.sk.suppli.braycurtis <- p.Cutibacterium.sk + p.Staphylococcus.sk + p.Corynebacterium.sk + p.Anaerococcus.sk
p.sk.suppli.braycurtis
ggsave2(filename = "../Suppl.figure/BC.rep.genus.Skin.pdf", height = 12, width = 8, dpi = 300)

ggsave(filename = "../Suppl.figure/PCOA/Corynebacterium.sk.pdf",p.Corynebacterium.sk, height = 7, width = 4, dpi = 300)
ggsave(filename = "../Suppl.figure/PCOA/Staphylococcus.sk.pdf",p.Staphylococcus.sk, height = 7, width = 4, dpi = 300)
ggsave(filename = "../Suppl.figure/PCOA/Cutibacterium.sk.pdf",p.Cutibacterium.sk, height = 7, width = 4, dpi = 300)
ggsave(filename = "../Suppl.figure/PCOA/Streptococcus.sk.pdf",p.Streptococcus.sk, height = 7, width = 4, dpi = 300)
ggsave(filename = "../Suppl.figure/PCOA/Anaerococcus.sk.pdf",p.Anaerococcus.sk, height = 7, width = 4, dpi = 300)
ggsave(filename = "../Suppl.figure/PCOA/Delftia.sk.pdf",p.Delftia.sk, height = 7, width = 4, dpi = 300)
ggsave(filename = "../Suppl.figure/PCOA/Peptoniphilus.sk.pdf",p.Peptoniphilus.sk, height = 7, width = 4, dpi = 300)
ggsave(filename = "../Suppl.figure/PCOA/Lysinibacillus.sk.pdf",p.Lysinibacillus.sk, height = 7, width = 4, dpi = 300)
ggsave(filename = "../Suppl.figure/PCOA/Prevotella.sk.pdf",p.Prevotella.sk, height = 7, width = 4, dpi = 300)
ggsave(filename = "../Suppl.figure/PCOA/Brucella.sk.pdf",p.Brucella.sk, height = 7, width = 4, dpi = 300)
ggsave(filename = "../Suppl.figure/PCOA/Moraxella.sk.pdf",p.Moraxella.sk, height = 7, width = 4, dpi = 300)

##############################<Oral Microbiome>####################################################################
physeq_C_OR_freq_filter <- filter_taxa(physeq_C_OR_freq, function(x) mean(x) > 0, TRUE)
otutable_OR <- vegan_otu(physeq_C_OR_freq_filter)
dim(otutable_OR)
#otutable_OR <- subset(otutable_OR, rownames(otutable_OR) != c("UBMkid429104718","UBMkid851105741")) 

decorana(otutable_OR)

#get sample metadata
hmp.sample.OR <- subset(hmp.sample, hmp.sample$RandomID %in% row.names(otutable_OR))

bray_OR <- vegdist(otutable_OR, method = "bray")
RDA_OR_all <- capscale(bray_OR ~ IRIS + FPG_Class + OTGG_Class + A1C_Class + Th17_Class, hmp.sample.OR)

plot(RDA_OR_all, type = "p", scaling = "sites")

temp3 <- vif.cca(RDA_OR_all)
temp3
select_para3 <- names(temp3[temp3 < 10])
select_para3

fit.or <- envfit(RDA_OR_all, hmp.sample.OR[,c("IRIS","FPG_Class", "OTGG_Class","A1C_Class","Th17_Class")], perm = 999, display = "lc", scaling = "sites")
fit.or$factors$pvals

plot(fit.or, col = "red", cex=1.2, axis=TRUE, p.max = 0.05)

summary(RDA_OR_all)
anova(RDA_OR_all)

scrs.or <- as.data.frame(scores(RDA_OR_all, display = "sites"))
scrs.or

spp.scrs.lib.or <- data.frame(fit.or$factors$centroids)
# spp.scrs.or <- filter(spp.scrs.lib.or, ! rownames(spp.scrs.lib.or) %in% c("OTGG_ClassUnknown", "IRISunknown"))
spp.scrs.or <- filter(spp.scrs.lib.or, rownames(spp.scrs.lib.or) %in% c("A1C_Class1.Normal", "A1C_Class2.NP","A1C_Class3.PN","A1C_Class4.VNP","A1C_Class5.P","A1C_Class6.VDP"))
# spp.scrs.or <- filter(spp.scrs.lib.or, rownames(spp.scrs.lib.or) %in% c("Th17_Class1", "Th17_Class2","Th17_Class3"))
# spp.scrs.or <- filter(spp.scrs.lib.or, rownames(spp.scrs.lib.or) %in% c("IRISIR", "IRISIS","IRISunknown"))
# spp.scrs.or <- filter(spp.scrs.lib.or, rownames(spp.scrs.lib.or) %in% c("OTGG_ClassNormal", "OTGG_ClassPrediabetes","OTGG_ClassDiabetes"))
 
spp.scrs.or

p.or <- ggplot(scrs.or) + geom_point(mapping = aes(x = CAP1, y = CAP2), alpha = 0.8, size = 1, color="#B24745FF") +
  geom_segment(data = spp.scrs.or, aes(x = 0, xend = CAP1*10, y = 0, yend = CAP2*10),arrow = arrow(length = unit(0.25, "cm")), colour = "blue", size=1) +
  geom_label_repel(data = spp.scrs.or, aes(x = CAP1*10, y = CAP2*10, label = row.names(spp.scrs.or)), size = 5, color="black") + theme_bw() +
  theme(axis.title=element_text(size=14), axis.text=element_text(size=14), legend.text=element_text(size=14)) +
  scale_color_manual(values =c("red", "blue", "black"), guide = guide_legend(ncol=1)) + theme(legend.position="right") 

p.or

#ggsave(filename = "../Suppl.figure/RDA_A1C_OR.pdf",p.or, width = 5, height = 5, dpi=300)

spp.scrs.or <- filter(spp.scrs.lib.or, rownames(spp.scrs.lib.or) %in% c("Th17_Class1", "Th17_Class2","Th17_Class3"))
spp.scrs.or

p.or <- ggplot(scrs.or) + geom_point(mapping = aes(x = CAP1, y = CAP2), alpha = 0.8, size = 1, color="#B24745FF") +
  geom_segment(data = spp.scrs.or, aes(x = 0, xend = CAP1*10, y = 0, yend = CAP2*10),arrow = arrow(length = unit(0.25, "cm")), colour = "blue", size=1) +
  geom_label_repel(data = spp.scrs.or, aes(x = CAP1*10, y = CAP2*10, label = row.names(spp.scrs.or)), size = 5, color="black") + theme_bw() +
  theme(axis.title=element_text(size=14), axis.text=element_text(size=14), legend.text=element_text(size=14)) +
  scale_color_manual(values =c("red", "blue", "black"), guide = guide_legend(ncol=1)) + theme(legend.position="right") 

p.or

#ggsave(filename = "../Suppl.figure/RDA_Th17_OR.pdf",p.or, width = 5, height = 5, dpi=300)

spp.scrs.or <- filter(spp.scrs.lib.or, rownames(spp.scrs.lib.or) %in% c("FPG_ClassNormal", "FPG_ClassPrediabetes","FPG_ClassDiabetes"))
spp.scrs.or

p.or <- ggplot(scrs.or) + geom_point(mapping = aes(x = CAP1, y = CAP2), alpha = 0.8, size = 1, color="#B24745FF") +
  geom_segment(data = spp.scrs.or, aes(x = 0, xend = CAP1*10, y = 0, yend = CAP2*10),arrow = arrow(length = unit(0.25, "cm")), colour = "blue", size=1) +
  geom_label_repel(data = spp.scrs.or, aes(x = CAP1*10, y = CAP2*10, label = row.names(spp.scrs.or)), size = 5, color="black") + theme_bw() +
  theme(axis.title=element_text(size=14), axis.text=element_text(size=14), legend.text=element_text(size=14)) +
  scale_color_manual(values =c("red", "blue", "black"), guide = guide_legend(ncol=1)) + theme(legend.position="right") 

p.or

#ggsave(filename = "../Suppl.figure/RDA_FPG_OR.pdf",p.or, width = 5, height = 5, dpi=300)

sort(colSums(otutable_OR)/1001,decreasing = T)[1:10]

p.Prevotella.or <- cbind(scrs.or,otutable_OR) %>%ggplot(aes(x=CAP1, y=CAP2, color=Prevotella)) + geom_point(size=0.3) + coord_equal() +
  theme_cowplot() + scale_colour_gradientn(values=c(1, .6, .4, .2, 0), colours=c("#770000","red", "orange", "grey", "grey"))+ base_theme + theme(legend.position="top")
p.Prevotella.or  

p.Neisseria.or <- cbind(scrs.or,otutable_OR) %>%ggplot(aes(x=CAP1, y=CAP2, color=Neisseria)) + geom_point(size=0.3) + coord_equal() +
  theme_cowplot() + scale_colour_gradientn(values=c(1, .3, .2, .1, 0), colours=c("#770000","red", "orange", "grey", "grey"))+ base_theme + theme(legend.position="top")
p.Neisseria.or  

p.Veillonella.or <- cbind(scrs.or,otutable_OR) %>%ggplot(aes(x=CAP1, y=CAP2, color=Veillonella)) + geom_point(size=0.3) + coord_equal() +
  theme_cowplot() + scale_colour_gradientn(values=c(1, .3, .2, .1, 0), colours=c("#770000","red", "orange", "grey", "grey"))+ base_theme + theme(legend.position="top")
p.Veillonella.or

p.Haemophilus.or <- cbind(scrs.or,otutable_OR) %>%ggplot(aes(x=CAP1, y=CAP2, color=Haemophilus)) + geom_point(size=0.3) + coord_equal() +
  theme_cowplot() + scale_colour_gradientn(values=c(1, .3, .2, .1, 0), colours=c("#770000","red", "orange", "grey", "grey"))+ base_theme + theme(legend.position="top")
p.Haemophilus.or  

p.Streptococcus.or <- cbind(scrs.or,otutable_OR) %>%ggplot(aes(x=CAP1, y=CAP2, color=Streptococcus)) + geom_point(size=0.3) + coord_equal() +
  theme_cowplot() + scale_colour_gradientn(values=c(1, .3, .2, .1, 0), colours=c("#770000","red", "orange", "grey", "grey"))+ base_theme + theme(legend.position="top")
p.Streptococcus.or  

p.Schaalia.or <- cbind(scrs.or,otutable_OR) %>%ggplot(aes(x=CAP1, y=CAP2, color=Schaalia)) + geom_point(size=0.3) + coord_equal() +
  theme_cowplot() + scale_colour_gradientn(values=c(1, .2, .1, .05, 0), colours=c("#770000","red", "orange", "grey", "grey"))+ base_theme + theme(legend.position="top")
p.Schaalia.or 

p.Leptotrichia.or <- cbind(scrs.or,otutable_OR) %>%ggplot(aes(x=CAP1, y=CAP2, color=Leptotrichia)) + geom_point(size=0.3) + coord_equal() +
  theme_cowplot() + scale_colour_gradientn(values=c(1, .3, .2, .1, 0), colours=c("#770000","red", "orange", "grey", "grey"))+ base_theme + theme(legend.position="top")
p.Leptotrichia.or 

p.Rothia.or <- cbind(scrs.or,otutable_OR) %>%ggplot(aes(x=CAP1, y=CAP2, color=Rothia)) + geom_point(size=0.3) + coord_equal() +
  theme_cowplot() + scale_colour_gradientn(values=c(1, .3, .2, .1, 0), colours=c("#770000","red", "orange", "grey", "grey"))+ base_theme + theme(legend.position="top")
p.Rothia.or 

p.Porphyromonas.or <- cbind(scrs.or,otutable_OR) %>%ggplot(aes(x=CAP1, y=CAP2, color=Porphyromonas)) + geom_point(size=0.3) + coord_equal() +
  theme_cowplot() + scale_colour_gradientn(values=c(1, .3, .2, .1, 0), colours=c("#770000","red", "orange", "grey", "grey"))+ base_theme + theme(legend.position="top")
p.Porphyromonas.or 

p.Campylobacter.or <- cbind(scrs.or,otutable_OR) %>%ggplot(aes(x=CAP1, y=CAP2, color=Campylobacter)) + geom_point(size=0.3) + coord_equal() +
  theme_cowplot() + scale_colour_gradientn(values=c(1, .3, .2, .1, 0), colours=c("#770000","red", "orange", "grey", "grey"))+ base_theme + theme(legend.position="top")
p.Campylobacter.or 

p.Alloprevotella.or <- cbind(scrs.or,otutable_OR) %>%ggplot(aes(x=CAP1, y=CAP2, color=Alloprevotella)) + geom_point(size=0.3) + coord_equal() +
  theme_cowplot() + scale_colour_gradientn(values=c(1, .3, .2, .1, 0), colours=c("#770000","red", "orange", "grey", "grey"))+ base_theme + theme(legend.position="top")
p.Alloprevotella.or 

p.Unclassified_Clostridia.or <- cbind(scrs.or,otutable_OR) %>%ggplot(aes(x=CAP1, y=CAP2, color=Unclassified_Clostridia)) + geom_point(size=0.3) + coord_equal() +
  theme_cowplot() + scale_colour_gradientn(values=c(1, .3, .2, .1, 0), colours=c("#770000","red", "orange", "grey", "grey"))+ base_theme + theme(legend.position="top")
p.Unclassified_Clostridia.or  

p.Unclassified.Selenomonadaceae.or <- cbind(scrs.or,otutable_OR) %>%ggplot(aes(x=CAP1, y=CAP2, color=Unclassified_Selenomonadaceae)) + geom_point(size=0.3) + coord_equal() +
  theme_cowplot() + scale_colour_gradientn(values=c(1, .3, .2, .1, 0), colours=c("#770000","red", "orange", "grey", "grey"))+ base_theme + theme(legend.position="top")
p.Unclassified.Selenomonadaceae.or 

p.or.suppli.braycurtis <- p.Prevotella.or  + p.Neisseria.or + p.Veillonella.or + p.Haemophilus.or + p.Streptococcus.or + p.Leptotrichia.or 
p.or.suppli.braycurtis
ggsave2(filename = "../Suppl.figure/BC.rep.genus.Oral.pdf", p.or.suppli.braycurtis,height = 5, width = 11, dpi = 300)
ggsave(filename = "../Suppl.figure/PCOA/Prevotella.or.pdf",p.Prevotella.or, height = 5.5, width = 4, dpi = 300)
ggsave(filename = "../Suppl.figure/PCOA/Neisseria.or.pdf",p.Neisseria.or, height = 5.5, width = 4, dpi = 300)
ggsave(filename = "../Suppl.figure/PCOA/Veillonella.or.pdf",p.Veillonella.or, height = 5.5, width = 4, dpi = 300)
ggsave(filename = "../Suppl.figure/PCOA/Haemophilus.or.pdf",p.Haemophilus.or, height = 5.5, width = 4, dpi = 300)
ggsave(filename = "../Suppl.figure/PCOA/Streptococcus.or.pdf",p.Streptococcus.or, height = 5.5, width = 4, dpi = 300)
ggsave(filename = "../Suppl.figure/PCOA/Schaalia.or.pdf",p.Schaalia.or, height = 5.5, width = 4, dpi = 300)
ggsave(filename = "../Suppl.figure/PCOA/Leptotrichia.or.pdf",p.Leptotrichia.or, height = 5.5, width = 4, dpi = 300)
ggsave(filename = "../Suppl.figure/PCOA/Rothia.or.pdf",p.Rothia.or, height = 5.5, width = 4, dpi = 300)
ggsave(filename = "../Suppl.figure/PCOA/Porphyromonas.or.pdf",p.Porphyromonas.or, height = 5.5, width = 4, dpi = 300)
ggsave(filename = "../Suppl.figure/PCOA/Campylobacter.or.pdf",p.Campylobacter.or, height = 5.5, width = 4, dpi = 300)
ggsave(filename = "../Suppl.figure/PCOA/Alloprevotella.or.pdf",p.Alloprevotella.or, height = 5.5, width = 4, dpi = 300)
ggsave(filename = "../Suppl.figure/PCOA/Unclassified_Clostridia.or.pdf",p.Unclassified_Clostridia.or, height = 5.5, width = 4, dpi = 300)
ggsave(filename = "../Suppl.figure/PCOA/Unclassified.Selenomonadaceae.or.pdf",p.Unclassified.Selenomonadaceae.or, height = 5.5, width = 4, dpi = 300)
ggsave(filename = "../Suppl.figure/PCOA/Unclassified.Selenomonadaceae.or.pdf",p.Unclassified.Selenomonadaceae.or, height = 5.5, width = 4, dpi = 300)

##############################<Nasal Microbiome>####################################################################
physeq_C_NS_freq_filter <- filter_taxa(physeq_C_NS_freq, function(x) mean(x) > 0, TRUE)
otutable_NS <- vegan_otu(physeq_C_NS_freq_filter)
dim(otutable_NS)

decorana(otutable_NS)

#get sample metadata
hmp.sample.NS <- subset(hmp.sample, hmp.sample$RandomID %in% row.names(otutable_NS))

bray_NS <- vegdist(otutable_NS, method = "bray")
RDA_NS_all <- capscale(bray_NS ~ IRIS + FPG_Class + OTGG_Class + A1C_Class + Th17_Class, hmp.sample.NS)

plot(RDA_NS_all, type = "p", scaling = "sites")

temp4 <- vif.cca(RDA_NS_all)
temp4
select_para4 <- names(temp4[temp4 < 10])
select_para4

fit.ns <- envfit(RDA_NS_all, hmp.sample.NS[,c("IRIS", "OTGG_Class","FPG_Class","A1C_Class","Th17_Class")], perm = 999, display = "lc", scaling = "sites")
fit.ns$factors$pvals

plot(fit.ns, col = "red", cex=1.2, axis=TRUE, p.max = 0.05)

summary(RDA_NS_all)
anova(RDA_NS_all)

scrs.ns <- as.data.frame(scores(RDA_NS_all, display = "sites"))
scrs.ns

spp.scrs.lib.ns <- data.frame(fit.ns$factors$centroids)
# spp.scrs.ns <- filter(spp.scrs.lib.ns, ! rownames(spp.scrs.lib.ns) %in% c("OTGG_ClassUnknown", "IRISunknown"))
spp.scrs.ns <- filter(spp.scrs.lib.ns, rownames(spp.scrs.lib.ns) %in% c("A1C_Class1.Normal", "A1C_Class2.NP","A1C_Class3.PN","A1C_Class4.VNP","A1C_Class5.P","A1C_Class6.VDP"))
# spp.scrs.ns <- filter(spp.scrs.lib.ns, rownames(spp.scrs.lib.ns) %in% c("Th17_Class1", "Th17_Class2","Th17_Class3"))
# spp.scrs.ns <- filter(spp.scrs.lib.ns, rownames(spp.scrs.lib.ns) %in% c("IRISIR", "IRISIS","IRISunknown"))
# spp.scrs.ns <- filter(spp.scrs.lib.ns, rownames(spp.scrs.lib.ns) %in% c("OTGG_ClassNormal", "OTGG_ClassPrediabetes","OTGG_ClassDiabetes"))
spp.scrs.ns

p.ns <- ggplot(scrs.ns) + geom_point(mapping = aes(x = CAP1, y = CAP2), alpha = 0.8, size = 1, color="#79AF97FF") +
  geom_segment(data = spp.scrs.ns, aes(x = 0, xend = CAP1*10, y = 0, yend = CAP2*10),arrow = arrow(length = unit(0.25, "cm")), colour = "blue", size=1) +
  geom_label_repel(data = spp.scrs.ns, aes(x = CAP1*10, y = CAP2*10, label = row.names(spp.scrs.ns)), size = 5, color="black", alpha=1) + theme_bw() +
  theme(axis.title=element_text(size=14), axis.text=element_text(size=14), legend.text=element_text(size=14)) +
  scale_color_manual(values =c("red", "blue", "black"), guide = guide_legend(ncol=1)) + theme(legend.position="right") 
p.ns

ggsave(filename = "../Suppl.figure/RDA_A1C_NS.pdf",p.ns, width = 5, height = 5, dpi=300)

spp.scrs.ns <- filter(spp.scrs.lib.ns, rownames(spp.scrs.lib.ns) %in% c("Th17_Class1", "Th17_Class2","Th17_Class3"))
spp.scrs.ns

p.ns <- ggplot(scrs.ns) + geom_point(mapping = aes(x = CAP1, y = CAP2), alpha = 0.8, size = 1, color="#79AF97FF") +
  geom_segment(data = spp.scrs.ns, aes(x = 0, xend = CAP1*10, y = 0, yend = CAP2*10),arrow = arrow(length = unit(0.25, "cm")), colour = "blue", size=1) +
  geom_label_repel(data = spp.scrs.ns, aes(x = CAP1*10, y = CAP2*10, label = row.names(spp.scrs.ns)), size = 5, color="black", alpha=1) + theme_bw() +
  theme(axis.title=element_text(size=14), axis.text=element_text(size=14), legend.text=element_text(size=14)) +
  scale_color_manual(values =c("red", "blue", "black"), guide = guide_legend(ncol=1)) + theme(legend.position="right") 
p.ns

spp.scrs.ns <- filter(spp.scrs.lib.ns, rownames(spp.scrs.lib.ns) %in% c("OTGG_ClassNormal", "OTGG_ClassPrediabetes","OTGG_ClassDiabetes"))
spp.scrs.ns
p.ns <- ggplot(scrs.ns) + geom_point(mapping = aes(x = CAP1, y = CAP2), alpha = 0.8, size = 1, color="#79AF97FF") +
  geom_segment(data = spp.scrs.ns, aes(x = 0, xend = CAP1*10, y = 0, yend = CAP2*10),arrow = arrow(length = unit(0.25, "cm")), colour = "blue", size=1) +
  geom_label_repel(data = spp.scrs.ns, aes(x = CAP1*10, y = CAP2*10, label = row.names(spp.scrs.ns)), size = 5, color="black", alpha=1) + theme_bw() +
  theme(axis.title=element_text(size=14), axis.text=element_text(size=14), legend.text=element_text(size=14)) +
  scale_color_manual(values =c("red", "blue", "black"), guide = guide_legend(ncol=1)) + theme(legend.position="right") 
p.ns

ggsave(filename = "../Suppl.figure/RDA_OGTT_NS.pdf",p.ns, width = 5, height = 5, dpi=300)

spp.scrs.ns <- filter(spp.scrs.lib.ns, rownames(spp.scrs.lib.ns) %in% c("FPG_ClassNormal", "FPG_ClassPrediabetes","FPG_ClassDiabetes"))
spp.scrs.ns
p.ns <- ggplot(scrs.ns) + geom_point(mapping = aes(x = CAP1, y = CAP2), alpha = 0.8, size = 1, color="#79AF97FF") +
  geom_segment(data = spp.scrs.ns, aes(x = 0, xend = CAP1*10, y = 0, yend = CAP2*10),arrow = arrow(length = unit(0.25, "cm")), colour = "blue", size=1) +
  geom_label_repel(data = spp.scrs.ns, aes(x = CAP1*10, y = CAP2*10, label = row.names(spp.scrs.ns)), size = 5, color="black", alpha=1) + theme_bw() +
  theme(axis.title=element_text(size=14), axis.text=element_text(size=14), legend.text=element_text(size=14)) +
  scale_color_manual(values =c("red", "blue", "black"), guide = guide_legend(ncol=1)) + theme(legend.position="right") 
p.ns

ggsave(filename = "../Suppl.figure/RDA_FPG_NS.pdf",p.ns, width = 5, height = 5, dpi=300)


sort(colSums(otutable_NS)/922,decreasing = T)[1:10]
sort(colSums(otutable_NS != 0),decreasing = T)[1:10]

p.Corynebacterium.ns <- cbind(scrs.ns,otutable_NS) %>%ggplot(aes(x=CAP1, y=CAP2, color=Corynebacterium)) + geom_point(size=0.3) + coord_equal() +
  theme_cowplot() + scale_colour_gradientn(values=c(1, .5, .2, .1, 0), colours=c("#770000","red", "orange", "grey", "grey"))+ base_theme+ theme(legend.position="top")
p.Corynebacterium.ns 

p.Staphylococcus.ns <- cbind(scrs.ns,otutable_NS) %>%ggplot(aes(x=CAP1, y=CAP2, color=Staphylococcus)) + geom_point(size=0.3) + coord_equal() +
  theme_cowplot() + scale_colour_gradientn(values=c(1, .5, .2, .1, 0), colours=c("#770000","red", "orange", "grey", "grey"))+ base_theme+ theme(legend.position="top")
p.Staphylococcus.ns 

p.Cutibacterium.ns <- cbind(scrs.ns,otutable_NS) %>%ggplot(aes(x=CAP1, y=CAP2, color=Cutibacterium)) + geom_point(size=0.3) + coord_equal() +
  theme_cowplot() + scale_colour_gradientn(values=c(1, .5, .2, .1, 0), colours=c("#770000","red", "orange", "grey", "grey"))+ base_theme+ theme(legend.position="top")
p.Cutibacterium.ns 

p.Dolosigranulum.ns <- cbind(scrs.ns,otutable_NS) %>%ggplot(aes(x=CAP1, y=CAP2, color=Dolosigranulum)) + geom_point(size=0.3) + coord_equal() +
  theme_cowplot() + scale_colour_gradientn(values=c(1, .5, .2, .1, 0), colours=c("#770000","red", "orange", "grey", "grey"))+ base_theme+ theme(legend.position="top")
p.Dolosigranulum.ns 

p.Unclassified_Betaproteobacteria.ns <- cbind(scrs.ns,otutable_NS) %>%ggplot(aes(x=CAP1, y=CAP2, color=Unclassified_Betaproteobacteria)) + geom_point(size=0.3) + coord_equal() +
  theme_cowplot() + scale_colour_gradientn(values=c(1, .5, .2, .1, 0), colours=c("#770000","red", "orange", "grey", "grey"))+ base_theme+ theme(legend.position="top")
p.Unclassified_Betaproteobacteria.ns 

p.Streptococcus.ns <- cbind(scrs.ns,otutable_NS) %>%ggplot(aes(x=CAP1, y=CAP2, color=Streptococcus)) + geom_point(size=0.3) + coord_equal() +
  theme_cowplot() + scale_colour_gradientn(values=c(1, .3, .2, .01, 0), colours=c("#770000","red", "orange", "grey", "grey"))+ base_theme+ theme(legend.position="top")
p.Streptococcus.ns 

p.Moraxella.ns <- cbind(scrs.ns,otutable_NS) %>%ggplot(aes(x=CAP1, y=CAP2, color=Moraxella)) + geom_point(size=0.3) + coord_equal() +
  theme_cowplot() + scale_colour_gradientn(values=c(1, .3, .2, .05, 0), colours=c("#770000","red", "orange", "grey", "grey"))+ base_theme+ theme(legend.position="top")
p.Moraxella.ns 

p.Anaerococcus.ns <- cbind(scrs.ns,otutable_NS) %>%ggplot(aes(x=CAP1, y=CAP2, color=Anaerococcus)) + geom_point(size=0.3) + coord_equal() +
  theme_cowplot() + scale_colour_gradientn(values=c(1, .3, .2, .05, 0), colours=c("#770000","red", "orange", "grey", "grey"))+ base_theme+ theme(legend.position="top")
p.Anaerococcus.ns 

p.Peptoniphilus.ns <- cbind(scrs.ns,otutable_NS) %>%ggplot(aes(x=CAP1, y=CAP2, color=Peptoniphilus)) + geom_point(size=0.3) + coord_equal() +
  theme_cowplot() + scale_colour_gradientn(values=c(1, .5, .2, .1, 0), colours=c("#770000","red", "orange", "grey", "grey"))+ base_theme+ theme(legend.position="top")
p.Peptoniphilus.ns 

p.Lawsonella.ns <- cbind(scrs.ns,otutable_NS) %>%ggplot(aes(x=CAP1, y=CAP2, color=Lawsonella)) + geom_point(size=0.3) + coord_equal() +
  theme_cowplot() + scale_colour_gradientn(values=c(1, .5, .2, .1, 0), colours=c("#770000","red", "orange", "grey", "grey"))+ base_theme+ theme(legend.position="top")
p.Lawsonella.ns 

p.ns.suppli.braycurtis  <- p.Corynebacterium.ns + p.Staphylococcus.ns + p.Cutibacterium.ns + p.Peptoniphilus.ns 
p.ns.suppli.braycurtis
ggsave2(filename = "../Suppl.figure/BC.rep.genus.Nasal.pdf", p.ns.suppli.braycurtis,height = 5, width = 8, dpi = 300)

pA1C <-  p.sk + p.ns + p.or + p.st 
ggsave(filename = "../Suppl.figure/A1Cgroup_BCdistance.pdf", pA1C, height = 10, width = 10, dpi=300)

ggsave(filename = "../Suppl.figure/PCOA/Corynebacterium.ns.pdf",p.Corynebacterium.ns, height = 6, width = 4, dpi = 300)
ggsave(filename = "../Suppl.figure/PCOA/Staphylococcus.ns.pdf",p.Staphylococcus.ns, height = 6, width = 4, dpi = 300)
ggsave(filename = "../Suppl.figure/PCOA/Cutibacterium.ns.pdf",p.Cutibacterium.ns, height = 6, width = 4, dpi = 300)
ggsave(filename = "../Suppl.figure/PCOA/Dolosigranulum.ns.pdf",p.Dolosigranulum.ns, height = 6, width = 4, dpi = 300)
ggsave(filename = "../Suppl.figure/PCOA/Unclassified_Betaproteobacteria.ns.pdf",p.Unclassified_Betaproteobacteria.ns, height = 6, width = 4, dpi = 300)
ggsave(filename = "../Suppl.figure/PCOA/Streptococcus.ns.pdf",p.Streptococcus.ns, height = 6, width = 4, dpi = 300)
ggsave(filename = "../Suppl.figure/PCOA/Moraxella.ns.pdf",p.Moraxella.ns, height = 6, width = 4, dpi = 300)
ggsave(filename = "../Suppl.figure/PCOA/Peptoniphilus.ns.pdf",p.Peptoniphilus.ns, height = 6, width = 4, dpi = 300)
ggsave(filename = "../Suppl.figure/PCOA/Lawsonella.ns.pdf",p.Lawsonella.ns, height = 6, width = 4, dpi = 300)
ggsave(filename = "../Suppl.figure/PCOA/Anaerococcus.ns.pdf",p.Anaerococcus.ns, height = 6, width = 4, dpi = 300)

####################Revision######
#stool
physeq_C_ST_freq_filter
physeq_C_ST_ord <- ordinate(physeq_C_ST_freq_filter,method = "PCoA", distance = "bray")

p1 <- plot_ordination(physeq_C_ST_freq_filter, physeq_C_ST_ord, type="sample", color="IRIS")
p1 <- p1 + ggtitle("The Stool Microbiome") +  theme_cowplot() + coord_equal() + base_theme + theme(legend.position="top") + scale_color_manual(values = iris_color)
p1
ggsave2("~/Desktop/PCoA/Stool.PCoA.pdf", p1, width = 5, height = 6, dpi = 300)

ST.pcoa.axis <- physeq_C_ST_ord$vectors[,1:2] %>% data.frame()

p.Bacteroides.st <- cbind(ST.pcoa.axis,otutable_ST) %>% ggplot(aes(x=Axis.1, y=Axis.2, color=Bacteroides)) + geom_point(size=0.3) + coord_equal() +
  theme_cowplot() + scale_colour_gradientn(values=c(1, .3, .2, .1, 0), colours=c("#770000","red", "orange", "grey", "grey")) + base_theme + theme(legend.position="top")
p.Bacteroides.st
ggsave2("~/Desktop/PCoA/p.Bacteroides.st.pdf", p.Bacteroides.st, width = 3, height = 3.5, dpi = 300)

p.Prevotella.st <- cbind(ST.pcoa.axis,otutable_ST) %>% ggplot(aes(x=Axis.1, y=Axis.2, color=Prevotella)) + geom_point(size=0.3) + coord_equal() +
  theme_cowplot() + scale_colour_gradientn(values=c(1, .3, .2, .1, 0), colours=c("#770000","red", "orange", "grey", "grey")) + base_theme+ theme(legend.position="top")
p.Prevotella.st
ggsave2("~/Desktop/PCoA/p.Prevotella.st.pdf", p.Prevotella.st, width = 3, height = 3.5, dpi = 300)

p.Phocaeicola.st <- cbind(ST.pcoa.axis,otutable_ST) %>%ggplot(aes(x=Axis.1, y=Axis.2, color=Phocaeicola)) + geom_point(size=0.3) + coord_equal() +
  theme_cowplot() + scale_colour_gradientn(values=c(1, .3, .2, .1, 0), colours=c("#770000","red", "orange", "grey", "grey")) + base_theme+ theme(legend.position="top")
p.Phocaeicola.st
ggsave2("~/Desktop/PCoA/p.Phocaeicola.st.pdf", p.Phocaeicola.st, width = 3, height = 3.5, dpi = 300)

p.Faecalibacterium.st <- cbind(ST.pcoa.axis,otutable_ST) %>%ggplot(aes(x=Axis.1, y=Axis.2, color=Faecalibacterium)) + geom_point(size=0.3) + coord_equal() +
  theme_cowplot() + scale_colour_gradientn(values=c(1, .3, .2, .1, 0), colours=c("#770000","red", "orange", "grey", "grey")) + base_theme+ theme(legend.position="top")
p.Faecalibacterium.st
ggsave2("~/Desktop/PCoA/p.Faecalibacterium.st.pdf", p.Faecalibacterium.st, width = 3, height = 3.5, dpi = 300)


#skin
physeq_C_SK_freq_filter
physeq_C_SK_ord <- ordinate(physeq_C_SK_freq_filter,method = "PCoA", distance = "bray")

p2 <- plot_ordination(physeq_C_SK_freq_filter, physeq_C_SK_ord, type="sample", color="IRIS")
p2 <- p2 + ggtitle("The Skin Microbiome") + theme_cowplot() + coord_equal() + base_theme + theme(legend.position="bottom") + scale_color_manual(values = iris_color)
p2
#ggsave2("~/Desktop/PCoA/Skin.PCoA.pdf", p2, width = 5, height = 6, dpi = 300)

SK.pcoa.axis <- physeq_C_SK_ord$vectors[,1:2] %>% data.frame()

p.Corynebacterium.sk <- cbind(SK.pcoa.axis,otutable_SK) %>%ggplot(aes(x=Axis.1, y=Axis.2, color=Corynebacterium)) + geom_point(size=0.3)  + coord_equal() +
  theme_cowplot() + scale_colour_gradientn(values=c(1, .3, .2, .1, 0), colours=c("#770000","red", "orange", "grey", "grey")) + base_theme + theme(legend.position="top")
p.Corynebacterium.sk  
ggsave2("~/Desktop/PCoA/p.Corynebacterium.sk.pdf", p.Corynebacterium.sk, width = 3, height = 3.5, dpi = 300)

p.Staphylococcus.sk <- cbind(SK.pcoa.axis,otutable_SK) %>%ggplot(aes(x=Axis.1, y=Axis.2, color=Staphylococcus)) + geom_point(size=0.3)  + coord_equal() +
  theme_cowplot() + scale_colour_gradientn(values=c(1, .3, .2, .1, 0), colours=c("#770000","red", "orange", "grey", "grey")) + base_theme + theme(legend.position="top")
p.Staphylococcus.sk  
ggsave2("~/Desktop/PCoA/p.Staphylococcus.sk.pdf", p.Staphylococcus.sk, width = 3, height = 3.5, dpi = 300)

p.Cutibacterium.sk <- cbind(SK.pcoa.axis,otutable_SK) %>%ggplot(aes(x=Axis.1, y=Axis.2, color=Cutibacterium)) + geom_point(size=0.3)  + coord_equal() +
  theme_cowplot() + scale_colour_gradientn(values=c(1, .5, .4, .2, 0), colours=c("#770000","red", "orange", "grey", "grey")) + base_theme + theme(legend.position="top")
p.Cutibacterium.sk  
ggsave2("~/Desktop/PCoA/p.Cutibacterium.sk.pdf", p.Cutibacterium.sk, width = 3, height = 3.5, dpi = 300)

p.Streptococcus.sk <- cbind(SK.pcoa.axis,otutable_SK) %>%ggplot(aes(x=Axis.1, y=Axis.2, color=Streptococcus)) + geom_point(size=0.3) + coord_equal() +
  theme_cowplot() + scale_colour_gradientn(values=c(1, .3, .2, .1, 0), colours=c("#770000","red", "orange", "grey", "grey")) + base_theme + theme(legend.position="top")
p.Streptococcus.sk
ggsave2("~/Desktop/PCoA/p.Streptococcus.sk.pdf", p.Streptococcus.sk, width = 3, height = 3.5, dpi = 300)


#oral
physeq_C_OR_freq_filter
physeq_C_OR_ord <- ordinate(physeq_C_OR_freq_filter,method = "PCoA", distance = "bray")

p3 <- plot_ordination(physeq_C_OR_freq_filter, physeq_C_OR_ord, type="sample", color="IRIS")
p3 <- p3 + ggtitle("The Oral Microbiome") + theme_cowplot() + coord_equal() + base_theme + theme(legend.position="bottom") + scale_color_manual(values = iris_color)
p3
ggsave2("~/Desktop/PCoA/Oral.PCoA.pdf", p3, width = 5, height = 6, dpi = 300)

OR.pcoa.axis <- physeq_C_OR_ord$vectors[,1:2] %>% data.frame()

p.Prevotella.or <- cbind(OR.pcoa.axis,otutable_OR) %>%ggplot(aes(x=Axis.1, y=Axis.2, color=Prevotella)) + geom_point(size=0.3) + coord_equal() +
  theme_cowplot() + scale_colour_gradientn(values=c(1, .6, .4, .2, 0), colours=c("#770000","red", "orange", "grey", "grey"))+ base_theme + theme(legend.position="top")
p.Prevotella.or  
ggsave2("~/Desktop/PCoA/p.Prevotella.or.pdf", p.Prevotella.or, width = 3, height = 3.5, dpi = 300)

p.Neisseria.or <- cbind(OR.pcoa.axis,otutable_OR) %>%ggplot(aes(x=Axis.1, y=Axis.2, color=Neisseria)) + geom_point(size=0.3) + coord_equal() +
  theme_cowplot() + scale_colour_gradientn(values=c(1, .3, .2, .1, 0), colours=c("#770000","red", "orange", "grey", "grey"))+ base_theme + theme(legend.position="top")
p.Neisseria.or 
ggsave2("~/Desktop/PCoA/p.Neisseria.or.pdf", p.Neisseria.or, width = 3, height = 3.5, dpi = 300)

p.Veillonella.or <- cbind(OR.pcoa.axis,otutable_OR) %>%ggplot(aes(x=Axis.1, y=Axis.2, color=Veillonella)) + geom_point(size=0.3) + coord_equal() +
  theme_cowplot() + scale_colour_gradientn(values=c(1, .3, .2, .1, 0), colours=c("#770000","red", "orange", "grey", "grey"))+ base_theme + theme(legend.position="top")
p.Veillonella.or
ggsave2("~/Desktop/PCoA/p.Veillonella.or.pdf", p.Veillonella.or, width = 3, height = 3.5, dpi = 300)

p.Haemophilus.or <- cbind(OR.pcoa.axis,otutable_OR) %>%ggplot(aes(x=Axis.1, y=Axis.2, color=Haemophilus)) + geom_point(size=0.3) + coord_equal() +
  theme_cowplot() + scale_colour_gradientn(values=c(1, .3, .2, .1, 0), colours=c("#770000","red", "orange", "grey", "grey"))+ base_theme + theme(legend.position="top")
p.Haemophilus.or  
ggsave2("~/Desktop/PCoA/p.Haemophilus.or.pdf", p.Haemophilus.or, width = 3, height = 3.5, dpi = 300)

#Nasal
physeq_C_NS_freq_filter
physeq_C_NS_ord <- ordinate(physeq_C_NS_freq_filter,method = "PCoA", distance = "bray")

p4 <- plot_ordination(physeq_C_NS_freq_filter, physeq_C_NS_ord, type="sample", color="IRIS")
p4 <- p4 + ggtitle("The Nasal Microbiome") + theme_cowplot() + coord_equal() + base_theme + theme(legend.position="bottom") + scale_color_manual(values = iris_color)
p4
ggsave2("~/Desktop/PCoA/Nasal.PCoA.pdf", p4, width = 5, height = 6, dpi = 300)

NS.pcoa.axis <- physeq_C_NS_ord$vectors[,1:2] %>% data.frame()

p.Corynebacterium.ns <- cbind(NS.pcoa.axis,otutable_NS) %>% ggplot(aes(x=Axis.1, y=Axis.2, color=Corynebacterium)) + geom_point(size=0.3) + coord_equal() +
  theme_cowplot() + scale_colour_gradientn(values=c(1, .5, .2, .1, 0), colours=c("#770000","red", "orange", "grey", "grey"))+ base_theme+ theme(legend.position="top")
p.Corynebacterium.ns 
ggsave2("~/Desktop/PCoA/p.Corynebacterium.ns.pdf", p.Corynebacterium.ns, width = 3, height = 3.5, dpi = 300)

p.Staphylococcus.ns <- cbind(NS.pcoa.axis,otutable_NS) %>% ggplot(aes(x=Axis.1, y=Axis.2, color=Staphylococcus)) + geom_point(size=0.3) + coord_equal() +
  theme_cowplot() + scale_colour_gradientn(values=c(1, .5, .2, .1, 0), colours=c("#770000","red", "orange", "grey", "grey"))+ base_theme+ theme(legend.position="top")
p.Staphylococcus.ns 
ggsave2("~/Desktop/PCoA/p.Staphylococcus.ns.pdf", p.Staphylococcus.ns, width = 3, height = 3.5, dpi = 300)

p.Cutibacterium.ns <- cbind(NS.pcoa.axis,otutable_NS) %>% ggplot(aes(x=Axis.1, y=Axis.2, color=Cutibacterium)) + geom_point(size=0.3) + coord_equal() +
  theme_cowplot() + scale_colour_gradientn(values=c(1, .5, .2, .1, 0), colours=c("#770000","red", "orange", "grey", "grey"))+ base_theme+ theme(legend.position="top")
p.Cutibacterium.ns 
ggsave2("~/Desktop/PCoA/p.Cutibacterium.ns.pdf", p.Cutibacterium.ns, width = 3, height = 3.5, dpi = 300)

p.Streptococcus.ns <- cbind(NS.pcoa.axis,otutable_NS) %>%ggplot(aes(x=Axis.1, y=Axis.2, color=Streptococcus)) + geom_point(size=0.3) + coord_equal() +
  theme_cowplot() + scale_colour_gradientn(values=c(1, .3, .2, .1, 0), colours=c("#770000","red", "orange", "grey", "grey")) + base_theme + theme(legend.position="top")
p.Streptococcus.ns
ggsave2("~/Desktop/PCoA/p.Streptococcus.ns.pdf", p.Streptococcus.ns, width = 3, height = 3.5, dpi = 300)

#########################<Combine Bodysites using sample ID for UMAP>#######################################################################################################################
physeq_C_ST_freq <- filter_taxa(physeq_C_ST_freq, function(x) mean(x) > 0, TRUE)
otutable_ST_C <- data.frame(t(otu_table(physeq_C_ST_freq)))
#tax_ST <- data.frame(tax_table(physeq_C_ST_freq))
#tax_ST$ASV <- rownames(tax_ST)
#colnames(otutable_ST_C) <- tax_ST$Genus[match(tax_ST$ASV,colnames(otutable_ST_C))]
colnames(otutable_ST_C) <- paste(colnames(otutable_ST_C), "_ST", sep="")
otutable_meta.ST <- hmp.sample.st %>% select(RandomID, SampleID) %>% merge(otutable_ST_C, by.x="RandomID", by.y ="row.names")
colnames(otutable_meta.ST)

physeq_C_SK_freq <- filter_taxa(physeq_C_SK_freq, function(x) mean(x) > 0, TRUE)
otutable_SK_C <- data.frame(t(otu_table(physeq_C_SK_freq)))
#tax_SK <- data.frame(tax_table(physeq_C_SK_freq))
#tax_SK$ASV <- rownames(tax_SK)
#colnames(otutable_SK_C) <- tax_SK$Genus[match(tax_SK$ASV,colnames(otutable_SK_C))]
colnames(otutable_SK_C) <- paste(colnames(otutable_SK_C), "_SK", sep="")
otutable_meta.SK <- hmp.sample.SK %>% select(RandomID, SampleID) %>% merge(otutable_SK_C, by.x="RandomID", by.y ="row.names")
colnames(otutable_meta.SK)

physeq_C_NS_freq <- filter_taxa(physeq_C_NS_freq, function(x) mean(x) > 0, TRUE)
otutable_NS_C <- data.frame(t(otu_table(physeq_C_NS_freq)))
#tax_NS <- data.frame(tax_table(physeq_C_NS_freq))
#tax_NS$ASV <- rownames(tax_NS)
#colnames(otutable_NS_C) <- tax_NS$Genus[match(tax_NS$ASV,colnames(otutable_NS_C))]
colnames(otutable_NS_C) <- paste(colnames(otutable_NS_C), "_NS", sep="")
otutable_meta.NS <- hmp.sample.NS %>% select(RandomID, SampleID) %>% merge(otutable_NS_C, by.x="RandomID", by.y ="row.names")
colnames(otutable_meta.NS)

physeq_C_OR_freq <- filter_taxa(physeq_C_OR_freq, function(x) mean(x) > 0, TRUE)
otutable_OR_C <- data.frame(t(otu_table(physeq_C_OR_freq)))
#tax_OR <- data.frame(tax_table(physeq_C_NS_freq))
#tax_OR$ASV <- rownames(tax_OR)
#colnames(otutable_OR_C) <- tax_OR$Genus[match(tax_OR$ASV,colnames(otutable_OR_C))]
colnames(otutable_OR_C) <- paste(colnames(otutable_OR_C), "_OR", sep="")
otutable_meta.OR <- hmp.sample.OR %>% select(RandomID, SampleID) %>% merge(otutable_OR_C, by.x="RandomID", by.y ="row.names")
colnames(otutable_meta.OR)

A1 <- merge(otutable_meta.ST, otutable_meta.SK, by="SampleID")
A2 <- merge(otutable_meta.OR, otutable_meta.NS, by="SampleID")

#make sure empty columns are already removed
table(select(otutable_meta.ST, -RandomID,-SampleID) %>% colSums() ==0)
table(select(otutable_meta.OR, -RandomID,-SampleID) %>% colSums() ==0)
table(select(otutable_meta.SK, -RandomID,-SampleID) %>% colSums() ==0)
table(select(otutable_meta.NS, -RandomID,-SampleID) %>% colSums() ==0)

OTU.Combined.Bodysite <- merge(A1, A2, by="SampleID")

row.names(OTU.Combined.Bodysite) <- OTU.Combined.Bodysite$SampleID
OTU.Combined.Bodysite <- select(OTU.Combined.Bodysite, -RandomID.x.x,-RandomID.x.y, -RandomID.y.x,-RandomID.y.y,-SampleID)
class(OTU.Combined.Bodysite)

OTU.Combined.Bodysite <- OTU.Combined.Bodysite[,colSums(OTU.Combined.Bodysite) != 0]

#save combined OTU table
#write.csv(file="../Genus Table/Combined_rel.abu.genus.csv", OTU.Combined.Bodysite)

OTU.Combined.Bodysite.su <- data.frame(t(OTU.Combined.Bodysite))
HMP_seurat_Combined <- CreateSeuratObject(OTU.Combined.Bodysite.su, project = "iHMP", min.cells = 1, min.features = 1)

hmp.sample.combined <- filter(hmp.sample.st, hmp.sample.st$SampleID %in% rownames(OTU.Combined.Bodysite))
hmp.sample.combined$names <- paste("X", hmp.sample.combined$SampleID, sep="")
hmp.sample.combined$names <- gsub("-",".",hmp.sample.combined$names)
rownames(hmp.sample.combined) <- hmp.sample.combined$names 
HMP_seurat_Combined <- AddMetaData(HMP_seurat_Combined, hmp.sample.combined)

HMP_seurat_Combined[["percent.unclassified"]] <- PercentageFeatureSet(HMP_seurat_Combined, pattern = "^Unclassified")
HMP_seurat_Combined

#transform to relative count per million
HMP_seurat_Combined <- NormalizeData(HMP_seurat_Combined, normalization.method = "RC", scale.factor = 1000000)

HMP_seurat_Combined <- FindVariableFeatures(HMP_seurat_Combined, selection.method = "vst", nfeatures = 500)
head(VariableFeatures(HMP_seurat_Combined), 10)

all.genes <- rownames(HMP_seurat_Combined)
HMP_seurat_Combined <- ScaleData(HMP_seurat_Combined, features = all.genes)

HMP_seurat_Combined <- RunPCA(HMP_seurat_Combined, features = VariableFeatures(object = HMP_seurat_Combined))
DimPlot(HMP_seurat_Combined, reduction = "pca")
ElbowPlot(HMP_seurat_Combined)

OTU_Combined_BCDist <- vegdist(OTU.Combined.Bodysite, method = "bray")
OTU_Combined_PCOA <- pcoa(OTU_Combined_BCDist)
OTU_Combined_PCOA$vectors[1:5, 1:5]
combined_bray <-  OTU_Combined_PCOA$vectors
row.names(combined_bray) <- paste("X", gsub("-", ".",row.names(combined_bray)), sep="")
HMP_seurat_Combined[["bray"]] <- CreateDimReducObject(embeddings = combined_bray, key = "Axis.", assay = DefaultAssay(HMP_seurat_Combined))
HMP_seurat_Combined <- RunUMAP(HMP_seurat_Combined, reduction="bray", dims = 1:40)
DimPlot(HMP_seurat_Combined, reduction = "umap",label = T)

HMP_seurat_Combined <- FindNeighbors(HMP_seurat_Combined, dims = 1:40)
HMP_seurat_Combined <- FindClusters(HMP_seurat_Combined, resolution = 0.5)
DimPlot(HMP_seurat_Combined)
DimPlot(HMP_seurat_Combined,group.by = "SubjectID")
DimPlot(HMP_seurat_Combined,split.by = "A1C_Class")
HMP_seurat_Combined[[]]
colSums(HMP_seurat_Combined)

DimPlot(HMP_seurat_Combined, pt.size = 0.2) + scale_color_npg()

GeneraList <- rownames(HMP_seurat_Combined)
GeneraList[str_detect(GeneraList, "Bacter")]

featurelist <- c("Corynebacterium","Staphylococcus","Cutibacterium","Bacteroides",
                 "Ruminococcus","Prevotella","Streptococcus","Veillonella",
                 "Haemophilus","Neisseria", "Moraxella","Unclassified-Clostridia")

featurelist2 <- c("Staphylococcus-ST", "Staphylococcus-SK", "Staphylococcus-OR","Staphylococcus-NS")
featurelist2 <- c("Prevotella-ST","Prevotella-SK", "Prevotella-OR","Prevotella-NS")
featurelist2 <- c("Cutibacterium-ST", "Cutibacterium-SK", "Cutibacterium-OR", "Cutibacterium-NS")

FeaturePlot(HMP_seurat_Combined, features=featurelist2, pt.size = 0.1)&theme(legend.position = "none")

df.subjectID <- data.frame(rownames(as.matrix(OTU_Combined_BCDist)))
colnames(df.subjectID)[1] <- "SampleID"
df.subjectID$SubjectID <- sub('^([^-]*-[^-]+).*', '\\1', df.subjectID$SampleID)

adonis(OTU_Combined_BCDist ~ SubjectID, data = df.subjectID)

metadata_ST <- data.frame(sample_data(physeq_C_ST_freq_filter))
adonis2(distance(physeq_C_ST_freq_filter, method="bray") ~ SubjectID + batch, data = metadata_ST)

p1.b <- plot_ordination(physeq_C_ST_freq_filter, physeq_C_ST_ord, type="sample", color="batch")
p1.b <- p1.b + geom_point(size = 0.1) + ggtitle("The Stool Microbiome by Batch") + theme_cowplot() + coord_equal() + base_theme + theme(legend.position="bottom")
p1.b

metadata_SK <- data.frame(sample_data(physeq_C_SK_freq_filter))
adonis2(distance(physeq_C_SK_freq_filter, method="bray") ~ SubjectID + batch, data = metadata_SK)

p2.b <- plot_ordination(physeq_C_SK_freq_filter, physeq_C_SK_ord, type="sample", color="batch")
p2.b <- p2.b +geom_point(size = 0.1) + ggtitle("The Skin Microbiome by Batch") + theme_cowplot() + coord_equal() + base_theme + theme(legend.position="bottom")
p2.b

metadata_OR <- data.frame(sample_data(physeq_C_OR_freq_filter))
adonis2(distance(physeq_C_OR_freq_filter, method="bray") ~ SubjectID + batch, data = metadata_OR)

p3.b <- plot_ordination(physeq_C_OR_freq_filter, physeq_C_OR_ord, type="sample", color="batch")
p3.b <- p3.b +geom_point(size = 0.1) + ggtitle("The Oral Microbiome by Batch") + theme_cowplot() + coord_equal() + base_theme + theme(legend.position="bottom")
p3.b

metadata_NS <- data.frame(sample_data(physeq_C_NS_freq_filter))
adonis2(distance(physeq_C_NS_freq_filter, method="bray") ~ SubjectID + batch, data = metadata_NS)

p4.b <- plot_ordination(physeq_C_NS_freq_filter, physeq_C_NS_ord, type="sample", color="batch")
p4.b <- p4.b +geom_point(size = 0.1) + ggtitle("The Nasal Microbiome by Batch") + theme_cowplot() + coord_equal() + base_theme + theme(legend.position="bottom")
p4.b

p1.b + p2.b + p3.b + p4.b

adonis2(distance(physeq_C_ST_freq_filter, method="bray") ~ IRIS + batch + SubjectID, data = metadata_ST)

adonis2(distance(physeq_C_SK_freq_filter, method="bray") ~ IRIS + batch, data = metadata_SK)
adonis2(distance(physeq_C_OR_freq_filter, method="bray") ~ SSPG, data = metadata_OR)
adonis2(distance(physeq_C_NS_freq_filter, method="bray") ~ SSPG, data = metadata_NS)





