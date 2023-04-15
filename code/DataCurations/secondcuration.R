library(phyloseq)
library(ggplot2)

setwd("~/Box/XinZhouFiles/Projects/ZXP2_HMP2_Core/Analysis/")

load("~/Box/XinZhouFiles/Projects/ZXP2_HMP2_Core/OLD HMP Files/Analysis/Analysis/HMP2_phyloseq_Objects.RData")
load("./Robject/PhyloseqObject.RData")

physeqord <- ordinate(physeqmean, "PCoA", dist = "bray")
p <- plot_ordination(physeqmean, physeqord  , type = "samples") #, color = "freqgroup" + scale_colour_manual(values = rhg_cols20) 
p <- p + geom_point(size = 0.5)
#p <- p + geom_text_repel(mapping = aes(label = SubjectID), size = 3,fontface = "bold")
p <- p + facet_wrap(~ V9,  drop=T) #ncol = 2,
#p <- p + scale_color_manual(values = c("blue", "gray","red","black")) #+ scale_alpha_manual(values = c(0.1,0.1,0.9))
#p <- p + stat_ellipse(mapping=aes(group=x),type="norm", linetype=2, show.legend=T)+ stat_ellipse(type="t")
p

sample_data()

physeq_UBM_skin <- subset_samples(physeq_UBM, SampleType=="Skin")
physeq_UBM_oral <- subset_samples(physeq_UBM, SampleType=="Mouth")
physeq_HMP_stool <- subset_samples(physeq_HMP, V3 == "ST")
physeq_HMP_nasal <- subset_samples(physeq_HMP, V3 == "NS")

physeq_HMP_stool_genux <- tax_glom(physeq_HMP_stool, taxrank = "Species")

physequbm_ord <- ordinate(physeq_HMP_stool_genux, "PCoA", dist = "bray")
p <- plot_ordination(physeq_HMP_stool_genux, physequbm_ord  , type = "samples") #, color = "freqgroup" + scale_colour_manual(values = rhg_cols20) 
p <- p + geom_point(size = 0.5)
#p <- p + geom_text_repel(mapping = aes(label = SubjectID), size = 3,fontface = "bold")
p <- p + facet_wrap(~ V9,  drop=T) #ncol = 2,
#p <- p + scale_color_manual(values = c("blue", "gray","red","black")) #+ scale_alpha_manual(values = c(0.1,0.1,0.9))
#p <- p + stat_ellipse(mapping=aes(group=x),type="norm", linetype=2, show.legend=T)+ stat_ellipse(type="t")
p

p1 <- plot_ordination(physeq_HMP_stool_genux, physequbm_ord, type="taxa", color="Phylum", title="taxa")
p1 + p
#save(physeq_HMP, physeq_UBM)

physeq_stool_69023 <- subset_samples(physeq_HMP_stool, SubjectID == "69-023")
physeq_stool_69023_genus <- tax_glom(physeq_stool_69023, taxrank = "Genus")

physeq_stool_69023_genus <- prune_taxa(names(sort(taxa_sums(physeq_stool_69023_genus),TRUE)[1:100]), physeq_stool_69023_genus)
sample_data(physeq_stool_69023_genus)
plot_heatmap(physeq_stool_69023_genus, taxa.label ="Genus", sample.label = "SampleID")

tax_table(physeq_stool_69023_genus) 

colnames(tax_table(physeq_stool_69023_genus))
plot_net(physeq_stool_69023_genus, maxdist=0.5)

Bacteroidetes_69023 <- subset_taxa(physeq_stool_69023, Genus == "Bacteroides")
tax_table(Bacteroidetes_69023) 
Bacteroidetes_69023_100 <- prune_taxa(names(sort(taxa_sums(Bacteroidetes_69023),TRUE)[1:100]), Bacteroidetes_69023)

plot_heatmap(Bacteroidetes_69023_100, taxa.label ="Species", sample.label = "SampleID")


subset(JAXmeta, SubjectID == "69-023")

subset(ls, SubjectID == "69-023")


subset(sample_data(physeq_UBM), SubjectID== "69-023")

sampledata_ALL 

#1remove Mikes sample beforeXXX date

sample_data(physeq)

sampleDataA <- data.frame(sample_data(physeq_UBM))
sampleDataB <- data.frame(sample_data(physeq_HMP))
 
sampledata_all <- merge(sampleDataA,sampleDataB,by="SampleID", all = T)
sampledata_all_ls <- merge(sampledata_all, ls, by = "SampleID", all = T)

p <- ggplot(sampleDataA,aes(x=Date, y=SubjectID)) + geom_point()
p

p1 <- ggplot(sampleDataB,aes(x=Date, y=SubjectID)) + geom_point()
p1

#write.csv(sampledata_all_ls, file = "~/Desktop/sampleDataALL.csv")

SampleDate <- read.csv("~/Desktop/SampleDate_Xin.csv",header = T)
SampleDate <- unique(SampleDate)
SampleDate$Date <- as.Date(SampleDate$Date , format = "%m/%d/%y")
SampleDate$Subject <- gsub("^([^-]*-[^-]*)-.*$", "\\1", SampleDate$SampleID)


write.csv(SampleDate, file = "./ori_meta_table/DateInfo.csv")





