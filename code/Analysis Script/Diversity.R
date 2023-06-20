#Diversity analysis
#This analysis is intended to understand the relationship between diversity and host physiology 

library(microbiome)
library(phyloseq)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(lme4)
library(sjPlot)
library(afex)
library(reshape)
library(patchwork)

base_theme = theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 13),
        panel.grid.minor = element_blank())

body_site_color = c(
  "Stool" = ggsci::pal_jama()(n=7)[2],
  "Skin" = ggsci::pal_jama()(n=7)[3],
  "Oral" = ggsci::pal_jama()(n=7)[4],
  "Nasal" = ggsci::pal_jama()(n=7)[5])


setwd("/Users/xzhou7/Library/CloudStorage/Box-Box/XinZhouFiles/Projects/ZXP2_HMP2_Core/Analysis/")
ls()
load("./Robject/DetailedPhyloseq.RData")
load("./Robject/Revision_MultiOmes_0509.RData")

sort(table(data.frame(tax_table(physeq_ST))$Genus),decreasing=T)[1:50]

#diversity estimate needs to use count data
physeq_ST <- filter_taxa(physeq_ST, function(x) sum(x) > 0, TRUE)
physeq_SK <- filter_taxa(physeq_SK, function(x) sum(x) > 0, TRUE)
physeq_OR <- filter_taxa(physeq_OR, function(x) sum(x) > 0, TRUE)
physeq_NS <- filter_taxa(physeq_NS, function(x) sum(x) > 0, TRUE)

#normalize to a minimum of 5000 case
set.seed(123)
physeq_ST_trim <- rarefy_even_depth(physeq_ST, sample.size = 5000, replace = T)
physeq_NS_trim <- rarefy_even_depth(physeq_NS, sample.size = 5000, replace = T)
physeq_SK_trim <- rarefy_even_depth(physeq_SK, sample.size = 5000, replace = T)
physeq_OR_trim <- rarefy_even_depth(physeq_OR, sample.size = 5000, replace = T)

#estimate diversity at genus level
#stool diversity by genus
Target_physeq <- physeq_ST_trim

ST.Genuslist <- unique(data.frame(tax_table(Target_physeq))$Genus)
richness_ST_byGenus <- data.frame()
diversity_ST_byGenus <- data.frame()

length(ST.Genuslist)

ntaxa(Target_physeq)


Target_physeq

for (i in 1:length(ST.Genuslist)){
  Genus.Target <- ST.Genuslist[i]
  print(Genus.Target)
  physeq_genus_temp <- subset_taxa(Target_physeq, Genus==Genus.Target )
  
  if(ntaxa(physeq_genus_temp) > 1){
    print("good")
    richness_temp <- estimate_richness(physeq_genus_temp, measures="Observed")
    diversity_temp <- estimate_richness(physeq_genus_temp, measures = "Shannon")
    
    richness_ST_byGenus[1:nsamples(Target_physeq), i] <- richness_temp
    colnames(richness_ST_byGenus)[i] <- Genus.Target
    
    diversity_ST_byGenus[1:nsamples(Target_physeq),i] <- diversity_temp
    colnames(diversity_ST_byGenus)[i] <- Genus.Target
  }
  
  else{
    print("bad")
    richness_ST_byGenus[1:nsamples(Target_physeq), i] <- 0
    colnames(richness_ST_byGenus)[i]<- Genus.Target
    
    diversity_ST_byGenus[1:nsamples(Target_physeq),i] <- 0
    colnames(diversity_ST_byGenus)[i] <- Genus.Target
  }
}

row.names(richness_ST_byGenus) <- row.names(richness_temp)
row.names(diversity_ST_byGenus) <- row.names(diversity_temp)

table(colSums(richness_ST_byGenus) == 0)
table(colSums(richness_ST_byGenus) == 1)

sort(colSums(richness_ST_byGenus), decreasing = T) [1:20]
sort(colMeans(richness_ST_byGenus), decreasing = T) [1:20]

table(colSums(diversity_ST_byGenus) == 0)

genus.list_ST <- names(sort(colMeans(diversity_ST_byGenus), decreasing = T)[1:20])

ST_Gene_Comp <- cbind(colMeans(richness_ST_byGenus),colMeans(diversity_ST_byGenus)) %>% data.frame()
colnames(ST_Gene_Comp) <- c("richness", "diversity")
ST_Gene_Comp <- filter(ST_Gene_Comp, diversity != 0 & richness != 0)

sort(ST_Gene_Comp$richness)

ST_Gene_Comp <- ST_Gene_Comp[order(ST_Gene_Comp$richness, decreasing = T),]
p <- ggplot(ST_Gene_Comp, aes(x=diversity, y=richness)) + geom_point() + ggtitle("Stool Microbiome Genetic Diversity")
p <- p + geom_text_repel(aes(label=ifelse(diversity > 1,rownames(ST_Gene_Comp),'')),hjust=0,vjust=0)
p

#ggsave(filename = "./CytokineRichness/Richness.stool.pdf", height=5, width=5, dpi=300)

taxa_ST <- data.frame(tax_table(Target_physeq)) %>% select(-Species) %>% unique()
sample_ST <- data.frame(sample_data(Target_physeq))

ST_diversity_cytokine <- merge(select(richness_ST_byGenus, genus.list_ST, Prevotella), sample_ST, by="row.names") %>% merge(ck.df, by="SampleID")
colnames(ST_diversity_cytokine)[33:94]
ST_diversity_cytokine_scaled <- ST_diversity_cytokine
ST_diversity_cytokine_scaled[,3:23] <- lapply(ST_diversity_cytokine[,3:23], scale) %>% data.frame()
ST_diversity_cytokine_scaled[,33:94] <- lapply(ST_diversity_cytokine[,33:94], scale) %>% data.frame()

ST_diversity_cytokine_scaled <- ST_diversity_cytokine_scaled %>% group_by(SubjectID.y) %>% mutate(Days = Date - min(Date))
  
# ST_diversity_cytokine$Unclassified_Lachnospiraceae_scaled <- scale(ST_diversity_cytokine$Unclassified_Lachnospiraceae)
#ST_diversity_cytokine$Unclassified_Ruminococcaceae_scaled <- scale(ST_diversity_cytokine$Unclassified_Ruminococcaceae)
# ST_diversity_cytokine$Bacteroides_scaled <- scale(ST_diversity_cytokine$Bacteroides)
# ST_diversity_cytokine$Prevotella_scaled <- scale(ST_diversity_cytokine$Prevotella)
# ST_diversity_cytokine$Phocaeicola_scaled <- scale(ST_diversity_cytokine$Phocaeicola)
# ST_diversity_cytokine$Alistipes_scaled <- scale(ST_diversity_cytokine$Alistipes)
# 
# ST_diversity_cytokine$LEPTIN_scaled <- scale(ST_diversity_cytokine$LEPTIN)
#ST_diversity_cytokine$IL17A_scaled <- scale(ST_diversity_cytokine$IL17A)
# ST_diversity_cytokine$IL17F_scaled <- scale(ST_diversity_cytokine$IL17F)
# ST_diversity_cytokine$IL22_scaled <- scale(ST_diversity_cytokine$IL22)
# ST_diversity_cytokine$IFNG_scaled <- scale(ST_diversity_cytokine$IFNG)
# 

#afex to operate the linear mixed model: https://cran.r-project.org/web/packages/afex/readme/README.html
genus_test_ST <- colnames(ST_diversity_cytokine)[3:23]
cytokine_list <- colnames(ST_diversity_cytokine)[33:94]

# p <- ggplot(data = ST_diversity_cytokine, aes(x=IL17F,y=Unclassified_Ruminococcaceae)) + geom_point() + geom_smooth(method=loess, se=FALSE)
# p
# 
# i <- "Unclassified_Ruminococcaceae_scaled"
# j <- "IL17A_scaled"

p.table.st <- data.frame()
beta.table.st <- data.frame()
for (i in 1:length(genus_test_ST)){
  for (j in 1:length(cytokine_list)){
  taxa <- genus_test_ST[i]
  print(taxa)
  cytokine <- cytokine_list[j]
  print(cytokine)
  basic.mix <- mixed(paste0(taxa, "~",cytokine, "+ Days + (1|SubjectID.y)"), data = ST_diversity_cytokine_scaled)
  
  p.table.st[j,i] <- basic.mix$anova_table$`Pr(>F)`[1]
  beta.table.st[j,i] <- basic.mix$full_model@beta[2]
  
  colnames(p.table.st)[i] <- taxa
  colnames(beta.table.st)[i] <- taxa
  
  row.names(p.table.st)[j] <- cytokine
  row.names(beta.table.st)[j] <- cytokine
  
  print(basic.mix$anova_table$`Pr(>F)`)[1]
  print(basic.mix$full_model@beta[2])
  }
}

p.table.st.melt<- melt(as.matrix(p.table.st))
beta.table.st.melt <- melt(as.matrix(beta.table.st))

stool_cytokine_richness <- p.table.st.melt

stool_cytokine_richness$beta <- beta.table.st.melt$value
colnames(stool_cytokine_richness) <- c("Cytokine","Richness","P","beta")
stool_cytokine_richness$adj.P <- stool_cytokine_richness$P * 21
stool_cytokine_richness$adj.P[stool_cytokine_richness$adj.P > 1] <- 1

# plot(density(stool_cytokine_richness$adj.P))
# pdf("./CytokineRichness/beta.density.stool.pdf", width = 5, height = 3)
# plot(density(stool_cytokine_richness$beta))
# dev.off()
# stool_cytokine_richness[(stool_cytokine_richness$adj.P < 0.05),]
# stool_cytokine_richness[(stool_cytokine_richness$adj.P < 0.2),]
# write.csv(file="./CytokineRichness/stool.richness.20.csv",stool_cytokine_richness)

# 
# mixed(Unclassified_Ruminococcaceae_scaled ~ IL17A_scaled + (1|SubjectID.y), data = ST_diversity_cytokine)
# basic.mix <- mixed(Unclassified_Ruminococcaceae_scaled ~ IL17A_scaled + (1|SubjectID.y), data = ST_diversity_cytokine)
# basic.mix$anova_table$`Pr(>F)`
# 
# em <- emmeans(basic.mix, "IL17A_scaled")
# plot(em)
# 
# basic.mix <- mixed(Alistipes_scaled ~ IFNG_scaled + (1|SubjectID.y), data = ST_diversity_cytokine)
# 
# basic.mix <- mixed(Phocaeicola_scaled ~ IL22_scaled + (1|SubjectID.y), data = ST_diversity_cytokine)
# 
# basic.mix <- mixed(Unclassified_Ruminococcaceae_scaled ~ IL22_scaled + (1|SubjectID.y), data = ST_diversity_cytokine)
# 
# basic.mix <- mixed(Phocaeicola_scaled ~ IL22_scaled + (1|SubjectID.y), data = ST_diversity_cytokine)
# 
# 
# lm.bs <- lm(Unclassified_Ruminococcaceae_scaled ~ IL17A_scaled, data = ST_diversity_cytokine)
# 
# summary(basic.mix)
# summary(lm.bs)
# 
# sjPlot::plot_model(basic.mix)
# 
# basiclm


#Skin 
Target_physeq <- physeq_SK_trim

Target_physeq

SK.Genuslist <- unique(data.frame(tax_table(Target_physeq))$Genus)

richness_SK_byGenus <- data.frame()
diversity_SK_byGenus <- data.frame()

for (i in 1:length(SK.Genuslist)){
  Genus.Target <- SK.Genuslist[i]
  print(Genus.Target)
  physeq_genus_temp <- subset_taxa(Target_physeq, Genus==Genus.Target )
  
  if(ntaxa(physeq_genus_temp) > 1){
    print("good")
    richness_temp <- estimate_richness(physeq_genus_temp, measures="ACE")
    diversity_temp <- estimate_richness(physeq_genus_temp, measures = "Shannon")
    
    richness_SK_byGenus[1:904, i] <- richness_temp
    colnames(richness_SK_byGenus)[i] <- Genus.Target
    
    diversity_SK_byGenus[1:904,i] <- diversity_temp
    colnames(diversity_SK_byGenus)[i] <- Genus.Target
  }
  
  else{
    print("bad")
    richness_SK_byGenus[1:904, i] <- 0
    colnames(richness_SK_byGenus)[i]<- Genus.Target
    
    diversity_SK_byGenus[1:904,i] <- 0
    colnames(diversity_SK_byGenus)[i] <- Genus.Target
  }
}

table(colSums(richness_SK_byGenus) == 0)
table(colSums(richness_SK_byGenus) == 1)

sort(colSums(richness_SK_byGenus), decreasing = T) [1:20]
sort(colMeans(richness_SK_byGenus), decreasing = T) [1:20]

table(colSums(diversity_SK_byGenus) == 0)

genus.list_SK <- names(sort(colMeans(diversity_SK_byGenus), decreasing = T)[1:20])

SK_Gene_Comp <- cbind(colMeans(richness_SK_byGenus),colMeans(diversity_SK_byGenus)) %>% data.frame()
colnames(SK_Gene_Comp) <- c("richness", "diversity")
SK_Gene_Comp <- filter(SK_Gene_Comp, diversity != 0 & richness != 0)

sort(SK_Gene_Comp$richness)

SK_Gene_Comp <- SK_Gene_Comp[order(SK_Gene_Comp$richness, decreasing = T),]
p.sk.di <- ggplot(SK_Gene_Comp, aes(x=diversity, y=richness)) + geom_point() + ggtitle("SKin Microbiome Genetic Diversity")
p.sk.di <- p.sk.di + geom_text_repel(aes(label=ifelse(diversity>0.1,rownames(SK_Gene_Comp),'')),hjust=0,vjust=0)
p.sk.di

#ggsave(filename = "./CytokineRichness/Richness.skin.pdf",p.sk.di, height=5, width=5, dpi=300)

taxa_SK <- data.frame(tax_table(Target_physeq)) %>% select(-Species) %>% unique()
sample_SK <- data.frame(sample_data(Target_physeq))

SK_diversity_cytokine <- merge(select(richness_SK_byGenus, genus.list_SK), sample_SK, by="row.names") %>% merge(ck.df, by="SampleID")
colnames(SK_diversity_cytokine)[30:91]
SK_diversity_cytokine_scaled <- SK_diversity_cytokine
SK_diversity_cytokine_scaled[,3:22] <- lapply(SK_diversity_cytokine[,3:22], scale) %>% data.frame()
SK_diversity_cytokine_scaled[,30:91] <- lapply(SK_diversity_cytokine[,30:91], scale) %>% data.frame()

SK_diversity_cytokine_scaled <- SK_diversity_cytokine_scaled %>% group_by(SubjectID.y) %>% mutate(Days = Date - min(Date))

#afex to operate the linear mixed model: https://cran.r-project.org/web/packages/afex/readme/README.html
genus_test_SK <- colnames(SK_diversity_cytokine)[3:22]
genus_test_SK
cytokine_list 

# p <- ggplot(data = SK_diversity_cytokine, aes(x=TNFA,y=Streptococcus)) + geom_point() + geom_smooth(method=loess, se=FALSE)
# p
# 
# i <- "Unclassified_Ruminococcaceae_scaled"
# j <- "IL17A_scaled"

p.table.sk <- data.frame()
beta.table.sk <- data.frame()
for (i in 1:length(genus_test_SK)){
  for (j in 1:length(cytokine_list)){
    taxa <- genus_test_SK[i]
    print(taxa)
    cytokine <- cytokine_list[j]
    print(cytokine)
    basic.mix <- mixed(paste0(taxa, "~",cytokine, "+ Days + (1|SubjectID.y)"), data = SK_diversity_cytokine_scaled)
    
    p.table.sk[j,i] <- basic.mix$anova_table$`Pr(>F)`[1]
    beta.table.sk[j,i] <- basic.mix$full_model@beta[2]
    
    colnames(p.table.sk)[i] <- taxa
    colnames(beta.table.sk)[i] <- taxa
    
    row.names(p.table.sk)[j] <- cytokine
    row.names(beta.table.sk)[j] <- cytokine
    
    print(basic.mix$anova_table$`Pr(>F)`)[1]
    print(basic.mix$full_model@beta[2])
  }
}

p.table.sk.melt<- melt(as.matrix(p.table.sk))
beta.table.sk.melt <- melt(as.matrix(beta.table.sk))

skin_cytokine_richness <- p.table.sk.melt

skin_cytokine_richness$beta <- beta.table.sk.melt$value
colnames(skin_cytokine_richness) <- c("Cytokine","Richness","P","beta")
skin_cytokine_richness$adj.P <- skin_cytokine_richness$P * 20
skin_cytokine_richness$adj.P[skin_cytokine_richness$adj.P > 1] <- 1

plot(density(skin_cytokine_richness$adj.P))
pdf("./CytokineRichness/beta.density.skin.pdf", width = 5, height = 3)
plot(density(skin_cytokine_richness$beta))
dev.off()
skin_cytokine_richness[(skin_cytokine_richness$adj.P < 0.05),]
skin_cytokine_richness[(skin_cytokine_richness$adj.P < 0.2),]

write.csv(file = "./CytokineRichness/skin.richness.20.csv",skin_cytokine_richness)

#Oral
Target_physeq <- physeq_OR_trim

Target_physeq

OR.Genuslist <- unique(data.frame(tax_table(Target_physeq))$Genus)

richness_OR_byGenus <- data.frame()
diversity_OR_byGenus <- data.frame()

for (i in 1:length(OR.Genuslist)){
  Genus.Target <- OR.Genuslist[i]
  print(Genus.Target)
  physeq_genus_temp <- subset_taxa(Target_physeq, Genus==Genus.Target )
  
  if(ntaxa(physeq_genus_temp) > 1){
    print("good")
    richness_temp <- estimate_richness(physeq_genus_temp, measures="ACE")
    diversity_temp <- estimate_richness(physeq_genus_temp, measures = "Shannon")
    
    richness_OR_byGenus[1:949, i] <- richness_temp
    colnames(richness_OR_byGenus)[i] <- Genus.Target
    
    diversity_OR_byGenus[1:949,i] <- diversity_temp
    colnames(diversity_OR_byGenus)[i] <- Genus.Target
  }
  
  else{
    print("bad")
    richness_OR_byGenus[1:949, i] <- 0
    colnames(richness_OR_byGenus)[i]<- Genus.Target
    
    diversity_OR_byGenus[1:949,i] <- 0
    colnames(diversity_OR_byGenus)[i] <- Genus.Target
  }
}

row.names(richness_OR_byGenus) <- row.names(richness_temp)
row.names(diversity_OR_byGenus) <- row.names(diversity_temp)

table(colSums(richness_OR_byGenus) == 0)
table(colSums(richness_OR_byGenus) == 1)

sort(colSums(richness_OR_byGenus), decreasing = T) [1:20]
sort(colMeans(richness_OR_byGenus), decreasing = T) [1:20]

table(colSums(diversity_OR_byGenus) == 0)

genus.list_OR <- names(sort(colMeans(diversity_OR_byGenus), decreasing = T)[1:20])

OR_Gene_Comp <- cbind(colMeans(richness_OR_byGenus),colMeans(diversity_OR_byGenus)) %>% data.frame()
colnames(OR_Gene_Comp) <- c("richness", "diversity")
OR_Gene_Comp <- filter(OR_Gene_Comp, diversity != 0 & richness != 0)

sort(OR_Gene_Comp$richness)

OR_Gene_Comp <- OR_Gene_Comp[order(OR_Gene_Comp$richness, decreasing = T),]
p.or.di <- ggplot(OR_Gene_Comp, aes(x=diversity, y=richness)) + geom_point() + ggtitle("Oral Microbiome Genetic Diversity")
p.or.di <- p.or.di + geom_text_repel(aes(label=ifelse(diversity>0.25,rownames(OR_Gene_Comp),'')),hjust=0,vjust=0)
p.or.di

ggsave(filename = "./CytokineRichness/Richness.oral.pdf",p.or.di, height=5, width=5, dpi=300)

taxa_OR <- data.frame(tax_table(Target_physeq)) %>% select(-Species) %>% unique()
sample_OR <- data.frame(sample_data(Target_physeq))

OR_diversity_cytokine <- merge(select(richness_OR_byGenus, genus.list_OR), sample_OR, by="row.names") %>% merge(ck.df, by="SampleID")

colnames(OR_diversity_cytokine)[30:91]
OR_diversity_cytokine_scaled <- OR_diversity_cytokine
OR_diversity_cytokine_scaled[,3:22] <- lapply(OR_diversity_cytokine[,3:22], scale) %>% data.frame()
OR_diversity_cytokine_scaled[,30:91] <- lapply(OR_diversity_cytokine[,30:91], scale) %>% data.frame()

OR_diversity_cytokine_scaled <- OR_diversity_cytokine_scaled %>% group_by(SubjectID.y) %>% mutate(Days = Date - min(Date))

#afex to operate the linear mixed model: https://cran.r-project.org/web/packages/afex/readme/README.html
genus_test_OR <- colnames(OR_diversity_cytokine)[3:22]
genus_test_OR
cytokine_list 

# p <- ggplot(data = SK_diversity_cytokine, aes(x=TNFA,y=Streptococcus)) + geom_point() + geom_smooth(method=loess, se=FALSE)
# p
# 
# i <- "Unclassified_Ruminococcaceae_scaled"
# j <- "IL17A_scaled"

p.table.or <- data.frame()
beta.table.or <- data.frame()
for (i in 1:length(genus_test_OR)){
  for (j in 1:length(cytokine_list)){
    taxa <- genus_test_OR[i]
    print(taxa)
    cytokine <- cytokine_list[j]
    print(cytokine)
    basic.mix <- mixed(paste0(taxa, "~",cytokine, "+ Days + (1|SubjectID.y)"), data = OR_diversity_cytokine_scaled)
    
    p.table.or[j,i] <- basic.mix$anova_table$`Pr(>F)`[1]
    beta.table.or[j,i] <- basic.mix$full_model@beta[2]
    
    colnames(p.table.or)[i] <- taxa
    colnames(beta.table.or)[i] <- taxa
    
    row.names(p.table.or)[j] <- cytokine
    row.names(beta.table.or)[j] <- cytokine
    
    print(basic.mix$anova_table$`Pr(>F)`[1])
    print(basic.mix$full_model@beta[2])
  }
}

p.table.or.melt<- melt(as.matrix(p.table.or))
beta.table.or.melt <- melt(as.matrix(beta.table.or))

oral_cytokine_richness <- p.table.or.melt

oral_cytokine_richness$beta <- beta.table.or.melt$value
colnames(oral_cytokine_richness) <- c("Cytokine","Richness","P","beta")
oral_cytokine_richness$adj.P <- oral_cytokine_richness$P * 20
oral_cytokine_richness$adj.P[oral_cytokine_richness$adj.P > 1] <- 1

plot(density(oral_cytokine_richness$adj.P))
pdf("./CytokineRichness/beta.density.oral.pdf", width = 5, height = 3)
plot(density(oral_cytokine_richness$beta))
dev.off()
oral_cytokine_richness[(oral_cytokine_richness$adj.P < 0.05),]
oral_cytokine_richness[(oral_cytokine_richness$adj.P < 0.2),]

#write.csv(file = "./CytokineRichness/oral.richness.20.csv",oral_cytokine_richness)

#Nasal
Target_physeq <- physeq_NS_trim

physeq_NS_trim

NS.Genuslist <- unique(data.frame(tax_table(Target_physeq))$Genus)

richness_NS_byGenus <- data.frame()
diversity_NS_byGenus <- data.frame()

for (i in 1:length(NS.Genuslist)){
  Genus.Target <- NS.Genuslist[i]
  print(Genus.Target)
  physeq_genus_temp <- subset_taxa(Target_physeq, Genus==Genus.Target)
  
  if(ntaxa(physeq_genus_temp) > 1){
    print("good")
    richness_temp <- estimate_richness(physeq_genus_temp, measures="ACE")
    diversity_temp <- estimate_richness(physeq_genus_temp, measures = "Shannon")
    
    richness_NS_byGenus[1:908, i] <- richness_temp
    colnames(richness_NS_byGenus)[i] <- Genus.Target
    
    diversity_NS_byGenus[1:908,i] <- diversity_temp
    colnames(diversity_NS_byGenus)[i] <- Genus.Target
  }
  
  else{
    print("bad")
    richness_NS_byGenus[1:908, i] <- 0
    colnames(richness_NS_byGenus)[i]<- Genus.Target
    
    diversity_NS_byGenus[1:908,i] <- 0
    colnames(diversity_NS_byGenus)[i] <- Genus.Target
  }
}

table(colSums(richness_NS_byGenus) == 0)
table(colSums(richness_NS_byGenus) == 1)

sort(colSums(richness_NS_byGenus), decreasing = T) [1:20]
sort(colMeans(richness_NS_byGenus), decreasing = T) [1:20]

table(colSums(diversity_NS_byGenus) == 0)

genus.list_NS <- names(sort(colMeans(diversity_NS_byGenus), decreasing = T)[1:20])

NS_Gene_Comp <- cbind(colMeans(richness_NS_byGenus),colMeans(diversity_NS_byGenus)) %>% data.frame()
colnames(NS_Gene_Comp) <- c("richness", "diversity")
NS_Gene_Comp <- filter(NS_Gene_Comp, diversity != 0 & richness != 0)

sort(NS_Gene_Comp$richness)

NS_Gene_Comp <- NS_Gene_Comp[order(NS_Gene_Comp$richness, decreasing = T),]
p.ns.di <- ggplot(NS_Gene_Comp, aes(x=diversity, y=richness)) + geom_point() + ggtitle("Nasal Microbiome Genetic Diversity")
p.ns.di <- p.ns.di + geom_text_repel(aes(label=ifelse(diversity>0.25,rownames(NS_Gene_Comp),'')),hjust=0,vjust=0)
p.ns.di

ggsave(filename = "./CytokineRichness/Richness.nasal.pdf", p.ns.di, height=5, width=5, dpi=300)

p.sk.di + p.ns.di

taxa_NS <- data.frame(tax_table(Target_physeq)) %>% select(-Species) %>% unique()
sample_NS <- data.frame(sample_data(Target_physeq))

NS_diversity_cytokine <- merge(select(richness_NS_byGenus, genus.list_NS), sample_NS, by="row.names") %>% merge(ck.df, by="SampleID")

colnames(NS_diversity_cytokine)[32:93]
NS_diversity_cytokine_scaled <- NS_diversity_cytokine
NS_diversity_cytokine_scaled[,3:22] <- lapply(NS_diversity_cytokine[,3:22], scale) %>% data.frame()
NS_diversity_cytokine_scaled[,32:93] <- lapply(NS_diversity_cytokine[,32:92], scale) %>% data.frame()

NS_diversity_cytokine_scaled <- NS_diversity_cytokine_scaled %>% group_by(SubjectID.y) %>% mutate(Days = Date - min(Date))

#afex to operate the linear mixed model: https://cran.r-project.org/web/packages/afex/readme/README.html
genus_test_NS <- colnames(NS_diversity_cytokine)[3:22]
genus_test_NS
cytokine_list 

# p <- ggplot(data = NS_diversity_cytokine, aes(x=IL17A,y=Cutibacterium)) + geom_point() + geom_smooth(method=loess, se=FALSE)
# p
# 
# i <- "Unclassified_Ruminococcaceae_scaled"
# j <- "IL17A_scaled"
# 

p.table.ns <- data.frame()
beta.table.ns <- data.frame()
for (i in 1:length(genus_test_NS)){
  for (j in 1:length(cytokine_list)){
    taxa <- genus_test_NS[i]
    print(taxa)
    cytokine <- cytokine_list[j]
    print(cytokine)
    basic.mix <- mixed(paste0(taxa, "~",cytokine, "+ Days + (1|SubjectID.y)"), data = NS_diversity_cytokine_scaled)
    
    p.table.ns[j,i] <- basic.mix$anova_table$`Pr(>F)`[1]
    beta.table.ns[j,i] <- basic.mix$full_model@beta[2]
    
    colnames(p.table.ns)[i] <- taxa
    colnames(beta.table.ns)[i] <- taxa
    
    row.names(p.table.ns)[j] <- cytokine
    row.names(beta.table.ns)[j] <- cytokine
    
    print(basic.mix$anova_table$`Pr(>F)`[1])
    print(basic.mix$full_model@beta[2])
  }
}

p.table.ns.melt<- melt(as.matrix(p.table.ns))
beta.table.ns.melt <- melt(as.matrix(beta.table.ns))

nasal_cytokine_richness <- p.table.ns.melt

nasal_cytokine_richness$beta <- beta.table.ns.melt$value
colnames(nasal_cytokine_richness) <- c("Cytokine","Richness","P","beta")
nasal_cytokine_richness$adj.P <- nasal_cytokine_richness$P * 20
nasal_cytokine_richness$adj.P[nasal_cytokine_richness$adj.P > 1] <- 1

plot(density(nasal_cytokine_richness$adj.P))
pdf("./CytokineRichness/beta.density.nasal.pdf", width = 5, height = 3)
plot(density(nasal_cytokine_richness$beta))
dev.off()
nasal_cytokine_richness[(nasal_cytokine_richness$adj.P < 0.05),]
nasal_cytokine_richness[(nasal_cytokine_richness$adj.P < 0.2),]

#write.csv(file = "./CytokineRichness/nasal.richness.20.csv",nasal_cytokine_richness)

stool_cytokine_richness[(stool_cytokine_richness$adj.P < 0.05),]
skin_cytokine_richness[(skin_cytokine_richness$adj.P < 0.05),]
oral_cytokine_richness[(oral_cytokine_richness$adj.P < 0.05),]
nasal_cytokine_richness[(nasal_cytokine_richness$adj.P < 0.05),]

getwd()

# save(richness_ST_byGenus, richness_SK_byGenus,richness_OR_byGenus, richness_NS_byGenus,
#         file = "./Robject/Richness_Bygenus.RData")

save(richness_ST_byGenus, richness_SK_byGenus,richness_OR_byGenus, richness_NS_byGenus,
     diversity_ST_byGenus, diversity_SK_byGenus, diversity_OR_byGenus, diversity_NS_byGenus,
     file = "./Robject/Diversity_Bygenus.RData")



######################
#Global Diversity
######################

#Stool
diversity.estimate.ST <- estimate_richness(physeq_ST_trim)
evenness.estimate.ST <- evenness(physeq_ST_trim, index = "all", zeroes = TRUE, detection = 0)
colnames(evenness.estimate.ST) <- paste(colnames(evenness.estimate.ST), "e", sep="_")
stool.diversity <- cbind(diversity.estimate.ST,evenness.estimate.ST)
colnames(stool.diversity)
write.csv(file = "./Diversity Table/Stool.Diversity.csv", stool.diversity)
stool.diversity$Bodysite <- "Stool"


#Skin
diversity.estimate.SK <- estimate_richness(physeq_SK_trim)
evenness.estimate.SK <- evenness(physeq_SK_trim, index = "all", zeroes = TRUE, detection = 0)
colnames(evenness.estimate.SK) <- paste(colnames(evenness.estimate.SK), "e", sep="_")
Skin.diversity <- cbind(diversity.estimate.SK, evenness.estimate.SK)
colnames(Skin.diversity)
write.csv(file = "./Diversity Table/Skin.Diversity.csv", Skin.diversity)
Skin.diversity$Bodysite <- "Skin"

#Oral
diversity.estimate.OR <- estimate_richness(physeq_OR_trim)
evenness.estimate.OR <- evenness(physeq_OR_trim, index = "all", zeroes = TRUE, detection = 0)
colnames(evenness.estimate.OR) <- paste(colnames(evenness.estimate.OR), "e", sep="_")
Oral.diversity <- cbind(diversity.estimate.OR,evenness.estimate.OR)
colnames(Oral.diversity)
write.csv(file = "./Diversity Table/Oral.Diversity.csv", Oral.diversity)
Oral.diversity$Bodysite <- "Oral"

#Nasal
diversity.estimate.NS <- estimate_richness(physeq_NS_trim)
evenness.estimate.NS <- evenness(physeq_NS_trim, index = "all", zeroes = TRUE, detection = 0)
colnames(evenness.estimate.NS) <- paste(colnames(evenness.estimate.NS), "e", sep="_")
Nasal.diversity <- cbind(diversity.estimate.NS,evenness.estimate.NS)
colnames(Nasal.diversity)
write.csv(file = "./Diversity Table/Nasal.Diversity.csv", Nasal.diversity)
Nasal.diversity$Bodysite <- "Nasal"

Diversity.Bysample <- rbind(stool.diversity, Skin.diversity, Oral.diversity, Nasal.diversity)

Diversity.Bysample
write.csv(file = "./Diversity Table/All_Diversity.csv",Diversity.Bysample)

colnames(Diversity.Bysample)

p.diver <- ggplot(Diversity.Bysample, aes(x=ACE, y=pielou_e, group=Bodysite, color=Bodysite)) + stat_density_2d(adjust=1.5)
p.diver <- p.diver + base_theme + scale_color_manual(values = body_site_color) + labs(title = "Ecological Characters of Microbiome", x="Richness(Observed)", y="Evenness(Pielou)")
p.diver <- p.diver + theme(
  legend.position = c(.95, .05),
  legend.justification = c("right", "bottom"),
  legend.box.just = "right",
  legend.margin = margin(6, 6, 6, 6),
  legend.box.background = element_rect(color="black", size=0.5),)
p.diver

ggsave(filename = "./Diversity Table/Diversity.pdf", p.diver, width = 4, height = 3, dpi = 300)





