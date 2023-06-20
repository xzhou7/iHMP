#cytokine summarize
library(x)

library(phyloseq)
library(ggplot2)
library(dplyr)
library(readxl)
library(patchwork)
library(ggrepel)
library(ggridges)
library(reshape2)
library(ggstatsplot)

setwd("~/Library/CloudStorage/Box-Box/XinZhouFiles/Projects/ZXP2_HMP2_Core/Analysis/")

base_theme = theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 13),
        panel.grid.minor = element_blank())
body_site_color = c(
  "Stool" = ggsci::pal_jama()(n=7)[2],
  "Skin" = ggsci::pal_jama()(n=7)[3],
  "Oral" = ggsci::pal_jama()(n=7)[4],
  "Nasal" = ggsci::pal_jama()(n=7)[5])

source("./Analysis Script/sxt.tools.R")

load("~/Library/CloudStorage/Box-Box/XinZhouFiles/Projects/ZXP2_HMP2_Core/Analysis/Robject/DetailedPhyloseq.RData")
load("~/Library/CloudStorage/Box-Box/XinZhouFiles/Projects/ZXP2_HMP2_Core/Analysis/Robject/Prevalance.RData")
taxa_ST <- data.frame(tax_table(physeq_ST))
taxa_SK <- data.frame(tax_table(physeq_SK))

stool.df <- read_xlsx("../../../../human_microbiome_project/data_analysis/brms/stool/result/stool_cytokine_beta.xlsx")
skin.df <- read_xlsx("../../../../human_microbiome_project/data_analysis/brms/skin/result/skin_cytokine_beta.xlsx")
oral.df <- read_xlsx("../../../../human_microbiome_project/data_analysis/brms/oral/result/oral_cytokine_beta.xlsx")
nasal.df <- read_xlsx("../../../../human_microbiome_project/data_analysis/brms/nasal/result/nasal_cytokine_beta.xlsx")

#get stool taxa distribution
stool.df$phylum <- "Xother"
stool.df$phylum[stool.df$genus %in% filter(taxa_ST, Phylum=="Bacteroidetes")$Genus] <- "Bacteroidetes"
stool.df$phylum[stool.df$genus %in% filter(taxa_ST, Phylum=="Firmicutes")$Genus] <- "Firmicutes"
#stool.df$phylum[stool.df$genus %in% filter(taxa_ST, Family=="Lachnospiraceae")$Genus] <- "Fam_Lachnospiraceae"
#stool.df$phylum[stool.df$genus %in% filter(taxa_ST, Family=="Ruminococcaceae")$Genus] <- "Fam_Ruminococcaceae"
stool.df$phylum[stool.df$genus %in% filter(taxa_ST, Phylum=="Proteobacteria")$Genus] <- "Proteobacteria"
stool.df$phylum[stool.df$genus %in% filter(taxa_ST, Phylum=="Actinobacteria")$Genus] <- "Actinobacteria"

stool.df$B_class <- "Xother"
stool.df$B_class[stool.df$genus %in% filter(taxa_ST, Class == "Clostridia")$Genus] <- "Clostridia"

tax.st <- prune_taxa(taxa_sums(physeqGenus_ST) > 0, physeqGenus_ST) %>% tax_table () %>% data.frame()
tax.st.total <- table(stool.df$phylum) %>% as.data.frame()

all.taxa.st <- table(tax.st$Phylum) %>% sort(decreasing = T) %>% data.frame()
sum(all.taxa.st$Freq) - sum((all.taxa.st$Freq)[1:4])

tax.st.total$total <- c(31,42,193,31,70)
colnames(tax.st.total) <- c("phylum", "Cytokine_Associated", "Total")

#write.csv(file = "../../../../human_microbiome_project/Figures/Figure4/stool.taxa.cytokine.csv",tax.st.total)

#get skin taxa distribution
skin.df$phylum <- "Xother"
skin.df$phylum[skin.df$genus %in% filter(taxa_SK, Phylum=="Bacteroidetes")$Genus] <- "Bacteroidetes"
skin.df$phylum[skin.df$genus %in% filter(taxa_SK, Phylum=="Firmicutes")$Genus] <- "Firmicutes"
skin.df$phylum[skin.df$genus %in% filter(taxa_SK, Phylum=="Proteobacteria")$Genus] <- "Proteobacteria"
skin.df$phylum[skin.df$genus %in% filter(taxa_SK, Phylum=="Actinobacteria")$Genus] <- "Actinobacteria"

skin.df$B_class <- "Xother"
skin.df$B_class[skin.df$genus %in% filter(taxa_SK, Class == "Clostridia")$Genus] <- "Clostridia"

tax.sk <- prune_taxa(taxa_sums(physeqGenus_SK) > 0, physeqGenus_SK) %>% tax_table () %>% data.frame()
tax.sk.total <- table(skin.df$phylum) %>% as.data.frame()

all.taxa.sk <- table(tax.sk$Phylum) %>% sort(decreasing = T) %>% data.frame()
sum(all.taxa.sk$Freq) - sum((all.taxa.sk$Freq)[1:4])

all.taxa.sk
tax.sk.total

tax.sk.total$total <- c(204,119,262,164,363)

colnames(tax.sk.total) <- c("phylum", "Cytokine_Associated", "Total")

#write.csv(file = "../../../../human_microbiome_project/Figures/Figure4/figure4d/skin.taxa.cytokine.csv",tax.sk.total)

#get oral taxa distribution, taxa table is the same as Skin microbiome
oral.df$phylum <- "Xother"
oral.df$phylum[oral.df$genus %in% filter(taxa_SK, Phylum=="Bacteroidetes")$Genus] <- "Bacteroidetes"
oral.df$phylum[oral.df$genus %in% filter(taxa_SK, Phylum=="Firmicutes")$Genus] <- "Firmicutes"
oral.df$phylum[oral.df$genus %in% filter(taxa_SK, Phylum=="Proteobacteria")$Genus] <- "Proteobacteria"
oral.df$phylum[oral.df$genus %in% filter(taxa_SK, Phylum=="Actinobacteria")$Genus] <- "Actinobacteria"

oral.df$B_class <- "Xother"
oral.df$B_class[oral.df$genus %in% filter(taxa_SK, Class == "Clostridia")$Genus] <- "Clostridia"

tax.or <- prune_taxa(taxa_sums(physeqGenus_OR) > 0, physeqGenus_OR) %>% tax_table () %>% data.frame()
tax.or.total <- table(oral.df$phylum) %>% as.data.frame()

all.taxa.or <- table(tax.or$Phylum) %>% sort(decreasing = T) %>% data.frame()
sum(all.taxa.or$Freq) - sum((all.taxa.or$Freq)[1:4])

all.taxa.or
tax.or.total

tax.or.total$total <- c(66,53,114,60,117)

colnames(tax.or.total) <- c("phylum", "Cytokine_Associated", "Total")

#write.csv(file = "../../../../human_microbiome_project/Figures/Figure4/figure4d/oral.taxa.cytokine.csv",tax.or.total)

#get nasal taxa distribution
nasal.df$phylum <- "Xother"
nasal.df$phylum[nasal.df$genus %in% filter(taxa_ST, Phylum=="Bacteroidetes")$Genus] <- "Bacteroidetes"
nasal.df$phylum[nasal.df$genus %in% filter(taxa_ST, Phylum=="Firmicutes")$Genus] <- "Firmicutes"
nasal.df$phylum[nasal.df$genus %in% filter(taxa_ST, Phylum=="Proteobacteria")$Genus] <- "Proteobacteria"
nasal.df$phylum[nasal.df$genus %in% filter(taxa_ST, Phylum=="Actinobacteria")$Genus] <- "Actinobacteria"

nasal.df$B_class <- "Xother"
nasal.df$B_class[nasal.df$genus %in% filter(taxa_ST, Class == "Clostridia")$Genus] <- "Clostridia"

tax.ns <- prune_taxa(taxa_sums(physeqGenus_NS) > 0, physeqGenus_NS) %>% tax_table () %>% data.frame()
tax.ns.total <- table(nasal.df$phylum) %>% as.data.frame()

all.taxa.ns <- table(tax.ns$Phylum) %>% sort(decreasing = T) %>% data.frame()
sum(all.taxa.ns$Freq) - sum((all.taxa.ns$Freq)[1:4])

all.taxa.ns
tax.ns.total

tax.ns.total$total <- c(224,121,252,165,376)

colnames(tax.ns.total) <- c("phylum", "Cytokine_Associated", "Total")

#write.csv(file = "../../../../human_microbiome_project/Figures/Figure4/figure4d/nasal.taxa.cytokine.csv",tax.ns.total)

cytokine.df <- rbind(stool.df, skin.df, oral.df, nasal.df)

head(cytokine.df)
dim(cytokine.df)
table(stool.df$significant)
colnames(stool.df)

summary(abs(stool.df$b))
summary(abs(skin.df$b))
summary(abs(oral.df$b))
summary(abs(nasal.df$b))

table(cytokine.df$cytokine) %>% sort()
table(stool.df$cytokine) %>% sort()
table(skin.df$cytokine) %>% sort()
table(oral.df$cytokine) %>% sort()
table(nasal.df$cytokine) %>% sort()

filter(stool.df, cytokine=="VCAM1")

p.density <- ggplot(cytokine.df, aes(x=abs(b),y=class, fill=class)) + geom_density_ridges2() + geom_vline(xintercept = 0.025,linetype="dotted")# +  geom_vline(xintercept = (-0.01),linetype="dotted")
p.density <- p.density + xlim(0, 0.1) + scale_fill_manual(values = body_site_color) + base_theme
p.density <- p.density + xlab("Absolute Value of Beta Coefficient") + ylab("Density")
p.density

cytokine.df$direction <- "Negative"
cytokine.df$direction[cytokine.df$b > 0] <- "Positive"
cytokine.df$bodysite.group <- paste(cytokine.df$class, cytokine.df$direction, sep = "_")
table(cytokine.df$bodysite.group)
p.density <- ggplot(filter(cytokine.df), aes(x=abs(b),y=bodysite.group)) + geom_boxplot(outlier.alpha = 0) + geom_jitter(aes(color=class), size=0.5)
p.density <- p.density + xlim(0, 0.1) + scale_color_manual(values = body_site_color) + base_theme
p.density <- p.density + xlab("Absolute Value of Beta Coefficient") + ylab("Group")
p.density

#ggsave(filename = "../../../../human_microbiome_project/Figures/Figure4/Figure4withbodysite.pdf", p.density,height = 3.5, width = 5, dpi=300)

#ggsave(filename = "../../../../human_microbiome_project/Figures/Figure4/Figure4C.pdf", height = 3.5, width = 5, dpi=300)

# this is to test the difference between positive and negative beta in each bodysite
wilcox.test(filter(cytokine.df, bodysite.group=="Stool_Positive")$b, -filter(cytokine.df, bodysite.group=="Stool_Negative")$b)
wilcox.test(filter(cytokine.df, bodysite.group=="Skin_Positive")$b, -filter(cytokine.df, bodysite.group=="Skin_Negative")$b)
wilcox.test(filter(cytokine.df, bodysite.group=="Oral_Positive")$b, -filter(cytokine.df, bodysite.group=="Oral_Negative")$b)
wilcox.test(filter(cytokine.df, bodysite.group=="Nasal_Positive")$b, -filter(cytokine.df, bodysite.group=="Nasal_Negative")$b)

################
tax.st.total$bodysite <- "Stool"
tax.sk.total$bodysite <- "Skin"
tax.or.total$bodysite <- "Oral"
tax.ns.total$bodysite <- "Nasal"
tax.total <- rbind(tax.st.total,tax.sk.total,tax.or.total,tax.ns.total)
tax.total.wide <- melt(tax.total, id=c("phylum", "bodysite"))
tax.total.wide$phylum <- as.character(tax.total.wide$phylum)
tax.total.wide$phylum[tax.total.wide$phylum == "Xother"] <- "Other"
tax.total.wide$phylum <- factor(tax.total.wide$phylum, levels = c("Actinobacteria", "Bacteroidetes","Firmicutes","Proteobacteria", "Others"))
tax.total.wide$bodysite <- factor(tax.total.wide$bodysite, levels = c("Stool", "Skin","Oral","Nasal"))
tax.total.wide$variable <- factor(tax.total.wide$variable, levels = c( "Total", "Cytokine_Associated"))

p.percent <- ggplot(tax.total.wide, aes(x=variable, y= value, fill=phylum)) + geom_bar(position="fill", stat="identity")
p.percent <- p.percent + facet_grid(.~bodysite) + base_theme + scale_fill_manual(values = phylum_color) + ylab("Percentage (%)") + scale_y_continuous(labels = scales::percent_format(accuracy = 1))
p.percent

################

#test based on Phylum
st.pr.mean <- colMeans(ST.Pr) %>% data.frame() %>% filter(. > 0)
sk.pr.mean <- colMeans(SK.Pr) %>% data.frame() %>% filter(. > 0)
or.pr.mean <- colMeans(OR.Pr) %>% data.frame() %>% filter(. > 0)
ns.pr.mean <- colMeans(NS.Pr) %>% data.frame() %>% filter(. > 0)

st.pr.mean$bodysite <- "Stool"
sk.pr.mean$bodysite <- "Skin"
or.pr.mean$bodysite <- "Oral"
ns.pr.mean$bodysite <- "Nasal"

st.pr.mean$genus <- row.names(st.pr.mean)
sk.pr.mean$genus <- row.names(sk.pr.mean)
or.pr.mean$genus <- row.names(or.pr.mean)
ns.pr.mean$genus <- row.names(ns.pr.mean)

stool.df.pr <- merge(st.pr.mean, stool.df, by="genus", all.x = T)
skinl.df.pr <- merge(sk.pr.mean, skin.df, by="genus", all.x = T)
oral.df.pr <- merge(or.pr.mean, oral.df, by="genus", all.x = T)
nasal.df.pr <- merge(ns.pr.mean, nasal.df, by = "genus", all.x = T)

colnames(stool.df.pr)[2] <- "prevalence"
colnames(skinl.df.pr)[2] <- "prevalence"
colnames(oral.df.pr)[2] <- "prevalence"
colnames(nasal.df.pr)[2] <- "prevalence"

all.df.pr <- rbind(stool.df.pr,skinl.df.pr,oral.df.pr,nasal.df.pr)
all.df.pr$bodysite <- factor(all.df.pr$bodysite, levels = c("Nasal", "Oral", "Skin", "Stool"))

all.df.clean <- all.df.pr[!is.na(all.df.pr$b),]

p.pre.den <- ggplot(all.df.pr, aes(y=bodysite, x=prevalence, fill=bodysite)) + geom_density_ridges() + geom_density_ridges(data=all.df.pr[!is.na(all.df.pr$b),],aes(y=bodysite, x=prevalence,fill=bodysite), alpha=0.4)
p.pre.den <- p.pre.den + geom_vline(xintercept = 0.2, linetype="dashed") +  geom_vline(xintercept = 0.8, linetype="dashed") + base_theme + scale_fill_manual(values=body_site_color)
p.pre.den
#ggsave(filename = "../../../../human_microbiome_project/Figures/Figure4/Prevalance_Cytokine.pdf",p.pre.den,height = 3.5, width = 5, dpi=300)

p.pre.den2 <- ggplot(all.df.pr, aes(y=b, x=prevalence, color=bodysite)) + geom_point()
p.pre.den2 <- p.pre.den2 + geom_vline(xintercept = 0.2, linetype="dashed") +  geom_vline(xintercept = 0.8, linetype="dashed") + base_theme + scale_color_manual(values=body_site_color)
p.pre.den2 <- p.pre.den2 + ggtitle("low-prevalance bacteria changes more dramatically with cytokines")
p.pre.den2
#ggsave(filename = "../../../../human_microbiome_project/Figures/Figure4/Prevalance_Cytokine_scatter.pdf",p.pre.den2,height = 3.5, width = 5, dpi=300)

p.taxa <- ggplot(data=all.df.pr[!is.na(all.df.pr$b),], aes(x=phylum, y=abs(b))) + geom_jitter() + geom_boxplot()
p.taxa <- p.taxa  + facet_grid(bodysite~.) + coord_trans(y="log2")
p.taxa

all.df.clean$bodysite <- factor(all.df.clean$bodysite, levels=c("Stool", "Skin", "Oral", "Nasal"))
p.taxa2 <- ggplot(all.df.clean, aes(x=bodysite, y=b)) + geom_jitter(aes(color=bodysite), size=0.5)# + geom_boxplot(outlier.alpha = 0, alpha=0.1)
p.taxa2 <- p.taxa2 + ylim(0,0.1)  + facet_grid(~phylum, scales = "free") + base_theme + scale_color_manual(values = body_site_color)
p.taxa2 <- p.taxa2 + ylab("Absolute Correlation Coefficients") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p.taxa2

p.taxa2.1 <- ggplot(all.df.clean, aes(x=bodysite, y=b)) + geom_jitter(aes(color=bodysite), size=0.5) #+ geom_boxplot(outlier.alpha = 0, alpha=0.1)
p.taxa2.1 <- p.taxa2.1 + ylim(-0.1,0)  + facet_grid(~phylum, scales = "free") + base_theme + scale_color_manual(values = body_site_color)
p.taxa2.1 <- p.taxa2.1 + ylab("Absolute Correlation Coefficients") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p.taxa2.1

p.taxa2/p.taxa2.1

#ggsave(filename = "../../../../human_microbiome_project/Figures/Figure4/bytaxa.cytokine.pdf", p.taxa2, height = 3.5, width = 10, dpi=300)
#ggsave(filename = "../../../../human_microbiome_project/Figures/Figure4/bytaxa.cytokine.2.pdf", p.taxa2/p.taxa2.1, height = 7, width = 10, dpi=300)

#write.csv(file = "../../../../human_microbiome_project/Figures/Figure4/phylum.beta.csv",all.df.clean)

all.df.clean$Pre <- "Middle"
all.df.clean$Pre[all.df.clean$prevalence > 0.8] <- "Core"
all.df.clean$Pre[all.df.clean$prevalence < 0.2] <- "Oppor"

p.taxa4 <- ggplot(all.df.clean, aes(x=Pre, y=abs(b))) + geom_jitter(aes(color=bodysite), size=0.5) + geom_boxplot(outlier.alpha = 0, alpha=0.1)
p.taxa4 <- p.taxa4 + ylim(0,0.1)  + facet_grid(.~bodysite, scales = "free") + base_theme + scale_color_manual(values = body_site_color) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p.taxa4 <- p.taxa4 + ylab("Correlation Coefficient")
p.taxa4
#ggsave(filename = "../../../../human_microbiome_project/Figures/Figure4/suppli_bytaxa.prev.cytokine2.pdf", p.taxa4, height = 3.5, width = 5, dpi=300)

all.df.clean$Pre <- factor(all.df.clean$Pre, levels = c("Oppor", "Middle", "Core"))
p.taxa3 <- ggplot(all.df.clean, aes(y=Pre, x=abs(b), fill=Pre)) + geom_density_ridges()
p.taxa3 <- p.taxa3 + xlim(0,0.1)  + facet_grid(.~bodysite, scales = "free") + base_theme + scale_fill_manual(values = c("#E64B35FF", "#48BBD5FF", "#91D1C2FF"))
p.taxa3 <- p.taxa3 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p.taxa3
#ggsave(filename = "../../../../human_microbiome_project/Figures/Figure4/Figure4D.pdf", p.taxa3, height = 3, width = 5, dpi=300)

#test core microbiome and oppor on their beta
wilcox.test(abs(filter(all.df.clean, bodysite == "Stool" & prevalence > 0.8)$b),abs(filter(all.df.clean, bodysite == "Stool" & prevalence < 0.2)$b))
wilcox.test(abs(filter(all.df.clean, bodysite == "Skin" & prevalence > 0.8)$b),abs(filter(all.df.clean, bodysite == "Skin" & prevalence < 0.2)$b))
wilcox.test(abs(filter(all.df.clean, bodysite == "Oral" & prevalence > 0.8)$b),abs(filter(all.df.clean, bodysite == "Oral" & prevalence < 0.2)$b))
wilcox.test(abs(filter(all.df.clean, bodysite == "Nasal" & prevalence > 0.8)$b),abs(filter(all.df.clean, bodysite == "Nasal" & prevalence < 0.2)$b))


table(all.df.clean$phylum)
#stool
wilcox.test(abs(filter(all.df.clean, bodysite == "Stool" & phylum == "Actinobacteria" & b>0)$b),
            abs(filter(all.df.clean, bodysite == "Stool" & phylum == "Actinobacteria" & b<0)$b))

wilcox.test(abs(filter(all.df.clean, bodysite == "Stool" & phylum == "Bacteroidetes" & b>0)$b),
            abs(filter(all.df.clean, bodysite == "Stool" & phylum == "Bacteroidetes" & b<0)$b))

wilcox.test(abs(filter(all.df.clean, bodysite == "Stool" & phylum == "Firmicutes" & b>0)$b),
            abs(filter(all.df.clean, bodysite == "Stool" & phylum == "Firmicutes" & b<0)$b))

wilcox.test(abs(filter(all.df.clean, bodysite == "Stool" & phylum == "Proteobacteria" & b>0)$b),
            abs(filter(all.df.clean, bodysite == "Stool" & phylum == "Proteobacteria" & b<0)$b))

wilcox.test(abs(filter(all.df.clean, bodysite == "Stool" & phylum == "Xother" & b>0)$b),
            abs(filter(all.df.clean, bodysite == "Stool" & phylum == "Xother" & b<0)$b))
#skin
wilcox.test(abs(filter(all.df.clean, bodysite == "Skin" & phylum == "Actinobacteria" & b>0)$b),
            abs(filter(all.df.clean, bodysite == "Skin" & phylum == "Actinobacteria" & b<0)$b))

wilcox.test(abs(filter(all.df.clean, bodysite == "Skin" & phylum == "Bacteroidetes" & b>0)$b),
            abs(filter(all.df.clean, bodysite == "Skin" & phylum == "Bacteroidetes" & b<0)$b))

wilcox.test(abs(filter(all.df.clean, bodysite == "Skin" & phylum == "Firmicutes" & b>0)$b),
            abs(filter(all.df.clean, bodysite == "Skin" & phylum == "Firmicutes" & b<0)$b))

wilcox.test(abs(filter(all.df.clean, bodysite == "Skin" & phylum == "Proteobacteria" & b>0)$b),
            abs(filter(all.df.clean, bodysite == "Skin" & phylum == "Proteobacteria" & b<0)$b))

wilcox.test(abs(filter(all.df.clean, bodysite == "Skin" & phylum == "Xother" & b>0)$b),
            abs(filter(all.df.clean, bodysite == "Skin" & phylum == "Xother" & b<0)$b))
#oral
wilcox.test(abs(filter(all.df.clean, bodysite == "Oral" & phylum == "Actinobacteria" & b>0)$b),
            abs(filter(all.df.clean, bodysite == "Oral" & phylum == "Actinobacteria" & b<0)$b))

wilcox.test(abs(filter(all.df.clean, bodysite == "Oral" & phylum == "Bacteroidetes" & b>0)$b),
            abs(filter(all.df.clean, bodysite == "Oral" & phylum == "Bacteroidetes" & b<0)$b))

wilcox.test(abs(filter(all.df.clean, bodysite == "Oral" & phylum == "Firmicutes" & b>0)$b),
            abs(filter(all.df.clean, bodysite == "Oral" & phylum == "Firmicutes" & b<0)$b))

wilcox.test(abs(filter(all.df.clean, bodysite == "Oral" & phylum == "Proteobacteria" & b>0)$b),
            abs(filter(all.df.clean, bodysite == "Oral" & phylum == "Proteobacteria" & b<0)$b))

wilcox.test(abs(filter(all.df.clean, bodysite == "Oral" & phylum == "Xother" & b>0)$b),
            abs(filter(all.df.clean, bodysite == "Oral" & phylum == "Xother" & b<0)$b))
#nasal
wilcox.test(abs(filter(all.df.clean, bodysite == "Nasal" & phylum == "Actinobacteria" & b>0)$b),
            abs(filter(all.df.clean, bodysite == "Nasal" & phylum == "Actinobacteria" & b<0)$b))

wilcox.test(abs(filter(all.df.clean, bodysite == "Nasal" & phylum == "Bacteroidetes" & b>0)$b),
            abs(filter(all.df.clean, bodysite == "Nasal" & phylum == "Bacteroidetes" & b<0)$b))

wilcox.test(abs(filter(all.df.clean, bodysite == "Nasal" & phylum == "Firmicutes" & b>0)$b),
            abs(filter(all.df.clean, bodysite == "Nasal" & phylum == "Firmicutes" & b<0)$b))

wilcox.test(abs(filter(all.df.clean, bodysite == "Nasal" & phylum == "Proteobacteria" & b>0)$b),
            abs(filter(all.df.clean, bodysite == "Nasal" & phylum == "Proteobacteria" & b<0)$b))

wilcox.test(abs(filter(all.df.clean, bodysite == "Nasal" & phylum == "Xother" & b>0)$b),
            abs(filter(all.df.clean, bodysite == "Nasal" & phylum == "Xother" & b<0)$b))



###############
#Chiq Square
###############
tax.st.total$bodysite <- "Stool"
tax.sk.total$bodysite <- "Skin"
tax.or.total$bodysite <- "Oral"
tax.ns.total$bodysite <- "Nasal"
tax.total <- rbind(tax.st.total,tax.sk.total,tax.or.total,tax.ns.total)
tax.total.wide <- melt(tax.total, id=c("phylum", "bodysite"))
tax.total.wide$phylum <- as.character(tax.total.wide$phylum)
tax.total.wide$phylum[tax.total.wide$phylum == "Xother"] <- "Other"
tax.total.wide$phylum <- factor(tax.total.wide$phylum, levels = c("Actinobacteria", "Bacteroidetes","Firmicutes","Proteobacteria", "Others"))
tax.total.wide$bodysite <- factor(tax.total.wide$bodysite, levels = c("Stool", "Skin","Oral","Nasal"))
tax.total.wide$variable <- factor(tax.total.wide$variable, levels = c( "Total", "Cytokine_Associated"))

p.percent <- ggplot(tax.total.wide, aes(x=variable, y= value, fill=phylum)) + geom_bar(position="fill", stat="identity")
p.percent <- p.percent + facet_grid(.~bodysite) + base_theme + scale_fill_manual(values = phylum_color) + ylab("Percentage (%)") + scale_y_continuous(labels = scales::percent_format(accuracy = 1))
p.percent

sum(all.taxa.st$Freq)
sum(all.taxa.sk$Freq)
sum(all.taxa.or$Freq)
sum(all.taxa.ns$Freq)

#Chiq for firmicutes
tax.st.total
st.chiq.firmicutes <- chisq.test(cbind(c(323, 154), c(193, 174)))
st.chiq.firmicutes

tax.sk.total
sk.chiq.firmicutes <- chisq.test(cbind(c(77, sum(tax.sk.total$Cytokine_Associated) - 77), c(262, sum(tax.sk.total$Total) - 262)))
sk.chiq.firmicutes

tax.or.total
or.chiq.firmicutes <- chisq.test(cbind(c(153, sum(tax.or.total$Cytokine_Associated) - 153), c(114, sum(tax.or.total$Total) - 114)))
or.chiq.firmicutes

tax.ns.total
ns.chiq.firmicutes <- chisq.test(cbind(c(89, sum(tax.ns.total$Cytokine_Associated) - 89), c(252, sum(tax.ns.total$Total) - 252)))
ns.chiq.firmicutes

filter(all.df.clean, B_class == "Clostridia")

####test if clostridia in stool is more immunogenic in firmicutes
#stool
table(filter(stool.df, phylum == "Firmicutes")$B_class)
table(filter(tax.st,Phylum == "Firmicutes")$Class)
sum(table(filter(tax.st,Phylum == "Firmicutes")$Class))
sum(table(filter(tax.st,Phylum == "Firmicutes")$Class)) - 129
chisq.test(cbind(c(220, 103), c(129, 64)))

#skin
table(filter(skin.df, phylum == "Firmicutes")$B_class)
table(filter(tax.sk,Phylum == "Firmicutes")$Class)
sum(table(filter(tax.sk,Phylum == "Firmicutes")$Class))
sum(table(filter(tax.sk,Phylum == "Firmicutes")$Class)) - 119
chisq.test(cbind(c(17, 60), c(119, 143)))

#oral
table(filter(oral.df, phylum == "Firmicutes")$B_class)
table(filter(tax.or,Phylum == "Firmicutes")$Class)
sum(table(filter(tax.or,Phylum == "Firmicutes")$Class))
sum(table(filter(tax.or,Phylum == "Firmicutes")$Class)) - 53
chisq.test(cbind(c(78, 75), c(53, 61)))

#nasal
table(filter(nasal.df, phylum == "Firmicutes")$B_class)
table(filter(tax.ns,Phylum == "Firmicutes")$Class)
sum(table(filter(tax.ns,Phylum == "Firmicutes")$Class))
sum(table(filter(tax.ns,Phylum == "Firmicutes")$Class)) - 121
chisq.test(cbind(c(31, 58), c(121, 131)))

# make heatmap for who is immunogenic
# library(ComplexHeatmap)
# library(circlize)
# library(dendextend)
# wide.cytokine.count.df <- table(cytokine.df$cytokine, cytokine.df$genus) %>% data.frame() %>% dcast(Var1 ~ Var2, value.var =  "Freq")
# dim(wide.cytokine.count.df)
# rownames(wide.cytokine.count.df) <- wide.cytokine.count.df$Var1
# wide.cytokine.count.df <- select(wide.cytokine.count.df, -Var1)
# 
# col_fun = colorRamp2(c(0, 3), c("white", "red"))
# row_dend = as.dendrogram(hclust(dist(wide.cytokine.count.df)))
# row_dend = color_branches(row_dend, k = 6)
# 
# heatmap.m.c<- Heatmap(wide.cytokine.count.df, col = col_fun,
#                          column_title = "Stool_Richness:Cytokine",
#                          column_km = 7,
#                          column_dend_side = "top",
#                          cluster_rows = row_dend,
#                          row_dend_reorder = TRUE,
#                          row_names_gp=grid::gpar(fontsize = 9,fontface = "bold"),
#                          column_names_gp = grid::gpar(fontsize = 7),
#                          column_names_rot = 45)
# heatmap.m.c

#check the most prevelant cytokine
table(cytokine.df$cytokine) %>% sort()

table(stool.df$cytokine) %>% sort() 
table(skin.df$cytokine) %>% sort()
table(oral.df$cytokine) %>% sort()
table(nasal.df$cytokine) %>% sort()

supplimentary.cytokine.table <- table(cytokine.df$cytokine, cytokine.df$bodysite.group) %>% data.frame()
#write.csv(file = "../../../../human_microbiome_project/Supplementary_data/Cytokine_number_withbodysites.csv",supplimentary.cytokine.table)

st.pre.cyto <- table(stool.df$cytokine) %>% data.frame()
sk.pre.cyto <- table(skin.df$cytokine) %>% data.frame()
or.pre.cyto <- table(oral.df$cytokine) %>% data.frame()
ns.pre.cyto <- table(nasal.df$cytokine) %>% data.frame()

#######which cytokine is the most potent
stool.df.bycytokine <- stool.df %>% group_by(cytokine) %>% mutate(mean_beta = mean(abs(b))) %>% select(cytokine,mean_beta,class) %>% unique()
skin.df.bycytokine <- skin.df %>% group_by(cytokine) %>% mutate(mean_beta = mean(abs(b))) %>% select(cytokine,mean_beta,class) %>% unique()
oral.df.bycytokine <- oral.df %>% group_by(cytokine) %>% mutate(mean_beta = mean(abs(b))) %>% select(cytokine,mean_beta,class) %>% unique()
nasal.df.bycytokine <- nasal.df %>% group_by(cytokine) %>% mutate(mean_beta = mean(abs(b))) %>% select(cytokine,mean_beta,class) %>% unique()

heatmap.df <- left_join(stool.df.bycytokine, skin.df.bycytokine, by="cytokine") %>% left_join(oral.df.bycytokine,by="cytokine") %>% left_join(nasal.df.bycytokine,by="cytokine") %>% 
  select(-class.x, -class.y, -class.x.x, -class.y.y) %>% rename(Stool = mean_beta.x) %>% rename(Skin = mean_beta.y) %>% rename(Oral = mean_beta.x.x) %>% rename(Nasal = mean_beta.y.y) %>% data.frame()
heatmap.df

library(ComplexHeatmap)
library(circlize)
library(dendextend)
library(dendsort)

rownames(heatmap.df) <- heatmap.df$cytokine
heatmap.df <- select(heatmap.df, -cytokine)

max(na.omit(heatmap.df))
min(na.omit(heatmap.df))

col_fun = colorRamp2(c(0, 0.1), c("white", "red"))

row_dend = as.dendrogram(hclust(dist(heatmap.df)))
row_dend = color_branches(row_dend, k = 1)

heatmap.cytokinecount <- Heatmap(as.matrix(heatmap.df), 
                                       col = col_fun,
                                       column_title = "Cytokine mean",
                                       column_km = 1,
                                       column_dend_side = "top",
                                       cluster_rows = row_dend,
                                       row_dend_reorder = TRUE,
                                       row_names_gp=grid::gpar(fontsize = 9,fontface = "bold"),
                                       column_names_gp = grid::gpar(fontsize = 7))

heatmap.cytokinecount
#to save it, save_as_pdf 3 * 12 inch


#https://jokergoo.github.io/2020/05/21/make-circular-heatmaps/

#save as pdf, on 5 * 5 inch for reproduce
circos.par(gap.after =  10)
circos.heatmap(heatmap.df, 
               col = col_fun, 
               rownames.side = "inside")

circos.clear()

circos.par(gap.after =  10)
circos.heatmap(heatmap.df, 
               col = col_fun, 
               dend.track.height = 0.7,
               dend.side = "inside")

circos.clear()



############how does cytokine sharp core microbiome
#Actinobacteria
p.act.pre <- filter(all.df.clean, phylum=="Actinobacteria") %>% filter(bodysite == "Stool") %>% 
  ggbetweenstats(x=Pre, y=b, title = "Actinobacteria from Stool Microbiome/cytokine correlation by Prevalence")
p.act.pre 
ggsave(filename = "./Suppl.figure/Figure4C_associated/ST_ACT.pdf", p.act.pre, height = 5, width = 7, dpi=300)

#Bacteroidetes
p.bac.pre <- filter(all.df.clean, phylum=="Bacteroidetes") %>% filter(bodysite == "Stool") %>% 
  ggbetweenstats(x=Pre, y=b, title = "Bacteroidetes from Stool Microbiome/cytokine correlation by Prevalence")
p.bac.pre 
ggsave(filename = "./Suppl.figure/Figure4C_associated/ST_BAC.pdf", p.bac.pre, height = 5, width = 7, dpi=300)

#Firmicutes
p.fir.pre <- filter(all.df.clean, phylum=="Firmicutes") %>% filter(bodysite == "Stool") %>% 
  ggbetweenstats(x=Pre, y=b, title = "Firmicutes from Stool Microbiome/cytokine correlation by Prevalence")
p.fir.pre 
ggsave(filename = "./Suppl.figure/Figure4C_associated/ST_FIR.pdf", p.fir.pre, height = 5, width = 7, dpi=300)

#Proteobacteria
p.pro.pre <- filter(all.df.clean, phylum=="Proteobacteria") %>% filter(bodysite == "Stool") %>% 
  ggbetweenstats(x=Pre, y=b, title = "Proteobacteria from Stool Microbiome/cytokine correlation by Prevalence")
p.pro.pre 
ggsave(filename = "./Suppl.figure/Figure4C_associated/ST_PRO.pdf", p.pro.pre, height = 5, width = 7, dpi=300)

#######Skin#####################
#Actinobacteria
p.act.pre.sk <- filter(all.df.clean, phylum=="Actinobacteria") %>% filter(bodysite == "Skin") %>% 
  ggbetweenstats(x=Pre, y=b, title = "Actinobacteria from Skin Microbiome/cytokine correlation by Prevalence")
p.act.pre.sk 
ggsave(filename = "./Suppl.figure/Figure4C_associated/SK_ACT.pdf", p.act.pre.sk, height = 5, width = 7, dpi=300)

#Bacteroidetes
p.bac.pre.sk <- filter(all.df.clean, phylum=="Bacteroidetes") %>% filter(bodysite == "Skin") %>% 
  ggbetweenstats(x=Pre, y=b, title = "Bacteroidetes from Skin Microbiome/cytokine correlation by Prevalence")
p.bac.pre.sk 
ggsave(filename = "./Suppl.figure/Figure4C_associated/SK_BAC.pdf", p.bac.pre.sk, height = 5, width = 7, dpi=300)

#Firmicutes
p.fir.pre.sk <- filter(all.df.clean, phylum=="Firmicutes") %>% filter(bodysite == "Skin") %>% 
  ggbetweenstats(x=Pre, y=b, title = "Firmicutes from Skin Microbiome/cytokine correlation by Prevalence")
p.fir.pre.sk
ggsave(filename = "./Suppl.figure/Figure4C_associated/SK_FIR.pdf", p.fir.pre.sk, height = 5, width = 7, dpi=300)

#Proteobacteria
p.pro.pre.sk <- filter(all.df.clean, phylum=="Proteobacteria") %>% filter(bodysite == "Skin") %>% 
  ggbetweenstats(x=Pre, y=b, title = "Proteobacteria from Skin Microbiome/cytokine correlation by Prevalence")
p.pro.pre.sk 
ggsave(filename = "./Suppl.figure/Figure4C_associated/SK_PRO.pdf", p.pro.pre.sk, height = 5, width = 7, dpi=300)

#######Oral#####################
#Actinobacteria
p.act.pre.or <- filter(all.df.clean, phylum=="Actinobacteria") %>% filter(bodysite == "Oral") %>% 
  ggbetweenstats(x=Pre, y=b, title = "Actinobacteria from Oral Microbiome/cytokine correlation by Prevalence")
p.act.pre.or 
ggsave(filename = "./Suppl.figure/Figure4C_associated/OR_ACT.pdf", p.act.pre.or, height = 5, width = 7, dpi=300)

#Bacteroidetes
p.bac.pre.or <- filter(all.df.clean, phylum=="Bacteroidetes") %>% filter(bodysite == "Oral") %>% 
  ggbetweenstats(x=Pre, y=b, title = "Bacteroidetes from Oral Microbiome/cytokine correlation by Prevalence")
p.bac.pre.or 
ggsave(filename = "./Suppl.figure/Figure4C_associated/OR_BAC.pdf", p.bac.pre.or, height = 5, width = 7, dpi=300)

#Firmicutes
p.fir.pre.or <- filter(all.df.clean, phylum=="Firmicutes") %>% filter(bodysite == "Oral") %>% 
  ggbetweenstats(x=Pre, y=b, title = "Firmicutes from Oral Microbiome/cytokine correlation by Prevalence")
p.fir.pre.or
ggsave(filename = "./Suppl.figure/Figure4C_associated/OR_FIR.pdf", p.fir.pre.or, height = 5, width = 7, dpi=300)

#Proteobacteria
p.pro.pre.or <- filter(all.df.clean, phylum=="Proteobacteria") %>% filter(bodysite == "Oral") %>% 
  ggbetweenstats(x=Pre, y=b, title = "Proteobacteria from Oral Microbiome/cytokine correlation by Prevalence")
p.pro.pre.or 
ggsave(filename = "./Suppl.figure/Figure4C_associated/OR_PRO.pdf", p.pro.pre.or, height = 5, width = 7, dpi=300)

#######Nasal#####################
#Actinobacteria
p.act.pre.ns <- filter(all.df.clean, phylum=="Actinobacteria") %>% filter(bodysite == "Nasal") %>% 
  ggbetweenstats(x=Pre, y=b, title = "Actinobacteria from Nasal Microbiome/cytokine correlation by Prevalence")
p.act.pre.ns 
ggsave(filename = "./Suppl.figure/Figure4C_associated/NS_ACT.pdf", p.act.pre.ns, height = 5, width = 7, dpi=300)

#Bacteroidetes
p.bac.pre.ns <- filter(all.df.clean, phylum=="Bacteroidetes") %>% filter(bodysite == "Nasal") %>% 
  ggbetweenstats(x=Pre, y=b, title = "Bacteroidetes from Nasal Microbiome/cytokine correlation by Prevalence")
p.bac.pre.ns 
ggsave(filename = "./Suppl.figure/Figure4C_associated/NS_BAC.pdf", p.bac.pre.ns, height = 5, width = 7, dpi=300)

#Firmicutes
p.fir.pre.ns <- filter(all.df.clean, phylum=="Firmicutes") %>% filter(bodysite == "Nasal") %>% 
  ggbetweenstats(x=Pre, y=b, title = "Firmicutes from Nasal Microbiome/cytokine correlation by Prevalence")
p.fir.pre.ns
ggsave(filename = "./Suppl.figure/Figure4C_associated/NS_FIR.pdf", p.fir.pre.ns, height = 5, width = 7, dpi=300)

#Proteobacteria
p.pro.pre.ns <- filter(all.df.clean, phylum=="Proteobacteria") %>% filter(bodysite == "Nasal") %>% 
  ggbetweenstats(x=Pre, y=b, title = "Proteobacteria from Nasal Microbiome/cytokine correlation by Prevalence")
p.pro.pre.ns 
ggsave(filename = "./Suppl.figure/Figure4C_associated/NS_PRO.pdf", p.pro.pre.ns, height = 5, width = 7, dpi=300)

