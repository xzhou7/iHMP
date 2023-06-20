#Plot strain replacement
library(x)
library(phyloseq)
library(dplyr)
library(patchwork)
library(ggplot2)
library(ggrepel)
library(ggridges)
library(ggstatsplot)
library(ggpubr)

base_theme = theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 13),
        panel.grid.minor = element_blank())
body_site_color = c(
  "Stool" = ggsci::pal_jama()(n=7)[2],
  "Skin" = ggsci::pal_jama()(n=7)[3],
  "Oral" = ggsci::pal_jama()(n=7)[4],
  "Nasal" = ggsci::pal_jama()(n=7)[5])

setwd("~/Library/CloudStorage/Box-Box/human_microbiome_project/Figures/Figure3/detailed/strainreplacement/")
load("~/Library/CloudStorage/Box-Box/XinZhouFiles/Projects/ZXP2_HMP2_Core/Analysis/Robject/Meta_MultiOmes0413.RData")
load("~/Library/CloudStorage/Box-Box/XinZhouFiles/Projects/ZXP2_HMP2_Core/Analysis/Robject/DetailedPhyloseq.RData")
load("~/Library/CloudStorage/Box-Box/XinZhouFiles/Projects/ZXP2_HMP2_Core/Analysis/Robject/Richness_Bygenus.RData")
load("~/Library/CloudStorage/Box-Box/XinZhouFiles/Projects/ZXP2_HMP2_Core/Analysis/Robject/Prevalance.RData")

load("../../../../data_analysis/microbiome_stability/stool_stability")

Stool.replacement <- read.csv("./stool.replacement_new.csv", header = T, row.names = 1)
Skin.replacement <- read.csv("./skin.replacement_new.csv", header = T, row.names = 1)
Skin.replacement <- filter(Skin.replacement, taxa != "Unclassified_Actinobacteria.1")
Oral.replacement <- read.csv("./oral.replacement_new.csv", header = T, row.names = 1)
Oral.replacement <- filter(Oral.replacement, taxa != "Unclassified_Actinobacteria.1")
Nasal.replacement <- read.csv("./nasal.replacement_new.csv", header = T, row.names = 1)
Nasal.replacement <- filter(Nasal.replacement, taxa != "Unclassified_Actinobacteria.1")

taxa_ST <- data.frame(tax_table(physeq_ST))
taxa_SK <- data.frame(tax_table(physeq_SK))

summary(Stool.replacement$dist)
summary(Skin.replacement$dist)
summary(Oral.replacement$dist)
summary(Nasal.replacement$dist)

dim(Stool.replacement)
dim(Skin.replacement)
dim(Oral.replacement)
dim(Nasal.replacement)

replacement.df <- rbind(Stool.replacement,Skin.replacement,Oral.replacement,Nasal.replacement)
replacement.df$Jdist <- 1 - replacement.df$jaccard
replacement.df$Interval <- "A.0-3month"
replacement.df$Interval[replacement.df$timediff > 120 & replacement.df$timediff < 360 | replacement.df$timediff == 120] <- "B.3-12month"
replacement.df$Interval[replacement.df$timediff > 360| replacement.df$timediff == 360] <- "C.12monthabove"

table(replacement.df$Interval) 

p.bc.jd <- ggscatterstats(replacement.df, x=dist, y=jaccard, title = "Bray Curtis Distance with Jaccard Distance of sample Pairs")
p.bc.jd
#ggsave(filename = "./Bray.Jaccard.pdf", p.bc.jd, width = 6, height = 5, dpi = 300)

density(replacement.df$asv_number)
table(replacement.df$taxa) %>% sort()
# 
# replacement.df %>% group_by(taxa,Bodysite, IRIS) %>% mutate(dist_mean = mean(jaccard)) %>% 
#   select(taxa,Bodysite,IRIS, dist_mean) %>% unique() %>% filter (Bodysite== "Nasal") %>% 
#   ggbetweenstats(x=IRIS, y=dist_mean)

#ggplot(replacement.df, aes(x=asv_number, y=dist)) + geom_point()

##############################here we need to set up distance matrix, either mean.dist=mean(dist) for BC distance, or mean.dist=mean(jaccard) for JSD
short.replacement.plot <- replacement.df %>% filter(Interval=="A.0-3month") %>% group_by(Bodysite,taxa) %>% mutate(mean.dist = mean(jaccard)) %>% select(taxa, Bodysite, mean.dist) %>% unique()
middle.replacement.plot <- replacement.df %>% filter(Interval=="B.3-12month") %>% group_by(Bodysite,taxa) %>% mutate(mean.dist = mean(jaccard)) %>% select(taxa, Bodysite, mean.dist) %>% unique()
long.replacement.plot <- replacement.df %>% filter(Interval=="C.12monthabove") %>% group_by(Bodysite,taxa) %>% mutate(mean.dist = mean(jaccard)) %>% select(taxa, Bodysite, mean.dist) %>% unique()

replacement.df.summary <- replacement.df %>% group_by(Interval,Bodysite,taxa) %>% mutate(mean.dist = mean(dist))%>% mutate(mean.jd = mean(Jdist)) %>% select(taxa, Interval,Bodysite, mean.dist,mean.jd) %>% unique()
replacement.df.summary

p.summary.bc <- ggplot(replacement.df.summary, aes(x=Bodysite, y=mean.jd)) + 
  geom_jitter(aes(color=Bodysite), size=0.5) +
  geom_boxplot(outlier.shape=NA, alpha=0.01) + 
  scale_color_manual(values=body_site_color) +
  base_theme + ggtitle("summary")

p.summary.bc

p.density.replacement <- ggplot(replacement.df.summary, aes(x = mean.jd, y = Bodysite,fill=Bodysite)) + geom_density_ridges2() + scale_fill_manual(values = body_site_color)
p.density.replacement

replacement.df.summary$Bodysite <- factor(replacement.df.summary$Bodysite, levels = c("Stool", "Skin", "Oral","Nasal"))
p.hist.replacement <- grouped_gghistostats(data = replacement.df.summary, 
                                           x = mean.jd,
                                           grouping.var = Bodysite, 
                                           normal.curve = TRUE,
                                           normal.curve.args = list(color = "red", size = 1),
                                           plotgrid.args  = list(nrow = 2))

p.hist.replacement
#ggsave(filename = "./p.hist.replacement.pdf", p.hist.replacement, width = 14, height = 10, dpi = 300 )

table(replacement.df.summary$mean.jd == 0,replacement.df.summary$Bodysite)
table(replacement.df.summary$mean.jd == 1,replacement.df.summary$Bodysite)

p.summary.jd <- ggplot(replacement.df.summary, aes(x=Bodysite, y=mean.jd)) + 
  geom_jitter(aes(color=Bodysite), size=0.5) +
  geom_boxplot(outlier.shape=NA, alpha=0.01) + 
  scale_color_manual(values=body_site_color) +
  base_theme + ggtitle("summary")

p.summary.jd

ggplot(replacement.df.summary, aes(x = mean.jd, y = Bodysite,fill=Bodysite)) + geom_density_ridges2() + scale_fill_manual(values = body_site_color) + base_theme

replacement.df.summary
# 
# replacement.3ormore <- table(replacement.df$taxa) %>% sort() %>% data.frame() %>% filter(Freq > 3)
# replacement.3ormore$Var1 <- as.character(replacement.3ormore$Var1)
# 
# p.time.dff1 <- replacement.df %>% filter(asv_number > 2) %>% filter(taxa %in% replacement.3ormore$Var1) %>% ggplot(aes(x=timediff, y=dist)) + geom_point(size=0.2) + facet_wrap(.~Bodysite) + geom_smooth() + base_theme
# p.time.dff1 <- p.time.dff1 + ggtitle("strain replacement rate not correlate with time") + xlab("time_interval") + ylab("BC Dissimilarity")
# p.time.dff1
# 
# p.time.dff <- replacement.df %>% filter(asv_number > 2) %>% filter(taxa %in% replacement.3ormore$Var1) %>% ggplot(aes(x=timediff, y=jaccard)) + geom_point(size=0.2) + facet_wrap(.~Bodysite) + geom_smooth() + base_theme
# p.time.dff <- p.time.dff + ggtitle("strain replacement rate not correlate with time") + xlab("time_interval") + ylab("Jaccard Similarity")
# p.time.dff
# 
# p.density <- ggplot(replacement.df, aes(x=timediff)) +  geom_density(aes(color=Bodysite)) + 
#   geom_vline(xintercept = 360, linetype="dotted") + geom_vline(xintercept = 120, linetype="dotted") +
#   scale_color_manual(values=body_site_color) + base_theme
# p.density
# 
# 
# p1 <- ggplot(short.replacement.plot, aes(x=Bodysite, y=mean.dist)) + 
#   geom_jitter(aes(color=Bodysite), size=0.5) +
#   geom_boxplot(outlier.shape=NA, alpha=0.01) + 
#   scale_color_manual(values=body_site_color) +
#   base_theme + ggtitle("Within 3 month")
# 
# p2 <- ggplot(middle.replacement.plot, aes(x=Bodysite, y=mean.dist)) + 
#   geom_jitter(aes(color=Bodysite), size=0.5) +
#   geom_boxplot(outlier.shape=NA, alpha=0.01) + 
#   scale_color_manual(values=body_site_color) +
#   base_theme + ggtitle("3~12 Month")
# 
# p3 <- ggplot(long.replacement.plot, aes(x=Bodysite, y=mean.dist)) + 
#   geom_jitter(aes(color=Bodysite), size=0.5) + 
#   geom_boxplot(outlier.shape=NA, alpha=0.01) + 
#   scale_color_manual(values=body_site_color) +
#   base_theme + ggtitle("More than 12 month")
# 
# p1 + p2 + p3 +  plot_layout(guides = 'collect')


short.replacement.plot$time_interval <- "Within 3 month"
middle.replacement.plot$time_interval <- "3~12 Month"
long.replacement.plot$time_interval <- "More than 12 month"

replacement.plot.time <- rbind(short.replacement.plot, middle.replacement.plot,long.replacement.plot)

replacement.plot.time$time_interval <- factor(replacement.plot.time$time_interval,levels = c("Within 3 month","3~12 Month", "More than 12 month"))

Stool.comp <- ggbetweenstats(data = filter(replacement.plot.time, Bodysite == "Stool"),
               x = time_interval, 
               y = mean.dist, 
               p.adjust.method = "BH",
               type = "nonparametric", 
               ggtheme = base_theme,
               xlab = "",
               ylab = "Recurrent Sample Similarity by Mean \n (1-Jaccard Distance)",
               title = "Stool Microbiome replacement by time"
               )

Stool.comp

Skin.comp <- ggbetweenstats(data = filter(replacement.plot.time, Bodysite == "Skin"),
                             x = time_interval, 
                             y = mean.dist, 
                             p.adjust.method = "BH",
                             type = "nonparametric", 
                             ggtheme = base_theme,
                             xlab = "",
                             ylab = "Recurrent Sample Similarity by Mean \n (1-Jaccard Distance)",
                             title = "Skin Microbiome replacement by time"
)

Skin.comp

Oral.comp <- ggbetweenstats(data = filter(replacement.plot.time, Bodysite == "Oral"),
                            x = time_interval, 
                            y = mean.dist, 
                            p.adjust.method = "BH",
                            type = "nonparametric", 
                            ggtheme = base_theme,
                            xlab = "",
                            ylab = "Recurrent Sample Similarity by Mean \n (1-Jaccard Distance)",
                            title = "Oral Microbiome replacement by time"
)

Oral.comp

Nasal.comp <- ggbetweenstats(data = filter(replacement.plot.time, Bodysite == "Nasal"),
                            x = time_interval, 
                            y = mean.dist, 
                            p.adjust.method = "BH",
                            type = "nonparametric", 
                            ggtheme = base_theme,
                            xlab = "",
                            ylab = "Recurrent Sample Similarity by Mean \n (1-Jaccard Distance)",
                            title = "Nasal Microbiome replacement by time"
)

Nasal.comp

replace_by_time <- Stool.comp + Skin.comp + Oral.comp + Nasal.comp
replace_by_time
#ggsave(filename = "./replacement.by.time.pdf", replace_by_time, width = 14, height = 10, dpi = 300)

####################################################################################
#IRIS comparasion
####################################################################################
replacement.df$IRIS <- "Unknown"
replacement.df$IRIS[replacement.df$SubjectID %in% filter(sc,IRIS =="IS")$SubjectID] <- "IS"
replacement.df$IRIS[replacement.df$SubjectID %in% filter(sc,IRIS =="IR")$SubjectID] <- "IR"
IR.df <- filter(replacement.df, IRIS == "IR")
IS.df <- filter(replacement.df, IRIS == "IS")

IR.replacement.plot <- IR.df %>% group_by(Bodysite,taxa) %>% mutate(mean.dist = mean(jaccard)) %>% select(taxa, Bodysite, mean.dist) %>% unique()
IS.replacement.plot <- IS.df %>% group_by(Bodysite,taxa) %>% mutate(mean.dist = mean(jaccard)) %>% select(taxa, Bodysite, mean.dist) %>% unique()

IS.replacement.plot$IRIS <- "IS"
IR.replacement.plot$IRIS <- "IR"

IRIS.replacement.plot <- replacement.df %>% group_by(Bodysite,taxa,IRIS,Interval) %>% mutate(mean.dist = mean(dist))%>%
  mutate(mean.jd = mean(jaccard)) %>% select(taxa, Bodysite, mean.dist,mean.jd,IRIS,Interval) %>% unique()

IRIS.replacement.plot

p1.1 <- ggplot(IS.replacement.plot, aes(x=Bodysite, y=mean.dist)) + 
  geom_jitter(aes(color=Bodysite), size=0.5) +
  geom_boxplot(outlier.shape=NA, alpha=0.01) + 
  scale_color_manual(values=body_site_color) +
  base_theme + ggtitle("insulin sensitive")
p1.1

p2.1 <- ggplot(IR.replacement.plot, aes(x=Bodysite, y=mean.dist)) + 
  geom_jitter(aes(color=Bodysite), size=0.5) + 
  geom_boxplot(outlier.shape=NA, alpha=0.01) + 
  scale_color_manual(values=body_site_color) +
  base_theme + ggtitle("insulin resistant")
p2.1

IS.replacement.plot %>% group_by(Bodysite) %>% mutate(median.jd = median(mean.dist)) %>% summarise(unique(median.jd))
IR.replacement.plot %>% group_by(Bodysite) %>% mutate(median.jd = median(mean.dist)) %>% summarise(unique(median.jd))

replacement.plot.IRIS <- rbind(IS.replacement.plot,IR.replacement.plot) 
table(replacement.plot.IRIS$IRIS)

Stool.comp.IRIS <- ggbetweenstats(data=filter(replacement.plot.IRIS, Bodysite == "Stool"),
                                  x= IRIS, 
                                  y=mean.dist, 
                                  p.adjust.method = "BH",
                                  type = "nonparametric", 
                                  ggtheme = base_theme,
                                  xlab = "",
                                  ylab = "Recurrent Sample Similarity by Mean \n (1-Jaccard Distance)",
                                  title = "Stool Microbiome replacement by IRIS")

Skin.comp.IRIS <- ggbetweenstats(data=filter(replacement.plot.IRIS, Bodysite == "Skin"),
                                  x= IRIS, 
                                  y=mean.dist, 
                                  p.adjust.method = "BH",
                                  type = "nonparametric", 
                                  ggtheme = base_theme,
                                  xlab = "",
                                  ylab = "Recurrent Sample Similarity by Mean \n (1-Jaccard Distance)",
                                  title = "Skin Microbiome replacement by IRIS")

Oral.comp.IRIS <- ggbetweenstats(data=filter(replacement.plot.IRIS, Bodysite == "Oral"),
                                  x= IRIS, 
                                  y=mean.dist, 
                                  p.adjust.method = "BH",
                                  type = "nonparametric", 
                                  ggtheme = base_theme,
                                  xlab = "",
                                  ylab = "Recurrent Sample Similarity by Mean \n (1-Jaccard Distance)",
                                  title = "Oral Microbiome replacement by IRIS")

Nasal.comp.IRIS <- ggbetweenstats(data=filter(replacement.plot.IRIS, Bodysite == "Nasal"),
                                  x= IRIS, 
                                  y=mean.dist, 
                                  p.adjust.method = "BH",
                                  type = "nonparametric", 
                                  ggtheme = base_theme,
                                  xlab = "",
                                  ylab = "Recurrent Sample Similarity by Mean \n (1-Jaccard Distance)",
                                  title = "Nasal Microbiome replacement by IRIS")

replacement.IRIS <- Stool.comp.IRIS + Skin.comp.IRIS + Oral.comp.IRIS + Nasal.comp.IRIS
replacement.IRIS
#ggsave(filename = "./replacement.by.IRIS.pdf", replacement.IRIS, width = 14, height = 10, dpi = 300)

####################################################################################
#by phylogenetic relationship 
####################################################################################
replacement.plot.3 <- replacement.df %>% group_by(Bodysite,taxa) %>% mutate(mean.dist = mean(dist)) %>% 
  mutate(mean.jd = mean(Jdist)) %>% mutate(meanrich = mean(asv_number)) %>% select(taxa, Bodysite, mean.dist,mean.jd,meanrich) %>% unique()
replacement.plot.3

filter(replacement.plot.3, meanrich > 1) %>% group_by(Bodysite) %>% arrange(mean.dist)

replacement.plot.rank <- replacement.plot.3 %>% group_by(Bodysite) %>% arrange(mean.dist)
replacement.plot.rank

replacement.plot.st <- filter(replacement.plot.3, Bodysite=="Stool")

#select(richness_ST_byGenus, Unclassified_Porphyromonadaceae, Elusimicrobium)

is.na(richness_ST_byGenus) <- richness_ST_byGenus == 0
richness.nonezero.st <- colMeans(richness_ST_byGenus, na.rm = T) %>% data.frame
colnames(richness.nonezero.st) <- "mean_richness"
richness.nonezero.st
# 
preva.st <- colMeans(ST.Pr) %>% data.frame 
preva.st
colnames(preva.st) <- "mean_preva"

stool.replacement.richness <- merge(replacement.plot.st,richness.nonezero.st, by.x="taxa", by.y="row.names")
stool.replacement.richness <- filter(stool.replacement.richness,mean_richness > 1)
stool.replacement.r.p <- merge(replacement.plot.st, preva.st, by.x = "taxa", by.y= "row.names")

stool.replacement.r.p <- stool.replacement.r.p %>% arrange(desc(meanrich))

table(taxa_ST$Phylum) %>% data.frame() %>% arrange(desc(Freq))

stool.replacement.r.p

stool.replacement.r.p$phylum <- "Other"
stool.replacement.r.p$phylum[stool.replacement.r.p$taxa %in% filter(taxa_ST, Phylum=="Bacteroidetes")$Genus] <- "Bacteroidetes"
stool.replacement.r.p$phylum[stool.replacement.r.p$taxa %in% filter(taxa_ST, Phylum=="Firmicutes")$Genus] <- "Firmicutes"
#stool.replacement.r.p$phylum[stool.replacement.r.p$taxa %in% filter(taxa_ST, Family=="Lachnospiraceae")$Genus] <- "Fam_Lachnospiraceae"
#stool.replacement.r.p$phylum[stool.replacement.r.p$taxa %in% filter(taxa_ST, Family=="Ruminococcaceae")$Genus] <- "Fam_Ruminococcaceae"
stool.replacement.r.p$phylum[stool.replacement.r.p$taxa %in% filter(taxa_ST, Phylum=="Proteobacteria")$Genus] <- "Proteobacteria"
stool.replacement.r.p$phylum[stool.replacement.r.p$taxa %in% filter(taxa_ST, Phylum=="Actinobacteria")$Genus] <- "Actinobacteria"

#stool.replacement.richness$phylum[stool.replacement.richness$taxa %in% filter(taxa_ST, Phylum=="Verrucomicrobia")$Genus] <- "Verrucomicrobia"

p.diver.pre <- ggplot(filter(stool.replacement.r.p), aes(x=meanrich, y=mean.dist)) + geom_jitter() + geom_smooth( method = "lm") + base_theme
p.diver.pre

p.diver.pre.st <- ggplot(stool.replacement.r.p, aes(x=meanrich, y=mean_preva,color=phylum)) + geom_jitter() + geom_smooth(method = "lm") + facet_wrap(.~phylum, scales = "free") + base_theme
p.diver.pre.st

#ggsave("../../../../XinZhouFiles/Projects/ZXP2_HMP2_Core/Analysis/Suppl.figure/rich_pre_stool.pdf",p.diver.pre, width = 6, height = 4)

stool.replacement.r.p$phylum2 <- "Other"
stool.replacement.r.p$phylum2[stool.replacement.r.p$taxa %in% filter(taxa_ST, Family=="Lachnospiraceae")$Genus] <- "Fam_Lachnospiraceae"
stool.replacement.r.p$phylum2[stool.replacement.r.p$taxa %in% filter(taxa_ST, Family=="Ruminococcaceae")$Genus] <- "Fam_Ruminococcaceae"

scatter_Fir.Stool <- ggscatterstats(data=filter(stool.replacement.r.p, phylum== "Firmicutes"), 
          x = meanrich, y=mean_preva, ylab="Prevalence (%)", xlab = "Observed Richness (Mean)" ,
          title= "Firmicutes/Stool"
            ) +  scale_y_continuous(labels = scales::percent)
scatter_Fir.Stool

ggscatterstats(data=filter(stool.replacement.r.p, phylum== "Actinobacteria"), 
               x = mean.dist, y=mean_preva, ylab="Prev (%)", xlab = "BC distance" ,
               title= "Actinobacteria/Stool") 

ggscatterstats(data=filter(stool.replacement.r.p, phylum== "Bacteroidetes"), 
               x = mean.dist, y=mean_preva, ylab="Richness (%)", xlab = "Jaccard Similarity (Mean)" ,
               title= "Bacteroidetes/Stool") 

ggscatterstats(data=filter(stool.replacement.r.p, phylum== "Firmicutes"), 
               x = mean.dist, y=mean_preva, ylab="Richness (%)", xlab = "Jaccard Similarity (Mean)" ,
               title= "Firmicutes/Stool") 

ggscatterstats(data=filter(stool.replacement.r.p, phylum== "Proteobacteria"), 
               x = mean.dist, y=mean_preva, ylab="Richness (%)", xlab = "Jaccard Similarity (Mean)" ,
               title= "Proteobacteria/Stool") 

ggscatterstats(data=filter(stool.replacement.r.p, phylum== "Other"), 
               x = mean.dist, y=mean_preva, ylab="Richness (%)", xlab = "Jaccard Similarity (Mean)" ,
               title= "Other/Stool") 

###########################################################################
#Skin microbiome
replacement.plot.sk <- filter(replacement.plot.3, Bodysite=="Skin")
replacement.plot.3 %>%  ggplot(aes(x=mean.dist, y=mean.jd)) + geom_point()

is.na(richness_SK_byGenus) <- richness_SK_byGenus == 0
richness.nonezero.sk <- colMeans(richness_SK_byGenus, na.rm = T) %>% data.frame
colnames(richness.nonezero.sk) <- "mean_richness"
richness.nonezero.sk

preva.sk <- colMeans(SK.Pr)%>% data.frame
preva.sk
colnames(preva.sk) <- "mean_preva"

skin.replacement.richness <- merge(replacement.plot.sk,richness.nonezero.sk, by.x="taxa", by.y="row.names")
skin.replacement.richness <- filter(skin.replacement.richness,mean_richness > 1)
skin.replacement.r.p <- merge(skin.replacement.richness, preva.sk, by.x = "taxa", by.y= "row.names")

skin.replacement.r.p <- skin.replacement.r.p %>% arrange(desc(mean_richness))
table(taxa_SK$Phylum) %>% sort()

#table(filter(taxa_SK, Phylum=="Firmicutes")$Class)
#table(filter(taxa_ST, Class=="Clostridia")$Family)

skin.replacement.r.p$phylum <- "Other"
skin.replacement.r.p$phylum[skin.replacement.r.p$taxa %in% filter(taxa_ST, Phylum=="Bacteroidetes")$Genus] <- "Bacteroidetes"
skin.replacement.r.p$phylum[skin.replacement.r.p$taxa %in% filter(taxa_ST, Phylum=="Firmicutes")$Genus] <- "Firmicutes"
skin.replacement.r.p$phylum[skin.replacement.r.p$taxa %in% filter(taxa_ST, Phylum=="Proteobacteria")$Genus] <- "Proteobacteria"
skin.replacement.r.p$phylum[skin.replacement.r.p$taxa %in% filter(taxa_ST, Phylum=="Actinobacteria")$Genus] <- "Actinobacteria"

#stool.replacement.richness$phylum[stool.replacement.richness$taxa %in% filter(taxa_ST, Phylum=="Verrucomicrobia")$Genus] <- "Verrucomicrobia"
ggplot(skin.replacement.r.p, aes(y=phylum, x=mean.dist)) + geom_jitter(color=ggsci::pal_jama()(n=7)[3])+ geom_boxplot(outlier.shape=NA, alpha=0.01) + base_theme

p.diver.pre.sk <- ggplot(skin.replacement.r.p, aes(x=mean_richness, y=mean_preva,color=phylum)) + geom_jitter() + geom_smooth(method = "lm") + facet_wrap(.~phylum, scales = "free") + base_theme
p.diver.pre.sk
#ggsave("../../../../XinZhouFiles/Projects/ZXP2_HMP2_Core/Analysis/Suppl.figure/rich_pre_skin.pdf",p.diver.pre.sk, width = 6, height = 4)

lm(mean_richness~mean_preva:phylum, skin.replacement.r.p) %>% summary()

###########################################################################
#Oral microbiome
replacement.plot.or <- filter(replacement.plot.3, Bodysite=="Oral")

is.na(richness_OR_byGenus) <- richness_OR_byGenus == 0
richness.nonezero.or <- colMeans(richness_OR_byGenus, na.rm = T) %>% data.frame
colnames(richness.nonezero.or) <- "mean_richness"
richness.nonezero.or

preva.or <- colMeans(OR.Pr)%>% data.frame
preva.or
colnames(preva.or) <- "mean_preva"

oral.replacement.richness <- merge(replacement.plot.or,richness.nonezero.or, by.x="taxa", by.y="row.names")
oral.replacement.richness <- filter(oral.replacement.richness,mean_richness > 1)
oral.replacement.r.p <- merge(oral.replacement.richness, preva.or, by.x = "taxa", by.y= "row.names")

oral.replacement.r.p <- oral.replacement.r.p %>% arrange(desc(mean_richness))
#oral and skin use same taxa table
table(taxa_SK$Phylum) %>% sort()

#table(filter(taxa_SK, Phylum=="Firmicutes")$Class)
#table(filter(taxa_ST, Class=="Clostridia")$Family)

oral.replacement.r.p$phylum <- "Other"
oral.replacement.r.p$phylum[oral.replacement.r.p$taxa %in% filter(taxa_SK, Phylum=="Bacteroidetes")$Genus] <- "Bacteroidetes"
oral.replacement.r.p$phylum[oral.replacement.r.p$taxa %in% filter(taxa_SK, Phylum=="Firmicutes")$Genus] <- "Firmicutes"
oral.replacement.r.p$phylum[oral.replacement.r.p$taxa %in% filter(taxa_SK, Phylum=="Proteobacteria")$Genus] <- "Proteobacteria"
oral.replacement.r.p$phylum[oral.replacement.r.p$taxa %in% filter(taxa_SK, Phylum=="Actinobacteria")$Genus] <- "Actinobacteria"

#stool.replacement.richness$phylum[stool.replacement.richness$taxa %in% filter(taxa_ST, Phylum=="Verrucomicrobia")$Genus] <- "Verrucomicrobia"
ggplot(oral.replacement.r.p, aes(y=phylum, x=mean.dist)) + geom_jitter(color=ggsci::pal_jama()(n=7)[4])+ geom_boxplot(outlier.shape=NA, alpha=0.01) + base_theme

p.diver.pre.or <- ggplot(oral.replacement.r.p, aes(x=mean_richness, y=mean_preva,color=phylum)) + geom_jitter() + geom_smooth(method = "lm") + facet_wrap(.~phylum, scales = "free") + base_theme
p.diver.pre.or

#ggsave("../../../../XinZhouFiles/Projects/ZXP2_HMP2_Core/Analysis/Suppl.figure/rich_pre_oral.pdf",p.diver.pre.or, width = 6, height = 4)

lm(mean_richness~mean_preva,stool.replacement.r.p %>% filter(phylum == "Proteobacteria")) %>% summary()

lm(mean_richness~mean_preva:phylum, oral.replacement.r.p) %>% summary()
###########################################################################
#Nasal microbiome
replacement.plot.ns <- filter(replacement.plot.3, Bodysite=="Nasal")

is.na(richness_NS_byGenus) <- richness_NS_byGenus == 0
richness.nonezero.ns <- colMeans(richness_NS_byGenus, na.rm = T) %>% data.frame
colnames(richness.nonezero.ns) <- "mean_richness"
richness.nonezero.ns

preva.ns <- colMeans(NS.Pr)%>% data.frame
preva.ns
colnames(preva.ns) <- "mean_preva"

nasal.replacement.richness <- merge(replacement.plot.ns,richness.nonezero.ns, by.x="taxa", by.y="row.names")
nasal.replacement.richness <- filter(nasal.replacement.richness,mean_richness > 1)
nasal.replacement.r.p <- merge(nasal.replacement.richness, preva.ns, by.x = "taxa", by.y= "row.names")

nasal.replacement.r.p <- nasal.replacement.r.p %>% arrange(desc(mean_richness))
#oral and skin use same taxa table
table(taxa_ST$Phylum) %>% sort()

#table(filter(taxa_SK, Phylum=="Firmicutes")$Class)
#table(filter(taxa_ST, Class=="Clostridia")$Family)

nasal.replacement.r.p$phylum <- "Other"
nasal.replacement.r.p$phylum[nasal.replacement.r.p$taxa %in% filter(taxa_ST, Phylum=="Bacteroidetes")$Genus] <- "Bacteroidetes"
nasal.replacement.r.p$phylum[nasal.replacement.r.p$taxa %in% filter(taxa_ST, Phylum=="Firmicutes")$Genus] <- "Firmicutes"
nasal.replacement.r.p$phylum[nasal.replacement.r.p$taxa %in% filter(taxa_ST, Phylum=="Proteobacteria")$Genus] <- "Proteobacteria"
nasal.replacement.r.p$phylum[nasal.replacement.r.p$taxa %in% filter(taxa_ST, Phylum=="Actinobacteria")$Genus] <- "Actinobacteria"

#stool.replacement.richness$phylum[stool.replacement.richness$taxa %in% filter(taxa_ST, Phylum=="Verrucomicrobia")$Genus] <- "Verrucomicrobia"
ggplot(nasal.replacement.r.p, aes(y=phylum, x=mean.dist)) + geom_jitter(color=ggsci::pal_jama()(n=7)[5])+ geom_boxplot(outlier.shape=NA, alpha=0.01) + base_theme

p.diver.pre.ns <- ggplot(nasal.replacement.r.p, aes(x=mean_richness, y=mean_preva,color=phylum)) + geom_jitter() + geom_smooth(method = "lm") + facet_wrap(.~phylum, scales = "free") + base_theme
p.diver.pre.ns

#ggsave("../../../../XinZhouFiles/Projects/ZXP2_HMP2_Core/Analysis/Suppl.figure/rich_pre_nasal.pdf",p.diver.pre.ns, width = 6, height = 4)

lm(mean_richness~mean_preva,nasal.replacement.r.p %>% filter(phylum == "Proteobacteria")) %>% summary()
lm(mean_richness~mean_preva:phylum, nasal.replacement.r.p) %>% summary()

stool.replacement.r.p %>% dim()
skin.replacement.r.p %>% dim()
oral.replacement.r.p %>% dim()
nasal.replacement.r.p %>% dim()

stool.replacement.r.p <- stool.replacement.r.p %>% mutate(phylum.bodysite = paste(phylum, Bodysite, sep="_"))
skin.replacement.r.p <- skin.replacement.r.p %>% mutate(phylum.bodysite = paste(phylum, Bodysite, sep="_"))
oral.replacement.r.p <- oral.replacement.r.p %>% mutate(phylum.bodysite = paste(phylum, Bodysite, sep="_"))
nasal.replacement.r.p <- nasal.replacement.r.p %>% mutate(phylum.bodysite = paste(phylum, Bodysite, sep="_"))

colnames(stool.replacement.r.p)
colnames(skin.replacement.r.p)
colnames(oral.replacement.r.p)
colnames(nasal.replacement.r.p)

stool.replacement.r.p.1 <- select(stool.replacement.r.p, taxa,Bodysite,mean.dist,mean.jd,meanrich,  
                                  mean_preva,phylum,phylum.bodysite)

skin.replacement.r.p.1 <- select(skin.replacement.r.p, taxa,Bodysite,mean.dist,mean.jd,meanrich,  
                                 mean_preva,phylum,phylum.bodysite)

oral.replacement.r.p.1 <- select(oral.replacement.r.p,taxa,Bodysite,mean.dist,mean.jd,meanrich,  
                                 mean_preva,phylum,phylum.bodysite)

nasal.replacement.r.p.1 <- select(nasal.replacement.r.p,taxa,Bodysite,mean.dist,mean.jd,meanrich,  
                                  mean_preva,phylum,phylum.bodysite)


combioned.replacement.r.p <- rbind(stool.replacement.r.p.1,
                                   skin.replacement.r.p.1,
                                   oral.replacement.r.p.1,
                                   nasal.replacement.r.p.1)
combioned.replacement.r.p %>% dim()

combined.plot <- ggplot(combioned.replacement.r.p, aes(y=phylum.bodysite, x=mean.dist)) + geom_jitter(aes(color=Bodysite))+ geom_boxplot(outlier.shape=NA, alpha=0.01) + base_theme
combined.plot <- combined.plot + facet_wrap(.~phylum, scales = "free") + scale_color_manual(values = body_site_color)
combined.plot

#############################################################################################
p.st <- stool.replacement.r.p %>% filter(meanrich > 1) %>% ggplot(aes(x=mean.jd,y=meanrich, label=taxa, color=phylum)) + geom_point() + geom_label_repel() + base_theme
p.st <- p.st + ggtitle("Stool")
p.st

p.sk <- skin.replacement.r.p %>% filter(mean_richness > 1.5) %>% ggplot(aes(x=mean.dist,y=mean_richness, label=taxa, color=phylum)) + geom_point() + geom_label_repel() + base_theme
p.sk <- p.sk + ggtitle("Skin")

p.or <- oral.replacement.r.p %>% filter(mean_richness > 1.5) %>% ggplot(aes(x=mean.dist,y=mean_richness, label=taxa, color=phylum)) + geom_point() + geom_label_repel() + base_theme
p.or <- p.or + ggtitle("Oral")

p.ns <- nasal.replacement.r.p %>% filter(mean_richness > 1.5) %>% ggplot(aes(x=mean.dist,y=mean_richness, label=taxa, color=phylum)) + geom_point() + geom_label_repel() + base_theme
p.ns <- p.ns +ggtitle("Nasal")

p.st + p.sk + p.or + p.ns

combioned.replacement.r.p %>% arrange(desc(mean.dist))

combioned.replacement.r.p %>% arrange((mean.dist))

colnames(combioned.replacement.r.p)

p.prev <- ggplot(combioned.replacement.r.p, aes(x=mean_preva, color=Bodysite)) + geom_density() + base_theme + scale_color_manual(values=body_site_color)
p.prev <- p.prev + labs(x="Prevalance", y="Density") 
p.prev

combioned.replacement.r.p$Bodysite <- factor(combioned.replacement.r.p$Bodysite, levels = c("Stool", "Skin", "Oral", "Nasal"))
p.prev.richness <- ggplot(combioned.replacement.r.p, aes(x=mean_richness, y=mean_preva, color=Bodysite)) + geom_jitter(size=0.7) + base_theme + scale_color_manual(values=body_site_color)
p.prev.richness <- p.prev.richness + labs(x="Mean_Richness", y="Prevalance")+ geom_smooth(method = "lm") + facet_wrap(.~Bodysite, scale="free")
p.prev.richness
#ggsave(filename = "../../Figure1/Diversity_Prevalance.pdf", width = 5, height = 4, dpi=300)


combioned.replacement.r.p

getwd()
load("../../../../data_analysis/combine_microbiome/distance/stool/personalized_score")
stool_personalized_score = personalized_score
stool_personalized_score$fc1 = stool_personalized_score$between_mean1 - stool_personalized_score$within_mean1

load("../../../../data_analysis/combine_microbiome/distance/skin/personalized_score")
skin_personalized_score = personalized_score
skin_personalized_score$fc1 = skin_personalized_score$between_mean1 - skin_personalized_score$within_mean1

load("../../../../data_analysis/combine_microbiome/distance/oral/personalized_score")
oral_personalized_score = personalized_score
oral_personalized_score$fc1 = oral_personalized_score$between_mean1 - oral_personalized_score$within_mean1

load("../../../../data_analysis/combine_microbiome/distance/nasal/personalized_score")
nasal_personalized_score = personalized_score
nasal_personalized_score$fc1 = nasal_personalized_score$between_mean1 - nasal_personalized_score$within_mean1

arrange(stool_personalized_score,desc(fc1))

dim(stool_personalized_score)
dim(skin_personalized_score)
dim(oral_personalized_score)
dim(nasal_personalized_score)

combined_personalizaed_score <- rbind(stool_personalized_score,skin_personalized_score,oral_personalized_score,nasal_personalized_score)
combined_personalizaed_score

dim(personalized_score)

temp <- data.frame(colMeans(ST.Pr)) 

stool.ps <- merge(stool_personalized_score,stool.replacement.r.p, by.x="genus", by.y="taxa", all.x = T)
stool.ps <- merge(stool.ps, data.frame(colMeans(ST.Pr)), by.x = "genus", by.y="row.names")
dim(stool.ps)
stool.ps

skin.ps <- merge(skin_personalized_score,skin.replacement.r.p, by.x="genus", by.y="taxa", all.x = T)
skin.ps <- merge(skin.ps, data.frame(colMeans(SK.Pr)), by.x = "genus", by.y="row.names")
dim(skin.ps)
skin.ps

oral.ps <- merge(oral_personalized_score,oral.replacement.r.p, by.x="genus", by.y="taxa", all.x = T)
oral.ps <- merge(oral.ps, data.frame(colMeans(OR.Pr)), by.x = "genus", by.y="row.names")
dim(oral.ps)
oral.ps

nasal.ps <- merge(nasal_personalized_score,nasal.replacement.r.p, by.x="genus", by.y="taxa", all.x = T)
nasal.ps <- merge(nasal.ps, data.frame(colMeans(NS.Pr)), by.x = "genus", by.y="row.names")
dim(nasal.ps)
nasal.ps

setdiff(stool_personalized_score$genus,stool.replacement.r.p$taxa)
setdiff(skin_personalized_score$genus,skin.replacement.r.p$taxa)
setdiff(oral_personalized_score$genus,oral.replacement.r.p$taxa)
setdiff(nasal_personalized_score$genus,nasal.replacement.r.p$taxa)

dim(stool.ps)
dim(skin.ps)
dim(oral.ps)
dim(nasal.ps)

psr.st <- ggplot(stool.ps, aes(x = mean.jd, y = fc1)) + geom_point(color=body_site_color[1]) + geom_smooth(method = "lm") + xlab("Jaccard Distance") + 
  ylab("DMI") + ggtitle("Stool") + base_theme + theme(axis.title=element_text(size=14)) +  theme(plot.title = element_text(size = 16, face = "bold"))
psr.sk <- ggplot(skin.ps, aes(x = mean.jd, y = fc1)) + geom_point(color=body_site_color[2]) + geom_smooth(method = "lm") + xlab("Jaccard Distance") + 
  ylab("DMI") + ggtitle("Skin") + base_theme + theme(axis.title=element_text(size=14)) +  theme(plot.title = element_text(size = 16, face = "bold"))
psr.or <- ggplot(oral.ps, aes(x = mean.jd, y = fc1)) + geom_point(color=body_site_color[3]) + geom_smooth(method = "lm") + xlab("Jaccard Distance") + 
  ylab("DMI") + ggtitle("Oral") + base_theme + theme(axis.title=element_text(size=14)) +  theme(plot.title = element_text(size = 16, face = "bold"))
psr.ns <- ggplot(nasal.ps, aes(x = mean.jd, y = fc1)) + geom_point(color=body_site_color[4]) + geom_smooth(method = "lm") + xlab("Jaccard Distance") + 
  ylab("DMI") + ggtitle("Nasal") + base_theme + theme(axis.title=element_text(size=14)) +  theme(plot.title = element_text(size = 16, face = "bold"))


psr <- psr.st + psr.sk + psr.or + psr.ns 
psr 




getwd()
load("../../../../data_analysis/combine_microbiome/distance/stool/personalized_score_permutation_trim")
stool.trim <- personalized_score_permutation_trim

load("../../../../data_analysis/combine_microbiome/distance/skin/personalized_score_permutation_trim")
skin.trim <- personalized_score_permutation_trim

load("../../../../data_analysis/combine_microbiome/distance/oral/personalized_score_permutation_trim")
oral.trim <- personalized_score_permutation_trim

load("../../../../data_analysis/combine_microbiome/distance/nasal/personalized_score_permutation_trim")
nasal.trim <- personalized_score_permutation_trim

colnames(stool.trim)
stool.trim$genus
colnames(stool.ps)
stool.ps$genus

# add bootstrap result
# load(file = "~/Desktop/1st.revision.CELL-S-22-04878/Figures/Figure3/detailed/revision.3a.data/Figure.Data.3A.RData")

# For stool dataset
stool.ps <- stool.ps %>% mutate(trim = "trim")
stool.ps$trim[stool.ps$genus %in% stool.trim$genus] <- "keep"
identical(stool.trim$genus, stool.ps$genus[stool.ps$trim == "keep"])

# For skin dataset
skin.ps <- skin.ps %>% mutate(trim = "trim")
skin.ps$trim[skin.ps$genus %in% skin.trim$genus] <- "keep"
identical(skin.trim$genus, skin.ps$genus[skin.ps$trim == "keep"])

# For oral dataset
oral.ps <- oral.ps %>% mutate(trim = "trim")
oral.ps$trim[oral.ps$genus %in% oral.trim$genus] <- "keep"
identical(oral.trim$genus, oral.ps$genus[oral.ps$trim == "keep"])

# For nasal dataset
nasal.ps <- nasal.ps %>% mutate(trim = "trim")
nasal.ps$trim[nasal.ps$genus %in% nasal.trim$genus] <- "keep"
identical(nasal.trim$genus, nasal.ps$genus[nasal.ps$trim == "keep"])

psr.st.1 <- filter(stool.ps, trim=="keep") %>% ggplot(aes(x = (1-mean.jd), y = fc1)) + geom_point(color=body_site_color[1]) + geom_smooth(method = "lm") + xlab("Taxa Recurrence Rate") + 
  ylab("DMI") + ggtitle("Stool") + base_theme + theme(axis.title=element_text(size=14)) +  theme(plot.title = element_text(size = 16, face = "bold"))
psr.sk.1 <- filter(skin.ps, trim=="keep") %>% ggplot(aes(x = (1-mean.jd), y = fc1)) + geom_point(color=body_site_color[2]) + geom_smooth(method = "lm") + xlab("Taxa Recurrence Rate") + 
  ylab("DMI") + ggtitle("Skin") + base_theme + theme(axis.title=element_text(size=14)) +  theme(plot.title = element_text(size = 16, face = "bold"))
psr.or.1 <- filter(oral.ps, trim=="keep") %>% ggplot(aes(x = (1-mean.jd), y = fc1)) + geom_point(color=body_site_color[3]) + geom_smooth(method = "lm") + xlab("Taxa Recurrence Rate") + 
  ylab("DMI") + ggtitle("Oral") + base_theme + theme(axis.title=element_text(size=14)) +  theme(plot.title = element_text(size = 16, face = "bold"))
psr.ns.1 <- filter(nasal.ps, trim=="keep") %>% ggplot(aes(x = (1-mean.jd), y = fc1)) + geom_point(color=body_site_color[4]) + geom_smooth(method = "lm") + xlab("Taxa Recurrence Rate") + 
  ylab("DMI") + ggtitle("Nasal") + base_theme + theme(axis.title=element_text(size=14)) +  theme(plot.title = element_text(size = 16, face = "bold"))

psr.1 <- psr.st.1 + psr.sk.1 + psr.or.1 + psr.ns.1 
psr.1 
#ggsave2("~/Desktop/1st.revision.CELL-S-22-04878/Figures/Figure3/detailed/Revision.fig3a.pdf",psr.1, height=5,width=7,dpi=300)
#save(stool.ps,skin.ps,oral.ps,nasal.ps,file = "~/Desktop/1st.revision.CELL-S-22-04878/Figures/Figure3/detailed/revision.3a.data/Figure.Data.3A.RData")


load("~/Desktop/1st.revision.CELL-S-22-04878/Figures/Figure3/detailed/revision.3a.data/Figure.Data.3A.RData")

stool.ps <- stool.ps %>%  mutate(trr= 1 - mean.jd)
skin.ps <- skin.ps %>%  mutate(trr= 1 - mean.jd)
oral.ps <- oral.ps %>%  mutate(trr= 1 - mean.jd)
nasal.ps <- nasal.ps %>%  mutate(trr= 1 - mean.jd)

test.stool <- cor.test(
  filter(stool.ps, trim == "keep")$trr,
  filter(stool.ps, trim == "keep")$fc1, method= "spearman")
print(test.stool)

table(stool.ps$trim)

test_skin <- cor.test(
  filter(skin.ps, trim == "keep")$trr,
  filter(skin.ps, trim == "keep")$fc1, method= "spearman")
print(test_skin)
test_skin

table(skin.ps$trim)

test_oral <- cor.test(
  filter(oral.ps, trim == "keep")$trr,
  filter(oral.ps, trim == "keep")$fc1, method= "spearman")
print(test_oral)

table(oral.ps$trim)

test_nasal <- cor.test(
  filter(nasal.ps, trim == "keep")$trr,
  filter(nasal.ps, trim == "keep")$fc1, method= "spearman")
print(test_nasal)

table(nasal.ps$trim)

test <- cor.test(skin.ps$mean.jd, skin.ps$fc1)
test$p.value

test <- cor.test(oral.ps$mean.jd, oral.ps$fc1)
test$p.value

test <- cor.test(nasal.ps$mean.jd, nasal.ps$fc1)
test$p.value

#make correlation based on different body sites
psr.st.2 <- ggscatterstats(stool.ps, x = mean.jd, y = fc1, type = "nonparametric", title = "Stool_JD")
psr.st.2
#ggsave(filename = "../Figure3a_suppl/Stool.corre.pdf",psr.st.2, width = 6, height = 5, dpi = 300)

psr.sk.2 <- ggscatterstats(skin.ps, x = mean.jd, y = fc1, type = "nonparametric", title = "Skin_JD")
psr.sk.2
#ggsave(filename = "../Figure3a_suppl/Skin.corre.pdf", psr.sk.2, width = 6, height = 5, dpi = 300)

psr.or.2 <- ggscatterstats(oral.ps, x = mean.jd, y = fc1, type = "nonparametric", title = "Oral_JD")
psr.or.2
#ggsave(filename = "../Figure3a_suppl/Oral.corre.pdf", psr.or.2 ,width = 6, height = 5, dpi = 300)

psr.ns.2 <- ggscatterstats(nasal.ps, x = mean.jd, y = fc1, type = "nonparametric", title = "Nasal_JD")
psr.ns.2
#ggsave(filename = "../Figure3a_suppl/Nasal.corre.pdf",psr.ns.2, width = 6, height = 5, dpi = 300)

test <- cor.test(stool.ps$mean.jd, stool.ps$fc1)
test$p.value

test <- cor.test(skin.ps$mean.jd, skin.ps$fc1)
test$p.value

test <- cor.test(oral.ps$mean.jd, oral.ps$fc1)
test$p.value

test <- cor.test(nasal.ps$mean.jd, nasal.ps$fc1)
test$p.value

#ggsave(filename = "../Figure3a.pdf", psr, width = 5, height = 5, dpi=300)

summary(lm(fc1 ~ mean.jd, data=stool.ps)) 
summary(lm(fc1 ~ mean.jd, data=skin.ps))
summary(lm(fc1 ~ mean.jd, data=oral.ps)) 
summary(lm(fc1 ~ mean.jd, data=nasal.ps)) 

colnames(stool.ps)
colnames(skin.ps)
colnames(oral.ps)
colnames(nasal.ps)

binded.phy<- rbind(select(stool.ps,genus,phylum,Bodysite,mean.jd),
      select(skin.ps,genus,phylum,Bodysite,mean.jd),
      select(oral.ps,genus,phylum,Bodysite,mean.jd),
      select(nasal.ps,genus,phylum,Bodysite,mean.jd))
 
binded.phy <- binded.phy %>% mutate(replace.rate = 1 - mean.jd)
 
binded.df <- filter(binded.phy, !is.na(phylum))
binded.df$phylum <- factor(binded.df$phylum, levels =c("Actinobacteria", "Bacteroidetes","Firmicutes","Proteobacteria","Other"))
binded.df$Bodysite <- factor(binded.df$Bodysite, levels=c("Stool", "Skin", "Oral", "Nasal"))
 
p.all <- binded.df %>% ggplot(aes(x=Bodysite, y=replace.rate)) + geom_jitter(size=0.1) + geom_boxplot(aes(color=Bodysite), alpha=0.7, outlier.alpha = 0) + facet_wrap(.~phylum)
p.all

############################################################################################################################################
Stool.dist.bysubject <- Stool.replacement %>% group_by(SubjectID) %>% 
  mutate(mean_timediff = mean(timediff)) %>% 
  mutate(mean_jd = mean(jaccard)) %>% mutate(mean_dist = mean(dist)) %>% 
  select(SubjectID,mean_timediff,mean_jd,mean_dist) %>% unique()

Skin.dist.bysubject <- Skin.replacement %>% group_by(SubjectID) %>% 
  mutate(mean_timediff = mean(timediff)) %>% 
  mutate(mean_jd = mean(jaccard)) %>% mutate(mean_dist = mean(dist)) %>% 
  select(SubjectID,mean_timediff,mean_jd,mean_dist) %>% unique()

Oral.dist.bysubject <- Oral.replacement %>% group_by(SubjectID) %>% 
  mutate(mean_timediff = mean(timediff)) %>% 
  mutate(mean_jd = mean(jaccard)) %>% mutate(mean_dist = mean(dist)) %>% 
  select(SubjectID,mean_timediff,mean_jd,mean_dist) %>% unique()

Nasal.dist.bysubject <- Nasal.replacement %>% group_by(SubjectID) %>% 
  mutate(mean_timediff = mean(timediff)) %>% 
  mutate(mean_jd = mean(jaccard)) %>% mutate(mean_dist = mean(dist)) %>% 
  select(SubjectID,mean_timediff,mean_jd,mean_dist) %>% unique()

Stool.dist.bysubject$SampleType <- "ST"
Skin.dist.bysubject$SampleType <- "Skin"
Oral.dist.bysubject$SampleType <- "Oral"
Nasal.dist.bysubject$SampleType <- "NS"
dist.bysubject <- rbind(Stool.dist.bysubject,Skin.dist.bysubject,Oral.dist.bysubject,Nasal.dist.bysubject)

load("/Users/xzhou7/Library/CloudStorage/Box-Box/XinZhouFiles/Projects/ZXP2_HMP2_Core/Analysis/Robject/Old_Diversity_Datatable.RData")
ID.list <- rbind(select(stool.nasal.ASV.diversity,RandomID,SampleID,SubjectID,SampleType),
      select(skin.oral.ASV.diversity,KitID,SampleID,SubjectID,SampleType) %>% mutate(RandomID=KitID) %>% select(-KitID))
dim(ID.list)
all.diversity <- read.csv(file = "../../../../../XinZhouFiles/Projects/ZXP2_HMP2_Core/Analysis/Diversity Table/All_Diversity.csv", header = T)
diversity_meta <- merge(all.diversity, ID.list, by.x="X", by.y ="RandomID")
dim(diversity_meta)

diversity_meta_mean <- diversity_meta %>% select(Observed, Shannon,pielou_e,SubjectID,SampleType) %>% group_by(SubjectID,SampleType) %>%
  mutate(mean_shannon = mean(Shannon)) %>% mutate(mean_observed = mean(Observed)) %>% mutate(mean_pielou_e = mean(pielou_e)) %>% select(-Observed, -Shannon,-pielou_e) %>% unique()

table(diversity_meta_mean$SampleType)
sc2 <- read.csv("../../../../../XinZhouFiles/Projects/ZXP2_HMP2_Core/Analysis/ori_meta_table/metadata.subject.csv", header = T)
sc2

diversitymean_meta <- merge(diversity_meta_mean, sc2, by="SubjectID")
diversity.DIstance_meta <- merge(diversitymean_meta, dist.bysubject, by=c("SubjectID","SampleType"))
diversity.DIstance_meta

p <- filter(diversity.DIstance_meta, IRIS != "Unknown") %>% ggplot(aes(x=IRIS, y=mean_jd)) + geom_boxplot()
p <- p + facet_wrap(.~SampleType) + stat_compare_means()
p

###check the mean value of BC distance within individuals and their relationship with this individuals' Shannon Diversity
diversity.DIstance_meta <- diversity.DIstance_meta %>% mutate(BC_Similarity = 1 - mean_dist)
diversity.DIstance_meta$IRIS <- factor(diversity.DIstance_meta$IRIS, levels = c("IS", "IR", "Unknown"))

pst.di.sta <- filter(diversity.DIstance_meta, SampleType == "ST") %>% ggscatterstats(mean_shannon,BC_Similarity,type="nonparametric", xlab="Shannon Diversity", ylab = "Similarity", title="Stool Microbiome")
psk.di.sta <- filter(diversity.DIstance_meta, SampleType == "Skin") %>% ggscatterstats(mean_shannon,BC_Similarity,type="nonparametric",xlab="Shannon Diversity", ylab = "Similarity", title="Skin Microbiome")
por.di.sta <- filter(diversity.DIstance_meta, SampleType == "Oral") %>% ggscatterstats(mean_shannon,BC_Similarity,type="nonparametric", xlab="Shannon Diversity", ylab = "Similarity", title="Oral Microbiome")
pns.di.sta <- filter(diversity.DIstance_meta, SampleType == "NS") %>% ggscatterstats(mean_shannon,BC_Similarity, type="nonparametric",xlab="Shannon Diversity", ylab = "Similarity", title="Nasal Microbiome")

p.di.sta <- pst.di.sta + psk.di.sta + por.di.sta + pns.di.sta
p.di.sta
#ggsave(filename = "./Stability_Diversity.pdf",p.di.sta, width = 14, height = 10, dpi = 300)


pst.jd.IRIS <- filter(diversity.DIstance_meta, SampleType == "ST" & IRIS != "Unknown") %>% ggbetweenstats(x=IRIS, y=mean_jd, xlab="", ylab = "Jaccard Similarity", title="Stool Microbiome")
pst.jd.IRIS
psk.jd.IRIS <- filter(diversity.DIstance_meta, SampleType == "Skin" & IRIS != "Unknown") %>% ggbetweenstats(x=IRIS, y=mean_jd, xlab="", ylab = "Jaccard Similarity", title="Skin Microbiome")
psk.jd.IRIS
por.jd.IRIS <- filter(diversity.DIstance_meta, SampleType == "Oral" & IRIS != "Unknown") %>% ggbetweenstats(x=IRIS, y=mean_jd, xlab="", ylab = "Jaccard Similarity", title="Oral Microbiome")
por.jd.IRIS
pns.jd.IRIS <- filter(diversity.DIstance_meta, SampleType == "NS" & IRIS != "Unknown") %>% ggbetweenstats(x=IRIS, y=mean_jd, xlab="", ylab = "Jaccard Similarity", title="Nasal Microbiome")
pns.jd.IRIS

p.jd.IRIS <- pst.jd.IRIS + psk.jd.IRIS + por.jd.IRIS + pns.jd.IRIS
p.jd.IRIS
#ggsave(filename = "./JD_IRIS_meanbysubject.pdf",p.jd.IRIS, width = 14, height = 10, dpi = 300)
################################################################################################################
mean(stool_personalized_score$fc1)
mean(skin_personalized_score$fc1)
mean(oral_personalized_score$fc1)
mean(nasal_personalized_score$fc1)

sd(stool_personalized_score$fc1)
sd(skin_personalized_score$fc1)
sd(oral_personalized_score$fc1)
sd(nasal_personalized_score$fc1)


colnames(stool_personalized_score)


#######power estimation
library("mediation")
library("dplyr")
library("ggplot2")

set.seed(777)

n_rep <- 1000  # Number of simulations
n <- 300       # Sample size
alpha <- 0.05  # Significance level
effect_size <- 0.3


simulate_data <- function(n, effect_size) {
  X <- rnorm(n)
  M <- effect_size * X + rnorm(n)
  Y <- effect_size * M + rnorm(n)
  data.frame(X = X, M = M, Y = Y)
}

data <- simulate_data(n, effect_size)

mediator_model <- lm(M ~ X, data = data)
outcome_model <- lm(Y ~ M + X, data = data)

med <- mediate(mediator_model,
               outcome_model,
               treat = "X",
               mediator = "M",
               sims = 500,
               conf.level = 1 - alpha)

power_mediation <- function(n_rep, n, effect_size, alpha) {
  results <- replicate(n_rep, {
    data <- simulate_data(n, effect_size)
    med <- mediate(lm(M ~ X, data = data),
                   lm(Y ~ M + X, data = data),
                   treat = "X",
                   mediator = "M",
                   sims = 500,
                   conf.level = 1 - alpha)
    as.numeric(med$d0 < alpha)
  })
  mean(results)
}

power <- power_mediation(n_rep, n, effect_size, alpha)
power
#power = 0.047


