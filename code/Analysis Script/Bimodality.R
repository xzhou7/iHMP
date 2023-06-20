#bimodality
#author Xin Zhou

library(microbiome)
library(dplyr)
library(plyr)
library(tidyr)
library(ggrepel)
library(gridExtra)
library(ggplot2)
library(cowplot)
library(patchwork)

setwd("~/Box/XinZhouFiles/Projects/ZXP2_HMP2_Core/")
load("./Analysis/Robject/DetailedPhyloseq.RData")
metadata <- read.csv("./Analysis/ori_meta_table/metadata.subject.csv", header = T, row.names = 1)

physeqmeta <- data.frame(sample_data(physeq_ST))
physeqmeta$IRIS <- "Unknown"
physeqmeta$IRIS[physeqmeta$SubjectID %in% metadata$SubjectID[metadata$IRIS=="IR"]] <- "IR"
physeqmeta$IRIS[physeqmeta$SubjectID %in% metadata$SubjectID[metadata$IRIS=="IS"]] <- "IS"

newsampleame <- sample_data(physeq_ST) %>% data.frame() %>%
  group_by(SubjectID) %>%
  mutate(time = Date - min(Date))

sample_data(physeq_ST)$time <- newsampleame$time
sample_data(physeq_ST)$subject <- sample_data(physeq_ST)$SubjectID

sample_data(physeq_ST)$IRIS <- physeqmeta$IRIS

physeqmeta5 <-subset(physeqmeta, SubjectID %in% names(which(table(physeqmeta$SubjectID) > 5)))
physeqmeta.less5 <-subset(physeqmeta, SubjectID %in% names(which(table(physeqmeta$SubjectID) < 6)))
set.seed(143)
physeqmeta5_subset <- ddply(physeqmeta5,.(SubjectID),function(x) x[sample(nrow(x),5),])
physeqmeta_random <- rbind(physeqmeta.less5, physeqmeta5_subset)
table(physeqmeta_random$SubjectID)

physeqmeta_random_IS <- subset(physeqmeta_random, IRIS=="IS")
physeqmeta_random_IR <- subset(physeqmeta_random, IRIS=="IR")

physeqis <- subset_samples(physeq_ST, SampleID %in% physeqmeta_random_IS$SampleID)
physeqir <- subset_samples(physeq_ST, SampleID %in% physeqmeta_random_IR$SampleID)

physeqis
physeqir

physeqis
physeqis <- microbiome::transform(physeqis, "compositional")
physeqis <- aggregate_rare(physeqis, level = "Genus", detection = .1/100, prevalence = 10/100)
physeqis.clr <- microbiome::transform(physeqis, "clr")

physeqir
physeqir <- microbiome::transform(physeqir, "compositional")
physeqir <- aggregate_rare(physeqir, level = "Genus", detection = .1/100, prevalence = 10/100)
physeqir.clr <- microbiome::transform(physeqir, "clr")

Sarle.finite.score.is <- bimodality(physeqis.clr, method = "Sarle.finite.sample", bs.iter = 1000, verbose=T)
Sarle.finite.score.ir <- bimodality(physeqir.clr, method = "Sarle.finite.sample", bs.iter = 1000, verbose=T)

Sarle.asymptotic.score.is <- bimodality(physeqis.clr, method = "Sarle.asymptotic", bs.iter = 1000, verbose=T)
Sarle.asymptotic.score.ir <- bimodality(physeqir.clr, method = "Sarle.asymptotic", bs.iter = 1000, verbose=T)

Bio.score.is <- bimodality(physeqis.clr, method = "potential_analysis", bs.iter = 1000, peak.threshold=20, min.density=6, verbose=T)
Bio.score.ir <- bimodality(physeqir.clr, method = "potential_analysis", bs.iter = 1000, peak.threshold=20, min.density=6, verbose=T)
####################################################################################################################################
df.loop <- data.frame()
taxa <- taxa(physeqis)
for (i in 1:10){
  print(i)
  Bio.score.is <- bimodality(physeqis.clr, method = "potential_analysis", bs.iter = 100, peak.threshold=20, min.density=(3+0.5*i), verbose=T)
  df.loop[1:82, i] <- Bio.score.is
}
row.names(df.loop) <- taxa
df.loop$taxa <- taxa
df.loop.gather1 <- df.loop %>% gather("key","value",-taxa)
df.loop.gather1$IRIS <- "IS"

df.loop2 <- data.frame()
taxa2 <- taxa(physeqir)
for (i in 1:10){
  print(i)
  Bio.score.ir <- bimodality(physeqir.clr, method = "potential_analysis", bs.iter = 100, peak.threshold=20, min.density= (3+0.5*i), verbose=T)
  df.loop2[1:72, i] <- Bio.score.ir
}
row.names(df.loop2) <- taxa2
df.loop2$taxa <- taxa2
df.loop.gather2 <- df.loop2 %>% gather(key = "key", value="value", -taxa)
df.loop.gather2$IRIS <- "IR"

df.loop.gather <- rbind(df.loop.gather1,df.loop.gather2)
df.loop.gather$key <- factor(df.loop.gather$key, levels = paste("V", 1:10,sep=""))

p.density.adjust <- ggplot(df.loop.gather, aes(x=key, y=value, group=IRIS, color=IRIS)) + geom_line()
p.density.adjust <- p.density.adjust + facet_wrap(.~taxa, ncol = 10, scales = "free")
p.density.adjust


colSums(df.loop==1)
colSums(df.loop2==1)
colSums(df.loop==0)
colSums(df.loop2==0)

df.loop_combined <- merge(df.loop,df.loop2, by="row.names")
View(df.loop_combined)
#df <- data.frame(group = taxa, Sarle.finite.is = Sarle.finite.score.is[taxa], Sarle.finite.ir = Sarle.finite.score.ir[taxa],
#                Sarle.asymptotic.is = Sarle.asymptotic.score.is[taxa], Sarle.asymptotic.ir = Sarle.asymptotic.score.ir[taxa], 
#                Bio.score.is = Bio.score.is[taxa], Bio.score.ir = Bio.score.ir[taxa])
####################################################################################################################################
tax <- "Unclassified_Bacteroidales"

x.is <- data.frame(abundances(microbiome::transform(physeqis, "clr"))[tax,])
x.is$group <- "IS"
colnames(x.is)[1] <- "abundance"
x.ir <- data.frame(abundances(microbiome::transform(physeqir, "clr"))[tax,])
x.ir$group <- "IR"
colnames(x.ir)[1] <- "abundance"
xdensity <- rbind(x.is, x.ir)

p <- ggplot(xdensity, aes(x=abundance, color = group)) + geom_density() + theme_classic() + ggtitle(tax)
p

#unimodal  <- names(sort(bimodality.score))[[1]]
#bimodal  <- rev(names(sort(bimodality.score)))[[1]]
#p1 <- plot_density(physeq.clr, variable = unimodal) 
#p2 <- plot_density(physeq.clr, variable = bimodal) 
#p1 + p2

#df.is <- df
#df.ir <- df

physeq_IS <- subset_samples(physeq_ST, IRIS == "IS")
physeq_IR <- subset_samples(physeq_ST, IRIS == "IR")

physeq_IS <- microbiome::transform(physeq_IS, "compositional")
physeq_IS <- aggregate_rare(physeq_IS, level = "Genus", detection = .1/100, prevalence = 10/100)
physeq_IS.clr <- microbiome::transform(physeq_IS, "clr")

physeq_IR <- microbiome::transform(physeq_IR, "compositional")
physeq_IR <- aggregate_rare(physeq_IR, level = "Genus", detection = .1/100, prevalence = 10/100)
physeq_IR.clr <- microbiome::transform(physeq_IR, "clr")

intermediate.stability_IS <- intermediate_stability(physeq_IS, output = "scores")
intermediate.stability_IS

intermediate.stability_IR <- intermediate_stability(physeq_IR, output = "scores")
intermediate.stability_IR

taxa <- unique(taxa(physeqis), taxa(physeqir))
df <- data.frame(group = taxa, Sarle.finite.is = Sarle.finite.score.is[taxa], Sarle.finite.ir = Sarle.finite.score.ir[taxa],
           Sarle.asymptotic.is = Sarle.asymptotic.score.is[taxa], Sarle.asymptotic.ir = Sarle.asymptotic.score.ir[taxa], 
           Bio.score.is = Bio.score.is[taxa], Bio.score.ir = Bio.score.ir[taxa], 
           Im_stab_IS = intermediate.stability_IS[taxa], Im_stab_IR = intermediate.stability_IR[taxa])

df$Sarle.finite_delta <- df$Sarle.finite.is - df$Sarle.finite.ir
df$Sarle.asymptotic_delta <- df$Sarle.asymptotic.is -df$Sarle.asymptotic.ir
df$Bio.score_delta <- df$Bio.score.is -df$Bio.score.ir
df$sta_delta <- df$Im_stab_IS - df$Im_stab_IR
df

listcore <- c("Coprococcus","Butyricicoccus","Barnesiella","Parasutterella","Butyrivibrio")
controllist <- c("Phocaeicola","Akkermansia", "Faecalibacterium")
meandiflist <- c("Bacteroides", "Oscillibacter" )

p1 <- ggplot(df, aes(x=Im_stab_IS, y=Bio.score.is)) + geom_point() + geom_smooth(method = "lm", se=F)
p1 <- p1 + xlim(-0.4, 0.4) + geom_text_repel(data=subset(df,Bio.score.is > 0.75 & Im_stab_IS < (-0.2)),aes(x=Im_stab_IS, y=Bio.score.is, label=group))
p1

p2 <- ggplot(df, aes(x=Im_stab_IR, y=Bio.score.ir)) + geom_point() + geom_smooth(method = "lm", se=F)
p2 <- p2 + xlim(-0.4, 0.4) + geom_text_repel(data=subset(df,Bio.score.ir > 0.75 & Im_stab_IR < (-0.2)),aes(x=Im_stab_IR, y=Bio.score.ir, label=group))
p2

p1+p2

p1 <- ggplot(df, aes(x=Im_stab_IS, y=Bio.score.is)) + geom_point() + geom_smooth(method = "lm", se=F)
p1 <- p1 + xlim(-0.4, 0.4) + geom_point(data=subset(df, group %in% listcore),aes(x=Im_stab_IS, y=Bio.score.is), color="red")
p1 <- p1 + geom_point(data=subset(df, group %in% controllist),aes(x=Im_stab_IS, y=Bio.score.is), color="blue")
p1 <- p1 + geom_point(data=subset(df, group %in% meandiflist),aes(x=Im_stab_IS, y=Bio.score.is), color="green")
p1

p2 <- ggplot(df, aes(x=Im_stab_IR, y=Bio.score.ir)) + geom_point() + geom_smooth(method = "lm", se=F)
p2 <- p2 + xlim(-0.4, 0.4) + geom_point(data=subset(df, group %in% listcore),aes(x=Im_stab_IR, y=Bio.score.ir), color="red")
p2 <- p2 + geom_point(data=subset(df, group %in% controllist),aes(x=Im_stab_IR, y=Bio.score.ir), color="blue")
p2 <- p2 + geom_point(data=subset(df, group %in% meandiflist),aes(x=Im_stab_IR, y=Bio.score.ir), color="green")
p2

p1+p2

p <- p + geom_text_repel(data=subset(df,, group %in% listcore),aes(x=Bio.score_delta, y=sta_delta, label=group)) 
p <- p + geom_point(data=subset(df, group %in% listcore),aes(x=Bio.score_delta, y=sta_delta), color="red")
p <- p + geom_point(data=subset(df, group %in% controllist),aes(x=Bio.score_delta, y=sta_delta), color="blue")
p <- p + geom_point(data=subset(df, group %in% meandiflist),aes(x=Bio.score_delta, y=sta_delta), color="green")
p

p <- ggplot(df, aes(x=Bio.score.is, y=Im_stab_IS)) + geom_point()
p <- p + geom_text_repel(data=subset(df,, group %in% listcore),aes(x=Bio.score_delta, y=sta_delta, label=group)) 
p <- p + geom_point(data=subset(df, group %in% listcore),aes(x=Bio.score_delta, y=sta_delta), color="red")
p <- p + geom_point(data=subset(df, group %in% controllist),aes(x=Bio.score_delta, y=sta_delta), color="blue")
p <- p + geom_point(data=subset(df, group %in% meandiflist),aes(x=Bio.score_delta, y=sta_delta), color="green")
p


listtaxa <- subset(df, !is.na(df$Sarle.finite.ir))

cor(df$Im_stab_IS, df$Bio.score.is,use = ,method = "pearson")

for (i in 1:67){
  print(i)
  print(listtaxa$group[i])
  tax <- (listtaxa$group)[i]
  x.is <- data.frame(abundances(microbiome::transform(physeqis, "clr"))[tax,])
  x.is$group <- "IS"
  colnames(x.is)[1] <- "abundance"
  x.ir <- data.frame(abundances(microbiome::transform(physeqir, "clr"))[tax,])
  x.ir$group <- "IR"
  colnames(x.ir)[1] <- "abundance"
  xdensity <- rbind(x.is, x.ir)
  ksresult <- ks.test(x.is$abundance, x.ir$abundance)
  
  p <- ggplot(xdensity, aes(x=abundance, color = group)) + geom_density() + theme_classic() + labs(title=tax, subtitle = p.adjust(ksresult$p.value, n=67))
  #path <- paste0("~/Box/XinZhouFiles/Projects/ZXP2_HMP2_Core/Analysis/Suppl.figure/biomodality/ST/", tax, ".pdf")
  #ggsave(filename = path, p, height = 1.5, width = 3, dpi = 300)
  listtaxa[i,14] <- ksresult$p.value
}
colnames(listtaxa)[14] <- "KS.Pvalue"
listtaxa$logp <- log(listtaxa$KS.Pvalue)

p <- ggplot(listtaxa, aes(x=Bio.score_delta,y=-logp,label=group)) + geom_point() + geom_text_repel()
p <- p + geom_point(data=subset(listtaxa, group %in% listcore),aes(x=Bio.score_delta, y=-logp), color="red")
p <- p + geom_point(data=subset(listtaxa, group %in% controllist),aes(x=Bio.score_delta, y=-logp), color="blue")
p <- p + geom_point(data=subset(listtaxa, group %in% meandiflist),aes(x=Bio.score_delta, y=-logp), color="green")
p


round(ksresult$p.value, digits = 9)

p.adjust(0.00247, n=67)

tax  <- names(which.max(bimodality.score))
bimodality.score
sort(bimodality.score)

tax <- "Prevotella"
x <- abundances(microbiome::transform(physeq_IS, "clr"))[tax,]
potential.minima <- potential_analysis(x, bs.iter = 100, min.density = 3)$minima
tipping.samples <- sapply(potential.minima, function (m) {names(which.min(abs(sort(x) - m)))})
tipping.point <- abundances(physeq_IS)[tax, tipping.samples]
print(tipping.point)
plot(density(x), main = tax)
abline(v = x[tipping.samples])

x <- abundances(microbiome::transform(physeq_IR, "clr"))[tax,]
potential.minima <- potential_analysis(x, bs.iter = 100,min.density = 3)$minima
tipping.samples <- sapply(potential.minima, function (m) {names(which.min(abs(sort(x) - m)))})
tipping.point <- abundances(physeq_IR)[tax, tipping.samples]
print(tipping.point)
plot(density(x), main = tax)
abline(v = x[tipping.samples])

save(df, listtaxa, file = "./Analysis/Suppl.figure/biomodality/STRData/StoolBimodality.RData")

