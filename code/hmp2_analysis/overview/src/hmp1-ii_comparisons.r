
library(plyr)
library(dplyr)

source("./common/disease_colors.r")
source("./common/load_bugs.r")

hmp12_bugs.pcl <- pcl.read(file.path(HMP2_data, "hmp1-ii/hmp1-II_metaphlan2-mtd-qcd.tsv"), metadata.rows=8) %>%
    pcl.filter.s(STSite == "Stool")

# Merge the datasets
hmp12_bugs.pcl$meta$diagnosis <- "HMP1-II"
hmp12_bugs.pcl$meta$VISNO[hmp12_bugs.pcl$meta$VISNO=="02S"] <- "2"
hmp12_bugs.pcl$meta$VISNO[hmp12_bugs.pcl$meta$VISNO=="03S"] <- "3"
hmp12_bugs.pcl$meta$VISNO <- as.numeric(as.character(hmp12_bugs.pcl$meta$VISNO))
merged.pcl <- pcl.merge(bugs.pcl, hmp12_bugs.pcl)

merged.pcl$meta$merged_subj <- as.character(merged.pcl$meta$subject)
merged.pcl$meta$merged_subj[is.na(merged.pcl$meta$merged_subj)] <-
    as.character(merged.pcl$meta$RANDSID)[is.na(merged.pcl$meta$merged_subj)]
merged.pcl$meta$merged_subj <- factor(merged.pcl$meta$merged_subj)

diag_subj <- merged.pcl$meta[match(levels(merged.pcl$meta$merged_subj), merged.pcl$meta$merged_subj),"diagnosis",drop=F]
rownames(diag_subj) <- levels(merged.pcl$meta$merged_subj)


# PERMANOVA to see if we can distinguish HMP1-II from non-IBD?
merged.pcl.healthy <- merged.pcl %>%
    pcl.filter.s(diagnosis == "nonIBD" || diagnosis == "HMP1-II") %>%
    pcl.only(rank="s") %>% pcl.nicenames %>%
    pcl.normalize

ad.naive <- adonis(merged.pcl.healthy$x ~ merged.pcl.healthy$meta$diagnosis, method="bray", permutations=999)
ad.naive
# Permutation: free
# Number of permutations: 999
#
# Terms added sequentially (first to last)
#
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# merged.pcl.healthy$meta$diagnosis   1      9.08  9.0804  33.105 0.03277  0.001 ***
# Residuals                         977    267.99  0.2743         0.96723
# Total                             978    277.06                 1.00000
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

library(vegan)
D <- vegdist(merged.pcl.healthy$x, method="bray");

ad <- PERMANOVA_repeat_measures(
    D, merged.pcl.healthy$meta[,c(),drop=F],
    factor(merged.pcl.healthy$meta$merged_subj, levels=rownames(diag_subj)),
    diag_subj, permutations=999)
ad
# Call:
# adonis(formula = D ~ ., data = mtdat[, metadata_order, drop = F],      permutations = 0)
#
# Permutation: free
# Number of permutations: 0
#
# Terms added sequentially (first to last)
#
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# diagnosis   1      9.08  9.0804  33.105 0.03277  0.001 ***
# Residuals 977    267.99  0.2743         0.96723  0.001 ***
# Total     978    277.06                 1.00000
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

library(vegan)
merged.pcl.healthy.adult <- pcl.filter.s(merged.pcl.healthy, adult | (diagnosis == "HMP1-II"))
Dadult <- vegdist(merged.pcl.healthy.adult$x, method="bray");

subjfact <- factor(merged.pcl.healthy.adult$meta$merged_subj)
ad <- PERMANOVA_repeat_measures(
    Dadult, merged.pcl.healthy.adult$meta[,c(),drop=F], subjfact,
    diag_subj[match(levels(subjfact), rownames(diag_subj)),,drop=F], permutations=999)
ad
# Call:
#     adonis(formula = D ~ ., data = mtdat[, metadata_order, drop = F],      permutations = 0)
#
# Permutation: free
# Number of permutations: 0
#
# Terms added sequentially (first to last)
#
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# diagnosis   1     6.796  6.7965  25.437 0.03198  0.001 ***
#     Residuals 770   205.736  0.2672         0.96802  0.001 ***
#     Total     771   212.533                 1.00000
# ---
#     Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Species-level joint ordination
joint.ord <- pcl.pcoa(merged.pcl %>% pcl.only(rank="s") %>% pcl.normalize)

pdf("./overview/joint_ordination_hmp1ii_alldiagnosis.pdf", 6.5, 5.5)
pcl.ordplot(merged.pcl, joint.ord, colour="diagnosis",
            colour_override=c(hmp2_disease_colors, c("HMP1-II" = "forestgreen")),
            colour_title = "Diagnosis", size_abs = 2)
dev.off()

# Species-level joint ordination (only non-IBD samples)
joint.nonibd.ord <- pcl.pcoa(merged.pcl %>% pcl.filter.s(diagnosis %in% c("HMP1-II", "nonIBD")) %>%
    pcl.only(rank="s") %>% pcl.normalize)

pdf("./overview/joint_ordination_hmp1ii_nonibd.pdf", 6.5, 5.5)
pcl.ordplot(merged.pcl, joint.nonibd.ord, colour="diagnosis",
            colour_override=c(hmp2_disease_colors, c("HMP1-II" = "forestgreen")),
            colour_title = "Diagnosis", size_abs = 2)
dev.off()

pdf("./overview/joint_ordination_hmp1ii_bacteroides.pdf", 6.5, 5.5)
pcl.ordplot(merged.pcl %>% pcl.nicenames,
            joint.nonibd.ord,
            colour="Bacteroides",
            size_abs = 3)
dev.off()

pdf("./overview/joint_ordination_hmp1ii_bacteroides_ovatus.pdf", 6.5, 5.5)
pcl.ordplot(merged.pcl %>% pcl.nicenames,
            joint.nonibd.ord,
            colour="Bacteroides ovatus",
            size_abs = 3)
dev.off()


pdf("./overview/joint_heatmap_hmp1ii.pdf", 12, 5)
merged.pcl %>%
    pcl.only(rank="s") %>%
    pcl.nicenames %>%
    pcl.top.f(mean(x), n=10) %>%
    pcl.heatmap(meta=c("diagnosis"), annotation_colors=list(diagnosis=c(hmp2_disease_colors, c("HMP1-II" = "forestgreen"))))
merged.pcl %>%
    pcl.only(rank="s") %>%
    pcl.nicenames %>%
    pcl.filter.s(diagnosis %in% c("HMP1-II", "nonIBD")) %>%
    pcl.top.f(mean(x), n=10) %>%
    pcl.heatmap(meta=c("diagnosis"), annotation_colors=list(diagnosis=c(hmp2_disease_colors, c("HMP1-II" = "forestgreen"))))
dev.off()




pdf("./overview/joint_ordination_hmp1ii_edfig.pdf", 1.15, 1.15)
ggp <- pcl.ordplot(merged.pcl, joint.ord, colour="diagnosis",
            colour_override=c(hmp2_disease_colors, c("HMP1-II" = "forestgreen")),
            colour_title = "Diagnosis", size_abs=1.4, outline_size=0.4) +
    theme_nature() +
    theme(axis.text.x=element_blank(), axis.text.y = element_blank())
print(ggp + guides(fill="none"))
print(ggp)
ggp <- pcl.ordplot(merged.pcl %>% pcl.nicenames, joint.ord, colour="Bacteroides ovatus",
        sortby="Bacteroides ovatus", colour_title = "B. ovatus", size_abs=1.4, outline_size=0.4) +
    theme_nature() +
    theme(axis.text.x=element_blank(), axis.text.y = element_blank())
print(ggp + guides(fill="none"))
print(ggp)
dev.off()



## Fraction with P-copri
pcopriperson <- sapply(split(pcl.nicenames(merged.pcl)$x[,"Prevotella copri"] > 0.1, merged.pcl$meta$merged_subj), any)
ct <- table(pcopriperson, diag_subj$diagnosis)
ct[2,] / colSums(ct)

# Likelihood ratio test based on binomial likelihoods for the two means being different
binoll <- function(x) {
    p <- mean(x)
    return (log(p) * sum(x) + log(1-p) * sum(!x))
}
ll_null <- binoll(pcopriperson[diag_subj$diagnosis %in% c("HMP1-II", "nonIBD")])
ll_alt <- binoll(pcopriperson[diag_subj$diagnosis == "HMP1-II"]) + binoll(pcopriperson[diag_subj$diagnosis == "nonIBD"])
D <- -2 * (ll_null - ll_alt)
p.value <- pchisq(D, 1, lower.tail=F)
p.value
# [1] 0.5497016



df <- merged.pcl.healthy$meta
df$pcopri <- merged.pcl.healthy$x[,"Prevotella copri"]

library(ggplot2)
library(cowplot)
ggplot(data=df, aes(x = VISNO, y = pcopri, color=factor(RANDSID))) +
    geom_line(aes(group=RANDSID)) +
    geom_point() +
    guides(color="none") +
    ylab("Relative abundance (Prevotella copri)") +
    xlab("Visit number") +
    #scale_y_log10() +
    theme_cowplot()



