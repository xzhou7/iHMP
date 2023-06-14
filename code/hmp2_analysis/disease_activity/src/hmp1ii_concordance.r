
library(plyr)
library(dplyr)

source("./common/disease_colors.r")
source("./common/load_bugs.r")
source("./common/disease_activity.r")

hmp12_bugs.pcl <- pcl.read(file.path(HMP2_root, "../HMP/data/metaphlan2/hmp1-II_metaphlan2-mtd-qcd.tsv"), metadata.rows=8) %>%
    pcl.filter.s(STSite == "Stool") %>%
    pcl.filter.s(any(x>0))

# Merge the datasets
hmp12_bugs.pcl$meta$diagnosis <- "HMP1-II"
hmp12_bugs.pcl$meta$subject <- as.character(hmp12_bugs.pcl$meta$RANDSID)
merged.pcl <- pcl.merge(bugs.pcl, hmp12_bugs.pcl) %>%
    pcl.only(rank="s")
merged.pcl$x[is.na(merged.pcl$x)] <- 0

identify_active_hmp1ii <- function(pcl) {
    library(vegan)
    D <- as.matrix(vegdist(pcl$x, method="bray"))

    # The "reference set" of HMP1-II
    ref_set <- pcl$meta$diagnosis == "HMP1-II"

    # Calculate the dysbiosis score
    dysbiosis_score <- sapply(seq_along(ref_set), function(i)
        median(D[i, ref_set & (pcl$meta$subject != pcl$meta$subject[i])]))
    names(dysbiosis_score) <- rownames(pcl$meta)

    # Threshold activity
    dysbiosis_threshold <- quantile(dysbiosis_score[pcl$meta$diagnosis=="nonIBD"], 0.9)
    dysbiotic <- dysbiosis_score >= dysbiosis_threshold

    return (list(active = dysbiotic, score = dysbiosis_score, thresh = dysbiosis_threshold))
}

dysbiotic_hmp1ii <- identify_active_hmp1ii(merged.pcl)
dysbiotic_hmp2nonibd <- list(active=merged.pcl$meta$active,
                             score=merged.pcl$meta$activity_index,
                             thresh=disease_activity_threshold)

library(ggplot2)

df <- data.frame(
    ishmp2 = merged.pcl$meta$diagnosis != "HMP1-II",
    hmp2_score = dysbiotic_hmp2nonibd$score,
    hmp2_dysbiotic = dysbiotic_hmp2nonibd$active,
    hmp1ii_score = dysbiotic_hmp1ii$score,
    hmp1ii_dysbiotic = dysbiotic_hmp1ii$active,
    diagnosis = merged.pcl$meta$diagnosis
)
df <- df[df$ishmp2,]


# Correlation
cor.test(df$hmp2_score, df$hmp1ii_score, method="spearman")
# Spearman's rank correlation rho
#
# data:  df$hmp2_score and df$hmp1ii_score
# S = 128930000, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
# rho
# 0.8180527

cor.test(df$hmp2_score, df$hmp1ii_score)
# Pearson's product-moment correlation
#
# data:  df$hmp2_score and df$hmp1ii_score
# t = 67.868, df = 1618, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.8470419 0.8724077
# sample estimates:
#       cor
# 0.8602562

# Confusion matrices
table(ifelse(df$hmp2_dysbiotic, "D2", "N2"), ifelse(df$hmp1ii_dysbiotic, "D1", "N1"))
#      D1   N1
# D2  215   57
# N2   84 1239


fisher.test(table(df$hmp2_dysbiotic, df$hmp1ii_dysbiotic))
# Fisher's Exact Test for Count Data
#
# data:  table(df$hmp2_dysbiotic, df$hmp1ii_dysbiotic)
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  38.66460 82.15876
# sample estimates:
# odds ratio
#   56.06572

# Subset to disease groups
grp <- df$diagnosis == "nonIBD"
table(df$hmp2_dysbiotic[grp], df$hmp1ii_dysbiotic[grp])
#         FALSE TRUE
#   FALSE   368   15
#   TRUE     15   28

grp <- df$diagnosis == "UC"
table(df$hmp2_dysbiotic[grp], df$hmp1ii_dysbiotic[grp])
#         FALSE TRUE
#   FALSE   361   32
#   TRUE     14   45

grp <- df$diagnosis == "CD"
table(df$hmp2_dysbiotic[grp], df$hmp1ii_dysbiotic[grp])
#         FALSE TRUE
#   FALSE   518   40
#   TRUE     29  155


# Plot scores vs each other
source("./common/disease_colors.r")
source("./common/theme_nature.r")
ggp <- ggplot(data=df) +
    geom_point(aes(x=hmp2_score, y=hmp1ii_score, color=diagnosis), size=0.5) +
    geom_hline(yintercept=dysbiotic_hmp1ii$thresh, size=.5) +
    geom_vline(xintercept=dysbiotic_hmp2nonibd$thresh, size=.5) +
    scale_color_manual(values=hmp2_disease_colors) +
    theme_nature() +
    xlab("Dysbiosis score (non-IBD reference)") +
    ylab("Dysbiosis score (HMP1-II reference)")
ggp

pdf("./disease_activity/hmp1ii_concordance.pdf", 2.4, 2.4)
print(ggp + guides(color="none"))
print(ggp)
dev.off()


# Fraction of HMP2 non-IBD which are "dysbiotic" according to HMP1-II definition
mean(df$hmp1ii_dysbiotic[df$diagnosis=="nonIBD"])
# [1] 0.100939
table(df$hmp1ii_dysbiotic[df$diagnosis=="nonIBD"])
# FALSE  TRUE
# 383    43
