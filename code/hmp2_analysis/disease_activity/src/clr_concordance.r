
library(plyr)
library(dplyr)

source("./common/disease_colors.r")
source("./common/load_bugs.r")
source("./common/disease_activity.r")

# The "reference set" of inactive samples
ref_set <- (bugs.pcl$meta$diagnosis == "nonIBD") & (bugs.pcl$meta$week_num >= 20)



#####  CLR-overall

# Load the CLR based distance matrix
load(file.path(HMP2_data, "mgx", "species_CLR_overall.RData"))
clrD <- as.matrix(dist(species_CLR_overall.pcl$x))

# Calculate the CLR-based activity index
bugs.pcl$meta$activity_index_clr <- sapply(seq_along(ref_set), function(i)
    median(clrD[i, ref_set & (bugs.pcl$meta$subject != bugs.pcl$meta$subject[i])]))

disease_activity_threshold_clr <- quantile(bugs.pcl$meta$activity_index_clr[bugs.pcl$meta$diagnosis=="nonIBD"], 0.9)

df <- data.frame(
    clr_score = bugs.pcl$meta$activity_index_clr,
    clr_dysbiotic = bugs.pcl$meta$activity_index_clr > disease_activity_threshold_clr,
    bc_score = bugs.pcl$meta$activity_index,
    bc_dysbiotic = bugs.pcl$meta$active,
    diagnosis = bugs.pcl$meta$diagnosis
)


# Correlation
cor.test(df$bc_score, df$clr_score, method="spearman")
# Spearman's rank correlation rho
#
# data:  df$bc_score and df$clr_score
# S = 366550000, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho
# 0.4579916

cor.test(df$bc_score, df$clr_score)
# Pearson's product-moment correlation
#
# data:  df$bc_score and df$clr_score
# t = 21.844, df = 1593, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.4414204 0.5170007
# sample estimates:
#       cor
# 0.4801011

# Confusion matrices
table(df$bc_dysbiotic, df$clr_dysbiotic)
# FALSE TRUE
# FALSE  1231   92
# TRUE    162  110

fisher.test(table(df$bc_dysbiotic, df$clr_dysbiotic))
# Fisher's Exact Test for Count Data
#
# data:  table(df$bc_dysbiotic, df$clr_dysbiotic)
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
# 6.495162 12.687648
# sample estimates:
# odds ratio
# 9.063362

# Subset to disease groups
grp <- df$diagnosis == "nonIBD"
table(df$bc_dysbiotic[grp], df$clr_dysbiotic[grp])
# FALSE TRUE
# FALSE   361   22
# TRUE     22   21

grp <- df$diagnosis == "UC"
table(df$bc_dysbiotic[grp], df$clr_dysbiotic[grp])
# FALSE TRUE
# FALSE   356   30
# TRUE     30   21

grp <- df$diagnosis == "CD"
table(df$bc_dysbiotic[grp], df$clr_dysbiotic[grp])
# FALSE TRUE
# FALSE   514   40
# TRUE    110   68


# Plot scores vs each other
source("./common/disease_colors.r")
source("./common/theme_nature.r")
ggp <- ggplot(data=df) +
    geom_point(aes(x=bc_score, y=clr_score, color=diagnosis), size=0.5) +
    geom_hline(yintercept=disease_activity_threshold_clr, size=.5) +
    geom_vline(xintercept=disease_activity_threshold, size=.5) +
    scale_color_manual(values=hmp2_disease_colors) +
    theme_nature() +
    xlab("Dysbiosis score (Bray-Curtis)") +
    ylab("Dysbiosis score (CLR)")

pdf("./disease_activity/clr_concordance_overall.pdf", 2.4, 2.4)
print(ggp + guides(color="none"))
print(ggp)
dev.off()



#####  CLR-perfeature

# Load the CLR based distance matrix
load(file.path(HMP2_data, "mgx", "species_CLR_perfeature.RData"))
clrD <- as.matrix(dist(species_CLR_perfeature.pcl$x))

# Calculate the CLR-based activity index
bugs.pcl$meta$activity_index_clr <- sapply(seq_along(ref_set), function(i)
    median(clrD[i, ref_set & (bugs.pcl$meta$subject != bugs.pcl$meta$subject[i])]))

disease_activity_threshold_clr <- quantile(bugs.pcl$meta$activity_index_clr[bugs.pcl$meta$diagnosis=="nonIBD"], 0.9)

df <- data.frame(
    clr_score = bugs.pcl$meta$activity_index_clr,
    clr_dysbiotic = bugs.pcl$meta$activity_index_clr > disease_activity_threshold_clr,
    bc_score = bugs.pcl$meta$activity_index,
    bc_dysbiotic = bugs.pcl$meta$active,
    diagnosis = bugs.pcl$meta$diagnosis
)


# Correlation
cor.test(df$bc_score, df$clr_score, method="spearman")
# Spearman's rank correlation rho
#
# data:  df$bc_score and df$clr_score
# S = 362740000, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#      rho
# 0.463635

cor.test(df$bc_score, df$clr_score)
# Pearson's product-moment correlation
#
# data:  df$bc_score and df$clr_score
# t = 22.623, df = 1593, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.4550421 0.5293811
# sample estimates:
#       cor
# 0.4931112

# Confusion matrices
table(df$bc_dysbiotic, df$clr_dysbiotic)
# FALSE TRUE
# FALSE  1231   92
# TRUE    160  112

fisher.test(table(df$bc_dysbiotic, df$clr_dysbiotic))
# Fisher's Exact Test for Count Data
#
# data:  table(df$bc_dysbiotic, df$clr_dysbiotic)
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   6.700169 13.075193
# sample estimates:
# odds ratio
#   9.347678

# Subset to disease groups
grp <- df$diagnosis == "nonIBD"
table(df$bc_dysbiotic[grp], df$clr_dysbiotic[grp])
# FALSE TRUE
# FALSE   361   22
# TRUE     22   21

grp <- df$diagnosis == "UC"
table(df$bc_dysbiotic[grp], df$clr_dysbiotic[grp])
# FALSE TRUE
# FALSE   355   31
# TRUE     29   22

grp <- df$diagnosis == "CD"
table(df$bc_dysbiotic[grp], df$clr_dysbiotic[grp])
# FALSE TRUE
# FALSE   515   39
# TRUE    109   69


# Plot scores vs each other
source("./common/disease_colors.r")
source("./common/theme_nature.r")
ggp <- ggplot(data=df) +
    geom_point(aes(x=bc_score, y=clr_score, color=diagnosis), size=0.5) +
    geom_hline(yintercept=disease_activity_threshold_clr, size=.5) +
    geom_vline(xintercept=disease_activity_threshold, size=.5) +
    scale_color_manual(values=hmp2_disease_colors) +
    theme_nature() +
    xlab("Dysbiosis score (Bray-Curtis)") +
    ylab("Dysbiosis score (CLR)")

pdf("./disease_activity/clr_concordance_perfeature.pdf", 2.4, 2.4)
print(ggp + guides(color="none"))
print(ggp)
dev.off()



#####  CLR-persample

# Load the CLR based distance matrix
load(file.path(HMP2_data, "mgx", "species_CLR_persample.RData"))
clrD <- as.matrix(dist(species_CLR_persample.pcl$x))

# Calculate the CLR-based activity index
bugs.pcl$meta$activity_index_clr <- sapply(seq_along(ref_set), function(i)
    median(clrD[i, ref_set & (bugs.pcl$meta$subject != bugs.pcl$meta$subject[i])]))

disease_activity_threshold_clr <- quantile(bugs.pcl$meta$activity_index_clr[bugs.pcl$meta$diagnosis=="nonIBD"], 0.9)

df <- data.frame(
    clr_score = bugs.pcl$meta$activity_index_clr,
    clr_dysbiotic = bugs.pcl$meta$activity_index_clr > disease_activity_threshold_clr,
    bc_score = bugs.pcl$meta$activity_index,
    bc_dysbiotic = bugs.pcl$meta$active,
    diagnosis = bugs.pcl$meta$diagnosis
)


# Correlation
cor.test(df$bc_score, df$clr_score, method="spearman")
# Spearman's rank correlation rho
#
# data:  df$bc_score and df$clr_score
# S = 454760000, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
# rho
# 0.3275614

cor.test(df$bc_score, df$clr_score)
# Pearson's product-moment correlation
#
# data:  df$bc_score and df$clr_score
# t = 15.394, df = 1593, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
# 0.3163640 0.4018436
# sample estimates:
# cor
# 0.3598587

# Confusion matrices
table(df$bc_dysbiotic, df$clr_dysbiotic)
# FALSE TRUE
# FALSE  1214  109
# TRUE    181   91

fisher.test(table(df$bc_dysbiotic, df$clr_dysbiotic))
# Fisher's Exact Test for Count Data
#
# data:  table(df$bc_dysbiotic, df$clr_dysbiotic)
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
# 4.010313 7.792087
# sample estimates:
# odds ratio
# 5.590791

# Subset to disease groups
grp <- df$diagnosis == "nonIBD"
table(df$bc_dysbiotic[grp], df$clr_dysbiotic[grp])
# FALSE TRUE
# FALSE   355   28
# TRUE     28   15

grp <- df$diagnosis == "UC"
table(df$bc_dysbiotic[grp], df$clr_dysbiotic[grp])
# FALSE TRUE
# FALSE   355   31
# TRUE     37   14

grp <- df$diagnosis == "CD"
table(df$bc_dysbiotic[grp], df$clr_dysbiotic[grp])
# FALSE TRUE
# FALSE   504   50
# TRUE    116   62

# What fraction of the different disease groups are dysbiotic
dia_dys <- table(df$clr_dysbiotic, df$diagnosis)
dia_dys_pct <- sweep(dia_dys, 2, colSums(dia_dys), FUN="/")
dia_dys_pct
# nonIBD        UC        CD
# FALSE 0.8990610 0.8970252 0.8469945
# TRUE  0.1009390 0.1029748 0.1530055

dia_dys <- table(df$bc_dysbiotic, df$diagnosis)
dia_dys_pct <- sweep(dia_dys, 2, colSums(dia_dys), FUN="/")
dia_dys_pct
# nonIBD        UC        CD
# FALSE 0.8990610 0.8832952 0.7568306
# TRUE  0.1009390 0.1167048 0.2431694

# Plot scores vs each other
source("./common/disease_colors.r")
source("./common/theme_nature.r")
ggp <- ggplot(data=df) +
    geom_point(aes(x=bc_score, y=clr_score, color=diagnosis), size=0.5) +
    geom_hline(yintercept=disease_activity_threshold_clr, size=.5) +
    geom_vline(xintercept=disease_activity_threshold, size=.5) +
    scale_color_manual(values=hmp2_disease_colors) +
    theme_nature() +
    xlab("Dysbiosis score (Bray-Curtis)") +
    ylab("Dysbiosis score (CLR)")

pdf("./disease_activity/clr_concordance_persample.pdf", 2.4, 2.4)
print(ggp + guides(color="none"))
print(ggp)
dev.off()



