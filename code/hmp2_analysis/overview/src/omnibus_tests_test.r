
library(plyr)
library(dplyr)

# Test dataset for PERMANOVA with repeat measures

# Dataset characteristics
Nsamp <- 1600 # Samples
Nsubj <- 100  # Subjects
Nfeat <- 15   # Features in the dataset
varmult <- 1  # Global variance multiplier
bugvar <- 0.2 # Variance of feature means
# Deeper explanation:
# varmult is multiplied into all the metadata variances. It does not alter the
# % variance explained, but can exacerbate the difficulties of measuring this
# given the limitation to the dynamic range of the differences between samples
# imposed by the exponentiation + TSS.
# bugvar is the variance between feature-level means. Higher variance = less
# features are "important" in the bug matrix

# Generate per-sample metadata
mtdat <- data.frame(
    # <name>_<# levels>rand<balance>
    a_2rand10 = factor((runif(Nsamp) < 0.1) + 1),
    b_2rand30 = factor((runif(Nsamp) < 0.3) + 1),
    c_2rand50 = factor((runif(Nsamp) < 0.5) + 1),
    d_2rand10 = factor((runif(Nsamp) < 0.1) + 1),
    e_2rand30 = factor((runif(Nsamp) < 0.3) + 1),
    f_2rand50 = factor((runif(Nsamp) < 0.5) + 1),
    subject_id = factor(sample(1:Nsubj, Nsamp, replace=T))
)

# Generate per-subject metadata
mtdat$g_subj_3rand <- factor(sample(1:3, Nsubj, replace=T)[mtdat$subject_id])
mtdat$h_subj_4rand <- factor(sample(1:4, Nsubj, replace=T)[mtdat$subject_id])
mtdat$i_subj_5rand <- factor(sample(1:5, Nsubj, replace=T)[mtdat$subject_id])

# "True" variance explained in the underlying basis space
# Set the ones that should come out as significant to positive values
variances <- data.frame(
    a_2rand10 = 0,
    b_2rand30 = 0,
    c_2rand50 = 0.1,
    d_2rand10 = 0,
    e_2rand30 = 0,
    f_2rand50 = 0,
    g_subj_3rand = 0,
    h_subj_4rand = 0.1,
    i_subj_5rand = 0,
    subject_id = 0.2, # Extra variation between subjects
    residual = 0.2    # Remaining variation between samples not explainable by the above
)

# Basis generation
mtdat <- mtdat[,colnames(variances)[1:ncol(variances)-1]]
test <- sqrt(variances$residual) * matrix(rnorm(Nfeat*Nsamp), Nsamp, Nfeat) %>%
    sweep(., 2, sqrt(bugvar) * matrix(rnorm(Nfeat), 1, Nfeat), FUN="+")
for (var in colnames(mtdat)) {
    test <- test + sqrt(variances[,var]) *
        matrix(rnorm(Nfeat*length(levels(mtdat[1,var]))),
               length(levels(mtdat[1,var])), Nfeat)[mtdat[,var],]
}

# Exponentiation and total sum scaling (TSS) to simulate compositional measurement
test <- exp(varmult * test)
test <- sweep(test, 1, rowSums(test), FUN="/")

# Visualization
#source("./common/pcl_utils.r")
#pcl.heatmap(pcl.make(test))

# Naive testing with adonis
# (p-values for per-subject metadata are anti-conservative)
library(vegan)
D <- vegdist(test, "bray")
adonis(D ~ ., data=mtdat, permutations=199)
# Per-subject metadata (*_subj_*) that should explain no variance
# are likely to come out significant.

# These are the "true" R2's
# Actual measures will differ slightly due to exponentiation and TSS
variances / sum(variances)

# PERMANOVA with repeat measure-aware permutations
# p-values for per-subject metadata are no longer anti-conservative
mtdat_subj <- mtdat[order(mtdat$subject_id),]
mtdat_subj <- mtdat_subj[!duplicated(mtdat_subj$subject_id),grepl("subj_", colnames(mtdat_subj)), drop=F]
rownames(mtdat_subj) <- mtdat_subj$subject_id
PERMANOVA_repeat_measures(
    D = D, permutations=199,
    permute_within=mtdat[,!grepl("subj", colnames(mtdat)), drop=F],
    blocks=mtdat$subject_id, block_data=mtdat_subj)

