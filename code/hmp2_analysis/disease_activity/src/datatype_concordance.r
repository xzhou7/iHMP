
library(ggplot2)

identify_active <- function(pcl) {

    library(vegan)
    D <- as.matrix(vegdist(pcl$x, method="bray"))

    # The "reference set" of inactive samples
    ref_set <- (pcl$meta$diagnosis == "nonIBD") & (pcl$meta$week_num >= 20)

    # Calculate the dysbiosis score
    dysbiosis_score <- sapply(seq_along(ref_set), function(i)
        median(D[i, ref_set & (pcl$meta$subject != pcl$meta$subject[i])]))
    names(dysbiosis_score) <- rownames(pcl$meta)

    # Threshold activity
    dysbiosis_threshold <- quantile(dysbiosis_score[pcl$meta$diagnosis=="nonIBD"], 0.9)
    dysbiotic <- dysbiosis_score >= dysbiosis_threshold

    return (dysbiotic)
}


source("./common/load_bugs.r")
source("./common/load_metabolites.r")

source("./common/match_datasets.r")

bugs.dys <- identify_active(bugs.pcl %>% pcl.only(rank="s") %>% pcl.normalize)
metabolites.dys <- identify_active(metabolites.pcl.nrm)


# Test MGX bugs - MBX
matching <- match_datasets(list(bugs.pcl, metabolites.pcl.nrm), lenience=0, matching=T)
ct <- table(ifelse(bugs.dys[matching[,1]], "YB", "NB"), ifelse(metabolites.dys[matching[,2]], "YM", "NM"))

ct
#     NM  YM
# NB 314  72
# YB  39  36

fisher.test(ct)
# Fisher's Exact Test for Count Data
#
# data:  ct
# p-value = 5.992e-09
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  2.674766 7.876497
# sample estimates:
# odds ratio
#   4.580315




