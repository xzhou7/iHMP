
library(plyr)
library(dplyr)

source("./common/match_datasets.r")
source("./common/load_bugs.r")

if (!exists("disease_activity_threshold")) {
    # Define active disease from the species-level Bray-Curtis dissimilarity matrix
    library(vegan)
    D <- as.matrix(vegdist((bugs.pcl %>%
        pcl.only(rank="s") %>% pcl.normalize)$x, method="bray"))

    # The "reference set" of inactive samples
    ref_set <- (bugs.pcl$meta$diagnosis == "nonIBD") & (bugs.pcl$meta$week_num >= 20)

    # Calculate the activity index
    bugs.pcl$meta$activity_index <- sapply(seq_along(ref_set), function(i)
        median(D[i, ref_set & (bugs.pcl$meta$subject != bugs.pcl$meta$subject[i])]))

    # Clean up
    rm(D, ref_set)

    # Threshold activity
    disease_activity_threshold <- quantile(bugs.pcl$meta$activity_index[bugs.pcl$meta$diagnosis=="nonIBD"], 0.9)
    eubiosis_lower_threshold <- quantile(bugs.pcl$meta$activity_index[bugs.pcl$meta$diagnosis=="nonIBD"], 0.1)
    bugs.pcl$meta$active <- bugs.pcl$meta$activity_index >= disease_activity_threshold

    # Function to merge activity from bugs.pcl into any other pcl
    merge_disease_activity <- function(pcl, lenience=0, week_offset=0) {
        mt <- match_datasets(list(pcl.filter.f(bugs.pcl, F), pcl), lenience, week_offset=week_offset)

        pcl$meta$activity_index <- NA
        pcl$meta$activity_index[match(rownames(mt[[2]]$x), rownames(pcl$meta))] <- mt[[1]]$meta$activity_index
        pcl$meta$active <- pcl$meta$activity_index >= disease_activity_threshold
        return (pcl)
    }
}
