
library(plyr)
library(dplyr)

nonredundant_samples <- function(pcl) {
    # Which samples are labelled as techreps?
    droptechrep <- grepl("_TR", rownames(pcl$x))

    # Which pilot samples are duplicates of other samples?
    sscpilot <- pcl$meta$site_sub_coll[pcl$meta$pilot, drop=F]
    sscpilot_dup <- sscpilot[sscpilot %in% pcl$meta$site_sub_coll[!pcl$meta$pilot, drop=F]]
    droppilot <- (pcl$meta$site_sub_coll %in% sscpilot_dup) & pcl$meta$pilot

    # Filter replicates out
    pcl <- pcl %>%
        pcl.filter.s(keep=!(droppilot | droptechrep))

    return (pcl)
}


