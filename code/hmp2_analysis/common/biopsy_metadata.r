
source("./common/merge_metadata.r")

merge_biopsy_metadata <- function(pcl) {
    # Pull in the metadata fields from host_transcriptomics since biopsy_16s seems to be incomplete at this time
    mt <- merge_metadata(pcl, datatype = "host_transcriptomics", c(
        loc = "biopsy_location",
        is_inflamed = "is_inflamed"))$meta

    if (all(is.na(mt$loc))) {
        # Not biopsy data type
        return (pcl)
    }

    # Non-inflamed samples are easy to identify
    pcl$meta$is_inflamed <- mt$is_inflamed

    # Reconstruct the biopsy location from the different fields [unnecessary now]
    # mt$loc <- as.character(mt$loc)
    # mt$loc[mt$loc==""] <- NA # make blanks NA
    # pcl$meta$biopsy_location <- mt$loc
    # mask <- !is.na(mt$loc) & mt$loc == "Non-inflamed"
    # pcl$meta$biopsy_location[mask] <- as.character(mt$noninflamed_loc[mask])
    # mask <- !is.na(mt$loc) & mt$loc == "Other Inflamed"
    # pcl$meta$biopsy_location[mask] <- as.character(mt$otherinflamed_loc[mask])
    # if (any(pcl$meta$biopsy_location == "", na.rm=T)) {
    #     mask <- !is.na(pcl$meta$biopsy_location) & pcl$meta$biopsy_location == ""
    #     warning(sprintf("Samples with biopsy_locations that ended up \"\": %s", paste(pcl$meta$site_sub_coll[mask])))
    #     pcl$meta$biopsy_location[mask] <- NA
    # }
    pcl$meta$biopsy_location <- mt$biopsy_location

    # HACK!
    pcl$meta$biopsy_location[pcl$meta$biopsy_location=="Non-inflamed"] <- NA

    # Simplified location (Ileum, Colon, Rectum)
    pcl$meta$biopsy_location_ilecolrec <- as.character(pcl$meta$biopsy_location)
    pcl$meta$biopsy_location_ilecolrec[grepl("[Ii]leum", pcl$meta$biopsy_location)] <- "Ileum"
    pcl$meta$biopsy_location_ilecolrec[grepl("[Cc]olon", pcl$meta$biopsy_location)] <- "Colon"
    pcl$meta$biopsy_location_ilecolrec[pcl$meta$biopsy_location == "Cecum"] <- "Colon"
    stopifnot(length(setdiff(unique(pcl$meta$biopsy_location_ilecolrec), c(NA, "Ileum", "Colon", "Rectum"))) == 0)

    # Factors
    pcl$meta$biopsy_location <- factor(pcl$meta$biopsy_location)
    pcl$meta$biopsy_location_ilecolrec <- factor(pcl$meta$biopsy_location_ilecolrec)

    return (pcl)
}
