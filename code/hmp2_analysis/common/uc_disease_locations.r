
source("./common/merge_metadata.r")

if (!exists("uc_disease_locations_master")) {
    mt <- hmp2_sample_metadata[(hmp2_sample_metadata$Location   != "") | (hmp2_sample_metadata$Location.1 != "") |
                               (hmp2_sample_metadata$Location.2 != "") | (hmp2_sample_metadata$Location.3 != "") |
                               (hmp2_sample_metadata$Location.4 != ""),]

    # Gather site-specific disease severity scores
    gather <- function(loc) {
        sev <- factor(rep("Normal or inactive disease", nrow(mt)), levels=c(
            "Normal or inactive disease",
            "Mild disease (erythema, decreased vascular pattern, mild friability)",
            "Moderate disease (marked erythema, absent vascular pattern, friability, erosions)",
            "Severe disease (spontaneous bleeding, ulceration)"))
        sev[mt$Location   == loc] <- mt$What.is.the.endoscopic.grading.of.severity.  [mt$Location   == loc]
        sev[mt$Location.1 == loc] <- mt$What.is.the.endoscopic.grading.of.severity..1[mt$Location.1 == loc]
        sev[mt$Location.2 == loc] <- mt$What.is.the.endoscopic.grading.of.severity..2[mt$Location.2 == loc]
        sev[mt$Location.3 == loc] <- mt$What.is.the.endoscopic.grading.of.severity..3[mt$Location.3 == loc]
        sev[mt$Location.4 == loc] <- mt$What.is.the.endoscopic.grading.of.severity..4[mt$Location.4 == loc]
        return (sev)
    }
    rectum <- gather("Rectum")
    lcolon <- gather("Left Colon")
    tcolon <- gather("Transverse Colon")
    rcolon <- gather("Right Colon")
    ileum  <- gather("Ileum")

    # Build the final location from the scores
    uc_disease_locations_master <- factor(rep(NA, nrow(mt)), levels=c("Proctosigmoiditis", "Left-sided colitis", "Pancolitis"))
    uc_disease_locations_master[!is.na(rcolon)] <- "Pancolitis"
    uc_disease_locations_master[(ileum == "Normal or inactive disease") & (rcolon == "Normal or inactive disease") &
               (tcolon == "Normal or inactive disease")] <- "Left-sided colitis"
    uc_disease_locations_master[(ileum == "Normal or inactive disease") & (rcolon == "Normal or inactive disease") &
               (tcolon == "Normal or inactive disease") & (lcolon == "Normal or inactive disease")] <- "Proctosigmoiditis"

    # Fix the CD guy who started diagnosed as UC
    uc_disease_locations_master[mt$diagnosis == "CD"] <- NA

    # Take the most complete set
    names(uc_disease_locations_master) <- mt$Participant.ID
    uc_disease_locations_master <- uc_disease_locations_master[!is.na(uc_disease_locations_master)]
    uc_disease_locations_master <- uc_disease_locations_master[!duplicated(names(uc_disease_locations_master))]

    # Clean up
    rm(mt, gather, rectum, lcolon, tcolon, rcolon, ileum)

    merge_uc_disease_locations <- function(pcl) {
        # Match with the dataset
        pcl$meta$uc_baseline_location <- uc_disease_locations_master[match(pcl$meta$subject, names(uc_disease_locations_master))]

        # Warn when not all get matched...
        if (any(is.na(pcl$meta$uc_baseline_location[pcl$meta$diagnosis=="UC" & !is.na(pcl$meta$diagnosis)]))) {
            ucid_fail <- unique(pcl$meta$subject[pcl$meta$diagnosis=="UC" & is.na(pcl$meta$uc_baseline_location)])
            warning(sprintf("Failed to find baseline UC locations for: %s",
                            paste(ucid_fail, collapse=" ")))
        }

        return (pcl)
    }
}
