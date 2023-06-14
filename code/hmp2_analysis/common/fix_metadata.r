
fix_metadata <- function(pcl) {
    # Fixes common problems in the metadata

    if (!("site_sub_coll" %in% colnames(pcl$meta))) {
        warning("No site_sub_coll! Assuming sample names are site_sub_coll...")
        pcl$meta$site_sub_coll <- rownames(pcl$meta)
    }

    # Check for NA SSC's
    if (any(is.na(pcl$meta$site_sub_coll))) {
        warning("Table contains NA site_sub_colls (will cause more warnings).")
    }

    # Fix metadata loaded as factors that should be numeric
    if (is.factor(pcl$meta$consent_age)) {
        warning("consent_age loaded as a factor. There are string elements in this column!")
        pcl$meta$consent_age <- as.numeric(levels(pcl$meta$consent_age))[as.numeric(pcl$meta$consent_age)]
    }
    if (is.factor(pcl$meta$hbi)) {
        warning("hbi loaded as a factor. There are string elements in this column!")
        pcl$meta$hbi <- as.numeric(levels(pcl$meta$hbi))[as.numeric(pcl$meta$hbi)]
    }
    if (is.factor(pcl$meta$sccai)) {
        warning("sccai loaded as a factor. There are string elements in this column!")
        pcl$meta$sccai <- as.numeric(levels(pcl$meta$sccai))[as.numeric(pcl$meta$sccai)]
    }
    if (is.factor(pcl$meta$week_num)) {
        warning("week_num loaded as a factor. There are string elements in this column!")
        pcl$meta$week_num <- as.numeric(levels(pcl$meta$week_num))[as.numeric(pcl$meta$week_num)]
    }
    if (any(is.na(pcl$meta$week_num))) {
        warning("NAs in week_num. Imputing 0..")
        pcl$meta$week_num[is.na(pcl$meta$week_num)] <- 0
    }

    # Make a subject ID for convenience
    pcl$meta$subject <- as.factor(gsub("^(\\w\\d\\d\\d\\d).*$", "\\1", pcl$meta$site_sub_coll))

    # Make a collection number for convenience
    coll <- suppressWarnings(as.numeric(gsub("^\\w\\d\\d\\d\\dC(\\d\\d?)$", "\\1", pcl$meta$site_sub_coll)))
    if (any(!is.na(coll))) {
        pcl$meta$collection <- coll
        if (!all(!is.na(coll))) {
            warning("Failed to extract collection numbers for all samples.")
        }
    } else {
        # Try a pseudo-collection number from baseline/follow-up SSC's
        coll <- suppressWarnings(as.numeric(gsub("^\\w\\d\\d\\d\\dC?(BL|FU)(\\d\\d?)$", "\\2", pcl$meta$site_sub_coll)))
        if (any(!is.na(coll))) {
            # Blood
            pcl$meta$collection <- ifelse(coll==1, 1, coll * 2)
            if (!all(!is.na(coll))) {
                warning("Failed to extract collection numbers for all samples.")
            }
        } else {
            # Biopsies
            pcl$meta$collection <- (pcl$meta$week_num + 2) / 2
        }
    }

    # Other -> nonIBD
    if ("diagnosis" %in% colnames(pcl$meta)) {
        if (!all(na.omit(pcl$meta$diagnosis) %in% c("nonIBD", "UC", "CD"))) {
            warning("Diagnosis contains elements not equal to nonIBD/UC/CD.")
            pcl$meta$diagnosis[pcl$meta$diagnosis=="Other"] <- "nonIBD"
        }

        # Order all plots as nonIBD -> UC -> CD
        pcl$meta$diagnosis <- factor(pcl$meta$diagnosis, levels=c("nonIBD", "UC", "CD"))
    }

    # Ensure sex is Male/Female
    if ("sex" %in% colnames(pcl$meta)) {
        if (!all(na.omit(pcl$meta$sex) %in% c("Male", "Female"))) {
            warning("Sex contains elements not equal to Male/Female.")
            pcl$meta$sex[!(pcl$meta$sex %in% c("Male", "Female"))] <- NA
        }
    }

    # Fix factors for any *counts fields (for read counts)
    for (field in colnames(pcl$meta)[grepl("counts$", colnames(pcl$meta)) | (colnames(pcl$meta) %in% c("Raw", "Mapped", "Unmapped"))]) {
        if (is.factor(pcl$meta[,field])) {
            warning(sprintf("%s loaded as a factor. There are string elements in this column!", field))
            pcl$meta[,field] <- as.numeric(levels(pcl$meta[,field]))[as.numeric(pcl$meta[,field])]

            # Counts should not have NA's
            if (any(is.na(pwy.dna.pcl$filtered_reads))) {
                warning(sprintf("There are NA %s.", field))
            }
        }
    }

    # Check for things that shouldn't be NA
    for (field in c("consent_age", "diagnosis", "week_num", "subject", "race", "site_name")) {
        if (field %in% colnames(pcl$meta) && any(is.na(pcl$meta[,field]))) {
            if (field == "consent_age" && all(pcl$meta$subject[is.na(pcl$meta[,field])] %in% c("C3031", "C3007", "M2059", "E5023", "M2082", "E5006"))) {
                # These subjects don't have a consent age, so suppress the warning...
            } else {
                warning(sprintf("There are NAs in %s.", field))
            }
        }
    }

    # Really high HBI scores are considered NA (coded as 999 in the data)
    if ("hbi" %in% colnames(pcl$meta)) {
        if (any(pcl$meta$hbi > 30, na.rm=T)) {
            warning(sprintf("Unusually high HBI scores of %d. Recoding as NA.", max(pcl$meta$hbi, na.rm=T)))
            pcl$meta$hbi[!is.na(pcl$meta$hbi) & (pcl$meta$hbi > 30)] <- NA
        }
    } else {
        warning("HBI is not in the metadata!")
    }

    # Really high SCCAI scores are considered NA (coded as 999 in the data)
    if ("sccai" %in% colnames(pcl$meta)) {
        if (any(pcl$meta$sccai > 100, na.rm=T)) {
            warning(sprintf("Unusually high SCCAI scores of %d. Recoding as NA.", max(pcl$meta$sccai, na.rm=T)))
            pcl$meta$sccai[!is.na(pcl$meta$sccai) & (pcl$meta$sccai > 100)] <- NA
        }
    } else {
        warning("SCCAI is not in the metadata!")
    }

    # Remove the age for the guy who's literally 1000 years old
    if ("consent_age" %in% colnames(pcl$meta)) {
        if (any(pcl$meta$consent_age > 100, na.rm=T)) {
            warning(sprintf("Someone is %d years old!! Recoding as NA.", max(pcl$meta$consent_age, na.rm=T)))
            pcl$meta$consent_age[pcl$meta$consent_age>100] <- NA
        }

        # Add an "adult" column
        pcl$meta$adult <- pcl$meta$consent_age >= 18
    }

    # Change "Massachusetts General Hospital" to "MGH"
    if ("site_name" %in% colnames(pcl$meta)) {
        if (any(levels(pcl$meta$site_name)=="Massachusetts General Hospital")) {
            warning("site_name coded with 'Massachusetts General Hospital'")
            levels(pcl$meta$site_name)[levels(pcl$meta$site_name)=="Massachusetts General Hospital"] <- "MGH"
        }
    }

    # Add a "pilot" field where pilot samples are
    pcl$meta$pilot <- grepl("^.*_P$", rownames(pcl$meta))

    return (pcl)
}

