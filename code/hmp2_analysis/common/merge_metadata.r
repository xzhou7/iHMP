
# Read in metadata
if (!exists("hmp2_sample_metadata")) {
    source("./env_config.r")
    hmp2_sample_metadata <- read.csv(file.path(HMP2_data, "metadata", "hmp2_metadata_2018-08-20.csv"))

    # Patients that answer "Yes" to the bowel surgery question at any timepoint have the "Yes"
    # filled in for all timepoints (some flip back and forth...)
    hmp2_sample_metadata <- hmp2_sample_metadata[order(hmp2_sample_metadata$week_num),]
    for (subj in levels(hmp2_sample_metadata$Participant.ID)) {
        samps <- hmp2_sample_metadata$Participant.ID == subj
        yes <- any(hmp2_sample_metadata$X6..Have.you.ever.had.bowel.surgery.[samps] == "Yes")
        hmp2_sample_metadata$X6..Have.you.ever.had.bowel.surgery.[samps] <- ifelse(yes, "Yes", "No")
        hmp2_sample_metadata$X6..Have.you.ever.had.bowel.surgery. <- factor(hmp2_sample_metadata$X6..Have.you.ever.had.bowel.surgery., levels=c("No", "Yes"))
    }
}

merge_metadata <- function(pcl, fields, datatype, week_offset=0) {
    # Subset to a smaller set of data
    if (!missing(datatype)) {
        subs_metadata <- hmp2_sample_metadata[hmp2_sample_metadata$data_type == datatype,,drop=F]
    } else {
        subs_metadata <- hmp2_sample_metadata
    }

    # Determine the cols to match
    merge_cols <- match(fields, colnames(subs_metadata))
    if (any(is.na(merge_cols))) {
        stop("Unknown metadata fields.")
    }

	if (week_offset != 0) {
		# Time-shifted match
	    merge_df <- subs_metadata[rep("", nrow(pcl$meta)), fields, drop=F]
	    rownames(merge_df) <- rownames(pcl$meta)
		for (i in seq_len(nrow(merge_df))) {
			person_metadata <- subs_metadata[subs_metadata$Participant.ID == as.character(pcl$meta$subject[i]),,drop=F]
			good <- person_metadata$week_num == pcl$meta$week_num[i] + week_offset
			if (any(good, na.rm=T)) {
			    sdf <- person_metadata[good, fields, drop=F]
			    if (nrow(sdf) > 1) {
			        sdf <- colwise(function(x) if(all(is.na(x))) {NA} else {na.omit(x)[1]})(sdf)
			    }
			    merge_df[i,] <- sdf
			}
		}
	} else {
		# Match by External ID
		samp_match <- match(rownames(pcl$meta), subs_metadata$External.ID)

		# Match my site/sub/coll if the External IDs fail
		if (any(is.na(samp_match))) {
			samp_match <- match(pcl$meta$site_sub_coll, subs_metadata$site_sub_coll)

			if (any(is.na(samp_match))) {
				warning(sprintf("Matching for %s incomplete, even using Site/Sub/Coll fallback!!", deparse(substitute(pcl))))
			} else {
				warning(sprintf("Matching for %s done by Site/Sub/Coll since at least some ExternalIDs failed to match.", deparse(substitute(pcl))))
			}
		}

		# Pull out metadata
		merge_df <- subs_metadata[samp_match, fields, drop=F]
	}

    # Merge fields
    if (!is.null(names(fields))) {
        named <- names(fields) != ""
        colnames(merge_df)[named] <- names(fields)[named]
    }
    for (field in colnames(merge_df)) {
        pcl$meta[,field] <- merge_df[,field]
    }

    return (pcl)
}

