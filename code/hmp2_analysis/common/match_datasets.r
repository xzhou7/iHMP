
source("./common/nonredundant_samples.r")

match_datasets <- function(pcl_list, lenience=0, week_offset=0, matching=F) {
    # pcl_list is a list of PCL structures to match
    # lenience is the maximum difference in week_num's to allow between matched samples
    # week_offset is a time-shift (in weeks) to apply to all datasets after the first,
    #             i.e. week x of dataset 1 is paired with week x+week_offset of the
    #             other datasets
    # if matching is T, then a matrix is returned, containing the sample names
    #                   from each dataset to use, otherwise a list of filtered
    #                   and reordered datasets is returned

    # Small number used to resolve duplicate week_nums
    small <- 0.001

    # Use only non-redundant samples from each dataset
    pcl_list <- lapply(pcl_list, nonredundant_samples)

    if (lenience == 0 && week_offset == 0) {
        # For a strict matching, match site_sub_colls directly
        sscs <- unique(pcl_list[[1]]$meta$site_sub_coll)
        for (i in seq_along(pcl_list))
            sscs <- intersect(sscs, unique(pcl_list[[i]]$meta$site_sub_coll))

        mt <- matrix("", ncol=length(pcl_list), nrow=length(sscs))
        for (i in seq_along(pcl_list))
            mt[,i] <- rownames(pcl_list[[i]]$meta)[match(sscs, pcl_list[[i]]$meta$site_sub_coll)]

        if (matching) {
            return (unname(mt))
        } else {
            ret <- list()
            for (i in seq_len(length(pcl_list)))
                ret[[i]] <- pcl.reorder.s(pcl_list[[i]], mt[,i])
            return (ret)
        }
    }

    # Use only samples which have week_num's
    pcl_list <- lapply(pcl_list, function(x)pcl.filter.s(x, !is.na(week_num)))

    # Get a list of all subjects with at least one data point for each pcl
    subjects <- unique(do.call(c, lapply(pcl_list, function(x) as.character(x$meta$subject))))
    subjects <- subjects[!is.na(subjects)]

    # Find the optimal matching per-subject
    mt <- matrix("", ncol=length(pcl_list), nrow=0)
    for (subj in subjects) {
        # Get the week_nums for each dataset
        week_nums <- lapply(pcl_list, function(pcl) {
            # Pull out the subject's samples
            I <- which(pcl$meta$subject == subj)

            # Order by collection number
            I <- I[order(pcl$meta$collection[I])]

            # Get week nums
            wk <- pcl$meta$week_num[I]

            # Offset later duplicate week nums by a tiny amount to eliminate ties
            wk_off <- wk + small*duplicated(wk)

            # Put the sample names in the names
            names(wk_off) <- rownames(pcl$meta)[I]

            return (wk_off)
        })

        week_nums[[1]] <- week_nums[[1]] + week_offset

        # Iteratively find the first point where all datatypes exist and
        # greedily add them to the set of global timepoints
        tp <- matrix("", ncol=length(pcl_list), nrow=0)
        while (all(sapply(week_nums, length) > 0)) {
            # Can we make a global timepoint out of the front set of week_nums?
            front <- sapply(week_nums, function(wk)wk[1])
            mf <- max(front)
            nf <- min(front)

            if (mf - nf <= lenience + 1.5*small) {
                # What samples might go into this timepoint?
                min_eligible <- mf + 0.5*small
                elig_wk <- lapply(week_nums, function(wk)wk[wk < min_eligible])

                # Select samples within the window such that we minimize the
                # std of the week_nums
                best_samps <- names(front)
                best_samp_sd <- sd(front)
                I <- rep(1, length(front))
                done <- F
                while (T) {
                    for (i in seq_along(I)) {
                        if (I[i] < length(elig_wk[[i]])) {
                            I[i] <- I[i] + 1
                            break
                        } else if (i == length(I)) {
                            done <- T
                            break
                        } else {
                            I[i] <- 1
                        }
                    }

                    if (done)
                        break

                    samp_wk <- sapply(seq_along(I), function(i) elig_wk[[i]][I[i]])
                    if (sd(samp_wk) < best_samp_sd) {
                        best_samp_sd <- sd(samp_wk)
                        best_samps <- names(samp_wk)
                    }
                }

                # Append the timepoint
                tp <- rbind(tp, best_samps)

            } else {
                # Discard samples up to the highest lenience boundary we can
                # accept based on the latest sample week
                min_eligible <- mf - (lenience + 1.5*small)
            }

            # Discard ineligible samples and continue
            week_nums <- lapply(week_nums, function(wk)wk[wk >= min_eligible])
        }

        # Append the samples to the list of matches
        mt <- rbind(mt, tp)
    }

    if (matching) {
        return (unname(mt))
    } else {
        ret <- list()
        for (i in seq_len(length(pcl_list)))
            ret[[i]] <- pcl.reorder.s(pcl_list[[i]], mt[,i])
        return (ret)
    }
}


# Test

#source("./common/load_bugs.r")
#mtt <- match_datasets(list(bugs.pcl, viruses.pcl), lenience=0, matching=T)
#head(mtt)
#mtt <- match_datasets(list(bugs.pcl, viruses.pcl), lenience=2, matching=T)
#head(mtt)
