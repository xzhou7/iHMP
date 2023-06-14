
source("./common/merge_metadata.r")


mdsum_site <- function(name, mt, sites, ord=NA) {
    # mdsum <- function(mt) {
    #     nna <- sum(is.na(mt))
    #     mt <- na.omit(mt)
    #
    #     if (is.logical(mt)) {
    #         s <- sprintf("%d", sum(mt))
    #     } else if (is.factor(mt)) {
    #         lvs <- as.character(unique(mt))
    #         if (length(ord)==1 && is.na(ord)) {
    #             lvs <- sort(lvs)
    #         } else {
    #             stopifnot(all(lvs %in% ord))
    #             lvs <- ord[ord %in% lvs]
    #         }
    #         s <- ""
    #         for (lv in lvs) {
    #             if (s != "") {
    #                 s <- sprintf("%s, ", s)
    #             }
    #             s <- sprintf("%s%s (%d)", s, lv, sum(mt==lv))
    #         }
    #
    #         if (nna > 0) {
    #             s <- sprintf("%s, N/A (%d)", s, nna)
    #         }
    #         return (s)
    #     } else if (length(mt) == 0) {
    #         s <- "Empty"
    #     } else if (all((mt %% 1 < 0.001) | (mt %% 1 > 0.999))) {
    #         s <- sprintf("%.2f ± %.2f [%d, %d]", mean(mt), sd(mt), min(mt), max(mt))
    #     } else {
    #         s <- sprintf("%.2f ± %.2f [%.2g, %.2g]", mean(mt), sd(mt), min(mt), max(mt))
    #     }
    #
    #     if (nna > 0) {
    #         s <- sprintf("%s, %d N/A", s, nna)
    #     }
    #     return (s)
    # }

    mdsum <- function(mt) {
        nna <- sum(is.na(mt))
        mt <- na.omit(mt)

        if (is.logical(mt)) {
            dat <- c(sum(mt), nna)
            names(dat) <- c(name, "N/A")
            return (dat)
        } else if (is.factor(mt)) {
            lvs <- levels(mt)
            if (length(ord)==1 && is.na(ord)) {
                lvs <- sort(lvs)
            } else {
                stopifnot(all(lvs %in% ord))
                lvs <- ord[ord %in% lvs]
            }
            dat <- rep(0, length(lvs)+1)
            for (i in seq_along(lvs)) {
                dat[i] <- sum(mt == lvs[i])
            }
            dat[length(lvs)+1] <- nna

            names(dat) <- c(lvs, "N/A")
            return (dat)
        } else if (length(mt) == 0) {
            dat <- "Empty"
            names(dat) <- name
        } else if (all((mt %% 1 < 0.001) | (mt %% 1 > 0.999))) {
            dat <- c(sprintf("%.1f ± %.1f", mean(mt), sd(mt)),
                     sprintf("[%d - %d]", min(mt), max(mt)))
            names(dat) <- c(name, "Range")
        } else {
            dat <- c(sprintf("%.1f ± %.1f", mean(mt), sd(mt)),
                     sprintf("[%.1f - %.1f]", min(mt), max(mt)))
            names(dat) <- c(name, "Range")
        }

        return (dat)
    }

    rdat <- data.frame(
        "Cedars-Sinai"   = mdsum(mt[sites=="Cedars-Sinai"]),
        Emory            = mdsum(mt[sites=="Emory"]),
        Cincinnati       = mdsum(mt[sites=="Cincinnati"]),
        MGH              = mdsum(mt[sites=="MGH"]),
        "MGH Pediatrics" = mdsum(mt[sites=="MGH Pediatrics"]),
        Total            = mdsum(mt)
        #row.names        = name
    )

    if (is.numeric(rdat$Total)) {
        rdat <- rdat[rowSums(rdat) > 0,,drop=F]
        rn <- rownames(rdat)
        rdat <- colwise(function(x) ifelse(x>0,as.character(x),""))(rdat)
        rownames(rdat) <- rn
    }

    return (rdat)
}

sub_mtdat <- hmp2_sample_metadata[hmp2_sample_metadata$data_type != "noproduct",]
sub_mtdat <- sub_mtdat[sub_mtdat$Participant.ID != "",]
sub_mtdat <- do.call(rbind, lapply(split(sub_mtdat, factor(sub_mtdat$Participant.ID)),
    function(x) colwise(function(q){
        q <- na.omit(q)
        q <- q[q != ""]
        if (length(q)>0) {
            return (q[1])
        } else {
            return (NA)
        }
    })(x)))
sites <- sub_mtdat$site_name
sub_mtdat$consent_age[sub_mtdat$consent_age>100] <- NA
sub_mtdat$Weight[sub_mtdat$Weight>500] <- NA
sub_mtdat$Height[sub_mtdat$Height>500] <- NA
separator <- data.frame("Cedars-Sinai" = "", Emory = "", Cincinnati = "", MGH = "", "MGH Pediatrics" = "", Total = "")
epi_table <- rbind(
    mdsum_site("Subjects", rep(T, length(sites)), sites),
    separator,
    mdsum_site("Diagnosis", sub_mtdat$diagnosis, sites, ord=c("nonIBD", "UC", "CD")),
    separator,
    mdsum_site("Sex", sub_mtdat$sex, sites),
    separator,
    mdsum_site("Ethnicity", sub_mtdat$race, sites),
    separator,
    mdsum_site("Age (years)", sub_mtdat$consent_age, sites),
    separator,
    mdsum_site("Height (cm)", sub_mtdat$Height, sites),
    separator,
    mdsum_site("Weight (kg)", sub_mtdat$Weight, sites),
    separator,
    mdsum_site("BMI", sub_mtdat$BMI, sites),
    separator,
    mdsum_site("CD Montreal Location", sub_mtdat$baseline_montreal_location[sub_mtdat$diagnosis=="CD"], sites[sub_mtdat$diagnosis=="CD"]),
    separator,
    mdsum_site("UC Extent", factor(sub_mtdat$Extent..E..1[sub_mtdat$diagnosis=="UC"]), sites[sub_mtdat$diagnosis=="UC"])
)

epi_table

write.table(epi_table, file="./overview/epi_table.tsv",
            sep = "\t", quote=F, row.names=T, col.names=T)

