
source("./common/disease_activity.r")
source("./common/merge_metadata.r")
source("./common/biopsy_metadata.r")
source("./common/merge_BMI.r")
source("./common/uc_disease_locations.r")

# Metadata merge function from omnibus_tests.r

merge_metadata_for_omnibus <- function(pcl, act_lenience=0, week_offset=0) {
    # Merge in/create the expected metadata
    pcl <- pcl %>% merge_metadata(c(
        abx = "Antibiotics",
        immsup = "Immunosuppressants..e.g..oral.corticosteroids.",
        bowelsurgery = "X6..Have.you.ever.had.bowel.surgery."),
        week_offset=week_offset) %>%
        merge_biopsy_metadata
    pcl$meta$ped_disease <- factor(sprintf("%s%s", ifelse(pcl$meta$diagnosis=="nonIBD", "", ifelse(pcl$meta$adult, "Adult", "Pediatric")), pcl$meta$diagnosis))

    if (!("is_inflamed" %in% colnames(pcl$meta))) {
        pcl$meta$is_inflamed <- F
    }
    if (!("biopsy_location_ilecolrec" %in% colnames(pcl$meta))) {
        pcl$meta$biopsy_location_ilecolrec <- "None"
    }

    pcl <- merge_BMI(pcl)
    pcl$meta$BMI_imputed <- pcl$meta$BMI
    pcl$meta$BMI_imputed[is.na(pcl$meta$BMI_imputed)] <- mean(pcl$meta$BMI_imputed[!duplicated(pcl$meta$subject)], na.rm=T)

    pcl$meta$montreal_l123 <- factor("None", levels=c("None", "L1", "L2", "L3"))
    pcl$meta$montreal_l123[grepl("L1", pcl$meta$baseline_montreal_location)] <- "L1"
    pcl$meta$montreal_l123[grepl("L2", pcl$meta$baseline_montreal_location)] <- "L2"
    pcl$meta$montreal_l123[grepl("L3", pcl$meta$baseline_montreal_location)] <- "L3"
    pcl$meta$montreal_l123[is.na(pcl$meta$baseline_montreal_location) & !is.na(pcl$meta$diagnosis) & (pcl$meta$diagnosis == "CD")] <- NA
    pcl$meta$montreal_l4 <- factor(ifelse(grepl("L4", pcl$meta$baseline_montreal_location), "L4", "None"), levels=c("None", "L4"))
    pcl$meta$montreal_l4[is.na(pcl$meta$baseline_montreal_location) & !is.na(pcl$meta$diagnosis) & (pcl$meta$diagnosis == "CD")] <- NA

    pcl <- merge_uc_disease_locations(pcl)
    pcl$meta$uc_baseline_location_nona <- factor(pcl$meta$uc_baseline_location, levels=c("None", levels(pcl$meta$uc_baseline_location)))
    pcl$meta$uc_baseline_location_nona[is.na(pcl$meta$uc_baseline_location_nona)] <- "None"
    #pcl$meta$baseline_uc_extent_nona <- factor(pcl$meta$baseline_uc_extent, levels=c("None", levels(pcl$meta$baseline_uc_extent)))
    #pcl$meta$baseline_uc_extent_nona[is.na(pcl$meta$baseline_uc_extent_nona)] <- "None"

    pcl <- merge_disease_activity(pcl, lenience=act_lenience, week_offset=week_offset)

    pcl$meta$race_simplified <- pcl$meta$race
    pcl$meta$race_simplified[pcl$meta$race_simplified=="American Indian or Alaska Native"] <- "Other"

    pcl$meta$bowelsurgery[pcl$meta$bowelsurgery==""] <- NA

    pcl$meta$ped_disease <- sprintf("%s_%s", pcl$meta$diagnosis, pcl$meta$adult)

    return (pcl)
}

covariates <- c(
    # "Name in heatmap" = c(list of metadata names to include as covariates)
    Age = "consent_age",
    Sex = "sex",
    BMI = "BMI_imputed",
    Race = "race_simplified",
    "Bowel Surgery" = "bowelsurgery",
    "Recruitment Site" = "site_name",
    "UC Extent" = "uc_baseline_location_nona",
    "Montreal Location (L123)" = "montreal_l123",
    "Montreal Location (L4)" = "montreal_l4",
    Antibiotics = "abx",
    Immunosuppressants = "immsup"
)

time_delays <- c(-4, -2, 0, 2, 4)

assoc_beta <- matrix(NA, nrow=length(time_delays), ncol=length(covariates))
colnames(assoc_beta) <- names(covariates)
rownames(assoc_beta) <- as.character(time_delays)
assoc_p <- assoc_beta
colnames(assoc_p) <- names(covariates)
rownames(assoc_p) <- as.character(time_delays)
convergence_failed <- matrix(F, nrow=length(time_delays), ncol=length(covariates))

library(lme4)
library(car)


for (j in seq_along(time_delays)) {
    bugsmt.pcl <- merge_metadata_for_omnibus(bugs.pcl, week_offset=time_delays[j])
    for (i in seq_along(covariates)) {
        bugsmt.pcl$meta$cv <- bugsmt.pcl$meta[,covariates[i]]

        mttest <- na.omit(bugsmt.pcl$meta[,c("active", "cv", "subject")])

        timevarying <- !all(sapply(split(mttest$cv, mttest$subject), function(x)all(x==x[1])))

        if (time_delays[j] == 0 || timevarying) {
            mdl <- glmer(active ~ cv + (1 | subject), data=mttest, family='binomial', nAGQ = 0)
            if (F) {
                convergence_failed[j,i] <- T

                if (is.numeric(mttest$cv)) {
                    mdl <- glm(active ~ cv, data=mttest, family='binomial')

                    ms <- summary(mdl)
                    est <- ms$coefficients["cv", "Estimate"]
                    p <- ms$coefficients["cv", "Pr(>|z|)"]
                } else {
                    fe <- fisher.test(mttest$cv, mttest$active, hybrid=length(levels(mttest$cv)) > 3)
                    if (is.null(fe$estimate)) {
                        est <- NA
                    } else {
                        est <- log(fe$estimate)
                    }
                    p <- fe$p.value
                }
            } else if (is.factor(mttest$cv) && length(levels(mttest$cv)) > 2) {
                ms <- Anova(mdl)
                p <- ms$`Pr(>Chisq)`
                est <- NA
            } else {
                ms <- summary(mdl)
                est <- ms$coefficients[2, "Estimate"]
                p <- ms$coefficients[2, "Pr(>|z|)"]
            }

            assoc_beta[j,i] <- est
            assoc_p[j,i] <- p
        }
    }
}

assoc_fdr <- matrix(p.adjust(assoc_p, method="fdr"), nrow=nrow(assoc_p))
rownames(assoc_fdr) <- as.character(time_delays)
colnames(assoc_fdr) <- names(covariates)

library(pheatmap)
library(RColorBrewer)
pdf("./disease_activity/correlations_with_metadata.pdf", 5, 3)
pheatmap(assoc_p, cluster_rows=F, cluster_cols=F, display_numbers=T,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "Oranges")[1:6]))(100))
dev.off()

