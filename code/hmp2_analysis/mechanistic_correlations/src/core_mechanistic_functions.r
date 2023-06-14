perfeature_model<- function(pcl, name) {

    # Merge disease activity
    if (name != 'bugs') {
        pcl<-merge_disease_activity(pcl, lenience = 0)
    }

    # Add antibiotics
    pcl <- merge_metadata(pcl, fields=c("Antibiotics"))

    # Filter out features with no variance or with >90% zeros
    pcl <- pcl %>% pcl.filter.f(any(is.finite(x)) && var(x, na.rm=T) > 0 && mean(!is.na(x) & (x>0)) >= 0.1)

    # Filter out features with little variance (except for metabolites)
    if (name!='metabolites'){
        sds <- pcl.apply.f(pcl, sd(x, na.rm=T))
        pcl <- pcl %>% pcl.filter.f(keep=sds > median(sds) / 2)
    }

    # Extract relevant metadata
    df <- pcl$meta
    df$active<-as.factor(df$active*1) # Convert to binary
    df$active <- relevel(df$active, ref = "0") # Set 'inactive' as reference

    # Allocate output matrices
    resid_full <- matrix(NA, pcl$ns, pcl$nf)
    rownames(resid_full) <- rownames(pcl$x)
    colnames(resid_full) <- colnames(pcl$x)
    resid_empty <- matrix(NA, pcl$ns, pcl$nf)
    rownames(resid_empty) <- rownames(pcl$x)
    colnames(resid_empty) <- colnames(pcl$x)
    usubjects <- unique(pcl$meta$subject)
    subjI <- match(usubjects, pcl$meta$subject)
    subj_full <- matrix(NA, length(subjI), pcl$nf)
    rownames(subj_full) <- usubjects
    colnames(subj_full) <- colnames(pcl$x)
    subj_empty <- matrix(NA, length(subjI), pcl$nf)
    rownames(subj_empty) <- usubjects
    colnames(subj_empty) <- colnames(pcl$x)

    # Remove Empty Column
    df<-df[,colnames(df)!=""]

    # Fit full model

    for (i in seq_along(colnames(pcl$x))) {
        df$y <- pcl$x[rownames(df),i]
        mdl_full <- try(lmer(y ~ diagnosis + diagnosis%in%active + Antibiotics  + consent_age + (1 | site_name) + (1 | subject), data=df, na.action = na.omit))
        mdl_empty <- try(lmer(y ~ diagnosis + Antibiotics  + consent_age + (1 | site_name) + (1 | subject), data=df, na.action = na.omit))

        if (class(mdl_full)=='try-error') {
            resid_full[, i]<-NA
            subj_full[,i]<-NA

        }
        else{

            # Save full residuals for downstream analysis
            r <- residuals(mdl_full, type="pearson")
            resid_full[names(r),i] <- r


            # Save full subject-specific effects for downstream analysis
            cf <- coef(mdl_full)$subject
            subj_full[,i] <- cf[match(usubjects, rownames(cf)),1]
        }

        if (class(mdl_empty)=='try-error') {
            resid_empty[, i]<-NA
            subj_empty[,i]<-NA
        }
        else{

            # Save empty residuals for downstream analysis
            r <- residuals(mdl_empty, type="pearson")
            resid_empty[names(r),i] <- r


            # Save empty subject-specific effects for downstream analysis
            cf <- coef(mdl_empty)$subject
            subj_empty[,i] <- cf[match(usubjects, rownames(cf)),1]
        }
    }
    return (list(residuals_full=resid_full,
                 residuals_empty=resid_empty,
                 subjects_full=subj_full,
                 subjects_empty=subj_empty
    ))
}


perfeature_model_binomial<- function(pcl, name) {

    # Merge disease activity
    if (name != 'bugs') {
        pcl<-merge_disease_activity(pcl, lenience=0)
    }

    # Add antibiotics
    pcl <- merge_metadata(pcl, fields=c("Antibiotics"))

    # Extract relevant metadata
    df <- pcl$meta
    df$active<-as.factor(df$active*1) # Convert to Binary
    df$active <- relevel(df$active, ref = "0") # Set 'inactive' as reference

    # Allocate output matrices
    resid_full <- matrix(NA, pcl$ns, pcl$nf)
    rownames(resid_full) <- rownames(pcl$x)
    colnames(resid_full) <- colnames(pcl$x)
    resid_empty <- matrix(NA, pcl$ns, pcl$nf)
    rownames(resid_empty) <- rownames(pcl$x)
    colnames(resid_empty) <- colnames(pcl$x)
    usubjects <- unique(pcl$meta$subject)
    subjI <- match(usubjects, pcl$meta$subject)
    subj_full <- matrix(NA, length(subjI), pcl$nf)
    rownames(subj_full) <- usubjects
    colnames(subj_full) <- colnames(pcl$x)
    subj_empty <- matrix(NA, length(subjI), pcl$nf)
    rownames(subj_empty) <- usubjects
    colnames(subj_empty) <- colnames(pcl$x)

    # Remove Empty Column
    df<-df[,colnames(df)!=""]

    # Fit full model

    for (i in seq_along(colnames(pcl$x))) {
        df$y <- pcl$x[rownames(df),i]
        mdl_full <- try(glmer(y ~ diagnosis + diagnosis%in%active + Antibiotics  + consent_age + (1 | site_name) + (1 | subject), data=df, na.action = na.omit,  family=binomial))
        mdl_empty <- try(glmer(y ~ diagnosis + Antibiotics  + consent_age + (1 | site_name) + (1 | subject), data=df, na.action = na.omit,  family=binomial))

        if (class(mdl_full)=='try-error') {
            resid_full[, i]<-NA
            subj_full[,i]<-NA

        }
        else{

            # Save full residuals for downstream analysis
            r <- residuals(mdl_full, type="pearson")
            resid_full[names(r),i] <- r


            # Save full subject-specific effects for downstream analysis
            cf <- coef(mdl_full)$subject
            subj_full[,i] <- cf[match(usubjects, rownames(cf)),1]
        }

        if (class(mdl_empty)=='try-error') {
            resid_empty[, i]<-NA
            subj_empty[,i]<-NA
        }
        else{

            # Save empty residuals for downstream analysis
            r <- residuals(mdl_empty, type="pearson")
            resid_empty[names(r),i] <- r


            # Save empty residuals for downstream analysis
            cf <- coef(mdl_empty)$subject
            subj_empty[,i] <- cf[match(usubjects, rownames(cf)),1]
        }
    }
    return (list(residuals_full=resid_full,
                 residuals_empty=resid_empty,
                 subjects_full=subj_full,
                 subjects_empty=subj_empty
    ))
}


extract_residuals <- function(pcl, name) {
    outfile <- paste(HMP2_root, sprintf("/analysis/mechanistic_correlations/residuals/residuals_%s.RData", name), sep='')


    if (!file.exists(outfile)) {

        # Extract residuals and subject random effects (full + empty models)
        if (name != 'viruses') {
            res <- perfeature_model(pcl, name)
        } else {
            res <- perfeature_model_binomial(pcl, name)
        }

        # Filter out subjects with too few timepoints (< 4)
        samps_per_subj <- sapply(split(pcl$meta$subject, pcl$meta$subject), length)
        keep_subj <- names(samps_per_subj)[samps_per_subj >= 4]
        pcl$meta<- pcl$meta[pcl$meta$subject %in% keep_subj,]
        res$residuals_full <- res$residuals_full[rownames(pcl$meta),]
        res$residuals_empty <- res$residuals_empty[rownames(pcl$meta),]
        res$subjects_full <- res$subjects_full[keep_subj,,drop=F]
        res$subjects_empty <- res$subjects_empty[keep_subj,,drop=F]

        # Add needed metadata to the results structure
        res$rmeta <- pcl$meta
        res$name <- name

        # Dump data
        save(res, file=outfile)

    } else {
        # Just read from disk
        res = LoadToEnvironment(outfile)
    }

    return (res)
}


match_residuals <- function(res1, res2) {

    # Intialize PCLs
    pcl.full1<- pcl.full2 <- pcl.empty1 <-pcl.empty2 <-list()

    # Extract full residuals and make PCLs
    pcl.full1$x <- res1$res$residuals_full
    pcl.full1$meta<-res1$res$rmeta

    pcl.full2$x <- res2$res$residuals_full
    pcl.full2$meta<-res2$res$rmeta

    pcl.empty1$x <- res1$res$residuals_empty
    pcl.empty1$meta<-res1$res$rmeta

    pcl.empty2$x <- res2$res$residuals_empty
    pcl.empty2$meta<-res2$res$rmeta

    # Match PCLs
    matched_pcl.full<-match_datasets(list(pcl.full1, pcl.full2), lenience=4, matching=F)
    matched_pcl.empty<-match_datasets(list(pcl.empty1, pcl.empty2), lenience=4, matching=F)

    # Extract matched residuals
    rdat1.full<-matched_pcl.full[[1]]$x
    rdat2.full<-matched_pcl.full[[2]]$x
    rdat1.empty<-matched_pcl.empty[[1]]$x
    rdat2.empty<-matched_pcl.empty[[2]]$x

    # Extract subject random effects
    common_subj.full <- intersect(rownames(res1$res$subjects_full), rownames(res2$res$subjects_full))
    sdat1.full <- res1$res$subjects_full[rownames(res1$res$subjects_full) %in% common_subj.full,, drop=F]
    sdat2.full <- res2$res$subjects_full[rownames(res2$res$subjects_full) %in% common_subj.full,, drop=F]

    common_subj.empty <- intersect(rownames(res1$res$subjects_empty), rownames(res2$res$subjects_empty))
    sdat1.empty <- res1$res$subjects_empty[rownames(res1$res$subjects_empty) %in% common_subj.empty,, drop=F]
    sdat2.empty <- res2$res$subjects_empty[rownames(res2$res$subjects_empty) %in% common_subj.empty,, drop=F]

    # HAllA-friendly outputs (features in rows, samples in columns, same column names pairwise)
    rrdat1.full<-as.data.frame(t(rdat1.full));
    colnames(rrdat1.full)<-paste('Sample', 1: ncol(rrdat1.full), sep='')
    rrdat1.full<-rownames_to_column(rrdat1.full, '#')


    rrdat2.full<-as.data.frame(t(rdat2.full));
    colnames(rrdat2.full)<-paste('Sample', 1: ncol(rrdat2.full), sep='')
    rrdat2.full<-rownames_to_column(rrdat2.full, '#')

    ssdat1.full<-as.data.frame(t(sdat1.full));
    ssdat1.full<-rownames_to_column(ssdat1.full, '#')

    ssdat2.full<-as.data.frame(t(sdat2.full));
    ssdat2.full<-rownames_to_column(ssdat2.full, '#')

    rrdat1.empty<-as.data.frame(t(rdat1.empty));
    colnames(rrdat1.empty)<-paste('Sample', 1: ncol(rrdat1.empty), sep='')
    rrdat1.empty<-rownames_to_column(rrdat1.empty, '#')

    rrdat2.empty<-as.data.frame(t(rdat2.empty));
    colnames(rrdat2.empty)<-paste('Sample', 1: ncol(rrdat2.empty), sep='')
    rrdat2.empty<-rownames_to_column(rrdat2.empty, '#')

    ssdat1.empty<-as.data.frame(t(sdat1.empty));
    ssdat1.empty<-rownames_to_column(ssdat1.empty, '#')

    ssdat2.empty<-as.data.frame(t(sdat2.empty));
    ssdat2.empty<-rownames_to_column(ssdat2.empty, '#')

    # Dump out the data tables for HAllA and AllA
    cmp_name <- sprintf("%s-%s", res1$res$name, res2$res$name)
    write.table(rrdat1.full,
                file=paste(HMP2_root,sprintf("/analysis/mechanistic_correlations/halla_input/%s_%s_residuals_full.tsv",
                                             cmp_name, res1$res$nam), sep=''),  sep = "\t", eol = "\n", quote = F, col.names = T, row.names = F)
    write.table(rrdat2.full,
                file=paste(HMP2_root,sprintf("/analysis/mechanistic_correlations/halla_input/%s_%s_residuals_full.tsv",
                                             cmp_name, res2$res$nam), sep=''),  sep = "\t", eol = "\n", quote = F, col.names = T, row.names = F)
    write.table(ssdat1.full,
                file=paste(HMP2_root,sprintf("/analysis/mechanistic_correlations/halla_input/%s_%s_subjects_full.tsv",
                                             cmp_name, res1$res$nam), sep=''),  sep = "\t", eol = "\n", quote = F, col.names = T, row.names = F)
    write.table(ssdat2.full,
                file=paste(HMP2_root,sprintf("/analysis/mechanistic_correlations/halla_input/%s_%s_subjects_full.tsv",
                                             cmp_name, res2$res$nam), sep=''),  sep = "\t", eol = "\n", quote = F, col.names = T, row.names = F)
    write.table(rrdat1.empty,
                file=paste(HMP2_root,sprintf("/analysis/mechanistic_correlations/halla_input/%s_%s_residuals_empty.tsv",
                                             cmp_name, res1$res$nam), sep=''),  sep = "\t", eol = "\n", quote = F, col.names =T, row.names = F)
    write.table(rrdat2.empty,
                file=paste(HMP2_root,sprintf("/analysis/mechanistic_correlations/halla_input/%s_%s_residuals_empty.tsv",
                                             cmp_name, res2$res$nam), sep=''),  sep = "\t", eol = "\n", quote = F, col.names = T, row.names = F)
    write.table(ssdat1.empty,
                file=paste(HMP2_root,sprintf("/analysis/mechanistic_correlations/halla_input/%s_%s_subjects_empty.tsv",
                                             cmp_name, res1$res$nam), sep=''),  sep = "\t", eol = "\n", quote = F, col.names = T, row.names = F)
    write.table(ssdat2.empty,
                file=paste(HMP2_root,sprintf("/analysis/mechanistic_correlations/halla_input/%s_%s_subjects_empty.tsv",
                                             cmp_name, res2$res$nam), sep=''),  sep = "\t", eol = "\n", quote = F, col.names = T, row.names = F)
}


# Load environment
LoadToEnvironment <- function(RData, env = new.env()){
    load(RData, env)
    return(env)
}

extract_residuals_fecalcal <- function() {

    outfile <- paste(HMP2_root, "/analysis/mechanistic_correlations/residuals/residuals_fecalcal.RData", sep="")


    if (!file.exists(outfile)) {


        pcl<-species.pcl.vst

        # Load general metadata table
        metadata<-hmp2_sample_metadata
        metadata<-metadata[metadata$data_type %in% c('metagenomics', 'noproduct'),]

        # Add antibiotics
        pcl <- merge_metadata(pcl, fields=c("Antibiotics"))

        # Extract relevant metadata
        df <- pcl$meta
        df$active<-as.factor(df$active*1) # Convert to Binary
        df$active <- relevel(df$active, ref = "0") # Set 'inactive' as reference

        # VST for fecalcal
        df$fecalcal<-log2(df$fecalcal+1)

        # Remove Empty Column
        df<-df[,colnames(df)!=""]

        # Fit full model
        mdl_full <- lmer(fecalcal ~ diagnosis + diagnosis%in%active + Antibiotics  + consent_age + (1 | site_name) + (1 | subject), data=df, na.action = na.omit)
        mdl_empty <- lmer(fecalcal ~ diagnosis + Antibiotics  + consent_age + (1 | site_name) + (1 | subject), data=df, na.action = na.omit)

        # Save full residuals for downstream analysis
        r_full <- as.data.frame(residuals(mdl_full, type="pearson"))
        names(r_full)<-'fecalcal'

        # Save full subject-specific effects for downstream analysis
        cf_full <- as.data.frame(coef(mdl_full)$subject)
        subj_full <- as.data.frame(cf_full[,"(Intercept)"])
        rownames(subj_full)<- rownames(cf_full)
        colnames(subj_full)<-'fecalcal'


        # Save empty residuals for downstream analysis
        r_empty <- as.data.frame(residuals(mdl_empty, type="pearson"))
        names(r_empty)<-'fecalcal'

        # Save empty subject-specific effects for downstream analysis
        cf_empty <- as.data.frame(coef(mdl_empty)$subject)
        subj_empty <- as.data.frame(cf_empty[,"(Intercept)"])
        rownames(subj_empty)<- rownames(cf_empty)
        colnames(subj_empty)<-'fecalcal'

        # Filter out subjects with too few timepoints (< 4)
        samps_per_subj <- sapply(split(pcl$meta$subject, pcl$meta$subject), length)
        keep_subj <- names(samps_per_subj)[samps_per_subj >= 4]
        pcl$meta<- pcl$meta[pcl$meta$subject %in% keep_subj,]
        r_full <- r_full[intersect(rownames(pcl$meta),rownames(r_full)), , drop=F]
        r_empty <- r_empty[intersect(rownames(pcl$meta),rownames(r_empty)), , drop=F]
        subj_full <- subj_full[intersect(keep_subj,rownames(subj_full)),,drop=F]
        subj_empty <- subj_empty[intersect(keep_subj,rownames(subj_empty)),,drop=F]

        # Add needed metadata to the results structure
        res<-list()
        res$residuals_full <- r_full
        res$residuals_empty <- r_empty
        res$subjects_full <- subj_full
        res$subjects_empty <- subj_empty
        res$rmeta <- pcl$meta[rownames(res$residuals_full),]
        res$name <- 'fecalcal'

        # Dump data
        save(res, file=outfile)

    } else {
        # Just read from disk
        res = LoadToEnvironment(outfile)
    }

    return (res)
}


perfeature_model_baseline<- function(pcl, name) {

    # Merge disease activity
    if (!name %in% c('bugs','serology', 'HTXrectum', 'HTXileum')) {
        pcl<-merge_disease_activity(pcl, lenience=0)
    }

    if (name %in% c('serology', 'HTXrectum', 'HTXileum')) {
        pcl<-merge_disease_activity(pcl, lenience=4)
    }

    # Add antibiotics
    pcl <- merge_metadata(pcl, fields=c("Antibiotics"))

    # Filter out features no variance or with >90% zeros
    pcl <- pcl %>% pcl.filter.f(any(is.finite(x)) && var(x, na.rm=T) > 0 && mean(!is.na(x) & (x>0)) >= 0.1)

    # Filter out features with little variance (not for metabolites)
    if (name!='metabolites'){
        sds <- pcl.apply.f(pcl, sd(x, na.rm=T))
        pcl <- pcl %>% pcl.filter.f(keep=sds > median(sds) / 2)
    }

    # Extract relevant metadata
    df0 <- pcl$meta
    df0$active<-as.factor(df0$active*1) # Convert to binary
    df0$active <- relevel(df0$active, ref = "0") # Set 'inactive' as reference

    # Extract baseline subjects
    df0<-df0[!is.na(df0$collection),]
    df0<-df0[order(df0$collection),]
    df1<-df0[!duplicated(df0$subject),]
    df<-df1[df1$collection<3,]


    # Allocate output matrices
    resid_full <- matrix(NA, nrow(df), pcl$nf)
    rownames(resid_full) <- rownames(df)
    colnames(resid_full) <- colnames(pcl$x)
    resid_empty <- matrix(NA, nrow(df), pcl$nf)
    rownames(resid_empty) <- rownames(df)
    colnames(resid_empty) <- colnames(pcl$x)

    # Remove Empty Column
    df<-df[,colnames(df)!=""]


    # Fit full model
    for (i in seq_along(colnames(pcl$x))) {
        df$y <- pcl$x[rownames(df),i]

        if (name %in% c('serology', 'HTXrectum', 'HTXileum')) {

            mdl_full <- try(lmer(y ~ diagnosis + diagnosis%in%active + consent_age + (1|site_name), data=df, na.action = na.omit))
            mdl_empty<- try(lmer(y ~ diagnosis + consent_age + (1|site_name), data=df, na.action = na.omit))
        }
        else{

            mdl_full <- try(lmer(y ~ diagnosis + diagnosis%in%active + Antibiotics + consent_age + (1|site_name), data=df, na.action = na.omit))
            mdl_empty<- try(lmer(y ~ diagnosis + Antibiotics + consent_age + (1|site_name), data=df, na.action = na.omit))
        }

        if (class(mdl_full)=='try-error') {
            resid_full[, i]<-NA

        }
        else{

            # Save full residuals for downstream analysis
            r <- residuals(mdl_full, type="pearson")
            resid_full[names(r),i] <- r

        }

        if (class(mdl_empty)=='try-error') {
            resid_empty[, i]<-NA
        }
        else{

            # Save empty residuals for downstream analysis
            r <- residuals(mdl_empty, type="pearson")
            resid_empty[names(r),i] <- r

        }
    }
    return (list(residuals_full=resid_full,
                 residuals_empty=resid_empty))
}

extract_residuals_baseline <- function(pcl, name) {
    outfile <- paste(HMP2_root, sprintf("/analysis/mechanistic_correlations/residuals/baseline_residuals_%s.RData", name), sep='')


    if (!file.exists(outfile)) {

        # Extract residuals (full + empty models)
        res <- perfeature_model_baseline(pcl, name)

        # Match rownames
        pcl$meta<-pcl$meta[rownames(res$residuals_full),]

        # Add needed metadata to the results structure
        res$rmeta <- pcl$meta
        res$name <- name

        # Dump data
        save(res, file=outfile)

    } else {
        # Just read from disk
        res = LoadToEnvironment(outfile)
    }

    return (res)
}



match_residuals_baseline <- function(res1, res2) {

    # Intialize PCLs
    pcl.full1<- pcl.full2 <- pcl.empty1 <-pcl.empty2 <-list()

    # Extract full residuals and make PCLs
    pcl.full1$x <- res1$res$residuals_full
    pcl.full1$meta<-res1$res$rmeta

    pcl.full2$x <- res2$res$residuals_full
    pcl.full2$meta<-res2$res$rmeta

    pcl.empty1$x <- res1$res$residuals_empty
    pcl.empty1$meta<-res1$res$rmeta

    pcl.empty2$x <- res2$res$residuals_empty
    pcl.empty2$meta<-res2$res$rmeta

    # Match PCLs
    matched_pcl.full<-match_datasets(list(pcl.full1, pcl.full2), lenience=4, matching=F)
    matched_pcl.empty<-match_datasets(list(pcl.empty1, pcl.empty2), lenience=4, matching=F)

    # Extract matched residuals
    rdat1.full<-matched_pcl.full[[1]]$x
    rdat2.full<-matched_pcl.full[[2]]$x
    rdat1.empty<-matched_pcl.empty[[1]]$x
    rdat2.empty<-matched_pcl.empty[[2]]$x

    # HAllA-friendly outputs (features in rows, samples in columns, same column names pairwise)
    rrdat1.full<-as.data.frame(t(rdat1.full));
    colnames(rrdat1.full)<-paste('Sample', 1: ncol(rrdat1.full), sep='')
    rrdat1.full<-rownames_to_column(rrdat1.full, '#')


    rrdat2.full<-as.data.frame(t(rdat2.full));
    colnames(rrdat2.full)<-paste('Sample', 1: ncol(rrdat2.full), sep='')
    rrdat2.full<-rownames_to_column(rrdat2.full, '#')

    rrdat1.empty<-as.data.frame(t(rdat1.empty));
    colnames(rrdat1.empty)<-paste('Sample', 1: ncol(rrdat1.empty), sep='')
    rrdat1.empty<-rownames_to_column(rrdat1.empty, '#')

    rrdat2.empty<-as.data.frame(t(rdat2.empty));
    colnames(rrdat2.empty)<-paste('Sample', 1: ncol(rrdat2.empty), sep='')
    rrdat2.empty<-rownames_to_column(rrdat2.empty, '#')


    # Dump out the data tables for HAllA and AllA
    cmp_name <- sprintf("%s-%s", res1$res$name, res2$res$name)
    write.table(rrdat1.full,
                file=paste(HMP2_root,sprintf("/analysis/mechanistic_correlations/halla_input/%s_%s_residuals_full.tsv",
                                             cmp_name, res1$res$nam), sep=''),  sep = "\t", eol = "\n", quote = F, col.names = T, row.names = F)
    write.table(rrdat2.full,
                file=paste(HMP2_root,sprintf("/analysis/mechanistic_correlations/halla_input/%s_%s_residuals_full.tsv",
                                             cmp_name, res2$res$nam), sep=''),  sep = "\t", eol = "\n", quote = F, col.names = T, row.names = F)
    write.table(rrdat1.empty,
                file=paste(HMP2_root,sprintf("/analysis/mechanistic_correlations/halla_input/%s_%s_residuals_empty.tsv",
                                             cmp_name, res1$res$nam), sep=''),  sep = "\t", eol = "\n", quote = F, col.names =T, row.names = F)
    write.table(rrdat2.empty,
                file=paste(HMP2_root,sprintf("/analysis/mechanistic_correlations/halla_input/%s_%s_residuals_empty.tsv",
                                             cmp_name, res2$res$nam), sep=''),  sep = "\t", eol = "\n", quote = F, col.names = T, row.names = F)
}
