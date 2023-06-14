run_DA<-function(pcl, name, transformation) {

    cat(sprintf('Running a combination of %s and %s...\n', name, transformation))

    # Merge disease activity
    if (name != 'bugs') {
        pcl<-merge_disease_activity(pcl, lenience = 0)
    }

    # Add antibiotics
    pcl <- merge_metadata(pcl, fields=c("Antibiotics"))

    # Filter out features with no variance or with >90% zeros
    pcl <- pcl %>% pcl.filter.f(any(is.finite(x)) && var(x, na.rm=T) > 0 && mean(!is.na(x) & (x>0)) >= 0.1)

    # Filter out features with little variance (except for metabolites)
    if (!name %in% c('metabolites', 'metabolites_Sena','metabolites_all','metabolites_nadine')){
        sds <- pcl.apply.f(pcl, sd(x, na.rm=T))
        pcl <- pcl %>% pcl.filter.f(keep=sds > median(sds) / 2)
    }

    # Save prevalence
    diag_sig <- data.frame(prevalence = pcl %>% pcl.apply.f(mean(x>0, na.rm=T)))

    # Transformation
    pcl<-pcl.transform(pcl, name, transformation)

    # Initialize
    diag_sig$coefCD <- NA
    diag_sig$coefUC <- NA
    diag_sig$coefActiveCD <- NA
    diag_sig$coefActiveUC <- NA
    diag_sig$coefActivenonIBD <- NA
    diag_sig$tvalCD <- NA
    diag_sig$tvalUC <- NA
    diag_sig$tvalActiveCD <- NA
    diag_sig$tvalActiveUC <- NA
    diag_sig$tvalActivenonIBD <- NA
    diag_sig$pvalCD <- NA
    diag_sig$pvalUC <- NA
    diag_sig$pvalActiveCD <- NA
    diag_sig$pvalActiveUC <- NA
    diag_sig$pvalActivenonIBD <- NA

    # Extract relevant metadata
    df <- pcl$meta
    df$active<-as.factor(df$active*1) # Convert to binary
    df$active <- relevel(df$active, ref = "0") # Set 'inactive' as reference

    # Remove Empty Column
    df<-df[,colnames(df)!=""]

    # Formula for fixed effects
    formula = as.formula(y ~ diagnosis + diagnosis%in%active + Antibiotics  + consent_age)

    # Fit full model
    for (i in seq_along(colnames(pcl$x))) {
        df$y <- pcl$x[rownames(df),i]
        mdl <- try(lme(formula , random= ~ 1 | site_name/subject, data=df, na.action = na.omit))

        if (class(mdl)=='try-error') {
            diag_sig$coefCD[i] <- NA
            diag_sig$coefUC[i] <- NA
            diag_sig$coefActiveCD[i] <- NA
            diag_sig$coefActiveUC[i] <- NA
            diag_sig$coefActivenonIBD[i] <- NA
            diag_sig$tvalCD[i] <- NA
            diag_sig$tvalUC[i] <- NA
            diag_sig$tvalActiveCD[i] <- NA
            diag_sig$tvalActiveUC[i] <- NA
            diag_sig$tvalActivenonIBD[i] <- NA
            diag_sig$pvalCD[i] <- NA
            diag_sig$pvalUC[i] <- NA
            diag_sig$pvalActiveCD[i] <- NA
            diag_sig$pvalActiveUC[i] <- NA
            diag_sig$pvalActivenonIBD[i] <- NA
        }
        else{

            # Extract relevant coefficients
            diag_sig$coefCD[i]  <- summary(mdl)$tTable['diagnosisCD', 'Value']
            diag_sig$coefUC[i] <- summary(mdl)$tTable['diagnosisUC', 'Value']
            diag_sig$coefActiveCD[i]<-summary(mdl)$tTable['diagnosisCD:active1', 'Value']
            diag_sig$coefActiveUC[i]<-summary(mdl)$tTable['diagnosisUC:active1', 'Value']
            diag_sig$coefActivenonIBD[i]<-summary(mdl)$tTable['diagnosisnonIBD:active1', 'Value']
            diag_sig$tvalCD[i] <-  summary(mdl)$tTable['diagnosisCD', 't-value']
            diag_sig$tvalUC[i] <-summary(mdl)$tTable['diagnosisUC', 't-value']
            diag_sig$tvalActiveCD[i] <- summary(mdl)$tTable['diagnosisCD:active1', 't-value']
            diag_sig$tvalActiveUC[i] <- summary(mdl)$tTable['diagnosisUC:active1', 't-value']
            diag_sig$tvalActivenonIBD[i] <- summary(mdl)$tTable['diagnosisnonIBD:active1', 't-value']
            diag_sig$pvalCD[i] <- summary(mdl)$tTable['diagnosisCD', 'p-value']
            diag_sig$pvalUC[i] <- summary(mdl)$tTable['diagnosisUC', 'p-value']
            diag_sig$pvalActiveCD[i] <- summary(mdl)$tTable['diagnosisCD:active1', 'p-value']
            diag_sig$pvalActiveUC[i] <- summary(mdl)$tTable['diagnosisUC:active1', 'p-value']
            diag_sig$pvalActivenonIBD[i] <- summary(mdl)$tTable['diagnosisnonIBD:active1', 'p-value']
        }
    }

    # Multiplicity correction
    diag_sig$qvalCD <- p.adjust(diag_sig$pvalCD, 'fdr')
    diag_sig$qvalUC <- p.adjust(diag_sig$pvalUC, 'fdr')
    diag_sig$qvalActiveCD <- p.adjust(diag_sig$pvalActiveCD, 'fdr')
    diag_sig$qvalActiveUC <- p.adjust(diag_sig$pvalActiveUC, 'fdr')
    diag_sig$qvalActivenonIBD <- p.adjust(diag_sig$pvalActivenonIBD, 'fdr')


    # Separate diagnosis and activity
    diagnosis_only<-c('coefCD', 'coefUC', 'tvalCD', 'tvalUC', 'pvalCD', 'pvalUC', 'qvalCD', 'qvalUC')
    diag_sig_diagnosis<-diag_sig[, c('prevalence', diagnosis_only)]
    diag_sig_activity<-diag_sig[, colnames(diag_sig)[!colnames(diag_sig) %in% diagnosis_only]]

    # Order by minimum q-value
    diag_sig_diagnosis$minQ<-pmin(diag_sig_diagnosis$qvalCD, diag_sig_diagnosis$qvalUC)
    diag_sig_diagnosis <- diag_sig_diagnosis[order(diag_sig_diagnosis$minQ),]

    diag_sig_activity$minQ<-pmin(diag_sig_activity$qvalActiveCD, diag_sig_activity$qvalActiveUC, diag_sig_activity$qvalActivenonIBD)
    diag_sig_activity <- diag_sig_activity[order(diag_sig_activity$minQ),]

    # Remove row names
    diag_sig_diagnosis<-rownames_to_column(diag_sig_diagnosis, '#')
    diag_sig_activity<-rownames_to_column(diag_sig_activity, '#')

    # Remove NA's
    library(tidyr)
    diag_sig_diagnosis<- diag_sig_diagnosis %>% drop_na()
    diag_sig_activity<- diag_sig_activity %>% drop_na()


    # Dump tables
    write.table(diag_sig_activity, file=paste(HMP2_root,
                                              sprintf("/analysis/differential_abundance/%s/active_vs_inactive/Active_IBD_relations_%s.tsv",  transformation, name), sep=''),
                row.names=F,col.names=T,quote=F,sep="\t")
    write.table(diag_sig_diagnosis, file=paste(HMP2_root,
                                               sprintf("/analysis/differential_abundance/%s/diagnosis/IBD_relations_%s.tsv",  transformation, name), sep=''),
                row.names=F,col.names=T,quote=F,sep="\t")
}


run_DA_binomial<-function(pcl, name, transformation) {

    cat(sprintf('Running a combination of %s and %s...\n', name, transformation))

    # Merge disease activity
    pcl<-merge_disease_activity(pcl, lenience = 0)

    # Add antibiotics
    pcl <- merge_metadata(pcl, fields=c("Antibiotics"))

    # Save prevalence
    diag_sig <- data.frame(prevalence = pcl %>% pcl.apply.f(mean(x>0, na.rm=T)))
    uniquenames<-make.names(colnames(pcl$x), unique=TRUE)
    rownames(diag_sig) <- uniquenames

    # Initialize
    diag_sig$coefCD <- NA
    diag_sig$coefUC <- NA
    diag_sig$coefActiveCD <- NA
    diag_sig$coefActiveUC <- NA
    diag_sig$coefActivenonIBD <- NA
    diag_sig$zvalCD <- NA
    diag_sig$zvalUC <- NA
    diag_sig$zvalActiveCD <- NA
    diag_sig$zvalActiveUC <- NA
    diag_sig$zvalActivenonIBD <- NA
    diag_sig$pvalCD <- NA
    diag_sig$pvalUC <- NA
    diag_sig$pvalActiveCD <- NA
    diag_sig$pvalActiveUC <- NA
    diag_sig$pvalActivenonIBD <- NA

    # Extract relevant metadata
    df <- pcl$meta
    df$active<-as.factor(df$active*1) # Convert to binary
    df$active <- relevel(df$active, ref = "0") # Set 'inactive' as reference

    # Remove Empty Column
    df<-df[,colnames(df)!=""]

    # Fit full model

    for (i in seq_along(colnames(pcl$x))) {
        df$y <- pcl$x[rownames(df),i]
        mdl <- try(glmer(y ~ diagnosis + diagnosis%in%active + Antibiotics + consent_age + ( 1 | site_name/subject), data=df, family=binomial))

        if (class(mdl)=='try-error') {
            diag_sig$coefCD[i] <- NA
            diag_sig$coefUC[i] <- NA
            diag_sig$coefActiveCD[i] <- NA
            diag_sig$coefActiveUC[i] <- NA
            diag_sig$coefActivenonIBD[i] <- NA
            diag_sig$zvalCD[i] <- NA
            diag_sig$zvalUC[i] <- NA
            diag_sig$zvalActiveCD[i] <- NA
            diag_sig$zvalActiveUC[i] <- NA
            diag_sig$zvalActivenonIBD[i] <- NA
            diag_sig$pvalCD[i] <- NA
            diag_sig$pvalUC[i] <- NA
            diag_sig$pvalActiveCD[i] <- NA
            diag_sig$pvalActiveUC[i] <- NA
            diag_sig$pvalActivenonIBD[i] <- NA
        }
        else{

            # Extract relevant coefficients
            diag_sig$coefCD[i]  <- summary(mdl)$coefficients['diagnosisCD', 1]
            diag_sig$coefUC[i] <- summary(mdl)$coefficients['diagnosisUC', 1]
            diag_sig$coefActiveCD[i] <- summary(mdl)$coefficients['diagnosisCD:active1', 1]
            diag_sig$coefActiveUC[i] <- summary(mdl)$coefficients['diagnosisUC:active1', 1]
            diag_sig$coefActivenonIBD[i] <-summary(mdl)$coefficients['diagnosisnonIBD:active1', 1]
            diag_sig$zvalCD[i] <-  summary(mdl)$coefficients['diagnosisCD', 3]
            diag_sig$zvalUC[i] <-summary(mdl)$coefficients['diagnosisUC', 3]
            diag_sig$zvalActiveCD[i] <- summary(mdl)$coefficients['diagnosisCD:active1', 3]
            diag_sig$zvalActiveUC[i] <- summary(mdl)$coefficients['diagnosisUC:active1', 3]
            diag_sig$zvalActivenonIBD[i] <- summary(mdl)$coefficients['diagnosisnonIBD:active1', 3]
            diag_sig$pvalCD[i] <- summary(mdl)$coefficients['diagnosisCD', 4]
            diag_sig$pvalUC[i] <- summary(mdl)$coefficients['diagnosisUC', 4]
            diag_sig$pvalActiveCD[i] <- summary(mdl)$coefficients['diagnosisCD:active1', 4]
            diag_sig$pvalActiveUC[i] <- summary(mdl)$coefficients['diagnosisUC:active1', 4]
            diag_sig$pvalActivenonIBD[i] <- summary(mdl)$coefficients['diagnosisnonIBD:active1', 4]
        }
    }

    # Multiplicity correction
    diag_sig$qvalCD <- p.adjust(diag_sig$pvalCD, 'fdr')
    diag_sig$qvalUC <- p.adjust(diag_sig$pvalUC, 'fdr')
    diag_sig$qvalActiveCD <- p.adjust(diag_sig$pvalActiveCD, 'fdr')
    diag_sig$qvalActiveUC <- p.adjust(diag_sig$pvalActiveUC, 'fdr')
    diag_sig$qvalActivenonIBD <- p.adjust(diag_sig$pvalActivenonIBD, 'fdr')


    # Separate diagnosis and activity
    diagnosis_only<-c( 'coefCD', 'coefUC', 'zvalCD', 'zvalUC', 'pvalCD', 'pvalUC', 'qvalCD', 'qvalUC')
    diag_sig_diagnosis<-diag_sig[, c('prevalence', diagnosis_only)]
    diag_sig_activity<-diag_sig[, colnames(diag_sig)[!colnames(diag_sig) %in% diagnosis_only]]

    # Order by minimum q-value
    diag_sig_diagnosis$minQ<-pmin(diag_sig_diagnosis$qvalCD, diag_sig_diagnosis$qvalUC)
    diag_sig_diagnosis <- diag_sig_diagnosis[order(diag_sig_diagnosis$minQ),]

    diag_sig_activity$minQ<-pmin(diag_sig_activity$qvalActiveCD, diag_sig_activity$qvalActiveUC, diag_sig_activity$qvalActivenonIBD)
    diag_sig_activity <- diag_sig_activity[order(diag_sig_activity$minQ),]

    # Remove row names
    diag_sig_diagnosis<-rownames_to_column(diag_sig_diagnosis, '#')
    diag_sig_activity<-rownames_to_column(diag_sig_activity, '#')

    # Remove NA's
    library(tidyr)
    diag_sig_diagnosis<- diag_sig_diagnosis %>% drop_na()
    diag_sig_activity<- diag_sig_activity %>% drop_na()

    # Dump tables
    write.table(diag_sig_activity, file=paste(HMP2_root,
                                              sprintf("/analysis/differential_abundance/%s/active_vs_inactive/Active_IBD_relations_%s.tsv",  transformation, name), sep=''),
                row.names=F,col.names=T,quote=F,sep="\t")
    write.table(diag_sig_diagnosis, file=paste(HMP2_root,
                                               sprintf("/analysis/differential_abundance/%s/diagnosis/IBD_relations_%s.tsv",  transformation, name), sep=''),
                row.names=F,col.names=T,quote=F,sep="\t")

}


run_DA_activity_index<-function(pcl, name, transformation) {

    cat(sprintf('Running a combination of %s and %s...\n', name, transformation))


    # Merge metadata
    if (name != 'bugs') {
        pcl<-merge_disease_activity(pcl, lenience = 0)
    }

    # Add antibiotics
    pcl <- merge_metadata(pcl, fields=c("Antibiotics"))

    # Filter out features with no variance or with >90% zeros
    pcl <- pcl %>% pcl.filter.f(any(is.finite(x)) && var(x, na.rm=T) > 0 && mean(!is.na(x) & (x>0)) >= 0.1)

    # Filter out features with little variance (except for metabolites)
    if (!name %in% c('metabolites', 'metabolites_Sena','metabolites_all','metabolites_nadine')){
        sds <- pcl.apply.f(pcl, sd(x, na.rm=T))
        pcl <- pcl %>% pcl.filter.f(keep=sds > median(sds) / 2)
    }

    # Save prevalence
    diag_sig <- data.frame(prevalence = pcl %>% pcl.apply.f(mean(x>0, na.rm=T)))

    # Transformation
    pcl<-pcl.transform(pcl, name, transformation)

    # Initialize
    diag_sig$coefCD <- NA
    diag_sig$coefUC <- NA
    diag_sig$coefActivityIndex <- NA
    diag_sig$tvalCD <- NA
    diag_sig$tvalUC <- NA
    diag_sig$tvalActivityIndex <- NA
    diag_sig$pvalCD <- NA
    diag_sig$pvalUC <- NA
    diag_sig$pvalActivityIndex <- NA

    # Extract metadata
    df <- pcl$meta

    # Remove Empty Column
    df<-df[,colnames(df)!=""]

    # Formula for fixed effects
    formula = as.formula(y ~ diagnosis + activity_index + Antibiotics  + consent_age)

    # Fit full model
    for (i in seq_along(colnames(pcl$x))) {
        df$y <- pcl$x[rownames(df),i]
        mdl <- try(lme(formula, random= ~ 1 | site_name/subject, data=df, na.action = na.omit))

        if (class(mdl)=='try-error') {
            diag_sig$coefCD[i] <- NA
            diag_sig$coefUC[i] <- NA
            diag_sig$coefActivityIndex[i] <- NA
            diag_sig$tvalCD[i] <- NA
            diag_sig$tvalUC[i] <- NA
            diag_sig$tvalActivityIndex[i] <- NA
            diag_sig$pvalCD[i] <- NA
            diag_sig$pvalUC[i] <- NA
            diag_sig$pvalActivityIndex[i] <- NA
        }
        else{

            # Extract relevant coefficients
            diag_sig$coefCD[i]  <- summary(mdl)$tTable['diagnosisCD', 'Value']
            diag_sig$coefUC[i] <- summary(mdl)$tTable['diagnosisUC', 'Value']
            diag_sig$coefActivityIndex[i] <-summary(mdl)$tTable['activity_index', 'Value']
            diag_sig$tvalCD[i] <-  summary(mdl)$tTable['diagnosisCD', 't-value']
            diag_sig$tvalUC[i] <-summary(mdl)$tTable['diagnosisUC', 't-value']
            diag_sig$tvalActivityIndex[i] <- summary(mdl)$tTable['activity_index', 't-value']
            diag_sig$pvalCD[i] <- summary(mdl)$tTable['diagnosisCD', 'p-value']
            diag_sig$pvalUC[i] <- summary(mdl)$tTable['diagnosisUC', 'p-value']
            diag_sig$pvalActivityIndex[i] <- summary(mdl)$tTable['activity_index', 'p-value']
        }
    }

    # Multiplicity correction
    diag_sig$qvalCD <- p.adjust(diag_sig$pvalCD, 'fdr')
    diag_sig$qvalUC <- p.adjust(diag_sig$pvalUC, 'fdr')
    diag_sig$qvalActivityIndex <- p.adjust(diag_sig$pvalActivityIndex, 'fdr')

    # Separate diagnosis and activity
    diagnosis_only<-c( 'coefCD', 'coefUC', 'tvalCD', 'tvalUC', 'pvalCD', 'pvalUC', 'qvalCD', 'qvalUC')
    diag_sig_diagnosis<-diag_sig[, c('prevalence', diagnosis_only)]
    diag_sig_activity<-diag_sig[, colnames(diag_sig)[!colnames(diag_sig) %in% diagnosis_only]]
    diag_sig_activity <- diag_sig_activity[order(diag_sig_activity$qvalActivityIndex),]

    # Order by minimum q-value
    diag_sig_diagnosis$minQ<-pmin(diag_sig_diagnosis$qvalCD, diag_sig_diagnosis$qvalUC)
    diag_sig_diagnosis <- diag_sig_diagnosis[order(diag_sig_diagnosis$minQ),]

    # Remove row names
    diag_sig_diagnosis<-rownames_to_column(diag_sig_diagnosis, '#')
    diag_sig_activity<-rownames_to_column(diag_sig_activity, '#')

    # Remove NA's
    library(tidyr)
    diag_sig_diagnosis<- diag_sig_diagnosis %>% drop_na()
    diag_sig_activity<- diag_sig_activity %>% drop_na()

    # Dump tables
    write.table(diag_sig_activity, file=paste(HMP2_root,
                                              sprintf("/analysis/differential_abundance/%s/active_vs_inactive/activity_index/Activity_Index_IBD_relations_%s.tsv",  transformation, name), sep=''),
                row.names=F,col.names=T,quote=F,sep="\t")
    write.table(diag_sig_diagnosis, file=paste(HMP2_root,
                                               sprintf("/analysis/differential_abundance/%s/diagnosis/activity_index/IBD_relations2_%s.tsv",  transformation, name), sep=''),
                row.names=F,col.names=T,quote=F,sep="\t")

}


run_DA_activity_index_binomial<-function(pcl, name, transformation) {

    cat(sprintf('Running a combination of %s and %s...\n', name, transformation))


    # Merge disease activity
    pcl<-merge_disease_activity(pcl, lenience = 0)

    # Add antibiotics
    pcl <- merge_metadata(pcl, fields=c("Antibiotics"))

    # Save prevalence
    diag_sig <- data.frame(prevalence = pcl %>% pcl.apply.f(mean(x>0, na.rm=T)))
    uniquenames<-make.names(colnames(pcl$x), unique=TRUE)
    rownames(diag_sig) <- uniquenames

    # Initialize
    diag_sig$coefCD <- NA
    diag_sig$coefUC <- NA
    diag_sig$coefActivityIndex <- NA
    diag_sig$zvalCD <- NA
    diag_sig$zvalUC <- NA
    diag_sig$zvalActivityIndex <- NA
    diag_sig$pvalCD <- NA
    diag_sig$pvalUC <- NA
    diag_sig$pvalActivityIndex <- NA

    # Extract relevant metadata
    df <- pcl$meta

    # Remove Empty Column
    df<-df[,colnames(df)!=""]

    # Fit full model
    for (i in seq_along(colnames(pcl$x))) {
        df$y <- pcl$x[rownames(df),i]
        mdl <- try(glmer(y ~ diagnosis + activity_index + Antibiotics + consent_age + ( 1 | site_name/subject), data=df, family=binomial))

        if (class(mdl)=='try-error') {
            diag_sig$coefCD[i] <- NA
            diag_sig$coefUC[i] <- NA
            diag_sig$coefActivityIndex[i] <- NA
            diag_sig$zvalCD[i] <- NA
            diag_sig$zvalUC[i] <- NA
            diag_sig$zvalActivityIndex[i] <- NA
            diag_sig$pvalCD[i] <- NA
            diag_sig$pvalUC[i] <- NA
            diag_sig$pvalActivityIndex[i] <- NA
        }
        else{

            # Extract relevant coefficients
            diag_sig$coefCD[i]  <- summary(mdl)$coefficients['diagnosisCD', 1]
            diag_sig$coefUC[i] <- summary(mdl)$coefficients['diagnosisUC', 1]
            diag_sig$coefActivityIndex[i] <-summary(mdl)$coefficients['activity_index', 1]
            diag_sig$zvalCD[i] <-  summary(mdl)$coefficients['diagnosisCD', 3]
            diag_sig$zvalUC[i] <-summary(mdl)$coefficients['diagnosisUC', 3]
            diag_sig$zvalActivityIndex[i] <- summary(mdl)$coefficients['activity_index', 3]
            diag_sig$pvalCD[i] <- summary(mdl)$coefficients['diagnosisCD', 4]
            diag_sig$pvalUC[i] <- summary(mdl)$coefficients['diagnosisUC', 4]
            diag_sig$pvalActivityIndex[i] <- summary(mdl)$coefficients['activity_index', 4]
        }
    }


    diag_sig$qvalCD <- p.adjust(diag_sig$pvalCD, 'fdr')
    diag_sig$qvalUC <- p.adjust(diag_sig$pvalUC, 'fdr')
    diag_sig$qvalActivityIndex <- p.adjust(diag_sig$pvalActivityIndex, 'fdr')


    # Separate diagnosis and activity
    diagnosis_only<-c( 'coefCD', 'coefUC', 'zvalCD', 'zvalUC', 'pvalCD', 'pvalUC', 'qvalCD', 'qvalUC')
    diag_sig_diagnosis<-diag_sig[, c('prevalence', diagnosis_only)]
    diag_sig_activity<-diag_sig[, colnames(diag_sig)[!colnames(diag_sig) %in% diagnosis_only]]
    diag_sig_activity <- diag_sig_activity[order(diag_sig_activity$qvalActivityIndex),]

    # Order by minimum q-value
    diag_sig_diagnosis$minQ<-pmin(diag_sig_diagnosis$qvalCD, diag_sig_diagnosis$qvalUC)
    diag_sig_diagnosis <- diag_sig_diagnosis[order(diag_sig_diagnosis$minQ),]

    # Remove row names
    diag_sig_diagnosis<-rownames_to_column(diag_sig_diagnosis, '#')
    diag_sig_activity<-rownames_to_column(diag_sig_activity, '#')

    # Remove NA's
    library(tidyr)
    diag_sig_diagnosis<- diag_sig_diagnosis %>% drop_na()
    diag_sig_activity<- diag_sig_activity %>% drop_na()

    # Dump tables
    write.table(diag_sig_activity, file=paste(HMP2_root,
                                              sprintf("/analysis/differential_abundance/%s/active_vs_inactive/activity_index/Activity_Index_IBD_relations_%s.tsv",  transformation, name), sep=''),
                row.names=F,col.names=T,quote=F,sep="\t")
    write.table(diag_sig_diagnosis, file=paste(HMP2_root,
                                               sprintf("/analysis/differential_abundance/%s/diagnosis/activity_index/IBD_relations2_%s.tsv",  transformation, name), sep=''),
                row.names=F,col.names=T,quote=F,sep="\t")
}


run_baseline_DA<-function(pcl, name, transformation) {

    cat(sprintf('Running a combination of %s and %s...\n', name, transformation))

    # Merge disease activity
    if (name != 'bugs') {
        pcl<-merge_disease_activity(pcl, lenience = 0)
    }

    # Add antibiotics
    pcl <- merge_metadata(pcl, fields=c("Antibiotics"))

    # Filter out features with no variance or with >90% zeros
    pcl <- pcl %>% pcl.filter.f(any(is.finite(x)) && var(x, na.rm=T) > 0 && mean(!is.na(x) & (x>0)) >= 0.1)

    # Filter out features with little variance (except for metabolites)
    if (!name %in% c('metabolites', 'metabolites_Sena','metabolites_all','metabolites_nadine')){
        sds <- pcl.apply.f(pcl, sd(x, na.rm=T))
        pcl <- pcl %>% pcl.filter.f(keep=sds > median(sds) / 2)
    }

    # Save prevalence
    diag_sig <- data.frame(prevalence = pcl %>% pcl.apply.f(mean(x>0, na.rm=T)))

    # Transformation
    pcl<-pcl.transform(pcl, name, transformation)

    # Initialize
    diag_sig$coefCD <- NA
    diag_sig$coefUC <- NA
    diag_sig$tvalCD <- NA
    diag_sig$tvalUC <- NA
    diag_sig$pvalCD <- NA
    diag_sig$pvalUC <- NA

    # Extract relevant metadata
    df0 <- pcl$meta
    df0<-df0[order(df0$collection),] # Order by collection to extract first sample

    # Extract baseline subjects
    if (!name %in% c('metabolites', 'metabolites_Sena','metabolites_all','metabolites_nadine')){
        df1<-df0[!duplicated(df0$subject),]
        df<-df1[df1$collection<3,]
    } else{
        df1<-df0[!duplicated(df0$subject), c("diagnosis", "Antibiotics", "consent_age", "site_name", "collection")]
        df<-df1[df1$collection<3,]
    }

    # Remove Empty Column
    df<-df[,colnames(df)!=""]

    # Formula for fixed effects
    formula = as.formula(y ~ diagnosis + Antibiotics  + consent_age)

    # Fit full model
    for (i in seq_along(colnames(pcl$x))) {
        df$y <- pcl$x[rownames(df),i]
        mdl <- try(lme(formula, random= ~ 1 | site_name, data=df, na.action = na.omit))

        if (class(mdl)=='try-error') {
            diag_sig$coefCD[i] <- NA
            diag_sig$coefUC[i] <- NA
            diag_sig$tvalCD[i] <- NA
            diag_sig$tvalUC[i] <- NA
            diag_sig$pvalCD[i] <- NA
            diag_sig$pvalUC[i] <- NA
        }
        else{

            if(all(c("diagnosisUC", "diagnosisCD") %in% rownames(summary(mdl)$tTable))){

                # Extract relevant coefficients
                diag_sig$coefCD[i]  <- summary(mdl)$tTable['diagnosisCD', 'Value']
                diag_sig$coefUC[i] <- summary(mdl)$tTable['diagnosisUC', 'Value']
                diag_sig$tvalCD[i] <-  summary(mdl)$tTable['diagnosisCD', 't-value']
                diag_sig$tvalUC[i] <- summary(mdl)$tTable['diagnosisUC', 't-value']
                diag_sig$pvalCD[i] <- summary(mdl)$tTable['diagnosisCD', 'p-value']
                diag_sig$pvalUC[i] <- summary(mdl)$tTable['diagnosisUC', 'p-value']

            }
            else{
                diag_sig$coefCD[i] <- NA
                diag_sig$coefUC[i] <- NA
                diag_sig$tvalCD[i] <- NA
                diag_sig$tvalUC[i] <- NA
                diag_sig$pvalCD[i] <- NA
                diag_sig$pvalUC[i] <- NA
            }
        }
    }

    # Multiplicity correction
    diag_sig$qvalCD <- p.adjust(diag_sig$pvalCD, 'fdr')
    diag_sig$qvalUC <- p.adjust(diag_sig$pvalUC, 'fdr')

    # Order by minimum q-value
    diag_sig$minQ<-pmin(diag_sig$qvalCD, diag_sig$qvalUC)
    diag_sig <- diag_sig[order(diag_sig$minQ),]

    # Remove row names
    diag_sig<-rownames_to_column(diag_sig, '#')

    # Remove NA's
    library(tidyr)
    diag_sig<- diag_sig %>% drop_na()

    # Dump tables
    write.table(diag_sig, file=paste(HMP2_root,
                                               sprintf("/analysis/differential_abundance/%s/diagnosis/baseline/Baseline_IBD_relations_%s.tsv",  transformation, name), sep=''),
                row.names=F,col.names=T,quote=F,sep="\t")

}

run_baseline_DA_binomial<-function(pcl, name, transformation) {

    cat(sprintf('Running a combination of %s and %s...\n', name, transformation))

    # Merge disease activity
    pcl<-merge_disease_activity(pcl, lenience = 0)

    # Add antibiotics
    pcl <- merge_metadata(pcl, fields=c("Antibiotics"))

    # Save prevalence
    diag_sig <- data.frame(prevalence = pcl %>% pcl.apply.f(mean(x>0, na.rm=T)))
    uniquenames<-make.names(colnames(pcl$x), unique=TRUE)
    rownames(diag_sig) <- uniquenames

    # Initialize
    diag_sig$coefCD <- NA
    diag_sig$coefUC <- NA
    diag_sig$zvalCD <- NA
    diag_sig$zvalUC <- NA
    diag_sig$pvalCD <- NA
    diag_sig$pvalUC <- NA

    # Extract relevant metadata
    df0 <- pcl$meta
    df0<-df0[order(df0$collection),] # Order by collection to extract first sample

    # Extract baseline subjects
    df1<-df0[!duplicated(df0$subject),]
    df<-df1[df1$collection<3,]

    # Remove Empty Column
    df<-df[,colnames(df)!=""]

    # Fit full model
    for (i in seq_along(colnames(pcl$x))) {
        df$y <- pcl$x[rownames(df),i]
        mdl <- try(glmer(y ~ diagnosis + Antibiotics + consent_age + (1 | site_name), data=df, family=binomial))

        if (class(mdl)=='try-error') {
            diag_sig$coefCD[i] <- NA
            diag_sig$coefUC[i] <- NA
            diag_sig$zvalCD[i] <- NA
            diag_sig$zvalUC[i] <- NA
            diag_sig$pvalCD[i] <- NA
            diag_sig$pvalUC[i] <- NA
        }
        else{

            # Extract relevant coefficients
            diag_sig$coefCD[i]  <- summary(mdl)$coefficients['diagnosisCD', 1]
            diag_sig$coefUC[i] <- summary(mdl)$coefficients['diagnosisUC', 1]
            diag_sig$zvalCD[i] <-  summary(mdl)$coefficients['diagnosisCD', 3]
            diag_sig$zvalUC[i] <-summary(mdl)$coefficients['diagnosisUC', 3]
            diag_sig$pvalCD[i] <- summary(mdl)$coefficients['diagnosisCD', 4]
            diag_sig$pvalUC[i] <- summary(mdl)$coefficients['diagnosisUC', 4]
        }
    }

    # Multiplicity correction
    diag_sig$qvalCD <- p.adjust(diag_sig$pvalCD, 'fdr')
    diag_sig$qvalUC <- p.adjust(diag_sig$pvalUC, 'fdr')

    # Order by minimum q-value
    diag_sig$minQ<-pmin(diag_sig$qvalCD, diag_sig$qvalUC)
    diag_sig <- diag_sig[order(diag_sig$minQ),]

    # Remove row Names
    diag_sig<-rownames_to_column(diag_sig, '#')

    # Remove NA's
    library(tidyr)
    diag_sig<- diag_sig %>% drop_na()

    # Dump tables
    write.table(diag_sig, file=paste(HMP2_root,
                                     sprintf("/analysis/differential_abundance/%s/diagnosis/baseline/Baseline_IBD_relations_%s.tsv",  transformation, name), sep=''),
                row.names=F,col.names=T,quote=F,sep="\t")

}


# Filter KOs (flag KOs highly correlated with a single bug)
filterKO<-function(KO, bugs, lenience=4, matching=F, correlationThreshold=0.6){



    # Basic filtering (KO)
    KO <- KO %>% pcl.filter.f(any(is.finite(x)) && var(x, na.rm=T) > 0 && mean(!is.na(x) & (x>0)) >= 0.1)

    # Filter out features with little variance (KO)
    sdsKO <- pcl.apply.f(KO, sd(x, na.rm=T))
    KO <- KO %>% pcl.filter.f(keep=sdsKO > median(sdsKO) / 2)

    # Basic filtering (bugs)
    bugs <- bugs %>% pcl.filter.f(any(is.finite(x)) && var(x, na.rm=T) > 0 && mean(!is.na(x) & (x>0)) >= 0.1)

    # Filter out features with little variance (bugs)
    sdsbugs <- pcl.apply.f(bugs, sd(x, na.rm=T))
    bugs <- bugs %>% pcl.filter.f(keep=sdsbugs > median(sdsKO) / 2)

    # Match datasets
    matchedKObugs<-match_datasets(list(KO, bugs), lenience=lenience, matching=matching)
    matchedKO<-matchedKObugs[[1]]$x
    matchedbugs<-matchedKObugs[[2]]$x


    # Pairwise correlations
    cat ('Creating filtered KO table...\n')
    correlationKObugs<-matrix(nrow=ncol(matchedKO), ncol=ncol(matchedbugs))
    for (i in 1:ncol(matchedKO)){
        for (j in 1:ncol(matchedbugs)){
            correlationKObugs[i,j]<-cor(as.numeric(as.vector(matchedKO[,i])), as.numeric(matchedbugs[,j]), method="spearman")
        }
        }
    cat ('Filtered KO table created\n')
    colnames(correlationKObugs)<-colnames(matchedbugs)
    rownames(correlationKObugs)<-colnames(matchedKO)

    # Flag KOs highly correlated with bug abundances
    col_sub<-apply(correlationKObugs, 1, function(row) any(abs(row) > correlationThreshold))
    KO$x<-KO$x[,col_sub]

    # Retain same PCL structure
    KO$ns<-dim(KO$x)[1]
    KO$nf<-dim(KO$x)[2]

    # Return filtered KO table
    return(KO)
}


# Create PCL for RNA/DNA ratios from KO bug roll ups
createPCL<-function(RNA, DNA, lenience=2, matching=F, prevalenceThreshold=0, abundanceThreshold=0.1){

    # Match samples
    matchedbugsKO<-match_datasets(list(RNA, DNA),lenience=lenience, matching=matching)
    matchedbugsKORNA.pcl<-matchedbugsKO[[1]]
    matchedbugsKODNA.pcl<-matchedbugsKO[[2]]

    # Match features after some basic trimming (10% non-zero)
    matchedbugsKORNA.pcl<-matchedbugsKORNA.pcl %>% pcl.filter.f(keep=colSums(matchedbugsKORNA.pcl$x > prevalenceThreshold) > nrow(matchedbugsKORNA.pcl$x )*abundanceThreshold)
    matchedbugsKODNA.pcl<-matchedbugsKODNA.pcl %>% pcl.filter.f(keep=colSums(matchedbugsKODNA.pcl$x > prevalenceThreshold) > nrow(matchedbugsKODNA.pcl$x )*abundanceThreshold)
    common_features<-intersect(colnames(matchedbugsKORNA.pcl$x), colnames(matchedbugsKODNA.pcl$x))
    matchedbugsKORNA.pcl<-matchedbugsKORNA.pcl %>% pcl.filter.f(keep=common_features)
    matchedbugsKODNA.pcl<-matchedbugsKODNA.pcl %>% pcl.filter.f(keep=common_features)

    # PCL file (retain same structure for DA testing)
    bugkoRNADNAs.pcl<-list()
    bugkoRNADNAs.pcl$x<-matchedbugsKORNA.pcl$x/matchedbugsKODNA.pcl$x
    bugkoRNADNAs.pcl$meta<-matchedbugsKORNA.pcl$meta
    bugkoRNADNAs.pcl$ns<-dim(bugkoRNADNAs.pcl$x)[1]
    bugkoRNADNAs.pcl$nf<-dim(bugkoRNADNAs.pcl$x)[2]

    # Return
    return(bugkoRNADNAs.pcl)
}



run_DA_PreActive<-function(pcl, name, PreActive, transformation) {

    cat(sprintf('Running a combination of %s and %s...\n', name, transformation))

    # Merge disease activity
    if (name != 'bugs') {
        pcl<-merge_disease_activity(pcl, lenience = 0)
    }

    # Add antibiotics
    pcl <- merge_metadata(pcl, fields=c("Antibiotics"))

    # Filter out features with no variance or with >90% zeros
    pcl <- pcl %>% pcl.filter.f(any(is.finite(x)) && var(x, na.rm=T) > 0 && mean(!is.na(x) & (x>0)) >= 0.1)

    # Filter out features with little variance (except for metabolites)
    if (!name %in% c('metabolites', 'metabolites_Sena','metabolites_all','metabolites_nadine')){
        sds <- pcl.apply.f(pcl, sd(x, na.rm=T))
        pcl <- pcl %>% pcl.filter.f(keep=sds > median(sds) / 2)
    }

    # Save prevalence
    diag_sig <- data.frame(prevalence = pcl %>% pcl.apply.f(mean(x>0, na.rm=T)))

    # Transformation
    pcl<-pcl.transform(pcl, name, transformation)

    # Initialize
    diag_sig$coefPreActiveCD <- NA
    diag_sig$coefPreActiveUC <- NA
    diag_sig$tvalPreActiveCD <- NA
    diag_sig$tvalPreActiveUC <- NA
    diag_sig$pvalPreActiveCD <- NA
    diag_sig$pvalPreActiveUC <- NA

    # Extract metadata
    df0 <- pcl$meta

    # Merge with Preactive/Inactive information
    Intersect<-intersect(rownames(df0), rownames(PreActive))
    df<-df0[Intersect,]
    PreActive<-PreActive[Intersect,]
    df$PreActive<-PreActive
    df$PreActive<-ifelse(df$PreActive=='PreActive',  '1', '0')
    df$PreActive<-as.factor(df$PreActive)
    df$PreActive <- relevel(df$PreActive, ref = "0")

    # Remove Empty Column
    df<-df[,colnames(df)!=""]

    # Fit model (CD)
    dfCD<-subset(df, df$diagnosis=='CD')

    for (i in seq_along(colnames(pcl$x))) {
        dfCD$y <- pcl$x[rownames(dfCD),i]
        mdl <- try(lme(y ~ PreActive + Antibiotics + consent_age, random= ~ 1 | site_name/subject, data=dfCD, na.action = na.omit))


        if (class(mdl)=='try-error') {
            diag_sig$coefPreActiveCD[i] <- NA
            diag_sig$tvalPreActiveCD[i] <- NA
            diag_sig$pvalPreActiveCD[i] <- NA
        }
        else{

            # How significant is the diagnosis difference?
            diag_sig$coefPreActiveCD[i] <- summary(mdl)$tTable['PreActive1', 'Value']
            diag_sig$tvalPreActiveCD[i] <- summary(mdl)$tTable['PreActive1', 't-value']
            diag_sig$pvalPreActiveCD[i] <- summary(mdl)$tTable['PreActive1', 'p-value']
        }
    }

    # Fit model (UC)
    dfUC<-subset(df, df$diagnosis=='UC') # Model Fitting

    for (i in seq_along(colnames(pcl$x))) {

        # Fit the model

        dfUC$y <- pcl$x[rownames(dfUC),i]
        mdl <- try(lme(y ~ PreActive + Antibiotics + consent_age, random= ~ 1 | site_name/subject, data=dfUC, na.action = na.omit))

        if (class(mdl)=='try-error') {
            diag_sig$coefPreActiveUC[i] <- NA
            diag_sig$tvalPreActiveUC[i] <- NA
            diag_sig$pvalPreActiveUC[i] <- NA
        }
        else{

            # How significant is the diagnosis difference?
            diag_sig$coefPreActiveUC[i] <- summary(mdl)$tTable['PreActive1', 'Value']
            diag_sig$tvalPreActiveUC[i] <- summary(mdl)$tTable['PreActive1', 't-value']
            diag_sig$pvalPreActiveUC[i] <- summary(mdl)$tTable['PreActive1', 'p-value']
        }
    }

    # Multiplicity correction
    diag_sig$qvalPreActiveCD <- p.adjust(diag_sig$pvalPreActiveCD, 'fdr')
    diag_sig$qvalPreActiveUC <- p.adjust(diag_sig$pvalPreActiveUC, 'fdr')

    # Order by minimum q-value
    diag_sig$minQ<-pmin(diag_sig$qvalPreActiveCD, diag_sig$qvalPreActiveUC)
    diag_sig <- diag_sig[order(diag_sig$minQ),]

    # Remove row names
    diag_sig<-rownames_to_column(diag_sig, '#')

    # Remove NA's
    library(tidyr)
    diag_sig<- diag_sig<- diag_sig %>% drop_na()

    # Save
    write.table(diag_sig, file=paste(HMP2_root,
                                     sprintf("/analysis/differential_abundance/%s/preactive_vs_inactive/PreActive_IBD_relations_%s.tsv",  transformation, name), sep=''),
                row.names=F,col.names=T,quote=F,sep="\t")

}

run_DA_PreActive_binomial<-function(pcl, name, PreActive, transformation) {

    cat(sprintf('Running a combination of %s and %s...\n', name, transformation))

    # Merge disease activity
    pcl<-merge_disease_activity(pcl, lenience = 0)

    # Add antibiotics
    pcl <- merge_metadata(pcl, fields=c("Antibiotics"))

    # Save prevalence
    diag_sig <- data.frame(prevalence = pcl %>% pcl.apply.f(mean(x>0, na.rm=T)))
    uniquenames<-make.names(colnames(pcl$x), unique=TRUE)
    rownames(diag_sig) <- uniquenames

    # Intialize
    diag_sig$coefPreActiveCD <- NA
    diag_sig$coefPreActiveUC <- NA
    diag_sig$zvalPreActiveCD <- NA
    diag_sig$zvalPreActiveUC <- NA
    diag_sig$pvalPreActiveCD <- NA
    diag_sig$pvalPreActiveUC <- NA

    # Extract relevant metadata
    df0 <- pcl$meta

    # Merge with Preactive/Inactive information
    Intersect<-intersect(rownames(df0), rownames(PreActive))
    df<-df0[Intersect,]
    PreActive<-PreActive[Intersect,]
    df$PreActive<-PreActive
    df$PreActive<-ifelse(df$PreActive=='PreActive',  '1', '0')
    df$PreActive<-as.factor(df$PreActive)
    df$PreActive <- relevel(df$PreActive, ref = "0")

    # Remove Empty Column
    df<-df[,colnames(df)!=""]

    # Fit model (CD)
    dfCD<-subset(df, df$diagnosis=='CD')

    for (i in seq_along(colnames(pcl$x))) {
        dfCD$y <- pcl$x[rownames(dfCD),i]
        mdl <- try(glmer(y ~ PreActive + Antibiotics + consent_age + ( 1 | site/subject), data=dfCD, family=binomial))

        if (class(mdl)=='try-error') {
            diag_sig$coefPreActiveCD[i] <- NA
            diag_sig$zvalPreActiveCD[i] <- NA
            diag_sig$pvalPreActiveCD[i] <- NA
        }
        else{

            # How significant is the diagnosis difference?
            diag_sig$coefPreActiveCD[i] <- summary(mdl)$coefficients['PreActive1', 1]
            diag_sig$zvalPreActiveCD[i] <- summary(mdl)$coefficients['PreActive1', 3]
            diag_sig$pvalPreActiveCD[i] <- summary(mdl)$coefficients['PreActive1', 4]
        }
    }

    # Fit model (UC)
    dfUC<-subset(df, df$diagnosis=='UC')  # Model Fitting

    for (i in seq_along(colnames(pcl$x))) {
        dfUC$y <- pcl$x[rownames(dfUC),i]
        mdl <- try(glmer(y ~ PreActive + Antibiotics + consent_age + ( 1 | site/subject), data=dfUC, family=binomial))

        if (class(mdl)=='try-error') {
            diag_sig$coefPreActiveUC[i] <- NA
            diag_sig$zvalPreActiveUC[i] <- NA
            diag_sig$pvalPreActiveUC[i] <- NA
        }
        else{

            # How significant is the diagnosis difference?
            diag_sig$coefPreActiveUC[i] <- summary(mdl)$coefficients['PreActive1', 1]
            diag_sig$zvalPreActiveUC[i] <- summary(mdl)$coefficients['PreActive1', 3]
            diag_sig$pvalPreActiveUC[i] <- summary(mdl)$coefficients['PreActive1', 4]
        }
    }

    # Multiplicity correction
    diag_sig$qvalPreActiveCD <- p.adjust(diag_sig$pvalPreActiveCD, 'fdr')
    diag_sig$qvalPreActiveUC <- p.adjust(diag_sig$pvalPreActiveUC, 'fdr')

    # Order by minimum q-value
    diag_sig$minQ<-pmin(diag_sig$qvalPreActiveCD, diag_sig$qvalPreActiveUC)
    diag_sig <- diag_sig[order(diag_sig$minQ),]

    # Remove row names
    diag_sig<-rownames_to_column(diag_sig, '#')

    # Remove NA's
    library(tidyr)
    diag_sig<- diag_sig %>% drop_na()

    # Save
    write.table(diag_sig, file=paste(HMP2_root,
                                     sprintf("/analysis/differential_abundance/%s/preactive_vs_inactive/PreActive_IBD_relations_%s.tsv",  transformation, name), sep=''),
                row.names=F,col.names=T,quote=F,sep="\t")

}

###################
# Transformations #
###################

pcl.transform<-function(pcl, name, transformation){

    if (transformation == 'CLR_perfeature'){
        if (name %in% c('bugs', 'ecDNAs', 'ecRNAs', 'pwyDNAs', 'pwyRNAs', 'koDNAs', 'koRNAs')){
            pcl<-  pcl %>% CLR_perfeature
        }
        if (name %in% c('metabolites', 'metabolites_Sena','metabolites_all','metabolites_nadine','proteins', 'proteinECs', 'proteinKOs')){
            pcl$x<- t(SpiecEasi::clr(pcl$x + 1, mar = 1, base = 2)) # counts
        }
        if (name %in% c('pwyRNADNAs','bugkoRNADNAs')){
            pcl$x <- log2(pcl$x) # ratios
            pcl$x[!is.finite(pcl$x)] <- NA
        }
    }

    if (transformation == 'CLR_persample'){
        if (name %in% c('bugs', 'ecDNAs', 'ecRNAs', 'pwyDNAs', 'pwyRNAs', 'koDNAs', 'koRNAs')){
            pcl<-  pcl %>% CLR_persample
        }
        if (name %in% c('metabolites', 'metabolites_Sena','metabolites_all','metabolites_nadine', 'proteins', 'proteinECs', 'proteinKOs')){
            pcl$x<- t(SpiecEasi::clr(pcl$x + 1, mar = 1, base = 2)) # counts
        }
        if (name %in% c('pwyRNADNAs','bugkoRNADNAs')){
            pcl$x <- log2(pcl$x) # ratios
            pcl$x[!is.finite(pcl$x)] <- NA
        }
    }

    if (transformation == 'CLR_overall'){
        if (name %in% c('bugs', 'ecDNAs', 'ecRNAs', 'pwyDNAs', 'pwyRNAs', 'koDNAs', 'koRNAs')){
            pcl<-  pcl %>% CLR_overall
        }
        if (name %in% c('metabolites', 'metabolites_Sena','metabolites_all','metabolites_nadine','proteins', 'proteinECs', 'proteinKOs')){
            pcl$x<- t(SpiecEasi::clr(pcl$x + 1, mar = 1, base = 2)) # counts
        }
        if (name %in% c('pwyRNADNAs','bugkoRNADNAs')){
            pcl$x <- log2(pcl$x) # ratios
            pcl$x[!is.finite(pcl$x)] <- NA
        }
    }
    if (transformation == 'Vanilla'){

        #######################################
        # Transformation (Ad hoc VST) Summary #
        #######################################
        # Counts -> log (base 2) with pseudocount 1
        # Proportions (between 0 and 1) -> Arc Sine Square Root
        # Ratios -> log (base 2)

        if (name %in% c('bugs', 'ecDNAs', 'ecRNAs', 'pwyDNAs', 'pwyRNAs', 'koDNAs', 'koRNAs')){
            pcl$x <- asin(sqrt(pcl$x)) # proportions

        }
        if (name %in% c('metabolites', 'metabolites_Sena','metabolites_all','metabolites_nadine','proteins', 'proteinECs', 'proteinKOs')){
            pcl$x <- log2(pcl$x + 1) # counts

        }
        if (name %in% c('pwyRNADNAs','bugkoRNADNAs')){
            pcl$x <- log2(pcl$x) # ratios
            pcl$x[!is.finite(pcl$x)] <- NA
        }
        }
    return(pcl)
}


# Function to fit full model for HMP2 and HMP1-II Comparisons
run_DA_HMP1_II<-function(pcl, transformation) {

    # Prepare output file and save prevalence
    diag_sig <- data.frame(prevalence = species.pcl %>% pcl.apply.f(mean(x>0, na.rm=T)))
    diag_sig$coefHMP1II <- NA
    diag_sig$tvalHMP1II <- NA
    diag_sig$pvalHMP1II <- NA

    if (transformation=='CLR_perfeature') {
        pcl<-  pcl %>% CLR_perfeature
    }

    if (transformation=='CLR_persample') {
        pcl<-  pcl %>% CLR_persample
    }

    if (transformation=='CLR_overall') {
        pcl<-  pcl %>% CLR_overall
    }

    if (transformation=='Vanilla') {
        pcl$x <- asin(sqrt(pcl$x))
    }

    for (i in seq_along(colnames(pcl$x))) {
        # Extract relevant metadata
        df <- pcl$meta
        df$y <- pcl$x[rownames(df),i]
        mdl <- try(lme(y ~ diagnosis, random= ~ 1 | merged_subj, data=df, na.action = na.omit))

        if (class(mdl)=='try-error') {
            diag_sig$coefHMP1II[i] <- NA
            diag_sig$tvalHMP1II[i] <- NA
            diag_sig$pvalHMP1II[i] <- NA

        }
        else{

            # Extract relevant coefficients
            diag_sig$coefHMP1II[i]  <- summary(mdl)$tTable['diagnosisHMP1-II', 'Value']
            diag_sig$tvalHMP1II[i] <- summary(mdl)$tTable['diagnosisHMP1-II', 't-value']
            diag_sig$pvalHMP1II[i]  <- summary(mdl)$tTable['diagnosisHMP1-II', 'p-value']
        }
    }

    # Multiplicity correction
    diag_sig$qvalHMP1II <- p.adjust(diag_sig$pvalHMP1II, 'fdr')

    # Order by minimum q-value
    diag_sig <- diag_sig[order(diag_sig$qvalHMP1II),]

    # Remove row names
    diag_sig<-rownames_to_column(diag_sig, '#')

    # Remove NA's
    library(tidyr)
    diag_sig<- diag_sig %>% drop_na()

    # Dump tables
    write.table(diag_sig, file=paste(HMP2_root,
                                     sprintf("/analysis/differential_abundance/HMP1-II/HMP1-II_%s.tsv",  transformation), sep=''),
                row.names=F,col.names=T,quote=F,sep="\t")
}



# Pseudocount-driven CLR Transformations
CLR_perfeature<-function(pcl, multiplier = 0.5){
    pseudocount_perfeature<- pcl %>% pcl.apply.f(min(x[x!=0], na.rm=T))
    pcl$x<- t(SpiecEasi::clr(pcl$x + rep(multiplier*pseudocount_perfeature, each = nrow(pcl$x)), base = 2, mar = 1))
    return(pcl)
}

CLR_persample<-function(pcl, multiplier = 0.5){
    pseudocount_persample<- pcl %>% pcl.apply.s(min(x[x!=0], na.rm=T))
    pcl$x<- t(SpiecEasi::clr(pcl$x + rep(multiplier*pseudocount_persample, ncol(pcl$x)), base = 2, mar = 1))
    return(pcl)
}

CLR_overall<-function(pcl, multiplier = 0.5){
    pseudocount_overall <- min(pcl$x[pcl$x!=0], na.rm=T)
    pcl$x<- t(SpiecEasi::clr(pcl$x + multiplier*pseudocount_overall, base = 2, mar = 1))
    return(pcl)
}

# DA with HBI and other covariates for CD patients only
run_DA_HBI<-function(pcl, name, transformation) {

    cat(sprintf('Running a combination of %s and %s...\n', name, transformation))

    # Merge disease activity
    if (name != 'bugs') {
        pcl<-merge_disease_activity(pcl, lenience = 0)
    }

    # Add antibiotics
    pcl <- merge_metadata(pcl, fields=c("Antibiotics"))

    # Filter out features with no variance or with >90% zeros
    pcl <- pcl %>% pcl.filter.f(any(is.finite(x)) && var(x, na.rm=T) > 0 && mean(!is.na(x) & (x>0)) >= 0.1)

    # Filter out features with little variance (except for metabolites)
    if (!name %in% c('metabolites', 'metabolites_Sena','metabolites_all','metabolites_nadine')){
        sds <- pcl.apply.f(pcl, sd(x, na.rm=T))
        pcl <- pcl %>% pcl.filter.f(keep=sds > median(sds) / 2)
    }

    # Save prevalence
    diag_sig <- data.frame(prevalence = pcl %>% pcl.apply.f(mean(x>0, na.rm=T)))

    # Transformation
    pcl<-pcl.transform(pcl, name, transformation)

    # Initialize
    diag_sig$coefHBI <- NA
    diag_sig$tvalHBI <- NA
    diag_sig$pvalHBI <- NA

    # Extract relevant metadata
    df <- pcl$meta
    df<-df[df$diagnosis=='CD',]

    # Remove Empty Column
    df<-df[,colnames(df)!=""]


    # Formula for fixed effects
    formula = as.formula(y ~ hbi + Antibiotics  + consent_age)

    # Fit full model
    for (i in seq_along(colnames(pcl$x))) {
        df$y <- pcl$x[rownames(df),i]
        mdl <- try(lme(formula , random= ~ 1 | site_name/subject, data=df, na.action = na.omit))

        if (class(mdl)=='try-error') {
            diag_sig$coefHBI[i] <- NA
            diag_sig$tvalHBI[i] <- NA
            diag_sig$pvalHBI[i] <- NA
        }
        else{

            # Extract relevant coefficients
            diag_sig$coefHBI[i]  <- summary(mdl)$tTable['hbi', 'Value']
            diag_sig$tvalHBI[i] <-  summary(mdl)$tTable['hbi', 't-value']
            diag_sig$pvalHBI[i] <- summary(mdl)$tTable['hbi', 'p-value']
        }
    }

    # Multiplicity correction and order
    diag_sig$qvalHBI <- p.adjust(diag_sig$pvalHBI, 'fdr')
    diag_sig<-diag_sig[order(diag_sig$qvalHBI, decreasing = FALSE),]

    # Remove row names
    diag_sig<-rownames_to_column(diag_sig, '#')

    # Remove NA's
    library(tidyr)
    diag_sig<- diag_sig %>% drop_na()

    # Dump tables
    write.table(diag_sig, file=paste(HMP2_root,
                                     sprintf("/analysis/differential_abundance/%s/HBI/HBI_relations_%s.tsv",  transformation, name), sep=''),
                row.names=F,col.names=T,quote=F,sep="\t")
}


run_DA_HBI_binomial<-function(pcl, name, transformation) {

    cat(sprintf('Running a combination of %s and %s...\n', name, transformation))

    # Merge disease activity
    pcl<-merge_disease_activity(pcl, lenience = 0)

    # Add antibiotics
    pcl <- merge_metadata(pcl, fields=c("Antibiotics"))

    # Save prevalence
    diag_sig <- data.frame(prevalence = pcl %>% pcl.apply.f(mean(x>0, na.rm=T)))
    uniquenames<-make.names(colnames(pcl$x), unique=TRUE)
    rownames(diag_sig) <- uniquenames

    # Initialize
    diag_sig$coefHBI <- NA
    diag_sig$zvalHBI <- NA
    diag_sig$pvalHBI <- NA

    # Extract relevant metadata
    df <- pcl$meta
    df<-df[df$diagnosis=='CD',]

    # Remove Empty Column
    df<-df[,colnames(df)!=""]


    # Fit full model
    for (i in seq_along(colnames(pcl$x))) {
        df$y <- pcl$x[rownames(df),i]
        mdl <- try(glmer(y ~ hbi + Antibiotics + consent_age + ( 1 | site_name/subject), data=df, family=binomial))

        if (class(mdl)=='try-error') {
            diag_sig$coefHBI[i] <- NA
            diag_sig$zvalHBI[i] <- NA
            diag_sig$pvalHBI[i] <- NA
        }
        else{

            # Extract relevant coefficients
            diag_sig$coefHBI[i]  <- summary(mdl)$coefficients['hbi', 1]
            diag_sig$zvalHBI[i] <-  summary(mdl)$coefficients['hbi', 3]
            diag_sig$pvalHBI[i] <- summary(mdl)$coefficients['hbi', 4]
        }
    }

    # Multiplicity correction and order
    diag_sig$qvalHBI <- p.adjust(diag_sig$pvalHBI, 'fdr')
    diag_sig<-diag_sig[order(diag_sig$qvalHBI, decreasing = FALSE),]

    # Remove row names
    diag_sig<-rownames_to_column(diag_sig, '#')

    # Remove NA's
    library(tidyr)
    diag_sig<- diag_sig %>% drop_na()

    # Dump tables
    write.table(diag_sig, file=paste(HMP2_root,
                                     sprintf("/analysis/differential_abundance/%s/HBI/HBI_relations_%s.tsv",  transformation, name), sep=''),
                row.names=F,col.names=T,quote=F,sep="\t")
}

# DA with HBIobjective and other covariates for CD patients only
run_DA_HBIobjective<-function(pcl, name, transformation) {

    cat(sprintf('Running a combination of %s and %s...\n', name, transformation))

    # Merge disease activity
    if (name != 'bugs') {
        pcl<-merge_disease_activity(pcl, lenience = 0)
    }

    # Add antibiotics
    pcl <- merge_metadata(pcl, fields=c("Antibiotics", "Number.of.liquid.or.very.soft.stools.in.the.past.24.hours."))

    # Filter out features with no variance or with >90% zeros
    pcl <- pcl %>% pcl.filter.f(any(is.finite(x)) && var(x, na.rm=T) > 0 && mean(!is.na(x) & (x>0)) >= 0.1)

    # Filter out features with little variance (except for metabolites)
    if (!name %in% c('metabolites', 'metabolites_Sena','metabolites_all','metabolites_nadine')){
        sds <- pcl.apply.f(pcl, sd(x, na.rm=T))
        pcl <- pcl %>% pcl.filter.f(keep=sds > median(sds) / 2)
    }

    # Save prevalence
    diag_sig <- data.frame(prevalence = pcl %>% pcl.apply.f(mean(x>0, na.rm=T)))

    # Transformation
    pcl<-pcl.transform(pcl, name, transformation)

    # Initialize
    diag_sig$coefHBIobjective <- NA
    diag_sig$tvalHBIobjective <- NA
    diag_sig$pvalHBIobjective <- NA

    # Extract relevant metadata
    df <- pcl$meta
    df<-df[df$diagnosis=='CD',]

    # Remove Empty Column
    df<-df[,colnames(df)!=""]


    # Formula for fixed effects
    formula = as.formula(y ~ Number.of.liquid.or.very.soft.stools.in.the.past.24.hours. + Antibiotics  + consent_age)

    # Fit full model
    for (i in seq_along(colnames(pcl$x))) {
        df$y <- pcl$x[rownames(df),i]
        mdl <- try(lme(formula , random= ~ 1 | site_name/subject, data=df, na.action = na.omit))

        if (class(mdl)=='try-error') {
            diag_sig$coefHBIobjective[i] <- NA
            diag_sig$tvalHBIobjective[i] <- NA
            diag_sig$pvalHBIobjective[i] <- NA
        }
        else{

            # Extract relevant coefficients
            diag_sig$coefHBIobjective[i]  <- summary(mdl)$tTable['Number.of.liquid.or.very.soft.stools.in.the.past.24.hours.', 'Value']
            diag_sig$tvalHBIobjective[i] <-  summary(mdl)$tTable['Number.of.liquid.or.very.soft.stools.in.the.past.24.hours.', 't-value']
            diag_sig$pvalHBIobjective[i] <- summary(mdl)$tTable['Number.of.liquid.or.very.soft.stools.in.the.past.24.hours.', 'p-value']
        }
    }

    # Multiplicity correction and order
    diag_sig$qvalHBIobjective <- p.adjust(diag_sig$pvalHBIobjective, 'fdr')
    diag_sig<-diag_sig[order(diag_sig$qvalHBIobjective, decreasing = FALSE),]

    # Remove row names
    diag_sig<-rownames_to_column(diag_sig, '#')

    # Remove NA's
    library(tidyr)
    diag_sig<- diag_sig %>% drop_na()

    # Dump tables
    write.table(diag_sig, file=paste(HMP2_root,
                                     sprintf("/analysis/differential_abundance/%s/HBIobjective/HBIobjective_relations_%s.tsv",  transformation, name), sep=''),
                row.names=F,col.names=T,quote=F,sep="\t")
}


run_DA_HBIobjective_binomial<-function(pcl, name, transformation) {

    cat(sprintf('Running a combination of %s and %s...\n', name, transformation))

    # Merge disease activity
    pcl<-merge_disease_activity(pcl, lenience = 0)

    # Add antibiotics
    pcl <- merge_metadata(pcl, fields=c("Antibiotics", "Number.of.liquid.or.very.soft.stools.in.the.past.24.hours."))

    # Save prevalence
    diag_sig <- data.frame(prevalence = pcl %>% pcl.apply.f(mean(x>0, na.rm=T)))
    uniquenames<-make.names(colnames(pcl$x), unique=TRUE)
    rownames(diag_sig) <- uniquenames

    # Initialize
    diag_sig$coefHBIobjective <- NA
    diag_sig$zvalHBIobjective <- NA
    diag_sig$pvalHBIobjective <- NA

    # Extract relevant metadata
    df <- pcl$meta
    df<-df[df$diagnosis=='CD',]

    # Remove Empty Column
    df<-df[,colnames(df)!=""]

    # Fit full model
    for (i in seq_along(colnames(pcl$x))) {
        df$y <- pcl$x[rownames(df),i]
        mdl <- try(glmer(y ~ Number.of.liquid.or.very.soft.stools.in.the.past.24.hours. + Antibiotics + consent_age + ( 1 | site_name/subject), data=df, family=binomial))

        if (class(mdl)=='try-error') {
            diag_sig$coefHBIobjective[i] <- NA
            diag_sig$zvalHBIobjective[i] <- NA
            diag_sig$pvalHBIobjective[i] <- NA
        }
        else{

            # Extract relevant coefficients
            diag_sig$coefHBIobjective[i]  <- summary(mdl)$coefficients['Number.of.liquid.or.very.soft.stools.in.the.past.24.hours.', 1]
            diag_sig$zvalHBIobjective[i] <-  summary(mdl)$coefficients['Number.of.liquid.or.very.soft.stools.in.the.past.24.hours.', 3]
            diag_sig$pvalHBIobjective[i] <- summary(mdl)$coefficients['Number.of.liquid.or.very.soft.stools.in.the.past.24.hours.', 4]
        }
    }

    # Multiplicity correction and order
    diag_sig$qvalHBIobjective <- p.adjust(diag_sig$pvalHBIobjective, 'fdr')
    diag_sig<-diag_sig[order(diag_sig$qvalHBIobjective, decreasing = FALSE),]

    # Remove row names
    diag_sig<-rownames_to_column(diag_sig, '#')

    # Remove NA's
    library(tidyr)
    diag_sig<- diag_sig %>% drop_na()

    # Dump tables
    write.table(diag_sig, file=paste(HMP2_root,
                                     sprintf("/analysis/differential_abundance/%s/HBIobjective/HBIobjective_relations_%s.tsv",  transformation, name), sep=''),
                row.names=F,col.names=T,quote=F,sep="\t")
}


# DA with SCCAI and other covariates for UC patients only
run_DA_SCCAI<-function(pcl, name, transformation) {

    cat(sprintf('Running a combination of %s and %s...\n', name, transformation))

    # Merge disease activity
    if (name != 'bugs') {
        pcl<-merge_disease_activity(pcl, lenience = 0)
    }

    # Add antibiotics
    pcl <- merge_metadata(pcl, fields=c("Antibiotics"))

    # Filter out features with no variance or with >90% zeros
    pcl <- pcl %>% pcl.filter.f(any(is.finite(x)) && var(x, na.rm=T) > 0 && mean(!is.na(x) & (x>0)) >= 0.1)

    # Filter out features with little variance (except for metabolites)
    if (!name %in% c('metabolites', 'metabolites_Sena','metabolites_all','metabolites_nadine')){
        sds <- pcl.apply.f(pcl, sd(x, na.rm=T))
        pcl <- pcl %>% pcl.filter.f(keep=sds > median(sds) / 2)
    }

    # Save prevalence
    diag_sig <- data.frame(prevalence = pcl %>% pcl.apply.f(mean(x>0, na.rm=T)))

    # Transformation
    pcl<-pcl.transform(pcl, name, transformation)

    # Initialize
    diag_sig$coefSCCAI <- NA
    diag_sig$tvalSCCAI <- NA
    diag_sig$pvalSCCAI <- NA

    # Extract relevant metadata
    df <- pcl$meta
    df<-df[df$diagnosis=='UC',]

    # Remove Empty Column
    df<-df[,colnames(df)!=""]

    # Formula for fixed effects
    formula = as.formula(y ~ sccai + Antibiotics  + consent_age)

    # Fit full model
    for (i in seq_along(colnames(pcl$x))) {
        df$y <- pcl$x[rownames(df),i]
        mdl <- try(lme(formula , random= ~ 1 | site_name/subject, data=df, na.action = na.omit))

        if (class(mdl)=='try-error') {
            diag_sig$coefSCCAI[i] <- NA
            diag_sig$tvalSCCAI[i] <- NA
            diag_sig$pvalSCCAI[i] <- NA
        }
        else{

            # Extract relevant coefficients
            diag_sig$coefSCCAI[i]  <- summary(mdl)$tTable['sccai', 'Value']
            diag_sig$tvalSCCAI[i] <-  summary(mdl)$tTable['sccai', 't-value']
            diag_sig$pvalSCCAI[i] <- summary(mdl)$tTable['sccai', 'p-value']
        }
    }

    # Multiplicity correction and order
    diag_sig$qvalSCCAI <- p.adjust(diag_sig$pvalSCCAI, 'fdr')
    diag_sig<-diag_sig[order(diag_sig$qvalSCCAI, decreasing = FALSE),]

    # Remove row names
    diag_sig<-rownames_to_column(diag_sig, '#')

    # Remove NA's
    library(tidyr)
    diag_sig<- diag_sig %>% drop_na()

    # Dump tables
    write.table(diag_sig, file=paste(HMP2_root,
                                     sprintf("/analysis/differential_abundance/%s/SCCAI/SCCAI_relations_%s.tsv",  transformation, name), sep=''),
                row.names=F,col.names=T,quote=F,sep="\t")
}


run_DA_SCCAI_binomial<-function(pcl, name, transformation) {

    cat(sprintf('Running a combination of %s and %s...\n', name, transformation))

    # Merge disease activity
    pcl<-merge_disease_activity(pcl, lenience = 0)

    # Add antibiotics
    pcl <- merge_metadata(pcl, fields=c("Antibiotics"))

    # Save prevalence
    diag_sig <- data.frame(prevalence = pcl %>% pcl.apply.f(mean(x>0, na.rm=T)))
    uniquenames<-make.names(colnames(pcl$x), unique=TRUE)
    rownames(diag_sig) <- uniquenames

    # Initialize
    diag_sig$coefSCCAI <- NA
    diag_sig$zvalSCCAI <- NA
    diag_sig$pvalSCCAI <- NA

    # Extract relevant metadata
    df <- pcl$meta
    df<-df[df$diagnosis=='UC',]

    # Remove Empty Column
    df<-df[,colnames(df)!=""]


    # Fit full model
    for (i in seq_along(colnames(pcl$x))) {
        df$y <- pcl$x[rownames(df),i]
        mdl <- try(glmer(y ~ sccai + Antibiotics + consent_age + ( 1 | site_name/subject), data=df, family=binomial))

        if (class(mdl)=='try-error') {
            diag_sig$coefSCCAI[i] <- NA
            diag_sig$zvalSCCAI[i] <- NA
            diag_sig$pvalSCCAI[i] <- NA
        }
        else{

            # Extract relevant coefficients
            diag_sig$coefSCCAI[i]  <- summary(mdl)$coefficients['sccai', 1]
            diag_sig$zvalSCCAI[i] <-  summary(mdl)$coefficients['sccai', 3]
            diag_sig$pvalSCCAI[i] <- summary(mdl)$coefficients['sccai', 4]
        }
    }

    # Multiplicity correction and order
    diag_sig$qvalSCCAI <- p.adjust(diag_sig$pvalSCCAI, 'fdr')
    diag_sig<-diag_sig[order(diag_sig$qvalSCCAI, decreasing = FALSE),]

    # Remove row names
    diag_sig<-rownames_to_column(diag_sig, '#')

    # Remove NA's
    library(tidyr)
    diag_sig<- diag_sig %>% drop_na()

    # Dump tables
    write.table(diag_sig, file=paste(HMP2_root,
                                     sprintf("/analysis/differential_abundance/%s/SCCAI/SCCAI_relations_%s.tsv",  transformation, name), sep=''),
                row.names=F,col.names=T,quote=F,sep="\t")
}





