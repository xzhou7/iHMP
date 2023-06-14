
source("./common/disease_colors.r")
source("./common/disease_activity.r")
source("./common/theme_nature.r")

library(ggplot2)
library(egg)
def_plot <- ggplot(data=bugs.pcl$meta, aes(x=activity_index, color=diagnosis, fill=diagnosis)) +
    geom_density(alpha=0.3, size=0.25) +
    geom_vline(xintercept=disease_activity_threshold, size=0.5) +
    scale_color_manual(values=hmp2_disease_colors, name="Diagnosis") +
    scale_fill_manual(values=hmp2_disease_colors, name="Diagnosis") +
    theme_nature() + xlab("Dysbiosis score") + ylab("Density") +
    guides(alpha="none", fill="none", color="none")
act_cor_plot <- ggplot(data=bugs.pcl$meta, aes(y=fecalcal, x=activity_index, color=diagnosis, fill=diagnosis)) +
    geom_point(size=0.2) + geom_smooth(method="lm", lwd=0.35) +
    theme_nature() + ylab("Calprotectin (ug/g)") + xlab("Dysbiosis Score") +
    scale_fill_manual(values=hmp2_disease_colors) +
    scale_color_manual(values=hmp2_disease_colors) +
    guides(color="none", fill="none")


# Paper plot

pdf("./disease_activity/definition_cor_paperfigure.pdf", 1.565, 2.0, onefile=F)
ggarrange(def_plot, act_cor_plot, ncol=1)
dev.off()

# What is the threshold?
disease_activity_threshold

# Fraction of samples at 90%
mean(bugs.pcl$meta$activity_index[bugs.pcl$meta$diagnosis=="CD"] > quantile(bugs.pcl$meta$activity_index[bugs.pcl$meta$diagnosis=="nonIBD"], 0.9))
mean(bugs.pcl$meta$activity_index[bugs.pcl$meta$diagnosis=="UC"] > quantile(bugs.pcl$meta$activity_index[bugs.pcl$meta$diagnosis=="nonIBD"], 0.9))

# Fraction of samples at 95%
mean(bugs.pcl$meta$activity_index[bugs.pcl$meta$diagnosis=="CD"] > quantile(bugs.pcl$meta$activity_index[bugs.pcl$meta$diagnosis=="nonIBD"], 0.95))
mean(bugs.pcl$meta$activity_index[bugs.pcl$meta$diagnosis=="UC"] > quantile(bugs.pcl$meta$activity_index[bugs.pcl$meta$diagnosis=="nonIBD"], 0.95))

# Ordinations
species.ord <- pcl.pcoa(bugs.pcl %>% pcl.only(rank="s"))
pdf("./disease_activity/active_disease_pcoa.pdf", 1.15, 1.15)
pcl.ordplot(bugs.pcl, species.ord, colour="active", size_abs=1.4, sortby="activity_index",
            colour_title="", colour_override=c("TRUE"="red", "FALSE"="dodgerblue"),
            colour_names=c("FALSE"="Eubiotic", "TRUE"="Dysbiotic"), outline_size=0.4) +
    theme_nature() + guides(fill="none") +
    theme(axis.text.x=element_blank(), axis.text.y=element_blank())
pcl.ordplot(bugs.pcl, species.ord, colour="active", size_abs=0.7, sortby="activity_index",
            colour_title="", colour_override=c("TRUE"="red", "FALSE"="dodgerblue"),
            colour_names=c("FALSE"="Eubiotic", "TRUE"="Dysbiotic"), outline_size=0.4) +
    theme_nature()
pcl.ordplot(bugs.pcl, species.ord, colour="activity_index", size_abs=1.4, sortby="activity_index",
            colour_title="Dysbiosis Score", outline_size=0.4) +
    theme_nature() + guides(fill="none") +
    theme(axis.text.x=element_blank(), axis.text.y=element_blank())
pcl.ordplot(bugs.pcl, species.ord, colour="activity_index", size_abs=0.7, sortby="activity_index",
            colour_title="Dysbiosis Score", outline_size=0.4) +
    theme_nature()
dev.off()
pdf("./disease_activity/active_disease_pcoa_vanilla.pdf", 1.6, 1.6)
pcl.ordplot(bugs.pcl, species.ord, colour="diagnosis", size_abs=1.4, sortby="activity_index",
            colour_title="Diagnosis", outline_size=0.4, colour_override=hmp2_disease_colors) +
    theme_nature() + guides(fill="none") +
    theme(axis.text.x=element_blank(), axis.text.y=element_blank())
dev.off()

# Marginal disease distribution (ED figure)
df <- cbind(pc1=species.ord$points[,1], pc2=species.ord$points[,2], bugs.pcl$meta)
ggp <- ggplot(data=df, aes(x=pc2, fill=diagnosis, color=diagnosis)) +
    geom_density(alpha=0.3, lwd=0.25) +
    scale_color_manual(values=hmp2_disease_colors, name="Diagnosis") +
    scale_fill_manual(values=hmp2_disease_colors, name="Diagnosis") +
    xlab("PCo 2") + ylab("Density") + theme_nature() +
    guides(color="none", fill="none")
pdf("./disease_activity/pco2_disease_enrichment.pdf", 1.6, 1)
print(ggp)
dev.off()


# How many dysbiotic samples total?
sum(bugs.pcl$meta$active)
# [1] 272
mean(bugs.pcl$meta$active)
# [1] 0.1705329
sum(bugs.pcl$meta$active[bugs.pcl$meta$diagnosis=="CD"])
# [1] 178
mean(bugs.pcl$meta$active[bugs.pcl$meta$diagnosis=="CD"])
# [1] 0.2431694
sum(bugs.pcl$meta$active[bugs.pcl$meta$diagnosis=="UC"])
# [1] 51
mean(bugs.pcl$meta$active[bugs.pcl$meta$diagnosis=="UC"])
# [1] 0.1167048


# Dump active/inactive
source("./common/load_ecs.r")
source("./common/load_proteins.r")
source("./common/load_metabolites.r")
write.table(data.frame(bugs.pcl$meta$active, row.names=rownames(bugs.pcl$meta)),
            file="./disease_activity/active_samples_mgx.tsv",
            sep="\t", quote=F, row.names=T, col.names=F)
metabolites.named.pcl.mt <- merge_disease_activity(metabolites.named.pcl, lenience=0)
write.table(data.frame(metabolites.named.pcl.mt$meta$active, row.names=rownames(metabolites.named.pcl.mt$meta)),
            file="./disease_activity/active_samples_mbx.tsv",
            sep="\t", quote=F, row.names=T, col.names=F)
proteins.pcl.mt <- merge_disease_activity(proteins.pcl, lenience=0)
write.table(data.frame(proteins.pcl.mt$meta$active, row.names=rownames(proteins.pcl.mt$meta)),
            file="./disease_activity/active_samples_mpx.tsv",
            sep="\t", quote=F, row.names=T, col.names=F)
ec.rna.unstrat.pcl.mt <- merge_disease_activity(ec.rna.unstrat.pcl, lenience=0)
write.table(data.frame(ec.rna.unstrat.pcl.mt$meta$active, row.names=rownames(ec.rna.unstrat.pcl.mt$meta)),
            file="./disease_activity/active_samples_mtx.tsv",
            sep="\t", quote=F, row.names=T, col.names=F)



# Dysbiosis stats
all_dysbiotic <- c()
dysbiosis_dt <- data.frame()
for (sub in levels(bugs.pcl$meta$subject)) {
    bugs.sub <- bugs.pcl %>% pcl.filter.s(subject==sub) %>% pcl.sort.s(collection)

    act_change <- diff(bugs.sub$meta$active)
    change_wk <- bugs.sub$meta$week_num[act_change!=0]
    change_delta <- diff(change_wk)

    dysbiosis_dt <- rbind(dysbiosis_dt, data.frame(
        active = (seq_along(change_delta) + bugs.sub$meta$active[1]) %% 2 == 1,
        delta = change_delta,
        censored = rep(F, length(change_delta)),
        diagnosis = rep(bugs.sub$meta$diagnosis[1], length(change_delta))
    ))
    if (length(change_wk) >= 1) {
        dysbiosis_dt <- rbind(dysbiosis_dt, data.frame(
            active = bugs.sub$meta$active[bugs.sub$ns],
            delta = bugs.sub$meta$week_num[bugs.sub$ns] - change_wk[length(change_wk)],
            censored = T,
            diagnosis = bugs.sub$meta$diagnosis[1]
        ))
    }

    if (all(bugs.sub$meta$active)) {
        all_dysbiotic <- c(all_dysbiotic, sub)
    }
}

# 0-length censored intervals contain no data
dysbiosis_dt <- dysbiosis_dt[!dysbiosis_dt$censored | dysbiosis_dt$delta>0,]


# Mean duration of dysbiotic and non-dysbiotic states
mean_dysbiosis_duration <- function(diagnosis, act) {
    dt_diag <- dysbiosis_dt$diagnosis == diagnosis
    dt_act <- dysbiosis_dt$active == act
    return (sum(dysbiosis_dt$delta[dt_diag & dt_act]) / sum(dt_diag & dt_act & !dysbiosis_dt$censored))
}
dys_uc_wks <- mean_dysbiosis_duration("UC", T)
dys_uc_wks
# [1] 4.076923
ndys_uc_wks<- mean_dysbiosis_duration("UC", F)
ndys_uc_wks
# [1] 17.16
dys_cd_wks <- mean_dysbiosis_duration("CD", T)
dys_cd_wks
# [1] 7.815789
ndys_cd_wks<- mean_dysbiosis_duration("CD", F)
ndys_cd_wks
# [1] 12.77419


km_estimator <- function(xx, delta, censored) {
    return (cumprod(sapply(xx, function(x)
        1 - sum((delta == x) & !censored) / sum((delta >= x) & ((delta > x) | !censored))
    )))
}

# Kaplan-Meier plot
xx <- seq(from=0, to=max(dysbiosis_dt$delta[(dysbiosis_dt$diagnosis != "nonIBD")&dysbiosis_dt$active]) + 1)
mask <- (dysbiosis_dt$diagnosis == "UC") & dysbiosis_dt$active
yy_uc <- km_estimator(xx, dysbiosis_dt$delta[mask], dysbiosis_dt$censored[mask])
mask <- (dysbiosis_dt$diagnosis == "CD") & dysbiosis_dt$active
yy_cd <- km_estimator(xx, dysbiosis_dt$delta[mask], dysbiosis_dt$censored[mask])
df_surv <- data.frame(x=c(xx, xx), y=c(yy_uc, yy_cd),
                      diagnosis=c(rep("UC", length(xx)), diagnosis=rep("CD", length(xx))))
df_surv$y[df_surv$y==0] <- NA

# Exponential fit
exp_fit_dys_x <- seq(from=0, to=max(xx), length=201)
exp_fit_dys_y_uc <- exp(-exp_fit_dys_x / dys_uc_wks)
exp_fit_dys_y_cd <- exp(-exp_fit_dys_x / dys_cd_wks)
df_fit <- data.frame(x=c(exp_fit_dys_x, exp_fit_dys_x), y=c(exp_fit_dys_y_uc, exp_fit_dys_y_cd),
                     diagnosis=c(rep("UC", length(exp_fit_dys_x)), diagnosis=rep("CD", length(exp_fit_dys_x))))

# Make the plot
ggp_dys <- ggplot() +
    geom_line(data=df_fit, aes(x=x, y=y, color=diagnosis, group=diagnosis), size=0.3, linetype=2, color="black") +
    geom_step(data=df_surv, aes(x=x, y=y, color=diagnosis, group=diagnosis)) +
    scale_color_manual(values=hmp2_disease_colors) +
    xlab("Duration (weeks)") + ylab("Survival fraction") +
    scale_y_continuous(expand=c(0, 0), limits=c(0, 1), breaks=c(0, 1)) + theme_nature()


# Kaplan-Meier plot
xx <- seq(from=0, to=max(dysbiosis_dt$delta[(dysbiosis_dt$diagnosis != "nonIBD")&!dysbiosis_dt$active]) + 1)
mask <- (dysbiosis_dt$diagnosis == "UC") & !dysbiosis_dt$active
yy_uc <- km_estimator(xx, dysbiosis_dt$delta[mask], dysbiosis_dt$censored[mask])
mask <- (dysbiosis_dt$diagnosis == "CD") & !dysbiosis_dt$active
yy_cd <- km_estimator(xx, dysbiosis_dt$delta[mask], dysbiosis_dt$censored[mask])
df_surv <- data.frame(x=c(xx, xx), y=c(yy_uc, yy_cd),
                      diagnosis=c(rep("UC", length(xx)), diagnosis=rep("CD", length(xx))))
df_surv$y[df_surv$y==0] <- NA

# Exponential fit
exp_fit_ndys_x <- seq(from=0, to=max(xx), length=201)
exp_fit_ndys_y_uc <- exp(-exp_fit_ndys_x / ndys_uc_wks)
exp_fit_ndys_y_cd <- exp(-exp_fit_ndys_x / ndys_cd_wks)
df_fit <- data.frame(x=c(exp_fit_ndys_x, exp_fit_ndys_x), y=c(exp_fit_ndys_y_uc, exp_fit_ndys_y_cd),
                     diagnosis=c(rep("UC", length(exp_fit_ndys_x)), diagnosis=rep("CD", length(exp_fit_ndys_x))))

# Make the plot
ggp_ndys <- ggplot() +
    geom_line(data=df_fit, aes(x=x, y=y, color=diagnosis, group=diagnosis), size=0.3, linetype=2, color="black") +
    geom_step(data=df_surv, aes(x=x, y=y, color=diagnosis, group=diagnosis)) +
    scale_color_manual(values=hmp2_disease_colors) +
    xlab("Interval (weeks)") + ylab(NULL) +
    scale_y_continuous(expand=c(0, 0), limits=c(0, 1), breaks=c(0, 1)) + theme_nature()


# Disease time distribution plots (duration and intervals)
pdf("./disease_activity/dysbiosis_time_distributions.pdf", 5, 4)
print(ggp_dys + guides(color="none"))
print(ggp_dys)
print(ggp_ndys + guides(color="none"))
print(ggp_ndys)
dev.off()

pdf("./disease_activity/dysbiosis_time_distributions_fig2.pdf", 1.6, .81, onefile=F)
library(egg)
ggarrange(ggp_dys + guides(color="none"), ggp_ndys + guides(color="none") + theme(axis.text.y = element_blank()), nrow=1)
dev.off()


# Correlation to other measures?
library(ggplot2)
library(cowplot)
pdf("./disease_activity/fecalcal_and_active_disease.pdf", 6, 4)
ggplot(data=bugs.pcl$meta, aes(x=diagnosis, y=fecalcal, fill=active)) +
    geom_boxplot() +
    geom_boxplot_n(10) +
    theme_cowplot()
dev.off()

source("./common/disease_colors.r")
pdf("./disease_activity/fecalcal_vs_activity.pdf", 6, 5)
ggplot(data=bugs.pcl$meta, aes(y=fecalcal, x=activity_index, color=diagnosis, fill=diagnosis)) +
    geom_point() +
    geom_smooth(method="lm") +
    ylab("Calprotectin (ug/g)") + xlab("Disease Activity") +
    theme_cowplot() + scale_fill_manual(values=hmp2_disease_colors) + scale_color_manual(values=hmp2_disease_colors)
dev.off()

pdf("./disease_activity/sccai_vs_activity.pdf", 5.2, 5)
ggplot(data=pcl.filter.s(bugs.pcl, diagnosis=="UC")$meta, aes(y=sccai, x=activity_index)) +
    geom_point() +
    geom_smooth(method="lm") +
    ylab("SCCAI") + xlab("Disease Activity") +
    theme_cowplot()
dev.off()
sccai_activity_cor <- with(pcl.filter.s(bugs.pcl, diagnosis=="UC")$meta, cor.test(activity_index, sccai))

pdf("./disease_activity/hbi_vs_activity.pdf", 5.2, 5)
ggplot(data=pcl.filter.s(bugs.pcl, diagnosis=="CD")$meta, aes(y=hbi, x=activity_index)) +
    geom_point() +
    geom_smooth(method="lm") +
    ylab("HBI") + xlab("Disease Activity") +
    theme_cowplot()
dev.off()
hbi_activity_cor <- with(pcl.filter.s(bugs.pcl, diagnosis=="CD")$meta, cor.test(activity_index, hbi))




# Define pre-active samples
source("./common/load_metabolites.r")
source("./common/load_pathways.r")
source("./common/load_proteins.r")
source("./common/match_datasets.r")
source("./common/disease_activity.r")

mgx_active_samples <- rownames(bugs.pcl$meta)[bugs.pcl$meta$active]
mbx_active_samples_lst <- match_datasets(list(bugs.pcl, metabolites.pcl))
mbx_active_samples <- rownames(mbx_active_samples_lst[[2]]$meta)[mbx_active_samples_lst[[1]]$meta$active]
rm(mbx_active_samples_lst)
mtx_active_samples_lst <- match_datasets(list(bugs.pcl, pwy.rna.unstrat.pcl))
mtx_active_samples <- rownames(mtx_active_samples_lst[[2]]$meta)[mtx_active_samples_lst[[1]]$meta$active]
rm(mtx_active_samples_lst)
mpx_active_samples_lst <- match_datasets(list(bugs.pcl, proteins.kos.pcl))
mpx_active_samples <- rownames(mpx_active_samples_lst[[2]]$meta)[mpx_active_samples_lst[[1]]$meta$active]
rm(mpx_active_samples_lst)

mgx_preactive_samples <- c()
mtx_preactive_samples <- c()
mbx_preactive_samples <- c()
mpx_preactive_samples <- c()
bugs.pcl.srt <- bugs.pcl %>% pcl.sort.s(collection)
for (subject in levels(bugs.pcl.srt$meta$subject)) {
    # Find collections where there were "activations"
    smt <- bugs.pcl.srt$meta[bugs.pcl.srt$meta$subject==subject,]
    activations <- smt$collection[which(diff(smt$active)>0)+1]

    mmt <- metabolites.pcl$meta[metabolites.pcl$meta$subject==subject,]
    tmt <- pwy.rna.unstrat.pcl$meta[pwy.rna.unstrat.pcl$meta$subject==subject,]
    pmt <- proteins.kos.pcl$meta[proteins.kos.pcl$meta$subject==subject,]
    for (actCol in activations) {
        # Find samples immediately before the collections where there was an
        # activation for each datatype
        cand <- which(smt$collection < actCol)
        if (length(cand)>0) {
            cand <- cand[smt$collection[cand] == max(smt$collection[cand])]
            mgx_preactive_samples <- c(mgx_preactive_samples, rownames(smt)[cand])
        }

        cand <- which(mmt$collection < actCol)
        if (length(cand)>0) {
            cand <- cand[mmt$collection[cand] == max(mmt$collection[cand])]
            mbx_preactive_samples <- c(mbx_preactive_samples, rownames(mmt)[cand])
        }

        cand <- which(tmt$collection < actCol)
        if (length(cand)>0) {
            cand <- cand[tmt$collection[cand] == max(tmt$collection[cand])]
            mtx_preactive_samples <- c(mtx_preactive_samples, rownames(tmt)[cand])
        }

        cand <- which(tmt$collection < actCol)
        if (length(cand)>0) {
            cand <- cand[pmt$collection[cand] == max(pmt$collection[cand])]
            mpx_preactive_samples <- c(mtx_preactive_samples, rownames(pmt)[cand])
        }
    }
}


length(setdiff(mgx_preactive_samples, mgx_active_samples))
# [1] 88

mgx_sample_activity <- rep("Inactive", bugs.pcl.srt$ns)
names(mgx_sample_activity) <- rownames((bugs.pcl.srt %>% pcl.sort.s(subject))$meta)
mgx_sample_activity[names(mgx_sample_activity) %in% mgx_preactive_samples] <- "PreActive"
mgx_sample_activity[names(mgx_sample_activity) %in% mgx_active_samples] <- "Active"
write.table(mgx_sample_activity,
            file="./disease_activity/preactive_samples_mgx.tsv",
            sep="\t", quote=F, row.names=T, col.names=F)


length(setdiff(mbx_preactive_samples, mbx_active_samples))
# [1] 56

mbx_sample_activity <- rep("Inactive", metabolites.pcl$ns)
names(mbx_sample_activity) <- rownames(metabolites.pcl$meta)
mbx_sample_activity[names(mbx_sample_activity) %in% mbx_preactive_samples] <- "PreActive"
mbx_sample_activity[names(mbx_sample_activity) %in% mbx_active_samples] <- "Active"
write.table(mbx_sample_activity,
            file="./disease_activity/preactive_samples_mbx.tsv",
            sep="\t", quote=F, row.names=T, col.names=F)


length(setdiff(mtx_preactive_samples, mtx_active_samples))
# [1] 64

mtx_sample_activity <- rep("Inactive", pwy.rna.unstrat.pcl$ns)
names(mtx_sample_activity) <- rownames(pwy.rna.unstrat.pcl$meta)
mtx_sample_activity[names(mtx_sample_activity) %in% mtx_preactive_samples] <- "PreActive"
mtx_sample_activity[names(mtx_sample_activity) %in% mtx_active_samples] <- "Active"
write.table(mtx_sample_activity,
            file="./disease_activity/preactive_samples_mtx.tsv",
            sep="\t", quote=F, row.names=T, col.names=F)

length(setdiff(mpx_preactive_samples, mpx_active_samples))
# [1] 70

mpx_sample_activity <- rep("Inactive", proteins.kos.pcl$ns)
names(mpx_sample_activity) <- rownames(proteins.kos.pcl$meta)
mpx_sample_activity[names(mpx_sample_activity) %in% mpx_preactive_samples] <- "PreActive"
mpx_sample_activity[names(mpx_sample_activity) %in% mpx_active_samples] <- "Active"
write.table(mpx_sample_activity,
            file="./disease_activity/preactive_samples_mpx.tsv",
            sep="\t", quote=F, row.names=T, col.names=F)


## Other definitions of activity
library(vegan)
D <- as.matrix(vegdist((bugs.pcl %>%
    #pcl.filter.f(!grepl("k__Viruses", Name)) %>%
    pcl.only(rank="s") %>% pcl.normalize)$x, method="bray"))

# The "reference set" of inactive samples
ref_set <- (bugs.pcl$meta$diagnosis == "nonIBD") & (bugs.pcl$meta$week_num >= 20)

# Calculate the activity index
quant <- 0.1
activity_index <- sapply(seq_len(length(ref_set)), function(i)
    quantile(D[i, ref_set & (bugs.pcl$meta$subject != bugs.pcl$meta$subject[i])], probs=quant))

ggplot(data=data.frame(diagnosis=bugs.pcl$meta$diagnosis, activity_index=activity_index), aes(x=activity_index, color=diagnosis, fill=diagnosis)) +
    geom_density(alpha=0.1) +
    geom_vline(xintercept=quantile(activity_index[bugs.pcl$meta$diagnosis=="nonIBD"], 0.95), size=1.5) +
    scale_color_manual(values=hmp2_disease_colors, name="Diagnosis") +
    scale_fill_manual(values=hmp2_disease_colors, name="Diagnosis") +
    xlab("Bray-Curtis dissimilarity to non-self nonIBD samples") +
    ggtitle(sprintf("Quantile = %.2f", quant)) +
    ylab(NULL)


mean(activity_index[bugs.pcl$meta$diagnosis=="CD"] > quantile(activity_index[bugs.pcl$meta$diagnosis=="nonIBD"], 0.95))
mean(activity_index[bugs.pcl$meta$diagnosis=="UC"] > quantile(activity_index[bugs.pcl$meta$diagnosis=="nonIBD"], 0.95))


# Really eubiotic ordination (<10%)
species.ord <- pcl.pcoa(bugs.pcl %>% pcl.only(rank="s") %>% pcl.normalize)
bugs.pcl$meta$active3 <- factor("Non-dysbiotic", levels=c("Dysbiotic", "Non-dysbiotic", "Eubiotic"))
bugs.pcl$meta$active3[bugs.pcl$meta$active] <- "Dysbiotic"
bugs.pcl$meta$active3[bugs.pcl$meta$activity_index < eubiosis_lower_threshold] <- "Eubiotic"
bugs.pcl$meta$active3_sort <- abs(bugs.pcl$meta$activity_index - (eubiosis_lower_threshold + disease_activity_threshold)/2)
pdf("./disease_activity/active3_disease_pcoa.pdf", 2.15, 2.15)
pcl.ordplot(bugs.pcl, species.ord, colour="active3", size_abs=1.0,, sortby = "active3_sort",
            colour_title="", colour_override=c("Dysbiotic"="red", "Non-dysbiotic"="grey70", "Eubiotic"="dodgerblue"), outline_size=0.3) +
    theme_nature() + guides(fill="none") +
    theme(axis.text.x=element_blank(), axis.text.y=element_blank())
dev.off()


# Longitudinal serology
source("./common/load_serology.r")
serology.eu.pcl <- merge_disease_activity(serology.eu.pcl, lenience=4)
serology.eu.pcl.cd <- pcl.filter.s(serology.eu.pcl, diagnosis == "CD")
df.cd <- cbind(serology.eu.pcl.cd$meta, serology.eu.pcl.cd$x)
serology.eu.pcl.uc <- pcl.filter.s(serology.eu.pcl, diagnosis == "UC")
df.uc <- cbind(serology.eu.pcl.cd$meta, serology.eu.pcl.cd$x)

library(ggplot2)
plot_active_titers <- function(df, var) {
    ggp <- ggplot(data=df, aes_string(x="week_num", y=sprintf("`%s`", var))) +
        geom_line() +
        geom_point(aes(color=active)) +
        facet_wrap(~subject) +
        scale_color_manual(values=c("TRUE"="red", "FALSE"="black"))
    return (ggp)
}

pdf("./disease_activity/serology_vs_activity.pdf", 10, 10)
for (v in colnames(serology.eu.pcl$x)) {
    print(plot_active_titers(df.uc, v) +
              ggtitle(sprintf("%s in %s", v, "UC")))
    print(plot_active_titers(df.cd, v) +
              ggtitle(sprintf("%s in %s", v, "CD")))
}
dev.off()


