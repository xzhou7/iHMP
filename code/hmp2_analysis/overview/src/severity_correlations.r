
source("./common/merge_metadata.r")
source("./common/disease_activity.r")
source("./common/disease_colors.r")
source("./common/theme_nature.r")

library(ggplot2)
library(GGally)
library(cowplot)
library(egg)

mt <- hmp2_sample_metadata
act_mtch <- match(mt$site_sub_coll, gsub("_P?T?R?$", "", bugs.pcl$meta$site_sub_coll))
mt$activity_index <- NA
mt$activity_index[!is.na(act_mtch)] <- bugs.pcl$meta$activity_index[na.omit(act_mtch)]
mt$frac_human_reads_mtx <- NA
mtt <- mt[mt$data_type == "metatranscriptomics",]
act_mtch <- match(mt$site_sub_coll, mtt$site_sub_coll)
mt$frac_human_reads_mtx[!is.na(act_mtch)] <- mtt$frac_human_reads[na.omit(act_mtch)]

mt <- mt[mt$data_type == "metagenomics",]
mt <- mt[!is.na(mt$fecalcal) | !is.na(mt$hbi) | !is.na(mt$sccai),]
mt <- mt[!duplicated(mt$site_sub_coll),]
mts <- mt[,c("hbi", "sccai", "fecalcal", "frac_human_reads", "frac_human_reads_mtx", "activity_index")]
colnames(mts) <- c("HBI", "SCCAI", "FecalCal", "Human Reads MGX", "Human Reads MTX", "Dysbiosis Score")

lowerFn <- function(data, mapping, ...) {
    x <- eval(mapping$x, data)
    y <- eval(mapping$y, data)

    ct <- cor.test(x, y, method="spearman")

    p <- ggally_text(sprintf("rho = %.2f\nP = %.3g", ct$estimate, ct$p.value)) +
        theme_cowplot()
    p
}

mt_sub <- mt
upperFn <- function(data, mapping, ...) {
    p <- ggplot(data = cbind(data, mt_sub[,"diagnosis", drop=F]), mapping = mapping) +
        geom_point(aes(fill=diagnosis), colour = "black", shape=21, size=1.5) +
        scale_fill_manual(values=hmp2_disease_colors) +
        geom_smooth(...) +
        theme_cowplot()
    p
}

diagFn <- function(data, mapping, ...) {
    dt <- eval(mapping$x, data)
    bins <- 30
    if (length(unique(dt)) < bins) {
        bins <- max(dt, na.rm=T) - min(dt, na.rm=T) + 1
    }
    p <- ggplot(data = data, mapping = mapping) +
        stat_bin(bins = bins, ...) +
        theme_cowplot()
    p
}

pdf("./overview/correlations_severity.pdf", 11, 11)
ggpairs(mts,
    upper = list(continuous = wrap(upperFn, method = "lm", color = "red")),
    diag = list(continuous = wrap(diagFn, colour = "blue", fill="slateblue1")),
    lower = list(continuous = wrap(lowerFn))
)
dev.off()

# Colonic CD correlations?
pdf("./overview/correlations_severity_colonic_cd.pdf", 11, 11)
mt_sub <- mt[grepl("L1", mt$baseline_montreal_location),]
ggpairs(mts[grepl("L1", mt$baseline_montreal_location),-2],
        lower = list(continuous = wrap(lowerFn, method = "lm", color = "red")),
        diag = list(continuous = wrap(diagFn, colour = "blue", fill="slateblue1")),
        upper = list(continuous = wrap(upperFn))
)
dev.off()


# Plots for Fig 2
doplot <- function(data, mapping) {
    mapping <-c(mapping, aes(color=diagnosis, fill=diagnosis))
    class(mapping) <- "uneval"
    ggplot(data = data, mapping) +
        geom_smooth(method="lm", size=0.25) +
        geom_point(size=1.1, stroke=0) +
        theme_nature() + scale_fill_manual(values=hmp2_disease_colors) +
        guides(fill="none", color="none") +
        scale_color_manual(values=hmp2_disease_colors)
}
pdf("./overview/fig2_severityassocs.pdf", 1.533, 0.815, onefile=F)
print(ggarrange(
    doplot(mt[mt$diagnosis=="CD",], aes(x=fecalcal, y=hbi)) +
        xlab("Calprotectin (ug/g)") + ylab("HBI") +
        scale_x_continuous(limits=c(0, 400), breaks=c(0, 100, 200, 300, 400)) +
        scale_y_continuous(limits=c(0, 10), breaks=c(0, 2, 4, 6, 8, 10)),
    doplot(mt[mt$diagnosis=="UC",], aes(x=fecalcal, y=sccai)) +
        xlab("Calprotectin (ug/g)") + ylab("SCCAI") +
        scale_x_continuous(limits=c(0, 500), breaks=c(0, 200, 400)) +
        scale_y_continuous(limits=c(0, 8), breaks=c(0, 2, 4, 6, 8)),
    nrow=1
))
dev.off()



# Time-delayed correlations

mtc <- c("hbi", "sccai", "fecalcal", "frac_human_reads", "frac_human_reads_mtx", "activity_index")
names(mtc) <- c("HBI", "SCCAI", "FecalCal", "Human Reads MGX", "Human Reads MTX", "Dysbiosis Score")

samesubject <- outer(mt$Participant.ID, mt$Participant.ID, FUN="==")
dw <- outer(mt$week_num, mt$week_num, FUN="-")
samples <- which((dw == 2) & (samesubject), arr.ind=T)

df <- data.frame(v1=c(), v2=c(), x=c(), y=c(), diagnosis=c())
for (i in seq_along(mtc)) {
    for (j in seq_along(mtc)) {
        df <- rbind(df, data.frame(
            v1 = rep(names(mtc)[i], nrow(samples)),
            v2 = rep(names(mtc)[j], nrow(samples)),
            x = mt[samples[,1], mtc[i]],
            y = mt[samples[,2], mtc[j]],
            diagnosis = mt$diagnosis[samples[,1]]
        ))
    }
}

df$v1 <- factor(df$v1, levels=names(mtc))
df$v2 <- factor(df$v2, levels=names(mtc))

ggplot(data=df, aes(x=x, y=y)) +
    geom_point(aes(fill=diagnosis), shape=21) +
    geom_smooth(method="lm", color="red") +
    facet_grid(v2 ~ v1, scales="free") +
    scale_fill_manual(values=hmp2_disease_colors)


library(RVAideMemoire)

lags <- seq(-8, 8, by=2)
samesubject <- outer(mt$Participant.ID, mt$Participant.ID, FUN="==")
dw <- outer(mt$week_num, mt$week_num, FUN="-")
df <- data.frame(v1=c(), v2=c(), lag=c(), rho=c(), rho_lci=c(), rho_uci=c(), diagnosis=c())
for (li in seq_along(lags)) {
    if (lags[li] == 0) {
        samples <- cbind(seq_len(nrow(mt)), seq_len(nrow(mt)))
    } else {
        samples <- which((dw == lags[li]) & (samesubject), arr.ind=T)
    }

    for (i in seq_along(mtc)) {
        for (j in seq_along(mtc)) {
            keep <- !is.na(mt[samples[,1], mtc[i]]) & !is.na(mt[samples[,2], mtc[j]])
            if (sum(keep) > 10) {
                rho <- spearman.ci(mt[samples[keep,1], mtc[i]], mt[samples[keep,2], mtc[j]])
                df <- rbind(df, data.frame(
                    v1 = names(mtc)[i],
                    v2 = names(mtc)[j],
                    lag = lags[li],
                    rho = rho$estimate,
                    rho_lci = rho$conf.int[1],
                    rho_uci = rho$conf.int[2],
                    diagnosis = "Overall"
                ))
            }
            dkeep <- mt$diagnosis == "UC"
            keep <- dkeep[samples[,1]] & dkeep[samples[,2]] & !is.na(mt[samples[,1], mtc[i]]) & !is.na(mt[samples[,2], mtc[j]])
            if (sum(keep) > 10) {
                rho <- spearman.ci(mt[samples[keep,1], mtc[i]], mt[samples[keep,2], mtc[j]])
                df <- rbind(df, data.frame(
                    v1 = names(mtc)[i],
                    v2 = names(mtc)[j],
                    lag = lags[li],
                    rho = rho$estimate,
                    rho_lci = rho$conf.int[1],
                    rho_uci = rho$conf.int[2],
                    diagnosis = "UC"
                ))
            }
            dkeep <- mt$diagnosis == "CD"
            keep <- dkeep[samples[,1]] & dkeep[samples[,2]] & !is.na(mt[samples[,1], mtc[i]]) & !is.na(mt[samples[,2], mtc[j]])
            if (sum(keep) > 10) {
                rho <- spearman.ci(mt[samples[keep,1], mtc[i]], mt[samples[keep,2], mtc[j]])
                df <- rbind(df, data.frame(
                    v1 = names(mtc)[i],
                    v2 = names(mtc)[j],
                    lag = lags[li],
                    rho = rho$estimate,
                    rho_lci = rho$conf.int[1],
                    rho_uci = rho$conf.int[2],
                    diagnosis = "CD"
                ))
            }
        }
    }
}

df$v1 <- factor(df$v1, levels=names(mtc))
df$v2 <- factor(df$v2, levels=names(mtc))


pdf("./disease_severity/lagged_correlations_95ci.pdf", 11, 11)
ggplot(data=df, aes(x=lag, y=rho, color=diagnosis)) +
    geom_vline(xintercept = 0, color="grey") +
    geom_hline(yintercept = 0) +
    geom_line() +
    geom_point() +
    geom_errorbar(aes(ymin=rho_lci, ymax=rho_uci)) +
    facet_grid(v2 ~ v1)
dev.off()

