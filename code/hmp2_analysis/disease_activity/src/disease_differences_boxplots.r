
library(plyr)
library(dplyr)
library(reshape2)
library(ggplot2)
library(ggridges)
library(cowplot)

source("./env_config.r")
source("./common/disease_colors.r")
source("./common/disease_activity.r")
source("./common/load_bugs.r")
source("./common/load_metabolites.r")
source("./common/theme_nature.r")

zeroed_logridgelineplot <- function(x, w, isz, df, extents, zscale, groups) {
    xs <- split(x, groups)
    ws <- split(w, groups)
    dfs <- split(df, groups)
    iszs <- split(isz, groups)

    dfnz <- data.frame()
    dfz <- data.frame()
    for (i in seq_along(xs)) {
        xi <- xs[[i]]
        wi <- ws[[i]]

        nz <- !iszs[[i]]
        if (sum(nz) > 1) {
            if (is.null(extents)) {
                dn <- density(xi[nz], weights=wi[nz] / sum(wi[nz]))
            } else {
                dn <- density(xi[nz], weights=wi[nz] / sum(wi[nz]), from=extents[1], to=extents[2])
            }
        } else if(is.null(extents)) {
            dn <- list(x=extents, y=c(0, 0))
        } else {
            dn <- list(x=c(0, 0), y=c(0, 0))
        }

        dfnzi <- dfs[[i]][rep(1, length(dn$x)),,drop=F]
        dfnzi$x <- dn$x
        dfnzi$height <- dn$y * sum(wi[nz]) / sum(wi)
        dfnz <- rbind(dfnz, dfnzi)

        dfzi <- dfs[[i]][c(1,1),,drop=F]
        dfzi$ymin <- dfzi$y[1]
        dfzi$ymax <- dfzi$ymin
        dfzi$ymax[2] <- dfzi$ymin[2] + sum(wi[!nz]) * zscale / sum(wi)
        dfzi$xrhs <- c(0, 1)

        dfz <- rbind(dfz, dfzi)
    }

    dfz <- dfz[order(dfz$ymin, decreasing=T),, drop=F]

    return (list(nz=dfnz, z=dfz))
}

make_boxplots <- function(pcl, to_show, active=F, rel=T, logspace=F, pwrspace=1, facet=F, italic=F, alpha=1,
                          reorder=T, limits=NULL, zerobar_width=0.035, sig=NULL, sig_name="FDR", sig_stars=F, label_right=T) {
    if (!("active" %in% colnames(pcl$meta))) {
        pcl <- merge_disease_activity(pcl, lenience=2)
    }

    x <- (pcl %>% pcl.reorder.f(match(to_show, colnames(pcl$x))))$x
    if (!is.null(sig)) {
        sig <- sig[match(to_show, rownames(sig)),,drop=F]
    }
    if (!is.null(names(to_show))) {
        colnames(x) <- names(to_show)[to_show %in% colnames(x)]
        if (!is.null(sig)) {
            rownames(sig) <- names(to_show)
        }
    }
    df <- melt(x)
    colnames(df) <- c("sample", "feature", "value")
    df$diagnosis <- pcl$meta$diagnosis[df$sample]
    if ("subject" %in% colnames(pcl$meta)) {
        df$subject <- pcl$meta$subject[df$sample]
        df$weight <- 1 / ave(rep(1, nrow(df)), df$subject, FUN=sum)
    } else {
        df$weight <- 1
    }

    if (logspace) {
        stopifnot(pwrspace==1)
        df$value <- log10(df$value)
        df$value[is.nan(df$value)] <- -Inf
        if (!is.null(limits)) {
            limits <- log10(limits)
        }
    }
    if (pwrspace!=1) {
        stopifnot(!logspace)
        df$value <- df$value^pwrspace
    }

    if (rel) {
        mask <- df$diagnosis == "nonIBD"
        med <- sapply(split(df$value[mask], df$feature[mask]), function(x)median(x[is.finite(x)]))
        df$value <- if (logspace) {
            df$value - med[df$feature]
        } else {
            df$value / med[df$feature]
        }
        rel_line <- geom_vline(xintercept=0, size=0.4,
                               color=colorRampPalette(c(hmp2_active_disease_colors["nonIBD"], "white"))(21)[5])
    } else {
        rel_line <- list()
    }

    if (facet) {
        faceting <- list(
            facet_wrap(~feature, ncol=1, scales="free_y"),
            theme(strip.background=element_blank(),
                  strip.text=element_text(face=ifelse(italic, "italic", "plain"), margin=c(0,0,0,0)),
                  axis.text.y=element_blank())
        )
    } else {
        faceting <- list(
            theme(
                axis.text.y=element_text(face=ifelse(italic, "italic", "plain"))
            )
        )
    }

    if (active) {
        df$active <- pcl$meta$active
        df$active[is.na(df$active)] <- F
        df$diag_act <- ifelse(df$diagnosis == "nonIBD", "nonIBD",
                              sprintf("%s %s", ifelse(df$active, "Active", "Inactive"), df$diagnosis))
        df$diag_act <- factor(df$diag_act, levels=rev(c("nonIBD", "Inactive UC", "Inactive CD", "Active UC",  "Active CD")))

        # Order by difference in active CD
        if (reorder) {
            mask <- df$diag_act == "Active CD"
            med_cd <- sapply(split(df$value[mask], df$feature[mask]), median)
            df$feature <- factor(df$feature, names(sort(med_cd, T)))
        }

        df$y <- as.numeric(df$feature)-1 + (as.numeric(df$diag_act)-1) * 0.3 / 4

        fillcolors <- hmp2_active_disease_colors_fill
        outlinecolors <- hmp2_active_disease_colors_outline
    } else {
        # Order by difference in CD
        if (reorder) {
            mask <- df$diagnosis == "CD"
            med_cd <- sapply(split(df$value[mask], df$feature[mask]), median)
            df$feature <- factor(df$feature, names(sort(med_cd, T)))
        }

        df$diagnosis <- factor(df$diagnosis, levels=rev(c("nonIBD", "UC", "CD")))
        df$diag_act <- df$diagnosis
        df$y <- as.numeric(df$feature)-1 + (as.numeric(df$diagnosis)-1) * 0.2 / 2

        fillcolors <- hmp2_disease_colors_fill
        outlinecolors <- hmp2_disease_colors_outline
    }

    # split the data into zero and non-zero components
    zll <- zeroed_logridgelineplot(df$value, df$weight, is.infinite(df$value) & df$value<0,
                                   df, extents=limits, 1, factor(df$y))

    # assign x values for the zero barplots
    if (is.null(limits)) {
        limits <- c(min(zll$nz$x), max(zll$nz$x))
    }
    if (any(zll$z$ymax > zll$z$ymin)) {
        zw <- zerobar_width * (limits[2] - limits[1])
        zll$z$xmin <- limits[1] + zw * (as.numeric(zll$z$diag_act) - 1) - zw * (length(levels(zll$z$diag_act)) + 1.35)
        zll$z$xmax <- limits[1] + zll$z$xrhs * (zw * (as.numeric(zll$z$diag_act) + 1) - zw * (length(levels(zll$z$diag_act)) + 1.35))

        zeroplots <- list(
            geom_rect(data=zll$z, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), lwd=0.25, alpha=alpha),
            geom_vline(xintercept=limits[1], color="black", size=0.3)
        )
    } else {
        zw <- 0
        zeroplots <- list()
    }

    # make a dataframe to show labels
    dflbl <- df[df$y == as.numeric(df$feature)-1,,drop=F]
    dflbl <- dflbl[!duplicated(dflbl$y),,drop=F]

    # make q-value labels
    if (!is.null(sig)) {
        dflbl$q <- sig$minQ[match(dflbl$feature, rownames(sig))]
        if (!sig_stars) {
            dflbl$qlbl <- sprintf("  %s %.3g", sig_name, dflbl$q)
        } else {
            dflbl$qlbl <- ifelse(dflbl$q < 0.001, "***",
                          ifelse(dflbl$q < 0.01,  "**",
                          ifelse(dflbl$q < 0.05,  "*", "")))
        }
        face <- if(sig_stars){"bold"}else{"plain"}

        qval_labels <- geom_text(data=dflbl, aes(y=y-0.14, label=qlbl), x=ifelse(label_right, limits[1], limits[2]),
                                 hjust=ifelse(label_right, 0, 1), color="black", size=1.785714, fontface=face)
    } else {
        qval_labels <- list()
    }

    # scale the heights to look nice
    zll$nz$height <- zll$nz$height * 0.9 / max(zll$nz$height)

    ggp <- ggplot(data=zll$nz, aes(group=factor(y), fill=diag_act, color=diag_act)) + rel_line +
        geom_ridgeline(aes(height=height, x=x, y=y), stat="identity", scale=1, lwd=0.25, alpha=alpha) +
        zeroplots +
        geom_text(data=dflbl, aes(y=y-0.14, label=sprintf("  %s ", feature)), x=ifelse(label_right, limits[2], limits[1]),
                  hjust=ifelse(label_right, 1, 0), color="black", fontface=ifelse(italic, "italic", "plain"), size=1.785714) +
        qval_labels +
        scale_fill_manual(values=fillcolors, name="Diagnosis") +
        scale_color_manual(values=outlinecolors, name="Diagnosis") +
        theme_nature() + xlab(NULL) + ylab("Density") + faceting +
        scale_y_continuous(expand=c(0, 0), limits=c(-0.35, 0.3+length(levels(df$feature)))) +
        scale_x_continuous(expand=c(0, 0), limits=c(min(limits[1], suppressWarnings(min(zll$z$xmin)-zw*0.35)), limits[2])) +
        theme(axis.text.y = element_blank())

    return (ggp)
}




# Significance
sig_mbx <- read.delim(header=T, file.path(HMP2_root, "analysis",
                                          "differential_abundance", "Vanilla", "active_vs_inactive", "Active_IBD_relations_metabolites.tsv"))
sig_mbx$ID <- sig_mbx$X.
rownames(sig_mbx) <- sig_mbx$ID

sig_bugs <- read.delim(header=T,
    file.path(HMP2_root, "analysis", "differential_abundance", "Vanilla",
              "active_vs_inactive", "Active_IBD_relations_bugs.tsv"))
sig_bugs$ID <- gsub("_", " ", sig_bugs$X.)
rownames(sig_bugs) <- sig_bugs$ID


strip_widths <- 1.2
ybplate <- 0.091
ysclbplate <- 0.65
entry_ht <- 0.28
full_height <- function(entries, ht) ybplate + (ysclbplate + entries) * entry_ht

# Paper plots
pdf("./disease_activity/act_taxonomy.pdf", strip_widths, full_height(7))
ggp <- make_boxplots(bugs.pcl %>% pcl.nicenames, active=T, logspace=T, italic=T,
                     limits=c(10^-3.5, 1e4),
              c(Escherichia = "Escherichia",
                "F. prausnitzii" = "Faecalibacterium",
                Subdoligranulum = "Subdoligranulum",
                Roseburia = "Roseburia",
                "R. gnavus" = "Ruminococcus gnavus",
                "R. torques" = "Ruminococcus torques",
                Alistipes = "Alistipes")) +
    guides(fill="none", color="none") + ylab(NULL) +
    theme(axis.ticks.y = element_blank())
plot(ggp)
dev.off()

pdf("./disease_activity/act_scfas.pdf", strip_widths, full_height(3))
ggp <- make_boxplots(metabolites.named.pcl, active=T, logspace=T,
                     reorder=F, limits=c(1e-2, 10^1.2), label_right=F,
              c("valerate/isovalerate", "propionate", "butyrate")) +
    guides(fill="none", color="none") + ylab(NULL) +
    theme(axis.ticks.y = element_blank())
plot(ggp)
dev.off()

pdf("./disease_activity/act_othermbx.pdf", strip_widths, full_height(3))
ggp <- make_boxplots(metabolites.named.pcl, active=T, logspace=T,
                     reorder=F, limits=c(1e-5, 1e1), label_right=F,
              c("urobilin", "indole-3-propionate", "adipate")) +
    guides(fill="none", color="none") + ylab(NULL) +
    theme(axis.ticks.y = element_blank())
plot(ggp)
dev.off()

pdf("./disease_activity/act_bileacids.pdf", strip_widths, full_height(7))
ggp <- make_boxplots(metabolites.named.pcl, active=T, logspace=T,
                     reorder=F, limits=c(10^-3.5, 10^3), label_right=F,
              c(lithocholate = "lithocholate",
                deoxycholate = "deoxycholate",
                cholate = "cholate",
                glycochenodeoxycholate = "glycochenodeoxycholate",
                taurochenodeoxycholate = "taurochenodeoxycholate",
                glycocholate = "glycocholate",
                taurocholate = "taurocholate")) +
    guides(fill="none", color="none") + ylab(NULL) +
    theme(axis.ticks.y = element_blank())
plot(ggp)
dev.off()


pdf("./disease_activity/act_xsectional.pdf", strip_widths, full_height(11))
ggp <- make_boxplots(metabolites.named.pcl, active=F, logspace=T,
                     reorder=F, limits=c(1e-2, 10^3),
                     c("nicotinate", "taurine", "putrescine", "porphobilinogen",
                       "arachidonate", "uridine", "hydroxycotinine", "adrenate", "nicotinuric acid",
                       "ethyl glucuronide")) +
    guides(fill="none", color="none") + ylab(NULL) +
    theme(axis.ticks.y = element_blank())
plot(ggp)
dev.off()


# ED figure for Carnitine
entry_ht <- 0.38
is_carnitine <- grepl("[Cc]arnitine", sig_mbx$ID)
include_carnitines <- head(which(is_carnitine), 10)
pdf("./disease_activity/act_carnitine.pdf", 2.25, full_height(10))
ggp <- make_boxplots(metabolites.named.pcl, active=T, logspace=T, reorder=F,
                     limits=c(10^-3, 10^4), zerobar_width=0.025, sig=sig_mbx, sig_name="FDR",
                     rev(sig_mbx$X.[include_carnitines])) +
    guides(fill="none", color="none") + ylab(NULL) +
    theme(axis.ticks.y = element_blank())
plot(ggp)
dev.off()


# ED figure for everything else
already_included <- c(
    as.character(sig_mbx$ID[is_carnitine]),
    "lithocholate", "deoxycholate", "cholate",
    "glycochenodeoxycholate", "taurochenodeoxycholate",
    "glycocholate", "taurocholate",
    "urobilin", "indole-3-propionate", "adipate",
    "valerate/isovalerate", "propionate", "butyrate"
)

pdf("./disease_activity/act_others1.pdf", 2.25, full_height(10))
other_stuff <- head(which(!(sig_mbx$ID %in% already_included)), 10)
ggp <- make_boxplots(metabolites.named.pcl, active=T, logspace=T, reorder=F,
                     limits=c(10^-3, 10^4), zerobar_width=0.025, sig=sig_mbx, sig_name="FDR",
                     rev(sig_mbx$ID[other_stuff])) +
    guides(fill="none", color="none") + ylab(NULL) +
    theme(axis.ticks.y = element_blank())
plot(ggp)
dev.off()

already_included <- sig_bugs$ID %in% c(
    "Escherichia coli",
    "Escherichia unclassified",
    "Faecalibacterium prausnitzii",
    "Subdoligranulum unclassified",
    "Roseburia hominis",
    "Roseburia intestinalis",
    "Roseburia inulinivorans",
    "Alistipes shahii",
    "Alistipes onderdonkii",
    "Alistipes putredinis",
    "Ruminococcus gnavus",
    "Ruminococcus torques",
    "Ruminococcus obeum"
)
low_prev <- sig_bugs$prevalence <= 0.5
keep <- !already_included & !low_prev

pdf("./disease_activity/act_bugs2.pdf", 2.25, full_height(7))
other_stuff <- head(which(keep), 7)
ggp <- make_boxplots(bugs.pcl %>% pcl.nicenames, active=T, logspace=T, reorder=F, italic=T,
                     limits=c(10^-3, 10^4), zerobar_width=0.025, sig=sig_bugs, sig_name="FDR",
                     rev(sig_bugs$ID[other_stuff])) +
    guides(fill="none", color="none") + ylab(NULL) +
    theme(axis.ticks.y = element_blank())
plot(ggp)
dev.off()

source("./common/load_kos.r")
source("./common/match_datasets.r")

pdf("./disease_activity/act_tax_mtx.pdf", strip_widths, full_height(3))
entry_ht <- 0.28
mt <- match_datasets(list(bugs.fromko.dna.pcl, bugs.fromko.rna.pcl))
toshow <- c("C. hathewayi" = "Clostridium hathewayi",
            "C. bolteae" = "Clostridium bolteae",
            "R. gnavus" = "Ruminococcus gnavus")
mt[[1]] <- pcl.reorder.f(mt[[1]] %>% pcl.nicenames, toshow)
mt[[2]] <- pcl.reorder.f(mt[[2]] %>% pcl.nicenames, toshow)
mt.bugs.fromko.rnadna.pcl <- mt[[2]]
mt.bugs.fromko.rnadna.pcl$x <- mt.bugs.fromko.rnadna.pcl$x / mt[[1]]$x
ggp <- make_boxplots(mt.bugs.fromko.rnadna.pcl, active=T, logspace=T, rev(toshow), reorder=F, italic=T,
                     limits=c(10^-1.5, 10^2), label_right=T) +
    guides(fill="none", color="none") + ylab(NULL) +
    theme(axis.ticks.y = element_blank())
plot(ggp)
dev.off()


pdf("./disease_activity/act_others1.pdf", 2.25, full_height(10))
other_stuff <- head(which(!(sig_mbx$ID %in% already_included)), 10)
ggp <- make_boxplots(metabolites.named.pcl, active=T, logspace=T, reorder=F,
                     limits=c(10^-3, 10^4), zerobar_width=0.025, sig=sig_mbx, sig_name="FDR",
                     rev(sig_mbx$ID[other_stuff])) +
    guides(fill="none", color="none") + ylab(NULL) +
    theme(axis.ticks.y = element_blank())
plot(ggp)
dev.off()





dump_all_features <- function(pcl, sigfile, outfile, ...) {
    sig <- read.delim(header=T, sigfile)
    sig$ID <- gsub("_", " ", sig$X.)
    rownames(sig) <- sig$ID
    ggp <- make_boxplots(pcl, rev(sig$ID), logspace=T, reorder=F,
                         limits=c(10^-4, 10^4), zerobar_width=0.025, sig=sig, sig_name="FDR", ...) +
        guides(fill="none", color="none") + ylab(NULL) +
        theme(axis.ticks.y = element_blank())

    pdf(outfile, strip_widths*1.8, full_height(nrow(sig))*1.4)
    plot(ggp)
    dev.off()
}


# Really tall plots of all features by model
for (model in c("CLR", "LOG")) {
    dump_all_features(bugs.pcl %>% pcl.nicenames,
                      file.path(HMP2_root, "analysis", "differential_abundance",
                                model, "diagnosis", "IBD_relations_bugs.tsv"),
                      sprintf("./disease_activity/diagnosis_taxa_%s.pdf", model),
                      italic=T)
    dump_all_features(bugs.pcl %>% pcl.nicenames,
                      file.path(HMP2_root, "analysis", "differential_abundance",
                                model, "active_vs_inactive", "Active_IBD_relations_bugs.tsv"),
                      sprintf("./disease_activity/dysbiosis_taxa_%s.pdf", model),
                      italic=T, active=T)
    dump_all_features(metabolites.named.pcl,
                      file.path(HMP2_root, "analysis", "differential_abundance",
                                model, "diagnosis", "IBD_relations_metabolites.tsv"),
                      sprintf("./disease_activity/diagnosis_mbx_%s.pdf", model))
    dump_all_features(metabolites.named.pcl,
                      file.path(HMP2_root, "analysis", "differential_abundance",
                                model, "active_vs_inactive", "Active_IBD_relations_metabolites.tsv"),
                      sprintf("./disease_activity/dysbiosis_mbx_%s.pdf", model),
                      active=T)
}

sig_bugs_clr <- read.delim(header=T,
    file.path(HMP2_root, "analysis", "differential_abundance",
        "diagnosis", "IBD_relations_bugs_CLR.tsv"))
sig_bugs_clr$ID <- gsub("_", " ", sig_bugs_clr$X.)
rownames(sig_bugs_clr) <- sig_bugs_clr$ID
pdf("./disease_activity/act_tax_clr.pdf", strip_widths*1.6, full_height(6)*1.4)
ggp <- make_boxplots(bugs.pcl %>% pcl.nicenames, active=F, logspace=T, reorder=F,
                     limits=c(10^-3, 10^4), zerobar_width=0.025, sig=sig_bugs_clr, sig_name="FDR",
                     rev(head(sig_bugs_clr$ID, 6))) +
    guides(fill="none", color="none") + ylab(NULL) +
    theme(axis.ticks.y = element_blank())
plot(ggp)
dev.off()
