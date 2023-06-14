
source("./common/bug_colors.r")
source("./common/merge_metadata.r")

# Aerotolerance theme
aerotolerance_colors <- c(obligate_anaerobe="palegreen",
                          aerotolerant="royalblue",
                          facultative_anaerobe="orangered",
                          unknown="gray60")

# Load the aerotolerance list
aerotolerance_tbl <- read.table(file.path(HMP2_data, "metadata/aerotolerance.txt"), sep="\t", header=F)
aerotolerance <- as.character(aerotolerance_tbl[,2])
names(aerotolerance) <- as.character(aerotolerance_tbl[,1])
rm(aerotolerance_tbl)

# Add genera that are in the bug table, but not the aerotolerance table as "unknown"
source("./common/load_bugs.r")
source("./common/disease_activity.r")
all_genera <- colnames((bugs.pcl %>% pcl.only(rank="g") %>% pcl.nicenames)$x) %>%
    gsub(" unclassified", "", .)
missing_genera <- setdiff(all_genera, names(aerotolerance))
missing_genera_unk <- rep("unknown", length(missing_genera))
names(missing_genera_unk) <- missing_genera
aerotolerance <- factor(c(aerotolerance, missing_genera_unk), levels=names(aerotolerance_colors))



ordination_plot <- function(pcl, ord, subject) {
    issubj <- pcl$meta$subject == subject
    if (!("active" %in% colnames(pcl$meta))) {
        pcl <- merge_disease_activity(pcl, lenience=2)
    }
    interesting <- ifelse(issubj, ifelse(pcl$meta$active, "SA", "SI"), "UI")
    subjpts <- which(issubj)[order(pcl$meta$collection[issubj])]

    ggp <- ggplot() +
        geom_point(data=data.frame(pc1=ord$points[!issubj,1], pc2=ord$points[!issubj,2], inter=interesting[!issubj]),
                   aes(pc1, pc2, fill=inter, size=inter), shape=21, stroke=0.13) +
        geom_path(data=data.frame(x=ord$points[subjpts,1], y=ord$points[subjpts,2]), aes(x,y), size=.9) +
        geom_point(data=data.frame(pc1=ord$points[subjpts,1], pc2=ord$points[subjpts,2], inter=interesting[subjpts]),
                   aes(pc1, pc2, fill=inter, size=inter), shape=21, stroke=0.13) +
        scale_fill_manual(values=c(UI="gray80", SI="dodgerblue", SA="firebrick")) +
        scale_size_manual(values=c(UI=1.5, SI=3, SA=3)) +
        xlab(ord$ordnames[1]) + ylab(ord$ordnames[2]) +
        theme_cowplot() + guides(fill="none", size="none") +
        theme(axis.text.x=element_blank(), axis.text.y=element_blank())

    return (ggp)
}

trimstring <- function(s, limit, ...) {
    # Helper to trim the length of a string at a specified width in real-world units (not just character length)
    if (length(s) > 1) {
        ret <- s
        for (i in seq_along(s)) {
            ret[i] <- trimstring(s[i], limit, ...)
        }
        return (ret)
    } else {
        s <- as.character(s)
        if (!is.finite(limit) || strwidth(s, ...) < limit) {
            return (s)
        } else {
            ellipses <- "..."
            cw <- sapply(1:nchar(s), function(i)strwidth(paste(substr(s, 1, i), ellipses, sep=""), ...))
            cutat <- min(which(cw > limit))
            return (paste(substr(s, 1, cutat-1), ellipses, sep=""))
        }
    }
}

heatmap_tracks <- function(df, colors=c(Yes="firebrick", No="gray90", Blank="gray30"), xax=NA) {
    library(reshape2)
    library(ggplot2)
    library(cowplot)

    df$ID <- rownames(df)
    df.melt <- suppressWarnings(melt(df, id.vars="ID")) # suppresses attributes are not identical across measure variables
    colnames(df.melt) <- c("Sample", "Feature", "val")
    df.melt$val[df.melt$val==""] <- "Blank"
    if (length(xax)>1 || !is.na(xax)) {
        df.melt$Sample <- factor(df.melt$Sample, levels=xax)
    }

    ggp <- ggplot(data=df.melt) +
        geom_tile(aes(x=Sample, y=Feature, fill=val)) +
        theme_cowplot() + ylab(NULL) +
        scale_fill_manual(values=colors, breaks=names(colors), name=NULL)

    return (ggp)
}

stacked_bar <- function(pcl, colors, maxn=10, norm=T, xax=NA) {
    library(reshape2)
    library(ggplot2)
    library(cowplot)

    if (pcl$ns == 0 || pcl$nf == 0) {
        return (ggplot())
    }

    pcl <- pcl %>% pcl.sort.s(collection)
    if (norm) {
        pcl <- pcl %>% pcl.normalize
    }

    if (missing(colors)) {
        namedx <- pcl$x

        library(RColorBrewer)
        color_themes <- scale_fill_manual(values=brewer.pal(12, "Set3"), name=NULL)
        other <- rep(0, pcl$ns)
    } else {
        namedCol <- colnames(pcl$x) %in% names(colors)
        other <- rowSums(pcl$x[, !namedCol, drop=F])
        namedx <- pcl$x[, namedCol, drop=F]

        color_themes <- scale_fill_manual(values=colors, name=NULL)
    }

    featX <- apply(namedx, 2, max, na.rm=T)
    keepCol <- featX >= sort(featX, T)[min(length(featX), maxn)]
    other <- other + rowSums(namedx[, !keepCol, drop=F])
    bestx <- namedx[, keepCol, drop=F]

    if (any(other>0)) {
        xnr <- cbind(bestx, other)
        colnames(xnr)[ncol(xnr)] <- "other"
    } else {
        xnr <- bestx
    }


    rownames(xnr) <- rownames(pcl$meta)
    df <- melt(xnr)
    colnames(df) <- c("Sample", "Feature", "y")
    df$x <- pcl$meta$collection[match(df$Sample, rownames(pcl$meta))]
    if (length(xax) > 1 || !is.na(xax)) {
        df$x <- factor(df$x, levels=xax)
        fill <- xax[!(xax %in% df$x)]
        df <- rbind(df, data.frame(
            Feature=rep(df$Feature[1], length(fill)),
            Sample=rep(df$Sample[1], length(fill)),
            y=rep(0, length(fill)),
            x=fill))
    }
    if (!missing(colors) && !is.null(names(colors))) {
        df$Feature <- factor(df$Feature, levels=names(colors))
    }
    #df$Sample <- factor(df$Sample, levels=rownames(pcl$meta)[order(pcl$meta$collection)])

    df$Feature <- trimstring(df$Feature, 3, units="inches")
    df$y[is.na(df$y)] <- 0
    ggp <- ggplot(data=df) +
        geom_bar(aes(x=Sample, fill=Feature, y=y), width=1, stat="identity") +
        color_themes + facet_grid(.~x, scales="free_x") +
        theme_cowplot() +
        theme(strip.background=element_blank(), strip.text=element_blank(),
              panel.spacing=unit(0, "pt"), strip.switch.pad.grid=unit(0, "pt"))

    return (ggp)
}

varstbl_lineplot <- function(pcl, nstable=5, nvariable=6, range=c(), xax=NA, norm=F, logy=F) {
    library(reshape2)
    library(ggplot2)
    library(cowplot)

    if (pcl$ns == 0 || pcl$nf == 0)
        return (ggplot())

    pcl <- pcl %>%
        pcl.sort.s(collection) %>%
        pcl.filter.f(any(x>0))

    if (norm)
        pcl <- pcl %>% pcl.normalize

    mu <- pcl.apply.f(pcl, mean(x))
    sigma2 <- pcl.apply.f(pcl, var(x))

    if (pcl$ns > 2) {
        keep <- union(names(head(sort(sigma2/mu^2), nstable)), names(head(sort(sigma2,T), nvariable)))
    } else {
        keep <- names(head(sort(mu, T), nstable+nvariable))
    }

    pcl.keep <- pcl %>% pcl.filter.f(keep=colnames(pcl$x) %in% keep)

    df <- melt(pcl.keep$x)
    colnames(df) <- c("Sample", "Feature", "y")
    df$x <- pcl$meta$collection[match(df$Sample, rownames(pcl$meta))]

    if (length(xax) > 1 || !is.na(xax)) {
        df$x <- factor(df$x, levels=xax)
    }

    df$Feature <- as.character(df$Feature)
    df$Feature <- trimstring(df$Feature, 3.5, units="inches")

    ggp <- ggplot(data=df, aes(x=x, color=Feature, y=y)) +
        geom_line(aes(group=Feature)) +
        geom_point() +
        theme_cowplot() + scale_color_brewer(palette="Paired", name=NULL) +
        theme(strip.background=element_blank(), strip.text=element_blank(),
              panel.spacing=unit(0, "pt"), strip.switch.pad.grid=unit(0, "pt")) +
        scale_x_discrete(drop=F)
    if (length(range) > 0) {
        ggp <- ggp + scale_y_continuous(limits=c(range))
    }
    if (logy) {
        ggp <- ggp + scale_y_log10()
    }
    ggp
    return (ggp)
}

just_lineplot <- function(x, y, range=c(), xax=NA, logy=F) {
    library(reshape2)
    library(ggplot2)
    library(cowplot)

    if (length(x) == 0 || all(is.na(y))) {
        return (ggplot())
    }

    df <- data.frame(y=y, x=x)
    df <- df[!is.na(df$y),]

    if (length(xax) > 1 || !is.na(xax)) {
        df$x <- factor(df$x, levels=xax)
    }

    ggp <- ggplot(data=df, aes(x=x, y=y)) +
        geom_line(aes(group=1)) +
        geom_point() +
        theme_cowplot() + scale_color_brewer(palette="Paired") +
        theme(strip.background=element_blank(), strip.text=element_blank(),
              panel.spacing=unit(0, "pt"), strip.switch.pad.grid=unit(0, "pt")) +
        scale_x_discrete(drop=F)
    if (length(range) > 0) {
        ggp <- ggp + scale_y_continuous(limits=c(range))
    }
    if (logy) {
        ggp <- ggp + scale_y_log10()
    }
    ggp
    return (ggp)
}

all_timeseries_for_subject <- function(subject) {
    plots <- list()
    plot_relheight <- c()

    # What timepoints to show?
    mto <- hmp2_sample_metadata[hmp2_sample_metadata$Participant.ID == subject,]
    mto <- mto[grepl("^\\w\\d\\d\\d\\dC\\d\\d?$", mto$site_sub_coll),]
    mto$collection <- as.numeric(gsub("^.*C(\\d+)$", "\\1", mto$site_sub_coll))
    mt <- mto[!duplicated(mto$site_sub_coll),]
    mt <- mt[order(mt$collection),]
    cols <- mt$collection

    # Overall information
    plt_title <- sprintf("%s: %.0f %s %s %s | %s %s", subject,
                         mt$consent_age[1], mt$sex[1], mt$race[1], mt$site_name[1],
                         mt$diagnosis[1], mt$baseline_montreal_location[1])

    # Dietary tracks
    mtd <- mt[,c("Soft.drinks..tea.or.coffee.with.sugar..corn.syrup..maple.syrup..cane.sugar..etc.",
                 "Diet.soft.drinks..tea.or.coffee.with.sugar..Stevia..Equal..Splenda.etc.",
                 "Fruit.juice..orange..apple..cranberry..prune.etc..",
                 "Water",
                 "Alcohol..beer..brandy..spirits..hard.liquor..wine..aperitif..etc..",
                 "Yogurt.or.other.foods.containing.active.bacterial.cultures..kefir..sauerkraut.",
                 "Dairy..milk..cream..ice.cream..cheese..cream.cheese.",
                 "Probiotic",
                 "Fruits..no.juice...Apples..raisins..bananas..oranges..strawberries..blueberries",
                 "Vegetables..salad..tomatoes..onions..greens..carrots..peppers..green.beans..etc.",
                 "Beans..tofu..soy..soy.burgers..lentils..Mexican.beans..lima.beans.etc.",
                 "Whole.grains..wheat..oats..brown.rice..rye..quinoa..wheat.bread..wheat.pasta.",
                 "Starch..white.rice..bread..pizza..potatoes..yams..cereals..pancakes..etc..",
                 "Eggs",
                 "Processed.meat..other.red.or.white.meat.such.as.lunch.meat..ham..salami..bologna",
                 "Red.meat..beef..hamburger..pork..lamb.",
                 "White.meat..chicken..turkey..etc..",
                 "Shellfish..shrimp..lobster..scallops..etc..",
                 "Fish..fish.nuggets..breaded.fish..fish.cakes..salmon..tuna..etc..",
                 "Sweets..pies..jam..chocolate..cake..cookies..etc..")]
    colnames(mtd) <- c("Sweet drinks", "Diet drinks", "Fruit juice", "Water", "Alcohol", "Yogurt", "Dairy", "Probiotic",
                       "Fruits", "Vegetables", "Beans", "Whole grains", "Starch",
                       "Eggs", "Proc. meat", "Red meat", "White meat", "Shellfish", "Fish", "Sweets")
    rownames(mtd) <- cols
    mtd <- colwise(function(x) {
        xf <- factor(x, levels=c("",
                                 "No, I did not consume these products in the last 7 days",
                                 "Within the past 4 to 7 days",
                                 "Within the past 2 to 3 days",
                                 "Yesterday, 1 to 2 times",
                                 "Yesterday, 3 or more times"))
        levels(xf) <- c("Blank", "None 7+ days", "Within 4-7 days", "Within 2-3 days", "Yesterday 1-2x", "Yesterday 3x+")
        return (xf)
    })(mtd)

    diet_tracks <- heatmap_tracks(mtd, xax=cols, colors=c(
        "Blank" = "gray60",
        "None 7+ days" = "#EDF8E9", # RColorBrewer brewer.pal("Greens", n=5)
        "Within 4-7 days" = "#BAE4B3",
        "Within 2-3 days" = "#74C476",
        "Yesterday 1-2x" = "#31A354",
        "Yesterday 3x+" = "#006D2C"
    ))
    plots <- c(plots, list(diet_tracks))
    plot_relheight <- c(plot_relheight, 1.9)

    # Medication
    mtv <- mt[,c("X5..In.the.past.2.weeks..have.you.been.hospitalized.",
                 "X1..In.the.past.2.weeks..have.you.received.any.of.the.following.medications",
                 "Antibiotics",
                 "Chemotherapy",
                 "Immunosuppressants..e.g..oral.corticosteroids.")]
    colnames(mtv) <- c("Hospitalized", "Medication", "Antibiotics", "Chemo", "Immunosupp.")
    rownames(mtv) <- cols

    isact <- unlist(lapply(split(mto$External.ID %in% rownames(bugs.pcl$x)[bugs.pcl$meta$active], mto$collection), any))
    isinmgx <- unlist(lapply(split(mto$External.ID %in% rownames(bugs.pcl$x), mto$collection), any))
    #mtv <- cbind(Dysbiotic=c("","No","Yes")[1 + (cols %in% names(which(isinmgx))) + (cols %in% names(which(isact)))], mtv)

    metadata_tracks <- heatmap_tracks(mtv, xax=cols)
    plots <- c(plots, list(metadata_tracks))
    plot_relheight <- c(plot_relheight, 0.6)

    ## TODO: Shift annotation?

    # Fecalcal
    fecalcal_line <- just_lineplot(cols, mt$fecalcal, xax=cols, range=c(0, 500)) +
        ylab("FecalCal")
    plots <- c(plots, list(fecalcal_line))
    plot_relheight <- c(plot_relheight, 0.5)

    # HBI/SCCAI
    if (mt$diagnosis[1] == "CD") {
        hbi_line <- just_lineplot(cols, mt$hbi, xax=cols, range=c(0, 20)) +
            ylab("HBI")
        plots <- c(plots, list(hbi_line))
    } else if (mt$diagnosis[1] == "UC") {
        sccai_line <- just_lineplot(cols, mt$sccai, xax=cols, range=c(0, 12)) +
            ylab("SCCAI")
        plots <- c(plots, list(sccai_line))
    } else {
        plots <- c(plots, list(ggplot()))
    }
    plot_relheight <- c(plot_relheight, 0.5)

    # Disease activity score
    bugs <- bugs.pcl %>%
        pcl.filter.s(keep=bugs.pcl$meta$subject == subject) %>%
        pcl.filter.f(any(x>0))
    act <- bugs$meta$activity_index
    onething <- data.frame(x=1)
    activity_line <- just_lineplot(bugs$meta$collection, act, xax=cols, range=c(0.5, 1)) +
        ylab("Dysbiosis") +
        geom_hline(data=onething, aes(yintercept=disease_activity_threshold, linetype="Dysbiosis"), color="red") +
        geom_hline(data=onething, aes(yintercept=eubiosis_lower_threshold, linetype="Eubiosis"), color="dodgerblue") +
        scale_linetype_manual(name="", values=c("dashed", "dashed"),
                              guide = guide_legend(override.aes = list(color = c("red", "dodgerblue"))))
    plots <- c(plots, list(activity_line))
    plot_relheight <- c(plot_relheight, 0.5)

    # Serology
    if (exists("serology.eu.pcl")) {
        ser <- serology.eu.pcl %>%
            pcl.filter.s(keep=!is.na(serology.eu.pcl$meta$subject) & serology.eu.pcl$meta$subject == subject)
        serology_lines <- varstbl_lineplot(ser, xax=cols, logy=F, range=c(0, 100)) +
            ylab("Serology")
        plots <- c(plots, list(serology_lines))
        plot_relheight <- c(plot_relheight, 0.5)
    }

    # MGX Sequencing depth
    if (exists("bugs.pcl")) {
        seqbugs <- merge_metadata(bugs, c("reads_qc_fail", "reads_human", "reads_filtered"), "metagenomics")
        if (any(is.na(seqbugs$meta$frac_human_reads))) {
            warning(sprintf("Samples with NA frac_human_reads: %s",
                            paste(rownames(seqbugs$meta)[is.na(seqbugs$meta$frac_human_reads)], collapse=" ")))
            seqbugs$meta$frac_human_reads[is.na(seqbugs$meta$frac_human_reads)] <- 0
        }
        seq_depth <- pcl.make(seqbugs$meta[,c("reads_human", "reads_filtered"),drop=F], meta=seqbugs$meta)
        colnames(seq_depth$x) <- c("Human", "Filtered")
        seqdepth_bars <- stacked_bar(seq_depth, xax=cols, norm=F,
            colors=c("Human" = "gray70", "Filtered" = "dodgerblue")) +
            ylab("MGX Reads")

        plots <- c(plots, list(seqdepth_bars))
        plot_relheight <- c(plot_relheight, 0.5)
    }

    # MTX Sequencing depth
    if (exists("ko.rna.unstrat.pcl")) {
        mtx.ko <- ko.rna.unstrat.pcl %>%
            pcl.filter.s(keep=ko.rna.unstrat.pcl$meta$subject == subject) %>%
            merge_metadata(c("reads_qc_fail", "reads_ribosomal", "reads_human", "reads_filtered"), "metatranscriptomics")
        if (any(is.na(mtx.ko$meta$frac_human_reads))) {
            warning(sprintf("Samples with NA frac_human_reads: %s",
                            paste(rownames(mtx.ko$meta)[is.na(mtx.ko$meta$frac_human_reads)], collapse=" ")))
            mtx.ko$meta$frac_human_reads[is.na(mtx.ko$meta$frac_human_reads)] <- 0
        }
        mtx_seq_depth <- pcl.make(mtx.ko$meta[,c("reads_human", "reads_filtered"),drop=F], meta=mtx.ko$meta)
        colnames(mtx_seq_depth$x) <- c("Human", "Filtered")
        mtx_seq_depth_bars <- stacked_bar(mtx_seq_depth, xax=cols, norm=F,
            colors=c("Human" = "gray70", "Filtered" = "dodgerblue")) +
            ylab("MTX Reads")

        plots <- c(plots, list(mtx_seq_depth_bars))
        plot_relheight <- c(plot_relheight, 0.5)
    }

    plots_left <- plots
    plot_left_relheight <- plot_relheight
    plots <- list()
    plot_relheight <- c()

    # Bugs
    species_bars <- stacked_bar(bugs %>% pcl.only(rank="s") %>% pcl.nicenames, hmp2_bug_colours, xax=cols) +
        ylab("Species") +
        theme(legend.text = element_text(face = "italic"))
    plots <- c(plots, list(species_bars))
    plot_relheight <- c(plot_relheight, 1)
    #genera_bars <- stacked_bar(bugs %>% pcl.only(rank="g") %>% pcl.nicenames, hmp2_bug_colours_genus, xax=cols) +
    #    ylab("Genus Abundance")
    # plots <- c(plots, list(genera_bars))
    # plot_relheight <- c(plot_relheight, c(1))

    # Bug transcription
    if (exists("bugs.fromko.rna.pcl")) {
        bug_transc <- bugs.fromko.rna.pcl %>%
            pcl.filter.s(keep=bugs.fromko.rna.pcl$meta$subject == subject) %>%
            pcl.filter.f(any(x>0))
        species_transc_bars <- stacked_bar(bug_transc, hmp2_bug_colours, xax=cols) +
            ylab("Transcription") +
            theme(legend.text = element_text(face = "italic"))

        plots <- c(plots, list(species_transc_bars))
        plot_relheight <- c(plot_relheight, 1)
    }

    # Aerotolerance
    aerotolerance_bars <- stacked_bar(bugs %>%
        pcl.only(rank="g") %>% pcl.nicenames %>% pcl.group.f(aerotolerance, mapping=T),
        aerotolerance_colors, xax=cols) + ylab("Aerotolerance")
    plots <- c(plots, list(aerotolerance_bars))
    plot_relheight <- c(plot_relheight, 0.5)

    # Viruses
    viruses <- viruses.pcl %>%
        pcl.filter.s(keep=viruses.pcl$meta$subject == subject) %>%
        pcl.filter.f(any(x>0))
    viromics_bars <- stacked_bar(viruses %>% pcl.nicenames, xax=cols, maxn=8) +
        ylab("Viruses")
    plots <- c(plots, list(viromics_bars))
    plot_relheight <- c(plot_relheight, 0.8)


    # MTX KO's
    if (exists("ko.rna.unstrat.pcl")) {
        mtx.ko <- ko.rna.unstrat.pcl %>%
            pcl.filter.s(keep=ko.rna.unstrat.pcl$meta$subject == subject)
        mtx_lines <- varstbl_lineplot(mtx.ko, xax=cols, logy=T) +
            ylab("Transcription")

        plots <- c(plots, list(mtx_lines))
        plot_relheight <- c(plot_relheight, c(1))
    }

    # Metabolites
    if (exists("metabolites.named.nonredundant.pcl")) {
        mtb <- metabolites.named.nonredundant.pcl %>%
            pcl.filter.s(keep=metabolites.named.nonredundant.pcl$meta$subject == subject)
        mbx_lines <- varstbl_lineplot(mtb, xax=cols, logy=T) +
            ylab("Metabolites")

        plots <- c(plots, list(mbx_lines))
        plot_relheight <- c(plot_relheight, c(1))
    }

    # Proteins
    if (exists("proteins.kos.pcl")) {
        ptx <- proteins.kos.pcl %>%
            pcl.filter.s(keep=proteins.pcl$meta$subject == subject)
        ptx_lines <- varstbl_lineplot(ptx, xax=cols, norm=T, logy=T) +
            ylab("Proteins")

        plots <- c(plots, list(ptx_lines))
        plot_relheight <- c(plot_relheight, c(1))
    }

    # # Pathways - DNA
    # if (exists("pwy.dna.unstrat.pcl")) {
    #     pwys_dna <- pwy.dna.unstrat.pcl %>%
    #         pcl.filter.s(keep=pwy.dna.unstrat.pcl$meta$subject == subject) %>%
    #         pcl.filter.f(any(x>0))
    #
    #     pwydna_lines <- varstbl_lineplot(pwys_dna, norm=T, logy=T, xax=cols) +
    #         ylab("Pathways (DNA)")
    #
    #     plots <- c(plots, list(pwydna_lines))
    #     plot_relheight <- c(plot_relheight, c(1))
    # }
    #
    # # Pathways - RNA
    # if (exists("pwy.rna.unstrat.pcl")) {
    #     pwys_rna <- pwy.rna.unstrat.pcl %>%
    #         pcl.filter.s(keep=pwy.rna.unstrat.pcl$meta$subject == subject) %>%
    #         pcl.filter.f(any(x>0))
    #
    #     pwyrna_lines <- varstbl_lineplot(pwys_rna, norm=T, logy=T, xax=cols) +
    #         ylab("Pathways (RNA)")
    #
    #     plots <- c(plots, list(pwyrna_lines))
    #     plot_relheight <- c(plot_relheight, c(1))
    # }

    # Ordinations
    tax_ord <- ordination_plot(species.pcl, species.ord, subject) + ggtitle("Taxonomy")
    mtx_ord <- ordination_plot(ko.rna.unstrat.pcl, mtx.ord, subject) + ggtitle("Transcripts")
    mpx_ord <- ordination_plot(proteins.kos.pcl, mpx.ord, subject) + ggtitle("Proteins")
    mbx_ord <- ordination_plot(metabolites.pcl.nrm, mbx.ord, subject) + ggtitle("Metabolites")
    plots_ord <- list(tax_ord, mtx_ord, mpx_ord, mbx_ord)

    # Put it all together
    library(egg)
    plots_left_clean <- lapply(plots_left, function(x) x + xlab(NULL) +
                              theme(axis.text.x=element_blank(), legend.key.size=unit(4.2, "mm"),
                                    legend.text=element_text(size=11)))
    plots_clean <- lapply(plots, function(x) x + xlab(NULL) +
                              theme(axis.text.x=element_blank(), legend.key.size=unit(4.2, "mm"),
                                    legend.text=element_text(size=11)))
    ggps_left <- ggarrange(plots=plots_left_clean, ncol=1, heights=c(plot_left_relheight), draw=F)
    ggps <- ggarrange(plots=plots_clean, ncol=1, heights=c(plot_relheight), draw=F)
    ggps_right <- ggarrange(plots=plots_ord, ncol=1, draw=F)
    ggpt <- ggdraw() + draw_label(plt_title)

    ggp_left <- plot_grid(ggpt, ggps_left, ncol=1, rel_heights=c(0.2, sum(plot_relheight)))
    ggp <- plot_grid(ggp_left, ggps, ggps_right, nrow=1, rel_widths=c(1, 1.22, 2.22*3.5/19))

    #ggp
    return (ggp)
}


source("./common/load_bugs.r")
source("./common/load_pathways.r")
source("./common/load_metabolites.r")
source("./common/load_proteins.r")
source("./common/load_serology.r")
source("./common/load_kos.r")

# Run the first part of define_disease_activity to include activity tracks

# metabolites.pcl.nrm.named <- metabolites.pcl.nrm %>%
#     pcl.filter.f(keep=metabolites.pcl.nrm$metaf$Metabolite!="")
# colnames(metabolites.pcl.nrm.named$x) <- sprintf("%s-%s",
#     metabolites.pcl.nrm.named$metaf$Method, metabolites.pcl.nrm.named$metaf$Metabolite)


# Generate ordinations
species.pcl <- bugs.pcl %>% pcl.only(rank="s") %>% pcl.normalize
species.ord <- pcl.pcoa(species.pcl)
mtx.ord <- pcl.pcoa(ko.rna.unstrat.pcl %>% pcl.normalize)
mpx.ord <- pcl.pcoa(proteins.kos.pcl %>% pcl.normalize)
mbx.ord <- pcl.pcoa(metabolites.pcl.nrm)

#ordination_plot(species.pcl, species.ord, "H4006")

pdf("./overview/subject_timeseries_test.pdf", 19, 12)
print(all_timeseries_for_subject("H4006"))
print(all_timeseries_for_subject("E5001"))
#print(all_timeseries_for_subject("C3019"))
dev.off()

pdf("./overview/subject_timeseries.pdf", 19, 12)
for (sub in levels(hmp2_sample_metadata$Participant.ID)) {
    if (sum(hmp2_sample_metadata$data_type[hmp2_sample_metadata$Participant.ID == sub] == "metagenomics") >= 2) {
        print(all_timeseries_for_subject(sub))
    }
}
dev.off()




