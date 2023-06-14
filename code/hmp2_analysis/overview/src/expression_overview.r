
source("./env_config.r")
source("./common/pcl_utils.r")
source("./common/disease_colors.r")
source("./common/bug_colors.r")

library(dplyr)


characterize_mtx_mgx_relationship <- function(name, mgx.pcl, mtx.pcl) {
    make_df <- function(mgx, mtx) {
        df <- data.frame(feature = names(mgx), mgx = mgx, mtx = mtx)
        return (df)
    }

    # Overall plot
    mtx.pcl <- pcl.match.f(mtx.pcl, mgx.pcl, unmatched.val=0)
    keep <- mgx.pcl$x > 0
    mgx.pcl$x[!keep] <- NA # code mgx zeros as NA so we remove them from both calculations
    mtx.pcl$x[!keep] <- NA
    df <- make_df(mgx.pcl %>% pcl.apply.f(mean(x, na.rm=T)), mtx.pcl %>% pcl.apply.f(mean(x, na.rm=T)))
    df$feat_lbl <- df$feature
    df$feat_lbl[!((df$mgx>2e-7) & (df$mtx>2e-7))] <- NA
    df$feat_lbl_sparse <- df$feature
    nzmtxmgx <- (df$mgx>0) & (df$mtx>0)
    df$ratio <- df$mtx / df$mgx
    df$feat_lbl_sparse[!((df$mgx>sort(df$mgx, T)[15]) | (df$mtx>sort(df$mtx, T)[15]) |
        (nzmtxmgx & ((df$ratio>=sort(df$ratio[nzmtxmgx], T)[15]) | (df$ratio<=sort(df$ratio[nzmtxmgx])[15]))))] <- NA

    library(ggplot2)
    library(ggrepel)
    ggp <- ggplot(data=df, aes(x=mgx, y=mtx)) +
        geom_abline(intercept=0, slope=1, colour="gray") +
        geom_point(shape=21, size=1.5, fill="royalblue") +
        scale_x_log10() + scale_y_log10()


    pdf(sprintf("./overview/expressionmean_%s_log.pdf", name), 7, 6)
    print(ggp + geom_text_repel(aes(label=feat_lbl_sparse), size=2, box.padding=unit(0.1, "lines")))
    dev.off()

    pdf(sprintf("./overview/expressionmean_%s_log_big.pdf", name), 15, 15)
    print(ggp + geom_text_repel(aes(label=feat_lbl), size=2, box.padding=unit(0.1, "lines")))
    dev.off()

    # Per-disease stats
    df.cd <- make_df(mgx.pcl %>% pcl.filter.s(diagnosis=="CD") %>% pcl.apply.f(mean(x, na.rm=T)),
                     mtx.pcl %>% pcl.filter.s(diagnosis=="CD") %>% pcl.apply.f(mean(x, na.rm=T)))
    df.cd$diagnosis <- "CD"
    df.uc <- make_df(mgx.pcl %>% pcl.filter.s(diagnosis=="UC") %>% pcl.apply.f(mean(x, na.rm=T)),
                     mtx.pcl %>% pcl.filter.s(diagnosis=="UC") %>% pcl.apply.f(mean(x, na.rm=T)))
    df.uc$diagnosis <- "UC"
    df.ni <- make_df(mgx.pcl %>% pcl.filter.s(diagnosis=="nonIBD") %>% pcl.apply.f(mean(x, na.rm=T)),
                     mtx.pcl %>% pcl.filter.s(diagnosis=="nonIBD") %>% pcl.apply.f(mean(x, na.rm=T)))
    df.ni$diagnosis <- "nonIBD"
    df <- rbind(df.cd, df.uc, df.ni)
    perimiter <-
        sqrt((log10(df.cd$mgx) - log10(df.uc$mgx))^2 + (log10(df.cd$mtx) - log10(df.uc$mtx))^2) +
        sqrt((log10(df.cd$mgx) - log10(df.ni$mgx))^2 + (log10(df.cd$mtx) - log10(df.ni$mtx))^2) +
        sqrt((log10(df.uc$mgx) - log10(df.ni$mgx))^2 + (log10(df.uc$mtx) - log10(df.ni$mtx))^2)
    names(perimiter) <- df.cd$feature
    df$feat_lbl <- df$feature
    df$feat_lbl[df$diagnosis!="nonIBD"] <- NA

    triangle_plot <- list(
        geom_abline(intercept=0, slope=1, colour="gray"),
        geom_polygon(aes(group=feature), alpha=0.3, fill="orange"),
        geom_point(aes(colour=diagnosis)), scale_colour_manual(values=hmp2_disease_colors),
        guides(alpha="none"),
        #geom_path(aes(group=feature)),
        scale_x_log10(), scale_y_log10()
    )

    thr <- 0
    keep <- (df.cd$mgx>thr) & (df.uc$mgx>thr) & (df.ni$mgx>thr) &
        (df.cd$mtx>thr) & (df.uc$mtx>thr) & (df.ni$mtx>thr)
    ggp <- ggplot(data=df[keep,], aes(x=mgx, y=mtx)) + triangle_plot

    pdf(sprintf("./overview/expressionmean_%s_log_disease.pdf", name), 7, 6)
    print(ggp)
    dev.off()
    pdf(sprintf("./overview/expressionmean_%s_log_disease_labels.pdf", name), 10, 9)
    print(ggp + geom_text_repel(aes(label=feat_lbl), size=2, box.padding=unit(0.1, "lines")))
    dev.off()

    # Dump perimiter statistics
    write.table(as.data.frame(sort(perimiter[is.finite(perimiter)], decreasing=T)), quote=F, sep="\t",
                file=sprintf("./overview/expressionmean_%s_perimiters.txt", name))

    # Triangle plot for largest changes
    ggp <- ggplot(aes(x=mgx, y=mtx),
        data=df[df$feature %in% names(head(sort(perimiter[is.finite(perimiter)], decreasing=T), 10)),]) +
        triangle_plot + geom_text_repel(aes(label=feat_lbl), size=2, box.padding=unit(0.1, "lines"))

    pdf(sprintf("./overview/expressionmean_%s_log_disease_10largest.pdf", name), 7, 6)
    print(ggp)
    dev.off()

    # Triangle plot for smallest changes
    ggp <- ggplot(aes(x=mgx, y=mtx),
        data=df[df$feature %in% names(head(sort(perimiter[is.finite(perimiter)]), 15)),]) +
        triangle_plot + geom_text_repel(aes(label=feat_lbl), size=2, box.padding=unit(0.1, "lines"))

    pdf(sprintf("./overview/expressionmean_%s_log_disease_15smallest.pdf", name), 7, 6)
    print(ggp + geom_text_repel(aes(label=feat_lbl), size=2, box.padding=unit(0.1, "lines")))
    dev.off()
}

# Characterize bug expression
source("./common/load_pathways.r")
#characterize_mtx_mgx_relationship("bugs",
#    pwy.dna.mt.strat.pcl %>% pcl.group.f(colnames(pcl.nicenames(pwy.dna.mt.strat.pcl)$x)),
#    pwy.rna.mt.strat.pcl %>% pcl.group.f(colnames(pcl.nicenames(pwy.rna.mt.strat.pcl)$x))
#)
source("./common/load_ecs.r")
characterize_mtx_mgx_relationship("bugs-ec", bug.dna.mt.pcl, bug.rna.mt.pcl)

# Characterize pathway expression
characterize_mtx_mgx_relationship("pathways", pwy.dna.mt.unstrat.pcl, pwy.rna.mt.unstrat.pcl)


# function to show stratified bug distributions
plot_stratified_distribution <- function(pwy, data, norm=F, sorting, show=NA,
                                         data_2ndsort, wt_2ndsort=0.0001, wt_abundsort=0.01,
                                         totdata, totdata_2ndsort,
                                         no.unclassified=norm, no.unaccounted=norm,
                                         sortmethod="tsp") {

    # pwy is a unique substring of the pathway name
    # data is the data pcl containing stratified relative abundances
    # norm is whether to normalize the data *among stratified abundances*
    # sorting is a vector of sample names used to impose an external sorting
    #         if not given, the sorting options below determing sample order
    # show is a vector of bug names used to choose what bugs to show
    #      if not given, all bugs with colours are displayed
    # data_2ndsort is a pcl object providing the secondary sorting of the samples
    # wt_2ndsort is the weight of the second sorting, relative to the first
    # wt_abundsort is the weight of the abundance sorting, relative to the first
    # totdata is a pcl of the unstratified pathway abundances
    # totdata_2ndsort is a pcl of the unstratified pathway abundances for the
    #                 2nd sort (for abundance sorting)
    # sortmethod is "tsp" or "pcoa", and determines the method used to
    #            determine sample ordering


    set.seed(666)

    # Put out the pathway's stratified abundances
    sub_pwy <- pcl.filter.f(data, keep=grepl(sprintf("%s|", pwy), colnames(data$x), fixed=T))
    # Clean up the names
    sub_pwy <- pcl.nicenames(sub_pwy)

    # Use the "total data" to create the "unaccounted" category
    if (!missing(totdata)) {
        sub_pwy.tot <- pcl.filter.f(totdata, keep=grepl(pwy, colnames(totdata$x), fixed=T))
        if (!no.unaccounted) {
            newcol <- matrix(sub_pwy.tot$x - pcl.apply.s(sub_pwy, sum(x)), sub_pwy.tot$ns, 1)
            rownames(newcol) <- rownames(sub_pwy.tot$x)
            colnames(newcol) <- "unaccounted"
            sub_pwy$x <- cbind(sub_pwy$x, newcol)
            sub_pwy$nf <- sub_pwy$nf + 1
        }
    }

    # Remove unclassified
    if (no.unclassified) {
        sub_pwy <- pcl.filter.f(sub_pwy, Name != "unclassified")
    }

    # Normalize the data if necessary
    sub_pwyNN <- sub_pwy
    if (norm) {
        sub_pwy <- pcl.normalize(sub_pwy)
    }

    # Determine and apply the sample sorting
    if (missing(sorting) || is.null(sorting)) {
        if (sub_pwy$nf == 0) {
        } else if (sortmethod == "tsp") {
            # get the inter-sample distance matrix
            library(labdsv)
            D <- as.matrix(dsvdis(asin(sqrt(sub_pwy$x)), "bray/curtis"))
            n <- dim(D)[1]

            # hack to force all-zero samples to be different
            totAb <- pcl.apply.s(sub_pwy, sum(x))
            ZeroSim <- matrix(kronecker(totAb==0, totAb==0, FUN="&"), n, n) & (diag(n)<1)
            D <- D + ZeroSim

            # abundance sort
            sub_pwyAB <- if (missing(totdata)) {
                pcl.apply.s(sub_pwyNN, sum(x))
            } else {
                sub_pwy.tot$x
            }
            sub_pwyAB <- sub_pwyAB / max(sub_pwyAB)
            Da <- wt_abundsort * matrix(kronecker(sub_pwyAB, sub_pwyAB, FUN="-"), n, n)^2
            D <- D + Da * wt_abundsort

            # secondary-sorting
            if (!missing(data_2ndsort)) {
                sub_pwy2 <- pcl.filter.f(data_2ndsort, keep=grepl(sprintf("%s|", pwy), colnames(data_2ndsort$x), fixed=T))
                sub_pwy2 <- pcl.nicenames(sub_pwy2)
                sub_pwy2 <- pcl.filter.f(sub_pwy2, Name != "unclassified")

                sub_pwyAB2 <- if (missing(totdata_2ndsort)) {
                    pcl.apply.s(sub_pwy2, sum(x))
                } else {
                    sub_pwy2.tot <- pcl.filter.f(totdata_2ndsort, keep=grepl(pwy, colnames(totdata_2ndsort$x), fixed=T))
                    sub_pwy2.tot$x
                }
                sub_pwyAB2 <- drop(sub_pwyAB2 / max(sub_pwyAB2))
                Da2 <- wt_abundsort * outer(sub_pwyAB2, sub_pwyAB2, FUN="-")^2
                D <- D + Da2 * wt_2ndsort

                sub_pwy2 <- pcl.normalize(sub_pwy2)
                D2 <- as.matrix(dsvdis(asin(sqrt(sub_pwy2$x)), "bray/curtis"))
                totAb2 <- pcl.apply.s(sub_pwy2, sum(x))
                ZeroSim2 <- outer(totAb2==0, totAb2==0, FUN="&") & (diag(n)<1)
                D2 <- D2 + ZeroSim2
                D <- D + D2 * wt_2ndsort
            }

            # find order by solving the TSP
            D <- rbind(cbind(D, rep(0, n)), rep(0, n+1)) # augment with fence
            library(TSP)
            tsp <- TSP(D, c(seq_along(rownames(sub_pwy$x)), 0))
            tour <- as.numeric(labels(solve_TSP(tsp, method="cheapest_insertion", reps=10)))
            library(wavethresh)
            tour <- guyrot(tour, -which(tour==0))[1:n] # rotate fence to correct position and remove

            # reorder based on the tour order
            sub_pwy <- pcl.reorder.s(sub_pwy, tour)

        } else if (sortmethod == "pcoa") {
            pca <- pcl.pcoa(sub_pwy)
            sub_pwy$x <- sub_pwy$x[order(pca$points[,1]),, drop=F]
        }
        sorting <- rownames(sub_pwy$x)
    } else {
        sub_pwy$x <- sub_pwy$x[match(sorting, rownames(sub_pwy$x)),, drop=F]
    }

    # Determine the set of bugs to show
    if (is.na(show) || is.numeric(show)) {
        show <- colnames(pcl.top.f(pcl.filter.f(sub_pwy, Name %in% names(hmp2_bug_colours)),
                                   mean(x), n=if(is.na(show)){29}else{show})$x)
        show <- names(hmp2_bug_colours)[names(hmp2_bug_colours) %in% show]
    }
    otherify_map <- colnames(sub_pwy$x)
    names(otherify_map) <- otherify_map
    otherify_map[!(otherify_map %in% show)] <- "other"
    sub_pwy <- pcl.group.f(sub_pwy, otherify_map)

    # Melt the data
    library(reshape2)
    library(RColorBrewer)
    sub_pwy.melt <- melt(sub_pwy$x)
    colnames(sub_pwy.melt) <- c("Sample", "Bug", "Abundance")
    # attach metadata
    sub_pwy.melt <- cbind(sub_pwy.melt, sub_pwy$meta[match(sub_pwy.melt$Sample, rownames(sub_pwy$meta)),])

    # Determine bug and color order
    buglevels <- show
    if ("other" %in% colnames(sub_pwy$x)) {
        buglevels <- c(buglevels, "other")
    }
    bugcolors = hmp2_bug_colours[buglevels]
    sub_pwy.melt$Bug <- factor(as.character(sub_pwy.melt$Bug), levels=names(bugcolors))

    # Actual plotting work
    library(ggplot2)
    ggp <- ggplot(sub_pwy.melt, aes(x=Sample, y=Abundance, fill=Bug)) +
        geom_bar(stat="identity", width=1) +
        scale_fill_manual(values=bugcolors, guide=guide_legend(ncol=1)) +
        # scale_fill_manual(guide=guide_legend(nrow = 1,
        #                                      direction = "horizontal",
        #                                      label.position = "bottom",
        #                                      label.hjust = 1, label.vjust = 0.5,
        #                                      title.theme = element_blank(),
        #                                      label.theme = element_text(angle = 90, face="italic")),
        # values=bugcolors) +
        theme(axis.text.x=element_blank(),
              legend.text = element_text(face = "italic"))

    return (list(ggp = ggp, sorting = sorting, show=show))
}

# helper to extract the legend grobs from a ggplot object
get_legend<-function(myggplot) {
    tmp <- ggplot_gtable(ggplot_build(myggplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
}

mgx_mtx_stratified_plot <- function(pwy) {

    legW <- 0.28
    rn <- plot_stratified_distribution(pwy, pwy.rna.mt.strat.pcl, norm=T, data_2ndsort=pwy.dna.mt.strat.pcl,
                                       totdata=pwy.rna.mt.unstrat.pcl, totdata_2ndsort=pwy.dna.mt.unstrat.pcl,
                                       sortmethod="tsp")
    dn <- plot_stratified_distribution(pwy, pwy.dna.mt.strat.pcl, sorting=rn$sorting, norm=T,
                                       totdata=pwy.dna.mt.unstrat.pcl)
    leg <- get_legend(dn$ggp)

    facets <- list(facet_grid(. ~ diagnosis, scales="free_x", space="free"))

    library(cowplot)
    print(ggdraw() +
        draw_plot(dn$ggp + ylab("DNA") + theme(legend.position="none") + facets,
                  x=0, y=0,   w=1-legW, h=0.48) +
        draw_plot(rn$ggp + ylab("RNA") + ggtitle(sprintf("%s", pwy)) + facets +
                  theme(legend.position="none"), x=0, y=0.48, w=1-legW, h=0.52) +
        draw_plot(leg, x=1-legW, y=0, w=legW, h=1))


    r <- plot_stratified_distribution(pwy, pwy.rna.mt.strat.pcl, sorting=rn$sorting, totdata=pwy.rna.mt.unstrat.pcl)
    d <- plot_stratified_distribution(pwy, pwy.dna.mt.strat.pcl, sorting=rn$sorting, totdata=pwy.dna.mt.unstrat.pcl)
    leg <- get_legend(d$ggp)

    library(cowplot)
    print(ggdraw() +
      draw_plot(d$ggp + ylab("DNA") + theme(legend.position="none") + facets,
                x=0, y=0,   w=1-legW, h=0.48) +
      draw_plot(r$ggp + ylab("RNA") + ggtitle(sprintf("%s", pwy)) + facets +
                    theme(legend.position="none"), x=0, y=0.48, w=1-legW, h=0.52) +
      draw_plot(leg, x=1-legW, y=0, w=legW, h=1))
}


pdf("./overview/pathway_rna_dna_barplots_tsp.pdf", 15, 7.5, onefile=T)
allpwys_strat <- intersect(colnames(pwy.rna.mt.strat.pcl$x), colnames(pwy.dna.mt.strat.pcl$x))
allpwys_strat_nounc <- allpwys_strat[!grepl("\\|unclassified$", allpwys_strat)]
for (pwy in names(sort(pwy.dna.mt.unstrat.pcl %>% pcl.apply.f(mean(x)), T))) {
    if (any(grepl(pwy, allpwys_strat_nounc, fixed=T)))
        mgx_mtx_stratified_plot(pwy)
}
dev.off()





