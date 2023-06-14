
source("./common/theme_nature.r")
source("./common/disease_colors.r")
source("./common/disease_activity.r")

source("./common/load_bugs.r")
source("./common/load_kos.r")
source("./common/load_metabolites.r")
source("./common/load_biopsies.r")
source("./common/load_proteins.r")

pcls <- list(
    "Taxonomy" = bugs.pcl %>% pcl.only(rank="s") %>% pcl.normalize,
    "MTX (KOs)" = ko.rna.unstrat.pcl %>% merge_disease_activity(lenience=2),
    "MPX (KOs)" = proteins.kos.pcl %>% pcl.normalize %>% merge_disease_activity(lenience=2),
    "MBX" = metabolites.pcl.nrm %>% merge_disease_activity(lenience=2),
    "Biopsy (16S)" = biopsy_16s.pcl %>% pcl.filter.s(!is.na(diagnosis)) %>% pcl.normalize %>% merge_disease_activity(lenience=4))

tsne.ords <- lapply(pcls, pcl.tsne)

make_plot <- function(pcl, tsne) {
    pcl.ordplot(pcl, tsne, colour="subject", size_abs=0.6*2,
                outline_size=0.1, connectwidth = 0.25,
                shape="diagnosis") +
        theme_nature() +
        theme(axis.ticks=element_blank(), axis.text=element_blank()) +
        scale_fill_manual(values=rainbow(length(unique(pcl$meta$subject)))) +
        guides(fill="none", shape="none")
}

library(egg)
plotlist <- list()
for (i in seq_along(tsne.ords)) {
    plotlist <- c(plotlist, list(make_plot(pcls[[i]], tsne.ords[[i]]) + ggtitle(names(pcls)[i])))
}
ggp <- ggarrange(plots=plotlist, ncol=3)

pdf("./overview/ed_tsne_plots_person.pdf", 7.20472, 5.22441)
print(ggp)
dev.off()


# Per-disease group t-SNEs

plotlist <- list()
for (i in seq_along(pcls)) {
    for (dg in c("nonIBD", "UC", "CD")) {
        pcl <- pcls[[i]] %>% pcl.filter.s(!is.na(diagnosis) && diagnosis==dg)
        pcl$meta$subject <- factor(pcl$meta$subject)
        tsne.ord <- pcl.tsne(pcl)
        ggp <- pcl.ordplot(pcl, tsne.ord, colour="subject", size_abs=0.6*2,
                           outline_size=0.1, connectwidth = 0.25, shape="active") +
            theme_nature() +
            theme(axis.ticks=element_blank(), axis.text=element_blank()) +
            scale_fill_manual(values=sample(rainbow(132))) +
            scale_shape_manual(values=c("FALSE"=21, "TRUE"=24)) +
            guides(fill="none", shape="none")

        if (i == 1) {
            ggp <- ggp + ggtitle(dg)
        }
        if (dg == "nonIBD") {
            ggp <- ggp + ylab(names(pcls)[i])
        }

        plotlist <- c(plotlist, list(ggp))
    }
}

library(egg)
ggp <- ggarrange(plots=plotlist, ncol=3)

pdf("./overview/ed_tsne_plots_person_perdisease.pdf", 7.20472, 5.22441*5/2)
print(ggp)
dev.off()


# ED 2 array of plots
library(cowplot)
plotlist <- list()
empty <- ggplot()
for (i in seq_along(pcls)) {
    pco.ord <- pcl.pcoa(pcls[[i]])
    if ("biopsy_location" %in% colnames(pcls[[i]]$meta)) {
        pcls[[i]]$meta$ilcolrect <- factor(rep("Colon", pcls[[i]]$ns), levels=c("Colon", "Ileum", "Rectum"))
        pcls[[i]]$meta$ilcolrect[pcls[[i]]$meta$biopsy_location=="Ileum"] <- "Ileum"
        pcls[[i]]$meta$ilcolrect[pcls[[i]]$meta$biopsy_location=="Rectum"] <- "Rectum"
        pcls[[i]]$meta$ilcolrect[is.na(pcls[[i]]$meta$biopsy_location)] <- NA
        pco.plot <- pcl.ordplot(pcls[[i]], pco.ord, size_abs=1.5, outline_size=0.12,
                    colour = "diagnosis", colour_override = hmp2_disease_colors,
                    size = "ilcolrect", size_override = c(Colon=1.5/2, Ileum=1.5/2, Rectum=1.5), # Size is funny due to Inkscape importer
                    shape = "ilcolrect", shape_override = c(Colon=22, Ileum=24, Rectum=21))
    } else {
        pco.plot <- pcl.ordplot(pcls[[i]], pco.ord, size_abs=1.5, outline_size=0.12,
            colour = "diagnosis", colour_override = hmp2_disease_colors)
    }
    pco.plot <- pco.plot + theme_nature() + theme(axis.text.x=element_blank(), axis.text.y=element_blank()) +
        guides(fill="none", shape="none", size="none")
    tsne.plot <- make_plot(pcls[[i]], tsne.ords[[i]])

    df <- cbind(data.frame(x=pco.ord$points[,1], y=pco.ord$points[,2]), pcls[[i]]$meta)
    dens_theme <- list(theme_nature(), theme(axis.text.x=element_blank(), axis.text.y=element_blank()))
    x_density <- ggplot(data=df, aes(x, fill=diagnosis)) +
        geom_density(alpha=0.3, size=0.25) + dens_theme +
        scale_fill_manual(values=hmp2_disease_colors_fill) +
        guides(fill="none") + xlab(NULL) + ylab("Density")
    y_density <- ggplot(data=df, aes(y, fill=diagnosis)) +
        geom_density(alpha=0.3, size=0.25) + coord_flip() + dens_theme +
        scale_fill_manual(values=hmp2_disease_colors_fill) +
        guides(fill="none") + xlab(NULL) + ylab("Density")


    plotlist <- c(plotlist, list(x_density, empty, empty, pco.plot, y_density, tsne.plot))
}

library(cowplot)
pdf("./overview/ed2_ordinations.pdf", 4.1, 11)
cowplot::plot_grid(plotlist=plotlist, rel_widths=c(1,.3,1), rel_heights=c(.3,1), ncol=3)
dev.off()

pdf("./overview/ed2_ordinations_transpose.pdf", 7.20472*5/4, 4.1*7.20472/11*5/4)
pl2 <- plotlist[unclass(t(matrix(seq_along(plotlist), nrow=3)))]
pl2[seq(2, 10, by=2)+9] <- pl2[seq(2, 10, by=2)]
pl2[seq(2, 10, by=2)] <- list(empty)
pl2[seq(21, 30, by=2)] <- pl2[seq(22, 30, by=2)]
pl2[seq(22, 30, by=2)] <- list(empty)
cowplot::plot_grid(plotlist=pl2, rel_heights=c(.3,1,1), rel_widths=c(1,.3), nrow=3)
dev.off()


