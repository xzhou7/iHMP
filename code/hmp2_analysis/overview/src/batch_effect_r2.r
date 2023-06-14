
source("./common/disease_colors.r")
source("./common/load_bugs.r")
source("./common/load_kos.r")
source("./common/merge_metadata.r")

batch_plots_and_tests <- function(pcl, name, nicename, datatype) {
    # Merge in batch info
    pcl <- merge_metadata(pcl, c(batch = "PDO.Number", date="date_of_receipt"), datatype=datatype)
    pcl$meta$date_since_010113 <- as.numeric(difftime(strptime(pcl$meta$date, "%m/%d/%y"),
                                                      strptime("1/1/13", "%m/%d/%y")), units="days")

    library(ggplot2)
    library(cowplot)
    # Sample selection-like plots
    pdf(sprintf("./overview/batch_%s_sampleselection.pdf", name), 8, 17)
    print(ggplot(pcl$meta, aes(week_num, subject, color=batch)) +
        geom_point() +
        theme_cowplot() +
        ylab(NULL) + xlab("Week") +
        ggtitle(sprintf("Batch distribution across samples in %s", nicename)))
    print(ggplot(pcl$meta, aes(date_since_010113, subject, color=batch)) +
        geom_point() +
        theme_cowplot() +
        ylab(NULL) + xlab("Date since Jan 1 2013") +
        ggtitle(sprintf("Batch distribution across samples in %s", nicename)))
    dev.off()

    pdf(sprintf("./overview/batch_%s_distributions.pdf", name), 10, 8)
    # Disease group by batch
    print(ggplot(pcl$meta, aes(x=diagnosis, fill=diagnosis)) +
              geom_bar() +
              scale_fill_manual(values=hmp2_disease_colors) +
              facet_wrap(~ batch))
    # Disease group by batch
    print(ggplot(pcl$meta, aes(x=site_name, fill=site_name)) +
              geom_bar() +
              facet_wrap(~ batch))
    dev.off()

    library(vegan)
    pcl <- pcl.filter.s(pcl, !is.na(site_sub_coll))
    D <- vegdist(pcl.normalize(pcl)$x, method="bray")
    sink(sprintf("./overview/batch_%s_permanova.txt", name))
    print(adonis(D ~ ., data=pcl$meta[,"batch",drop=F],
                 permutations=how(nperm=499, blocks=pcl$meta$subject)))
    print(adonis(D ~ ., data=pcl$meta[,c("site_name","batch"),drop=F],
                 permutations=how(nperm=499, blocks=pcl$meta$subject)))
    print(adonis(D ~ ., data=pcl$meta[,c("subject","batch"),drop=F],
                 permutations=499))
    print(adonis(D ~ ., data=pcl$meta[,c("subject","batch"),drop=F],
                 permutations=how(nperm=499, blocks=pcl$meta$subject)))
    sink()


}


batch_plots_and_tests(bugs.pcl %>% pcl.only(rank="s"), "species", "MGX Species", "metagenomics")
batch_plots_and_tests(ko.dna.unstrat.pcl, "ko_mgx", "MGX KOs", "metagenomics")
batch_plots_and_tests(ko.rna.unstrat.pcl, "ko_mtx", "MTX KOs", "metatranscriptomics")
