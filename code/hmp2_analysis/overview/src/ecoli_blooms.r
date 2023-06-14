
source("./common/pcl_utils.r")

# Assumes bugs.pcl is loaded

df <- data.frame(entero=bugs.pcl.nn$x[,"Enterobacteriaceae"],
                 week=bugs.pcl.nn$meta$week_num,
                 subject=bugs.pcl.nn$meta$subject,
                 diagnosis=bugs.pcl.nn$meta$diagnosis)

library(RColorBrewer)
ggp <- ggplot(data=df, aes(x=week, y=entero, colour=subject)) +
    geom_line(aes(group=subject)) + geom_point() +
    scale_color_manual(values=colorRampPalette(rev(brewer.pal(n = 7, name = "Dark2")))(length(levels(df$subject)))) +
    guides(colour="none") +
    facet_grid(diagnosis ~ .) +
    ylab("Enterobacteriaceae") + xlab("Week")

pdf("./overview/enterobacteriaceae_blooms.pdf", 7, 7)
print(ggp)
dev.off()


characterize_blooms <- function(bloom_thresh) {
    braycurtis <- function(A, B) 2*sum(pmin(A,B)) / (sum(A) + sum(B))

    n_blooms <- rep(0, length(levels(df$diagnosis)))
    names(n_blooms) <- levels(df$diagnosis)
    bloom_time <- n_blooms
    bloom_bc <- data.frame()
    for (sub in levels(bugs.pcl.nn$meta$subject)) {
        subpcl <- pcl.filter.s(bugs.pcl, keep = bugs.pcl.nn$meta$subject == sub)
        subpcl <- pcl.reorder.s(subpcl, order(subpcl$meta$week))

        dg <- subpcl$meta$diagnosis[1]

        bloom <- pcl.nicenames(subpcl)$x[,"Enterobacteriaceae"] >= bloom_thresh
        blooms <- which(diff(bloom) > 0) + 1
        n_blooms[dg] <- n_blooms[dg] + length(blooms)
        bloom_time[dg] <- bloom_time[dg] + max(subpcl$meta$week) - min(subpcl$meta$week)

        # Species-only PCL with Enterobacteriaceae removed and renormalized
        sp_noentero <- subpcl %>%
            pcl.filter.f(keep=!grepl("Enterobacteriaceae", colnames(subpcl$x))) %>%
            pcl.only(rank="s") %>%
            pcl.normalize

        for (bli in blooms) {
            # find the end of the bloom
            bloom_end <- bli + 1
            while (bloom_end <= length(bloom) && bloom[bloom_end])
                bloom_end <- bloom_end + 1

            # measure BC's
            bloom_bc <- rbind(bloom_bc, data.frame(
                diagnosis = dg,
                before_during = braycurtis(sp_noentero$x[bli,], sp_noentero$x[bli-1,]),
                before_after = if (bloom_end > length(bloom)) { NA } else {
                    braycurtis(sp_noentero$x[bli-1,], sp_noentero$x[bloom_end,]) },
                during_after = if (bloom_end > length(bloom)) { NA } else {
                    braycurtis(sp_noentero$x[bli,], sp_noentero$x[bloom_end,]) })
            )
        }
    }

    # Overall bloom rate
    blooms_per_year <- 52 * n_blooms / bloom_time

    sink(sprintf("./overview/enterobacteriaceae_bloom_stats_%.2f.txt", bloom_thresh))
    cat("Automatically generated.\n\n")
    cat(sprintf("Bloom Threshold: %.2f\n", bloom_thresh))
    cat("\nBlooms observed\n")
    print(n_blooms)
    cat("\nTime observed (weeks)\n")
    print(bloom_time)
    cat("\nBloom Rate (blooms / 52weeks)\n")
    print(blooms_per_year)
    cat("\n\nAll Bray-Curtis measures are species-level and with Enterobacteriaceae\nremoved and the table renormalized.\n")
    cat("\nbefore: The sample immediately before the bloom")
    cat("\nduring: The first sample in the bloom")
    cat("\nafter:  The first sample after the bloom\n\n")
    print(split(bloom_bc, bloom_bc$diagnosis))
    cat("\nSummaries:\n")
    print(sapply(split(bloom_bc, bloom_bc$diagnosis), summary))
    sink()
}

sapply(seq(0.1, 0.5, by=0.05), characterize_blooms)

