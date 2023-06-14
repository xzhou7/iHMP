
source("./common/load_metabolites.r")
source("./common/load_ecs.r")
source("./common/load_proteins.r")
source("./common/match_datasets.r")
source("./common/disease_colors.r")

read_attributevalue <- function(filename, ignore_fields = c("DBLINKS")) {
    data <- list()
    entry <- list()

    con <- file(filename, "r")
    while (T) {
        line <- readLines(con, n = 1)
        if (length(line) == 0) {
            break
        }
        if (line == "//") {
            # New entry
            if (length(entry) > 0) {
                data[[entry[["UNIQUE-ID"]]]] <- entry
                if (length(data) %% 2000 == 0) {
                    cat(sprintf("Reading %s: %d entries\n", filename, length(data)))
                    flush.console()
                }
            }
            entry <- list()
        } else if (startsWith(line, "#")) {
            # Comment... ignore
        } else {
            m <- regmatches(line, regexec("^([A-Z\\-]+) \\- (.*)$", line))[[1]]

            if (length(m) > 1 && !(m[2] %in% ignore_fields)) {
                if (m[2] %in% names(entry)) {
                    entry[[m[2]]] <- c(entry[[m[2]]], m[3])
                } else {
                    entry[[m[2]]] <- m[3]
                }
            }
        }
    }
    close(con)

    cat(sprintf("Finished reading %s: %d entries\n", filename, length(data)))
    return (data)
}

# Assumes these MetaCyc files are in ./data:
reactions.dat <- read_attributevalue("./data/reactions.dat")
compounds.dat <- read_attributevalue("./data/compounds.dat")
enzrxns.dat <- read_attributevalue("./data/enzrxns.dat")


# Build a compound name -> Metacyc UNIQUE-ID mapping
compoundmap <- do.call(c, unname(lapply(compounds.dat, function(x) {
    q <- tolower(c(x[["COMMON-NAME"]], x[["SYNONYMS"]]))
    q <- gsub("^l-", "", q)
    q <- gsub("</?[a-z]+>", "", q)
    q <- gsub("ic acid", "ate", q)
    #q <- gsub("[^a-z0-9 ]", ".", q)
    qid <- rep(x[["UNIQUE-ID"]], length(q))
    names(qid) <- q
    return (qid)
}))) %>% bmatch_sort

# Map metabolites to UNIQUE-IDs
metabolites.named.pcl$metaf$uniqueid <- unname(compoundmap[match(
    #gsub("[^a-z0-9 ]", ".", tolower(metabolites.named.pcl$metaf$Metabolite)),
    tolower(metabolites.named.pcl$metaf$Metabolite),
    names(compoundmap))])

# How many compounds match?
sum(!is.na(metabolites.named.pcl$metaf$uniqueid))
# [1] 181
mean(!is.na(metabolites.named.pcl$metaf$uniqueid))
# [1] 0.3284936

# These metabolites failed to match
head(metabolites.named.pcl$metaf$Metabolite[is.na(metabolites.named.pcl$metaf$uniqueid)], 20)


# Build an enzyme -> reaction ID mapping
enzrxnmap <- do.call(c, unname(lapply(enzrxns.dat, function(x) {
    q <- tolower(c(x[["COMMON-NAME"]], x[["SYNONYMS"]]))
    q <- gsub("</?[a-z]+>", "", q)
    q <- q[nchar(q) > 0]
    #q <- gsub("ic acid", "ate", q)
    qid <- rep(x[["REACTION"]], length(q))
    names(qid) <- q
    return (qid)
}))) %>% bmatch_sort


# MTX/MTX table
mt.rnadna <- match_datasets(list(ec.rna.unstrat.pcl, ec.dna.unstrat.pcl))
ec.common <- intersect(colnames(mt.rnadna[[1]]$x), colnames(mt.rnadna[[2]]$x))
ec.rnadna.unstrat.pcl <- pcl.reorder.f(mt.rnadna[[1]], ec.common)
ec.rnadna.unstrat.pcl$x <- ec.rnadna.unstrat.pcl$x / pcl.reorder.f(mt.rnadna[[2]], ec.common)$x



# Get a list of all enzymes in the data
ecs.all <- union(union(colnames(ec.dna.unstrat.pcl$x), colnames(ec.rna.unstrat.pcl$x)), colnames(proteins.ecs.pcl$x))

# Match the enzymes with their reactions
rxmt <- unname(enzrxnmap[match(tolower(gsub("^.*: ", "", ecs.all)), names(enzrxnmap))])

# Get sample matchings between all pairs of datasets
datasets <- list(metabolites.named.pcl, pcl.normalize(proteins.ecs.pcl), ec.rna.unstrat.pcl, ec.dna.unstrat.pcl)
names(datasets) <- c("MBX", "MPX", "MTX", "MGX")
matchings <- list()
for (i in seq_along(datasets)) {
    for (j in seq_along(datasets)) {
        if (i != j) {
            matchings[[i*length(datasets)+j]] <- match_datasets(datasets[c(i,j)], matching=T)
        }
    }
}
matching_rnadna_mbx <- match_datasets(list(ec.rnadna.unstrat.pcl, metabolites.named.pcl), matching=T)


library(ggplot2)
library(cowplot)
library(egg)

plot_multiomic_mbx_panel <- function(eci, cmpi) {
    cat(sprintf("\n\n%s -- %s\n", ecs.all[eci], metabolites.named.pcl$metaf$Metabolite[cmpi]))


    # Get appropriate scale for each axis
    scales <- lapply(datasets, function(pcl) c(min(pcl$x[pcl$x>0], na.rm=T), max(pcl$x, na.rm=T)))
    rnadna_scale <- c(min(ec.rnadna.unstrat.pcl$x[ec.rnadna.unstrat.pcl$x>0], na.rm=T), max(ec.rnadna.unstrat.pcl$x[is.finite(ec.rnadna.unstrat.pcl$x)], na.rm=T))

    mxcor <- 0
    dsi <- c(cmpi, match(ecs.all[eci], colnames(datasets[[2]]$x)), match(ecs.all[eci], colnames(datasets[[3]]$x)), match(ecs.all[eci], colnames(datasets[[4]]$x)))

    maxprev <- max(sapply(seq_along(datasets), function(i) mean(datasets[[i]]$x[,dsi[i]] > 0, na.rm=T))[-1], na.rm=T)
    maxprev <- min(maxprev, mean(datasets[[1]]$x[,dsi[1]] > 0))

    plt <- list()
    for (di in rev(seq_along(datasets))) {
        for (dj in rev(seq_along(datasets))) {
            if (di < dj) {
                mt <- matchings[[di*length(datasets)+dj]]
                df <- cbind(datasets[[di]]$meta[mt[,1],], data.frame(
                    y = if(is.na(dsi[di])) {0} else {unname(datasets[[di]]$x[mt[,1], dsi[di]])},
                    x = if(is.na(dsi[dj])) {0} else {unname(datasets[[dj]]$x[mt[,2], dsi[dj]])}))
                if (di == 1) {
                    mxcor <- max(mxcor, abs(cor(df$x, df$y, method="spearman", use="na.or.complete")))
                    cat(sprintf("%s:%s\n", names(datasets)[dj], names(datasets)[di]))
                    print(cor.test(df$x, df$y, method="spearman"))
                }
                ggp <- ggplot(data=df, aes(x=x, y=y)) +
                    geom_point(aes(fill=diagnosis), shape=21, size=1.6, stroke=0.1) +
                    scale_fill_manual(values=hmp2_disease_colors) +
                    guides(fill="none") +
                    scale_x_log10(limits=scales[[dj]]) +
                    scale_y_log10(limits=scales[[di]]) +
                    xlab(if (di == 1) {sprintf("EC (%s)", names(datasets)[dj])} else {NULL}) +
                    ylab(if (dj == length(datasets)) {
                        if (di == 1) {metabolites.named.pcl$metaf$Metabolite[cmpi]}
                        else {sprintf("EC (%s)", names(datasets)[di])}} else {NULL})
                if (di != 1) {
                    ggp <- ggp + theme(axis.text.x=element_blank())
                }
                if (dj != 4) {
                    ggp <- ggp + theme(axis.text.y=element_blank())
                }
                plt <- c(plt, list(ggp))
            }
        }
        if (di > 1 && di < length(datasets)) {
            plt <- c(plt, rep(list(ggplot()), di))
        }
    }

    mt <- matching_rnadna_mbx
    rnadna_mt <- match(ecs.all[eci], colnames(ec.rnadna.unstrat.pcl$x))
    df <- cbind(ec.rnadna.unstrat.pcl$meta[mt[,1],], data.frame(
        x = if (is.na(rnadna_mt)) {0} else {unname(ec.rnadna.unstrat.pcl$x[mt[,1], rnadna_mt])},
        y = unname(metabolites.named.pcl$x[mt[,2], dsi[1]])))
    if (sum(!is.nan(df$x) & !is.nan(df$y)) > 1) {
        mxcor <- max(mxcor, abs(cor(df$x, df$y, method="spearman", use="na.or.complete")))
        cat(sprintf("%s:%s\n", "RNADNA", "MBX"))
        print(cor.test(df$x, df$y, method="spearman"))
    }
    rnadna_plt <- ggplot(data=df, aes(x=x, y=y)) +
        geom_point(aes(fill=diagnosis), shape=21, size=1.6, stroke=0.1) +
        scale_fill_manual(values=hmp2_disease_colors) +
        guides(fill="none") +
        xlab("EC (MTX/MGX)") +
        scale_x_log10(limits=rnadna_scale) +
        scale_y_log10(limits=scales[[1]]) +
        theme(axis.text.y=element_blank()) +
        #ylab(metabolites.named.pcl$metaf$Metabolite[cmpi])
        ylab(NULL)
    plt <- c(plt, list(rnadna_plt))

    dtype_grid <- ggarrange(plots=plt, nrow=length(datasets)-1, draw=F)
    dtype_grid

    ggp <- ggdraw() +
        draw_plot(dtype_grid, 0, 0, 1, 1) +
        #draw_plot(rnadna_plt, 3/4, 0, 1/4, 1/3) +
        draw_label(sprintf("%s\n[%s]\nvs.\n%s [%s]", ecs.all[eci], rxmt[eci], metabolites.named.pcl$metaf$Metabolite[cmpi], metabolites.named.pcl$metaf$uniqueid[cmpi]),
                   1/4 + .5*3/4, 2/3 + .5*1/3, hjust=.5, vjust=.5, fontface="bold", size=15)

    return (list(ggp=ggp, mxcor=mxcor, maxprev=maxprev))
}


# For each EC...
plots <- list()
sortf <- c()
sink("./mechanistic_correlations/enzyme_mbx_correlations.txt")
for (i in seq_along(ecs.all)) {
    # .. if there's a matching MetaCyc reaction ..
    if (!is.na(rxmt[i])) {
        rx <- reactions.dat[[rxmt[i]]]

        compounds <- c(rx$LEFT, rx$RIGHT)
        #if (grepl("carnitine", ecs.all[i])) {
        #    compounds <- c(compounds, metabolites.named.pcl$metaf$uniqueid[grepl("carnitine", metabolites.named.pcl$metaf$Metabolite)])
        #}

        # .. check its substrate and products ...
        for (cmp in compounds) {
            cmpi <- match(cmp, metabolites.named.pcl$metaf$uniqueid)
            # .. if they're in our dataset ...
            if (!is.na(cmpi)) {
                # .. then plot the relationships between abundance in the different measurements
                pl <- plot_multiomic_mbx_panel(i, cmpi)

                if (pl$maxprev >= 0.1) {
                    sortf <- c(sortf, pl$mxcor)
                    plots <- c(plots, list(pl$ggp))
                }
            }
        }
    }
}
sink()

cat("Printing plots...\n")
pdf("./mechanistic_correlations/enzyme_mbx_scatters.pdf", 11+0.3, 11*3/4)
for (i in order(sortf, decreasing = T)) {
    print(plots[[i]])
}
dev.off()


sig_dnambx <- read.delim(file.path(HMP2_root, "analysis", "mechanistic_correlations", "halla_output", "ecDNAs-metabolites_residuals_full_q05_alla", "associations.txt"), header=T, sep="\t", stringsAsFactors=F)
sig_rnambx <- read.delim(file.path(HMP2_root, "analysis", "mechanistic_correlations", "halla_output", "ecRNAs-metabolites_residuals_full_q05_alla", "associations.txt"), header=T, sep="\t", stringsAsFactors=F)
sig_mpxmbx <- read.delim(file.path(HMP2_root, "analysis", "mechanistic_correlations", "halla_output", "proteinECs-metabolites_residuals_full_q05_alla", "associations.txt"), header=T, sep="\t", stringsAsFactors=F)


sig <- unique(rbind(
    sig_dnambx[sig_dnambx$qvalue<0.001,c("cluster1", "cluster2")],
    sig_rnambx[sig_rnambx$qvalue<0.001,c("cluster1", "cluster2")],
    sig_mpxmbx[sig_mpxmbx$qvalue<0.001,c("cluster1", "cluster2")]))
dim(sig)
# [1] 3838    2

# attach significance
sig$qvalue <- NA
for (i in seq_len(nrow(sig))) {
    c1 <- sig$cluster1[i]
    c2 <- sig$cluster2[i]
    minq <- min(sig_dnambx$qvalue[(sig_dnambx$cluster1==c1) & (sig_dnambx$cluster2==c2)],
                sig_rnambx$qvalue[(sig_rnambx$cluster1==c1) & (sig_rnambx$cluster2==c2)],
                sig_mpxmbx$qvalue[(sig_mpxmbx$cluster1==c1) & (sig_mpxmbx$cluster2==c2)])
    sig$qvalue[i] <- minq
}

# sort by significance
sig <- sig[order(sig$qvalue),]

plots <- list()
sortf <- c()
sink("./mechanistic_correlations/enzyme_mbx_correlations_mechcorr.txt")
for (i in 1:10) {#seq_len(nrow(sig))) {
    eci <- match(sig$cluster1[i], ecs.all)
    cmpi <- match(sig$cluster2[i], metabolites.named.pcl$metaf$Metabolite)
    if (!is.na(eci) && !is.na(cmpi)) {
        pl <- plot_multiomic_mbx_panel(eci, cmpi)

        qstr <- ""
        q <- sig_dnambx$qvalue[(sig_dnambx$cluster1==sig$cluster1[i]) & (sig_dnambx$cluster2==sig$cluster2[i])]
        qstr <- if (length(q) > 0) {sprintf("%sMGX:MBX q = %.2g\n", qstr, q)} else {qstr}
        q <- sig_rnambx$qvalue[(sig_rnambx$cluster1==sig$cluster1[i]) & (sig_rnambx$cluster2==sig$cluster2[i])]
        qstr <- if (length(q) > 0) {sprintf("%sMTX:MBX q = %.2g\n", qstr, q)} else {qstr}
        q <- sig_mpxmbx$qvalue[(sig_mpxmbx$cluster1==sig$cluster1[i]) & (sig_mpxmbx$cluster2==sig$cluster2[i])]
        qstr <- if (length(q) > 0) {sprintf("%sMGX:MBXo q = %.2g\n", qstr, q)} else {qstr}

        ggp <- pl$ggp +
            draw_label(qstr, 3/4, 1/2)
        sortf <- c(sortf, pl$mxcor)
        plots <- c(plots, list(ggp))
    }
}
sink()

cat("Printing plots...\n")
pdf("./mechanistic_correlations/enzyme_mbx_scatters_mechcorr.pdf", 11+0.3, 11*3/4)
for (i in order(sortf, decreasing = T)) {
    print(plots[[i]])
}
dev.off()
