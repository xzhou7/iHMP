
# Load metabolite tables
#   metabolites.pcl: Raw median-scaled table of all metabolites
#   metabolites.pcl.nrm: All metabolites, abundances normalized within each method
#   metabolites.named.pcl: Median-scaled table subsetted subsetted to named
#       metabolies by taking the representative with minimal CV2 in the
#       pooled runs.

if (!exists("metabolites.pcl")) {
    source("./env_config.r")
    source("./common/pcl_utils.r")
    source("./common/fix_metadata.r")

    library(plyr)
    library(dplyr)

    # Load in the giant table
    file <- gzfile(file.path(HMP2_data, "mbx", "iHMP_metabolomics.pcl.csv.gz"), "r")
    raw <- read.table(file, sep=",", header=F, stringsAsFactors=F, quote="\"")
    close(file)

    # Figure out how much metadata there is
    nmetaf <- 0
    while (raw[1, nmetaf+1] == "")
        nmetaf <- nmetaf + 1
    nmetas <- 0
    while (raw[nmetas+1, 1] == "")
        nmetas <- nmetas + 1

    # Cut out the metadata and actual tables
    metas <- as.data.frame(t(raw[1:nmetas, (nmetaf+2):ncol(raw)]))
    colnames(metas) <- raw[1:nmetas, nmetaf+1]
    rownames(metas) <- raw[nmetas+1, (nmetaf+2):ncol(raw)]
    metaf <- as.data.frame(raw[(nmetas+2):nrow(raw), 1:nmetaf])
    colnames(metaf) <- raw[nmetas+1, 1:nmetaf]
    rownames(metaf) <- raw[(nmetas+2):nrow(raw), nmetaf+1]
    data <- as.matrix(raw[(nmetas+2):nrow(raw), (nmetaf+2):ncol(raw)]) %>% t
    mode(data) <- "numeric"
    data[is.na(data)] <- 0
    colnames(data) <- rownames(metaf)
    rownames(data) <- rownames(metas)

    # Make a PCL object out of it all
    metabolites.pcl <- pcl.make(data, metas, metaf)

    # Clean up
    rm(raw, nmetas, nmetaf, metas, metaf, data)

    # Apply common metadata fixes
    metabolites.pcl <- fix_metadata(metabolites.pcl)

    # Normalize within-method
    metabolites.pcl.nrm <- metabolites.pcl
    imputed <- apply(metabolites.pcl$x, 2, median) # what to impute for crazy values?
    nfixed <- 0
    for (i in seq_along(rownames(metabolites.pcl.nrm$x))) {
        prf <- metabolites.pcl.nrm$x[i,]
        nprf <- prf / ave(prf, metabolites.pcl.nrm$metaf$Method, FUN=sum)
        outlier <- nprf > 0.5
        while (any(outlier)) {
            # Fix crazy outliers
            nfixed <- nfixed + sum(outlier)
            prf[outlier] <- imputed[outlier]
            nprf <- prf / ave(prf, metabolites.pcl.nrm$metaf$Method, FUN=sum)
            outlier <- nprf > 0.5
        }
        metabolites.pcl.nrm$x[i,] <- nprf
    }
    if (nfixed > 0) {
        warning(sprintf("Imputed %d unusually large values in the normalized metabolomics data.", nfixed))
    }
    rm(imputed, prf, nprf, outlier, nfixed)

    # Subset to named metabolites
    metabolites.named.rd.pcl <- pcl.filter.f(metabolites.pcl, keep=metabolites.pcl$metaf$Metabolite != "")

    # Subset again to non-redundant measurements of named metabolites
    unq_metabolites <- sort(unique(metabolites.named.rd.pcl$metaf$Metabolite))
    unq_metabolites_keepidx <- rep(NA, length(unq_metabolites))
    for (i in seq_along(unq_metabolites_keepidx)) {
        mtI <- which(metabolites.named.rd.pcl$metaf$Metabolite == unq_metabolites[i])
        unq_metabolites_keepidx[i] <- mtI[which.min(metabolites.named.rd.pcl$metaf$`Pooled QC sample CV`[mtI])]
    }
    metabolites.named.pcl <- pcl.reorder.f(metabolites.named.rd.pcl, unq_metabolites_keepidx)
    colnames(metabolites.named.pcl$x) <- metabolites.named.pcl$metaf$Metabolite
    rownames(metabolites.named.pcl$metaf) <- metabolites.named.pcl$metaf$Metabolite

    # Clean up
    rm(unq_metabolites, unq_metabolites_keepidx, mtI, metabolites.named.rd.pcl)
}






