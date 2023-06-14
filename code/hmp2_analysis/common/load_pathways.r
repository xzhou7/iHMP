
# Load pathway tables (relative abundances except for the ratio tables)
#   pwy.dna.unstrat.pcl: Unstratified DNA pathway table (no pathways with |)
#   pwy.rna.unstrat.pcl: Unstratified RNA pathway table (no pathways with |)
#   pwy.dna.strat.pcl: Stratified DNA pathway table (only pathways with |)
#   pwy.rna.strat.pcl: Stratified RNA pathway table (only pathways with |)
#   pwy.rnadna.unstrat.pcl: Unstratified RNA/DNA table
#   pwy.rnadna.strat.pcl: Stratified RNA/DNA table

# Additionally produces pwy.[dr]na.mt.(un)?strat.pcl: tables which have been
#   subset to a common set of pathways and samples for RNA and DNA.

if (!exists("pwy.dna.strat.pcl")) {
    source("./env_config.r")
    source("./common/pcl_utils.r")
    source("./common/fix_metadata.r")

    library(dplyr)

    # Load data
    pwy.dna.pcl <- pcl.read(file.path(HMP2_data, "mgx", "pathabundance_relab.pcl.tsv.gz"), metadata.rows=" ")
    pwy.rna.pcl <- pcl.read(file.path(HMP2_data, "mtx", "pathabundance_relab.pcl.tsv.gz"), metadata.rows=" ")

    # Filter out samples with low seq depth
    stopifnot(all(!is.na(pwy.dna.pcl$meta$reads_filtered)))
    stopifnot(all(!is.na(pwy.rna.pcl$meta$reads_filtered)))
    pwy.dna.pcl <- pwy.dna.pcl %>% fix_metadata %>%
        pcl.filter.s(reads_filtered >= 1e6)
    pwy.rna.pcl <- pwy.rna.pcl %>% fix_metadata %>%
        pcl.filter.s(reads_filtered >= 1e6)

    # Remove bizarre samples
    bizarre <- union(names(which((pwy.dna.pcl %>% pcl.apply.s(sum(x))) == 0)),
                     names(which((pwy.rna.pcl %>% pcl.apply.s(sum(x))) == 0)))
    if (length(bizarre) > 0) {
        warning(sprintf("These samples have 0 abundance in the RNA or DNA pwy tables: %s", do.call(paste, c(as.list(bizarre), list(sep=" ")))))
        pwy.dna.pcl <- pcl.filter.s(pwy.dna.pcl, keep=!(rownames(pwy.dna.pcl$x) %in% bizarre))
        pwy.rna.pcl <- pcl.filter.s(pwy.rna.pcl, keep=!(rownames(pwy.rna.pcl$x) %in% bizarre))
    }

    # Remove UNINTEGRATED and UNMAPPED as they're not consistently there in all parts of the dataset
    if (any(grepl("UNINTEGRATED", colnames(pwy.dna.pcl$x), fixed=T) | colnames(pwy.dna.pcl$x) == "UNMAPPED")) {
        warning("UNINTEGRATED/UNMAPPED rows in DNA pathway table. Removing.")
        pwy.dna.pcl <- pwy.dna.pcl %>% pcl.filter.f(keep=!grepl("UNINTEGRATED", colnames(pwy.dna.pcl$x), fixed=T) & colnames(pwy.dna.pcl$x) != "UNMAPPED")
    }
    if (any(grepl("UNINTEGRATED", colnames(pwy.rna.pcl$x), fixed=T) | colnames(pwy.rna.pcl$x) == "UNMAPPED")) {
        warning("UNINTEGRATED/UNMAPPED rows in RNA pathway table. Removing.")
        pwy.rna.pcl <- pwy.rna.pcl %>% pcl.filter.f(keep=!grepl("UNINTEGRATED", colnames(pwy.rna.pcl$x), fixed=T) & colnames(pwy.rna.pcl$x) != "UNMAPPED")
    }

    # Subset to common samples
    samps <- intersect(rownames(pwy.dna.pcl$x), rownames(pwy.rna.pcl$x))
    pwy.dna.mt.pcl <- pwy.dna.pcl %>% pcl.reorder.s(samps)
    pwy.rna.mt.pcl <- pwy.rna.pcl %>% pcl.reorder.s(samps)

    # Match features
    feats <- union(colnames(pwy.dna.mt.pcl$x), colnames(pwy.dna.mt.pcl$x))
    pwy.dna.mt.pcl <- pwy.dna.mt.pcl %>% pcl.reorder.f(feats, na.val=0)
    pwy.rna.mt.pcl <- pwy.rna.mt.pcl %>% pcl.reorder.f(feats, na.val=0)

    # Stratified/unstratified tables
    # Costs ~300MB to keep all tables in RAM currently (2017-08-09)
    pwy.dna.unstrat.pcl    <- pwy.dna.pcl    %>% pcl.filter.f(keep=!grepl("|", fixed=T, colnames(pwy.dna.pcl$x)))
    pwy.rna.unstrat.pcl    <- pwy.rna.pcl    %>% pcl.filter.f(keep=!grepl("|", fixed=T, colnames(pwy.rna.pcl$x)))
    pwy.dna.mt.unstrat.pcl <- pwy.dna.mt.pcl %>% pcl.filter.f(keep=!grepl("|", fixed=T, colnames(pwy.dna.mt.pcl$x)))
    pwy.rna.mt.unstrat.pcl <- pwy.rna.mt.pcl %>% pcl.filter.f(keep=!grepl("|", fixed=T, colnames(pwy.rna.mt.pcl$x)))

    pwy.dna.strat.pcl      <- pwy.dna.pcl    %>% pcl.filter.f(keep= grepl("|", fixed=T, colnames(pwy.dna.pcl$x)))
    pwy.rna.strat.pcl      <- pwy.rna.pcl    %>% pcl.filter.f(keep= grepl("|", fixed=T, colnames(pwy.rna.pcl$x)))
    pwy.dna.mt.strat.pcl   <- pwy.dna.mt.pcl %>% pcl.filter.f(keep= grepl("|", fixed=T, colnames(pwy.dna.mt.pcl$x)))
    pwy.rna.mt.strat.pcl   <- pwy.rna.mt.pcl %>% pcl.filter.f(keep= grepl("|", fixed=T, colnames(pwy.rna.mt.pcl$x)))

    # Transcription ratio tables
    pwy.rnadna.unstrat.pcl <- pwy.rna.mt.unstrat.pcl
    pwy.rnadna.unstrat.pcl$x <- pwy.rna.mt.unstrat.pcl$x / pwy.dna.mt.unstrat.pcl$x
    pwy.rnadna.strat.pcl   <- pwy.rna.mt.strat.pcl
    pwy.rnadna.strat.pcl$x <- pwy.rna.mt.strat.pcl$x / pwy.dna.mt.strat.pcl$x

    # Cleanup
    rm(samps, feats, pwy.dna.pcl, pwy.rna.pcl, pwy.dna.mt.pcl, pwy.rna.mt.pcl, bizarre)
}
