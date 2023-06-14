
# Loads EC tables (in relative abundances)
#   ec.dna.strat.pcl: Stratified EC abundances in DNA (| in features)
#   ec.rna.strat.pcl: Stratified EC abundances in RNA (| in features)

if (!exists("ec.dna.strat.pcl")) {
    source("./env_config.r")
    source("./common/pcl_utils.r")
    source("./common/fix_metadata.r")

    library(dplyr)

    # Load data
    ec.dna.strat.pcl <- pcl.read(file.path(HMP2_data, "mgx", "ecs_relab.pcl.tsv.gz"), metadata.rows=" ")
    ec.rna.strat.pcl <- pcl.read(file.path(HMP2_data, "mtx", "ecs_relab.pcl.tsv.gz"), metadata.rows=" ")

    # Filter out samples with low seq depth
    stopifnot(all(!is.na(ec.dna.strat.pcl$meta$reads_filtered)))
    stopifnot(all(!is.na(ec.rna.strat.pcl$meta$reads_filtered)))
    ec.dna.strat.pcl <- ec.dna.strat.pcl %>% fix_metadata %>%
        pcl.filter.s(reads_filtered >= 1e6)
    ec.rna.strat.pcl <- ec.rna.strat.pcl %>% fix_metadata %>%
        pcl.filter.s(reads_filtered >= 1e6)
}
