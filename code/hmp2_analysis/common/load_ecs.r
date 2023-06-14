
# Loads EC tables (in relative abundances)
#   ec.dna.unstrat.pcl: Unstratified EC abundances in DNA (NO | in features)
#   ec.rna.unstrat.pcl: Unstratified EC abundances in RNA (NO | in features)

if (!exists("ec.dna.unstrat.pcl")) {
    source("./env_config.r")
    source("./common/pcl_utils.r")
    source("./common/fix_metadata.r")

    library(dplyr)

    # Load data
    ec.dna.unstrat.pcl <- pcl.read(file.path(HMP2_data, "mgx", "ecs_relab.slim.pcl.tsv.gz"), metadata.rows=" ")
    ec.rna.unstrat.pcl <- pcl.read(file.path(HMP2_data, "mtx", "ecs_relab.slim.pcl.tsv.gz"), metadata.rows=" ")

    # Filter out samples with low seq depth
    stopifnot(all(!is.na(ec.dna.unstrat.pcl$meta$reads_filtered)))
    stopifnot(all(!is.na(ec.rna.unstrat.pcl$meta$reads_filtered)))
    ec.dna.unstrat.pcl <- ec.dna.unstrat.pcl %>% fix_metadata %>%
        pcl.filter.s(reads_filtered >= 1e6)
    ec.rna.unstrat.pcl <- ec.rna.unstrat.pcl %>% fix_metadata %>%
        pcl.filter.s(reads_filtered >= 1e6)
}
