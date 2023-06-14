
# Loads KO tables (in relative abundances)
#   ko.dna.unstrat.pcl: Unstratified KO abundances in DNA (NO | in features)
#   ko.rna.unstrat.pcl: Unstratified KO abundances in RNA (NO | in features)
#   bugs.fromko.dna.pcl: Bug-level rollup of the stratified KO DNA table
#   bugs.fromko.rna.pcl: Bug-level rollup of the stratified KO RNA table

if (!exists("ko.dna.unstrat.pcl")) {
    source("./env_config.r")
    source("./common/pcl_utils.r")
    source("./common/fix_metadata.r")

    library(dplyr)

    # Load data
    ko.dna.unstrat.pcl <- pcl.read(file.path(HMP2_data, "mgx", "kos_relab.slim.pcl.tsv.gz"), metadata.rows=" ")
    ko.rna.unstrat.pcl <- pcl.read(file.path(HMP2_data, "mtx", "kos_relab.slim.pcl.tsv.gz"), metadata.rows=" ")

    # Filter out samples with low seq depth
    stopifnot(all(!is.na(ko.dna.unstrat.pcl$meta$reads_filtered)))
    stopifnot(all(!is.na(ko.rna.unstrat.pcl$meta$reads_filtered)))
    ko.dna.unstrat.pcl <- ko.dna.unstrat.pcl %>% fix_metadata %>%
        pcl.filter.s(reads_filtered >= 1e6)
    ko.rna.unstrat.pcl <- ko.rna.unstrat.pcl %>% fix_metadata %>%
        pcl.filter.s(reads_filtered >= 1e6)

    # Load bug-level summaries of KO expression data
    # These files are generated in overview/src/bug_expression_rollup.r
    bugs.fromko.dna.pcl <- pcl.read(file.path(HMP2_data, "mgx", "kos_relab_bugsummary.pcl.tsv"), metadata.rows=" ") %>%
        fix_metadata
    bugs.fromko.rna.pcl <- pcl.read(file.path(HMP2_data, "mtx", "kos_relab_bugsummary.pcl.tsv"), metadata.rows=" ") %>%
        fix_metadata
}
