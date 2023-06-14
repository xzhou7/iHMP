

# Load biopsy-related tables
#   biopsy_16s.pcl: OTU count matrix
#   biopsy_htx.counts.pcl: HTX count matrix

if (!exists("biopsy_16s.pcl")) {
    source("./env_config.r")
    source("./common/pcl_utils.r")
    source("./common/fix_metadata.r")

    library(plyr)
    library(dplyr)

    # Read the biospy 16S table
    biopsy_16s.pcl <- pcl.read(file.path(HMP2_data, "16s", "biopsy_16s_otus_Rfriendly.pcl.tsv"), metadata.rows=" ") %>%
        pcl.filter.s(sum(x) > 3000) %>% # read depth filter
        fix_metadata

    # Read the HTX table
    biopsy_htx.counts.pcl <- pcl.read(file.path(HMP2_data, "htx", "host_tx_counts.pcl.tsv"), metadata.rows=" ") %>%
        pcl.filter.s(sum(x) > 1e7) %>%
        fix_metadata
}


