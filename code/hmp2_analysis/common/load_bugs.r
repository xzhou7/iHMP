

# Load Bugs and Virus tables
#   bugs.pcl: Full MetaPhlAn2 table (all taxonomic levels), in relative abundances
#   viruses.pcl: VirMAP viral abundance table, in counts

if (!exists("bugs.pcl")) {
    source("./env_config.r")
    source("./common/pcl_utils.r")
    source("./common/fix_metadata.r")

    library(dplyr)

    # Read in raw table
    bugs.pcl <- pcl.read(file.path(HMP2_data, "mgx", "taxonomic_profiles.pcl.tsv"), metadata.rows=" ")

    # Filter out low-read depth samples and fix metadata
    stopifnot(all(!is.na(bugs.pcl$meta$reads_filtered)))
    bugs.pcl <- bugs.pcl %>% fix_metadata %>%
        pcl.filter.s(reads_filtered >= 1e6)

    # Remove bizarre samples
    bizarre <- names(which(!pcl.apply.s(bugs.pcl, any((x>0) & (x<100)))))
    if (length(bizarre) > 0) {
        warning(sprintf("These samples have only one bug: %s", do.call(paste, c(as.list(bizarre), list(sep=" ")))))
        bugs.pcl <- pcl.filter.s(bugs.pcl, keep=!(rownames(bugs.pcl$x) %in% bizarre))
    }


    # Read in the viromics table
    viruses.pcl <- pcl.read(file.path(HMP2_data, "mvx", "HMP2.Virome.VirMAP.pcl.tsv"), metadata.rows=" ")

    # Filter samples with low reads and fix metadata
    viruses.pcl <- viruses.pcl %>%
        pcl.filter.s(keep=viruses.pcl$meta$reads_viral >= 10) %>%
        # Makes the virus names look more like metaphlan
        pcl.map.fnames(chartr(";, ", "|__", gsub("(s__.*?);[^;_]*$", "\\1",
            gsub("(super)?(\\w)\\w*=([^;]+)", "\\2__\\3", gsub(";taxId=\\d+", "", Name))))) %>%
        fix_metadata

    # Cleanup
    rm(bizarre)
}


