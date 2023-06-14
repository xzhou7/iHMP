
# Load protein tables
#   proteins.pcl: Raw counts table, taken from the 1pep-pro <.05 FDR table
#   proteins.ecs.pcl: EC-level rollup of the protein counts
#   proteins.kos.pcl: KO-level rollup of the protein counts

source("./env_config.r")
source("./common/pcl_utils.r")
source("./common/fix_metadata.r")

library(dplyr)

if (!exists("proteins.pcl")) {
    # Load in giant table
    proteins.pcl <- pcl.read(file.path(HMP2_data, "mpx", "iHMP_all_proteomics_1pep-pro_1p_FDR_5ppm.pcl.tsv.gz"), metadata.rows=" ")

    # Filter out samples with no metadata
    proteins.pcl <- proteins.pcl %>%
        pcl.filter.s(!is.na(diagnosis))

    # NA = 0
    proteins.pcl$x[is.na(proteins.pcl$x)] <- 0

    # Apply common metadata fixes and normalize the data
    proteins.pcl <- fix_metadata(proteins.pcl)


    # Load up more readable versions
    proteins.ecs.pcl <- pcl.read(file.path(HMP2_data, "mpx", "hmp2_MPX_not_normalized_ecs.names.pcl.tsv"), metadata.rows=" ")
    proteins.kos.pcl <- pcl.read(file.path(HMP2_data, "mpx", "hmp2_MPX_not_normalized_kos.names.pcl.tsv"), metadata.rows=" ")

    # Apply common metadata fixes
    proteins.ecs.pcl <- fix_metadata(proteins.ecs.pcl)
    proteins.kos.pcl <- fix_metadata(proteins.kos.pcl)
}






