
# Load Serology tables
#   serology.pcl: Raw mixed count/boolean table
#   serology.eu.pcl: Filtered to the actual measurements
#   serology.pos.pcl: Filtered to positive/negative calls

if (!exists("serology.pos.pcl")) {
    source("./env_config.r")
    source("./common/pcl_utils.r")
    source("./common/fix_metadata.r")

    library(dplyr)

    # Read in raw table
    serology.pcl <- pcl.read(file.path(HMP2_data, "serology", "hmp2_serology_Compiled_ELISA_Data.pcl.tsv"), metadata.rows=" ") %>%
        fix_metadata

    # Filter to actual measurements only
    serology.eu.pcl <- serology.pcl %>%
        pcl.filter.f(grepl("EU", Name))

    # Filter to nonredundant Pos only
    serology.pos.pcl <- serology.pcl %>%
        pcl.filter.f(grepl("Pos", Name))
}


