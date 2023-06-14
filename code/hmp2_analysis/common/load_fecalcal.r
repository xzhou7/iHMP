
# Load fecalcal measurements
#   fecalcal.pcl: Fecalcal measurements

if (!exists("fecalcal.pcl")) {
    library(plyr)
    library(dplyr)

    source("./env_config.r")
    source("./common/pcl_utils.r")
    source("./common/fix_metadata.r")
    source("./common/merge_metadata.r")

    # Read the file
    fc_data <- read.table(file.path(HMP2_data, "fc", "fecalcal_studytrax_measures.csv"), sep=",", header=T)

    # Split into metadata and data
    fc_meta <- fc_data[,"site_sub_coll", drop=F]
    fc_data <- as.matrix(fc_data[,"fecalcal_ng_ml", drop=F])
    rownames(fc_meta) <- fc_meta$site_sub_coll
    rownames(fc_data) <- fc_meta$site_sub_coll

    # Make it a pcl structure and apply common metadata fixes
    fecalcal.pcl <- pcl.make(fc_data, meta=fc_meta) %>%
        fix_metadata %>%
        merge_metadata("week_num")

    # There are samples here that are not in the main metdata table, so merge diagnosis separately
    fecalcal.pcl$meta$diagnosis <- hmp2_sample_metadata$diagnosis[match(fecalcal.pcl$meta$subject, hmp2_sample_metadata$Participant.ID)]

    # Clean up
    rm(fc_data, fc_meta)
}





