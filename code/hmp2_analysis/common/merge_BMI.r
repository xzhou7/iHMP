
source("./common/merge_metadata.r")

merge_BMI <- function(pcl) {
    # Merge in BMI based on subject label

    mtbmi <- hmp2_sample_metadata$BMI
    names(mtbmi) <- hmp2_sample_metadata$Participant.ID

    mtbmi <- mtbmi[!is.na(mtbmi)]

    pcl$meta$BMI <- mtbmi[match(pcl$meta$subject, names(mtbmi))]

    return (pcl)
}

