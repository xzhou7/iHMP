
library(dplyr)
library(ggplot2)

source("./common/merge_metadata.r")
source("./common/theme_nature.r")

## ========================================================================= ##
##                                                                 Load data ##

# Load the overall sample metadata file
metadata <- hmp2_sample_metadata

### Merge in noproducts
aliquot_data <- read.delim(file.path(HMP2_data, "aliquots", "IBDMDB Aliquot Tracking.tsv"))
# Calculate week_nums
aliquot_data$Actual.Date.of.Receipt <- as.Date(aliquot_data$Actual.Date.of.Receipt, format="%m-%d-%y")
aliquot_data <- aliquot_data[!is.na(aliquot_data$Actual.Date.of.Receipt),]
subject_first_stool <- aliquot_data$Actual.Date.of.Receipt[aliquot_data$Collection..==1]
names(subject_first_stool) <- aliquot_data$Subject[aliquot_data$Collection..==1]
any(!(aliquot_data$Subject %in% names(subject_first_stool)))
# [1] FALSE  # Check if any subjects in the main aliquot list don't have a first stool according collection #1
aliquot_data$first_stool_date <- subject_first_stool[match(aliquot_data$Subject, names(subject_first_stool))]
aliquot_data$week_num <- floor((aliquot_data$Actual.Date.of.Receipt-aliquot_data$first_stool_date)/7)
# Remove samples with all-"missing"/"destroyed"/"N/A"/etc..
has_aliquot <- as.matrix(aliquot_data[,5:ncol(aliquot_data)])
has_aliquot <- rowSums(matrix(grepl("^SM\\-", has_aliquot), ncol=ncol(has_aliquot))) > 0
aliquot_data <- aliquot_data[has_aliquot,] # Remove rows with no data
# Regenerate site-sub-coll
aliquot_data$site_sub_coll <- NA
aliquot_data$Participant.ID <- NA
site_abbrev <- c("MGH Peds" = "P", "MGH" = "M", "Emory" = "E", "CSMC" = "C", "CCHMC" = "H")
for (i in seq_len(nrow(aliquot_data))) {
    aliquot_data$Participant.ID[i] <- sprintf("%s%d",
        site_abbrev[as.character(aliquot_data$Sheet.Name[i])], aliquot_data$Subject[i])
    aliquot_data$site_sub_coll[i] <- sprintf("%sC%d",
        aliquot_data$Participant.ID[i], aliquot_data$Collection..[i])
}
# Add noproducts to the metadata table for site_sub_colls that don't match
noproduct_aliquots <- aliquot_data[!(aliquot_data$site_sub_coll %in% as.character(metadata$site_sub_coll)),
                                   c("site_sub_coll", "week_num", "Participant.ID")]
noproduct_aliquots$data_type <- "noproduct"
noproduct_aliquots$diagnosis <- metadata$diagnosis[match(noproduct_aliquots$Participant.ID, metadata$Participant.ID)]
metadata <- as.data.frame(data.table::rbindlist(list(metadata, noproduct_aliquots), fill=T, use.names=T))

### Merge in fecalcal
fc_data <- read.table(file.path(HMP2_data, "fc", "fecalcal_studytrax_measures.csv"), sep=",", header=T)
sum(!(fc_data$site_sub_coll %in% metadata$site_sub_coll))
# [1] 0  # Check how many fecalcals we can't match by site/sub/coll
metadata$fecalcal <- fc_data$fecalcal_ng_ml[match(metadata$site_sub_coll, fc_data$site_sub_coll)]

# Filter out subjects with only noproduct
keep_subj <- names(which(sapply(split(metadata, metadata$Participant.ID), function(df)
    any(df$data_type != "noproduct"))))
metadata <- metadata[metadata$Participant.ID %in% keep_subj,]

# Dataframe for plotting
samples.df <- data.frame(ssc=unique(metadata$site_sub_coll))
# filter to non-biopsy/blood
samples.df <- samples.df[grepl("^\\w\\d\\d\\d\\dC\\d\\d?$", samples.df$ssc),,drop=F]
get_minimal_NA <- function(field) {
    metadata.nna <- metadata[!is.na(metadata[,field]),]
    return (metadata.nna[match(samples.df$ssc, metadata.nna$site_sub_coll), field])
}
samples.df$week_num <- get_minimal_NA("week_num")
samples.df$hbi <- get_minimal_NA("hbi")
if (is.factor(samples.df$hbi))
    samples.df$hbi <- as.numeric(levels(samples.df$hbi))[as.numeric(samples.df$hbi)]
samples.df$hbi[samples.df$hbi>20] <- NA
samples.df$sccai <- get_minimal_NA("sccai")
samples.df$fecalcal <- get_minimal_NA("fecalcal")
samples.df$diagnosis <- get_minimal_NA("diagnosis")
samples.df$diagnosis[samples.df$diagnosis=="Other"] <- "nonIBD"
samples.df$diagnosis <- factor(samples.df$diagnosis, levels=c("CD", "UC", "nonIBD"))
samples.df$consent_age <- get_minimal_NA("consent_age")
samples.df$sex <- get_minimal_NA("sex")
samples.df$race <- get_minimal_NA("race")
samples.df$site_name <- get_minimal_NA("site_name")
samples.df$baseline_montreal_location <- get_minimal_NA("baseline_montreal_location")
samples.df$participant <- factor(gsub("^\\w(\\d\\d\\d\\d)C.*$", "\\1", samples.df$ssc))

append_if <- function(base, what, cond, sep=" ") {
    app <- base
    for (i in which(cond)) {
        if (app[i] == "") {
            app[i] <- what
        } else {
            app[i] <- sprintf("%s%s%s", app[i], sep, what)
        }
    }
    return (app)
}

# All data type combinations
samples.df$datatypes <- rep("", nrow(samples.df)) %>%
    append_if("MGX", samples.df$ssc %in% unique(metadata$site_sub_coll[(metadata$data_type == "metagenomics") | (metadata$data_type == "metatranscriptomics")])) %>%
    append_if("MTX", samples.df$ssc %in% unique(metadata$site_sub_coll[ metadata$data_type == "metatranscriptomics"])) %>%
    append_if( "VX", samples.df$ssc %in% unique(metadata$site_sub_coll[ metadata$data_type == "viromics"])) %>%
    append_if("MBX", samples.df$ssc %in% unique(metadata$site_sub_coll[ metadata$data_type == "metabolomics"])) %>%
    append_if("MPX", samples.df$ssc %in% unique(metadata$site_sub_coll[ metadata$data_type == "proteomics"])) %>%
    gsub("^$", "None", .) %>%
    factor

# Data type aesthetics split by MTX/MBX/MPX and MGX
samples.df$datatypes_tbp <- rep("", nrow(samples.df)) %>%
    append_if("MTX", samples.df$ssc %in% unique(metadata$site_sub_coll[metadata$data_type == "metatranscriptomics"])) %>%
    append_if("MBX", samples.df$ssc %in% unique(metadata$site_sub_coll[metadata$data_type == "metabolomics"])) %>%
    append_if("MPX", samples.df$ssc %in% unique(metadata$site_sub_coll[metadata$data_type == "proteomics"])) %>%
    gsub("^$", "None", .) %>%
    factor
samples.df$hasmgx <- samples.df$ssc %in% unique(metadata$site_sub_coll[metadata$data_type == "metagenomics"])

# Data type aesthetics - only MGX + MTX
samples.df$datatypes_mt <- rep("", nrow(samples.df)) %>%
    append_if("MGX", samples.df$ssc %in% unique(metadata$site_sub_coll[metadata$data_type == "metagenomics"])) %>%
    append_if("MTX", samples.df$ssc %in% unique(metadata$site_sub_coll[metadata$data_type == "metatranscriptomics"])) %>%
    gsub("^$", "None", .) %>%
    factor

# Summary of data type combinations
samples.df %>%
    group_by(datatypes) %>%
    summarise(count = length(datatypes))

# Reprioritize participants by number of data types generated for them
participant_weight <- samples.df %>%
    group_by(participant) %>%
    summarise(weight = length(participant))

samples.df$participant_wt <- factor(samples.df$participant,
    levels=as.character(participant_weight$participant)[order(participant_weight$weight)])
samples.df$participant_samps <- participant_weight$weight[as.numeric(samples.df$participant)]


# Make a dataframe with participant info
participants.df <- data.frame(id=unique(samples.df$participant))
get_minimal_NA_participant <- function(field) {
    samples.df.nna <- samples.df[!is.na(samples.df[,field]),]
    return (samples.df.nna[match(participants.df$id, samples.df.nna$participant), field])
}
participants.df$consent_age <- get_minimal_NA_participant("consent_age")
participants.df$diagnosis <- get_minimal_NA_participant("diagnosis")
participants.df$sex <- get_minimal_NA_participant("sex")
participants.df$race <- get_minimal_NA_participant("race")
participants.df$site_name <- get_minimal_NA_participant("site_name")
participants.df$baseline_montreal_location <- get_minimal_NA_participant("baseline_montreal_location")
participants.df$nsamples <- sapply(split(samples.df$participant, samples.df$participant), length)[participants.df$id]

## ========================================================================= ##
##          Giant set of tables for different breakdowns of the participants ##

sink("./overview/participant_breakdowns.txt")
cat("Participant metadata breakdown\nAll participants have min 1 sample with min 1 datatype\n\n")
cat("total_samples = Number of samples for which min 1 datatype has been generated.\n\n")
dump_breakdown <- function(groups) {
    if (length(groups) > 1) {
        cat(sprintf(" #### %s distribution by %s\n", groups[length(groups)], do.call(paste, c(as.list(groups[-length(groups)]), list(sep=", ")))))
    } else if (length(groups) > 0) {
        cat(sprintf(" #### %s distribution\n", groups))
    }
    participants.df %>%
        group_by_at(vars(groups)) %>%
        summarise(count=length(id),
                  total_samples=sum(nsamples)) %>%
        print(n=nrow(.), width=Inf)
    cat("\n")
}

breakdown_by <- c("site_name", "race", "sex", "diagnosis", "baseline_montreal_location")
cat(sprintf("Breakdown order: %s\n\n", do.call(paste, c(as.list(breakdown_by), list(sep=", ")))))

cat(" #### Summary\n")
participants.df %>%
    summarise(count=length(id),
              total_samples=sum(nsamples)) %>%
    print(width=Inf)
cat("\n")

for (i in seq_along(breakdown_by)) {
    combs <- combn(seq_along(breakdown_by), i)
    for (j in seq_along(combs[1,])) {
        dump_breakdown(breakdown_by[combs[,j]])
    }
}
sink()


## ========================================================================= ##
##                                                    Samples over time plot ##

hasmgxshape <- c("TRUE" = 21, "FALSE" = 25)
tbp_color <- c(
    "MTX"         = "forestgreen",
    "MBX"         = "mediumblue",
    "MPX"         = "red",
    "MTX MBX"     = "deepskyblue1",
    "MTX MPX"     = "goldenrod1",
    "MBX MPX"     = "darkorchid1",
    "MTX MBX MPX" = "white",
    "None"        = "gray50"
)
ggp <- ggplot(data=samples.df) +
    geom_point(aes(x=week_num, y=participant_wt, fill=datatypes_tbp, shape=hasmgx), size=2) +
    scale_fill_manual(values=tbp_color) +
    scale_shape_manual(values=hasmgxshape) +
    facet_grid(diagnosis ~ ., scales="free", space="free") +
    theme(axis.text.y=element_text(size=8)) +
    guides(fill=guide_legend(title="Data Types", override.aes=list(size=4, shape=21)),
           shape=guide_legend(title="Has MGX?", override.aes=list(size=3))) +
    ylab("Participant") + xlab("Week")
pdf("./overview/sample_datatypes.pdf", 8, 14)
print(ggp)
dev.off()

ggp <- ggplot(data=samples.df) +
    geom_point(aes(x=week_num, y=participant_wt, fill=datatypes_tbp, shape=hasmgx, size=hasmgx), stroke=0.25) +
    scale_fill_manual(values=tbp_color) +
    scale_shape_manual(values=hasmgxshape) +
    scale_size_manual(values=c("TRUE"=1.2, "FALSE"=0.7)) +
    facet_grid(diagnosis ~ ., scales="free", space="free") +
    guides(fill="none", shape="none", size="none") +
    ylab("Participant") + xlab("Week") + theme_nature() +
    theme(axis.text.y=element_text(size=4.8))
pdf("./overview/sample_datatypes_edfig.pdf", 2.566, 8.3)
print(ggp)
dev.off()

## ========================================================================= ##
##                                               MGX|MTX only over time plot ##

# Print-ready size

gtn_color <- c(
    "MGX"     = "gray70",
    "MGX MTX" = "gray10"
)
min_datatypes <- 5
library(cowplot)
ggp <- ggplot(data=samples.df[(samples.df$participant_samps >= min_datatypes) & samples.df$hasmgx,]) +
    theme_cowplot() +
    geom_vline(data=data.frame(x=c(0, 2, 4, 8, 16, 24, 36)), aes(xintercept=x), color="gray50", size=0.25) +
    geom_point(aes(x=week_num, y=participant_wt, color=datatypes_mt), size=0.3, stroke=0.25) +
    scale_color_manual(values=gtn_color) +
    scale_x_continuous(limits = c(-1, max(samples.df$week_num)+1), expand = c(0, 0)) +
    facet_grid(diagnosis ~ ., scales="free", space="free") +
    ylab("Participant") + xlab("Week") +
    theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.ticks.x=element_blank(),
          axis.text.x=element_text(size=5, margin=margin(1, 1, 1, 1)),
          strip.text.y=element_text(size=5, margin=margin(1, 1, 1, 1)),
          axis.title.x=element_text(size=6, margin=margin(1, 1, 1, 1)),
          axis.title.y=element_text(size=6, margin=margin(1, 1, 1, 1)),
          panel.spacing=unit(0.1, "lines"),
          axis.line.x=element_line(size=0.25),
          axis.line.y=element_line(size=0.25),
          plot.margin=margin(2, 2, 2, 2)) +
    guides(color="none", fill="none")
pdf("./overview/samples_mtxmgx_over_time.pdf", 1.2, 1.96)
print(ggp)
dev.off()


## ========================================================================= ##
##                                                   Disease Severity Scores ##

samples.df.cd <- samples.df[samples.df$diagnosis=="CD",]
dsev_theme <- list(
    scale_fill_manual(values=tbp_color),
    scale_shape_manual(values=hasmgxshape),
    facet_wrap(~participant, ncol=6),
    theme(strip.text.x = element_text(size=8, margin = margin(.1, 0, .1, 0, "cm"))),
    guides(fill=guide_legend(title="Data Types", override.aes=list(size=4, shape=21)),
           shape=guide_legend(title="Has MGX?", override.aes=list(size=3))),
    xlab("Week")
)

ggp <- ggplot(data=samples.df.cd[!is.na(samples.df.cd$hbi),], aes(x=week_num, y=hbi)) +
    geom_line(aes(group=participant)) +
    geom_point(aes(fill=datatypes_tbp, shape=hasmgx), size=2) +
    ylab("HBI") + dsev_theme
pdf("./overview/timeseries_cd_hbi.pdf", 11, 15)
print(ggp)
dev.off()

ggp <- ggplot(data=samples.df.cd[!is.na(samples.df.cd$fecalcal),], aes(x=week_num, y=fecalcal)) +
    geom_line(aes(group=participant)) +
    geom_point(aes(fill=datatypes_tbp, shape=hasmgx), size=2) +
    ylab("FecalCal") + dsev_theme
pdf("./overview/timeseries_cd_fecalcal.pdf", 11, 11)
print(ggp)
dev.off()


samples.df.uc <- samples.df[samples.df$diagnosis=="UC",]
ggp <- ggplot(data=samples.df.uc[!is.na(samples.df.uc$sccai),], aes(x=week_num, y=sccai)) +
    geom_line(aes(group=participant)) +
    geom_point(aes(fill=datatypes_tbp, shape=hasmgx), size=2) +
    ylab("SCCAI") + dsev_theme
pdf("./overview/timeseries_uc_sccai.pdf", 11, 9)
print(ggp)
dev.off()

ggp <- ggplot(data=samples.df.uc[!is.na(samples.df.uc$fecalcal),], aes(x=week_num, y=fecalcal)) +
    geom_line(aes(group=participant)) +
    geom_point(aes(fill=datatypes_tbp, shape=hasmgx), size=2) +
    ylab("FecalCal") + dsev_theme
pdf("./overview/timeseries_uc_fecalcal.pdf", 11, 7)
print(ggp)
dev.off()


## ========================================================================= ##
##                                            Distribution of Age at Consent ##

sink("./overview/age_distribution_na.txt")
cat(sprintf("These %d participants have NA consent_age:\n", sum(is.na(participants.df$consent_age))))
as.character(participants.df$id[is.na(participants.df$consent_age)])
sink()

ggp <- ggplot(data=participants.df) + theme_gray() +
    geom_histogram(aes(x=consent_age), color="black", fill="lightskyblue") +
    xlab("Age at Consent") + ylab("Count")
pdf("./overview/age_distribution_overall.pdf", 3.5, 3.5)
print(ggp)
dev.off()

pdf("./overview/age_distribution_diagnosis.pdf", 9, 3.5)
print(ggp + facet_grid(. ~ diagnosis))
dev.off()

pdf("./overview/age_distribution_site.pdf", 10, 3)
print(ggp + facet_grid(. ~ site_name))
dev.off()

pdf("./overview/age_distribution_diagnosis_site.pdf", 9, 2.5*3.5)
print(ggp + facet_grid(site_name ~ diagnosis))
dev.off()



