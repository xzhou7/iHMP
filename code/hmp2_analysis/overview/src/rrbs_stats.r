
source("./common/merge_metadata.r")

rrbs_stats <- read.delim(file.path(HMP2_data, "rrbs", "RRBS Picard Report for 504 Samples 17Aug2017.txt"),
                         sep="\t", stringsAsFactors = F)

summary(rrbs_stats)

colwise(sd)(rrbs_stats)

sum(rrbs_stats$Total.Reads > 10e6)
# [1] 448
mean(rrbs_stats$Total.Reads > 10e6)
# [1] 0.8888889
plot(ecdf(rrbs_stats$Total.Reads))

# Alignment rates
plot(rrbs_stats$Total.Reads, rrbs_stats$Reads.Aligned..RRBS.)

summary(rrbs_stats$Reads.Aligned..RRBS. / rrbs_stats$Total.Reads)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.7064  0.9465  0.9588  0.9516  0.9634  0.9744

