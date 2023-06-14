
### RUN taxonomy_overview.r first to get bugs.pcl.unfiltered ###
source("./overview/src/taxonomy_overview.r")

# Species-level BC dissimilarity matrix
library(labdsv)
bc <- dsvdis(pcl.only(bugs.pcl.unfiltered, rank="s")$x, "bray/curtis")

# Histogram of distances
pdf("./overview/ecological_qc_bc_distribution.pdf", 6, 5)
hist(bc, xlab="Bray-Curtis", main="Distribution of all Bray-Curtis dissimilarities")
dev.off()

# Median dissimilarities
mdd <- apply(as.matrix(bc), 1, median)

# Upper inner fence = upper quartile + 1.5*IQR
uif <- quantile(mdd, 0.75) + 1.5*IQR(mdd)
uif
# 75%
# 0.9880927

# How many and which samples get dropped?
fail_qc <- mdd>=uif
sum(fail_qc)
# [1] 38

# Dump out the table of failed sample IDs
write.table(file="./overview/hmp1-ii_ecological_qc_fail.txt", quote=F, sep="\t", row.names=F,
            data.frame(externalid=names(mdd)[fail_qc],
                       site_sub_coll=as.character(bugs.pcl.unfiltered$meta$site_sub_coll[fail_qc])))


# Read depths of failed samples?
pdf("./overview/hmp1-ii_ecological_qc_failed_read_depths.pdf", 6, 5)
hist(bugs.pcl.unfiltered$meta$filtered_reads[fail_qc],
     xlab="Filtered Read Count",
     main="Read depth of samples that failed ecological QC")
dev.off()


low_readdepth <- bugs.pcl.unfiltered$meta$filtered_reads < 1e6
sum(low_readdepth & fail_qc, na.rm=T)
# [1] 6
sum(is.na(low_readdepth) & fail_qc)
# [1] 2


# Per-person ecological QC
subjects <- levels(bugs.pcl.unfiltered$meta$subject)
fail_qc <- rep(F, bugs.pcl.unfiltered$ns)
names(fail_qc) <- rownames(bugs.pcl.unfiltered$x)
for (sub in subjects) {
    subj.only <- pcl.filter.s(bugs.pcl.unfiltered, keep=bugs.pcl.unfiltered$meta$subject == sub)
    bc <- dsvdis(pcl.only(subj.only, rank="s")$x, "bray/curtis")

    # Ecological QC (same as above but per-subject)
    mdd <- apply(as.matrix(bc), 1, median)
    uif <- quantile(mdd, 0.75) + 1.5*IQR(mdd)
    fail_qc[rownames(subj.only$x)] <- mdd>=uif
}

sum(fail_qc)
# [1] 130

# Dump out the table of failed sample IDs
write.table(file="./overview/persubject_ecological_qc_fail.txt", quote=F, sep="\t", row.names=F,
            data.frame(externalid=names(fail_qc)[fail_qc],
                       site_sub_coll=as.character(bugs.pcl.unfiltered$meta$site_sub_coll[fail_qc])))




