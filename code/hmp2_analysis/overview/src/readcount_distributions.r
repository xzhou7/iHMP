
source("./common/merge_metadata.r")
source("./common/load_ecs.r")

library(ggplot2)

mgx_meta <- hmp2_sample_metadata[hmp2_sample_metadata$data_type=="metagenomics",]
summary(mgx_meta$reads_raw)
max(mgx_meta$reads_raw)/min(mgx_meta$reads_raw)
sd(mgx_meta$reads_raw)/mean(mgx_meta$reads_raw)

pdf("./read_mgx.pdf", 4, 3)
ggplot(data=mgx_meta) +
    geom_histogram(aes(x=reads_raw))
dev.off()


mtx_meta <- hmp2_sample_metadata[hmp2_sample_metadata$data_type=="metatranscriptomics",]
summary(mtx_meta$reads_raw)
max(mtx_meta$reads_raw)/min(mtx_meta$reads_raw)
sd(mtx_meta$reads_raw)/mean(mtx_meta$reads_raw)

pdf("./read_mtx.pdf", 4, 3)
ggplot(data=mtx_meta) +
    geom_histogram(aes(x=reads_raw))
dev.off()


hmp_bugs.pcl <- pcl.read(file.path(HMP2_data, "hmp1-ii", "hmp1-II_metaphlan2-mtd-qcd.tsv"))


max(bugs.pcl$meta$reads_filtered)/min(bugs.pcl$meta$reads_filtered)
max(ec.rna.unstrat.pcl$meta$reads_filtered)/min(ec.rna.unstrat.pcl$meta$reads_filtered)
pdf("./read_mgxmtxqc.pdf", 4, 3)
ggplot(data=ec.dna.unstrat.pcl$meta) +
    geom_histogram(aes(x=reads_filtered))
ggplot(data=ec.rna.unstrat.pcl$meta) +
    geom_histogram(aes(x=reads_filtered))
dev.off()


library(truncnorm)
MuLibSize = 10.04278
SDLibSize = 1.112657
viThreshold = 3
iLibSize = exp(truncnorm::rtruncnorm(n=20, mean = MuLibSize, sd = SDLibSize, b = MuLibSize + viThreshold * SDLibSize))
max(iLibSize) / min(iLibSize)

pdf("./read_sparsedossa.pdf", 4, 3)
ggplot(data=data.frame(reads=iLibSize)) +
    geom_histogram(aes(x=reads), bins=30)
dev.off()
