
source("./env_config.r")
source("./common/merge_metadata.r")
source("./common/disease_colors.r")
source("./common/theme_nature.r")

mds <- read.table(file.path(HMP2_data, "exome", "PCs", "HMP.and1kG.MDS.mds"), header=T)
pop <- read.table(file.path(HMP2_data, "exome", "PCs", "kgp.pop"), header=T)

mds <- merge(mds, pop, by="IID", all=T)
mds$hmp2 <- is.na(mds$POP)
mds <- mds[order(mds$hmp2),]
mds$diagnosis <- hmp2_sample_metadata$diagnosis[match(mds$IID, hmp2_sample_metadata$External.ID)]
mds$ANC <- as.character(mds$ANC)
mds$ANC[mds$hmp2] <- "HMP2"

library(ggplot2)
library(cowplot)
pdf("./overview/genetics_ordination.pdf", 1.29921 * (20.523 / 29.612), 0.826772)
ggp <- ggplot(data=mds, aes(x=C1, y=C2)) +
    geom_point(aes(color=ANC, shape=hmp2, size=hmp2, alpha=hmp2), fill="gray90", stroke=0.25) +
    scale_shape_manual(values=c("TRUE"=21, "FALSE"=3)) +
    scale_size_manual(values=c("TRUE"=0.65, "FALSE"=0.4)) +
    scale_alpha_manual(values=c("TRUE"=1, "FALSE"=0.25)) +
    #scale_fill_manual(values=c(hmp2_disease_colors, c("NA"="black"))) +
    scale_color_manual(values=c(EUR="#e41a1c", EAS="#ff7f00", AMR="#377eb8", SAS="#984ea3", AFR="#4daf4a", HMP2="black")) +
    theme_nature() +
    theme(axis.ticks=element_blank(), axis.text=element_blank()) +
    xlab(NULL) + ylab(NULL)
print(ggp + guides(color="none", shape="none", size="none", alpha="none"))
print(ggp)
dev.off()


