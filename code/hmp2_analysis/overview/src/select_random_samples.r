
## Script to select some random samples which evenly tile the bug PCoA

library(dplyr)
source("./common/load_bugs.r")


bug.ord <- pcl.pcoa(bugs.pcl %>% pcl.only(rank="s"))

bugs.pcl$meta$picked <- F
avail <- seq_len(nrow(bugs.pcl$meta))
rad <- 0.22
while (length(avail) > 0) {
    samp <- sample(avail, 1)
    bugs.pcl$meta$picked[samp] <- T
    d2 <- (bug.ord$points[avail,1]-bug.ord$points[samp,1])^2 + (bug.ord$points[avail,2]-bug.ord$points[samp,2])^2
    avail <- avail[d2>rad^2]
}

sum(bugs.pcl$meta$picked)

pcl.ordplot(bugs.pcl, bug.ord, colour="picked", colour_override = c("TRUE"="red", "FALSE"="black"), size_abs=2)

rownames(bugs.pcl$meta)[bugs.pcl$meta$picked]

