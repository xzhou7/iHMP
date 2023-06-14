
source("./common/pcl_utils.r")
source("./env_config.r")

mgx.bug.ko <- pcl.read(file.path(HMP2_data, "mgx/kos_relab.pcl.tsv.gz"), metadata.rows = " ")
mtx.bug.ko <- pcl.read(file.path(HMP2_data, "mtx/kos_relab.pcl.tsv.gz"), metadata.rows = " ")

mgx.bug.ko.strat <- mgx.bug.ko %>%
    pcl.filter.f(keep=grepl("\\|", colnames(mgx.bug.ko$x)))
mgx.bug.ko.bug <- mgx.bug.ko.strat %>%
    pcl.group.f(colnames(pcl.nicenames(mgx.bug.ko.strat)$x), avg=F)

mtx.bug.ko.strat <- mtx.bug.ko %>%
    pcl.filter.f(keep=grepl("\\|", colnames(mtx.bug.ko$x)))
mtx.bug.ko.bug <- mtx.bug.ko.strat %>%
    pcl.group.f(colnames(pcl.nicenames(mtx.bug.ko.strat)$x), avg=F)

mgx.bug.ko.bug %>% pcl.top.f(mean(x), n=20) %>% pcl.heatmap
mtx.bug.ko.bug %>% pcl.top.f(mean(x), n=20) %>% pcl.heatmap

pcl.write(mgx.bug.ko.bug, file.path(HMP2_data, "mgx", "kos_relab_bugsummary.pcl.tsv"))
pcl.write(mtx.bug.ko.bug, file.path(HMP2_data, "mtx", "kos_relab_bugsummary.pcl.tsv"))

# these files need to have the blank line between metadata and data inserted manually
