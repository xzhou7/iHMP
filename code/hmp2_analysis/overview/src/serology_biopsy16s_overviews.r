

source("./common/disease_colors.r")

library(viridis)
make_heatmap <- function(pcl, logspace=T, ...) {
    pcl.heatmap(pcl,
                meta=c("diagnosis"), annotation_colors=list(diagnosis=hmp2_disease_colors),
                show_rownames=T, color=viridis(100), logspace=logspace, minx=1e-4, zerospecial=NA, ...)
}


# --------------------------------------------
# Serology

source("./common/load_serology.r")

# Merge in active disease data
source("./common/load_bugs.r")
source("./common/merge_metadata.r")
source("./common/disease_activity.r")

#serology.eu.pcl <- merge_metadata(serology.eu.pcl, "week_num", "serology")
#serology.pos.pcl <- merge_metadata(serology.pos.pcl, "week_num", "serology")

serology.eu.pcl <- merge_disease_activity(serology.eu.pcl, lenience=2)
serology.pos.pcl <- merge_disease_activity(serology.pos.pcl, lenience=2)


library(reshape2)
library(ggplot2)
library(cowplot)
library(ggsignif)
library(ggpubr)
library(metap)
serology.eu.pcl$meta$active_diagnosis <- sprintf("%s %s", ifelse(serology.eu.pcl$meta$active, "Active", "Inactive"), serology.eu.pcl$meta$diagnosis)
serology.eu.pcl$meta$active_diagnosis[serology.eu.pcl$meta$diagnosis=="nonIBD"] <- "nonIBD"
serology.eu.pcl$meta$active_diagnosis <- factor(serology.eu.pcl$meta$active_diagnosis,
    levels=c("nonIBD", "Inactive UC", "Inactive CD", "Active UC", "Active CD"))
df <- melt(cbind(serology.eu.pcl$meta[,c("diagnosis", "active_diagnosis", "active"),drop=F], serology.eu.pcl$x))
df <- df[!is.na(df$diagnosis),,drop=F]
pdf("./overview/serology_boxplots.pdf", 10, 6)
common_theme <- list(
    scale_fill_manual(values=hmp2_disease_colors, name="Diagnosis"),
    facet_wrap(~variable, drop=T, scales="free", ncol=5),
    theme_cowplot(), theme(strip.background = element_blank(), strip.text.x = element_blank()),
    xlab(NULL)
)
common_theme_active <- list(
    scale_fill_manual(values=hmp2_active_disease_colors, name="Diagnosis"),
    facet_wrap(~variable, drop=T, scales="free", ncol=5),
    theme_cowplot(), theme(strip.background = element_blank(), strip.text.x = element_blank()),
    xlab(NULL)
)
ggplot(data=df, aes(x=variable, y=value, fill=diagnosis)) +
    geom_boxplot() + geom_boxplot_n(offset=4) + common_theme +
    ggtitle("Serology by Diagnosis") + ylab("Titer")
ggplot(data=df, aes(x=variable, y=value, fill=diagnosis)) +
    geom_violin(scale="width") + common_theme +
    ggtitle("Serology by Diagnosis") + ylab("Titer")
ggplot(data=df, aes(x=variable, y=value, fill=diagnosis)) +
    geom_boxplot() + geom_boxplot_n(offset=4) + common_theme +
    scale_y_log10() +
    ggtitle("Serology by Diagnosis") + ylab("Titer (log)")
ggplot(data=df, aes(x=variable, y=value, fill=diagnosis)) +
    geom_violin(scale="width") + common_theme +
    scale_y_log10() +
    ggtitle("Serology by Diagnosis") + ylab("Titer (log)")

ggplot(data=df, aes(x=variable, y=value)) +
    geom_boxplot() + geom_boxplot_n(offset=4) + common_theme +
    ggtitle("Serology overall") + ylab("Titer")
ggplot(data=df, aes(x=variable, y=value)) +
    geom_violin(scale="width") + common_theme +
    ggtitle("Serology overall") + ylab("Titer")
ggplot(data=df, aes(x=variable, y=value)) +
    geom_boxplot() + scale_y_log10() + geom_boxplot_n(offset=0.2) + common_theme +
    ggtitle("Serology overall (log)") + ylab("Titer (log)")
ggplot(data=df, aes(x=variable, y=value)) +
    geom_violin(scale="width") + scale_y_log10() + common_theme +
    ggtitle("Serology overall (log)") + ylab("Titer (log)")

ggplot(data=df[!is.na(df$active),,drop=F], aes(x=variable, y=value, fill=active_diagnosis)) +
    geom_boxplot() + geom_boxplot_n(offset=4) +
    common_theme_active +
    ggtitle("Serology vs Disease Activity") + ylab("Titer")
ggplot(data=df[!is.na(df$active),,drop=F], aes(x=variable, y=value, fill=active_diagnosis)) +
    geom_violin(scale="width") + common_theme_active +
    ggtitle("Serology vs Disease Activity") + ylab("Titer")
ggplot(data=df[!is.na(df$active),,drop=F], aes(x=variable, y=value, fill=active_diagnosis)) +
    geom_boxplot() + geom_boxplot_n(offset=0.2) +
    common_theme_active +
    scale_y_log10() +
    ggtitle("Serology vs Disease Activity (log)") + ylab("Titer (log)")
ggplot(data=df[!is.na(df$active),,drop=F], aes(x=variable, y=value, fill=active_diagnosis)) +
    geom_violin(scale="width") +
    common_theme_active +
    scale_y_log10() +
    ggtitle("Serology vs Disease Activity (log)") + ylab("Titer (log)")
dev.off()

### Figure for paper
paper_theme <- list(
    scale_fill_manual(values=hmp2_active_disease_colors_fill, name="Diagnosis"),
    scale_color_manual(values=hmp2_active_disease_colors_outline, name="Diagnosis"),
    facet_wrap(~variable, drop=T, scales="free", ncol=2),
    theme_nature(), theme(strip.background = element_blank(), strip.text.x = element_blank()),
    xlab(NULL)
)
levels(df$variable) <- gsub(" EU$", "", levels(df$variable))
levels(df$variable)[levels(df$variable)=="IgA ASCA"] <- "ASCA (IgA)"
levels(df$variable)[levels(df$variable)=="IgG ASCA"] <- "ASCA (IgG)"
pdf("./overview/serology_panel_fig2.pdf", 1.3, 2.15)
ggplot(data=df[!is.na(df$active),,drop=F], aes(x=variable, y=value, color=active_diagnosis, fill=active_diagnosis)) +
    geom_boxplot(lwd=0.25, outlier.size=0.35, outlier.stroke=0.25) +
    geom_boxplot_n(offset = 3, size=1.5) +
    paper_theme + guides(fill="none", color="none") +
    ylab("Titer (EU/mL)")
ggplot(data=df[!is.na(df$active),,drop=F], aes(x=active_diagnosis, y=value, color=active_diagnosis)) +
    geom_point(size=0.35, position=position_jitter(width=0.1)) +
    paper_theme + guides(fill="none", color="none") +
    ylab("Titer (EU/mL)")
dev.off()


# Wilcoxon tests of differences between active/inactive
ks.test.p <- matrix(NA, 6, 2)
rownames(ks.test.p) <- c(levels(df$variable), "Combined_Fisher")
colnames(ks.test.p) <- c("UC", "CD")
wk.test.p <- ks.test.p
t.test.p <- ks.test.p
df$value <- df$value + 0.01*rnorm(length(df$value))
for (feat in levels(df$variable)) {
    ks.test.p[feat, "UC"] <- ks.test(
        df$value[(df$variable == feat) & (df$diagnosis == "UC") & (!is.na(df$active) & df$active)],
        df$value[(df$variable == feat) & (df$diagnosis == "UC") & (!is.na(df$active) &!df$active)])$p.value
    ks.test.p[feat, "CD"] <- ks.test(
        df$value[(df$variable == feat) & (df$diagnosis == "CD") & (!is.na(df$active) & df$active)],
        df$value[(df$variable == feat) & (df$diagnosis == "CD") & (!is.na(df$active) &!df$active)])$p.value
    wk.test.p[feat, "UC"] <- wilcox.test(
        df$value[(df$variable == feat) & (df$diagnosis == "UC") & (!is.na(df$active) & df$active)],
        df$value[(df$variable == feat) & (df$diagnosis == "UC") & (!is.na(df$active) &!df$active)])$p.value
    wk.test.p[feat, "CD"] <- wilcox.test(
        df$value[(df$variable == feat) & (df$diagnosis == "CD") & (!is.na(df$active) & df$active)],
        df$value[(df$variable == feat) & (df$diagnosis == "CD") & (!is.na(df$active) &!df$active)])$p.value
    t.test.p[feat, "UC"] <- t.test(
        df$value[(df$variable == feat) & (df$diagnosis == "UC") & (!is.na(df$active) & df$active)],
        df$value[(df$variable == feat) & (df$diagnosis == "UC") & (!is.na(df$active) &!df$active)])$p.value
    t.test.p[feat, "CD"] <- t.test(
        df$value[(df$variable == feat) & (df$diagnosis == "CD") & (!is.na(df$active) & df$active)],
        df$value[(df$variable == feat) & (df$diagnosis == "CD") & (!is.na(df$active) &!df$active)])$p.value
}
combine_pvalues <- function(m)apply(m, 2, function(x)sumlog(x)$p)
ks.test.p["Combined_Fisher",] <- combine_pvalues(ks.test.p[1:5,])
wk.test.p["Combined_Fisher",] <- combine_pvalues(wk.test.p[1:5,])
t.test.p["Combined_Fisher",] <- combine_pvalues(t.test.p[1:5,])

write.table(ks.test.p, "./overview/serology_activity_ks_test.txt", quote=F)
write.table(wk.test.p, "./overview/serology_activity_wilcox_test.txt", quote=F)
write.table(t.test.p, "./overview/serology_activity_t_test.txt", quote=F)

pdf("./overview/serology_ecdf_activity_lenience2.pdf", 10, 6)
ggplot(data=df, aes(x=value, color=active)) +
    stat_ecdf() + facet_grid(diagnosis ~ variable, scale="free") +
    scale_color_discrete(name="Active") +
    xlab("Titer") + ylab(NULL)
dev.off()
pdf("./overview/serology_boxpot_split_activity_lenience2_wilcox.pdf", 5, 12)
ggplot(data=df[!is.na(df$active) & df$diagnosis!="nonIBD",,drop=F], aes(y=value, x=active, fill=active)) +
    geom_boxplot() + geom_boxplot_n(offset=4) +
    facet_grid(variable ~ diagnosis, scale="free") +
    stat_compare_means(method="wilcox.test") +
    scale_fill_discrete(guide="none") +
    xlab("Active?") + ylab(NULL)
dev.off()


serology.pca.ord <- pcl.pcoa(serology.eu.pcl, D=dist(serology.eu.pcl$x, method="euclidean")^2)
pdf("./overview/serology_pca.pdf", 5.5, 5)
pcl.ordplot(serology.eu.pcl, serology.pca.ord,
            colour="diagnosis", colour_override=hmp2_disease_colors)
dev.off()

serology.pcoa.ord <- pcl.pcoa(serology.eu.pcl %>% pcl.normalize)
pdf("./overview/serology_pcoa_bc.pdf", 5.5, 5)
pcl.ordplot(serology.eu.pcl, serology.pcoa.ord,
            colour="diagnosis", colour_override=hmp2_disease_colors)
dev.off()

serology.tsne.ord <- pcl.tsne(serology.eu.pcl %>% pcl.normalize)
pdf("./overview/serology_tsne_bc.pdf", 5.5, 5)
pcl.ordplot(serology.eu.pcl, serology.tsne.ord,
            colour="diagnosis", colour_override=hmp2_disease_colors)
dev.off()

pdf("./overview/serology_heatmap.pdf", 15, 5, onefile=F)
make_heatmap(serology.eu.pcl, logspace=F)
dev.off()


# Active disease plots
pdf("./overview/serology_pcoa_bc_coloractive.pdf", 5.5, 5)
pcl.ordplot(serology.eu.pcl, serology.pcoa.ord,
            colour="active", shape="diagnosis")
dev.off()
pdf("./overview/serology_pcoa_bc_shapeactive.pdf", 5.5, 5)
pcl.ordplot(serology.eu.pcl, serology.pcoa.ord,
            colour="diagnosis", shape="active", colour_override=hmp2_disease_colors)
dev.off()

pdf("./overview/serology_pos_heatmap.pdf", 15, 5, onefile=F)
pcl.heatmap(serology.pos.pcl,
            meta=c("diagnosis", "active"), annotation_colors=list(diagnosis=hmp2_disease_colors),
            show_rownames=T, color=viridis(100), logspace=F, minx=1e-4, zerospecial=NA)
dev.off()


# --------------------------------------------
# Biopsy 16S

source("./common/load_biopsies.r")

biopsy_16s.pcl$meta$biopsy_location_ilcolrec <- NA
biopsy_16s.pcl$meta$biopsy_location_ilcolrec[!is.na(biopsy_16s.pcl$meta$biopsy_location)] <- "Colon"
biopsy_16s.pcl$meta$biopsy_location_ilcolrec[biopsy_16s.pcl$meta$biopsy_location=="Ileum"] <- "Ileum"
biopsy_16s.pcl$meta$biopsy_location_ilcolrec[biopsy_16s.pcl$meta$biopsy_location=="Rectum"] <- "Rectum"

biopsy_16s.pcoa.ord <- pcl.pcoa(biopsy_16s.pcl)
pdf("./overview/biopsy_16s_pcoa_bc.pdf", 5.5, 5)
pcl.ordplot(biopsy_16s.pcl, biopsy_16s.pcoa.ord,
            colour="diagnosis", colour_override=hmp2_disease_colors)
pcl.ordplot(biopsy_16s.pcl, biopsy_16s.pcoa.ord,
            colour="diagnosis", colour_override=hmp2_disease_colors,
            shape="biopsy_location_ilcolrec", shape_override=c(Colon=22, Ileum=24, Rectum=21))
dev.off()

biopsy_16s.tsne.ord <- pcl.tsne(biopsy_16s.pcl)
pdf("./overview/biopsy_16s_tsne_bc.pdf", 5.5, 5)
pcl.ordplot(biopsy_16s.pcl, biopsy_16s.tsne.ord,
            colour="diagnosis", colour_override=hmp2_disease_colors)
dev.off()

pdf("./overview/biopsy_16s_heatmap.pdf", 22, 8, onefile=F)
make_heatmap(biopsy_16s.pcl %>% pcl.top.f(mean(x>0), n=20))
dev.off()





