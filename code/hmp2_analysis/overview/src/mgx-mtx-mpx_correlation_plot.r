
# load the data
source("./env_config.r")
df_Summary_Corrs <- readRDS(file.path(HMP2_root, "analysis", "mgx-mtx-mpx", "df_Summary_Corrs.rds"))

# reorder factors
df_Summary_Corrs$between <- factor(df_Summary_Corrs$between,
                                   levels=c("MGX_MTX", "MTX_MPX", "MGX_MPX"))
levels(df_Summary_Corrs$between) <- gsub("_", "-", levels(df_Summary_Corrs$between))

# make the plot
source("./common/theme_nature.r")
library(ggplot2)
ggp <- ggplot(df_Summary_Corrs, aes(x=Correlation, fill=between)) +
    geom_density(alpha=0.35, size=0.25) +
    guides(fill=guide_legend(title=NULL)) +
    scale_x_continuous(limits=c(-0.2, 0.8)) +
    xlab("Spearman Correlation") + ylab("Density") +
    theme_nature() +
    theme(legend.justification = c(1, 1), legend.position = c(1, 1))

pdf("./overview/correlation_persample_density_mgx-mtx-mpx.pdf", 2.37, 1.8)
print(ggp)
dev.off()

# summary statistics

summary(df_Summary_Corrs$Correlation[df_Summary_Corrs$between=="MGX-MTX"])
sd(df_Summary_Corrs$Correlation[df_Summary_Corrs$between=="MGX-MTX"])
IQR(df_Summary_Corrs$Correlation[df_Summary_Corrs$between=="MGX-MTX"])

summary(df_Summary_Corrs$Correlation[df_Summary_Corrs$between=="MGX-MPX"])
sd(df_Summary_Corrs$Correlation[df_Summary_Corrs$between=="MGX-MPX"])
IQR(df_Summary_Corrs$Correlation[df_Summary_Corrs$between=="MGX-MPX"])

summary(df_Summary_Corrs$Correlation[df_Summary_Corrs$between=="MTX-MPX"])
sd(df_Summary_Corrs$Correlation[df_Summary_Corrs$between=="MTX-MPX"])
IQR(df_Summary_Corrs$Correlation[df_Summary_Corrs$between=="MTX-MPX"])

