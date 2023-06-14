#### Remove All the Variables from the Workspace ####
#rm(list = ls())

#### Setting Directory ####
#setwd('/Users/hmallick/Dropbox (Personal)/Repos/hmp2_analysis')

###### The Required Package List ######
library(ggplot2); library(dplyr); library(cowplot); library(plyr); library(nlme)


# Load ALL Utility Functions At Once
# pkgmaker::source_files('common', '*.r')

# Load Particular Loaders
source("./common/load_bugs.r")
source("./env_config.r")
source("./common/disease_colors.r")
source("./common/disease_activity.r")
source("./common/theme_nature.r")

# Significance Testing with ANOVA
species.pcl<-pcl.only(bugs.pcl, rank="s")
gs <- pcl.apply.s(species.pcl, 1-sum(x^2))
is <- 1/(1-gs)
sh <- pcl.apply.s(species.pcl, -sum(x[x>0]*log2(x[x>0])))
df <- cbind.data.frame(species.pcl$meta, data.frame(gs=gs, is=is, sh=sh))
df$newfactor <- with(df, interaction(diagnosis,  active), drop = TRUE )
df <- transform(df, newfactor=revalue(newfactor,c("nonIBD.FALSE"="nonIBD",
                                                  "nonIBD.TRUE"="nonIBD",
                                                  "CD.TRUE"="Active CD",
                                                  "CD.FALSE"="Inactive CD",
                                                  "UC.TRUE"="Active UC",
                                                  "UC.FALSE"="Inactive UC")))

df$newfactor <- factor(df$newfactor, c("nonIBD", "Inactive UC", "Inactive CD", "Active UC", "Active CD"))

# 5-group comparison
aov1 <- aov(gs ~ newfactor, data = df); summary(aov1) # ANOVA, <2e-16
lme1<-lme(gs ~ newfactor, data = df, random = ~ 1 | subject); anova.lme(lme1)
summary(lme1)$tTable

# 3-group comparison
aov2 <- aov(gs ~ diagnosis, data = df); summary(aov2) # ANOVA, 1.86e-08
lme2<-lme(gs ~ diagnosis, data = df, random = ~ 1 | subject); anova.lme(lme2)
summary(lme2)$tTable

# Alpha diversity plots

plot_alpha_divs <- function(pcl, name) {
  gs <- pcl.apply.s(pcl, 1-sum(x^2))
  is <- 1/(1-gs)
  sh <- pcl.apply.s(pcl, -sum(x[x>0]*log2(x[x>0])))

  df <- cbind(pcl$meta, data.frame(gs=gs, is=is, sh=sh))
  df$newfactor <- with(df, interaction(diagnosis,  active), drop = TRUE )
  df <- transform(df, newfactor=revalue(newfactor,c("nonIBD.FALSE"="nonIBD",
                                                    "nonIBD.TRUE"="nonIBD",
                                                    "CD.TRUE"="Active CD",
                                                    "CD.FALSE"="Inactive CD",
                                                    "UC.TRUE"="Active UC",
                                                    "UC.FALSE"="Inactive UC")))

  df$newfactor <- factor(df$newfactor, c("nonIBD", "Inactive UC", "Inactive CD", "Active UC", "Active CD"))
  ggp_base <- ggplot(data=df, aes(x=newfactor, gs, fill=newfactor, color=newfactor)) +
    theme_cowplot() + xlab(NULL) +
      scale_color_manual(values=hmp2_active_disease_colors_outline,
                         labels = c("non-IBD",
                                    "Non-dysbiotic UC",
                                    "Non-dysbiotic CD",
                                    "Dysbiotic UC",
                                    "Dysbiotic CD")) +
      scale_fill_manual(values=hmp2_active_disease_colors_fill,
                        labels = c("non-IBD",
                                   "Non-dysbiotic UC",
                                   "Non-dysbiotic CD",
                                   "Dysbiotic UC",
                                   "Dysbiotic CD")) +
      guides(fill=guide_legend(title=""), color=guide_legend(title="")) +
      theme_nature() + theme(axis.title.x=element_blank(),
                             axis.text.x=element_blank(),
                             axis.ticks.x=element_blank()) +

  pdf(sprintf("./disease_activity/alphadiversity_%s.pdf", name), 2.3, 1.1)
  print(ggp_base + geom_boxplot(size=0.25, outlier.size=0.6, outlier.stroke=0) +
            geom_boxplot_n() +
          ylab("Gini-Simpson"))
  # print(ggp_base + geom_boxplot(aes(y=is), size=0.5, outlier.size=1) +
  #         ylab("Inverse Simpson"))
  # print(ggp_base + geom_boxplot(aes(y=sh), size=0.5, outlier.size=1) +
  #         ylab("Shannon Index (bits)"))
  dev.off()

  ggp <- ggplot(df, aes(x=activity_index, y=gs, fill=diagnosis)) +
      scale_fill_manual(values=hmp2_disease_colors, name="") +
      geom_vline(xintercept = disease_activity_threshold, size=0.7, color="black") +
      geom_point(size=0.7, shape=21, stroke=0.1) +
      ylab("Gini-Simpson") + xlab("Dysbiosis Score") +
      theme_nature()
  pdf(sprintf("./disease_activity/alphadiversity_%s_scatter.pdf", name), 2.3, 1.1)
  print(ggp)
  dev.off()
}

plot_alpha_divs(pcl.only(bugs.pcl, rank="s"), "dysbiosis_species")
plot_alpha_divs(pcl.only(bugs.pcl, rank="g"), "dysbiosis_genus")




