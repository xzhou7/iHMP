
source("./common/disease_activity.r")

# Make L123/L4 versions of the Montreal location
bugs.pcl$meta$baseline_montreal_l123 <- bugs.pcl$meta$baseline_montreal_location
bugs.pcl$meta$baseline_montreal_l123[bugs.pcl$meta$baseline_montreal_l123=="L1+L4"] <- "L1"
bugs.pcl$meta$baseline_montreal_l123[bugs.pcl$meta$baseline_montreal_l123=="L2+L4"] <- "L2"
bugs.pcl$meta$baseline_montreal_l123[bugs.pcl$meta$baseline_montreal_l123=="L3+L4"] <- "L3"
bugs.pcl$meta$baseline_montreal_l123 <- factor(bugs.pcl$meta$baseline_montreal_l123, levels=c("L1", "L2", "L3"))

bugs.pcl$meta$baseline_montreal_l4 <- grepl("L4", bugs.pcl$meta$baseline_montreal_location)


meta_cd <- bugs.pcl$meta[bugs.pcl$meta$diagnosis == "CD",,drop=F]

# Test overall activity -- Montreal location
with(meta_cd, table(active, baseline_montreal_location))
with(meta_cd, fisher.test(active, baseline_montreal_location))

# Test activity -- L123 location
with(meta_cd, table(active, baseline_montreal_l123))
with(meta_cd, fisher.test(active, baseline_montreal_l123))

# Test activity -- L4 location
with(meta_cd, table(active, baseline_montreal_l4))
with(meta_cd, fisher.test(active, baseline_montreal_l4))

# GLMER
library(lme4)
anova(glmer(active ~ baseline_montreal_location + ( 1 | subject), data=meta_cd, family=binomial)) #P=0.4233
anova(glmer(active ~ baseline_montreal_l123 + ( 1 | subject), data=meta_cd, family=binomial)) # P=0.1112
anova(glmer(active ~ baseline_montreal_l4 + ( 1 | subject), data=meta_cd, family=binomial)) # P=0.2543



