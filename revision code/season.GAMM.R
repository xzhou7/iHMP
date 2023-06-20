#Season Effect Using GAMM

library(mgcv)
library(tidyverse)
library(lubridate)

setwd("~/Library/CloudStorage/Box-Box/XinZhouFiles/Projects/ZXP2_HMP2_Core/Analysis/")

#########
#define two functions
check_asv_dates <- function(data, asv_col) {
  asv_data <- data %>% filter(!!sym(asv_col) != 0) %>% select(Date)
  return(n_distinct(asv_data$Date) > 10)
}

#Zscore for ASV
zscore_normalize <- function(x) {
  x_mean <- mean(x, na.rm = TRUE)
  x_sd <- sd(x, na.rm = TRUE)
  
  if (length(x) < 2 || is.na(x_sd) || x_sd == 0) {
    return(ifelse(is.na(x), NA, 0))
  } else {
    return((x - x_mean) / x_sd)}
}

# define a function for z-score normalization of diversity
z_score <- function(x) {
  (x - mean(x)) / sd(x)
}

####################################################################################################################################
#######
#load data
date.info <- read.csv(file = "./ori_meta_table/DateInfo.csv", header = T)
participant.info <- read.csv(file = './ori_meta_table/metadata.subject.csv', header = T)

stool.diversity <- read.csv("./Diversity Table/Stool.Diversity.csv", header = T)
skin.diversity <- read.csv("./Diversity Table/Skin.Diversity.csv",header = T)
oral.diversity <- read.csv("./Diversity Table/Oral.Diversity.csv",header = T)
nasal.diversity <- read.csv("./Diversity Table/Nasal.Diversity.csv",header = T)

####################################################################################################################################
#####stool
stool.otu <- otu_table(physeqGenus_ST) %>% data.frame()
stool.otu[1:5,1:5]

stool.sample <- sample_data(physeqGenus_ST) %>% data.frame()
stool.sample

stool.taxa <- tax_table(physeqGenus_ST) %>% data.frame()
stool.taxa$ASV <- rownames(stool.taxa)

get_stool_genus <- function(ASV) {
  genus <- stool.taxa$Genus[which(stool.taxa$ASV == ASV)]
  return(genus)
}

output_genus <- get_stool_genus("ASV1")
output_genus

identical(rownames(stool.otu),stool.sample$RandomID)

stool.otu$Subject <- stool.sample$SubjectID
stool.otu <- select(stool.otu, Subject, everything())
stool.otu[1:5,1:5]

filter(stool.otu, Subject == "69-118") %>% dim()

Z.stool.rb <- stool.otu %>%
  group_by(Subject) %>%
  mutate(across(starts_with("ASV"), zscore_normalize))

Z.stool.rb$SampleID <- stool.sample$SampleID

Z.stool.rb.date <- Z.stool.rb %>%
  left_join(date.info[, c("SampleID", "Date")], by = "SampleID")

Z.stool.rb.date <- Z.stool.rb.date %>% select(SampleID,Date, everything())
Z.stool.rb.date
Z.stool.rb.date$Date = as.Date(Z.stool.rb.date$Date, format="%m/%d/%y")
Z.stool.rb.date$Time = yday(Z.stool.rb.date$Date)

Z.stool.rb.date

asv_columns <- colnames(Z.stool.rb.date) %>% str_subset("^ASV")
filtered_asv_columns <- asv_columns %>% keep(~check_asv_dates(Z.stool.rb.date, .))

removed_asv_columns <- setdiff(asv_columns, filtered_asv_columns)

filtered_Z_stool_rb_date <- Z.stool.rb.date %>% select(c("SampleID", "Date", "Subject","Time", all_of(filtered_asv_columns)))
filtered_Z_stool_rb_date$IRIS <- participant.info$IRIS[match(filtered_Z_stool_rb_date$Subject,participant.info$SubjectID)]
#how many taxa are removed
length(removed_asv_columns)

ctrl <- lmeControl(opt='optim')

stool.gene.list <- filtered_Z_stool_rb_date %>% select(starts_with("ASV")) %>% colnames()
# 
# set.seed(77777)
# p.vals = c()
# coeffs = c()
# 

# 
# mod  <- gamm(ASV2420 ~ IRIS + s(Time, bs = "cc"),
#              data=filtered_Z_stool_rb_date,
#              method = "REML",
#              control=ctrl,
#              random = list(Subject = ~ 1),
#              knots = list(TimeOfYear = c(0, 366)))
# 
# plot(mod$gam, shade=T,shade.col="#9a9a00",xlab = "Day of the year", ylab = get_stool_genus("ASV1"),cex.lab=1.5,font.lab=2, cex.axis=1.5,font=2, main=get_stool_genus("ASV1"), font.main=2,cex.main=2.5)
# 
# summary(mod$gam)
# 
# mod$gam$coefficients

stool.gene.names <- setdiff(stool.gene.list, "Subject")
stool.gene.names <- setdiff(stool.gene.names, "ASV2420")
stool.p.vals = c()
stool.IRIS = c()
stool.coeffs = c()

for (gene in stool.gene.names) {
  path <- paste0("~/Desktop/season/stool/plots/",gene,".png")
  mod <- gamm(as.formula(paste(gene, " ~ IRIS + s(Time, bs = \"cc\")")),
              data = filtered_Z_stool_rb_date, 
              method = "REML",
              control=ctrl,
              random = list(Subject = ~ 1),
              knots = list(TimeOfYear = c(0, 366)))
  
  png(path)
  
  plot(mod$gam, shade=T,shade.col="#9a9a00",xlab = "Day of the year", ylab = get_stool_genus(gene),cex.lab=1.5,font.lab=2, cex.axis=1.5,font=2, main=get_stool_genus(gene), font.main=2,cex.main=2.5)
  
  
  dev.off()
  
  stool.p.vals = rbind(stool.p.vals, summary(mod$gam)$s.table[4])
  stool.IRIS = rbind(stool.IRIS, summary(mod$gam)$pTerms.table[3])
  stool.coeffs = rbind(stool.coeffs, mod$gam$coefficients)
  
  print(gene)
}

table(stool.p.vals < 0.05)
stool.coeffs

#######perform the same analysis for diversity
stool.diversity$SampleID <- stool.sample$SampleID[match(stool.diversity$X, stool.sample$RandomID)]
stool.diversity.date <- merge(stool.diversity, select(filtered_Z_stool_rb_date, SampleID:Time), by = "SampleID")
stool.diversity.date$IRIS <- participant.info$IRIS[match(stool.diversity.date$Subject,participant.info$SubjectID)]
stool.diversity.date

stool.diversity.list <- select(stool.diversity.date, Observed:bulla_e) %>% colnames()
stool.diversity.list

#loop starts here
stool.diversity.p.vals = c()
stool.diversity.IRIS.p =c()
stool.diversity.coeffs = c()

for (gene in stool.diversity.list) {
  path <- paste0("~/Desktop/season/stool/plots/Diversity_",gene,".png")
  mod.1 <- gamm(as.formula(paste(gene, " ~ IRIS + s(Time, bs = \"cc\")")),
              data = stool.diversity.date, 
              method = "REML",
              control=ctrl,
              random = list(Subject = ~ 1),
              knots = list(TimeOfYear = c(0, 366)))
  
  png(path)
  
  plot(mod.1$gam, shade=T,shade.col="#9a9a00",xlab = "Day of the year", ylab = paste(gene),cex.lab=1.5,font.lab=2, cex.axis=1.5,font=2, main=paste(gene), font.main=2,cex.main=2.5)
  
  dev.off()
  
  stool.diversity.p.vals = rbind(stool.diversity.p.vals, summary(mod.1$gam)$s.table[4])
  stool.diversity.IRIS.p = rbind(stool.diversity.IRIS.p, summary(mod.1$gam)$pTerms.table[3])
  stool.diversity.coeffs = rbind(stool.diversity.coeffs, mod.1$gam$coefficients)
  
  print(gene)
}

stool.diversity.p.vals
stool.diversity.IRIS.p
stool.diversity_gamm <- cbind(stool.diversity.p.vals,stool.diversity.list) %>% data.frame()

####################################################################################################################################
#skin
skin.otu <- otu_table(physeqGenus_SK) %>% data.frame()
skin.otu[1:5,1:5]

skin.sample <- sample_data(physeqGenus_SK) %>% data.frame()
skin.sample

skin.taxa <- tax_table(physeqGenus_SK) %>% data.frame()
skin.taxa$ASV <- rownames(skin.taxa)

get_skin_genus <- function(ASV) {
  genus <- skin.taxa$Genus[which(skin.taxa$ASV == ASV)]
  return(genus)
}

output_genus <- get_skin_genus("OTU_1")
output_genus

identical(rownames(skin.otu), skin.sample$KitID)

skin.otu$Subject <- skin.sample$SubjectID
skin.otu <- select(skin.otu, Subject, everything())
skin.otu[1:5, 1:5]

filter(skin.otu, Subject == "69-118") %>% dim()

Z.skin.rb <- skin.otu %>%
  group_by(Subject) %>%
  mutate(across(starts_with("OTU"), zscore_normalize))

Z.skin.rb$SampleID <- skin.sample$SampleID

Z.skin.rb.date <- Z.skin.rb %>%
  left_join(date.info[, c("SampleID", "Date")], by = "SampleID")

Z.skin.rb.date <- Z.skin.rb.date %>% select(SampleID, Date, everything())
Z.skin.rb.date
Z.skin.rb.date$Date = as.Date(Z.skin.rb.date$Date, format="%m/%d/%y")
Z.skin.rb.date$Time = yday(Z.skin.rb.date$Date)

Z.skin.rb.date

asv_columns <- colnames(Z.skin.rb.date) %>% str_subset("^OTU")
filtered_asv_columns <- asv_columns %>% keep(~check_asv_dates(Z.skin.rb.date, .))

removed_asv_columns <- setdiff(asv_columns, filtered_asv_columns)

filtered_Z_skin_rb_date <- Z.skin.rb.date %>% select(c("SampleID", "Date", "Subject", "Time", all_of(filtered_asv_columns)))
filtered_Z_skin_rb_date$IRIS <- participant.info$IRIS[match(filtered_Z_skin_rb_date$Subject, participant.info$SubjectID)]

length(removed_asv_columns)

skin.gene.list <- filtered_Z_skin_rb_date %>% select(starts_with("OTU")) %>% colnames()

set.seed(77777)

skin.gene.names <- setdiff(skin.gene.list, "Subject")
skin.gene.names <- setdiff(skin.gene.names, "ASV2420")
skin.p.vals = c()
skin.IRIS = c()
skin.coeffs = c()

for (gene in skin.gene.names) {
  path <- paste0("~/Desktop/season/skin/plots/", gene, ".png")
  mod <- gamm(as.formula(paste(gene, " ~ IRIS + s(Time, bs = \"cc\")")),
              data = filtered_Z_skin_rb_date, 
              method = "REML",
              control = ctrl,
              random = list(Subject = ~ 1),
              knots = list(TimeOfYear = c(0, 366)))
  
  png(path)
  
  plot(mod$gam, shade = T, shade.col = "#9a9a00", xlab = "Day of the year", ylab = get_skin_genus(gene), cex.lab = 1.5, font.lab = 2, cex.axis = 1.5, font = 2, main = get_skin_genus(gene), font.main = 2, cex.main = 2.5)
  
  dev.off()
  
  skin.p.vals = rbind(skin.p.vals, summary(mod$gam)$s.table[4])
  skin.IRIS = rbind(skin.IRIS, summary(mod$gam)$pTerms.table[3])
  skin.coeffs = rbind(skin.coeffs, mod$gam$coefficients)
  
  print(gene)
}

table(skin.p.vals < 0.05)
skin.p.vals
skin.coeffs

#######perform the same analysis for diversity
skin.diversity$SampleID <- skin.sample$SampleID[match(skin.diversity$X, skin.sample$KitID)]
skin.diversity.date <- merge(skin.diversity, select(filtered_Z_skin_rb_date, SampleID:Time), by = "SampleID")
skin.diversity.date$IRIS <- participant.info$IRIS[match(skin.diversity.date$Subject, participant.info$SubjectID)]
skin.diversity.date

skin.diversity.list <- select(skin.diversity.date, Observed:bulla_e) %>% colnames()
skin.diversity.list

#loop starts here
skin.diversity.p.vals = c()
skin.diversity.IRIS.p = c()
skin.diversity.coeffs = c()

for (gene in skin.diversity.list) {
  path <- paste0("~/Desktop/season/skin/plots/Diversity_", gene, ".png")
  mod.1 <- gamm(as.formula(paste(gene, " ~ IRIS + s(Time, bs = \"cc\")")),
                data = skin.diversity.date, 
                method = "REML",
                control = ctrl,
                random = list(Subject = ~ 1),
                knots = list(TimeOfYear = c(0, 366)))
  
  png(path)
  
  plot(mod.1$gam, shade = T, shade.col = "#9a9a00", xlab = "Day of the year", ylab = paste(gene), cex.lab = 1.5, font.lab = 2, cex.axis = 1.5, font = 2, main = paste(gene), font.main = 2, cex.main = 2.5)
  
  dev.off()
  
  skin.diversity.p.vals = rbind(skin.diversity.p.vals, summary(mod.1$gam)$s.table[4])
  skin.diversity.IRIS.p = rbind(skin.diversity.IRIS.p, summary(mod.1$gam)$pTerms.table[3])
  skin.diversity.coeffs = rbind(skin.diversity.coeffs, mod.1$gam$coefficients)
  
  print(gene)
}

skin.diversity.p.vals
skin.diversity.IRIS.p
skin.diversity_gamm <- cbind(skin.diversity.p.vals,skin.diversity.list) %>% data.frame()

####################################################################################################################################
#oral
oral.otu <- otu_table(physeqGenus_OR) %>% data.frame()
oral.otu[1:5,1:5]

oral.sample <- sample_data(physeqGenus_OR) %>% data.frame()
oral.sample

oral.taxa <- tax_table(physeqGenus_OR) %>% data.frame()
oral.taxa$ASV <- rownames(oral.taxa)

get_oral_genus <- function(ASV) {
  genus <- oral.taxa$Genus[which(oral.taxa$ASV == ASV)]
  return(genus)
}

output_genus <- get_oral_genus("OTU_1")
output_genus

identical(rownames(oral.otu), oral.sample$KitID)

oral.otu$Subject <- oral.sample$SubjectID
oral.otu <- select(oral.otu, Subject, everything())
oral.otu[1:5, 1:5]

filter(oral.otu, Subject == "69-118") %>% dim()

Z.oral.rb <- oral.otu %>%
  group_by(Subject) %>%
  mutate(across(starts_with("OTU"), zscore_normalize))

Z.oral.rb$SampleID <- oral.sample$SampleID

Z.oral.rb.date <- Z.oral.rb %>%
  left_join(date.info[, c("SampleID", "Date")], by = "SampleID")

Z.oral.rb.date <- Z.oral.rb.date %>% select(SampleID, Date, everything())
Z.oral.rb.date
Z.oral.rb.date$Date = as.Date(Z.oral.rb.date$Date, format="%m/%d/%y")
Z.oral.rb.date$Time = yday(Z.oral.rb.date$Date)

Z.oral.rb.date

asv_columns <- colnames(Z.oral.rb.date) %>% str_subset("^OTU")
filtered_asv_columns <- asv_columns %>% keep(~check_asv_dates(Z.oral.rb.date, .))

removed_asv_columns <- setdiff(asv_columns, filtered_asv_columns)

filtered_Z_oral_rb_date <- Z.oral.rb.date %>% select(c("SampleID", "Date", "Subject", "Time", all_of(filtered_asv_columns)))
filtered_Z_oral_rb_date$IRIS <- participant.info$IRIS[match(filtered_Z_oral_rb_date$Subject, participant.info$SubjectID)]

length(removed_asv_columns)

oral.gene.list <- filtered_Z_oral_rb_date %>% select(starts_with("OTU")) %>% colnames()

set.seed(77777)

oral.gene.names <- setdiff(oral.gene.list, "Subject")
oral.gene.names <- setdiff(oral.gene.names, "OTU_207")
oral.p.vals = c()
oral.IRIS = c()
oral.coeffs = c()

for (gene in oral.gene.names) {
  path <- paste0("~/Desktop/season/oral/plots/", gene, ".png")
  mod <- gamm(as.formula(paste(gene, " ~ IRIS + s(Time, bs = \"cc\")")),
              data = filtered_Z_oral_rb_date, 
              method = "REML",
              control = ctrl,
              random = list(Subject = ~ 1),
              knots = list(TimeOfYear = c(0, 366)))
  
  png(path)
  
  plot(mod$gam, shade = T, shade.col = "#9a9a00", xlab = "Day of the year", ylab = get_oral_genus(gene), cex.lab = 1.5, font.lab = 2, cex.axis = 1.5, font = 2, main = get_oral_genus(gene), font.main = 2, cex.main = 2.5)
  
  dev.off()
  
  oral.p.vals = rbind(oral.p.vals, summary(mod$gam)$s.table[4])
  oral.IRIS = rbind(oral.IRIS, summary(mod$gam)$pTerms.table[3])
  oral.coeffs = rbind(oral.coeffs, mod$gam$coefficients)
  
  print(gene)
}

table(oral.p.vals < 0.05)
oral.coeffs

#######perform the same analysis for diversity
oral.diversity$SampleID <- oral.sample$SampleID[match(oral.diversity$X, oral.sample$KitID)]
oral.diversity.date <- merge(oral.diversity, select(filtered_Z_oral_rb_date, SampleID:Time), by = "SampleID")
oral.diversity.date$IRIS <- participant.info$IRIS[match(oral.diversity.date$Subject, participant.info$SubjectID)]
oral.diversity.date

oral.diversity.list <- select(oral.diversity.date, Observed:bulla_e) %>% colnames()
oral.diversity.list

#loop starts here
oral.diversity.p.vals = c()
oral.diversity.IRIS.p = c()
oral.diversity.coeffs = c()

for (gene in oral.diversity.list) {
  path <- paste0("~/Desktop/season/oral/plots/Diversity_", gene, ".png")
  mod.1 <- gamm(as.formula(paste(gene, " ~ IRIS + s(Time, bs = \"cc\")")),
                data = oral.diversity.date, 
                method = "REML",
                control = ctrl,
                random = list(Subject = ~ 1),
                knots = list(TimeOfYear = c(0, 366)))
  
  png(path)
  
  plot(mod.1$gam, shade = T, shade.col = "#9a9a00", xlab = "Day of the year", ylab = get_oral_genus(gene), cex.lab = 1.5, font.lab = 2, cex.axis = 1.5, font = 2, main = get_oral_genus(gene), font.main = 2, cex.main = 2.5)
  
  dev.off()
  
  oral.diversity.p.vals = rbind(oral.diversity.p.vals, summary(mod.1$gam)$s.table[4])
  oral.diversity.IRIS.p = rbind(oral.diversity.IRIS.p, summary(mod.1$gam)$pTerms.table[3])
  oral.diversity.coeffs = rbind(oral.diversity.coeffs, mod.1$gam$coefficients)
  
  print(gene)
}

oral.diversity.p.vals
oral.diversity.IRIS.p
oral.diversity_gamm <- cbind(oral.diversity.p.vals,oral.diversity.list) %>% data.frame()

####################################################################################################################################
#nasal
nasal.otu <- otu_table(physeqGenus_NS) %>% data.frame()
nasal.otu[1:5,1:5]

nasal.sample <- sample_data(physeqGenus_NS) %>% data.frame()
nasal.sample

nasal.taxa <- tax_table(physeqGenus_NS) %>% data.frame()
nasal.taxa$ASV <- rownames(nasal.taxa)

get_nasal_genus <- function(ASV) {
  genus <- nasal.taxa$Genus[which(nasal.taxa$ASV == ASV)]
  return(genus)
}

output_genus <- get_nasal_genus("ASV1")
output_genus

identical(rownames(nasal.otu), nasal.sample$RandomID)

nasal.otu$Subject <- nasal.sample$SubjectID
nasal.otu <- select(nasal.otu, Subject, everything())
nasal.otu[1:5, 1:5]

filter(nasal.otu, Subject == "69-118") %>% dim()

Z.nasal.rb <- nasal.otu %>%
  group_by(Subject) %>%
  mutate(across(starts_with("ASV"), zscore_normalize))

Z.nasal.rb$SampleID <- nasal.sample$SampleID

Z.nasal.rb.date <- Z.nasal.rb %>%
  left_join(date.info[, c("SampleID", "Date")], by = "SampleID")

Z.nasal.rb.date <- Z.nasal.rb.date %>% select(SampleID, Date, everything())
Z.nasal.rb.date
Z.nasal.rb.date$Date = as.Date(Z.nasal.rb.date$Date, format="%m/%d/%y")
Z.nasal.rb.date$Time = yday(Z.nasal.rb.date$Date)

Z.nasal.rb.date

asv_columns <- colnames(Z.nasal.rb.date) %>% str_subset("^ASV")
filtered_asv_columns <- asv_columns %>% keep(~check_asv_dates(Z.nasal.rb.date, .))

removed_asv_columns <- setdiff(asv_columns, filtered_asv_columns)

filtered_Z_nasal_rb_date <- Z.nasal.rb.date %>% select(c("SampleID", "Date", "Subject", "Time", all_of(filtered_asv_columns)))
filtered_Z_nasal_rb_date$IRIS <- participant.info$IRIS[match(filtered_Z_nasal_rb_date$Subject, participant.info$SubjectID)]

#how many taxa are removed
length(removed_asv_columns)

nasal.gene.list <- filtered_Z_nasal_rb_date %>% select(starts_with("ASV")) %>% colnames()

nasal.gene.names <- setdiff(nasal.gene.list, "Subject")
nasal.gene.names <- setdiff(nasal.gene.names, "ASV42971")
nasal.p.vals = c()
nasal.IRIS = c()
nasal.coeffs = c()

for (gene in nasal.gene.names) {
  path <- paste0("~/Desktop/season/nasal/plots/", gene, ".png")
  mod <- gamm(as.formula(paste(gene, " ~ IRIS + s(Time, bs = \"cc\")")),
              data = filtered_Z_nasal_rb_date, 
              method = "REML",
              control = ctrl,
              random = list(Subject = ~ 1),
              knots = list(TimeOfYear = c(0, 366)))
  
  png(path)
  
  plot(mod$gam, shade = T, shade.col = "#9a9a00", xlab = "Day of the year", ylab = get_nasal_genus(gene), cex.lab = 1.5, font.lab = 2, cex.axis = 1.5, font = 2, main = get_nasal_genus(gene), font.main = 2, cex.main = 2.5)
  
  dev.off()
  
  nasal.p.vals = rbind(nasal.p.vals, summary(mod$gam)$s.table[4])
  nasal.IRIS = rbind(nasal.IRIS, summary(mod$gam)$pTerms.table[3])
  nasal.coeffs = rbind(nasal.coeffs, mod$gam$coefficients)
  
  print(gene)
}

table(nasal.p.vals < 0.05)
nasal.coeffs

#######perform the same analysis for diversity
nasal.diversity$SampleID <- nasal.sample$SampleID[match(nasal.diversity$X, nasal.sample$RandomID)]
nasal.diversity.date <- merge(nasal.diversity, select(filtered_Z_nasal_rb_date, SampleID:Time), by = "SampleID")
nasal.diversity.date$IRIS <- participant.info$IRIS[match(nasal.diversity.date$Subject,participant.info$SubjectID)]
nasal.diversity.date

nasal.diversity.list <- select(nasal.diversity.date, Observed:bulla_e) %>% colnames()
nasal.diversity.list

#loop starts here
nasal.diversity.p.vals = c()
nasal.diversity.IRIS.p = c()
nasal.diversity.coeffs = c()

for (gene in nasal.diversity.list) {
  path <- paste0("~/Desktop/season/nasal/plots/Diversity_", gene, ".png")
  mod.1 <- gamm(as.formula(paste(gene, " ~ IRIS + s(Time, bs = \"cc\")")),
                data = nasal.diversity.date, 
                method = "REML",
                control = ctrl,
                random = list(Subject = ~ 1),
                knots = list(TimeOfYear = c(0, 366)))
  
  png(path)
  
  plot(mod.1$gam, shade = T, shade.col = "#9a9a00", xlab = "Day of the year", ylab = paste(gene), cex.lab = 1.5, font.lab = 2, cex.axis = 1.5, font = 2, main = paste(gene), font.main = 2, cex.main = 2.5)
  
  dev.off()
  
  nasal.diversity.p.vals = rbind(nasal.diversity.p.vals, summary(mod.1$gam)$s.table[4])
  nasal.diversity.IRIS.p = rbind(nasal.diversity.IRIS.p, summary(mod.1$gam)$pTerms.table[3])
  nasal.diversity.coeffs = rbind(nasal.diversity.coeffs, mod.1$gam$coefficients)
  
  print(gene)
}

nasal.diversity.p.vals
nasal.diversity.IRIS.p
nasal.diversity_gamm <- cbind(nasal.diversity.p.vals, nasal.diversity.list) %>% data.frame()


colnames(stool.diversity_gamm)[1] <- "Stool"
colnames(skin.diversity_gamm)[1] <- "Skin"
colnames(oral.diversity_gamm)[1] <- "Oral"
colnames(nasal.diversity_gamm)[1] <- "Nasal"

GAMM_diversity_table <- cbind(stool.diversity_gamm,skin.diversity_gamm,oral.diversity_gamm,nasal.diversity_gamm)
write.csv(file = "~/Desktop/season/Diversity.P.csv",GAMM_diversity_table)

stool.p.result <- cbind(stool.gene.names,stool.p.vals) %>% data.frame()
stool.p.result$Genus <- stool.taxa$Genus[match(stool.p.result$stool.gene.names, stool.taxa$ASV)]
stool.p.result$bodysite <- "Stool"

stool.p.result$p.adj <- p.adjust(stool.p.result$V2, method = "BH", n = length(stool.p.vals))

table(stool.p.result$p.adj < 0.05)

skin.p.result <- cbind(skin.gene.names, skin.p.vals) %>% data.frame()
skin.p.result$Genus <- skin.taxa$Genus[match(skin.p.result$skin.gene.names, skin.taxa$ASV)]
skin.p.result$bodysite <- "Skin"

skin.p.result$p.adj <- p.adjust(skin.p.result$V2, method = "BH", n = length(skin.p.vals))

table(skin.p.result$p.adj < 0.05)

oral.p.result <- cbind(oral.gene.names, oral.p.vals) %>% data.frame()
oral.p.result$Genus <- oral.taxa$Genus[match(oral.p.result$oral.gene.names, oral.taxa$ASV)]
oral.p.result$bodysite <- "Oral"

oral.p.result$p.adj <- p.adjust(oral.p.result$V2, method = "BH", n = length(oral.p.vals))

table(oral.p.result$p.adj < 0.05)

nasal.p.result <- cbind(nasal.gene.names, nasal.p.vals) %>% data.frame()
nasal.p.result$Genus <- nasal.taxa$Genus[match(nasal.p.result$nasal.gene.names, nasal.taxa$ASV)]
nasal.p.result$bodysite <- "Nasal"

nasal.p.result$p.adj <- p.adjust(nasal.p.result$V2, method = "BH", n = length(nasal.p.vals))

table(nasal.p.result$p.adj < 0.05)

names(stool.p.result)[1] <- "gene.names"
names(skin.p.result)[1] <- "gene.names"
names(oral.p.result)[1] <- "gene.names"
names(nasal.p.result)[1] <- "gene.names"

combined_GAMM.p <- rbind(stool.p.result, skin.p.result, oral.p.result, nasal.p.result)

combined_GAMM.p.csv <- select(combined_GAMM.p, bodysite,Genus,V2,p.adj) %>% filter(V2 < 0.05)
combined_GAMM.p.csv

table(combined_GAMM.p.csv$bodysite)

#remove non-bacteria signal
combined_GAMM.p.csv <- combined_GAMM.p.csv %>% filter(Genus != "Zea")

#write.csv(file = "~/Desktop/season/GAMM.pvalue.csv",combined_GAMM.p.csv)

################
#Analysis of infection in the same way

infection.meta <- read.csv(file = "./ori_meta_table/hmp2_infection_metadata_1027.csv", header = T)
infection.meta

infection.meta_filtered <- filter(infection.meta, InfectionState_jax %in% c("pre", "early", "late", "recovery", "post"))

# Create a vector of mapping values
mapping <- c("pre" = 1, "early" = 2, "late" = 3, "recovery" = 4, "post" = 5)

# Apply the mapping to create the new column
infection.meta_filtered$Time <- mapping[infection.meta_filtered$InfectionState_jax]

# Initialize the event counter
event_counter <- 1

# Initialize the current_subject
current_subject <- infection.meta_filtered$Subject[1]

# Initialize the event column with NA
infection.meta_filtered$event <- NA

# Loop through each row
for (i in 1:nrow(infection.meta_filtered)) {
  # Check if the subject has changed
  if (infection.meta_filtered$Subject[i] != current_subject) {
    # Update the current subject
    current_subject <- infection.meta_filtered$Subject[i]
    # Reset the event counter
    event_counter <- 1
  }
  # Check if the value in "Time" for the current row is smaller than the above row (if not the first row)
  if (i > 1 && infection.meta_filtered$Subject[i-1] == current_subject && 
      infection.meta_filtered$Time[i] < infection.meta_filtered$Time[i-1]) {
    # Increment the event counter
    event_counter <- event_counter + 1
  }
  # Assign the current event counter value to the "event" column
  infection.meta_filtered$event[i] <- paste(infection.meta_filtered$Subject[i], event_counter, sep = "_")
}

event_counter
infection.meta_filtered

write.csv(file = "~/Desktop/infection.meta_filtered.csv", infection.meta_filtered)
# 
# # Create a vector of mapping values
# mapping <- c("pre" = 1, "early" = 2, "late" = 3, "recovery" = 4, "post" = 5)
# 
# # Apply the mapping to create the new column
# infection.meta_filtered$Time <- mapping[infection.meta_filtered$InfectionState_jax]
# 
# # Initialize the event counter
# event_counter <- 1
# 
# # Initialize the event column with NA
# infection.meta_filtered$event <- NA
# 
# # Loop through each row starting from the second row
# for (i in 2:nrow(infection.meta_filtered)) {
#   # Check if the value in "Time" for the current row is smaller than the above row
#   if (infection.meta_filtered$Time[i] < infection.meta_filtered$Time[i-1]) {
#     # Increment the event counter
#     event_counter <- event_counter + 1
#   }
#   # Assign the current event counter value to the "event" column
#   infection.meta_filtered$event[i] <- event_counter
# }
# infection.meta_filtered$event[1] <- 1

infection.z.stool <- merge(infection.meta_filtered, filtered_Z_stool_rb_date, by="SampleID")
infection.z.skin <- merge(infection.meta_filtered, filtered_Z_skin_rb_date, by="SampleID")
infection.z.oral <- merge(infection.meta_filtered, filtered_Z_oral_rb_date, by="SampleID")
infection.z.nasal <- merge(infection.meta_filtered, filtered_Z_nasal_rb_date, by="SampleID")


infection.diversity.stool <- merge(infection.meta_filtered, stool.diversity.date, by="SampleID") 
infection.diversity.skin <- merge(infection.meta_filtered, skin.diversity.date, by="SampleID") 
infection.diversity.oral <- merge(infection.meta_filtered, oral.diversity.date, by="SampleID") 
infection.diversity.nasal <- merge(infection.meta_filtered, nasal.diversity.date, by="SampleID") 


# apply the function to the desired columns
infection.diversity.z.stool <- infection.diversity.stool %>%
  group_by(Subject.x) %>%
  mutate(across(Observed:bulla_e, z_score)) 

infection.diversity.z.skin <- infection.diversity.skin %>%
  group_by(Subject.x) %>%
  mutate(across(Observed:bulla_e, z_score)) 

infection.diversity.z.oral <- infection.diversity.oral %>%
  group_by(Subject.x) %>%
  mutate(across(Observed:bulla_e, z_score)) 

infection.diversity.z.nasal <- infection.diversity.nasal %>%
  group_by(Subject.x) %>%
  mutate(across(Observed:bulla_e, z_score)) 

stool.diversity.list <- infection.diversity.z.stool %>% select(Observed:bulla_e) %>% colnames() %>% as.character()

ir.list <- subset(participant.info, IRIS == "IR")
is.list <- subset(participant.info, IRIS == "IS")

infection.diversity.stool

#loop starts here
stool.diversity.p.vals.IF = c()
stool.diversity.IRIS.p.IF = c()
stool.diversity.coeffs.IF = c()

stool.diversity.list <-  stool.diversity.list[stool.diversity.list != "Subject.x"]

for (gene in stool.diversity.list) {
  path <- paste0("~/Desktop/infection/stool/Diversity_", gene, ".png")
  mod.1 <- gamm(as.formula(paste(gene, " ~  IRIS + s(Time.x, bs = \"cc\", k=5)")),
                data = infection.diversity.z.stool, 
                method = "REML",
                control = ctrl,
                random = list(event = ~ 1))
  
  png(path)
  
  plot(mod.1$gam, shade = T, shade.col = "#DF8F44FF", xlab = "Infection Stage", ylab = paste(gene), cex.lab = 1.5, font.lab = 2, cex.axis = 1.5, font = 2, main = paste(gene), font.main = 2, cex.main = 2.5)
  
  dev.off()
  
  stool.diversity.p.vals.IF = rbind(stool.diversity.p.vals.IF, summary(mod.1$gam)$s.table[4])
  stool.diversity.IRIS.p.IF = rbind(stool.diversity.IRIS.p.IF, summary(mod.1$gam)$pTerms.table[3])
  stool.diversity.coeffs.IF = rbind(stool.diversity.coeffs.IF, mod.1$gam$coefficients)
  
  print(gene)
}

stool.infection.diversity <- cbind(stool.diversity.list,stool.diversity.p.vals.IF,stool.diversity.IRIS.p.IF)
colnames(stool.infection.diversity) <- c("measurement", "p-value", "IRIS_p")
stool.infection.diversity

summary(mod.1$gam)

# Skin Loop starts here
skin.diversity.p.vals.IF <- c()
skin.diversity.IRIS.p.IF <- c()
skin.diversity.coeffs.IF <- c()

for (gene in skin.diversity.list) {
  path <- paste0("~/Desktop/infection/skin/Diversity_", gene, ".png")
  mod.1 <- gamm(as.formula(paste(gene, " ~ IRIS + s(Time.x, bs = \"cc\", k=5)")),
                data = infection.diversity.skin, 
                method = "REML",
                control = ctrl,
                random = list(event = ~ 1))
  
  png(path)
  
  plot(mod.1$gam, shade = T, shade.col = "#00A1D5FF", xlab = "Infection Stage", ylab = paste(gene), cex.lab = 1.5, font.lab = 2, cex.axis = 1.5, font = 2, main = paste(gene), font.main = 2, cex.main = 2.5)
  
  dev.off()
  
  skin.diversity.p.vals.IF <- rbind(skin.diversity.p.vals.IF, summary(mod.1$gam)$s.table[4])
  skin.diversity.IRIS.p.IF <- rbind(skin.diversity.IRIS.p.IF, summary(mod.1$gam)$pTerms.table[3])
  skin.diversity.coeffs.IF <- rbind(skin.diversity.coeffs.IF, mod.1$gam$coefficients)
  
  print(gene)
}

skin.infection.diversity <- cbind(skin.diversity.list, skin.diversity.p.vals.IF, skin.diversity.IRIS.p.IF)
colnames(skin.infection.diversity) <- c("measurement", "p-value", "IRIS_p")
skin.infection.diversity


# Loop starts here
oral.diversity.p.vals.IF <- c()
oral.diversity.IRIS.p.IF <- c()
oral.diversity.coeffs.IF <- c()

for (gene in oral.diversity.list) {
  path <- paste0("~/Desktop/infection/oral/Diversity_", gene, ".png")
  mod.1 <- gamm(as.formula(paste(gene, " ~ IRIS + s(Time.x, bs = \"cc\", k=5)")),
                data = infection.diversity.oral, 
                method = "REML",
                control = ctrl,
                random = list(event = ~ 1))
  
  png(path)
  
  plot(mod.1$gam, shade = T, shade.col = "#B24745FF", xlab = "Infection Stage", ylab = paste(gene), cex.lab = 1.5, font.lab = 2, cex.axis = 1.5, font = 2, main = paste(gene), font.main = 2, cex.main = 2.5)
  
  dev.off()
  
  oral.diversity.p.vals.IF <- rbind(oral.diversity.p.vals.IF, summary(mod.1$gam)$s.table[4])
  oral.diversity.IRIS.p.IF <- rbind(oral.diversity.IRIS.p.IF, summary(mod.1$gam)$pTerms.table[3])
  oral.diversity.coeffs.IF <- rbind(oral.diversity.coeffs.IF, mod.1$gam$coefficients)
  
  print(gene)
}

oral.infection.diversity <- cbind(oral.diversity.list, oral.diversity.p.vals.IF, oral.diversity.IRIS.p.IF)
colnames(oral.infection.diversity) <- c("measurement", "p-value", "IRIS_p")
oral.infection.diversity


# Loop starts here
nasal.diversity.p.vals.IF <- c()
nasal.diversity.IRIS.p.IF <- c()
nasal.diversity.coeffs.IF <- c()

nasal.diversity.list <- nasal.diversity.list[nasal.diversity.list != "simpson_e"]
nasal.diversity.list <- nasal.diversity.list[nasal.diversity.list != "evar_e"]

for (gene in nasal.diversity.list) {
  path <- paste0("~/Desktop/infection/nasal/Diversity_", gene, ".png")
  mod.1 <- gamm(as.formula(paste(gene, " ~ IRIS + s(Time.x, bs = \"cc\", k=5)")),
                data = infection.diversity.nasal, 
                method = "REML",
                control = ctrl,
                random = list(event = ~ 1))
  
  png(path)
  
  plot(mod.1$gam, shade = T, shade.col = "#79AF97FF", xlab = "Infection Stage", ylab = paste(gene), cex.lab = 1.5, font.lab = 2, cex.axis = 1.5, font = 2, main = paste(gene), font.main = 2, cex.main = 2.5)
  
  dev.off()
  
  nasal.diversity.p.vals.IF <- rbind(nasal.diversity.p.vals.IF, summary(mod.1$gam)$s.table[4])
  nasal.diversity.IRIS.p.IF <- rbind(nasal.diversity.IRIS.p.IF, summary(mod.1$gam)$pTerms.table[3])
  nasal.diversity.coeffs.IF <- rbind(nasal.diversity.coeffs.IF, mod.1$gam$coefficients)
  
  print(gene)
}

nasal.infection.diversity <- cbind(nasal.diversity.list, nasal.diversity.p.vals.IF, nasal.diversity.IRIS.p.IF)
colnames(nasal.infection.diversity) <- c("measurement", "p-value", "IRIS_p")
nasal.infection.diversity


#loop start here
stool.gene.names.infection <- colnames(select(stool.otu, -Subject))[colMeans(select(stool.otu, -Subject)) > 0.0001]

stool.gene.names.infection <- setdiff(stool.gene.names.infection, "ASV1172")
stool.gene.names.infection <- setdiff(stool.gene.names.infection, "ASV1601")
stool.gene.names.infection <- setdiff(stool.gene.names.infection, "ASV1992")
stool.gene.names.infection <- setdiff(stool.gene.names.infection, "ASV2706")
stool.gene.names.infection <- setdiff(stool.gene.names.infection, "ASV4396")
#stool.gene.names.infection <- setdiff(stool.gene.names.infection, "ASV4673")

stool.p.vals.IF = c()
stool.IRIS.p.IF = c()
stool.coeffs.IF = c()

for (gene in stool.gene.names.infection) {
  path <- paste0("~/Desktop/infection/stool/ASV_", gene, ".png")
  mod.1 <- gamm(as.formula(paste(gene, " ~ IRIS + s(Time.x, bs = \"cc\", k=5)")),
                data = infection.z.stool, 
                method = "REML",
                control = ctrl,
                random = list(event = ~ 1))
  
  png(path)
  
  plot(mod.1$gam, shade = T, shade.col = "#DF8F44FF", xlab = "Infection Stage", ylab = get_stool_genus(gene), cex.lab = 1.5, font.lab = 2, cex.axis = 1.5, font = 2, main = get_stool_genus(gene), font.main = 2, cex.main = 2.5)
  
  dev.off()
  
  stool.p.vals.IF = rbind(stool.p.vals.IF, summary(mod.1$gam)$s.table[4])
  stool.IRIS.p.IF = rbind(stool.IRIS.p.IF, summary(mod.1$gam)$pTerms.table[3])
  stool.coeffs.IF = rbind(stool.coeffs.IF, mod.1$gam$coefficients)
  
  print(gene)
}

table(stool.p.vals.IF < 0.05)
table(stool.IRIS.p.IF< 0.05)
stool.coeffs.IF

stool.genus.infection.p <- cbind(stool.gene.names.infection,stool.p.vals.IF,stool.IRIS.p.IF) %>% data.frame()
colnames(stool.genus.infection.p) <- c("ASV", "infection_P", "IRIS_P")
stool.genus.infection.p$bodysite <- "Stool"
stool.genus.infection.p$Genus <- stool.taxa$Genus[match(stool.genus.infection.p$ASV, stool.taxa$ASV)]

stool.genus.infection.p[stool.genus.infection.p$infection_P < 0.05,]

# Skin loop starts here
skin.gene.names.infection <- colnames(select(skin.otu, -Subject))[colMeans(select(skin.otu, -Subject)) > 0.0001]

skin.gene.names.infection <- setdiff(skin.gene.names.infection, "OTU_70")
skin.gene.names.infection <- setdiff(skin.gene.names.infection, "OTU_154")
skin.gene.names.infection <- setdiff(skin.gene.names.infection, "OTU_408")
skin.gene.names.infection <- setdiff(skin.gene.names.infection, "OTU_411")
skin.gene.names.infection <- setdiff(skin.gene.names.infection, "OTU_764")
skin.gene.names.infection <- setdiff(skin.gene.names.infection, "OTU_1674")
skin.gene.names.infection <- setdiff(skin.gene.names.infection, "OTU_2608")
skin.gene.names.infection <- setdiff(skin.gene.names.infection, "OTU_2697")
skin.gene.names.infection <- setdiff(skin.gene.names.infection, "OTU_5256")


skin.p.vals.IF <- c()
skin.IRIS.p.IF <- c()
skin.coeffs.IF <- c()

for (gene in skin.gene.names.infection) {
  path <- paste0("~/Desktop/infection/skin/ASV_", gene, ".png")
  mod.1 <- gamm(as.formula(paste(gene, " ~ IRIS + s(Time.x, bs = \"cc\", k=5)")),
                data = infection.z.skin,
                method = "REML",
                control = ctrl,
                random = list(event = ~ 1))
  
  png(path)
  
  plot(mod.1$gam, shade = T, shade.col = "#00A1D5FF", xlab = "Infection Stage", ylab = get_skin_genus(gene), cex.lab = 1.5, font.lab = 2, cex.axis = 1.5, font = 2, main = get_skin_genus(gene), font.main = 2, cex.main = 2.5)
  
  dev.off()
  
  skin.p.vals.IF <- rbind(skin.p.vals.IF, summary(mod.1$gam)$s.table[4])
  skin.IRIS.p.IF <- rbind(skin.IRIS.p.IF, summary(mod.1$gam)$pTerms.table[3])
  skin.coeffs.IF <- rbind(skin.coeffs.IF, mod.1$gam$coefficients)
  
  print(gene)
}

table(skin.p.vals.IF < 0.05)
table(skin.IRIS.p.IF < 0.05)
skin.coeffs.IF

skin.genus.infection.p <- cbind(skin.gene.names.infection, skin.p.vals.IF, skin.IRIS.p.IF) %>% data.frame()
colnames(skin.genus.infection.p) <- c("ASV", "infection_P", "IRIS_P")
skin.genus.infection.p$bodysite <- "Skin"
skin.genus.infection.p$Genus <- skin.taxa$Genus[match(skin.genus.infection.p$ASV, skin.taxa$ASV)]

# Oral loop starts here
oral.gene.names.infection <- colnames(select(oral.otu, -Subject))[colMeans(select(oral.otu, -Subject)) > 0.0001]

oral.gene.names.infection <- setdiff(oral.gene.names.infection, "OTU_40")

oral.p.vals.IF <- c()
oral.IRIS.p.IF <- c()
oral.coeffs.IF <- c()

for (gene in oral.gene.names.infection) {
  path <- paste0("~/Desktop/infection/oral/ASV_", gene, ".png")
  mod.1 <- gamm(as.formula(paste(gene, " ~ IRIS + s(Time.x, bs = \"cc\", k=5)")),
                data = infection.z.oral,
                method = "REML",
                control = ctrl,
                random = list(event = ~ 1))
  
  png(path)
  
  plot(mod.1$gam, shade = T, shade.col = "#B24745FF", xlab = "Infection Stage", ylab = get_oral_genus(gene), cex.lab = 1.5, font.lab = 2, cex.axis = 1.5, font = 2, main = get_oral_genus(gene), font.main = 2, cex.main = 2.5)
  
  dev.off()
  
  oral.p.vals.IF <- rbind(oral.p.vals.IF, summary(mod.1$gam)$s.table[4])
  oral.IRIS.p.IF <- rbind(oral.IRIS.p.IF, summary(mod.1$gam)$pTerms.table[3])
  oral.coeffs.IF <- rbind(oral.coeffs.IF, mod.1$gam$coefficients)
  
  print(gene)
}

table(oral.p.vals.IF < 0.05)
table(oral.IRIS.p.IF < 0.05)
oral.coeffs.IF

oral.genus.infection.p <- cbind(oral.gene.names.infection, oral.p.vals.IF, oral.IRIS.p.IF) %>% data.frame()
colnames(oral.genus.infection.p) <- c("ASV", "infection_P", "IRIS_P")
oral.genus.infection.p$bodysite <- "Oral"
oral.genus.infection.p$Genus <- oral.taxa$Genus[match(oral.genus.infection.p$ASV, oral.taxa$ASV)]

# Nasal loop starts here
nasal.gene.names.infection <- colnames(select(nasal.otu, -Subject))[colMeans(select(nasal.otu, -Subject)) > 0.0001]

nasal.gene.names.infection <- setdiff(nasal.gene.names.infection, "ASV4316")
nasal.gene.names.infection <- setdiff(nasal.gene.names.infection, "ASV4557")
nasal.gene.names.infection <- setdiff(nasal.gene.names.infection, "ASV5718")

nasal.p.vals.IF <- c()
nasal.IRIS.p.IF <- c()
nasal.coeffs.IF <- c()

for (gene in nasal.gene.names.infection) {
  path <- paste0("~/Desktop/infection/nasal/ASV_", gene, ".png")
  mod.1 <- gamm(as.formula(paste(gene, " ~ IRIS + s(Time.x, bs = \"cc\", k=5)")),
                data = infection.z.nasal,
                method = "REML",
                control = ctrl,
                random = list(event = ~ 1))
  
  png(path)
  
  plot(mod.1$gam, shade = T, shade.col = "#79AF97FF", xlab = "Infection Stage", ylab = get_nasal_genus(gene), cex.lab = 1.5, font.lab = 2, cex.axis = 1.5, font = 2, main = get_nasal_genus(gene), font.main = 2, cex.main = 2.5)
  
  dev.off()
  
  nasal.p.vals.IF <- rbind(nasal.p.vals.IF, summary(mod.1$gam)$s.table[4])
  nasal.IRIS.p.IF <- rbind(nasal.IRIS.p.IF, summary(mod.1$gam)$pTerms.table[3])
  nasal.coeffs.IF <- rbind(nasal.coeffs.IF, mod.1$gam$coefficients)
  
  print(gene)
}

table(nasal.p.vals.IF < 0.05)
table(nasal.IRIS.p.IF < 0.05)
nasal.coeffs.IF

nasal.genus.infection.p <- cbind(nasal.gene.names.infection, nasal.p.vals.IF, nasal.IRIS.p.IF) %>% data.frame()
colnames(nasal.genus.infection.p) <- c("ASV", "infection_P", "IRIS_P")
nasal.genus.infection.p$bodysite <- "Nasal"
nasal.genus.infection.p$Genus <- nasal.taxa$Genus[match(nasal.genus.infection.p$ASV, nasal.taxa$ASV)]

stool.genus.infection.p$infection.adj <- p.adjust(stool.genus.infection.p$infection_P, method = "BH", n = length(stool.genus.infection.p$infection_P))
stool.genus.infection.p$IRIS.adj <- p.adjust(stool.genus.infection.p$IRIS_P, method = "BH", n = length(stool.genus.infection.p$IRIS_P))

skin.genus.infection.p$infection.adj <- p.adjust(skin.genus.infection.p$infection_P, method = "BH", n = length(skin.genus.infection.p$infection_P))
skin.genus.infection.p$IRIS.adj <- p.adjust(skin.genus.infection.p$IRIS_P, method = "BH", n = length(skin.genus.infection.p$IRIS_P))

oral.genus.infection.p$infection.adj <- p.adjust(oral.genus.infection.p$infection_P, method = "BH", n = length(oral.genus.infection.p$infection_P))
oral.genus.infection.p$IRIS.adj <- p.adjust(oral.genus.infection.p$IRIS_P, method = "BH", n = length(oral.genus.infection.p$IRIS_P))

nasal.genus.infection.p$infection.adj <- p.adjust(nasal.genus.infection.p$infection_P, method = "BH", n = length(nasal.genus.infection.p$infection_P))
nasal.genus.infection.p$IRIS.adj <- p.adjust(nasal.genus.infection.p$IRIS_P, method = "BH", n = length(nasal.genus.infection.p$IRIS_P))

infection.pvalue <- rbind(stool.genus.infection.p,skin.genus.infection.p,oral.genus.infection.p,nasal.genus.infection.p) %>% filter(infection_P < 0.05)
infection.pvalue 


rbind(stool.genus.infection.p,skin.genus.infection.p,oral.genus.infection.p,nasal.genus.infection.p) %>% filter(IRIS_P < 0.05)

write.csv(file = "~/Desktop/infection/pva.csv",infection.pvalue)


#loop starts here
stool.diversity.p.vals.IF = c()
stool.diversity.IRIS.p.IF = c()
stool.diversity.coeffs.IF = c()

infection.diversity.stool.ir <- filter(infection.diversity.stool, Subject.x %in% ir.list$SubjectID)
infection.diversity.stool.is <- filter(infection.diversity.stool, Subject.x %in% is.list$SubjectID)

infection.diversity.skin.ir <- filter(infection.diversity.skin, Subject.x %in% ir.list$SubjectID)
infection.diversity.skin.is <- filter(infection.diversity.skin, Subject.x %in% is.list$SubjectID)

infection.diversity.oral.ir <- filter(infection.diversity.oral, Subject.x %in% ir.list$SubjectID)
infection.diversity.oral.is <- filter(infection.diversity.oral, Subject.x %in% is.list$SubjectID)

infection.diversity.nasal.ir <- filter(infection.diversity.nasal, Subject.x %in% ir.list$SubjectID)
infection.diversity.nasal.is <- filter(infection.diversity.nasal, Subject.x %in% is.list$SubjectID)


for (gene in stool.diversity.list) {
  path <- paste0("~/Desktop/infection/stool/Diversity_IR", gene, ".png")
  mod.1 <- gamm(as.formula(paste(gene, " ~ s(Time.x, bs = \"cc\", k=5)")),
                data = infection.diversity.stool.ir, 
                method = "REML",
                control = ctrl,
                random = list(event = ~ 1))
  
  png(path)
  
  plot(mod.1$gam, shade = T, shade.col = "#DF8F44FF", xlab = "Infection Stage", ylab = paste(gene), cex.lab = 1.5, font.lab = 2, cex.axis = 1.5, font = 2, main = paste(gene), font.main = 2, cex.main = 2.5)
  
  dev.off()
  
  stool.diversity.p.vals.IF = rbind(stool.diversity.p.vals.IF, summary(mod.1$gam)$s.table[4])
  stool.diversity.IRIS.p.IF = rbind(stool.diversity.IRIS.p.IF, summary(mod.1$gam)$pTerms.table[3])
  stool.diversity.coeffs.IF = rbind(stool.diversity.coeffs.IF, mod.1$gam$coefficients)
  
  print(gene)
}

stool.infection.diversity <- cbind(stool.diversity.list,stool.diversity.p.vals.IF,stool.diversity.IRIS.p.IF)
colnames(stool.infection.diversity) <- c("measurement", "p-value")
stool.infection.diversity

summary(mod.1$gam)

stool.diversity.p.vals.IF = c()
stool.diversity.IRIS.p.IF = c()
stool.diversity.coeffs.IF = c()

stool.diversity.list <- setdiff(stool.diversity.list, "pielou_e")
for (gene in stool.diversity.list) {
  path <- paste0("~/Desktop/infection/stool/Diversity_IS", gene, ".png")
  mod.1 <- gamm(as.formula(paste(gene, " ~ s(Time.x, bs = \"cc\", k=5)")),
                data = infection.diversity.stool.is, 
                method = "REML",
                control = ctrl,
                random = list(event = ~ 1))
  
  png(path)
  
  plot(mod.1$gam, shade = T, shade.col = "#DF8F44FF", xlab = "Infection Stage", ylab = paste(gene), cex.lab = 1.5, font.lab = 2, cex.axis = 1.5, font = 2, main = paste(gene), font.main = 2, cex.main = 2.5)
  
  dev.off()
  
  stool.diversity.p.vals.IF = rbind(stool.diversity.p.vals.IF, summary(mod.1$gam)$s.table[4])
  stool.diversity.IRIS.p.IF = rbind(stool.diversity.IRIS.p.IF, summary(mod.1$gam)$pTerms.table[3])
  stool.diversity.coeffs.IF = rbind(stool.diversity.coeffs.IF, mod.1$gam$coefficients)
  
  print(gene)
}

stool.infection.diversity <- cbind(stool.diversity.list,stool.diversity.p.vals.IF,stool.diversity.IRIS.p.IF)
colnames(stool.infection.diversity) <- c("measurement", "p-value")
stool.infection.diversity

summary(mod.1$gam)

#seperate IR IS: Stool infection
infection.diversity.stool.ir <- filter(infection.diversity.z.stool, Subject.x %in% ir.list$SubjectID)
infection.diversity.stool.is <- filter(infection.diversity.z.stool, Subject.x %in% is.list$SubjectID)

IS.Shannon <- gamm(Shannon ~ s(Time.x, bs = "cc", k=5),
                   data = infection.diversity.stool.is, 
                   method = "REML",
                   control = ctrl,
                   random = list(event = ~ 1))

summary(IS.Shannon$gam)

pdf("~/Desktop/infection/Diversity_IRIS/Stool.IS.Shannon.pdf", width = 5, height = 4)
plot(IS.Shannon$gam, shade = T, shade.col = "#DF8F44FF", xlab = "Infection Stage", ylab = "Diversity", cex.lab = 1.5, font.lab = 2, cex.axis = 1.5, font = 2, main = "IS", font.main = 2, cex.main = 2.5)
dev.off()

IR.Shannon <- gamm(Shannon ~ s(Time.x, bs = "cc", k=5),
                   data = infection.diversity.stool.ir, 
                   method = "REML",
                   control = ctrl,
                   random = list(event = ~ 1))

summary(IR.Shannon$gam)

pdf("~/Desktop/infection/Diversity_IRIS/Stool.IR.Shannon.pdf", width = 5, height = 4)
plot(IR.Shannon$gam, shade = T, shade.col = "#DF8F44FF", xlab = "Infection Stage", ylab = "Diversity", cex.lab = 1.5, font.lab = 2, cex.axis = 1.5, font = 2, main = "IR", font.main = 2, cex.main = 2.5)
dev.off()
anova.gam(IS.Shannon$gam,IR.Shannon$gam)


#seperate IR IS: Skin infection
infection.diversity.skin.ir <- filter(infection.diversity.z.skin, Subject.x %in% ir.list$SubjectID)
infection.diversity.skin.is <- filter(infection.diversity.z.skin, Subject.x %in% is.list$SubjectID)

IS.Shannon <- gamm(Simpson ~ s(Time.x, bs = "cc", k=5),
                   data = infection.diversity.skin.is, 
                   method = "REML",
                   control = ctrl,
                   random = list(event = ~ 1))

summary(IS.Shannon$gam)
pdf("~/Desktop/infection/Diversity_IRIS/Skin.IS.Shannon.pdf", width = 5, height = 4)
plot(IS.Shannon$gam, shade = T, shade.col = "#00A1D5FF", xlab = "Infection Stage", ylab = "Diversity", cex.lab = 1.5, font.lab = 2, cex.axis = 1.5, font = 2, main = "IS", font.main = 2, cex.main = 2.5)
dev.off()

IR.Shannon <- gamm(Simpson ~ s(Time.x, bs = "cc", k=5),
                   data = infection.diversity.skin.ir, 
                   method = "REML",
                   control = ctrl,
                   random = list(event = ~ 1))

summary(IR.Shannon$gam)
pdf("~/Desktop/infection/Diversity_IRIS/Skin.IR.Shannon.pdf", width = 5, height = 4)
plot(IR.Shannon$gam, shade = T, shade.col = "#00A1D5FF", xlab = "Infection Stage", ylab = "Diversity", cex.lab = 1.5, font.lab = 2, cex.axis = 1.5, font = 2, main = "IR", font.main = 2, cex.main = 2.5)
dev.off()

anova.gam(IS.Shannon$gam,IR.Shannon$gam)


#seperate IR IS: oral microbiome
infection.diversity.oral.ir <- filter(infection.diversity.z.oral, Subject.x %in% ir.list$SubjectID)
infection.diversity.oral.is <- filter(infection.diversity.z.oral, Subject.x %in% is.list$SubjectID)

IS.Shannon <- gamm(Shannon ~ s(Time.x, bs = "cc", k=5),
                   data = infection.diversity.oral.is, 
                   method = "REML",
                   control = ctrl,
                   random = list(event = ~ 1))

summary(IS.Shannon$gam)

pdf("~/Desktop/infection/Diversity_IRIS/Oral.IS.Shannon.pdf", width = 5, height = 4)
plot(IS.Shannon$gam, shade = T, shade.col = "#B24745FF", xlab = "Infection Stage", ylab = "Diversity", cex.lab = 1.5, font.lab = 2, cex.axis = 1.5, font = 2, main = "IS", font.main = 2, cex.main = 2.5)
dev.off()

IR.Shannon <- gamm(Shannon ~ s(Time.x, bs = "cc", k=5),
                   data = infection.diversity.oral.ir, 
                   method = "REML",
                   control = ctrl,
                   random = list(event = ~ 1))

summary(IR.Shannon$gam)

pdf("~/Desktop/infection/Diversity_IRIS/Oral.IR.Shannon.pdf", width = 5, height = 4)
plot(IR.Shannon$gam, shade = T, shade.col = "#B24745FF", xlab = "Infection Stage", ylab = "Diversity", cex.lab = 1.5, font.lab = 2, cex.axis = 1.5, font = 2, main = "IR", font.main = 2, cex.main = 2.5)
dev.off()

anova.gam(IS.Shannon$gam,IR.Shannon$gam)



#seperate IR IS: Nasal Microbiome
infection.diversity.nasal.ir <- filter(infection.diversity.z.nasal, Subject.x %in% ir.list$SubjectID)
infection.diversity.nasal.is <- filter(infection.diversity.z.nasal, Subject.x %in% is.list$SubjectID)

IS.Shannon <- gamm(Shannon~ s(Time.x, bs = "cc", k=5),
                   data = infection.diversity.nasal.is, 
                   method = "REML",
                   control = ctrl,
                   random = list(event = ~ 1))

summary(IS.Shannon$gam)

pdf("~/Desktop/infection/Diversity_IRIS/Nasal.IS.Shannon.pdf", width = 5, height = 4)
plot(IS.Shannon$gam, shade = T, shade.col = "#79AF97FF", xlab = "Infection Stage", ylab = "Diversity", cex.lab = 1.5, font.lab = 2, cex.axis = 1.5, font = 2, main = "IS", font.main = 2, cex.main = 2.5)
dev.off()

IR.Shannon <- gamm(Shannon ~ s(Time.x, bs = "cr", k=5),
                   data = infection.diversity.nasal.ir, 
                   method = "REML",
                   control = ctrl,
                   random = list(event = ~ 1))

summary(IR.Shannon$gam)

pdf("~/Desktop/infection/Diversity_IRIS/Nasal.IR.Shannon.pdf", width = 5, height = 4)
plot(IR.Shannon$gam, shade = T, shade.col = "#79AF97FF", xlab = "Infection Stage", ylab = "Diversity", cex.lab = 1.5, font.lab = 2, cex.axis = 1.5, font = 2, main = "IR", font.main = 2, cex.main = 2.5)
dev.off()

anova.gam(IS.Shannon$gam,IR.Shannon$gam)



#######perform same analysis for evenness
IS.Pielou <- gamm(pielou_e ~ s(Time.x, bs = "cc", k=5),
                  data = infection.diversity.stool.is, 
                  method = "REML",
                  control = ctrl,
                  random = list(event = ~ 1))

summary(IS.Pielou$gam)

pdf("~/Desktop/infection/evenness_IRIS/Stool.IS.Pielou.pdf", width = 5, height = 4)
plot(IS.Pielou$gam, shade = T, shade.col = "#DF8F44FF", xlab = "Infection Stage", ylab = "Evenness", cex.lab = 1.5, font.lab = 2, cex.axis = 1.5, font = 2, main = "IS", font.main = 2, cex.main = 2.5)
dev.off()

IR.Pielou <- gamm(pielou_e ~ s(Time.x, bs = "cc", k=5),
                  data = infection.diversity.stool.ir, 
                  method = "REML",
                  control = ctrl,
                  random = list(event = ~ 1))

summary(IR.Pielou$gam)

pdf("~/Desktop/infection/evenness_IRIS/Stool.IR.Pielou.pdf", width = 5, height = 4)
plot(IR.Pielou$gam, shade = T, shade.col = "#DF8F44FF", xlab = "Infection Stage", ylab = "Evenness", cex.lab = 1.5, font.lab = 2, cex.axis = 1.5, font = 2, main = "IR", font.main = 2, cex.main = 2.5)
dev.off()

anova.gam(IS.Pielou$gam,IR.Pielou$gam)

###skin evenness
IS.Shannon <- gamm(pielou_e ~ s(Time.x, bs = "cc", k=5),
                   data = infection.diversity.skin.is, 
                   method = "REML",
                   control = ctrl,
                   random = list(event = ~ 1))

summary(IS.Shannon$gam)

pdf("~/Desktop/infection/evenness_IRIS/skin.IS.pielou_e.pdf", width = 5, height = 4)
plot(IS.Shannon$gam, shade = T, shade.col = "#00A1D5FF", xlab = "Infection Stage", ylab = "Evenness", cex.lab = 1.5, font.lab = 2, cex.axis = 1.5, font = 2, main = "IS", font.main = 2, cex.main = 2.5)
dev.off()

IR.Shannon <- gamm(pielou_e ~ s(Time.x, bs = "cc", k=5),
                   data = infection.diversity.skin.ir, 
                   method = "REML",
                   control = ctrl,
                   random = list(event = ~ 1))

summary(IR.Shannon$gam)

pdf("~/Desktop/infection/evenness_IRIS/skin.IR.pielou_e.pdf", width = 5, height = 4)
plot(IR.Shannon$gam, shade = T, shade.col = "#00A1D5FF", xlab = "Infection Stage", ylab = "Evenness", cex.lab = 1.5, font.lab = 2, cex.axis = 1.5, font = 2, main = "IR", font.main = 2, cex.main = 2.5)
dev.off()
anova.gam(IS.Shannon$gam,IR.Shannon$gam)

###oral 
IS.Shannon <- gamm(pielou_e ~ s(Time.x, bs = "cc", k=5),
                   data = infection.diversity.oral.is, 
                   method = "REML",
                   control = ctrl,
                   random = list(event = ~ 1))

summary(IS.Shannon$gam)

pdf("~/Desktop/infection/evenness_IRIS/oral.IS.pielou_e.pdf", width = 5, height = 4)
plot(IS.Shannon$gam, shade = T, shade.col = "#B24745FF", xlab = "Infection Stage", ylab = "Evenness", cex.lab = 1.5, font.lab = 2, cex.axis = 1.5, font = 2, main = "IS", font.main = 2, cex.main = 2.5)
dev.off()

IR.Shannon <- gamm(pielou_e ~ s(Time.x, bs = "cc", k=5),
                   data = infection.diversity.oral.ir, 
                   method = "REML",
                   control = ctrl,
                   random = list(event = ~ 1))

summary(IR.Shannon$gam)

pdf("~/Desktop/infection/evenness_IRIS/oral.IR.pielou_e.pdf", width = 5, height = 4)
plot(IR.Shannon$gam, shade = T, shade.col = "#B24745FF", xlab = "Infection Stage", ylab = "Evenness", cex.lab = 1.5, font.lab = 2, cex.axis = 1.5, font = 2, main = "IR", font.main = 2, cex.main = 2.5)
dev.off()
anova.gam(IS.Shannon$gam,IR.Shannon$gam)

###Nasal
IS.Shannon <- gamm(pielou_e ~ s(Time.x, bs = "cc", k=5),
                   data = infection.diversity.nasal.is, 
                   method = "REML",
                   control = ctrl,
                   random = list(event = ~ 1))

summary(IS.Shannon$gam)

pdf("~/Desktop/infection/evenness_IRIS/nasal.IS.pielou_e.pdf", width = 5, height = 4)
plot(IS.Shannon$gam, shade = T, shade.col = "#79AF97FF", xlab = "Infection Stage", ylab = "Evenness", cex.lab = 1.5, font.lab = 2, cex.axis = 1.5, font = 2, main = "IS", font.main = 2, cex.main = 2.5)
dev.off()

IR.Shannon <- gamm(pielou_e ~ s(Time.x, bs = "cc", k=5),
                   data = infection.diversity.nasal.ir, 
                   method = "REML",
                   control = ctrl,
                   random = list(event = ~ 1))

summary(IR.Shannon$gam)

pdf("~/Desktop/infection/evenness_IRIS/nasal.IR.pielou_e.pdf", width = 5, height = 4)
plot(IR.Shannon$gam, shade = T, shade.col = "#79AF97FF", xlab = "Infection Stage", ylab = "Evenness", cex.lab = 1.5, font.lab = 2, cex.axis = 1.5, font = 2, main = "IR", font.main = 2, cex.main = 2.5)
dev.off()

anova.gam(IS.Shannon$gam,IR.Shannon$gam)


stool.diveristy.gamm1 <- gamm(Shannon ~  IRIS + s(Time.x, bs = "cc", k=5),
     data = filter(infection.diversity.z.stool), 
     method = "REML",
     control = ctrl,
     random = list(event = ~ 1))

stool.diveristy.gamm2 <- gamm(Shannon ~  s(Time.x, bs = "cc", k=5),
                              data = filter(infection.diversity.z.stool, IRIS == "IS"), 
                              method = "REML",
                              control = ctrl,
                              random = list(event = ~ 1))

summary(stool.diveristy.gamm$gam)


stool.diveristy.gamm$lme


AIC(stool.diveristy.gamm1$lme)
AIC(stool.diveristy.gamm2$lme)

stool.diveristy.gamm3 <- gamm(pielou_e ~  IRIS + s(Time.x, bs = "cc", k=5),
                              data = filter(infection.diversity.z.stool), 
                              method = "REML",
                              control = ctrl,
                              random = list(event = ~ 1))

stool.diveristy.gamm4 <- gamm(pielou_e ~  s(Time.x, bs = "cc", k=5),
                              data = filter(infection.diversity.z.stool, IRIS == "IS"), 
                              method = "REML",
                              control = ctrl,
                              random = list(event = ~ 1))

AIC(stool.diveristy.gamm3$lme)
AIC(stool.diveristy.gamm4$lme)


stool.taxa[str_detect(stool.taxa$Genus, "Finego"),]
skin.taxa[str_detect(skin.taxa$Genus, "Finego"),]


infection.meta_filtered$Date <- as.Date(infection.meta_filtered$Date)
event_duration <- tapply(infection.meta_filtered$Date, infection.meta_filtered$event, function(x) diff(range(x)))
max(event_duration)
median(event_duration)
mean(event_duration)

sd(event_duration)
(88+72)/30

length(unique(infection.meta_filtered$event))

length(unique(infection.meta_filtered$Subject))

infection.subject <- table(infection.meta_filtered$Subject) %>% data.frame()
intersect(as.character(infection.subject$Var1), ir.list$SubjectID)
intersect(as.character(infection.subject$Var1), is.list$SubjectID)
32-19

