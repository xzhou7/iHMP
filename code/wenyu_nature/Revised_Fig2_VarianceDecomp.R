

###===Variance decomposition of healthy baseline, using linear mixed effects===
#Importing revised data
load("Box Sync/PHI_Protected/Shared_WZ_HMP/Revision/Data/Revision_MultiOmes.RData")

##Variance decomposition
#First build a custom function with a dataframe containing variables, A1C, SubjectID, CollectionDate and CL4
#Providing the function with ds and n (1:m for variables)
library(lme4)
  myVD <- function(ds, m, trans){
    htb = {}
    for (subject in unique(ds$SubjectID)) {
      hset = ds[(ds$SubjectID == subject & ds$CL4 == "Healthy"), ]
      hset$CollectionDate <- as.Date(hset$CollectionDate, "%m/%d/%y")
      hset = hset[order(hset$CollectionDate),]
      hset$Days = difftime(hset$CollectionDate, hset$CollectionDate[1], units = "days")
      htb = rbind(htb, hset)
    }
    htb[, 1:m] <- apply(htb[, 1:m], 2, function(x) as.numeric(x))
    htb$Days <- as.numeric(htb$Days)
    vd = {}
    for (i in 1:m) {
      if (trans == "Original") {model <- tryCatch(lmer(scale(htb[, i]) ~ 1 + Days + A1C + SSPG + FPG + (1|SubjectID), data = htb, REML = FALSE), error=function(err) NA)}
      if (trans == "Log10") {model <- tryCatch(lmer(scale(log10(htb[, i])) ~ 1 + Days + A1C + SSPG + FPG + (1|SubjectID), data = htb, REML = FALSE), error=function(err) NA)}
      if (trans == "Arcsin") {model <- tryCatch(lmer(scale(asin(sqrt(htb[, i]))) ~ 1 + Days + A1C + SSPG + FPG + (1|SubjectID), data = htb, REML = FALSE), error=function(err) NA)}
      if (!is.na(model)) {
        v.var = as.data.frame(VarCorr(model), comp=c("Variance"))[, 4] #Extracting random effect variances
        v.var <- c(v.var, as.data.frame(anova(model))[, 2]) #Extracting Sum Sq for fixed effects
      } else {
        #for variable that has too many zeros or for any other error situations
        v.var = rep(NA, 6)
      }
      vd <- rbind(vd, v.var)
    }
    rownames(vd) <- colnames(htb)[1:m]
    vd <- as.data.frame(vd)
    colnames(vd) <- c("random_Subject", "random_Residual", "fix_Days", "fix_A1C", "fix_SSPG", "fix_FPG")
    return(vd)
  } #End of function
    
  #Applying function to Clinics
  clinic.sc = merge(clinic.df, sc, by = "SubjectID")
  clinic.sc[, 62:63] <- apply(clinic.sc[, 62:63], 2, function(x) as.numeric(x))
  vd.clinic = myVD(clinic.sc[, c(4:53, 3, 1, 54, 58, 62, 63)], 50, "Original") #variables, A1C, SubjectID, CollectionDate, CL4, SSPG and FPG
  a1c.df = clinic.sc[, c(2,3,62,63)]
  #Cytokines
  ck.sc = merge(ck.df, a1c.df, by = "SampleID")
  vd.ck = myVD(ck.sc[, c(-1, -2, -3)], 62, "Original")
  #Metabolites
  metb.sc = merge(metbcr.df, a1c.df, by = "SampleID")
  vd.metb = myVD(metb.sc[, -1], 724, "Log10")
  #Proteins
  prot.sc = merge(swathprot.PCR.df, a1c.df, by = "SampleID")
  vd.prot = myVD(prot.sc[, -1], 302, "Original")
  #16S_st
  st.sc = merge(st.df, a1c.df, by.x = "HostSampleID", by.y = "SampleID") #740 obs for overlapping IDs
  vd.st = myVD(st.sc[, -1], 96, "Arcsin")
  stKO.sc = merge(stKO.df, a1c.df, by.x = "HostSampleID", by.y = "SampleID") #740 obs for overlapping IDs
  stKO.sc[2:362] <- apply(stKO.sc[2:362], 2, function(x) x/1000)  #Converting to dicimal 
  vd.stKO = myVD(stKO.sc[, -1], 361, "Arcsin")
  #16S_ns
  ns.sc = merge(ns.df, a1c.df, by.x = "HostSampleID", by.y = "SampleID") #777 obs for overlapping IDs
  vd.ns = myVD(ns.sc[, -1], 80, "Arcsin")
  nsKO.sc = merge(nsKO.df, a1c.df, by.x = "HostSampleID", by.y = "SampleID") #735 obs for overlapping IDs
  nsKO.sc[2:358] <- apply(nsKO.sc[2:358], 2, function(x) x/1000)
  vd.nsKO = myVD(nsKO.sc[, -1], 357, "Arcsin")
  #RNAseq
  rna.sc = merge(rnaseq.log.df, a1c.df, by = "SampleID")
  vd.rna = myVD(rna.sc[, -1], 10343, "Original")
  #Combine all vd tables
  vd.comb = rbind(vd.clinic, vd.ck, vd.metb, vd.rna, vd.prot, vd.st, vd.stKO, vd.ns, vd.nsKO)
  vd.comb$dt = c(rep("ClinicLabs", 50), rep("Cytokines", 62), rep("Metabolites", 724), rep("Transcripts", 10343), rep("Proteins", 302), rep("Gut_Microbes", 96), rep("Gut_MicrobialGenes", 361), rep("Nasal_Microbes", 80), rep("Nasal_MicrobialGenes", 357))
  vd.comb$ICC = vd.comb$random_Subject/(vd.comb$random_Subject+vd.comb$random_Residual)*100
  write.table(vd.comb, file = "Desktop/Revised_Fig2A_VarianceDecomposition.txt", sep = "\t") #12375 obs. 7 variables
  save(vd.comb, file = "Box Sync/WenyuEx/Rcodes/Revised_VarianceDecomposition.RData") 

#Visualization
  load("Box Sync/WenyuEx/Rcodes/Revised_VarianceDecomposition.RData")
  
  library(reshape2)
  pd = vd.comb[, 7:8]
  pd$Variables = rownames(pd)
  pd <- melt(pd)
  
  #Boxplot
  ggplot(data=pd, aes(x = dt, y = value)) +
    geom_jitter(data=pd[pd$dt == "Cytokines",], size = 1.6, alpha =0.6, color = "grey50") +
    geom_jitter(data=pd[pd$dt == "ClinicLabs",], size = 1.6, alpha =0.6, color = "grey50") +
    geom_jitter(data=pd[pd$dt == "Metabolites",], size = 0.6, alpha =0.4, color = "grey50") +
    geom_jitter(data=pd[pd$dt == "Transcripts",], size = 0.1, alpha =0.1, color = "grey50") +
    geom_jitter(data=pd[pd$dt == "Proteins",], size = 0.8, alpha =0.5, color = "grey50") +
    geom_jitter(data=pd[pd$dt == "Gut_Microbes",], size = 1.5, alpha =0.5, color = "grey50") +
    geom_jitter(data=pd[pd$dt == "Gut_MicrobialGenes",], size = 0.8, alpha =0.5, color = "grey50") +
    geom_jitter(data=pd[pd$dt == "Nasal_Microbes",], size = 1.5, alpha =0.5, color = "grey50") +
    geom_jitter(data=pd[pd$dt == "Nasal_MicrobialGenes",], size = 0.8, alpha =0.5, color = "grey50") +
    geom_boxplot(color = "red", size = 0.5, outlier.shape = NA, alpha = 0.1) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text.x = element_text(size = 8, angle = 30, hjust = 1), axis.title.y = element_text(size=10)) +
    scale_x_discrete(limits=c("ClinicLabs", "Cytokines", "Metabolites", "Transcripts","Proteins", "Gut_Microbes", "Gut_MicrobialGenes", "Nasal_Microbes", "Nasal_MicrobialGenes")) +
    labs(x = "", y = "ICC of Normalized Variables")
  ggsave("Desktop/Revised_Fig2A_ICC_Healthy.pdf", width =6, height = 3)
  
  #Scatterplot
  ggplot(data=pd, aes(x=All, y = Individual_Mean, color = ICC)) + geom_point(data=pd[pd$ICC < 50, ], alpha = 0.5) + geom_point(data=pd[pd$ICC > 50, ], alpha = 0.5) +
    scale_colour_gradient(low = "grey80", high = "black") + labs(x="IQR across all subjects", y = "Mean of individual's IQR") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text = element_text(size = 10)) + ylim(c(0,2))
  ggsave("Desktop/Revised_Fig2_HealthyIQR.pdf", width = 6, height = 4)

