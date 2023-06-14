

###===Indiviudal variance at healthy baseline
#Importing revised data
load("Box Sync/PHI_Protected/Shared_WZ_HMP/Revision/Data/Revision_MultiOmes.RData")

  #Converting microbial KO tables
  stKO.df[2:362] <- apply(stKO.df[2:362], 2, function(x) x/1000)  #Converting to dicimal 
  nsKO.df[2:358] <- apply(nsKO.df[2:358], 2, function(x) x/1000)
  
##Variance calculation using interquantile range (IQR)
#Build custom function to calculate IQR given a scaled variable
#Proving ds, 1: m variable, then SubjectID, CL4
  library(Hmisc)
  myIQR <- function(ds, m, trans){
    iqrtb = {}
    if (trans == "none") {hscaled = data.frame(apply(ds[ds$CL4 == "Healthy", 1:m], 2, function(x) scale(as.numeric(x), center = TRUE, scale = TRUE)))}
    if (trans == "log10") {hscaled = data.frame(apply(ds[ds$CL4 == "Healthy", 1:m], 2, function(x) scale(log10(as.numeric(x)), center = TRUE, scale = TRUE)))}
    if (trans == "arcsin") {hscaled = data.frame(apply(ds[ds$CL4 == "Healthy", 1:m], 2, function(x) scale(asin(sqrt(as.numeric(x))), center = TRUE, scale = TRUE)))}
    hscaled$SubjectID = ds[ds$CL4 == "Healthy", "SubjectID"]
    all.iqr = apply(hscaled[, 1:m], 2, function(x) IQR(as.numeric(x), na.rm = TRUE, type = 8))
    all.mean = apply(ds[ds$CL4 == "Healthy", 1:m], 2, function(x) mean(as.numeric(x), na.rm = TRUE))
    iqrtb = rbind(iqrtb, c("All", all.iqr), c("Expression_Mean", all.mean))
    for (subject in sc$SubjectID) {
      subset = hscaled[(hscaled$SubjectID == subject), 1:m]
      sub.iqr = apply(subset, 2, function(x) IQR(as.numeric(x), na.rm = TRUE, type = 8))
      iqrtb = rbind(iqrtb, c(subject, sub.iqr))
    }
    iqrtb[iqrtb == 0] <- NA
    n = nrow(iqrtb)
    ind.mean = apply(iqrtb[3:n, 2:(m+1)], 2, function(x) mean(as.numeric(as.character(x)), na.rm = TRUE))
    ind.sd = apply(iqrtb[3:n, 2:(m+1)], 2, function(x) sd(as.numeric(as.character(x)), na.rm = TRUE))
    ind.n = ind.who = c()
    for (i in 1:m){
      col.num = as.numeric(as.character(iqrtb[3:n, i+1]))
      col.eva = rowSums(cbind(col.num > (ind.mean[i]+2*ind.sd[i]), col.num < (ind.mean[i]-2*ind.sd[i])))
      col.eva[is.na(col.eva)] <- 0
      col.n = sum(col.eva, na.rm = TRUE)
      col.who = iqrtb[3:n, 1][as.logical(col.eva)]
      ind.n = c(ind.n, col.n)
      ind.who = c(ind.who, paste(col.who, collapse = ","))
    }
    iqrtb <- rbind(iqrtb, c("Individual_Mean", ind.mean))
    iqrtb <- rbind(iqrtb, c("Individual_SD", ind.sd))
    iqrtb <- rbind(iqrtb, c("Num_Outlier_byIQR", ind.n))
    iqrtb <- rbind(iqrtb, c("Outlier_Subject_byIQR", ind.who))
    return(data.frame(iqrtb))
  }
  
  #Host
  iqr.ck = myIQR(ck.df[, c(4:65, 70, 75)], 62, "none") 
  iqr.clinic = myIQR(clinic.df[, c(2:53, 58)], 51, "none") 
  iqr.metb = myIQR(metbcr.df[, c(2:726, 731)], 724, "log10") 
  iqr.rna = myIQR(rnaseq.log.df[, c(2:10348, 10353)], 10346, "none") 
  iqr.prot = myIQR(swathprot.PCR.df[, c(2:303, 304, 309)], 302, "none") 
  #Microbiome
  iqr.st = myIQR(st.df[, c(2:97, 103, 102)], 96, "arcsin") 
  iqr.stKO = myIQR(stKO.df[, c(2:362, 368, 367)], 361, "arcsin") 
  iqr.ns = myIQR(ns.df[, c(2:81, 87, 86)], 80, "arcsin") 
  iqr.nsKO = myIQR(nsKO.df[, c(2:358, 364, 363)], 357, "arcsin") 
  
  iqr.tb = cbind(iqr.ck, iqr.clinic[, -1], iqr.metb[, -1], iqr.st[, -1], iqr.ns[, -1], iqr.rna[, -1], iqr.prot[, -1], iqr.stKO[, -1], iqr.nsKO[, -1])
  iqr.tb = cbind(t(iqr.tb), c("Type", rep("Cytokines", 62), rep("ClinicLabs", 51), rep("Metabolites", 724), rep("Gut_Microbes", 96), rep("Nasal_Microbes", 80), 
                              rep("Transcripts", 10346), rep("Proteins", 302), rep("Gut_MicrobialGenes", 361), rep("Nasal_MicrobialGenes", 357)))
  
  save(iqr.tb, file = "Box Sync/WenyuEx/Rcodes/Revised_MultiOmes_HealthyIQR.RData")
  write.table(iqr.tb, file = "Desktop/Revised_Fig2_Table_MultiOmes_HealthyIQR.txt", sep = "\t")
  
#Visualization
  library(ggplot2)
  library(RColorBrewer)
  
##Density plot for all IQR
  pd = data.frame(IQR=iqr.tb[-1, 1])
  pd$Variable = rownames(iqr.tb)[-1]
  pd$Type = c(rep("Cytokines", 62), rep("ClinicLabs", 51), rep("Metabolites", 724), rep("Gut_Microbes", 96), rep("Nasal_Microbes", 80), 
              rep("Transcripts", 10346), rep("Proteins", 302), rep("Gut_MicrobialGenes", 361), rep("Nasal_MicrobialGenes", 357))
  
  ggplot(data=pd, aes(x = as.numeric(as.character(IQR)), fill = Type, color = Type)) +
    geom_density(data=pd[pd$Type == "Cytokines", ], color = "#e41a1c", size = 0.5, alpha = 0.3) + geom_density(data=pd[pd$Type == "Nasal_Microbes", ], color = "#377eb8", size = 0.5, alpha = 0.3) + geom_density(data=pd[pd$Type == "Gut_MicrobialGenes", ], color = "#ff7f00", size = 0.5, alpha = 0.3) +
    geom_density(data=pd[pd$Type == "Nasal_MicrobialGenes", ], color = "#4daf4a", size = 0.5, alpha = 0.3) + geom_density(data=pd[pd$Type == "Gut_Microbes", ], color = "#984ea3", size = 0.5, alpha = 0.3) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text = element_text(size = 12)) + labs(x="Interquantile Range (IQR)", y="Density", fill="Measurements")
  ggsave("Desktop/Revised_Fig2_HealthyIQR_Bimodality.pdf", width = 6, height = 4)
  
  ggplot(data=pd, aes(x = as.numeric(as.character(IQR)), fill = Type, color = Type)) +
    geom_density(data=pd[pd$Type == "Metabolites", ], color = "#ffff33", size = 0.5, alpha = 0.3) + geom_density(data=pd[pd$Type == "ClinicLabs", ], color = "#a65628", size = 0.5, alpha = 0.3) + 
    geom_density(data=pd[pd$Type == "Transcripts", ], color = "#f781bf", size = 0.5, alpha = 0.3) + geom_density(data=pd[pd$Type == "Proteins", ], color = "#999999", size = 0.5, alpha = 0.3) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text = element_text(size = 12)) + labs(x="Interquantile Range (IQR)", y="Density", fill="Measurements")
  ggsave("Desktop/Revised_Fig2_HealthyIQR_unimodality.pdf", width = 6, height = 4)
  
##Examples showing low (IL6) and high (LEPTIN) IQR
  ds = ck.df[ck.df$CL4 == "Healthy", c("IL6", "TNFB", "IL1A", "RESISTIN", "PAI1", "LEPTIN")]
  ds <- data.frame(apply(ds, 2, function(x) scale(x)))
  pd = melt(ds)
  
  ggplot(data=pd, aes(x = variable, y = value)) + geom_jitter(size = 0.8, alpha =0.8, color = "grey40") + geom_boxplot(color = "red", size = 0.5, outlier.shape = NA, alpha = 0.1) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text.x = element_text(size = 10), axis.title.y = element_text(size=10)) +
   labs(x = "", y = "Normalized Value")
  ggsave("Desktop/Revised_Fig2_IQRexamples.pdf", width =6, height = 4)
  
##Scattor plot for IQR versus expression, color by ICC
  pd = data.frame(IQR=iqr.tb[-1, c(1:2, 112)])
  pd$Variable = rownames(iqr.tb)[-1]
  colnames(pd)[1:3] <- c("IQR", "Expression_mean", "dt")
  pd[, 1:2] <- apply(pd[, 1:2], 2, function(x) as.numeric(as.character(x)))
  
  #Transform for visuliazation
  pd[1:62, 2] <- log2(pd[1:62, 2])  #Cytokines
  pd[114:837, 2] <- log10(pd[114:837, 2])  #Metabolites
  pd <- merge(pd, vd.comb, by.x = "Variable", by.y = "row.names")
  pd$dt.x <- factor(pd$dt.x, levels = c("ClinicLabs", "Cytokines", "Transcripts", "Metabolites", "Proteins", "Gut_Microbes", "Gut_MicrobialGenes", "Nasal_Microbes", "Nasal_MicrobialGenes"))
 
  ggplot(data=pd, aes(x=IQR, y = Expression_mean, color = ICC)) + geom_point(alpha = 0.5) + facet_wrap(~dt.x, nrow = 3, scales = "free") + 
    scale_colour_gradient(low = "grey80", high = "black") + labs(x="IQR across all subjects", y = "Expression Mean") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"), strip.text = element_text(size=12),
          axis.text = element_text(size = 12), strip.text.x = element_text(size = 8))
  ggsave("Desktop/Revised_Fig2_HealthyIQR_facet.pdf", width = 6, height = 6)
  
##Investigating microbes ICC and expression
  pd <- pd[pd$dt.x %in% c("Gut_Microbes", "Gut_MicrobialGenes", "Nasal_Microbes", "Nasal_MicrobialGenes"), ]
  
  ggplot(data=pd, aes(x=Expression_mean, y=ICC)) + geom_point(alpha = 0.5) + facet_wrap(~dt.x, nrow = 2, scales = "free_x") + 
    scale_colour_gradient(low = "grey80", high = "black") + labs(x="Mean Abundance", y = "Inter-class Correlation") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"), strip.text = element_text(size=12),
          axis.text = element_text(size = 8), strip.text.x = element_text(size = 10)) + ylim(c(0, 80))
  ggsave("Desktop/Revised_Fig2_ICC_exp_16s.pdf", width = 4, height = 4)
  
  
  
  
  
  
  
  
  
  
  
  
  
##Other plots
  #Demonstrate the ICC and individual IQR outlier
  #Using MIG, LIF as examples
  ds = ck.df[, c(4:65, 70, 75)]
  m = 62
  hscaled = data.frame(apply(ds[ds$CL4 == "Healthy", 1:m], 2, function(x) scale(as.numeric(x), center = TRUE, scale = TRUE)))
  hscaled$SubjectID = ds[ds$CL4 == "Healthy", "SubjectID"]
  ggplot(data=hscaled, aes(x=factor(SubjectID), y=MIG)) + geom_boxplot(size = 0.3, color = 'red', outlier.shape = NA) + geom_jitter(size = 0.05) + 
    theme_bw() + theme(axis.text.x = element_text(angle = 45, size = 4, hjust = 1)) + labs(x = "SubjectID")
  ggsave("Desktop/PF_HealthyIQR_MIG.pdf", width = 6, height = 4)

  #load RData
  load("Box Sync/WenyuEx/Rcodes/Revised_MultiOmes_HealthyIQR.RData")
  load("Box Sync/WenyuEx/Rcodes/Revised_VarianceDecomposition.RData")
  
  pd = as.data.frame(iqr.tb[, c(1, 108)])
  colnames(pd) <- iqr.tb[1, c(1, 108)]
  pd <- merge(pd[-1, ],vd.comb[, c(7,8)], by = "row.names")
  pd[, 2:3] <- apply(pd[, 2:3], 2, function(x) as.numeric(as.character(x)))
  #Scatterplot
  ggplot(data=pd, aes(x=All, y = Individual_Mean, color = ICC)) + geom_point(data=pd[pd$ICC < 50, ], alpha = 0.5) + geom_point(data=pd[pd$ICC > 50, ], alpha = 0.5) +
    scale_colour_gradient(low = "grey80", high = "black") + labs(x="IQR across all subjects", y = "Mean of individual's IQR") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text = element_text(size = 10)) + ylim(c(0,2))
  ggsave("Desktop/Revised_Fig2_HealthyIQR.pdf", width = 6, height = 4)
  #Facets
  pd$dt <- factor(pd$dt, levels = c("ClinicLabs", "Cytokines", "Transcripts", "Metabolites", "Proteins", "Gut_Microbes", "Gut_MicrobialGenes", "Nasal_Microbes", "Nasal_MicrobialGenes"))
  ggplot(data=pd, aes(x=All, y = Individual_Mean, color = ICC)) + geom_point(alpha = 0.5) + facet_wrap(~dt, nrow = 3) + 
    scale_colour_gradient(low = "grey80", high = "black") + labs(x="IQR across all subjects", y = "Mean of individual's IQR") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"), strip.text = element_text(size=12),
          axis.text = element_text(size = 12)) + ylim(c(0,2))
  ggsave("Desktop/Revised_Fig2_HealthyIQR_facet.pdf", width = 6, height = 6)
  
  #Inspecting the ratio of individual/all in relation to expression_mean
  pd = as.data.frame(iqr.tb[-1, c(1:2, 108, 112)])
  colnames(pd) <- iqr.tb[1, c(1:2, 108, 112)]
  pd[, 1:3] <- apply(pd[, 1:3], 2, function(x) as.numeric(as.character(x)))
  pd$ratio = pd[, 3]/pd[, 1]
  
  #Facets
  pd$Type <- factor(pd$Type, levels = c("ClinicLabs", "Cytokines", "Transcripts", "Metabolites", "Proteins", "Gut_Microbes", "Gut_MicrobialGenes", "Nasal_Microbes", "Nasal_MicrobialGenes"))
  ggplot(data=pd[(pd$ratio < 20 & pd$ratio >0 & !is.na(pd$Type)), ], aes(x=ratio, y = Expression_Mean)) + geom_point(alpha = 0.5) + facet_wrap(~Type, nrow = 4, scales = "free") + 
    labs(x="IQR ratio individual mean/all subjects", y = "Expression Mean") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"), 
          strip.text = element_text(size=12), axis.text = element_text(size = 12)) 
  ggsave("Desktop/Revised_Fig2_HealthyIQR_ratio_expression_facet.pdf", width = 6, height = 6)
  
  #Boxplot
  pd = data.frame(IQR=iqr.tb[-1, 1])
  pd$Variable = rownames(iqr.tb)[-1]
  pd$Type = c(rep("Cytokines", 62), rep("ClinicLabs", 51), rep("Metabolites", 724), rep("Gut_Microbes", 96), rep("Nasal_Microbes", 80), 
              rep("Transcripts", 10346), rep("Proteins", 302), rep("Gut_MicrobialGenes", 361), rep("Nasal_MicrobialGenes", 357))
  ggplot(data=pd, aes(x = Type, y = as.numeric(as.character(IQR))))  + geom_boxplot(color = "red", size = 0.8, outlier.shape = NA) +
    geom_jitter(data=pd[pd$Type == "Cytokines",], size = 1.6, alpha =0.6, color = "grey50") +
    geom_jitter(data=pd[pd$Type == "ClinicLabs",], size = 1.6, alpha =0.6, color = "grey50") +
    geom_jitter(data=pd[pd$Type == "Metabolites",], size = 0.8, alpha =0.4, color = "grey50") +
    geom_jitter(data=pd[pd$Type == "Transcripts",], size = 0.1, alpha =0.1, color = "grey50") +
    geom_jitter(data=pd[pd$Type == "Proteins",], size = 1, alpha =0.5, color = "grey50") +
    geom_jitter(data=pd[pd$Type == "Gut_Microbes",], size = 1.5, alpha =0.5, color = "grey50") +
    geom_jitter(data=pd[pd$Type == "Gut_MicrobialGenes",], size = 1, alpha =0.5, color = "grey50") +
    geom_jitter(data=pd[pd$Type == "Nasal_Microbes",], size = 1.5, alpha =0.5, color = "grey50") +
    geom_jitter(data=pd[pd$Type == "Nasal_MicrobialGenes",], size = 1, alpha =0.5, color = "grey50") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text.x = element_text(size = 10, angle = 30, hjust = 1)) +
    scale_x_discrete(limits=c("ClinicLabs", "Cytokines", "Metabolites", "Transcripts","Proteins", "Gut_Microbes", "Gut_MicrobialGenes", "Nasal_Microbes", "Nasal_MicrobialGenes")) +
    labs(x = "", y = "ICC of Normalized Variables")
 
  