
##===Variance calculation and variance decomposition of healthy baseline===
#Importing data
load("Box Sync/PHI_Protected/Shared_WZ_HMP/OmicsData/HMP_MultiOmics.RData")
load("Box Sync/PHI_Protected/Shared_WZ_HMP/OmicsData/HMP_Microbiome_Perc.RData")

#Preparing data
  metbcr.df <- metb.df[, c(1, c(1:1798)[colnames(metb.df) %in% metb.curated$Compounds_ID], 1793:1798)] #Only 724 curated metabolites
  #Nasal
  pr = apply(r16s.ns.perc.df[, c(2:954)], 2, function(x) sum(x == 0))
  ns.df <- r16s.ns.perc.df[, c(1, c(2:954)[pr < 400],3807:3812)] #Narrowing to 79 taxa levels that were present in more than half sample size
  pr = apply(r16sKO.ns.perc.df[, c(2:6910)], 2, function(x) sum(x < 1))
  nsKO.df <- r16sKO.ns.perc.df[, c(1, c(2:6910)[pr < 400],6911:6916)] #Narrowing to 357 genes that have more than 1% in more than half sample size
  #Stool
  pr = apply(r16s.st.perc.df[, c(2:322)], 2, function(x) sum(x == 0))
  st.df <- r16s.st.perc.df[, c(1, c(2:322)[pr < 430],2276:2281)] #Narrowing to 96 taxa levels
  pr = apply(r16sKO.st.perc.df[, c(2:6910)], 2, function(x) sum(x < 1))
  stKO.df <- r16sKO.st.perc.df[, c(1, c(2:6910)[pr < 400],6911:6916)] #Narrowing to 361 genes that have more than 1% in more than half sample size
  

##Variance calculation using interquantile range (IQR)
  #Build custom function to calculate IQR given a scaled variable
  #Proving ds, 1: m variable, then SubjectID, CL4
  library(Hmisc)
  myIQR <- function(ds, m){
    iqrtb = {}
    hscaled = data.frame(apply(ds[ds$CL4 == "Healthy", 1:m], 2, function(x) scale(as.numeric(x), center = TRUE, scale = TRUE)))
    hscaled$SubjectID = ds[ds$CL4 == "Healthy", "SubjectID"]
    all.iqr = apply(hscaled[, 1:m], 2, function(x) IQR(as.numeric(x), na.rm = TRUE, type = 8))
    all.mean = apply(ds[ds$CL4 == "Healthy", 1:m], 2, function(x) mean(as.numeric(x), na.rm = TRUE))
    iqrtb = rbind(iqrtb, c("All", all.iqr), c("Expression_Mean", all.mean))
    for (subject in sc$SubjectID) {
      subset = hscaled[(hscaled$SubjectID == subject), 1:m]
      sub.iqr = apply(subset, 2, function(x) IQR(as.numeric(x), na.rm = TRUE, type = 8))
      iqrtb = rbind(iqrtb, c(subject, sub.iqr))
    }
    return(data.frame(iqrtb))
  }
   
  #Host
  iqr.ck = myIQR(ck.df[, c(4:65, 70, 75)], 62) 
  iqr.clinic = myIQR(clinic.df[, c(2:53, 58)], 51) 
  iqr.metb = myIQR(metbcr.df[, c(2:726, 731)], 724) 
  iqr.rna = myIQR(rnaseq.log.df[, c(2:10348, 10353)], 10346) 
  iqr.prot = myIQR(swathprot.df[, c(2:307, 309, 314)], 306) 
  #Microbiome
  iqr.st = myIQR(st.df[, c(2:97, 103, 102)], 96) 
  iqr.stKO = myIQR(stKO.df[, c(2:362, 368, 367)], 361) 
  iqr.ns = myIQR(ns.df[, c(2:81, 87, 86)], 80) 
  iqr.nsKO = myIQR(nsKO.df[, c(2:358, 364, 363)], 357) 
  
  iqr.tb = cbind(iqr.ck, iqr.clinic[, -1], iqr.metb[, -1], iqr.st[, -1], iqr.ns[, -1], iqr.rna[, -1], iqr.prot[, -1], iqr.stKO[, -1], iqr.nsKO[, -1])
  iqr.tb = cbind(t(iqr.tb), c("Type", rep("Cytokines", 62), rep("ClinicLabs", 51), rep("Metabolites", 724), rep("Gut_Microbes", 96), rep("Nasal_Microbes", 80), 
                              rep("Transcripts", 10346), rep("Proteins", 306), rep("Gut_MicrobialGenes", 361), rep("Nasal_MicrobialGenes", 357)))
  
  save(iqr.tb, file = "Box Sync/WenyuEx/Rcodes/MultiOmes_HealthyIQR.RData")
  write.table(iqr.tb, file = "Desktop/Fig2C_Table_MultiOmes_HealthyIQR.txt", sep = "\t")
  
  #To calculate quantiles
  quantile(as.numeric(as.character(unlist(iqr.ck[1, -1]))))
  #ck: 0%       25%       50%       75%      100% 
  #0.1437546 0.3187255 0.4861053 0.8720697 1.5254708 

  #IQR distribution boxplot
  library(ggplot2)
  iqr.tb = data.frame(iqr.tb[-1, c(1, 108)])
  colnames(iqr.tb) = c("iqr.value", "Type")
  iqr.tb$iqr.value = as.numeric(as.character(iqr.tb$iqr.value))
  ggplot(data=iqr.tb, aes(x = Type, y = iqr.value))  + geom_boxplot(color = "red", size = 0.8, outlier.shape = NA) +
    geom_jitter(data=iqr.tb[iqr.tb$Type == "Cytokines",], size = 1.6, alpha =0.6, color = "grey50") +
    geom_jitter(data=iqr.tb[iqr.tb$Type == "ClinicLabs",], size = 1.6, alpha =0.6, color = "grey50") +
    geom_jitter(data=iqr.tb[iqr.tb$Type == "Metabolites",], size = 0.8, alpha =0.4, color = "grey50") +
    geom_jitter(data=iqr.tb[iqr.tb$Type == "Transcripts",], size = 0.1, alpha =0.1, color = "grey50") +
    geom_jitter(data=iqr.tb[iqr.tb$Type == "Proteins",], size = 1, alpha =0.5, color = "grey50") +
    geom_jitter(data=iqr.tb[iqr.tb$Type == "Gut_Microbes",], size = 1.5, alpha =0.5, color = "grey50") +
    geom_jitter(data=iqr.tb[iqr.tb$Type == "Gut_MicrobialGenes",], size = 1, alpha =0.5, color = "grey50") +
    geom_jitter(data=iqr.tb[iqr.tb$Type == "Nasal_Microbes",], size = 1.5, alpha =0.5, color = "grey50") +
    geom_jitter(data=iqr.tb[iqr.tb$Type == "Nasal_MicrobialGenes",], size = 1, alpha =0.5, color = "grey50") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text.x = element_text(size = 10, angle = 30, hjust = 1)) +
    scale_x_discrete(limits=c("ClinicLabs", "Cytokines", "Metabolites", "Transcripts","Proteins", "Gut_Microbes", "Gut_MicrobialGenes", "Nasal_Microbes", "Nasal_MicrobialGenes")) +
    labs(x = "", y = "Interquartile Range of Normalized Variables")
  ggsave("Desktop/Fig2A_IQR_Healthy.pdf", width =6, height = 4)
  
  #IQR versus mean plot for supplemental figure
  iqr.tb$Type <- factor(iqr.tb$Type, levels = c("ClinicLabs", "Cytokines", "Metabolites", "Transcripts", "Proteins", "Gut_Microbes", "Gut_MicrobialGenes", "Nasal_Microbes", "Nasal_MicrobialGenes"))
  ggplot(data=iqr.tb, aes(x = log10(Raw_Mean), y = iqr.value)) + geom_point(size = 0.5, alpha =0.5) + facet_wrap(~Type, scales = "free", nrow =3) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"), 
          strip.text = element_text(size=8, face="bold"), axis.text.x = element_text(size = 4)) +
    labs(x = "Mean of Non-normalized Variables (log10))", y = "Interquartile Range of Normalized Variables")
  ggsave("Desktop/Fig2A_suppl_IQR_versus_mean.pdf", width =4, height = 6)

##Variance decomposition
#First build a custom function with a dataframe containing SampleID, variables, A1C, SubjectID, CollectionDate and CL4
  #Providing the function with ds and n (2:n for variables)
  myVD <- function(ds, n, trans){
    htb = {}
    for (subject in unique(ds$SubjectID)) {
      subset = ds[ds$SubjectID == subject, ]
      hset = subset[subset$CL4 == "Healthy", ]
      hset$CollectionDate <- as.Date(hset$CollectionDate, "%m/%d/%y")
      hset = hset[order(hset$CollectionDate),]
      hset$Days = difftime(hset$CollectionDate, hset$CollectionDate[1], units = "days")
      htb = rbind(htb, hset)
    }
    vd = {}
    for (i in 2:n) {
      if (trans == "Original") {model <- tryCatch(lm(as.numeric(htb[,i]) ~ 1 + as.numeric(Days) + as.numeric(A1C) +as.numeric(SSPG) + as.numeric(FPG) + factor(SubjectID), data = htb), error=function(err) NA)}
      if (trans == "Log10") {model <- tryCatch(lm(log10(as.numeric(htb[,i])) ~ 1 + as.numeric(Days) + as.numeric(A1C) +as.numeric(SSPG) + as.numeric(FPG) + factor(SubjectID), data = htb), error=function(err) NA)}
      if (!is.na(model)) {sq = summary(aov(model))[[1]][, 'Sum Sq']} else {sq = 0}
      if (sum(sq) != 0) {
        rr = c(summary(model)$r.squared, sq / sum(sq) * 100 ) 
      } else {
        #for variable that has too many zeros or for any other error situations
        rr = rep(NA, 7)
      }
      rr = c(colnames(htb)[i], rr)
      vd = rbind(vd, rr)
    }
    colnames(vd) <- c("Variable", "Rsq",  "Days", "A1C", "SSPG", "FPG", "Personal", "Other")
    vd_new = data.frame(Variable=vd[,1], 
                        Days=as.numeric(vd[,3]), A1C=as.numeric(vd[,4]), SSPG=as.numeric(vd[,5]), FPG=as.numeric(vd[,6]),
                        Personal=as.numeric(vd[,7]), Other=as.numeric(vd[,8]),  
                        Experiment=as.numeric(vd[,3])+as.numeric(vd[,4])+as.numeric(vd[,5])+as.numeric(vd[,6]))
    return(vd_new)
  }
  
  #Applying function to Clinics
  clinic.sc = merge(clinic.df, sc, by = "SubjectID")
  vd.clinic = myVD(clinic.sc[, c(2, 4:53, 3, 1, 54, 58, 60, 61)], 51, "Original") #SampleID, variables, A1C, SubjectID, CollectionDate, CL4, SSPG and FPG
  a1c.df = clinic.sc[, c(2,3,60,61)]
  #Cytokines
  ck.sc = merge(ck.df, a1c.df, by = "SampleID")
  vd.ck = myVD(ck.sc[, c(-2, -3)], 63, "Original")
  #Metabolites
  metb.sc = merge(metbcr.df, a1c.df, by = "SampleID")
  vd.metb = myVD(metb.sc, 725, "Log10")
  #Proteins
  prot.sc = merge(swathprot.df, a1c.df, by = "SampleID")
  vd.prot = myVD(prot.sc, 307, "Original")
  #16S_st
  st.sc = merge(st.df, clinic.sc[, c(2,3,60,61,54)], by.x = "HostSampleID", by.y = "SampleID") #740 obs for overlapping IDs
  vd.st = myVD(st.sc, 97, "Original")
  stKO.sc = merge(stKO.df, clinic.sc[, c(2,3,60,61,54)], by.x = "HostSampleID", by.y = "SampleID") #740 obs for overlapping IDs
  vd.stKO = myVD(stKO.sc, 362, "Original")
  #16S_ns
  ns.sc = merge(ns.df, clinic.sc[, c(2,3,60,61,54)], by.x = "HostSampleID", by.y = "SampleID") #777 obs for overlapping IDs
  vd.ns = myVD(ns.sc, 80, "Original")
  nsKO.sc = merge(nsKO.df, clinic.sc[, c(2,3,60,61,54)], by.x = "HostSampleID", by.y = "SampleID") #735 obs for overlapping IDs
  vd.nsKO = myVD(nsKO.sc, 358, "Original")
  #RNAseq
  rna.sc = merge(rnaseq.log.df, a1c.df, by = "SampleID")
  vd.rna = myVD(rna.sc, 10344, "Original")
  #Combine all vd tables
  vd.comb = rbind(vd.clinic, vd.ck, vd.metb, vd.rna, vd.prot, vd.st, vd.stKO, vd.ns, vd.nsKO)
  vd.comb$dt = c(rep("ClinicLabs", 50), rep("Cytokines", 62), rep("Metabolites", 724), rep("Transcripts", 10343), rep("Proteins", 306), rep("Gut_Microbes", 96), rep("Gut_MicrobialGenes", 361), rep("Nasal_Microbes", 79), rep("Nasal_MicrobialGenes", 357))
  write.table(vd.comb, file = "Desktop/Fig2A_Table_MultiOmes_VarianceDecomposition.txt", sep = "\t") #11354 obs. 9 variables
  save(vd.comb, file = "Box Sync/WenyuEx/Rcodes/MultiOmes_VarianceDecomposition.RData") 
  
  #Visualization
  library(reshape2)
  pd = melt(vd.comb[, c("Variable", "Personal", "Experiment", "Other", "dt")], id.vars = c("Variable", "dt"))
  pd$variable <- factor(pd$variable, levels = c("Experiment", "Personal", "Other"))
  ggplot(data=pd, aes(x = dt, y = value, color = variable)) + 
    geom_point(data=pd[pd$dt == "Cytokines",], size = 1.6, alpha = 0.6, position=position_jitterdodge(dodge.width=0.6)) +
    geom_point(data=pd[pd$dt == "ClinicLabs",], size = 1.6, alpha = 0.6, position=position_jitterdodge(dodge.width=0.6)) +
    geom_point(data=pd[pd$dt == "Metabolites",], size = 0.4, alpha = 0.3, position=position_jitterdodge(dodge.width=0.6)) +
    geom_point(data=pd[pd$dt == "Transcripts",], size = 0.1, alpha = 0.1, position=position_jitterdodge(dodge.width=0.6)) +
    geom_point(data=pd[pd$dt == "Proteins",], size = 1, alpha = 0.4, position=position_jitterdodge(dodge.width=0.6)) +
    geom_point(data=pd[pd$dt == "Gut_Microbes",], size = 1.5, alpha = 0.5, position=position_jitterdodge(dodge.width=0.6)) +
    geom_point(data=pd[pd$dt == "Gut_MicrobialGenes",], size = 1, alpha = 0.4, position=position_jitterdodge(dodge.width=0.6)) +
    geom_point(data=pd[pd$dt == "Nasal_Microbes",], size = 1.5, alpha = 0.5, position=position_jitterdodge(dodge.width=0.6)) +
    geom_point(data=pd[pd$dt == "Nasal_MicrobialGenes",], size = 1, alpha = 0.4, position=position_jitterdodge(dodge.width=0.6)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text.x = element_text(angle=45, vjust = 0.6, size=10)) +
    scale_x_discrete(limits=c("ClinicLabs", "Cytokines", "Metabolites", "Transcripts", "Proteins", "Gut_Microbes", "Gut_MicrobialGenes", "Nasal_Microbes", "Nasal_MicrobialGenes")) +
    labs(x = "", y = "Variance Contributing Percentage (%)") + scale_color_discrete(name="Explained By") 
    
  ggsave("Desktop/Fig2A_MultiOmes_VarianceDecomp.pdf", width =7, height = 4)
  
  
