

###===For revision II: IR/IS correlation using linear mixed models with interaction term
    ##Healthy visits, IR and IS
    #Loading dataset
    load("Box Sync/PHI_Protected/Shared_WZ_HMP/Revision/Data/Revision_MultiOmes.RData") #Revised version shared on 04/18/18
    
    #Build custom fuction, uses linear mixed models to find significance for interaction term
    #Provding a data frame with 1:m variable, SubjectID and CL4
    #Return tables with correlation IS first and IR next
    library(nlme)
    #myHC <- function(ds, m, trans){ #For intra-omes
    myHC <- function(ds, m, n, trans){ #For inter-omes
      h_IRIS = {}
      for (subject in unique(sc$SubjectID)) { 
        if (sc[sc$SubjectID == subject, "IRIS"] != "Unknown") {
          hset = ds[(ds$SubjectID == subject & ds$CL4 == "Healthy"), ]
          h_IRIS <- rbind(h_IRIS, hset)
        }
      }
      h_IRIS$IRIS = sapply(h_IRIS$SubjectID, function(x) sc$IRIS[match(x, sc$SubjectID)])
      h_IRIS$SubjectID <- factor(h_IRIS$SubjectID)
      h_IRIS$IRIS <- factor(h_IRIS$IRIS)
      #Transform
      if (trans == "none") {h_IRIS[, 1:m] <- apply(h_IRIS[, 1:m], 2, function(x) as.numeric(x))} 
      if (trans == "log2") {h_IRIS[, 1:m] <- apply(h_IRIS[, 1:m], 2, function(x) log2(as.numeric(x)))} 
      if (trans == "arcsin") {h_IRIS[, 1:m] <- apply(h_IRIS[, 1:m], 2, function(x) asin(sqrt(as.numeric(x))))}
      if (trans == "arcsin_log2") {
        h_IRIS[, 1:m] <- apply(h_IRIS[, 1:m], 2, function(x) asin(sqrt(as.numeric(x))))
        h_IRIS[, (m+1):n] <- apply(h_IRIS[, (m+1):n], 2, function(x) log2(as.numeric(x)))
      }
      cor_char = pvalue = IRIS.pvalue = c()
      #for (i in 1:(m-1)) { #for intra-omes
        #for (j in (i+1):m) { #for intra-omes
      for (i in 1:m) { #for inter-omes
        for (j in (m+1):n) { #for inter-omes
          model = tryCatch(lme(as.formula(paste(colnames(h_IRIS)[i], paste(colnames(h_IRIS)[j],"*IRIS", sep = ""), sep = "~")),random = ~1|SubjectID, data=h_IRIS, na.action=na.omit), error=function(err) NA)
          p.value = tryCatch(anova(model)[2,4], error=function(err) NA)
          int.term.p = tryCatch(anova(model)[4,4], error=function(err) NA)
          cor_char = c(cor_char, paste(colnames(h_IRIS)[i], colnames(h_IRIS)[j], sep = " vs "))
          pvalue = c(pvalue, p.value)
          IRIS.pvalue = c(IRIS.pvalue, int.term.p)
         }
      }
      tb_ISIR = data.frame(cor_char = cor_char, LMM.p.value = pvalue, Interaction.p.value = IRIS.pvalue)
      return(tb_ISIR)
    }
    
#Intra-omes:
    LMM.ck = myHC(ck.df[, c(4:65, 70, 75)], 62, "none")
    LMM.clinic = myHC(clinic.df[, c(2:28, 30:50, 52, 53, 58)], 49, "none") #No 29, 51
    LMM.metb = myHC(metbcr.df[, c(2:726, 731)], 724, "log2")
    LMM.prot = myHC(swathprot.PCR.df[, c(2:303, 304, 309)], 302, "none")
    LMM.st = myHC(st.df[, c(2:97, 103, 102)], 96, "arcsin")
    
    save(LMM.ck, LMM.clinic, LMM.metb, LMM.prot, LMM.st, file = "Box Sync/WenyuEx/Rcodes/Revised2_IRIS_signTest_IntraOmes.RData")
    
#Inter-omes to gut microbes
    st_metb.df <- merge(st.df[, 1:97], metbcr.df, by.x = "HostSampleID", by.y = "SampleID")
    LMM.st.metb = myHC(st_metb.df[, c(2:821, 822, 827)], 96, 820, "arcsin_log2")  #2:97 st 98:821 met
    
    st_prot.df <- merge(st.df[, 1:97], swathprot.PCR.df, by.x = "HostSampleID", by.y = "SampleID")
    LMM.st.prot = myHC(st_prot.df[, c(2:399, 400, 405)], 96, 398, "arcsin") #2:97 st 98:399 proteins
    
    st_ck.df <- merge(st.df[, 1:97], ck.df[, c(1,4:75)], by.x = "HostSampleID", by.y = "SampleID")
    LMM.st.ck = myHC(st_ck.df[, c(2:159, 164, 169)], 96, 158, "arcsin") #2:97 st 98:159 cytokines
    
    clinic.clean = apply(clinic.df[, c(2:28, 30:50, 52)], 2, function(x) as.numeric(x)) #No col 29 and 51
    clinic.clean = cbind(clinic.clean, clinic.df[, c(1, 53, 58)])
    st_clinic.df <- merge(st.df[, 1:97], clinic.clean, by.x = "HostSampleID", by.y = "SampleID")
    LMM.st.clinic = myHC(st_clinic.df[, c(2:146, 147, 148)], 96, 145, "arcsin")
    
    st_rna.df <- merge(st.df[, 1:97], rnaseq.log.df, by.x = "HostSampleID", by.y = "SampleID")
    LMM.st.rna = myHC(st_rna.df[, c(2:10443, 10444, 10449)], 96, 10442, "arcsin") #2:97 st 98:10443 rna_genes
    
    save(LMM.st.metb, LMM.st.prot, LMM.st.ck, LMM.st.clinic, LMM.st.rna, file = "Box Sync/WenyuEx/Rcodes/Revised2_IRIS_signTest_InterOmes.RData")
