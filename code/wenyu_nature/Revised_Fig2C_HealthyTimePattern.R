###Revised version at April 20, 2018
###===For finding time associating variables within healthy baselines
##Healthy visits, longitudinal time pattern, by correlation of delta change verus delta time 
#Loading dataset
load("Box Sync/PHI_Protected/Shared_WZ_HMP/Revision/Data/Revision_MultiOmes.RData") #Revised version shared on 04/18/18

##Subjects with lengthy collection > 900 in all dataset
SubjectID.d900 = c("69-001","69-010","69-012","69-022","69-023","69-026","69-027","69-032","69-033","69-034",
                   "69-035","69-036","69-037","69-039","69-045","69-047","69-048","70-1001","70-1002","70-1003",
                   "70-1004","70-1005","70-1006","70-1008","70-1012", "70-1014", "70-1010")  #27 subjects verified

#Build custom fuction (Healthy Ordered, HO) to align time points for individuals having at least THREE longitudinal visits
  #Provding a data frame with SubjectID, CollectionDate, CL4, 1:m variable
  #Return tables with significant(all < p after bonferroni MHC, started with variables with smaller p) variables with the delta values and DiffDay
  library(Hmisc)
  library(rmcorr)
  myHO <- function(ds, m, p){
    ho.nord = {}
    #for (subject in unique(ds$SubjectID)) { #for all
    for (subject in SubjectID.d900) { # for SubjectID that have more than 900 days' collection
        hset = ds[(ds$SubjectID == subject & ds$CL4 == "Healthy"), ]
        if (nrow(hset) > 2 ){
          hset$CollectionDate <- as.Date(hset$CollectionDate, "%m/%d/%y")
          hset = hset[order(hset$CollectionDate), ]
          # Need a model with random intercept and fixed slope here
          T1.Date = hset[, "CollectionDate"][1]
          hset$DiffDay = lapply(hset[, "CollectionDate"], function(x) difftime(x, T1.Date, units = "days"))
          ho.nord = rbind(ho.nord, hset)
      }
    }
    ho.nord[, 1:m] <- apply(ho.nord[, 1:m], 2, function(x) as.numeric(x))
    ho.nord$SubjectID <- factor(ho.nord$SubjectID)
    ho.nord$DiffDay <- as.numeric(ho.nord$DiffDay)
    cor_char = cor_r = cor_p = CI_left = CI_right = c()
    for (i in 1:m) {
      rmc.dset <- tryCatch(rmcorr(participant = SubjectID, measure1 = DiffDay, measure2 = colnames(ho.nord)[i], dataset = ho.nord), error=function(err) NA)
      if (!is.na(rmc.dset)) {
      cor_char <- c(cor_char, colnames(ho.nord)[i])
      cor_r <- c(cor_r, unlist(rmc.dset[1]))
      cor_p <- c(cor_p, unlist(rmc.dset[3]))
      CI_left <- c(CI_left, unlist(rmc.dset[4])[1])
      CI_right <- c(CI_right, unlist(rmc.dset[4])[2])
      }
    }
    cor_time = data.frame(cor_char, cor_r, cor_p, CI_left, CI_right)
    cor_time$p.adj = p.adjust(cor_time$cor_p, method = "fdr")
    cor_time <- cor_time[order(cor_time$p.adj), ]
    sigvar = cor_time[cor_time$p.adj < p & !is.na(cor_time$p.adj), ] # leaving out variables that didn't pass the significance cut-off
    #return(sigvar) 
    return(cor_time)
  } #End of custom function
    #the custom function doesn't work with "rmcorr" inside. Need to excute the following by defining ds, m, p in separate lines first
  
#Provding a data frame with SubjectID, CollectionDate, CL4, and variable in column 1:m 
  #Host
  ho.ck = myHO(ck.df[, c(4:65, 70, 71, 75)], 62, 0.2)  #5 hits for 20% FDR, 2 under 5%. All: 23
  ho.clinic = myHO(clinic.df[, c(2:54, 58)], 51, 0.05) #22 hits for 5% FDR. All: 19
  ho.metb = myHO(metbcr.df[, c(2:727, 731)], 724, 0.05) #191 hits for 5% FDR. All: 158
  ho.prot = myHO(swathprot.PCR.df[, c(2:305, 309)], 302, 0.05) #57 hit for 5% FDR. All: 75
  ho.rna = myHO(rnaseq.log.df[, c(2:10349, 10353)], 10346, 0.2) #none
  #Microbiome
  ho.st = myHO(st.df[, c(2:97, 102:104)], 96, 0.05) #27 hits. All: 34
  ho.ns = myHO(ns.df[, c(2:81, 86:88)], 80, 0.05) #10 hits. All: 11
  ho.stKO = myHO(stKO.df[, c(2:362, 367:369)], 361, 0.05) #1 hits. All: 1
  ho.nsKO = myHO(nsKO.df[, c(2:358, 363:365)], 357, 0.05) #42 hits. All: 22
  #When returning rho and p.adj tables, 355 variables for D900, 343 variables for all
  write.table(rbind(ho.ck, ho.clinic, ho.metb, ho.prot, ho.st, ho.ns, ho.stKO, ho.nsKO), file = "Desktop/Revised_Fig2C_TimeAssociatingVariables.txt", sep = "\t", row.names = FALSE)
  
  #Merge tables
  cor_D900 = read.table("Box Sync/WenyuEx/Research_HMP/Results_Tables/Revision/Revised_Fig2B_D900_TimeAssociatingVariables.txt", header = TRUE, sep = "\t")
  cor_D900$DataType = c(rep("Cytokines", 5), rep("ClinicLabs", 22), rep("Metabolites", 191), rep("Proteins", 57), rep("GutMicrobes", 27), rep("NasalMicrobes", 10), rep("GutKOs", 1), rep("NasalKO", 42))
  cor_all = read.table("Box Sync/WenyuEx/Research_HMP/Results_Tables/Revision/Revised_Fig2B_ALL_TimeAssociatingVariables.txt", header = TRUE, sep = "\t")
  cor_all$DataType = c(rep("Cytokines", 23), rep("ClinicLabs", 19), rep("Metabolites", 158), rep("Proteins", 75), rep("GutMicrobes", 34), rep("NasalMicrobes", 11), rep("GutKOs", 1), rep("NasalKO", 22))
  
  cor_comb = merge(cor_all, cor_D900, by = "cor_char", all = TRUE)
  write.table(cor_comb, file = "Desktop/Revised_Fig2B_Combined.txt", sep = "\t", row.names = FALSE)
  
  
#Visualization
  #For ALKP
  char = "ALKP"
  rmc.dset <- rmcorr(participant = SubjectID, measure1 = DiffDay, measure2 = char, dataset = ho.nord)
  plot(rmc.dset, xlab = "Days since the onset", ylab = paste(char), cex = 0.5)
  
  #For nRPLC_395.1897_10, Pregnenolone sulfate, with r -0.06826777, p-value 0.1450883, 95% CI -0.1591974 0.02381033
  char = "nRPLC_395.1897_10" 
  rmc.dset <- rmcorr(participant = SubjectID, measure1 = DiffDay, measure2 = char, dataset = ho.nord)
  plot(rmc.dset, xlab = "Days since the onset", ylab = paste(char), cex = 0.5)
  