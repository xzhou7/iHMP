

###===For revision: correlation between individuals, using median values
##Healthy visits
#Loading dataset
load("Box Sync/PHI_Protected/Shared_WZ_HMP/Revision/Data/Revision_MultiOmes.RData") #Revised version shared on 04/18/18

  sc[c(7, 18,29), 10] <- NA

#Build custom fuction, uses rmcorr to find 
#Provding a data frame with 1:m variable, SubjectID and CL4
#Return tables with correlation IS first and IR next
#Turn on/off intra-omes or inter-omes clause, there are seven clauses needed accordingly
  library(ppcor)
  #myHPC <- function(ds, m, trans){ #CLAUSE #1: For intra-omes 
  myHPC <- function(ds, m, n, trans){ #For inter-omes
    h_all = h_IS = h_IR = {}
    for (subject in unique(ds$SubjectID)) { 
      if (subject %in% sc$SubjectID) {
        hset = ds[(ds$SubjectID == subject & ds$CL4 == "Healthy"), ]
        #h_median = apply(hset[, 1:m], 2, function(x) median(as.numeric(x), na.rm = TRUE)) #CLAUSE #2:For intra-omes
        h_median = apply(hset[, 1:(n-1)], 2, function(x) median(as.numeric(x), na.rm = TRUE)) #for inter-omes
        if ( sc[sc$SubjectID == subject, "IRIS"] == "IS") {
          h_IS <- rbind(h_IS, c(h_median, subject))
        }
        if ( sc[sc$SubjectID == subject, "IRIS"] == "IR") {
          h_IR <- rbind(h_IR, c(h_median, subject))
        }
        h_all <- rbind(h_all, c(h_median, subject))
      }
    }
    #cor_final = as.data.frame(matrix(, nrow = m*(m-1)/2, ncol = 0)) #CLAUSE #3: for intra-omes
    cor_final = as.data.frame(matrix(, nrow = m*(n-1-m), ncol = 0))  #for inter-omes
    for (df in c("h_all", "h_IS", "h_IR")) {
      hdf = as.data.frame(get(df))
      #colnames(hdf)[m+1] <- "SubjectID" #CLAUSE #4: Intra-omes
      colnames(hdf)[n] <- "SubjectID" #Inter-omes
      #Transform
      if (trans == "none") {hdf[, 1:m] <- apply(hdf[, 1:m], 2, function(x) as.numeric(as.character(x)))} 
      if (trans == "log2") {hdf[, 1:m] <- apply(hdf[, 1:m], 2, function(x) log2(as.numeric(as.character(x))))} 
      if (trans == "arcsin") {hdf[, 1:m] <- apply(hdf[, 1:m], 2, function(x) asin(sqrt(as.numeric(as.character(x)))))}
      if (trans == "arcsin_log2") {
        hdf[, 1:m] <- apply(hdf[, 1:m], 2, function(x) asin(sqrt(as.numeric(as.character(x)))))
        hdf[, (m+1):(n-1)] <- apply(hdf[, (m+1):(n-1)], 2, function(x) log2(as.numeric(as.character(x))))
      }
      hdf <- merge(hdf, sc[, c(1, 11:12)], by = "SubjectID") #Categorial variable doesn't work for pcor.test
      cor_char = cor_r = cor_p = c()
      #for (i in 2:m) { #CLAUSE #5: for intra-omes
        #for (j in (i+1):(m+1)) { #for intra-omes
      for (i in 2:(m+1)) { #CLAUSE #6: for inter-omes
        for (j in (m+2):n) { #for inter-omes
              #m.ds = hdf[, c(i, j, (m+2), (m+3))]  #CLAUSE #7: For intra-omes
              m.ds = hdf[, c(i, j, (n+1), (n+2))]  #For inter-omes
              del.row = unique(which(is.na(m.ds), arr.ind = TRUE)[, 1])
              if (length(del.row) != 0) {m.ds <- m.ds[-del.row, ]}  #deleting rows having NAs
              model <- pcor.test(m.ds[, 1], m.ds[, 2], m.ds[, 3:4])
              cor_char <- c(cor_char, paste(colnames(hdf)[i], colnames(hdf)[j], sep = " vs "))
              cor_r <- c(cor_r, model[, 1])
              cor_p <- c(cor_p, model[, 2])
        }
      }
      cor_df = data.frame(cor_char, cor_r, cor_p)
      cor_df$adj.p = p.adjust(cor_p, method = "fdr")
      cor_final = cbind(cor_final, cor_df)
    }
    return(cor_final)
  }
    
  #Intra-omes:
  HPC.ck = myHPC(ck.df[, c(4:65, 70, 75)], 62, "none")
  HPC.clinic = myHPC(clinic.df[, c(2:28, 30:50, 52, 53, 58)], 49, "none")
  HPC.metb = myHPC(metbcr.df[, c(2:726, 731)], 724, "log2")
  HPC.prot = myHPC(swathprot.PCR.df[, c(2:303, 304, 309)], 302, "none")
  HPC.st = myHPC(st.df[, c(2:97, 103, 102)], 96, "arcsin")
  
  #Inter-omes
  st_ck.df <- merge(st.df[, 1:97], ck.df[, c(1,4:75)], by.x = "HostSampleID", by.y = "SampleID")
  HPC.st.ck = myHPC(st_ck.df[, c(2:159, 164, 169)], 96, 159, "arcsin")
  
  clinic.clean = apply(clinic.df[, c(2:28, 30:50, 52)], 2, function(x) as.numeric(x)) #No col 29 and 51
  clinic.clean = cbind(clinic.clean, clinic.df[, c(1, 53, 58)])
  st_clinic.df <- merge(st.df[, 1:97], clinic.clean, by.x = "HostSampleID", by.y = "SampleID")
  HPC.st.clinic = myHPC(st_clinic.df[, c(2:146, 147, 148)], 96, 146, "arcsin") #2:97 st 98:146 clinic labs

  st_metb.df <- merge(st.df[, 1:97], metbcr.df, by.x = "HostSampleID", by.y = "SampleID")
  HPC.st.metb = myHPC(st_metb.df[, c(2:821, 822, 827)], 96, 821,"arcsin_log2")  #2:97 st 98:821 metb
  
  st_prot.df <- merge(st.df[, 1:97], swathprot.PCR.df, by.x = "HostSampleID", by.y = "SampleID")
  HPC.st.prot = myHPC(st_prot.df[, c(2:399, 400, 405)], 96, 399,"arcsin") #2:97 st 98:399 proteins
  
  save(HPC.ck, HPC.clinic, HPC.metb, HPC.prot, HPC.st, HPC.st.ck, HPC.st.clinic, HPC.st.metb, HPC.st.prot, file = "Box Sync/WenyuEx/Rcodes/Revised_InterCor.RData")
  
##Visualization
  #Barplot with error bars, showing the average edge per analyte (epa)
  #Preparing matrix
  epa = ave.epa = ratio = {}
  for (df in c("HPC.st.ck", "HPC.st.clinic", "HPC.st.metb", "HPC.st.prot")){
  #for (df in c("HPC.ck", "HPC.clinic", "HPC.metb", "HPC.prot", "HPC.st", "HPC.st.ck", "HPC.st.clinic", "HPC.st.metb", "HPC.st.prot")){
    ds = get(df)
    #For all
    ds.sig = as.character(ds[(ds[, 4] < 0.05), 1])
    print(paste(c(df, "Combined", length(ds.sig)), collapse = " "))
    epa.all = data.frame(table(unlist(strsplit(ds.sig, split = " vs "))))
    epa.all <- data.frame(epa.all, ISIR=rep("Combined", nrow(epa.all)), Type=rep(df, nrow(epa.all)))
    ave.all = c(df, "Combined", mean(epa.all$Freq), sd(epa.all$Freq))
    #For IS
    ds.sig = as.character(ds[(ds[, 8] < 0.05), 1])
    print(paste(c(df, "IS", length(ds.sig)), collapse = " "))
    epaIS = data.frame(table(unlist(strsplit(ds.sig, split = " vs "))))
    epaIS <- data.frame(epaIS, ISIR=rep("IS", nrow(epaIS)), Type=rep(df, nrow(epaIS)))
    ave.IS = c(df, "IS", mean(epaIS$Freq), sd(epaIS$Freq))
    #For IR
    ds.sig = as.character(ds[(ds[, 12] < 0.05), 1])
    print(paste(c(df, "IR", length(ds.sig)), collapse = " "))
    epaIR = data.frame(table(unlist(strsplit(ds.sig, split = " vs "))))
    epaIR <- data.frame(epaIR, ISIR=rep("IR", nrow(epaIR)), Type=rep(df, nrow(epaIR)))
    ave.IR = c(df, "IR", mean(epaIR$Freq), sd(epaIR$Freq))
    if (nrow(epa.all) == 0 |nrow(epaIS) == 0 |nrow(epaIR) == 0){
      ISIR = {}
    } else {
      ISIR = merge(merge(epaIS[, 1:3], epaIR, by = "Var1", all = TRUE), epa.all[, 1:3], by = "Var1", all = TRUE)
    }
    epa = rbind(epa, epa.all, epaIS, epaIR)
    ave.epa = rbind(ave.epa, ave.all, ave.IS, ave.IR)
  }
  ave.epa <- as.data.frame(ave.epa)
  colnames(ave.epa) <- c("Type", "ISIR", "Freq_mean", "Freq_sd")
  ave.epa[, 3:4] <- apply(ave.epa[, 3:4], 2, function(x) as.numeric(as.character(x))) 
  
  ave.epa$ISIR <- factor(ave.epa$ISIR, levels = c("Combined", "IS", "IR"))
  ggplot(data = ave.epa, aes(x = Type, y = Freq_mean, fill = ISIR)) + geom_bar(position = "dodge", stat = "identity") +
    geom_errorbar(aes(ymin=Freq_mean, ymax=Freq_mean+Freq_sd), width=.1, position=position_dodge(.9)) + 
    scale_fill_manual(breaks = c("Combined", "IS", "IR"), values=c("grey80", "#90AA3C", "#EF6125")) +
    #scale_x_discrete(breaks=c("HPC.st.ck", "HPC.st.clinic", "HPC.st.metb", "HPC.st.prot"),labels=c("GutMicrobes\nCytokines", "GutMicrobe\nClinicLabs", "GutMicrobe\nMetabolites", "GutMicrobes\nProteins")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text.x = element_text(size = 10, angle = 45, hjust = 1)) + 
    labs(x="Omics Association", y="Average number of edges per analyte", fill='')
  ggsave("Desktop/Revised_Fig5_InterCor_all.pdf", height = 4, width = 6)
  
#Comparing with rmcorr results
  load("Box Sync/WenyuEx/Rcodes/Revised_IntraOmes_rmcorr.RData")
  diff.int = diff.rmc = {}
  tp = "clinic"  #for ck: Identify IL2 (col 9) vs MIG (col 50) as one of examples
  df_inter = get(paste("HPC", tp, sep = "."))
  df_rmcorr = get(paste("HC", tp, sep = "."))
  int.ind = cbind(df_inter[, 8] < 0.05, df_inter[, 12] < 0.05)
  rmc.ind = cbind(df_rmcorr[, 6] < 0.05, df_rmcorr[, 12] < 0.05)
  for (i in 1: nrow(df_inter)) {
    if (sum(int.ind[i, ]) ==2 & sum(rmc.ind[i, ]) == 0) {
      diff.int = rbind(diff.int, df_inter[i, ])
      diff.rmc = rbind(diff.rmc, df_rmcorr[i, ])
    }
  }
  
  library(RColorBrewer)
  IS_color = colorRampPalette(brewer.pal(12, 'Paired'))
  IR_color = colorRampPalette(brewer.pal(9, 'Set1'))
  
  ds = clinic.df[, c(2:28, 30:50, 52, 53, 58)]
  m = 49
  #For rmcorr 29 IS with 209 obs, 34 IR with 199 obs.
  h_IS = h_IR = {}
  for (subject in unique(ds$SubjectID)) { 
    if (subject %in% sc$SubjectID) {
      hset = ds[(ds$SubjectID == subject & ds$CL4 == "Healthy"), ]
      if ( sc[sc$SubjectID == subject, "IRIS"] == "IS") {
        h_IS <- rbind(h_IS, hset)
        }
      if ( sc[sc$SubjectID == subject, "IRIS"] == "IR") {
        h_IR <- rbind(h_IR, hset)
        }
      }
    }
  h_IS[, 1:m] <- apply(h_IS[, 1:m], 2, function(x) as.numeric(x))
  h_IS$SubjectID <- factor(h_IS$SubjectID)
  rmc.dset <- rmcorr(participant = SubjectID, measure1 = colnames(h_IS)[1], measure2 = colnames(h_IS)[21], dataset = h_IS)
  plot(rmc.dset, palette = IS_color, xlab="A1C levels", ylab="GLU levels")
  
  #For InterCor matrixes, first run the first part of the function
  #29 IS, 34 IR
  h_IS = as.data.frame(h_IS)
  colnames(h_IS)[50] <- "SubjectID"
  h_IS <- h_IS[h_IS$SubjectID != "69-049", ]  #for IR: h_IR <- h_IR[h_IR$SubjectID != "69-030", ]
  h_IS$SubjectID <- factor(h_IS$SubjectID)
  ggplot(data=h_IS, aes(x=as.numeric(as.character(A1C)), y=as.numeric(as.character(GLU)), color=SubjectID)) + geom_point(size = 3) +
    stat_smooth(method = "lm", col = "red") + scale_color_manual(values=IS_color(29), guide=FALSE) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text.x = element_text(size = 12)) + labs(x="A1C levels", y="GLU levels")
  ggsave("Desktop/Revised_CorComp_inter_IS.pdf", width = 6, height = 4)
  
  #Export significant associations
  HPC.tb = rbind(HPC.ck, HPC.clinic, HPC.metb, HPC.prot, HPC.st, HPC.st.ck, HPC.st.clinic, HPC.st.metb, HPC.st.prot)
  HPC.sig = HPC.tb[apply(cbind(HPC.tb[, 4]<0.05, HPC.tb[, 8]<0.05, HPC.tb[, 12]<0.05), 1, any), ]
  write.table(HPC.sig, file = "Desktop/Between_Subject_Associations.txt", sep = "\t", row.names = FALSE)
