


###===For revision: correlation using rmcorr, considering repeated sampling from subjects
##Healthy visits, IR and IS
#Loading dataset
load("Box Sync/PHI_Protected/Shared_WZ_HMP/Revision/Data/Revision_MultiOmes.RData") #Revised version shared on 04/18/18

#Build custom fuction, uses rmcorr to find 
#Provding a data frame with 1:m variable, SubjectID and CL4
#Return tables with correlation IS first and IR next
library(rmcorr)
  myHC <- function(ds, m){
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
    cor_ISIR = as.data.frame(matrix(, nrow = m*(m-1)/2, ncol = 0))
    #for IS
    h_IS$SubjectID <- factor(h_IS$SubjectID)
    #Transform
    if (trans == "none") {h_IS[, 1:m] <- apply(h_IS[, 1:m], 2, function(x) as.numeric(x))} 
    if (trans == "log2") {h_IS[, 1:m] <- apply(h_IS[, 1:m], 2, function(x) log2(as.numeric(x)))} 
    if (trans == "arcsin") {h_IS[, 1:m] <- apply(h_IS[, 1:m], 2, function(x) asin(sqrt(as.numeric(x))))}
    if (trans == "arcsin_log2") {
      h_IS[, 1:m] <- apply(h_IS[, 1:m], 2, function(x) asin(sqrt(as.numeric(x))))
      h_IS[, (m+1):n] <- apply(h_IS[, (m+1):n], 2, function(x) log2(as.numeric(x)))
      }
    cor_char = cor_r = cor_p = CI_left = CI_right = c()
    for (i in 1:(m-1)) { #for intra-omes
      for (j in (i+1):m) { #for intra-omes
    #for (i in 1:m) { #for inter-omes
      #for (j in (m+1):n) { #for inter-omes
        rmc.dset <- rmcorr(participant = SubjectID, measure1 = colnames(h_IS)[i], measure2 = colnames(h_IS)[j], dataset = h_IS) #allows uneven sampling
        cor_char <- c(cor_char, paste(colnames(h_IS)[i], colnames(h_IS)[j], sep = " vs "))
        cor_r <- c(cor_r, unlist(rmc.dset[1]))
        cor_p <- c(cor_p, unlist(rmc.dset[3]))
        CI_left <- c(CI_left, unlist(rmc.dset[4])[1])
        CI_right <- c(CI_right, unlist(rmc.dset[4])[2])
      }
    }
    cor_IS = data.frame(cor_char, cor_r, cor_p, CI_left, CI_right)
    cor_IS$adj.p = p.adjust(cor_p, method = "fdr")
    #for IR
    h_IR$SubjectID <- factor(h_IR$SubjectID)
    #Transform
    if (trans == "none") {h_IR[, 1:m] <- apply(h_IR[, 1:m], 2, function(x) as.numeric(x))} 
    if (trans == "log2") {h_IR[, 1:m] <- apply(h_IR[, 1:m], 2, function(x) log2(as.numeric(x)))} 
    if (trans == "arcsin") {h_IR[, 1:m] <- apply(h_IR[, 1:m], 2, function(x) asin(sqrt(as.numeric(x))))} 
    if (trans == "arcsin_log2") {
      h_IR[, 1:m] <- apply(h_IR[, 1:m], 2, function(x) asin(sqrt(as.numeric(x))))
      h_IR[, (m+1):n] <- apply(h_IR[, (m+1):n], 2, function(x) log2(as.numeric(x)))
      }
    cor_char = cor_r = cor_p = CI_left = CI_right = c()
    for (i in 1:(m-1)) { #for intra-omes
      for (j in (i+1):m) { #for intra-omes
    #for (i in 1:m) { #for inter-omes
      #for (j in (m+1):n) { #for inter-omes
        rmc.dset <- rmcorr(participant = SubjectID, measure1 = h_IR[, i], measure2 = h_IR[, j], dataset = h_IR) #allows uneven sampling
        cor_char <- c(cor_char, paste(colnames(h_IR)[i], colnames(h_IR)[j], sep = " vs "))
        cor_r <- c(cor_r, unlist(rmc.dset[1]))
        cor_p <- c(cor_p, unlist(rmc.dset[3]))
        CI_left <- c(CI_left, unlist(rmc.dset[4])[1])
        CI_right <- c(CI_right, unlist(rmc.dset[4])[2])
      }
    }
    cor_IR = data.frame(cor_char, cor_r, cor_p, CI_left, CI_right)
    cor_IR$adj.p = p.adjust(cor_p, method = "fdr")
    cor_ISIR <- cbind(cor_IS, cor_IR) #IS first
    cor_ISIR$Sig_num = rowSums(cbind(cor_ISIR[, 6] < 0.05, cor_ISIR[, 12] < 0.05))
    return(cor_ISIR)
  } #rmcorr doesn't work properly within a function
  
  #Intra-omes:
  ds = ck.df[, c(4:65, 70, 75)]  
  m = 62
  trans = "none"
  HC.ck = cor_ISIR
  ds = clinic.df[, c(2:28, 30:50, 52, 53, 58)] #No 29, 51
  m = 49
  trans = "none"
  HC.clinic = cor_ISIR
  ds = metbcr.df[, c(2:726, 731)]
  m = 724
  trans = "log2"
  HC.metb = cor_ISIR
  ds = swathprot.PCR.df[, c(2:303, 304, 309)]
  m = 302
  trans = "none"
  HC.prot = cor_ISIR
  ds = st.df[, c(2:97, 103, 102)]
  m = 96
  trans = "arcsin"
  HC.st = cor_ISIR
  save(HC.ck, HC.clinic, HC.metb, HC.prot, HC.st, file = "Box Sync/WenyuEx/Rcodes/Revised_IntraOmes_rmcorr.RData")
  
  
  #Inter-omes to gut microbes
  st_metb.df <- merge(st.df[, 1:97], metbcr.df, by.x = "HostSampleID", by.y = "SampleID")
  ds = st_metb.df[, c(2:821, 822, 827)] #2:97 st 98:821 metb
  m = 96
  n = 820
  trans = "arcsin_log2"
  HC.st.metb = cor_ISIR
  
  st_prot.df <- merge(st.df[, 1:97], swathprot.PCR.df, by.x = "HostSampleID", by.y = "SampleID")
  ds = st_prot.df[, c(2:399, 400, 405)] #2:97 st 98:399 proteins
  m = 96
  n = 398
  trans = "arcsin"
  HC.st.prot = cor_ISIR
  
  st_ck.df <- merge(st.df[, 1:97], ck.df[, c(1,4:75)], by.x = "HostSampleID", by.y = "SampleID")
  ds = st_ck.df[, c(2:159, 164, 169)] #2:97 st 98:159 cytokines
  m = 96
  n = 158
  trans = "arcsin"
  HC.st.ck = cor_ISIR
  
  clinic.clean = apply(clinic.df[, c(2:28, 30:50, 52)], 2, function(x) as.numeric(x)) #No col 29 and 51
  clinic.clean = cbind(clinic.clean, clinic.df[, c(1, 53, 58)])
  st_clinic.df <- merge(st.df[, 1:97], clinic.clean, by.x = "HostSampleID", by.y = "SampleID")
  ds = st_clinic.df[, c(2:146, 147, 148)] #2:97 st 98:146 clinic labs
  m = 96
  n = 145
  trans = "arcsin"
  HC.st.clinic = cor_ISIR
  
  st_rna.df <- merge(st.df[, 1:97], rnaseq.log.df, by.x = "HostSampleID", by.y = "SampleID")
  ds = st_rna.df[, c(2:10443, 10444, 10449)] #2:97 st 98:10443 rna_genes
  m = 96
  n = 10442
  trans = "arcsin"
  HC.st.rna = cor_ISIR
  
  save(HC.st.metb, HC.st.prot, HC.st.ck, HC.st.clinic, HC.st.rna, file = "Box Sync/WenyuEx/Rcodes/Revised_InterOmes_rmcorr.RData")

#Inter-omes to A1C, GLU, HSCRP
  setlab.df = clinic.df[, c(1, 2, 22, 26)]
  setlab.df$HSCRP <- as.numeric(setlab.df$HSCRP)
  lab_metb.df <- merge(setlab.df, metbcr.df, by = "SampleID")
  ds = lab_metb.df[, c(5:728, 2:4, 729, 734)] #5:728 metb 
  m = 724
  n = 727
  trans = "log2"
  HC.labs.metb = cor_ISIR
  
  lab_prot.df <- merge(setlab.df, swathprot.PCR.df, by = "SampleID")
  ds = lab_prot.df[, c(5:306, 2:4, 307, 312)] #5:306 proteins
  m = 302
  n = 305
  trans = "none"
  HC.labs.prot = cor_ISIR
  
  lab_ck.df <- merge(setlab.df, ck.df[, c(1,4:75)], by = "SampleID")
  ds = lab_ck.df[, c(5:66, 2:4, 71, 76)] #5:66 cytokines
  m = 62
  n = 65
  trans = "none"
  HC.labs.ck = cor_ISIR
  
  lab_rna.df <- merge(setlab.df, rnaseq.log.df, by = "SampleID")
  ds = lab_rna.df[, c(5:10350, 2:4, 10351, 10356)] #5:10350 proteins
  m = 10346
  n = 10349
  trans = "none"
  HC.labs.rna = cor_ISIR
  
  ds = merge(setlab.df, clinic.df[, c(1, 3:21, 23:25, 27:28, 30:50, 52, 53, 58)], by = "SampleID") #No 29, 51
  ds = ds[, c(5:50, 2:4,51,52)]
  m = 46
  n = 49
  trans = "none"
  HC.labs = cor_ISIR
  
  lab_st.df <- merge(setlab.df, st.df, by.x = "SampleID", by.y = "HostSampleID")
  ds = lab_st.df[, c(5:100, 2:4, 106, 105)] #5:100 gut micorbes
  m = 96
  n = 99
  trans = "arcsin"
  HC.labs.st = cor_ISIR
  
  save(HC.labs, HC.labs.ck, HC.labs.metb, HC.labs.prot, HC.labs.rna, HC.labs.st, file = "Box Sync/WenyuEx/Rcodes/Revised_rmcorr_labs.RData")
  
#Visualization    
  #Plotting examples in IR and IS
  #IL17F and IFNB: both significant, IL5 and VEGFD: IS sig, TGFB and IL15: IR sig
  char1 = "IL17F" #or: char1="IL5" char1="TGFB"
  char2 = "IFNB"  #or: char2="VEGFD" char2="IL15"
  rmc.IS = rmcorr(participant = SubjectID, measure1 = char1, measure2 = char2, dataset = h_IS) 
  plot(rmc.IS, xlab = paste(char1), ylab = paste(char2))
  rmc.IR = rmcorr(participant = SubjectID, measure1 = char1, measure2 = char2, dataset = h_IR) 
  plot(rmc.IR, xlab = paste(char1), ylab = paste(char2))
  
  #Barplot with error bars, showing the average edge per analyte (epa)
  #Preparing matrix
  epa = ave.epa = ratio = cor.tb = {}
  #for (df in c("HC.st.ck", "HC.st.clinic", "HC.st.metb", "HC.st.prot")){
  for (df in c("HC.ck", "HC.clinic", "HC.metb", "HC.prot", "HC.st", "HC.st.ck", "HC.st.clinic", "HC.st.metb", "HC.st.prot", "HC.st.rna")){
    ds = get(df)
    #For IS, 5% FDR
    ds.sig = as.character(ds[(ds[, 6] < 0.05), 1])
    ds_IS = ds[(ds[, 6] < 0.05), 1:6]
    print(paste(c(df, "IS", length(ds.sig)), collapse = " "))
    epaIS = data.frame(table(unlist(strsplit(ds.sig, split = " vs "))))
    epaIS <- data.frame(epaIS, ISIR=rep("IS", nrow(epaIS)), Type=rep(df, nrow(epaIS)))
    ave.IS = c(df, "IS", mean(epaIS$Freq), sd(epaIS$Freq))
    #For IR, 5% FDR
    ds.sig = as.character(ds[(ds[, 12] < 0.05), 1])
    ds_IR = ds[(ds[, 12] < 0.05), 7:12]
    print(paste(c(df, "IR", length(ds.sig)), collapse = " "))
    epaIR = data.frame(table(unlist(strsplit(ds.sig, split = " vs "))))
    epaIR <- data.frame(epaIR, ISIR=rep("IR", nrow(epaIR)), Type=rep(df, nrow(epaIR)))
    ave.IR = c(df, "IR", mean(epaIR$Freq), sd(epaIR$Freq))
    if (nrow(epaIS) == 0 |nrow(epaIR) == 0){
      ISIR = {}
    } else {
      ISIR = merge(epaIS[, 1:3], epaIR, by = "Var1", all = TRUE)
      ISIR$ratio = ISIR[, 2]/ISIR[, 4]
      cor_sig = merge(ds_IS, ds_IR, by = "cor_char", all = TRUE)
      cor_sig$Type = rep(df, nrow(cor_sig))
    }
    epa = rbind(epa, epaIS, epaIR)
    ave.epa = rbind(ave.epa, ave.IS, ave.IR)
    ratio = rbind(ratio, ISIR)
    cor.tb = rbind(cor.tb, cor_sig)
  }
  ave.epa <- as.data.frame(ave.epa)
  colnames(ave.epa) <- c("Type", "ISIR", "Freq_mean", "Freq_sd")
  ave.epa[, 3:4] <- apply(ave.epa[, 3:4], 2, function(x) as.numeric(as.character(x))) 
  
  ave.epa$ISIR <- factor(ave.epa$ISIR, levels = c("IS", "IR"))
  ggplot(data = ave.epa, aes(x = Type, y = Freq_mean, fill = ISIR)) + geom_bar(position = "dodge", stat = "identity") +
    geom_errorbar(aes(ymin=Freq_mean, ymax=Freq_mean+Freq_sd), width=.1, position=position_dodge(.9)) + 
    scale_fill_manual(breaks = c("IS", "IR"), values=c("#90AA3C", "#EF6125")) +
    #scale_x_discrete(breaks=c("HC.st.ck", "HC.st.clinic", "HC.st.metb", "HC.st.prot"),labels=c("GutMicrobes\nCytokines", "GutMicrobe\nClinicLabs", "GutMicrobe\nMetabolites", "GutMicrobes\nProteins")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text.x = element_text(size = 10, angle = 45, hjust = 1)) + 
    labs(x="Inter-omes Association", y="Average number of edges per analyte", fill='')
  ggsave("Desktop/Revised_Fig5_Interome_cor.pdf", height = 4, width = 6)
  
#Inspecting individual associations
  taxa = "genus_Veillonella"  # taxa = "genus_unclassified_Firmicutes"
  HC.subset = HC.st.metb[grepl(taxa, as.character(HC.st.metb[, 1])) & HC.st.metb[, 6] < 0.05, ]
  HC.subset <- rbind(HC.subset, HC.st.metb[grepl(taxa, as.character(HC.st.metb[, 1])) & HC.st.metb[, 12] < 0.05, ])
  HC.subset$Compounds_ID = gsub(".*vs ", "", as.character(HC.subset[, 1]))
  HC.subset <- merge(HC.subset, metb.curated[, c(2:5, 14, 15)], by = "Compounds_ID")
  
  #Preparing dataset for plotting
  #Ribbon plot
  pd = rbind(HC.subset[, c(15, 3:7)], HC.subset[, c(15, 9:13)])
  pd$ISIR = c(rep("IS", nrow(HC.subset)), rep("IR", nrow(HC.subset)))
  pd$ISIR <- factor(pd$ISIR, levels = c("IS", "IR"))
  pd$Metabolite <- factor(pd$Metabolite, levels = as.character(HC.subset[order(HC.subset[, 7], decreasing = TRUE), "Metabolite"]))
    
  ggplot(data=pd, aes(x=Metabolite, y=cor_r, fill=ISIR)) + geom_ribbon(aes(ymin=CI_left, ymax=CI_right, color=ISIR), size = 4, position=position_dodge(0.6)) + 
    geom_point(position = position_dodge(width=0.6), size = 4) + geom_hline(yintercept = 0, size = 0.3) +
    coord_flip() + scale_color_manual(breaks = c("IS", "IR"), values=c("#cddb9d", "#ffc2a8")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line.x = element_line(colour = "black", size =1),
          axis.line.y = element_blank(), axis.text = element_text(size = 12)) + labs(x="",y="Correlation Coefficient r", color="") + guides(fill=FALSE)
  ggsave("Desktop/Revised_Fig5_cor_Veillonella.pdf", height = 6, width = 10, useDingbats=FALSE)
  
  #p-value heatmap
  pd = HC.subset[, c(15, 7, 13)]
  pd[, 2:3] <- apply(pd[, 2:3], 2, function(x) -log10(x))
  colnames(pd)[2:3] <- c("IS", "IR")
  pd <- melt(pd)
  pd$Metabolite <- factor(pd$Metabolite, levels = as.character(HC.subset[order(HC.subset[, 7], decreasing = TRUE), "Metabolite"]))
  
  ggplot(pd, aes(variable, Metabolite)) + geom_tile(aes(fill = value), colour = "black") +
    scale_fill_gradient(low = "beige",high = "red") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.text.x = element_text(size = 12),
          axis.text.y = element_blank()) + labs(x="", y="", fill="-log10\n(p-value)")
  ggsave("Desktop/Revised_Fig5_cor_Veillonella_pvalue.pdf", height = 6, width = 2)
  
#IS vs IR specific visualization
  IR.spec = cor.tb[(is.na(cor.tb[, 2]) & !is.na(cor.tb[, 7])), ]  #IS unique 46483, IR unique 42690, all 138707
  IS.spec = cor.tb[(!is.na(cor.tb[, 2]) & is.na(cor.tb[, 7])), ] 
  IR.perc.type = table(IR.spec$Type)/c(nrow(HC.ck), nrow(HC.clinic), nrow(HC.metb), nrow(HC.prot), nrow(HC.st), nrow(HC.st.clinic), nrow(HC.st.metb), nrow(HC.st.prot))
  IS.perc.type = table(IS.spec$Type)/c(nrow(HC.ck), nrow(HC.clinic), nrow(HC.metb), nrow(HC.prot), nrow(HC.st), nrow(HC.st.clinic), nrow(HC.st.metb), nrow(HC.st.prot))
  perc.tb = data.frame(IR=IR.perc.type/sum(IR.perc.type), IS=IS.perc.type/sum(IS.perc.type))
  
  ggplot(data=perc.tb[, 1:2], aes(x="", y=IR.Freq, fill=IR.Var1)) + geom_bar(stat = "identity") + coord_polar("y") + labs(fill="Associations") +
    theme_minimal() + theme(axis.title.x = element_blank(),axis.title.y = element_blank(),panel.border = element_blank(),axis.text.x=element_blank(), 
                            panel.grid=element_blank(),axis.ticks = element_blank(),plot.title=element_text(size=14, face="bold")) +
    scale_fill_hue(labels=c("Cytokines", "ClinicLabs", "Metabolites", "Proteins", "GutMicrobes", "GutMicrobes_ClinicLabs", "GutMicrobes_Metabolites", "GutMicrobes_Proteins"))
  ggsave("Desktop/Revised_Fig5_CorSpec_IR.pdf", height = 6, width = 6)
  
#IS and IR specific gut genus associations
  #IR
  IR.st = IR.spec[IR.spec$Type == "HC.st", c(1, 7:11)]
  IR.genus = IR.st[grepl("genus_.*genus_.*", as.character(IR.st$cor_char)), ] #At genus level IR 222/2, IS 332/2
  IR.genus$node1 = gsub(" vs genus_.*", "", as.character(IR.genus[, 1]))
  IR.genus$node2 = gsub("genus_.*vs ", "", as.character(IR.genus[, 1]))
  IR.list = as.data.frame(table(unlist(strsplit(as.character(IR.genus[, 1]), split = " vs "))))
  #IS
  IS.st = IS.spec[IS.spec$Type == "HC.st", 1:6]
  IS.genus = IS.st[grepl("genus_.*genus_.*", as.character(IS.st$cor_char)), ]
  IS.genus$node1 = gsub(" vs genus_.*", "", as.character(IS.genus[, 1]))
  IS.genus$node2 = gsub("genus_.*vs ", "", as.character(IS.genus[, 1]))
  IS.list = as.data.frame(table(unlist(strsplit(as.character(IS.genus[, 1]), split = " vs "))))
  
  #Barplot for each genus, each pair was counted twice due a bug?
  list.genus = merge(IS.list, IR.list, by = "Var1", all = TRUE)
  colnames(list.genus)[2:3] <-c("IS", "IR")
  list.genus[, 2:3] <- list.genus[, 2:3]/2
  pd = melt(list.genus)
  ggplot(data=pd, aes(x=as.character(Var1), y=value, fill=variable)) + geom_bar(position = "dodge", stat = "identity") +
    scale_fill_manual(breaks = c("IS", "IR"), values=c("#90AA3C", "#EF6125")) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text.x = element_blank(), plot.background = element_blank(), axis.title.y = element_text(size=5)) + scale_y_continuous(position = "right") +
    labs(x="", y="Number of Associations", fill="")
  ggsave("Desktop/Revised_CorSpec_st_nodes.pdf", width = 6, height = 1) 
  
  #Dotplot for coefficient of each genus, each pair was counted twice due a bug?
  r.genus = as.data.frame(rbind(as.matrix(IS.genus[1:166, c(7, 2)]), as.matrix(IS.genus[1:166, c(8, 2)]), 
                  as.matrix(IR.genus[1:111, c(7, 2)]), as.matrix(IR.genus[1:111, c(8, 2)])))
  r.genus$ISIR = c(rep("IS", nrow(IS.genus)), rep("IR", nrow(IR.genus)))
  r.genus$ISIR <- factor(r.genus$ISIR, levels=c("IS", "IR"))
  ggplot(data=r.genus, aes(x=as.character(node1), y=as.numeric(as.character(cor_r.x)), fill=ISIR)) + 
    geom_point(size=1.5, alpha = 0.8, colour="black", shape=21, stroke = 0.2, position=position_jitterdodge(dodge.width=0.6)) +
    scale_fill_manual(breaks = c("IS", "IR"), values=c("#90AA3C", "#EF6125")) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_line(colour = "grey60"),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text.x = element_text(angle=45, vjust=1, hjust=1, size=6)) +
  labs(x="", y="Correlation Coefficient", color="")
  ggsave("Desktop/Revised_CorSpec_st_cor.pdf", width = 6, height = 4)
  
  #Tile plots to visualize all interactions using geom_tile
  HC.st.genus = HC.st[grepl("genus_.*genus_.*", as.character(HC.st$cor_char)), ]
  HC.st.genus$node1 = gsub(" vs genus_.*", "", as.character(HC.st.genus[, 1]))
  HC.st.genus$node2 = gsub("genus_.*vs ", "", as.character(HC.st.genus[, 1]))
  
  HC.st.IS = HC.st.genus[, c(2,6,14,15)]
  ggplot(data=HC.st.IS, aes(x=node1,y=node2)) + geom_tile(fill = "white", color = "black") +
    geom_tile(data=HC.st.IS[HC.st.IS$adj.p < 0.05, ], aes(fill=cor_r)) + theme_minimal() +
    scale_fill_gradient2(low = "blue4", mid = "white",high = "red") +
    theme(axis.text.x = element_text(angle = 30, size = 8, hjust = 1)) + labs(x="", y="", fill="Correlation\nCoefficient r")
  ggsave("Desktop/Revsied_st_cor_IStile.pdf", width = 12, height = 8)
  HC.st.IR = HC.st.genus[, c(8,12,14,15)]
  
#IS specific gut-metb
  #IS
  IS.st.metb = IS.spec[IS.spec$Type == "HC.st.metb", 1:6]
  IS.st.metb$node1 = gsub(" vs .*", "", as.character(IS.st.metb[, 1]))
  IS.st.metb$node2 = gsub(".*vs ", "", as.character(IS.st.metb[, 1]))
  IS.st.metb <- merge(IS.st.metb, metb.curated[, c(2:5, 14, 15)], by.x = "node2", by.y = "Compounds_ID")
  #IR
  #IS
  IR.st.metb = IR.spec[IR.spec$Type == "HC.st.metb",  c(1, 7:11)]
  IR.st.metb$node1 = gsub(" vs .*", "", as.character(IR.st.metb[, 1]))
  IR.st.metb$node2 = gsub(".*vs ", "", as.character(IR.st.metb[, 1]))
  IR.st.metb <- merge(IR.st.metb, metb.curated[, c(2:5, 14, 15)], by.x = "node2", by.y = "Compounds_ID")
  
  write.table(rbind(as.matrix(IS.st.metb), as.matrix(IR.st.metb)), file = "Desktop/Revised_Fig5d_st_metb_cor.txt", sep = "\t", row.names = FALSE)
  
  #Export significant associations
  HC.tb = rbind(HC.ck, HC.clinic, HC.metb, HC.prot, HC.st, HC.st.ck, HC.st.clinic, HC.st.metb, HC.st.prot, HC.st.rna)
  HC.sig = HC.tb[apply(cbind(HC.tb[, 6]<0.05, HC.tb[, 12]<0.05), 1, any), ]
  write.table(HC.sig, file = "Desktop/Within_Subject_Associations.txt", sep = "\t", row.names = FALSE)
  