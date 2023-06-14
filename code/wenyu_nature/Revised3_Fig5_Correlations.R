###Re-generate Fig5a and 5b using python sparCC data
#microbiome to microbiome: sparCC only at the same rank
#Excuted 04/07/19 with added modifications
  #IS
  epaIS.all = tb.IS = {}
  n.epaIS = 0
  for (L in c("phylum", "class", "order", "family", "genus")){ 
    fname = paste(c("Box Sync/PHI_Protected/Shared_WZ_HMP/Revision3/Ahmed/sparcc_pvalue_ipop/output/", L, "_IS_cor_sparcc.txt"), collapse = "")
    fname.p = paste(c("Box Sync/PHI_Protected/Shared_WZ_HMP/Revision3/Ahmed/sparcc_pvalue_ipop/pvalue/", L, "_IS_one_sided.txt"), collapse = "")
    #For IS, abs(rho) > 0.2 cutoff
    ds.IS = read.table(fname)
    ds.IS.p = read.table(fname.p)
    pr.name = c()
    for (i in 1:(nrow(ds.IS)-1)) {pr.name <- c(pr.name, paste(colnames(ds.IS)[i], rownames(ds.IS)[(i+1):nrow(ds.IS)], sep = " vs "))}
    ds.IS.transf = data.frame(cor_char=pr.name, rho=ds.IS[lower.tri(ds.IS)], p=ds.IS.p[lower.tri(ds.IS.p)])
    ds.IS.sig = ds.IS.transf[abs(ds.IS.transf$rho) > 0.2 & ds.IS.transf$p < 0.05, ]
    epaIS = data.frame(table(unlist(strsplit(as.character(ds.IS.sig[, 1]), " vs "))))
    epaIS <- data.frame(epaIS, ISIR=rep("IS", nrow(epaIS)), Type=rep(L, nrow(epaIS)))
    epaIS.all <- rbind(epaIS.all, epaIS)
    n.epaIS <- n.epaIS+nrow(ds.IS.sig)
    tb.IS = rbind(tb.IS, ds.IS.sig)
  }
  ave.IS = c("sparCC.st", "IS", mean(epaIS$Freq), sd(epaIS$Freq)) #Only at genus level
  print(c("IS", n.epaIS))
  #IR
  epaIR.all = tb.IR = {}
  n.epaIR = 0
  for (L in c("phylum", "class", "order", "family", "genus")){ 
    fname = paste(c("Box Sync/PHI_Protected/Shared_WZ_HMP/Revision3/Ahmed/sparcc_pvalue_ipop/output/", L, "_IR_cor_sparcc.txt"), collapse = "")
    fname.p = paste(c("Box Sync/PHI_Protected/Shared_WZ_HMP/Revision3/Ahmed/sparcc_pvalue_ipop/pvalue/", L, "_IR_one_sided.txt"), collapse = "")
    #For IS, abs(rho) > 0.2 cutoff
    ds.IR = read.table(fname)
    ds.IR.p = read.table(fname.p)
    pr.name = c()
    for (i in 1:(nrow(ds.IR)-1)) {pr.name <- c(pr.name, paste(colnames(ds.IR)[i], rownames(ds.IR)[(i+1):nrow(ds.IR)], sep = " vs "))}
    ds.IR.transf = data.frame(cor_char=pr.name, rho=ds.IR[lower.tri(ds.IR)], p=ds.IR.p[lower.tri(ds.IR.p)])
    ds.IR.sig = ds.IR.transf[abs(ds.IR.transf$rho) > 0.2 & ds.IR.transf$p < 0.05, ]
    epaIR = data.frame(table(unlist(strsplit(as.character(ds.IR.sig[, 1]), " vs "))))
    epaIR <- data.frame(epaIR, ISIR=rep("IR", nrow(epaIR)), Type=rep(L, nrow(epaIR)))
    epaIR.all <- rbind(epaIR.all, epaIR)
    n.epaIR <- n.epaIR+nrow(ds.IR.sig)
    tb.IR = rbind(tb.IR, ds.IR.sig)
  }
  ave.IR = c("sparCC.st", "IR", mean(epaIR$Freq), sd(epaIR$Freq))
  print(c("IR", n.epaIR))
  
  write.table(merge(tb.IS, tb.IR, by = "cor_char", all = TRUE), file = "Desktop/Revised3_InterMicrobiome.txt", sep = "\t", row.names = FALSE)
  
  #Barplot for each genus, each pair was counted twice due a bug?
  list.genus = merge(epaIS.all[44:88, ], epaIR.all[51:95, ], by = "Var1", all = TRUE)
  colnames(list.genus)[c(2,5)] <-c("IS", "IR")
  
  pd = melt(list.genus[, c(1,2,5)])
  ggplot(data=pd, aes(x=as.character(Var1), y=value, fill=variable)) + geom_bar(position = "dodge", stat = "identity") +
    scale_fill_manual(breaks = c("IS", "IR"), values=c("#90AA3C", "#EF6125")) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text.x = element_blank(), plot.background = element_blank(), axis.title.y = element_text(size=5)) + scale_y_continuous(position = "right") +
    labs(x="", y="Number of Associations", fill="")
  ggsave("Desktop/Revised3_CorSpec_st_nodes.pdf", width = 6, height = 1) 
  
  #Dotplot for coefficient of each genus
  #IR
  IR.genus = tb.IR[84:278, ] #At genus level IR 195
  IR.genus$node1 = gsub(" vs genus_.*", "", as.character(IR.genus[, 1]))
  IR.genus$node2 = gsub("genus_.*vs ", "", as.character(IR.genus[, 1]))
  IR.list = as.data.frame(table(unlist(strsplit(as.character(IR.genus[, 1]), split = " vs "))))
  #IS
  IS.genus = tb.IS[74:317, ] #At genus level IS 244
  IS.genus$node1 = gsub(" vs genus_.*", "", as.character(IS.genus[, 1]))
  IS.genus$node2 = gsub("genus_.*vs ", "", as.character(IS.genus[, 1]))
  IS.list = as.data.frame(table(unlist(strsplit(as.character(IS.genus[, 1]), split = " vs "))))
  
  
  r.genus = as.data.frame(rbind(as.matrix(IS.genus[, c(4, 2)]), as.matrix(IS.genus[, c(5, 2)]), 
                                as.matrix(IR.genus[, c(4, 2)]), as.matrix(IR.genus[, c(5, 2)])))
  r.genus$ISIR = c(rep("IS", 2*nrow(IS.genus)), rep("IR", 2*nrow(IR.genus)))
  r.genus$ISIR <- factor(r.genus$ISIR, levels=c("IS", "IR"))
  ggplot(data=r.genus, aes(x=as.character(node1), y=as.numeric(as.character(rho)), fill=ISIR)) + 
    geom_point(size=1.5, alpha = 0.8, colour="black", shape=21, stroke = 0.2, position=position_jitterdodge(dodge.width=0.6)) +
    scale_fill_manual(breaks = c("IS", "IR"), values=c("#90AA3C", "#EF6125")) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_line(colour = "grey60"),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text.x = element_text(angle=45, vjust=1, hjust=1, size=6)) +
    labs(x="", y="Correlation Coefficient", color="")
  ggsave("Desktop/Revised_CorSpec_st_cor.pdf", width = 6, height = 4)
  
  #Tile plots to visualize interactions using geom_tile for abs(rho) > 0.2
  ggplot(data=IS.genus, aes(x=node1,y=node2)) + geom_tile(fill = "white", color = "black") +
    geom_tile(data=, aes(fill=rho)) + theme_minimal() +
    scale_fill_gradient2(low = "blue4", mid = "white",high = "red") +
    theme(axis.text.x = element_text(angle = 30, size = 8, hjust = 1)) + labs(x="", y="", fill="Correlation\nCoefficient r")
  ggsave("Desktop/Revsied_st_cor_IStile.pdf", width = 12, height = 8)
  #Then use data=IR.genus
  
###Sent on 02/10 for python sparCC results
##Compare Python sparCC and R SpiecEasi package
load("Box Sync/PHI_Protected/Shared_WZ_HMP/Revision3/Ahmed/IntraOme_microbiome_sparcc_IR_IS_Healthy.RData")
  rse.IS = genus_sparcc$IS_sparcc
  rse.IR = genus_sparcc$IR_sparcc
  
  py.IS = read.table("Box Sync/PHI_Protected/Shared_WZ_HMP/Revision3/Ahmed/sparcc_pvalue_ipop/output/genus_IS_cor_sparcc.txt")
  py.IR = read.table("Box Sync/PHI_Protected/Shared_WZ_HMP/Revision3/Ahmed/sparcc_pvalue_ipop/output/genus_IR_cor_sparcc.txt")
  
  #Build data.frame for all
  pr.name = c()
  for (i in 1:(nrow(rse.IS)-1)) {pr.name <- c(pr.name, paste(colnames(rse.IS)[i], rownames(rse.IS)[(i+1):nrow(rse.IS)]))}
  
  pd = data.frame(cor_char = pr.name, rse_IS = rse.IS[lower.tri(rse.IS)], rse_IR = rse.IR[lower.tri(rse.IR)], 
                                      py_IS = py.IS[lower.tri(py.IS)], py_IR = py.IR[lower.tri(py.IR)])
   
  ggplot(data=pd, aes(x=rse_IS, y=py_IS)) + geom_point(color = "gray50", alpha = 0.8) + labs(x="SpiecEasi rho", y = "Python sparCC rho") +
    theme(plot.background = element_blank(),panel.background = element_blank(),
          axis.text.x = element_text(size = 12, hjust = 1),axis.title = element_text(size = 16),
          axis.line = element_line(color="black", size = 1), panel.border = element_blank())
  ggsave("Desktop/Python_SpiecEasi_sparcc_IS.pdf", width = 6, height = 4)
  

###Sent on 02/09 for microbiome to cytokines
load("Box Sync/PHI_Protected/Shared_WZ_HMP/Revision3/st_ck.df.RData")
load("Box Sync/WenyuEx/Research_HMP/OmicsData/Revision_MultiOmes.RData")

  #Extracting features with significant correlations
  tb = st_ck.df[, c(1, 41, 58, 64, 74, 105, 108, 119, 125, 164, 169)]
  
  tb.IR = tb.IS = {}
  for (subject in unique(tb$SubjectID)) {
    sset = tb[tb$SubjectID == subject & tb$CL4 == "Healthy", ]
    if (nrow(sc[sc$SubjectID == subject, ]) != 0) {
      if (sc[sc$SubjectID == subject, "IRIS"] == "IS") {
        tb.IS = rbind(tb.IS, sset)
      }
      if (sc[sc$SubjectID == subject, "IRIS"] == "IR") {
        tb.IR = rbind(tb.IR, sset)
      }
    }
  }
  
  library(rmcorr)
  tb.IS$SubjectID <- as.factor(tb.IS$SubjectID)
  tb.IS[, 2:5] <- apply(tb.IS[, 2:5], 2, function(x) as.numeric(as.character(x)))
  tb.IS[, 6:9] <- apply(tb.IS[, 6:9], 2, function(x) log10(x))
  
  rmc.dset <- rmcorr(participant = SubjectID, measure1 = IL1B, measure2 = genus_Barnesiella, dataset = tb.IS)
  plot(rmc.dset,  cex = 0.5)
  
  rmc.dset <- rmcorr(participant = SubjectID, measure1 = TNFA, measure2 = genus_Faecalibacterium, dataset = tb.IS)
  plot(rmc.dset,  cex = 0.5)
  #Same for tb.IR
  
  #Inspecting individual associations
  taxa = "genus_Veillonella"  # taxa = "genus_unclassified_Firmicutes"
  HC.subset = HC.st.metb[grepl(taxa, as.character(HC.st.metb[, 1])) & HC.st.metb[, 6] < 0.05, ]
  HC.subset <- rbind(HC.subset, HC.st.metb[grepl(taxa, as.character(HC.st.metb[, 1])) & HC.st.metb[, 12] < 0.05, ])
  HC.subset$Compounds_ID = gsub(".*vs ", "", as.character(HC.subset[, 1]))
  HC.subset <- merge(HC.subset, metb.curated[, c(2:5, 14, 15)], by = "Compounds_ID")
  
  #Preparing dataset for plotting
  clr.st.metb = read.csv("Box Sync/PHI_Protected/Shared_WZ_HMP/Revision3/HC.st.metb.csv", stringsAsFactors = FALSE)[, -1]
  #Ribbon plot
  clr.st.metb.sig = clr.st.metb[clr.st.metb$Sig_num > 0, ]
  pd = as.data.frame(rbind(as.matrix(clr.st.metb.sig[5:8, 2:6]), as.matrix(clr.st.metb.sig[5:8, 8:12])))  #Only for genus Butyricimonas
  pd$ISIR = c(rep("IS", 4), rep("IR", 4))
  pd$ISIR <- factor(pd$ISIR, levels = c("IS", "IR"))
  pd$Metabolite <- rep(c("LysoPC(17:0)", "LysoPE(22:0)", "C19H32O2", "Androsterone sulfate"), 2)
  pd$Metabolite <- factor(pd$Metabolite)
  
  ggplot(data=pd, aes(x=Metabolite, y=cor_r, fill=ISIR)) + geom_ribbon(aes(ymin=CI_left, ymax=CI_right, color=ISIR), size = 6, position=position_dodge(0.6)) + 
    geom_point(position = position_dodge(width=0.6), size = 6) + geom_hline(yintercept = 0, size = 0.3) +
    coord_flip() + scale_color_manual(breaks = c("IS", "IR"), values=c("#cddb9d", "#ffc2a8")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          plot.background = element_rect(fill = "transparent"), panel.background = element_blank(),
          axis.line.x = element_line(colour = "black", size =1),
          axis.line.y = element_blank(), axis.text = element_text(size = 14), 
          axis.title = element_text(size = 16)) + labs(x="",y="Correlation Coefficient", color="") + guides(fill=FALSE)
  ggsave("Desktop/Revised3_Fig5_cor_Butyricimonas.pdf", height = 4, width = 6, useDingbats=FALSE)
  
  #p-value heatmap
  pd = clr.st.metb.sig[5:8, c(6, 12)]
  pd <- as.data.frame(apply(pd, 2, function(x) -log10(x)))
  colnames(pd) <- c("IS", "IR")
  pd$Metabolite <- c("LysoPC(17:0)", "LysoPE(22:0)", "C19H32O2", "Androsterone sulfate")
  pd$Metabolite <- factor(pd$Metabolite)
  pd <- melt(pd)
  
  ggplot(pd, aes(variable, Metabolite)) + geom_tile(aes(fill = value), colour = "black") +
    scale_fill_gradient(low = "beige",high = "red") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.text.x = element_text(size = 16),
          axis.text.y = element_blank()) + labs(x="", y="", fill="-log10\n(p-value)")
  ggsave("Desktop/Revised3_Fig5_cor_Butyricimonas_pvalue.pdf", height = 4, width = 2)
  
###Sent on 02/08 for microbiome to host omes
clr.st.ck = read.csv("Box Sync/PHI_Protected/Shared_WZ_HMP/Revision3/HC.st.ck.csv", stringsAsFactors = FALSE)[, -1]
clr.st.clinic = read.csv("Box Sync/PHI_Protected/Shared_WZ_HMP/Revision3/HC.st.clinic.csv", stringsAsFactors = FALSE)[, -1]
clr.st.metb = read.csv("Box Sync/PHI_Protected/Shared_WZ_HMP/Revision3/HC.st.metb.csv", stringsAsFactors = FALSE)[, -1]
clr.st.prot = read.csv("Box Sync/PHI_Protected/Shared_WZ_HMP/Revision3/HC.st.prot.csv", stringsAsFactors = FALSE)[, -1]
load("Box Sync/PHI_Protected/Shared_WZ_HMP/Revision3/HC.st.rna_combined.RData")

write.table(rbind(clr.st.ck[clr.st.ck$Sig_num > 0, ], clr.st.clinic[clr.st.clinic$Sig_num > 0, ], 
                  clr.st.metb[clr.st.metb$Sig_num > 0, ], clr.st.prot[clr.st.prot$Sig_num > 0, ], 
                  HC.st.rna_combined[HC.st.rna_combined$Sig_num > 0, ]), file = "Desktop/Revised3_within_Intraomes.txt", sep = "\t", row.names = FALSE)

  #Visualization only for intra-omes, 
 
  ###Executed on 04/07/19 with new modifications
  #Preparing matrix
  load("Box Sync/WenyuEx/Rcodes/Revised_IntraOmes_rmcorr.RData")
  #Barplot with error bars, showing the average edge per analyte (epa)
  ave.epa = epa = {}
  for (df in c("HC.ck", "HC.clinic", "HC.metb", "HC.prot")){ 
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
    #rbind together
    ave.epa = rbind(ave.epa, ave.IS, ave.IR)
    epa = rbind(epa, epaIS, epaIR)
  }
  ave.epa <- as.data.frame(ave.epa)
  colnames(ave.epa) <- c("Type", "ISIR", "Freq_mean", "Freq_sd")
  ave.epa[, 3:4] <- apply(ave.epa[, 3:4], 2, function(x) as.numeric(as.character(x))) 
  ave.epa[, 1] <- as.character(ave.epa[, 1])

  #microbiome to microbiome: sparCC only at the same rank
  #IS
  epaIS.all = tb.IS = {}
  n.epaIS = 0
  for (L in c("phylum_sparcc", "class_sparcc", "order_sparcc", "family_sparcc", "genus_sparcc")){ 
    ds = get(L)
    #For IS, abs(rho) > 0.2 cutoff
    ds.IS = ds$IS_sparcc
    pr.name = c()
    for (i in 1:(nrow(ds.IS)-1)) {pr.name <- c(pr.name, paste(colnames(ds.IS)[i], rownames(ds.IS)[(i+1):nrow(ds.IS)]))}
    ds.IS.transf = data.frame(cor_char=pr.name, rho=ds.IS[lower.tri(ds.IS)])
    ds.IS.sig = ds.IS.transf[abs(ds.IS.transf$rho) > 0.2, ]
    epaIS = data.frame(table(unlist(strsplit(as.character(ds.IS.sig[, 1]), " "))))
    epaIS <- data.frame(epaIS, ISIR=rep("IS", nrow(epaIS)), Type=rep(L, nrow(epaIS)))
    epaIS.all <- rbind(epaIS.all, epaIS)
    n.epaIS <- n.epaIS+nrow(ds.IS.sig)
    tb.IS = rbind(tb.IS, ds.IS.sig)
  }
  ave.IS = c("sparCC.st", "IS", mean(epaIS$Freq), sd(epaIS$Freq))
  print(c("IS", n.epaIS))
  #IR
  epaIR.all = tb.IR = {}
  n.epaIR = 0
  for (L in c("phylum_sparcc", "class_sparcc", "order_sparcc", "family_sparcc", "genus_sparcc")){ 
    ds = get(L)
    #For IR, abs(rho) > 0.2 cutoff
    ds.IR = ds$IR_sparcc
    pr.name = c()
    for (i in 1:(nrow(ds.IR)-1)) {pr.name <- c(pr.name, paste(colnames(ds.IR)[i], rownames(ds.IR)[(i+1):nrow(ds.IR)]))}
    ds.IR.transf = data.frame(cor_char=pr.name, rho=ds.IR[lower.tri(ds.IR)])
    ds.IR.sig = ds.IR.transf[abs(ds.IR.transf$rho) > 0.2, ]
    epaIR = data.frame(table(unlist(strsplit(as.character(ds.IR.sig[, 1]), " "))))
    epaIR <- data.frame(epaIR, ISIR=rep("IR", nrow(epaIR)), Type=rep(L, nrow(epaIR)))
    epaIR.all <- rbind(epaIR.all, epaIR)
    n.epaIR <- n.epaIR+nrow(ds.IR.sig)
    tb.IR = rbind(tb.IR, ds.IR.sig)
  }
  ave.IR = c("sparCC.st", "IR", mean(epaIR$Freq), sd(epaIR$Freq))
  print(c("IR", n.epaIR))
  
  write.table(merge(tb.IS, tb.IR, by = "cor_char", all = TRUE), file = "Desktop/Revised3_InterMicrobiome.txt", sep = "\t", row.names = FALSE)
  
  ave.epa <- rbind(ave.epa, ave.IS, ave.IR)
  
  ave.epa[, 3:4] <- apply(ave.epa[, 3:4], 2, function(x) as.numeric(as.character(x))) 
  ave.epa$ISIR <- factor(ave.epa$ISIR, levels = c("IS", "IR"))
  ggplot(data = ave.epa, aes(x = Type, y = Freq_mean, fill = ISIR)) + geom_bar(position = "dodge", stat = "identity") +
    geom_errorbar(aes(ymin=Freq_mean, ymax=Freq_mean+Freq_sd), width=.1, position=position_dodge(.9)) + 
    scale_fill_manual(breaks = c("IS", "IR"), values=c("#90AA3C", "#EF6125")) +
    scale_x_discrete(breaks=c("HC.ck", "HC.clinic", "HC.metb", "HC.prot", "sparCC.st"),labels=c("Cytokines", "ClinicLabs", "Metabolites", "Proteins", "Gut_Microbiome")) +
    scale_y_log10() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text.x = element_text(size = 12, face = "bold", angle = 45, hjust = 1), 
          plot.title = element_text(hjust = 0.5, size = 14)) + 
    labs(x="", y="Average number of edges per analyte", fill='', title="Intra-omes Significant Associations")
  ggsave("Desktop/Revised3_Fig5_Intraome_cor.pdf", height = 4, width = 6)

###Sent on 02/08 for microbiome-microbiome sparCC
load("Box Sync/PHI_Protected/Shared_WZ_HMP/Revision3/Revised_IntraOmes_rmcorr_clr_abundant.RData")
load("Box Sync/PHI_Protected/Shared_WZ_HMP/Revision3/IntraOme_microbiome_sparcc_IR_IS_Healthy.RData")

  HC.st.sub = HC.st[grepl("genus_.*genus_.*", as.character(HC.st$cor_char)), ]  

  gs_IS = genus_sparcc$IS_sparcc
  gs_IR = genus_sparcc$IR_sparcc
 
  gs.IS = data.frame(cor_char=HC.st.sub$cor_char, rho=gs_IS[lower.tri(gs_IS)])
  gs.IR = data.frame(cor_char=HC.st.sub$cor_char, rho=gs_IR[lower.tri(gs_IR)])
  
  #IS and IR specific gut genus associations
  #IR
  IR.genus = gs.IR[abs(gs.IR$rho) > 0.2, ] #At genus level IR 203
  IR.genus$node1 = gsub(" vs genus_.*", "", as.character(IR.genus[, 1]))
  IR.genus$node2 = gsub("genus_.*vs ", "", as.character(IR.genus[, 1]))
  IR.list = as.data.frame(table(unlist(strsplit(as.character(IR.genus[, 1]), split = " vs "))))
  #IS
  IS.genus = gs.IS[abs(gs.IS$rho) > 0.2, ] #At genus level IS 244
  IS.genus$node1 = gsub(" vs genus_.*", "", as.character(IS.genus[, 1]))
  IS.genus$node2 = gsub("genus_.*vs ", "", as.character(IS.genus[, 1]))
  IS.list = as.data.frame(table(unlist(strsplit(as.character(IS.genus[, 1]), split = " vs "))))
  
  #Barplot for each genus, each pair was counted twice due a bug?
  list.genus = merge(IS.list, IR.list, by = "Var1", all = TRUE)
  colnames(list.genus)[2:3] <-c("IS", "IR")
  ##list.genus[, 2:3] <- list.genus[, 2:3]/2  no need, a bug in old scripts
  pd = melt(list.genus)
  ggplot(data=pd, aes(x=as.character(Var1), y=value, fill=variable)) + geom_bar(position = "dodge", stat = "identity") +
    scale_fill_manual(breaks = c("IS", "IR"), values=c("#90AA3C", "#EF6125")) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text.x = element_blank(), plot.background = element_blank(), axis.title.y = element_text(size=5)) + scale_y_continuous(position = "right") +
    labs(x="", y="Number of Associations", fill="")
  ggsave("Desktop/Revised3_CorSpec_st_nodes.pdf", width = 6, height = 1) 
  
  #Dotplot for coefficient of each genus, each pair was counted twice due a bug?
  r.genus = as.data.frame(rbind(as.matrix(IS.genus[, c(3, 2)]), as.matrix(IS.genus[, c(4, 2)]), 
                                as.matrix(IR.genus[, c(3, 2)]), as.matrix(IR.genus[, c(4, 2)])))
  r.genus$ISIR = c(rep("IS", 2*nrow(IS.genus)), rep("IR", 2*nrow(IR.genus)))
  r.genus$ISIR <- factor(r.genus$ISIR, levels=c("IS", "IR"))
  ggplot(data=r.genus, aes(x=as.character(node1), y=as.numeric(as.character(rho)), fill=ISIR)) + 
    geom_point(size=1.5, alpha = 0.8, colour="black", shape=21, stroke = 0.2, position=position_jitterdodge(dodge.width=0.6)) +
    scale_fill_manual(breaks = c("IS", "IR"), values=c("#90AA3C", "#EF6125")) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_line(colour = "grey60"),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text.x = element_text(angle=45, vjust=1, hjust=1, size=6)) +
    labs(x="", y="Correlation Coefficient", color="")
  ggsave("Desktop/Revised_CorSpec_st_cor.pdf", width = 6, height = 4)
  
  #Tile plots to visualize interactions using geom_tile for abs(rho) > 0.2
  ggplot(data=IS.genus, aes(x=node1,y=node2)) + geom_tile(fill = "white", color = "black") +
    geom_tile(data=, aes(fill=rho)) + theme_minimal() +
    scale_fill_gradient2(low = "blue4", mid = "white",high = "red") +
    theme(axis.text.x = element_text(angle = 30, size = 8, hjust = 1)) + labs(x="", y="", fill="Correlation\nCoefficient r")
  ggsave("Desktop/Revsied_st_cor_IStile.pdf", width = 12, height = 8)
  #Then use data=IR.genus

  #Comparing sparCC and clr_abundant, obtain IR.genus and IS.genus from sparcc matrixes above
  HC.st_clr.sub = HC.st_clr_abundant[grepl("genus_.*genus_.*", as.character(HC.st_clr_abundant$cor_char)), ]
  
  #plots of sparCC rho and clr rho for all
  pd = data.frame(CLR=HC.st_clr.sub[, 2], sparCC=gs.IS$rho) #IS
  ggplot(data=pd, aes(x=CLR, y=sparCC)) + geom_point(color = "gray50", alpha = 0.8) + labs(x="Centered Log Ratio") +
    theme(plot.background = element_blank(),panel.background = element_blank(),
          axis.text.x = element_text(size = 12, hjust = 1),axis.title = element_text(size = 16),
          axis.line = element_line(color="black", size = 1), panel.border = element_blank())
  ggsave("Desktop/Revsied_sparcc_clr_IR.pdf", width = 6, height = 4)
  pd = data.frame(CLR=HC.st_clr.sub[, 8], sparCC=gs.IR$rho) #IR
  
  #After filtering based on p or rho
  clr.IS = HC.st_clr.sub[HC.st_clr.sub[, 6] < 0.05, ]  #130 obs.
  clr.IR = HC.st_clr.sub[HC.st_clr.sub[, 12] < 0.05, ]  #97 obs.
  
  sum(clr.IS$cor_char %in% IS.genus$cor_char)  #67 overlaps of IS
  sum(clr.IR$cor_char %in% IR.genus$cor_char)  #45 overlaps of IR

###Sent on 02/07
load("Downloads/Revised_IntraOmes_rmcorr_clr_abundant.RData")
load("Downloads/IntraOme_microbiome_sparcc.RData")

  #Only genus pairwise: 990
  HC.st_clr.sub = HC.st_clr_abundant[grepl("genus_.*genus_.*", as.character(HC.st_clr_abundant$cor_char)), ]
  HC.st.clr.IS = HC.st_clr.sub[HC.st_clr.sub[, 6] < 0.05, ]  #130 obs.
  
  HC.st.sub = HC.st[grepl("genus_.*genus_.*", as.character(HC.st$cor_char)), ]
  HC.st.IS = HC.st.sub[HC.st.sub[, 6] < 0.05, ]  #312 obs.
  
  #From sparCC result
  spc.st.sub = data.frame(cor_char=HC.st.sub$cor_char, rho=genus_features_sparcc_cor[lower.tri(genus_features_sparcc_cor)])
  spc.st.sig = spc.st.sub[abs(spc.st.sub[, 2]) > 0.2, ]  #140 obs.
  
  sum(spc.st.sig$cor_char %in% HC.st.clr.IS$cor_char)  #53 out of 140
  sum(spc.st.sig$cor_char %in% HC.st.IS$cor_char)  #58 out of 312
  sum(HC.st.clr.IS$cor_char %in% HC.st.IS$cor_char)  #73 out of 312
  
###Sent on 02/06
  load("Downloads/Revised_IntraOmes_rmcorr_clr.RData")  #13,861 obs.

  st.sig.IS = HC.st_clr[HC.st_clr[, 6] < 0.05, 1:6]  #7,644 obs.
  st.sig.IR = HC.st_clr[HC.st_clr[, 12] < 0.05, 7:12]  #7,183 obs.
  
  st.IS.sub = st.sig.IS[st.sig.IS$cor_char %in% HC.st$cor_char, ]  #397 obs. compared to 166 IS.genus from old method
  