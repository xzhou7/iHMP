

###Inspecting A1C association difference between IR/IS
load("Box Sync/WenyuEx/Rcodes/Revised_rmcorr_labs.RData")

#Calculating InterCor for A1C, GLU and HSCRP
  setlab.df = clinic.df[, c(1, 2, 22, 26)]
  setlab.df$HSCRP <- as.numeric(setlab.df$HSCRP)
  lab_metb.df <- merge(setlab.df, metbcr.df, by = "SampleID")
  HPC.labs.metb = myHPC(lab_metb.df[, c(5:728, 2:4, 729, 734)], 724, 728, "log2") #5:728 metb 
  
  lab_prot.df <- merge(setlab.df, swathprot.PCR.df, by = "SampleID")
  HPC.labs.prot = myHPC(lab_prot.df[, c(5:306, 2:4, 307, 312)], 302, 306, "none") #5:306 proteins
 
  lab_ck.df <- merge(setlab.df, ck.df[, c(1,4:75)], by = "SampleID")
  HPC.labs.ck = myHPC(lab_ck.df[, c(5:66, 2:4, 71, 76)], 62, 66, "none") #5:66 cytokines
  
  lab_rna.df <- merge(setlab.df, rnaseq.log.df, by = "SampleID")
  HPC.labs.rna = myHPC(lab_rna.df[, c(5:10350, 2:4, 10351, 10356)], 10346, 10350, "none") #5:10350 proteins
  
  ds = merge(setlab.df, clinic.df[, c(1, 3:21, 23:25, 27:28, 30:50, 52, 53, 58)], by = "SampleID") #No 29, 51
  HPC.labs = myHPC(ds[, c(5:50, 2:4,51,52)], 46, 50, "none")
  
  lab_st.df <- merge(setlab.df, st.df, by.x = "SampleID", by.y = "HostSampleID")
  HPC.labs.st = myHPC(lab_st.df[, c(5:100, 2:4, 106, 105)], 96, 100, "arcsin") #5:100 gut micorbes
  
  save(HPC.labs, HPC.labs.ck, HPC.labs.metb, HPC.labs.prot, HPC.labs.rna, HPC.labs.st, file = "Box Sync/WenyuEx/Rcodes/Revised_InterCor_labs.RData")
  

#Barplot with error bars, showing the average edge per analyte (epa)
  #Preparing matrix for IS or IR specific ones
     var = "A1C"
     IS.sig = IR.sig = {}
     for (df in c("HC.labs", "HC.labs.ck", "HC.labs.metb", "HC.labs.prot", "HC.labs.rna", "HC.labs.st")){
       ds = get(df)
       #For IS
       ds.sig = ds[(ds[, 6] < 0.05 & ds[, 12] > 0.05 & grepl(var, ds[, 1])), 1:12]
       IS.sig <- rbind(IS.sig, ds.sig)
       #For IR
       ds.sig = ds[(ds[, 12] < 0.05 & ds[, 6] > 0.05 & grepl(var, ds[, 7])), 1:12]
       IR.sig <- rbind(IR.sig, ds.sig)
     }
     tb.sig = rbind(IS.sig, IR.sig)
     tb.sig$Variable = gsub(paste(" vs ", var, sep=""), "", tb.sig$cor_char)
     tb.sig <- merge(tb.sig, metb.curated[, c(2:5, 14, 15)], by.x = "Variable", by.y = "Compounds_ID", all.x = TRUE)
     tb.sig$Variable[2:14] <- tb.sig$Metabolite[2:14]  #For A1C #custome for differnt var
    #For GLU: tb.sig$Variable[2:6] <- tb.sig$Metabolite[2:6] 
    #         tb.sig <- tb.sig[-3, ]  #duplicated Hexosamine
    #For HSCRP: tb.sig$Variable[12:15] <- tb.sig$Metabolite[12:15] 
     
#Visualization for IS or IR specific ones
     pd = as.data.frame(rbind(as.matrix(tb.sig[, c(1, 3:7)]), as.matrix(tb.sig[, c(1, 9:13)])))
     pd$ISIR = c(rep("IS", nrow(tb.sig)), rep("IR", nrow(tb.sig)))
     pd$ISIR <- factor(pd$ISIR, levels = c("IS", "IR"))
     pd$Variable <- factor(pd$Variable, levels = as.character(tb.sig[order(tb.sig[, 7], decreasing = TRUE), "Variable"]))
     pd[, 2:6] <- apply(pd[, 2:6], 2, function(x) as.numeric(as.character(x)))
     
     ggplot(data=pd, aes(x=Variable, y=cor_r, fill=ISIR)) + geom_ribbon(aes(ymin=CI_left, ymax=CI_right, color=ISIR), size = 3, position=position_dodge(0.6)) + 
       geom_point(position = position_dodge(width=0.6), size = 3) +
       geom_hline(yintercept = 0, size = 0.3) + coord_flip() + scale_color_manual(breaks = c("IS", "IR"), values=c("#cddb9d", "#ffc2a8")) +
       theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             panel.background = element_blank(), axis.line.x = element_line(colour = "black", size =1),
             axis.line.y = element_blank(), axis.text = element_text(size = 12)) + labs(x="",y="Correlation Coefficient r", color="") + guides(fill=FALSE)
     ggsave("Desktop/Revised_Fig5_cor_GLU.pdf", height = 6, width = 6)
     
     #For InterCor graph, first matching variables in tb.sig
     #no associations with RNAs
     tb.HPC = {}
     for (df in c("HPC.labs", "HPC.labs.ck", "HPC.labs.metb", "HPC.labs.prot", "HPC.labs.st")){
       ds = get(df)
       match = ds[(ds[, 1] %in% tb.sig$cor_char), ]
       tb.HPC = rbind(tb.HPC, match)
     }
     #Prepare dataset
     pd_match = rbind(tb.HPC[, c(1,2,4)], tb.HPC[, c(5,6,8)], tb.HPC[, c(9,10,12)])
     pd_match$ISIR = c(rep("Combined", nrow(tb.HPC)), rep("IS", nrow(tb.HPC)), rep("IR", nrow(tb.HPC)))
     pd_match$ISIR <- factor(pd_match$ISIR, levels=c("Combined", "IS", "IR"))
     pd_match$Variable = gsub(" vs HSCRP", "", as.character(pd_match$cor_char))
     pd_match <- merge(pd_match, metb.curated[, c(2:5, 14, 15)], by.x = "Variable", by.y = "Compounds_ID", all.x = TRUE)
     pd_match$Variable[34:45] <- pd_match$Metabolite[34:45] 
     pd_match$Variable <- factor(pd_match$Variable, levels = levels(pd$Variable))
       
     ggplot(data=pd_match, aes(x=Variable, y=cor_r, fill = ISIR)) + geom_bar(position = "dodge", stat = "identity") +
       scale_fill_manual(breaks = c("Combined", "IS", "IR"), values=c("grey80", "#90AA3C", "#EF6125")) + coord_flip() + geom_hline(yintercept = 0, size = 0.3) +
       theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             panel.background = element_blank(), axis.line.x = element_line(colour = "black", size =1),
             axis.line.y = element_blank(), axis.text = element_text(size = 12)) + 
       labs(x="", y="Correlation Coefficient r", fill='')
     ggsave("Desktop/Revised_HSCRP_InterCor.pdf", height = 6, width = 6)
      