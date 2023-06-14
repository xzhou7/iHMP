
##====== Revised on 02/07/19: Plotting heatmaps of variables associated with SSPG after correcting BMI and age==========
#Importing data
load("Box Sync/PHI_Protected/Shared_WZ_HMP/Revision/Data/Revision_MultiOmes.RData") #Revised version shared on 04/18/18

##Adding possible confounding factors for SSPG correlation, TGL, TGL/HDL
  ds = clinic.df[, c(24, 48, 53, 58)]
  conff = {}
  for (subject in unique(ds$SubjectID)) {
    hset = ds[(ds$SubjectID == subject & ds$CL4 == "Healthy"), ]
    if (nrow(hset) != 0 ) {
      hvalue = apply(hset[, 1:2], 2, function(x) median(as.numeric(x), na.rm = TRUE))
      hset$THratio = as.numeric(hset[, 2])/as.numeric(hset[, 1])  #Calculate TGL/HDL ratio for each visit and then take the median
      ratio.m = median(hset$THratio, na.rm = TRUE)
      hrow = data.frame(subject, t(hvalue), Ratio=ratio.m)
      conff = rbind(conff, hrow)
    }
  }
  ##Adding medication records as the conff factors
  med = read.table("Box Sync/WenyuEx/Research_HMP/OmicsData/iPOP_MedicationUse.txt", stringsAsFactors = FALSE, sep = "\t", header = TRUE)
  med[, 2:3] <- apply(med[, 2:3], 2, function(x) sub("No", 0, sub("Yes", 1, x)))
  
  conff <- merge(conff, med, by.x = "subject", by.y = "SubjectID", all.x = TRUE)
  conff[, 5:6] <- apply(conff[, 5:6], 2, function(x) as.numeric(x))

#Build costom function (HMT, Healthy Median Test) to extract median value of healthy baselines per subject per variable
  #Proving ds, 1: m variable, then SubjectID, CL4
  library(ppcor)
  myHMT <- function(ds, m, fdr, trans){
    htb = {}
    for (subject in unique(ds$SubjectID)) {
      hset = ds[(ds$SubjectID == subject & ds$CL4 == "Healthy"), ]
      if (nrow(hset) != 0 ) {
        hvalue = apply(hset[, 1:m], 2, function(x) median(as.numeric(x), na.rm = TRUE))
        hrow = data.frame(subject, t(hvalue))
        htb = rbind(htb, hrow)
      }
    }
    htb.cl = merge(merge(htb, sc[, c(1,5,11,12)], by.x = "subject", by.y = "SubjectID"), conff, by = "subject") #Add classifications and confounding matrixes
    hdf = htb.cl[!is.na(htb.cl$SSPG), ]
    #Transform
    hdf$SSPG <- as.numeric(hdf$SSPG)
    if (trans == "none") {hdf[, 2:(m+1)] <- apply(hdf[, 2:(m+1)], 2, function(x) as.numeric(as.character(x)))} 
    if (trans == "log2") {hdf[, 2:(m+1)] <- apply(hdf[, 2:(m+1)], 2, function(x) log2(as.numeric(as.character(x))))} 
    if (trans == "arcsin") {hdf[, 2:(m+1)] <- apply(hdf[, 2:(m+1)], 2, function(x) asin(sqrt(as.numeric(as.character(x)))))}
    if (trans == "arcsin/1000") {hdf[, 2:(m+1)] <- apply(hdf[, 2:(m+1)], 2, function(x) asin(sqrt(as.numeric(as.character(x))/1000)))}
    #Starting correlating
    cor_char = cor_r = cor_p = c()
    conf_p = {}
    for (i in 2:(m+1)) { 
        m.ds = hdf[, c(i, (m+2), (m+3), (m+4), (m+5), (m+6), (m+7), (m+8), (m+9))] 
        del.row = unique(which(is.na(m.ds), arr.ind = TRUE)[, 1])
        if (length(del.row) != 0) {m.ds <- m.ds[-del.row, ]}  #deleting rows having NAs
        model <- pcor.test(m.ds[, 1], m.ds[, 2], m.ds[, 3:4]) #SSPG
        m.hdl <- tryCatch(pcor.test(m.ds[, 1], m.ds[, 5], m.ds[, 3:4]), error = function(e) NA) #HDL
        m.tgl <- tryCatch(pcor.test(m.ds[, 1], m.ds[, 6], m.ds[, 3:4]), error = function(e) NA) #TGL
        m.ratio <- pcor.test(m.ds[, 1], m.ds[, 7], m.ds[, 3:4]) #TGL/HDL ratio
        m.stat <- tryCatch(pcor.test(m.ds[, 1], m.ds[, 8], m.ds[, 3:4]), error = function(e) NA) #Statin
        m.glc <- tryCatch(pcor.test(m.ds[, 1], m.ds[, 9], m.ds[, 3:4]), error = function(e) NA) #Glucose medication
        cor_char <- c(cor_char, colnames(hdf)[i])
        cor_r <- c(cor_r, model[, 1])
        cor_p <- c(cor_p, model[, 2])
        conf_p <- rbind(conf_p, data.frame(HDL_p=tryCatch(m.hdl[, 2], error = function(e) NA), TGL_p=tryCatch(m.tgl[, 2], error = function(e) NA), ratio_P=m.ratio[, 2],
                        Statin_p=tryCatch(m.stat[, 2], error = function(e) NA), GLC_p=tryCatch(m.glc[, 2], error = function(e) NA)))
    }
    conf_p_adj = apply(conf_p, 2, function(x) p.adjust(x, method = "fdr"))
    cor_df = data.frame(cor_char, cor_r, cor_p, conf_p_adj)  
    cor_df$adj.p = p.adjust(cor_p, method = "fdr")
    #Filtering significance
    sig.var = as.character(cor_df$cor_char)[!is.na(cor_df$adj.p) & (cor_df$adj.p < fdr)]   # Using FDR to filter
    htb.sig = data.frame(SubjectID=hdf$subject, hdf[, sig.var], SSPG=hdf$SSPG)
    if (length(sig.var) == 1) {colnames(htb.sig)[2] <- sig.var}
    #Below are adding p-value and rho for each sig.var
    htb.sig$SubjectID <- as.character(htb.sig$SubjectID)
    htb.sig[nrow(htb.sig)+1, ] <- c("p_ad", cor_df$adj.p[cor_df$cor_char %in% sig.var], "NA") #p.ad values
    htb.sig[nrow(htb.sig)+1, ] <- c("r", cor_df$cor_r[cor_df$cor_char %in% sig.var], "NA") #rho values
    htb.sig[c((nrow(htb.sig)+1):(nrow(htb.sig)+5)), ] <- cbind(c("HDL_p", "TGL_p", "ratio_p", "Statin_p", "GLC_p"), t(cor_df[cor_df$cor_char %in% sig.var, 4:8]), rep("NA", 5))
    return(htb.sig)
  } #End of definiting custom function
  
  #Host 11/06
  hmt.ck = myHMT(ck.df[, c(4:65, 70, 75)], 62, 0.1, "none") #1 at 10%
  hmt.clinic = myHMT(clinic.df[, c(2:53, 58)], 51, 0.1, "none") #4 for 10% FDR
  hmt.metb = myHMT(metbcr.df[, c(2:726, 731)], 724, 0.1, "log2")  #5 for 10% FDR
  hmt.rna = myHMT(rnaseq.log.df[, c(2:10348, 10353)], 10346, 0.2, "none") # 0 for 20% FDR
  hmt.prot = myHMT(swathprot.PCR.df[, c(2:303, 304, 309)], 302, 0.2, "none") #6 for 10% FDR
  #Microbiome
  hmt.st = myHMT(st.df[, c(2:97, 103, 102)], 96, 0.1, "arcsin") #for 10% FDR, 27 taxonomic levels
  hmt.ns = myHMT(ns.df[, c(2:81, 87, 86)], 80, 0.2, "arcsin") #0 for 20% FDR
  hmt.stKO = myHMT(stKO.df[, c(2:362, 368, 367)], 361, 0.1, "arcsin/1000") #for 10% FDR, 49 KOs
  hmt.nsKO = myHMT(nsKO.df[, c(2:358, 364, 363)], 357, 0.2, "arcsin/1000") #for 20% FDR, none
  
  #Visualzing SSPG association at baseline
  #Merge first, clinic, metb, ns have 62 subjects, st have 59 subjects
  hmt.tb = merge(merge(merge(merge(merge(hmt.clinic, hmt.ck[, -3], by = "SubjectID", all = TRUE), hmt.metb[, -6], by = "SubjectID", all = TRUE), hmt.st[, -29], by = "SubjectID", all = TRUE), hmt.prot[, -8], by = "SubjectID", all = TRUE), hmt.stKO[, -44], by = "SubjectID", all = TRUE)
  hmt.tb <- hmt.tb[order(as.numeric(as.character(hmt.tb$SSPG))), ]
  write.table(t(hmt.tb), file = "Desktop/Revised_SSPG_HealthyMedian.txt", sep = '\t')
  
  #Taking gut microbes at genus level and no stKO
  pd = t(apply(hmt.tb[1:63, c(2:6, 8:12, 26:45)], 2, function(x) scale(as.numeric(x))))
  
  #Replacing metabolite name
  for (i in 1:63) {
    if (grepl("pHILIC|nHILIC|pRPLC|nRPLC", rownames(pd)[i])) {rownames(pd)[i] <- metb.curated[metb.curated$Compounds_ID == rownames(pd)[i], "Metabolite"]}
  }
  
  var_annotation = c(rep(c("pink"), times=5), rep(c("red"), times=1), rep(c("grey"), times=4), rep(c("darkgreen"), times=14), rep(c("orange"), times=6))
  
  library("gplots")
  library("RColorBrewer")
  mycol = rev(colorRampPalette(brewer.pal(11, "RdBu"))(53))
  
  heatmap.2(pd, Rowv = TRUE, Colv = FALSE, key = TRUE, keysize = 1.5,
            scale = "none", trace = "none", dendrogram = "row",  
            distfun = function(x) dist(x,method = "manhattan"), hclustfun = function(x) hclust(x,method = "complete"),
            RowSideColors=var_annotation, margins = c(1,14),
            col = mycol, cexRow = 1, labCol = FALSE, offsetRow = 0, breaks = seq (-5.3, 5.3, by = 0.2),
            lwid = c(0.3, 4), lhei = c(0.05, 4))
  
  library(ggplot2)
  pd.sspg = data.frame(hmt.tb[1:63,c("SubjectID","SSPG")])
  pd.sspg$SubjectID <- factor(pd.sspg$SubjectID, levels =  pd.sspg$SubjectID)
  ggplot(data=pd.sspg, aes(x = SubjectID, y = as.numeric(as.character(SSPG)))) + geom_point() +
    theme(axis.text.x = element_blank(), panel.grid.major = element_blank(), 
          panel.background = element_blank(), axis.line.y = element_line(colour = "black"))+
    labs(x="", y = "SSPG") + scale_y_continuous(position = "right")
  ggsave("Desktop/Revised3_SSPG_HealthyMedian_sspg.pdf", width = 8, height = 1)
  