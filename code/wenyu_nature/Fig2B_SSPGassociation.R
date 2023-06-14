
##====== Plotting heatmaps of variables associated with SSPG==========
#Importing data
load("Box Sync/PHI_Protected/Shared_WZ_HMP/OmicsData/HMP_MultiOmics.RData")
load("Box Sync/PHI_Protected/Shared_WZ_HMP/OmicsData/HMP_Microbiome_Perc.RData")

  #Preparing data
  metbcr.df <- metb.df[, c(1, c(1:1798)[colnames(metb.df) %in% metb.curated$Compounds_ID], 1793:1798)] #Only 724 curated metabolites
  #Nasal
  pr = apply(r16s.ns.perc.df[, c(2:954)], 2, function(x) sum(x == 0))
  ns.df <- r16s.ns.perc.df[, c(1, c(2:954)[pr < 400],3807:3812)] #Narrowing to 80 taxa levels that were present in more than half sample size
  pr = apply(r16sKO.ns.perc.df[, c(2:6910)], 2, function(x) sum(x < 1))
  nsKO.df <- r16sKO.ns.perc.df[, c(1, c(2:6910)[pr < 400],6911:6916)] #Narrowing to 357 genes that have more than 1% in more than half sample size
  #Stool
  pr = apply(r16s.st.perc.df[, c(2:322)], 2, function(x) sum(x == 0))
  st.df <- r16s.st.perc.df[, c(1, c(2:322)[pr < 430],2276:2281)] #Narrowing to 96 taxa levels
  pr = apply(r16sKO.st.perc.df[, c(2:6910)], 2, function(x) sum(x < 1))
  stKO.df <- r16sKO.st.perc.df[, c(1, c(2:6910)[pr < 400],6911:6916)] #Narrowing to 361 genes that have more than 1% in more than half sample size

  #Build costom function (HMT, Healthy Median Test) to extract median value of healthy baselines per subject per variable
  #Proving ds, 1: m variable, then SubjectID, CL4
  library(Hmisc)
  myHMT <- function(ds, m, fdr){
    htb = {}
    for (subject in unique(ds$SubjectID)) {
      hset = ds[(ds$SubjectID == subject & ds$CL4 == "Healthy"), ]
      if (nrow(hset) != 0 ) {
        hvalue = apply(hset[, 1:m], 2, function(x) median(as.numeric(x), na.rm = TRUE))
        hrow = data.frame(subject, t(hvalue))
        htb = rbind(htb, hrow)
      }
    }
    htb.cl = merge(htb, sc[, c(1,3)], by.x = "subject", by.y = "SubjectID") #Add classifications
    htb.cl.sub = htb.cl[!is.na(htb.cl$SSPG), ]
    #Spearman correlations
    n = m+1 
    p.val = rho = {}
    for (i in 2:n) { 
      p.val[colnames(htb.cl.sub)[i]] = rcorr(htb.cl.sub[, i], htb.cl.sub[, "SSPG"], type = "spearman")$P[1,2]
      rho[colnames(htb.cl.sub)[i]] = rcorr(htb.cl.sub[, i], htb.cl.sub[, "SSPG"], type = "spearman")$r[1,2]
    }
    p.ad <- p.adjust(p.val, method = "fdr")
    sig.var = names(p.ad)[!is.na(p.ad) & p.ad < fdr]   # Using FDR to filter
    htb.sig = data.frame(SubjectID=htb.cl.sub$subject, htb.cl.sub[, sig.var], SSPG=htb.cl.sub$SSPG)
    #Below are adding p-value and rho for each sig.var
    htb.sig$SubjectID <- as.character(htb.sig$SubjectID)
    htb.sig[nrow(htb.sig)+1, ] <- c("p_ad", p.ad[sig.var], "NA") #p.ad values
    htb.sig[nrow(htb.sig)+1, ] <- c("rho", rho[sig.var], "NA") #rho values
    return(htb.sig)
  } #End of definiting custom function
  
  #Host
  hmt.ck = myHMT(ck.df[, c(4:65, 70, 75)], 62, 0.2) #0 for 20% FDR
  hmt.clinic = myHMT(clinic.df[, c(2:53, 58)], 51, 0.1) #8 for 10% FDR
  hmt.metb = myHMT(metbcr.df[, c(2:726, 731)], 724, 0.1)  #36 for 10% FDR
  hmt.rna = myHMT(rnaseq.log.df[, c(2:10348, 10353)], 10346, 0.2) # 0 for 20% FDR
  hmt.prot = myHMT(swathprot.df[, c(2:307, 309, 314)], 306, 0.2) #0 for 20% FDR
  #Microbiome
  hmt.st = myHMT(st.df[, c(2:97, 103, 102)], 96, 0.01) #for 1% FDR, 18 taxonomic levels
  hmt.stKO = myHMT(stKO.df[, c(2:362, 368, 367)], 361, 0.01) #for 1% FDR, 27 KOs
  hmt.ns = myHMT(ns.df[, c(2:81, 87, 86)], 80, 0.1) #for 10% FDR, 10 taxonomic levels
  hmt.nsKO = myHMT(nsKO.df[, c(2:358, 364, 363)], 357, 0.2) #for 20% FDR, none
  
  #Visualzing SSPG association at baseline
  #Merge first, clinic, metb, ns have 62 subjects, st have 59 subjects
  hmt.tb = merge(merge(merge(merge(hmt.clinic, hmt.ns[, -12], by = "SubjectID"), hmt.metb[, -38], by = "SubjectID", all.x = TRUE), hmt.st[, -20], by = "SubjectID", all.x = TRUE), hmt.stKO[, -29], by = "SubjectID", all.x = TRUE)
  hmt.tb <- hmt.tb[order(as.numeric(as.character(hmt.tb$SSPG))), ]
  write.table(t(hmt.tb), file = "Desktop/Fig2B_SSPG_HealthyMedian.txt", sep = '\t')
  
  pd = t(apply(hmt.tb[1:62, c(2:9, 11:101)], 2, function(x) scale(as.numeric(x))))
  rownames(pd) = colnames(hmt.tb)[c(2:9, 11:101)]
  #Replacing metabolite name
  for (i in 1:99) {
    if (grepl("pHILIC|nHILIC|pRPLC|nRPLC", rownames(pd)[i])) {rownames(pd)[i] <- metb.curated[metb.curated$Compounds_ID == rownames(pd)[i], "Metabolite"]}
  }
  
  var_annotation = c(rep(c("pink"), times=8), rep(c("orange"), times=10), rep(c("grey"), times=36), rep(c("darkgreen"), times=18), rep(c("lightgreen"), times=27))
  
  library("gplots")
  library("RColorBrewer")
  mycol = rev(colorRampPalette(brewer.pal(11, "RdBu"))(69))
  
  heatmap.2(pd, Rowv = TRUE, Colv = FALSE, key = TRUE, keysize = 1.5,
            scale = "none", trace = "none", dendrogram = "row",  
            distfun = function(x) dist(x,method = "manhattan"), hclustfun = function(x) hclust(x,method = "complete"),
            RowSideColors=var_annotation,
            col = mycol, cexRow = 0.4, labCol = FALSE, offsetRow = 0, breaks = seq ( -6.9, 6.9, by = 0.2 ),
            lwid = c(0.1, 4), lhei = c(0.05, 4))
  
  pd.sspg = data.frame(hmt.tb[1:62,c("SubjectID","SSPG")])
  pd.sspg$SubjectID <- factor(pd.sspg$SubjectID, levels =  pd.sspg$SubjectID)
  ggplot(data=pd.sspg, aes(x = SubjectID, y = as.numeric(as.character(SSPG)))) + geom_point() +
    theme(axis.text.x = element_blank(), panel.grid.major = element_blank(), 
          panel.background = element_blank(), axis.line.y = element_line(colour = "black"))+
    labs(x="", y = "SSPG") + scale_y_continuous(position = "right")
  ggsave("Desktop/Fig2B_SSPG_HealthyMedian_sspg.pdf", width = 7, height = 1)
  
##====== Plotting heatmaps of variables associated with IS/IR
  #Build costom function (HMT, Healthy Median Test) to extract median value of healthy baselines per subject per variable
  #Proving ds, 1: m variable, then SubjectID, CL4
  library(Hmisc)
  myHMT <- function(ds, m, fdr){
    htb = {}
    for (subject in unique(ds$SubjectID)) {
      hset = ds[(ds$SubjectID == subject & ds$CL4 == "Healthy"), ]
      if (nrow(hset) != 0 ) {
        hvalue = apply(hset[, 1:m], 2, function(x) median(as.numeric(x), na.rm = TRUE))
        hrow = data.frame(subject, t(hvalue))
        htb = rbind(htb, hrow)
      }
    }
    htb.cl = merge(htb, sc[, c(1,2)], by.x = "subject", by.y = "SubjectID") #Add IR/IS classifications
    htb.cl.sub = htb.cl[htb.cl$IRIS != "Unknown", ]
    IR.ind = which(htb.cl.sub$IRIS == "IR")
    IS.ind = which(htb.cl.sub$IRIS == "IS")
    #Wilcoxon signed-rank test
    n = m+1 
    p.val = rho = {}
    for (i in 2:n) { 
      p.val[colnames(htb.cl.sub)[i]] = wilcox.test(htb.cl.sub[, i][IR.ind], htb.cl.sub[, i][IS.ind])$p.value
      }
    p.ad <- p.adjust(p.val, method = "fdr")
    sig.var = names(p.ad)[!is.na(p.ad) & p.ad < fdr]   # Using FDR to filter
    htb.sig = data.frame(SubjectID=htb.cl.sub$subject, htb.cl.sub[, sig.var], IRIS=htb.cl.sub$IRIS)
    #Below are adding p-value and rho for each sig.var
    htb.sig$SubjectID <- as.character(htb.sig$SubjectID)
    htb.sig[nrow(htb.sig)+1, ] <- c("p_ad", p.ad[sig.var], "NA") #p.ad values
    return(htb.sig)
  } #End of definiting custom function
  
  #Host
  hmt.ck = myHMT(ck.df[, c(4:65, 70, 75)], 62, 0.2) #0 for 20% FDR
  hmt.clinic = myHMT(clinic.df[, c(2:53, 58)], 51, 0.2) #13 for 20% FDR
  hmt.metb = myHMT(metbcr.df[, c(2:726, 731)], 724, 0.2)  #44 for 10% FDR
  hmt.rna = myHMT(rnaseq.log.df[, c(2:10348, 10353)], 10346, 0.2) # 0 for 20% FDR
  hmt.prot = myHMT(swathprot.df[, c(2:307, 309, 314)], 306, 0.2) #0 for 20% FDR
  #Microbiome
  hmt.st = myHMT(st.df[, c(2:97, 103, 102)], 96, 0.05) #for 1% FDR, 21 taxonomic levels
  hmt.stKO = myHMT(stKO.df[, c(2:362, 368, 367)], 361, 0.05) #for 1% FDR, 28 KOs
  hmt.ns = myHMT(ns.df[, c(2:81, 87, 86)], 80, 0.2) #for 10% FDR, 11 taxonomic levels
  hmt.nsKO = myHMT(nsKO.df[, c(2:358, 364, 363)], 357, 0.2) #for 20% FDR, none
  
  #Merge first, clinic, metb, ns have 62 subjects, st have 59 subjects
  hmt.tb = merge(merge(merge(merge(hmt.clinic, hmt.ns[, -13], by = "SubjectID"), hmt.metb[, -46], by = "SubjectID", all.x = TRUE), hmt.st[, -23], by = "SubjectID", all.x = TRUE), hmt.stKO[, -30], by = "SubjectID", all.x = TRUE)
  write.table(t(hmt.tb), file = "Desktop/Fig2d_IRIS_HealthyMedian.txt", sep = '\t')
  
  pd = t(apply(hmt.tb[1:62, c(2:14, 16:119)], 2, function(x) scale(as.numeric(x))))
  #Replacing metabolite name
  for (i in 1:117) {
    if (grepl("pHILIC|nHILIC|pRPLC|nRPLC", rownames(pd)[i])) {rownames(pd)[i] <- metb.curated[metb.curated$Compounds_ID == rownames(pd)[i], "Metabolite"]}
  }
  #Clinic:pink, nasal: orange, metb: grey, st:darkgreen, stKO:lightgreen
  var_annotation = c(rep(c("pink"), times=13), rep(c("orange"), times=11), rep(c("grey"), times=44), rep(c("darkgreen"), times=21), rep(c("lightgreen"), times=28))
  IRIS_annotation = gsub("IS", "#90AA3C", gsub("IR", "#EF6125", hmt.tb$IRIS[1:62]))
  
  library("gplots")
  library("RColorBrewer")
  mycol = rev(colorRampPalette(brewer.pal(11, "RdBu"))(69))
  
  heatmap.2(pd, Rowv = TRUE, Colv = TRUE, key = TRUE, keysize = 1.5,
            scale = "none", trace = "none", dendrogram = "row",  
            distfun = function(x) dist(x,method = "euclidean"), hclustfun = function(x) hclust(x,method = "complete"),
            RowSideColors=var_annotation, ColSideColors=IRIS_annotation,
            col = mycol, cexRow = 0.3, labCol = FALSE, offsetRow = 0, breaks = seq ( -6.9, 6.9, by = 0.2 ),
            lwid = c(0.1, 4), lhei = c(0.05, 4))
  
  