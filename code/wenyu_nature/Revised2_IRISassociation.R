

##====== Revised 2nd round: Plotting heatmaps of variables associated with IRIS after correcting BMI, age and gender==========
  #Importing data
  load("Box Sync/PHI_Protected/Shared_WZ_HMP/Revision/Data/Revision_MultiOmes.RData") #Revised version shared on 04/18/18
  
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
    htb.cl = merge(htb, sc[, c(1,4,11,12)], by.x = "subject", by.y = "SubjectID") #Add classifications
    hdf = htb.cl[htb.cl$IRIS != "Unknown", ]
    hdf$IRIS <- as.numeric(gsub("IR", 1, gsub("IS", 0, hdf$IRIS))) #Converting IR/IS to dummy variables, IS 0
    #Transform
    if (trans == "none") {hdf[, 2:(m+1)] <- apply(hdf[, 2:(m+1)], 2, function(x) as.numeric(as.character(x)))} 
    if (trans == "log2") {hdf[, 2:(m+1)] <- apply(hdf[, 2:(m+1)], 2, function(x) log2(as.numeric(as.character(x))))} 
    if (trans == "arcsin") {hdf[, 2:(m+1)] <- apply(hdf[, 2:(m+1)], 2, function(x) asin(sqrt(as.numeric(as.character(x)))))}
    if (trans == "arcsin/1000") {hdf[, 2:(m+1)] <- apply(hdf[, 2:(m+1)], 2, function(x) asin(sqrt(as.numeric(as.character(x))/1000)))}
    #Starting correlating
    cor_char = cor_r = cor_p = c()
    for (i in 2:(m+1)) { 
      m.ds = hdf[, c(i, (m+2), (m+3), (m+4))] 
      del.row = unique(which(is.na(m.ds), arr.ind = TRUE)[, 1])
      if (length(del.row) != 0) {m.ds <- m.ds[-del.row, ]}  #deleting rows having NAs
      model <- pcor.test(m.ds[, 1], m.ds[, 2], m.ds[, 3:4])
      cor_char <- c(cor_char, colnames(hdf)[i])
      cor_r <- c(cor_r, model[, 1])
      cor_p <- c(cor_p, model[, 2])
    }
    cor_df = data.frame(cor_char, cor_r, cor_p)
    cor_df$adj.p = p.adjust(cor_p, method = "fdr")
    #Filtering significance
    sig.var = as.character(cor_df$cor_char)[!is.na(cor_df$adj.p) & (cor_df$adj.p < fdr)]   # Using FDR to filter
    htb.sig = data.frame(SubjectID=hdf$subject, hdf[, sig.var], IRIS=hdf$IRIS)
    if (length(sig.var) == 1) {colnames(htb.sig)[2] <- sig.var}
    #Below are adding p-value and rho for each sig.var
    htb.sig$SubjectID <- as.character(htb.sig$SubjectID)
    htb.sig[nrow(htb.sig)+1, ] <- c("p_ad", cor_df$adj.p[cor_df$cor_char %in% sig.var], "NA") #p.ad values
    htb.sig[nrow(htb.sig)+1, ] <- c("r", cor_df$cor_r[cor_df$cor_char %in% sig.var], "NA") #rho values
    return(htb.sig)
  } #End of definiting custom function
  
  #Host
  hmt.ck = myHMT(ck.df[, c(4:65, 70, 75)], 62, 0.2, "none") #0
  hmt.clinic = myHMT(clinic.df[, c(2:53, 58)], 51, 0.2, "none") #14 under 20% FDR, 6 for 10% FDR
  hmt.metb = myHMT(metbcr.df[, c(2:726, 731)], 724, 0.2, "log2")  #10 under 20% FDR, 5 for 10% FDR
  hmt.rna = myHMT(rnaseq.log.df[, c(2:10348, 10353)], 10346, 0.2, "none") # 0 for 20% FDR
  hmt.prot = myHMT(swathprot.PCR.df[, c(2:303, 304, 309)], 302, 0.2, "none") #3 under 20% FDR
  #Microbiome
  hmt.st = myHMT(st.df[, c(2:97, 103, 102)], 96, 0.2, "arcsin") #for 31 under 20% FDR, 10% FDR, 14 taxonomic levels
  hmt.ns = myHMT(ns.df[, c(2:81, 87, 86)], 80, 0.2, "arcsin") #0 for 20% FDR
  hmt.stKO = myHMT(stKO.df[, c(2:362, 368, 367)], 361, 0.2, "arcsin/1000") #94 under 20% FDR, 10% FDR, 47 KOs
  hmt.nsKO = myHMT(nsKO.df[, c(2:358, 364, 363)], 357, 0.2, "arcsin/1000") #for 20% FDR, none
  
  #Visualzing IRIS association at baseline
  #Merge first, clinic, metb, ns have 62 subjects, st have 59 subjects
  hmt.tb = merge(merge(merge(merge(hmt.clinic, hmt.metb[, -12], by = "SubjectID", all = TRUE), hmt.st[, -33], by = "SubjectID", all = TRUE), hmt.prot[, -5], by = "SubjectID", all = TRUE), hmt.stKO[, -96], by = "SubjectID", all = TRUE)

  write.table(t(hmt.tb), file = "Desktop/Revised2_IRIS_HealthyMedian.txt", sep = '\t')
  
  #Taking gut microbes at genus level and no stKO
  pd = t(apply(hmt.tb[1:65, c(2:15, 17:26, 41:51)], 2, function(x) scale(as.numeric(x))))
  #Replacing metabolite name
  for (i in 1:65) {
    if (grepl("pHILIC|nHILIC|pRPLC|nRPLC", rownames(pd)[i])) {rownames(pd)[i] <- metb.curated[metb.curated$Compounds_ID == rownames(pd)[i], "Metabolite"]}
  }
  
  IRIS_annotation = sub("0", "#90AA3C", sub("1", "#EF6125", hmt.tb$IRIS[1:65]))[c(-6, -39)]
  var_annotation = c(rep(c("pink"), times=14), rep(c("grey"), times=10), rep(c("darkgreen"), times=11))
  
  library("gplots")
  library("RColorBrewer")
  mycol = rev(colorRampPalette(brewer.pal(11, "RdBu"))(50))
  
  heatmap.2(pd[, c(-6, -39)], Rowv = TRUE, Colv = TRUE, key = TRUE, keysize = 1.5,
            scale = "none", trace = "none", dendrogram = "none",  
            distfun = function(x) dist(x,method = "manhattan"), hclustfun = function(x) hclust(x,method = "complete"),
            ColSideColors=IRIS_annotation, RowSideColors=var_annotation,
            col = mycol, cexRow = 0.8, labCol = FALSE, offsetRow = 0, breaks = seq ( -5, 5, by = 0.2 ),
            lwid = c(0.02, 4), lhei = c(0.05, 2.5))
  
  
 