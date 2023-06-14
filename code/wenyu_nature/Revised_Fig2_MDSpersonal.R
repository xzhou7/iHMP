

##====== Plotting MDS graphs to access personal separation==========
#Loading dataset
load("Box Sync/PHI_Protected/Shared_WZ_HMP/Revision/Data/Revision_MultiOmes.RData") #Revised version shared on 04/18/18

#Build a custom function (HBC, Healthy Baseline Clutering) that prepare data for MDS plot
#Proving a ds, with 1:m variables, SubjectID and CL4, limiting subjects with more than n baseline visits
library(vegan)
library(MASS)
  myHBC <- function(ds, m, n, trans){
    tb = {}
    for (subject in unique(ds$SubjectID)) {
      subset = ds[(ds$CL4 == "Healthy" & ds$SubjectID == subject), ]
      if (nrow(subset) > n) {
        tb = rbind(tb, subset)
      }
    }
    if (trans == "none") {tb.scl = scale(apply(tb[, 1:m], 2, function(x) as.numeric(x)))}
    if (trans == "log10") {tb.scl = scale(apply(tb[, 1:m], 2, function(x) log10(as.numeric(x))))}
    if (trans == "arcsin") {tb.scl = scale(apply(tb[, 1:m], 2, function(x) asin(sqrt(as.numeric(x)))))}
    if (trans == "arcsin/1000") {tb.scl = scale(apply(tb[, 1:m], 2, function(x) asin(sqrt(as.numeric(x)/1000))))}
    X <- metaMDSdist(tb.scl, distance = "manhattan", autotransform = FALSE, noshare = FALSE, trace = 1, zerodist = "add", na.rm = TRUE)
    mds <- metaMDS(X, k = 2, try = 50, trymax = 200, wascores = FALSE, expand = FALSE, trace = 1, plot = FALSE)
    mdsdata = data.frame(mds$points, tb$SubjectID)
    return(mdsdata)
  } #End of custom function
  
#Host
  mds.ck = myHBC(ck.df[, c(4:65,70, 75)], 62, 9, "none") #11 subjects, solution reached
  mds.metb = myHBC(metbcr.df[, c(2:726, 731)], 724, 9, "log10") #10 subjects, no convergence
  mds.rna = myHBC(rnaseq.log.df[, c(2:10348, 1, 10353)], 10346, 8, "none") #10 subjects, no convergence
  mds.clinic = myHBC(clinic.df[, c(2:53, 58)], 51, 9, "none") #12 subjects, solution reached
  mds.prot = myHBC(swathprot.PCR.df[, c(2:303,304, 309)], 302, 9, "none") #11 subjects, no convergence
#Microbiome
  mds.st = myHBC(st.df[, c(2:97, 102:103)], 96, 9, "arcsin") #11 subjects, no convergence
  mds.ns = myHBC(ns.df[, c(2:81, 86:87)], 80, 8, "arcsin") #14 subjects, solution reached
  mds.stKO = myHBC(stKO.df[, c(2:362, 367:368)], 361, 9, "arcsin/1000") #11 subjects, no convergence
  mds.nsKO = myHBC(nsKO.df[, c(2:358, 363:364)], 357, 8, "arcsin/1000") #13 subjects, no convergence

#Top 30, 100, 300 and 1000 most personal variables
  load("Box Sync/WenyuEx/Rcodes/Revised_VarianceDecomposition.RData")
  
  #Top 30
  dcomp = vd.comb[order(vd.comb$ICC, decreasing = TRUE), ]
  pstb = merge(merge(merge(clinic.df, ck.df[, c(1, 4:65)], by = "SampleID"), metbcr.df[, 1:725], by = "SampleID"), swathprot.PCR.df[, 1:303],  by = "SampleID") #846 obs. of 1146 variables
  pstb.sub = pstb[, c(1, (1:1146)[colnames(pstb) %in% rownames(dcomp)[1:30]], 53, 58)]
  pstb.sub[, 22:30] <- apply(pstb.sub[, 22:30], 2, function(x) log10(x))
  mds.top30 = myHBC(pstb.sub[, 2:33], 30, 6, "none") #31 subjects by length(unique(mds.top30$tb.SubjectID)), solution reached

  #Top 100
  pstb = merge(merge(merge(merge(merge(clinic.df, ck.df[, c(1, 4:65)], by = "SampleID"), metbcr.df[, 1:725], by = "SampleID"), swathprot.PCR.df[, 1:303],  by = "SampleID"), ns.df[, 1:81], by.x = "SampleID", by.y="HostSampleID"), nsKO.df[, 1:358], by.x = "SampleID", by.y="HostSampleID") #659 obs. 1583 variables
  pstb.sub = pstb[, c(1, (1:1583)[colnames(pstb) %in% rownames(dcomp)[1:100]], 53, 58)]
  pstb.sub[, 48:92] <- apply(pstb.sub[, 48:92], 2, function(x) log10(x))
  pstb.sub[, 97:99] <- apply(pstb.sub[, 97:99], 2, function(x) asin(sqrt(x)))
  mds.top100 = myHBC(pstb.sub[, 2:101], 98, 5, "none") #34 subjects, no convergence
  
  #Top 300
  pstb = merge(merge(merge(merge(merge(merge(merge(clinic.df, ck.df[, c(1, 4:65)], by = "SampleID"), metbcr.df[, 1:725], by = "SampleID"), swathprot.PCR.df[, 1:303],  by = "SampleID"), ns.df[, 1:81], by.x = "SampleID", by.y="HostSampleID"), nsKO.df[, 1:358], by.x = "SampleID", by.y="HostSampleID"), st.df[, 1:97], by.x = "SampleID", by.y="HostSampleID"), stKO.df[, 1:362], by.x = "SampleID", by.y="HostSampleID") #533 obs. 2040 variables
  pstb.sub = pstb[, c(1, (1:2040)[colnames(pstb) %in% rownames(dcomp)[1:300]], 53, 58)]
  pstb.sub[, 89:212] <- apply(pstb.sub[, 89:212], 2, function(x) log10(x))
  pstb.sub[, 220:234] <- apply(pstb.sub[, 220:234], 2, function(x) asin(sqrt(x))) #for 16S
  pstb.sub[, 257:265] <- apply(pstb.sub[, 257:265], 2, function(x) asin(sqrt(x)))
  pstb.sub[, 235:256] <- apply(pstb.sub[, 235:256], 2, function(x) asin(sqrt(x/1000))) #for 16s genes
  pstb.sub[, 266] <- sapply(pstb.sub[, 266], function(x) asin(sqrt(x/1000)))
  mds.top300 = myHBC(pstb.sub[, 2:268], 265, 5, "none") #27 subjects, no convergence
  
  save(mds.ck, mds.clinic, mds.metb, mds.prot, mds.st, mds.stKO, mds.ns, mds.nsKO, mds.top30, mds.top100, mds.top300, 
       file = "Box Sync/WenyuEx/Rcodes/Revised_MDS.RData")
  
#Visualizing
library(RColorBrewer)
library(plyr)
  #For top ones
  mycol = sample(colorRampPalette(brewer.pal(12,"Paired"))(34), 34) #Random order, so different color for neighboring IDs
  names(mycol) = unique(mds.top100$tb.SubjectID) #Name each color so to match for subset
  
  #For top 30
  hulls = ddply(mds.top30, "tb.SubjectID", function(df) df[chull(df$MDS1, df$MDS2), ])
  ggplot(mds.top30, aes(x = MDS1, y = MDS2, color = tb.SubjectID, fill = tb.SubjectID)) + geom_point(size = 1, alpha = 0.8) + 
     geom_polygon(data = hulls, alpha = 0.3) + scale_colour_manual(values=mycol) + scale_fill_manual(values=mycol) +
     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
           panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1),
           axis.text = element_text(size = 10), legend.position="none") + ylim(c(-1.5, 4.3)) + xlim(c(-3.5,2))
  ggsave("Desktop/Revised_Fig2D_PersonalMDS_top30.pdf", width = 4, height = 4)

  #For different omes
  mycol = sample(colorRampPalette(brewer.pal(12,"Paired"))(14), 14) #Random order, so different color for neighboring IDs
  names(mycol) = unique(mds.ns$tb.SubjectID) #Name each color so to match for subset
  
  ggplot(mds.ck, aes(x = MDS1, y = MDS2, color = tb.SubjectID)) + geom_point(size = 3, alpha = 0.6) + scale_colour_manual(values=mycol) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1),
          axis.text = element_text(size = 14, face = "bold"), axis.title = element_blank(), legend.position="none") + xlim(c(-20, 20))
  ggsave("Desktop/Suppl_PersonalMDS_ck.pdf", width = 6, height = 4)


  #Same for all others
  #mds.metb: xlim(c(-500, 500)), mds.stKO: xlim(c(-200,200)), mds.ns: xlim(c(-50,50)), mds.nsKO: xlim(c(-100,100))
  #mds.prot resulted from too little variance, thus weird plotting
