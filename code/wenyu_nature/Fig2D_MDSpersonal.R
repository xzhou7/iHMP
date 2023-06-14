
##====== Plotting MDS graphs to access personal difference==========
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
  
#Build a custom function (HBC, Healthy Baseline Clutering) that prepare data for MDS plot
  #Proving a ds, with 1:m variables, SubjectID and CL4, limiting subjects with more than n baseline visits
  library(vegan)
  library(MASS)
  myHBC <- function(ds, m, n){
    tb = {}
    for (subject in unique(ds$SubjectID)) {
      subset = ds[(ds$CL4 == "Healthy" & ds$SubjectID == subject), ]
      if (nrow(subset) > n) {
        tb = rbind(tb, subset)
      }
    }
    tb.scl = scale(apply(tb[, 1:m], 2, function(x) as.numeric(x)))
    X <- metaMDSdist(tb.scl, distance = "manhattan", autotransform = FALSE, noshare = FALSE, trace = 1, zerodist = "add", na.rm = TRUE)
    mds <- metaMDS(X, k = 2, try = 50, trymax = 50, wascores = FALSE, expand = FALSE, trace = 1, plot = FALSE)
    mdsdata = data.frame(mds$points, tb$SubjectID)
    return(mdsdata)
  } #End of custom function
  
  #Top 50 most personal variables
  dcomp = read.table("Box Sync/PHI_Protected/Shared_WZ_HMP/PublicationDraft/Tables_Final/Fig2A_Table_MultiOmes_VarianceDecomposition.txt", header = TRUE)
  dcomp <- dcomp[order(dcomp$Personal, decreasing = TRUE), ]
  pstb = merge(merge(clinic.df, ck.df[, c(1, 4:65)], by = "SampleID"), metbcr.df[, 1:725], by = "SampleID")
  pstb.sub = pstb[, c(1, (1:844)[colnames(pstb) %in% as.character(dcomp$Variable[1:50])], 53, 58)]
  mds.top = myHBC(pstb.sub[, 2:53], 50, 8) #14 subjects, solution reached
  
  #Host
  mds.ck = myHBC(ck.df[, c(4:65,70, 75)], 62, 9) #11 subjects
  mds.metb = myHBC(metbcr.df[, c(2:726, 731)], 724, 9) #10 subjects
  mds.rna = myHBC(rnaseq.log.df[, c(2:10348, 1, 10353)], 10346, 8) #10 subjects
  mds.clinic = myHBC(clinic.df[, c(2:53, 58)], 51, 9) #12 subjects, solution reached
  mds.prot = myHBC(swathprot.df[, c(2:307,309, 314)], 306, 9) #11 subjects, solution reacher, Warning message:Stress is (nearly) zero - you may have insufficient data
  #Microbiome
  mds.st = myHBC(st.df[, c(2:97, 102:103)], 96, 9) #11 subjects
  mds.ns = myHBC(ns.df[, c(2:81, 86:87)], 80, 8) #14 subjects
  mds.stKO = myHBC(stKO.df[, c(2:362, 367:368)], 361, 9) #11 subjects, solution reached
  mds.nsKO = myHBC(nsKO.df[, c(2:358, 363:364)], 357, 8) #13 subjects, solution reached
  
#Visualizing
  library(RColorBrewer)
  mycol = sample(colorRampPalette(brewer.pal(12,"Paired"))(14), 14) #Random order, so different color for neighboring IDs
  names(mycol) = unique(mds.top$tb.SubjectID) #Name each color so to match for subset
  
  ggplot(mds.top, aes(x = MDS1, y = MDS2)) + geom_point(size = 3, aes(color = tb.SubjectID), alpha = 0.8) + 
    theme(legend.position="none") + scale_colour_manual(values=mycol) 
  ggsave("Desktop/Fig2D_PersonalMDS_top50.pdf", width = 6, height = 4)
  
  ggplot(mds.ck, aes(x = MDS1, y = MDS2)) + geom_point(size = 3, aes(color = tb.SubjectID), alpha = 0.8) + 
    theme(legend.position="none") + scale_colour_manual(values=mycol) +
    xlim(c(-20, 20))
  ggsave("Desktop/Fig2D_suppl_PersonalMDS_ck.pdf", width = 6, height = 4)
  
  ggplot(mds.ns, aes(x = MDS1, y = MDS2)) + geom_point(size = 3, aes(color = tb.SubjectID), alpha = 0.8) + 
    theme(legend.position="none") + scale_colour_manual(values=mycol) 
  ggsave("Desktop/Fig2D_suppl_PersonalMDS_ns.pdf", width = 6, height = 4)
  #Same for all others
  #mds.prot resulted from too little variance, thus weird plotting
  