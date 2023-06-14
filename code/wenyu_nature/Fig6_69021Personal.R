
###===Plotting 69-021, or ZNED4XZ, personal characters as compared to the rest of the cohort and over the time
##Using healthy visits, calculate cumulative Z-score for each variable
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

##Build a function (CZS, cumulative Z-score) to return the cumutative z-score for each SampleID
  #Providing ds with 1:m, SampleID, CL4
  myCZS <- function(ds, m){
    ds = ds[ds$CL4 == "Healthy", ]
    ds.scaled = apply(ds[, 1:m], 2, function(x) scale(as.numeric(x), center = TRUE, scale = TRUE)) #Calculating Z-score for each Sample
    ds.cum = apply(ds.scaled, 1, function(x) sum(abs(x), na.rm = TRUE))
    ds.f = data.frame(SampleID=ds$SampleID, cumulative=ds.cum)
    ##If ranked:
    #ds.rk = rank(ds.cum, ties.method = "average")
    #ds.f = data.frame(SampleID=ds$SampleID, cumulative=ds.rk)
    return(ds.f)
  } #End of custom function
  
  #Host
  csz.ck = myCZS(ck.df[, c(4:65, 1, 70, 75)], 62)
  csz.clinic = myCZS(clinic.df[, c(2:54, 1, 58)], 51) 
  csz.metb = myCZS(metbcr.df[, c(2:726, 1, 731)], 724) 
  csz.rna = myCZS(rnaseq.log.df[, c(2:10349, 1, 10353)], 10346) 
  csz.prot = myCZS(swathprot.df[, c(2:307,1, 314)], 306) 
  
  #Microbiome
  colnames(st.df)[1] <- "SampleID" 
  colnames(ns.df)[1] <- "SampleID" 
  colnames(stKO.df)[1] <- "SampleID" 
  colnames(nsKO.df)[1] <- "SampleID"
  csz.st = myCZS(st.df[, c(2:97, 1, 102)], 96) 
  csz.ns = myCZS(ns.df[, c(2:81, 1, 86)], 80) 
  csz.stKO = myCZS(stKO.df[, c(2:362, 1, 367)], 361) 
  csz.nsKO = myCZS(nsKO.df[, c(2:358, 1, 363)], 357) 
  #Merge all
  csz.tb = merge(merge(merge(merge(csz.ck, merge(csz.clinic, merge(merge(merge(csz.metb, csz.rna, by = "SampleID", all = TRUE), csz.prot, by = "SampleID", all = TRUE),
                csz.st, by = "SampleID", all = TRUE), by = "SampleID", all = TRUE), by = "SampleID", all = TRUE), 
                csz.stKO, by = "SampleID", all = TRUE), csz.ns, by = "SampleID", all = TRUE), csz.nsKO, by = "SampleID", all = TRUE)
                                          
  colnames(csz.tb) <- c("SampleID", "Cytokines", "ClinicLabs", "Metabolites", "Transcripts", "Proteins", "Gut_Microbes", "Gut_MicrobialGenes", "Nasal_Microbes", "Nasal_MicrobialGenes")
  write.table(csz.tb, file = "Desktop/Fig6_Healthy_Outlier_CumZscore.txt", sep = "\t", row.names = FALSE)

#Visualizing
  library(ggplot2)
  csz.tb$SubjectID = gsub("-\\d+$", "", csz.tb$SampleID) #With minor errors for special SampleID
  #For 69-021
  ggplot(data = csz.tb, aes(x = log10(Metabolites), y = log10(Cytokines))) + geom_point(color = "grey60", alpha = 0.5) + 
    geom_point(data = csz.tb[csz.tb$SubjectID == "69-021", ], aes(x = log10(Metabolites), y = log10(Cytokines)), color = "tomato", size = 2) +
    geom_text(data = csz.tb[csz.tb$SubjectID == "69-021", ], aes(label = gsub("69-021-2011", "H6", gsub("69-021-4011", "H5", gsub("69-021-0", "H", SampleID)))), parse = TRUE, hjust = 1, size = 2) +
    theme_bw() + labs(x = "Cumulative Z-score of Metabolites", y = "Cumulative Z-score of Cytokines")
  ggsave("Desktop/Fig6_69021_MetbCyto.pdf", width = 6, height = 4)
  
##Plotting immune panel for 69-021 longitudinal collection
#Importing immune genes
  InputFile = "Box Sync/WenyuEx/Rcodes/HMP_ImmunePanels - Final.tsv"
  Imp = unlist(read.table(InputFile, header = TRUE, sep="\t",quote="",stringsAsFactors=FALSE))
  Imlist = c(Imp, "IL1RA", "TNFA", "IL1A", "TNFB") #Adding four more from cytokine set
  
  #Preparing data with cytokines and HSCRP in ClinicLabs
  ds = ck.df[ck.df$SubjectID == "69-021", c("SampleID", intersect(colnames(ck.df), Imlist))]
  ds <- merge(ds, clinic.df[clinic.df$SubjectID == "69-021", c("SampleID", "HSCRP", "CollectionDate", "CL4")], by = "SampleID")
  
  ds$CollectionDate = as.Date(ds$CollectionDate, "%m/%d/%y")
  ds <- ds[order(ds$CollectionDate), ]
  
  #Heatmap
  library("gplots")
  library("RColorBrewer")
  
  pd = scale(apply(ds[, 2:13], 2, function(x) as.numeric(as.character(x))))
  
  mycol = rev(colorRampPalette(brewer.pal(11, "RdYlBu"))(80))
  mybk = seq(-0.5,1.6,length=81)
  
  heatmap.2(t(pd), Rowv = TRUE, Colv = FALSE, key = TRUE,
            symm = FALSE,symkey = FALSE,
            scale = "none", trace = "none", dendrogram = "row",  
            distfun = function(x) dist(x,method = "euclidean"), hclustfun = function(x) hclust(x,method = "median"),
            col = mycol)
  
##Finding correlation with IL1RA and HSCRP with only 69-021 data
#Build a function (PPC, personal pearson correlation) to return a table with variables correlated with IL1RA or HSCRP in 69-021's data
  #Providing ds with 1:m, SampleID and SubjectID
  myPPC <- function(ds, m, target, fdr){
    pds = ds[ds$SubjectID == "69-021", ]
    if (target == "IL1RA"){
      tvalue = ck.df[ck.df$SubjectID == "69-021", c("SampleID", target)]
    }
    if (target == "HSCRP"){
      tvalue = clinic.df[clinic.df$SubjectID == "69-021", c("SampleID", target)]
    }
    #Pearson correlations
    p.val = r = {}
    pds <- merge(pds, tvalue, by = "SampleID", all = TRUE)
    n = m+1
    for (i in 2:n) { 
      p.val[colnames(pds)[i]] = rcorr(as.numeric(pds[, i]), as.numeric(pds[, target]), type = "pearson")$P[1,2]
      r[colnames(pds)[i]] = rcorr(as.numeric(pds[, i]), as.numeric(pds[, target]), type = "pearson")$r[1,2]
    }
    p.ad <- p.adjust(p.val, method = "fdr")
    sig.var = names(p.ad)[!is.na(p.ad) & p.ad < fdr]   # Using FDR to filter
    pcor.sig = data.frame(SampleID = pds$SampleID, pds[, sig.var], target=pds[, target])
    #Below are adding p-value and rho for each sig.var
    pcor.sig$SampleID <- as.character(pcor.sig$SampleID)
    pcor.sig$target <- as.character(pcor.sig$target)
    pcor.sig[nrow(pcor.sig)+1, ] <- c("p-adj", p.ad[sig.var], target) #p.ad values
    pcor.sig[nrow(pcor.sig)+1, ] <- c("r", r[sig.var], target) #rho values
    return(pcor.sig)
  } #End of custom function
  
  #Host
  ppc.ck = myPPC(ck.df[, c(4:25, 27:65, 1, 70)], 61, "IL1RA", 0.2) # minus "IL1RA"(col 26), none
  ppc.clinic = myPPC(clinic.df[, c(2:52, 1, 53)], 51, "IL1RA", 0.2) # 1 hit at 0.1, 2 hit at 0.2, two have too much missing data
  ppc.metb = myPPC(metbcr.df[, c(2:725, 1, 726)], 724, "IL1RA", 0.2) # 10 hits at 0.1 fdr, 14, at 0.2 fdr
  ppc.rna = myPPC(rnaseq.log.df[, c(2:10347, 1, 10348)], 10346, "IL1RA", 0.2) 
  rnazero <- ppc.rna[6, ] == 0 #Filtering genes that are zero for "69-021-09"
  ppc.rna <- ppc.rna[, (!rnazero)] #236 hits at 0.05 fdr, 704 at 0.2 fdr
  ppc.prot = myPPC(swathprot.df[, c(2:307,1, 309)], 306, "IL1RA", 0.2) #15 hits at 0.1 fdr, 20 at 02 fdr
  
  #Microbiome #none
  colnames(st.df)[1] <- "SampleID" 
  colnames(ns.df)[1] <- "SampleID" 
  colnames(stKO.df)[1] <- "SampleID" 
  colnames(nsKO.df)[1] <- "SampleID"
  ppc.st = myPPC(st.df[, c(2:97, 1, 103)], 96, "IL1RA", 0.2) 
  ppc.ns = myPPC(ns.df[, c(2:81, 1, 87)], 80, "IL1RA", 0.2) 
  ppc.stKO = myPPC(stKO.df[, c(2:362, 1, 368)], 361, "IL1RA", 0.2) 
  ppc.nsKO = myPPC(nsKO.df[, c(2:358, 1, 364)], 357, "IL1RA", 0.2) 
  #Merge all
  ppc.IL1RA = cbind(ppc.clinic[, c(1, 3, 5)], ppc.metb[, 2:11], ppc.rna[, 2:237], ppc.prot[, 2:16])
  #ppc.IL1RA.fdr20 = cbind(ppc.clinic[, c(1, 3, 4, 6)], ppc.metb[, 2:15], ppc.rna[, 2:703], ppc.prot[, 2:21])
  write.table(t(ppc.IL1RA), file = "Desktop/Fig6_IL1RAcorrelation.txt", row.names = TRUE, sep = "\t")
  
  #Visualize, the graph is not interesting to present as a separate figure
  InputFile = "Box Sync/WenyuEx/Research_HMP/Results_Tables/Tables_Final/Fig6_IL1RAcorrelation.txt"
  tb.cor = read.table(InputFile, header = TRUE, sep="\t",quote="",stringsAsFactors=FALSE)
  
  pd = tb.cor[c(2:264), 2:15]
  colnames(pd) = tb.cor[1, 2:15]
  pd.ord <- pd[, ds$SampleID] #Order by CollectionDate
  rownames(pd.ord) = tb.cor[c(2:264), 1]
  rownames(pd.ord)[which(tb.cor$X.1 != "")-1] <- tb.cor$X.1[which(tb.cor$X.1 != "")] #Replacing metb and prot's annotation
  rownames(pd.ord)[2] <- "IL1RA"
  
  pd.sc = t(scale(t(apply(pd.ord, 2, function(x) as.numeric(x)))))
  rownames(pd.sc) = rownames(pd.ord)
  
  mycol = rev(colorRampPalette(brewer.pal(11, "RdYlBu"))(80))
  
  heatmap.2(pd.sc, Rowv = TRUE, Colv = FALSE, key = TRUE,
            symm = FALSE,symkey = FALSE,
            scale = "none", trace = "none", dendrogram = "row",  
            distfun = function(x) dist(x,method = "euclidean"), hclustfun = function(x) hclust(x,method = "median"),
            col = mycol)
  
#Correlation among variables
  colnames(ppc.IL1RA)[3] <- "IL1RA"
  cor.all = rcorr(apply(ppc.IL1RA[, 2:264], 2, function(x) as.numeric(x)), type = "pearson")
  #Turning into a data frame
  ind <- which(upper.tri(cor.all$r,diag=F) , arr.ind = TRUE)
  cor.all.tb = data.frame(col = dimnames(cor.all$r)[[2]][ind[,2]] ,
                          row = dimnames(cor.all$r)[[1]][ind[,1]] ,
                          r = cor.all$r[ ind ], pval = cor.all$P[ ind ])
  cor.all.tb <- cor.all.tb[order(cor.all.tb$r, decreasing = TRUE), ]
  write.table(cor.all.tb, file = "Desktop/Fig6_IL1RAcorrelation_Cor.txt", row.names = FALSE, sep = "\t")
  
#For "HSCRP"
  #Host
  ppc.ck = myPPC(ck.df[, c(4:65, 1, 70)], 61, "HSCRP", 0.2) # 38 hits at 0.2 fdr, 0 at 0.1 fdr
  ppc.clinic = myPPC(clinic.df[, c(2:25, 27:52, 1, 53)], 51, "HSCRP", 0.2) # minus "HSCRP"(col 26), 5 hit, there were seven, but two were due to missing data
  ppc.metb = myPPC(metbcr.df[, c(2:725, 1, 726)], 724, "HSCRP", 0.2) #10 hits
  ppc.rna = myPPC(rnaseq.log.df[, c(2:10347, 1, 10348)], 10346, "HSCRP", 0.2) 
  rnazero <- ppc.rna[6, ] == 0 #Filtering genes that are zero for "69-021-09"
  ppc.rna <- ppc.rna[, (!rnazero)] #20 hits
  ppc.prot = myPPC(swathprot.df[, c(2:307,1, 309)], 306, "HSCRP", 0.2) #2 hits
  
  #Microbiome #none
  colnames(st.df)[1] <- "SampleID" 
  colnames(ns.df)[1] <- "SampleID" 
  colnames(stKO.df)[1] <- "SampleID" 
  colnames(nsKO.df)[1] <- "SampleID"
  ppc.st = myPPC(st.df[, c(2:97, 1, 103)], 96, "HSCRP", 0.2) #1 hit, but only two data points, therefore remove
  ppc.ns = myPPC(ns.df[, c(2:81, 1, 87)], 80, "HSCRP", 0.2) 
  ppc.stKO = myPPC(stKO.df[, c(2:362, 1, 368)], 361, "HSCRP", 0.2) 
  ppc.nsKO = myPPC(nsKO.df[, c(2:358, 1, 364)], 357, "HSCRP", 0.2) 
  #Merge all
  ppc.HSCRP = cbind(ppc.ck, ppc.clinic[, c(3:6, 8)], ppc.metb[, 2:11], ppc.rna[, 2:21], ppc.prot[, 2:3])
  write.table(t(ppc.HSCRP), file = "Desktop/Fig6_HSCRPcorrelation.txt", row.names = TRUE, sep = "\t")
  
  