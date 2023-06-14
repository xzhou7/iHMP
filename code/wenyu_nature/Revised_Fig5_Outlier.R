

###=====2018 June 30 Outlier analysis========
#Using ealthy visits, calculate cumulative Z-score for each variable
#Loading dataset
load("Box Sync/PHI_Protected/Shared_WZ_HMP/Revision/Data/Revision_MultiOmes.RData") #Revised version shared on 04/18/18

##Build a function (CZS, cumulative Z-score) to return the cumutative z-score for each SampleID
#Providing ds with 1:m, SampleID, CL4
  myCZS <- function(ds, m, trans){
    ds = ds[ds$CL4 == "Healthy", ]
    if (trans == "none") {ds.scaled = apply(ds[, 1:m], 2, function(x) scale(as.numeric(x), center = TRUE, scale = TRUE))}
    if (trans == "log2") {ds.scaled = apply(ds[, 1:m], 2, function(x) scale(log2(as.numeric(x)), center = TRUE, scale = TRUE))}
    if (trans == "arcsin") {ds.scaled = apply(ds[, 1:m], 2, function(x) scale(asin(sqrt(as.numeric(x))), center = TRUE, scale = TRUE))}
    ds.cum = apply(ds.scaled, 1, function(x) sum(abs(x), na.rm = TRUE))
    #ranked:
    ds.rk = rank(ds.cum, ties.method = "average")
    ds.f = data.frame(SampleID=ds$SampleID, Cumulative=ds.cum, Rank=ds.rk)
    return(ds.f)
  }

  csz.ck = myCZS(ck.df[, c(4:65, 1, 70, 75)], 62, "none")
  csz.clinic = myCZS(clinic.df[, c(2:54, 1, 58)], 51, "none") 
  csz.metb = myCZS(metbcr.df[, c(2:726,1,731)], 724, "log2") 
  csz.rna = myCZS(rnaseq.log.df[, c(2:10349, 1, 10353)], 10346, "none") 
  csz.prot = myCZS(swathprot.PCR.df[, c(2:303, 1, 304, 309)], 302, "none")
  colnames(st.df)[1] <- "SampleID"
  csz.16s = myCZS(st.df[, c(2:97, 1, 103, 102)], 96, "arcsin") 
  colnames(stKO.df)[1] <- "SampleID"
  csz.16sko = myCZS(stKO.df[, c(2:362, 1, 367,368)], 361, "arcsin") 
  csz.tb = merge(merge(csz.16sko, merge(csz.16s, merge(merge(merge(csz.ck, csz.clinic, by = "SampleID", all = TRUE), csz.metb, by = "SampleID", all = TRUE),
                                                 csz.rna, by = "SampleID", all = TRUE), by = "SampleID", all = TRUE), by = "SampleID", all = TRUE), csz.prot, by = "SampleID", all = TRUE) 
  colnames(csz.tb) <- c("SampleID", "16sko_CumZ", "16sko_Rank", "16s_CumZ", "16s_Rank", "Cytokine_CumZ", "Cytokine_Rank", 
                        "Clinic_CumZ", "Clinic_Rank", "Metb_CumZ", "Metb_Rank", "RNAseq_CumZ", "RNAseq_Rank", "Protein_CumZ", "Protein_Rank")
  write.table(csz.tb, file = "Desktop/Healthy_Outlier_CumZscore_ranked.txt", sep = "\t", row.names = FALSE)

library(ggplot2)
csz.tb$SubjectID = gsub("-\\d+$", "", csz.tb$SampleID) #With minor errors for special SampleID
#For 69-021
ggplot(data = csz.tb, aes(x = log10(Metb_CumZ), y = log10(Cytokine_CumZ))) + geom_point(color = "grey30", alpha = 0.6) + 
  geom_point(data = csz.tb[csz.tb$SubjectID == "69-021", ], aes(x = log10(Metb_CumZ), y = log10(Cytokine_CumZ)), color = "tomato", size = 2) +
  theme_bw() + labs(x = "Cumulative Z-score of Metabolites", y = "Cumulative Z-score of Cytokines") + xlim(c(2.46,2.95)) + ylim(c(1.16, 2.4))
ggsave("Desktop/69021_MetbCyto.pdf", width = 6, height = 4)
