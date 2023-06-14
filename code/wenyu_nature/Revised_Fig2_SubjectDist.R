

##====== Calculating distance between subjects==========
#Loading dataset
load("Box Sync/PHI_Protected/Shared_WZ_HMP/Revision/Data/Revision_MultiOmes.RData") #Revised version shared on 04/18/18

 #Creating merged dataset
 colnames(ck.df)[4:65] <- paste("Cytokines", colnames(ck.df[, 4:65]), sep = "_")
 colnames(clinic.df)[2:52] <- paste("ClinicLabs", colnames(clinic.df[, 2:52]), sep = "_")
 colnames(metbcr.df)[2:725] <- paste("Metabolites", colnames(metbcr.df[, 2:725]), sep = "_")
 colnames(rnaseq.log.df)[2:10347] <- paste("Transcripts", colnames(rnaseq.log.df[, 2:10347]), sep = "_")
 colnames(swathprot.PCR.df)[2:303] <- paste("Proteins", colnames(swathprot.PCR.df[, 2:303]), sep = "_")
 
 colnames(st.df)[2:97] <- paste("Gut_Microbes", colnames(st.df[, 2:97]), sep = "_")
 colnames(stKO.df)[2:362] <- paste("Gut_MicrobialGenes", colnames(stKO.df[, 2:362]), sep = "_")
 colnames(ns.df)[2:81] <- paste("Nasal_Microbes", colnames(ns.df[, 2:81]), sep = "_")
 colnames(nsKO.df)[2:358] <- paste("Nasal_MicrobialGenes", colnames(nsKO.df[, 2:358]), sep = "_")
 
 merged.df = merge(merge(merge(merge(merge(ls, ck.df[, c(1, 4:65)], by = "SampleID", all = TRUE), clinic.df[, 1:52],  by = "SampleID", all = TRUE), 
                   metbcr.df[, 1:725], by = "SampleID", all = TRUE), swathprot.PCR.df[, 1:303],  by = "SampleID", all = TRUE), 
                   rnaseq.log.df[, 1:10347], by = "SampleID", all = TRUE)
 merged.mb.df = merge(merge(merge(merge(st.df[, 1:97], merged.df, by.x = "HostSampleID", by.y = "SampleID", all.y = TRUE), stKO.df[, 1:362], by = "HostSampleID", all.x = TRUE),
                      ns.df[, 1:81], by = "HostSampleID", all.x = TRUE), nsKO.df[, 1:358], by = "HostSampleID", all.x = TRUE)
 #2:97: st, 98:103:ls, 104:165:ck, 166:216:Clinic, 217:940:metb, 941:1242:prot, 1243:11588:rnas
 #11589:11949:stKO, 11950:12029:ns, 12030:12386:nsKO
 
 #Linearly transform
 #metb
 merged.mb.df[, 217:940] <- apply(merged.mb.df[, 217:940], 2, function(x) log10(x))
 #st and ns
 merged.mb.df[, 2:97] <- apply(merged.mb.df[, 2:97], 2, function(x) asin(sqrt(x)))
 merged.mb.df[, 11950:12029] <- apply(merged.mb.df[, 11950:12029], 2, function(x) asin(sqrt(x)))
 #stKO and nsKO
 merged.mb.df[, 11589:11949] <- apply(merged.mb.df[, 11589:11949], 2, function(x) asin(sqrt(x/1000)))
 merged.mb.df[, 12030:12386] <- apply(merged.mb.df[, 12030:12386], 2, function(x) asin(sqrt(x/1000)))
 
#Custom function to calculate subject distance (SDist) using median values, given top n variables
 mySDist <- function(n) {
  tgt = dcomp$Var[1:n]
 #mySDist <- function(assay){
   #tgt = dcomp[dcomp$dt == assay, "Var"]
   #n = length(tgt)
   ds = merged.mb.df[merged.mb.df$CL4 == "Healthy", c(1, 98:103, (1:12386)[colnames(merged.mb.df) %in% tgt])]
   ds[, 8:dim(ds)[2]] <- apply(ds[, 8:dim(ds)[2]], 2, function(x) scale(as.numeric(x)))  #scale to sd 1
   ds.median = {}
   for (subject in unique(ds$SubjectID)){
     hset = ds[ds$SubjectID == subject, ]
     if (nrow(hset) > 2) {
       hset.median = apply(hset[, 8:dim(ds)[2]], 2, function(x) median(as.numeric(x), na.rm = TRUE))
       ds.median <- rbind(ds.median, c(subject, hset.median))
     }
   }
   rownames(ds.median) <- ds.median[, 1]
   dm = as.matrix(dist(apply(ds.median[, 2:dim(ds.median)[2]], 2, function(x) as.numeric(as.character(x)))))
   dm.median = median(dm[lower.tri(dm)], na.rm = TRUE)/n
   dm.ave = mean(dm[lower.tri(dm)], na.rm = TRUE)/n
   return(c(n, dm.median, dm.ave, 
            sum(grepl("Cytokines", tgt))/n, sum(grepl("ClinicLabs", tgt))/n, sum(grepl("Metabolites", tgt))/n, sum(grepl("Transcripts", tgt))/n, sum(grepl("Proteins", tgt))/n,
            sum(grepl("Gut_Microbes", tgt))/n, sum(grepl("Gut_MicrobialGenes", tgt))/n, sum(grepl("Nasal_Microbes", tgt))/n, sum(grepl("Nasal_MicrobialGenes", tgt))/n))
   #return(c(assay, n, dm.median, dm.ave))
   }
 
 dcomp = vd.comb[order(vd.comb$ICC, decreasing = TRUE), ]
 dcomp$Var = paste(dcomp$dt, rownames(dcomp), sep = "_")
 
 sdist.top = t(sapply(seq(30,3000, by =30), function(x) mySDist(x)))
 sdist.assay = t(sapply(c("Cytokines", "ClinicLabs", "Metabolites", "Proteins", "Gut_Microbes", "Gut_MicrobialGenes", "Nasal_Microbes", "Nasal_MicrobialGenes", "Transcripts"), function(x) mySDist(x)))
 
 write.table(sdist.top, file = "Desktop/Revised_Fig2_SubjectDistance_tops.txt", row.names = FALSE, sep = "\t")
 write.table(sdist.assay, file = "Desktop/Revised_Fig2_SubjectDistance_byomes.txt", row.names = FALSE, sep = "\t")
 