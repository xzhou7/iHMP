
###===For microbiome data: calculate Shannon diversity and Chao richness
load("Box Sync/PHI_Protected/Shared_WZ_HMP/OmicsData/HMP_Microbiome_Perc.RData")

  library(vegan)
  StartDate = as.Date("03/01/13", format = "%m/%d/%y")
  #For st
  st.dv.df = data.frame(HostSampleID=r16s.st.perc.df[, 1],
                        phylum_chao=specpool(r16s.st.perc.df[, 2:15], r16s.st.perc.df[, 1])$chao, 
                        class_chao=specpool(r16s.st.perc.df[, 16:41], r16s.st.perc.df[, 1])$chao,
                        order_chao=specpool(r16s.st.perc.df[, 42:80], r16s.st.perc.df[, 1])$chao,
                        family_chao=specpool(r16s.st.perc.df[, 81:155], r16s.st.perc.df[, 1])$chao,
                        genus_chao=specpool(r16s.st.perc.df[, 156:322], r16s.st.perc.df[, 1])$chao, 
                        phylum_shannon=diversity(r16s.st.perc.df[, 2:15]), 
                        class_shannon=diversity(r16s.st.perc.df[, 16:41]),
                        order_shannon=diversity(r16s.st.perc.df[, 42:80]),
                        family_shannon=diversity(r16s.st.perc.df[, 81:155]),
                        genus_shannon=diversity(r16s.st.perc.df[, 156:322]), 
                        r16s.st.perc.df[, 2276:2281])
  st.dv.df$CollectionDate = StartDate + st.dv.df$collection.date
  stKO.dv.df = data.frame(HostSampleID=r16sKO.st.perc.df[, 1],
                          KO_chao=specpool(r16sKO.st.perc.df[, 2:6910], r16sKO.st.perc.df[, 1])$chao,
                          KO_shannon=diversity(r16sKO.st.perc.df[, 2:6910]))
  stKO.dv.df$CollectionDate = StartDate + stKO.dv.df$collection.date
  #For ns
  ns.dv.df = data.frame(HostSampleID=r16s.ns.perc.df[, 1],
                        phylum_chao=specpool(r16s.ns.perc.df[, 2:25], r16s.ns.perc.df[, 1])$chao, 
                        class_chao=specpool(r16s.ns.perc.df[, 26:78], r16s.ns.perc.df[, 1])$chao,
                        order_chao=specpool(r16s.ns.perc.df[, 79:171], r16s.ns.perc.df[, 1])$chao,
                        family_chao=specpool(r16s.ns.perc.df[, 172:378], r16s.ns.perc.df[, 1])$chao,
                        genus_chao=specpool(r16s.ns.perc.df[, 379:954], r16s.ns.perc.df[, 1])$chao, 
                        phylum_shannon=diversity(r16s.ns.perc.df[, 2:25]), 
                        class_shannon=diversity(r16s.ns.perc.df[, 26:78]),
                        order_shannon=diversity(r16s.ns.perc.df[, 79:171]),
                        family_shannon=diversity(r16s.ns.perc.df[, 172:378]),
                        genus_shannon=diversity(r16s.ns.perc.df[, 379:954]), 
                        r16s.ns.perc.df[, 3807:3812])
  ns.dv.df$CollectionDate = StartDate + ns.dv.df$Collection_data
  
  save(st.dv.df, ns.dv.df, file="Box Sync/WenyuEx/Research_HMP/OmicsData/HMP_Microbiome_Diversity.RData")
   
###===For revision: list sample collections
#Loading dataset
load("Box Sync/PHI_Protected/Shared_WZ_HMP/Revision/Data/Revision_MultiOmes.RData") #Revised version shared on 04/18/18

#Generating a list with available omics mapped, 1091 visits with 923 had at least six profiles
  ls$Cytokines = ifelse(ls$SampleID %in% ck.df$SampleID, 1, 0)
  ls$ClinicLabs = ifelse(ls$SampleID %in% clinic.df$SampleID, 1, 0)
  ls$Metabolites = ifelse(ls$SampleID %in% metbcr.df$SampleID, 1, 0)
  ls$Proteins = ifelse(ls$SampleID %in% swathprot.PCR.df$SampleID, 1, 0)
  ls$Transcripts = ifelse(ls$SampleID %in% rnaseq.log.df$SampleID, 1, 0)
  ls$Gut_16S = ifelse(ls$SampleID %in% st.df$HostSampleID, 1, 0)
  ls$Nasal_16S = ifelse(ls$SampleID %in% ns.df$HostSampleID, 1, 0)
  ls$Num_Type = rowSums(ls[, 8:14]) 
  
  write.table(ls, file = "Desktop/SampleCollectionProfileList.txt", sep = "\t", row.names = FALSE)
  write.table(sc, file = "Desktop/Revised_SubjectList.txt", sep = "\t", row.names = FALSE)
  
#Generating timespan and # of visits per subject
  sc.tb = {}
  for (subject in sc$SubjectID){
    sset = ls[ls$SubjectID == subject & ls$CollectionDate != "", ]
    sset$CollectionDate <- as.Date(sset$CollectionDate, "%m/%d/%y")
    sset = sset[order(sset$CollectionDate), ]
    span = difftime(sset$CollectionDate[nrow(sset)], sset$CollectionDate[1], units = "days")
    nvisit = nrow(sset)
    hset = sset[sset$CL4 == "Healthy", ]
    h_visit = nrow(hset)
    sc.tb=rbind(sc.tb, c(subject, nvisit, span, h_visit))
  }
  write.table(sc.tb, file = "Desktop/Revised_Subject_n_all_visit.txt", sep = "\t", row.names = FALSE)
  
  