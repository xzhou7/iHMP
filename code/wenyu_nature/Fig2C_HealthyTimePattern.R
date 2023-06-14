
###===For finding time associating variables within healthy baselines
##Healthy visits, longitudinal time pattern, by correlation of delta change verus delta time 
#Loading dataset
load("Box Sync/PHI_Protected/Shared_WZ_HMP/OmicsData/HMP_MultiOmics.RData")
load("Box Sync/PHI_Protected/Shared_WZ_HMP/OmicsData/HMP_Microbiome_Perc.RData")

#Preparing data
  metbcr.df <- metb.df[, c(1, c(1:1798)[colnames(metb.df) %in% metb.curated$Compounds_ID], 1793:1798)] #Only 724 curated metabolites
  
  StartDate = as.Date("03/01/13", format = "%m/%d/%y")
  #Nasal
  pr = apply(r16s.ns.perc.df[, c(2:954)], 2, function(x) sum(x == 0))
  ns.df <- r16s.ns.perc.df[, c(1, c(2:954)[pr < 400],3807:3812)] #Narrowing to 80 taxa levels that were present in more than half sample size
  ns.df$CollectionDate = StartDate + ns.df$Collection_data
  pr = apply(r16sKO.ns.perc.df[, c(2:6910)], 2, function(x) sum(x < 1))
  nsKO.df <- r16sKO.ns.perc.df[, c(1, c(2:6910)[pr < 400],6911:6916)] #Narrowing to 357 genes that have more than 1% in more than half sample size
  nsKO.df$CollectionDate = StartDate + nsKO.df$Collection_data
  #Stool
  pr = apply(r16s.st.perc.df[, c(2:322)], 2, function(x) sum(x == 0))
  st.df <- r16s.st.perc.df[, c(1, c(2:322)[pr < 430],2276:2281)] #Narrowing to 96 taxa levels
  st.df$CollectionDate = StartDate + st.df$collection.date
  pr = apply(r16sKO.st.perc.df[, c(2:6910)], 2, function(x) sum(x < 1))
  stKO.df <- r16sKO.st.perc.df[, c(1, c(2:6910)[pr < 400],6911:6916)] #Narrowing to 361 genes that have more than 1% in more than half sample size
  stKO.df$CollectionDate = StartDate + stKO.df$collection.date
  
#Build custom fuction (Healthy Ordered, HO) to align time points for individuals having at least THREE longitudinal visits
  #Provding a data frame with SubjectID, CollectionDate, CL4, 1:m variable
  #Return tables with significant(all < p after bonferroni MHC, started with variables with smaller p) variables with the delta values and DiffDay
  library(Hmisc)
  myHO <- function(ds, m, p){
    ho.nord = {}
    #for (subject in unique(ds$SubjectID)) { #for all
    for (subject in SubjectID.d900) { # for SubjectID that have more than 900 days' collection
        hset = ds[(ds$SubjectID == subject & ds$CL4 == "Healthy"), ]
        if (nrow(hset) > 2 ){
          hset$CollectionDate <- as.Date(hset$CollectionDate, "%m/%d/%y")
          hset = hset[order(hset$CollectionDate), ]
          #add delta day
          T1.value = as.numeric(unlist(hset[1, 1:m]))
          T1.Date = hset[, "CollectionDate"][1]
          nor.sbt = as.data.frame(t(apply(hset[, 1:m], 1, function(x) as.numeric(x)-T1.value))) #delta value
          colnames(nor.sbt) = colnames(hset)[1:m]
          nor.sbt$DiffDay = lapply(hset[, "CollectionDate"], function(x) difftime(x, T1.Date, units = "days"))
          nor.sbt$SubjectID = rep(subject, nrow(hset))
          ho.nord = rbind(ho.nord, nor.sbt)
      }
    }
    ho.nord <- ho.nord[ho.nord$DiffDay < 1200, ] #Limiting to DiffDay to major subject pool
    rvalue = apply(ho.nord[, 1:m], 2, function(x) rcorr(x, as.numeric(unlist(ho.nord[, "DiffDay"])), type = "pearson")$r[1,2])
    p.value = apply(ho.nord[, 1:m], 2, function(x) rcorr(x, as.numeric(unlist(ho.nord[, "DiffDay"])), type = "pearson")$P[1,2])
    p.adjvalue = p.adjust(p.value, method = "fdr")
    ho.ls = data.frame(Variable=colnames(ho.nord)[1:m], rho=rvalue, p.adj=p.adjvalue)
    ho.ls <- ho.ls[order(ho.ls$p.adj),]
    sigvar = ho.ls[ho.ls$p.adj < p & !is.na(ho.ls$p.adj), "Variable"] # leaving out variables that didn't pass the significance cut-off
    #return(ho.nord[, c(as.character(sigvar), "DiffDay", "SubjectID")])
    #return(ho.ls[ho.ls$p.adj < p & !is.na(ho.ls$p.adj), ]) #Or returning rho and p.adj tables
    return(ho.nord) #Or return the ho.nord table regardless significance
  } #End of custom function
  
#Provding a data frame with SubjectID, CollectionDate, CL4, and variable in column 1:m 
  #Host
  ho.ck = myHO(ck.df[, c(4:65, 70, 71, 75)], 62, 0.2)  #none
  ho.clinic = myHO(clinic.df[, c(2:54, 58)], 51, 0.05) #16 hits for 5% FDR
  ho.metb = myHO(metbcr.df[, c(2:727, 731)], 724, 0.05) #62 hits for 5% FDR
  ho.prot = myHO(swathprot.df[, c(2:307, 309, 310, 314)], 306, 0.05) #1 hit for 5% FDR
  ho.rna = myHO(rnaseq.log.df[, c(2:10349, 10353)], 10346) #none
  #Microbiome
  ho.st = myHO(st.df[, c(2:97, 102:104)], 96, 0.05) #36 hits
  ho.ns = myHO(ns.df[, c(2:81, 86:88)], 80, 0.05) #none
  ho.stKO = myHO(stKO.df[, c(2:362, 367:369)], 361, 0.05) #11 hits
  ho.nsKO = myHO(nsKO.df[, c(2:358, 363:365)], 357, 0.05) #48 hits
  #When returning rho and p.adj tables, 174 variables:
  write.table(rbind(ho.ck, ho.clinic, ho.metb, ho.prot, ho.st, ho.ns, ho.stKO, ho.nsKO), file = "Desktop/Fig2C_TimeAssociatingVariables.txt", sep = "\t", row.names = FALSE)
  
##Subjects with lengthy collection > 900 in all dataset
  SubjectID.d900 = c("69-001","69-010","69-012","69-022","69-023","69-026","69-027","69-032","69-033","69-034",
                     "69-035","69-036","69-037","69-039","69-045","69-047","69-048","70-1001","70-1002","70-1003",
                     "70-1004","70-1005","70-1006","70-1008","70-1012", "70-1014", "70-1010")  #27 subjects verified
  #Provding a data frame with SubjectID, CollectionDate, CL4, and variable in column 1:m 
  #Host
  ho.ck.d900 = myHO(ck.df[, c(4:65, 70, 71, 75)], 62, 0.2)  #none
  ho.clinic.d900 = myHO(clinic.df[, c(2:54, 58)], 51, 0.05) #11 hits
  ho.metb.d900 = myHO(metbcr.df[, c(2:727, 731)], 724, 0.05) #12 hits
  ho.prot.d900 = myHO(swathprot.df[, c(2:307, 309, 310, 314)], 306, 0.05) #1 hit 
  ho.rna.d900 = myHO(rnaseq.log.df[, c(2:10349, 10353)], 10346) #none
  #Microbiome
  ho.st.d900 = myHO(st.df[, c(2:97, 102:104)], 96, 0.05) #13 hits
  ho.ns.d900 = myHO(ns.df[, c(2:81, 86:88)], 80, 0.05) #none
  ho.stKO.d900 = myHO(stKO.df[, c(2:362, 367:369)], 361, 0.05) #1 hits
  ho.nsKO.d900 = myHO(nsKO.df[, c(2:358, 363:365)], 357, 0.05) #56 hits
  #When returning rho and p.adj tables, 94 variables:
  write.table(rbind(ho.ck.d900, ho.clinic.d900, ho.metb.d900, ho.prot.d900, ho.st.d900, ho.ns.d900, ho.stKO.d900, ho.nsKO.d900), file = "Desktop/Fig2C_TimeAssociatingVariables_D900.txt", sep = "\t", row.names = FALSE)
  
#Visualization
  library(RColorBrewer)
  mycol = sample(colorRampPalette(brewer.pal(12,"Set3"))(79), 79) #Random order, so different color for neighboring IDs
  names(mycol) = unique(ho.clinic$SubjectID) #Name each color so to match for subset
  ggplot(data=ho.clinic, aes(x=unlist(DiffDay), y=ALKP, color = SubjectID)) + geom_point() +
    stat_smooth(method = "lm", col = "grey") +
    scale_colour_manual(values = mycol) +
    labs(x = "Days lapsed", y = "Change in ALKP") +
    theme_bw() + theme(axis.line.x = element_line(color="black", size = 0.5),
                       axis.line.y = element_line(color="black", size = 0.5),
                       panel.grid.minor = element_blank(),
                       panel.border = element_blank(), 
                       legend.position="none")
  ggsave("Desktop/Fig1C_DeltaChange_Daylapse_ALKP_all.pdf", width = 6, height = 4)
  
  #SubjectID.d900
  ggplot(data=ho.clinic.d900, aes(x=unlist(DiffDay), y=ALKP, color = SubjectID)) + scale_colour_manual(values = mycol) +
    geom_point() +
    stat_smooth(method = "lm", col = "grey") +
    labs(x = "Days lapsed", y = "Change in ALKP") +
    theme_bw() + theme(axis.line.x = element_line(color="black", size = 0.5),
                       axis.line.y = element_line(color="black", size = 0.5),
                       panel.grid.minor = element_blank(),
                       panel.border = element_blank(), 
                       legend.position="none")
  ggsave("Desktop/Fig1C_DeltaChange_Daylapse_ALKP_D900.pdf", width = 6, height = 4)
  
  #nHILIC_181.0506_6.2, C03672 HMDB00755 Amino Acid Tyrosine Metabolism
  #For metabolite example that loose significance in SubjectID.d900
  ggplot(data=ho.metb, aes(x=unlist(DiffDay), y=nHILIC_181.0506_6.2, color = SubjectID)) + geom_point() +
    stat_smooth(method = "lm", col = "grey") +
    scale_colour_manual(values = mycol) +
    labs(x = "Days lapsed", y = "Change in Hydroxyphenyllactic acid") +
    theme_bw() + theme(axis.line.x = element_line(color="black", size = 0.5),
                       axis.line.y = element_line(color="black", size = 0.5),
                       panel.grid.minor = element_blank(),
                       panel.border = element_blank(), 
                       legend.position="none")
  ggsave("Desktop/Fig1C_DeltaChange_Daylapse_Hydroxyphenyllactic acid_all.pdf", width = 6, height = 4)
  
  ggplot(data=ho.metb.d900, aes(x=unlist(DiffDay), y=nHILIC_181.0506_6.2, color = SubjectID)) + scale_colour_manual(values = mycol) +
    geom_point() +
    labs(x = "Days lapsed", y = "Change in Hydroxyphenyllactic acid") +
    theme_bw() + theme(axis.line.x = element_line(color="black", size = 0.5),
                       axis.line.y = element_line(color="black", size = 0.5),
                       panel.grid.minor = element_blank(),
                       panel.border = element_blank(), 
                       legend.position="none")
  ggsave("Desktop/Fig1C_DeltaChange_Daylapse_Hydroxyphenyllactic acid_D900.pdf", width = 6, height = 4)

