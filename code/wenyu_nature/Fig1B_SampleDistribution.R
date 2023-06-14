

##===For plotting Sample collection distributions in Fig 1B
  #First loading data
  load("Box Sync/PHI_Protected/Shared_WZ_HMP/OmicsData/HMP_MultiOmics.RData")

  #Clean up rows that should be left out, and add subject classifications to sample list
  ls.clean <- ls[ls$CollectionDate != "" & ls$CL4 != "Fiber" & ls$CL4 != "", ] #1046 obs.
  ls.clean <- merge(ls.clean, sc[, c(1,2,6)], by = "SubjectID") #Adding IR/IS, and Subject ClassA #1044 obs.

  library(ggplot2)
  #Plotting subject group distribution (Fig1B left)
  sc$ClassA <- factor(sc$ClassA, levels = c("Diabetic", "Crossover", "Prediabetic", "Control"))
  mycolor <- c("IR" = "#EF6125", "IS" = "#90AA3C", "Unknown" = "grey70")
  ggplot(data=sc, aes(x = ClassA)) + geom_bar(aes(fill = IRIS)) + 
    theme(plot.background = element_blank(),panel.background = element_blank(),
          axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
          axis.line.x = element_line(color="black", size = 0.5),
          axis.line.y = element_line(color="black", size = 0.5),panel.border = element_blank()) +
    scale_fill_manual(values = mycolor) + labs(x = "", y = "Number of Subjects") + geom_text(stat="count", aes(label=..count..),vjust=-1)
  
  ggsave("Desktop/Fig1B_SubjectGroup.pdf", width = 3, height = 4)

  #Plotting visit distribution (Fig1B right)
  ls.clean$CLM <- sub("Allergy|Ant_L/Infection_L|Colonoscopy|Colonoscopy_L|Post-Travel|Post-Travel_L|Special|Stress|Stress_L|Surgery", "Others", ls.clean$CL4)  #Lumping together all minor visits types into "Others" category
  ls.clean$CLM <- sub("Infection|Infection_L", "Infection", ls.clean$CLM) #Lumping Infection together
  ls.clean$CLM <- sub("Imz|Imz_L", "Immun", ls.clean$CLM) #Lumping Imz together
  ls.clean$CLM <- sub("Ant|Ant_L", "Abx", ls.clean$CLM) #Lumping Antibotics together
  ls.clean$CLM <- sub("Weight-gain|Weight-loss", "Overfeed", ls.clean$CLM) #Lumping overfeeding together
  ls.clean$CLM <- factor(ls.clean$CLM, levels = c("Healthy", "Infection", "Immun", "Overfeed", "Abx", "Others"))
  ggplot(data=ls.clean, aes(x = CLM)) + geom_bar(aes(fill = IRIS)) +
    theme(plot.background = element_blank(),panel.background = element_blank(),
          axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
          axis.line.x = element_line(color="black", size = 0.5),
          axis.line.y = element_line(color="black", size = 0.5),panel.border = element_blank()) +
    scale_fill_manual(values = mycolor) + labs(x = "", y = "Number of Visits") +
    geom_text(stat="count", aes(label=..count..),vjust=-1)
  
  ggsave("Desktop/Fig1B_VisitGroup.pdf", width = 4, height = 6)
