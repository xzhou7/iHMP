#### Remove All the Variables from the Workspace ####
rm(list = ls())

#### Setting Directory ####
setwd('/Users/hmallick/Dropbox (Huttenhower Lab)/HMP2/analysis/differential_abundance')

#### Load libraries ####
library(data.table)
library(stringr)
library(ggplot2)
library(gridExtra)
library(tibble)
library(scales)

# Parameters
FDR_Threshold<-0.05

# Define all  transformations
transMethods<-c('CLR_perfeature', 'CLR_persample', 'CLR_overall', 'Vanilla')

# Load utility function
findFiles<-function(x){
    subdir1<-sprintf('./%s/diagnosis/', x)
    subdir2<-sprintf('./%s/active_vs_inactive/', x)
    foo <-list.files(subdir1, pattern="*.tsv")
    files1<-paste(subdir1, foo, sep='')
    bar <-list.files(subdir2, pattern="*.tsv")
    files2<-paste(subdir2, bar, sep='')
    files<-c(files1, files2)
    return(files)
}


# Extract file names
file_names<-list()
for (i in 1:length(transMethods)){
    file_names[[i]]<-findFiles(transMethods[i])
}

# Define main looping variables
names<-c('bugs',
         'metabolites',
         'pwyDNAs',
         'pwyRNAs',
         'ecDNAs',
         'ecRNAs',
         'koDNAs',
         'koRNAs',
         'proteins',
         'proteinECs',
         'proteinKOs')

names_display<-c('Species',
         'Metabolites',
         'Pwy DNAs',
         'Pwy RNAs',
         'EC DNAs',
         'EC RNAs',
         'KO DNAs',
         'KO RNAs',
         'Proteins',
         'Protein ECs',
         'Protein KOs')

strings<-c('diagnosis/IBD_relations','Active_IBD_relations')

strings2<-c('diagnosis','active_vs_inactive')

strings_display<-c('Diagnosis', 'Dysbiosis')

axes_display<- list(c("Significant", "Not Significant"), c("Significant", "Not Significant"))

# Delete metabolites_Sena
file_names[[4]]<-file_names[[4]][!str_detect(file_names[[4]], 'metabolites_Sena')]


########################################
# Contingency Tables and Scatter Plots #
########################################

mcnemarTest<-matrix(nrow=length(names), ncol=length(strings))
fisherTest<-matrix(nrow=length(names), ncol=length(strings))
contingency_mcnemar.list<-matrix(list(), length(names),length(strings))
contingency_fisher.list<-matrix(list(), length(names),length(strings))
scatterplot.list<-matrix(list(), length(names),length(strings))

# Loop Over Each Pairwise Combinations
for (i in 1:(length(transMethods)-1)) { # Number of cross-method comparisons
    for (m in (i+1):length(transMethods)) { # Number of cross-method comparisons
        for (j in 1:length(names)) {
            for (k in 1:length(strings)){

            # Extract Relevant File Names
            files_temp_1<-file_names[[i]][str_detect(file_names[[i]], paste(strings[k], names[j], sep='_'))]
            files_temp_2<-file_names[[m]][str_detect(file_names[[m]], paste(strings[k], names[j], sep='_'))]

            # Read the corresponding data frames
            dat_1 <- data.table::fread(files_temp_1, header = TRUE)
            dat_2 <- data.table::fread(files_temp_2, header = TRUE)

            # Subset to ID and Qval
            df_1<-dat_1[, c('#', 'prevalence', 'minQ')]
            colnames(df_1)<-c('ID', 'prevalence', 'qval_1')
            df_2<-dat_2[, c('#', 'prevalence', 'minQ')]
            colnames(df_2)<-c('ID', 'prevalence', 'qval_2')

            # Merge
            df<-merge(df_1, df_2, c('ID', 'prevalence'))

            # Create Table
            nonVanilla_Significant<-df$ID[df$qval_1<FDR_Threshold]
            nonVanilla_nonSignificant<-df$ID[df$qval_1>=FDR_Threshold]
            Vanilla_Significant<-df$ID[df$qval_2<FDR_Threshold]
            Vanilla_nonSignificant<-df$ID[df$qval_2>=FDR_Threshold]


            # Contingency
            A<-length(intersect(Vanilla_Significant, nonVanilla_Significant))
            B<-length(intersect(Vanilla_Significant, nonVanilla_nonSignificant))
            C<-length(intersect(Vanilla_nonSignificant, nonVanilla_Significant))
            D<-length(intersect(Vanilla_nonSignificant, nonVanilla_nonSignificant))

            # McNemar's Exact Test - Calculate P value
            Mat<-matrix(c(A, B, C, D),nrow=2,ncol=2)
            rownames(Mat)<-c('Significant', 'Not Significant')
            colnames(Mat)<-c('Significant', 'Not Significant')
            mcnemarTest[j, k]<-ifelse(is.finite(round(mcnemar.test(Mat)$p.value,4)), round(mcnemar.test(Mat)$p.value,4), 1)
            fisherTest[j, k]<-round(fisher.test(Mat)$estimate,2)


            # Create Data for Plot
            Mat0<-as.data.frame(Mat)
            Mat0<-rownames_to_column(Mat0, 'ID')
            Mat_Final<-melt(Mat0)
            Mat_Final$variable<-as.factor(Mat_Final$variable)
            Mat_Final$ID<-as.factor(Mat_Final$ID)
            Mat_Final$percent <- Mat_Final$value / sum(Mat_Final$value)
            Mat_Final$label<-ifelse(Mat_Final$percent<1e-04, paste0(comma(Mat_Final$value), ' (', round(100 * Mat_Final$percent, 4), '%)'), paste0(comma(Mat_Final$value), ' (', round(100 * Mat_Final$percent, 2), '%)'))

            # Contingency Plot (McNemar)
            contingency_mcnemar.list[[j,k]]<-ggplot(Mat_Final, aes(variable, ID, fill=percent, label=label)) +
                geom_tile(colour = "black") +
                geom_text(aes(label = label), colour = "white") +
                scale_fill_gradient(low = "#6bb6ff", high = "#006ad1",  name="Percent", breaks = 0.1*0:9, labels = percent(0.1*0:9)) +
                theme_classic() +
                labs(title = "") +
                labs(subtitle = paste(names_display[j], ':', strings_display[k],
                                      ", Mcnemar P = ",  mcnemarTest[j, k],
                                      ", OR = ",  fisherTest[j, k],
                                      sep="")) +
                theme(plot.title = element_text(hjust = 0.5)) +
                theme(plot.subtitle = element_text(hjust = 0.5)) +
                theme(axis.ticks = element_blank(),
                      text = element_text(size = 16),
                      axis.line.x = element_blank(),
                      axis.line.y = element_blank(),
                      legend.background = element_blank(),
                      legend.text = element_text( size = 8 ),
                      axis.title = element_text( size = 16 )) +
                xlab(transMethods[m]) + ylab(transMethods[i]) + guides(fill=FALSE)

            # Scatter Plot (Colored by Prevalence)
            corr <- cor(as.numeric(as.vector(df$qval_1)), as.numeric(df$qval_2), method="spearman")
            scatterplot.list[[j,k]] <- ggplot(df, aes(qval_1, qval_2, color = df$prevalence)) +
                ggtitle(paste(names[j], '/', strings2[k], ": Spearman ", round(corr, 2), sep='')) +
                geom_point(aes_string(color = df$prevalence), size = 6) +
                scale_color_gradient( low = 'red', high = 'blue', space = "Lab", guide_legend( title = 'Prevalence')) +
                xlab(sprintf('Min Q value (%s)', transMethods[i])) +
                ylab(sprintf('Min Q value (%s)', transMethods[m])) +
                theme_classic() +
                theme( panel.grid = element_blank(),
                       axis.title.y = element_text(size = 10),
                       axis.title.x = element_text(size = 10),
                       axis.text.x = element_text(size = 10),
                       plot.title = element_text( size = 10))
            }
        }

    # Save and Plot
    rownames(mcnemarTest)<-names
    colnames(mcnemarTest)<-c('diagnosis', 'active_vs_inactive')

    rownames(fisherTest)<-names
    colnames(fisherTest)<-c('diagnosis', 'active_vs_inactive')

    # Print
    ml1<-marrangeGrob(contingency_mcnemar.list, ncol=1, nrow = 1,  top = NULL)
    ggsave(sprintf('./Benchmarking/%s_vs_%s_Contingency_Plots_McNemar.pdf', transMethods[i], transMethods[m]), ml1)
    dev.off()

    # Print
    ml2<-marrangeGrob(scatterplot.list, ncol=1, nrow = 1,  top = NULL)
    ggsave(sprintf('./Benchmarking/%s_vs_%s_Scatter_Plots.pdf', transMethods[i], transMethods[m]), ml2)
    dev.off()

    }
}

###################
# Response Figure #
###################

#############
# Diagnosis #
#############

diagnosis_CLR_file<-file_names[[3]][str_detect(file_names[[3]], paste(strings[1], names[1], sep='_'))]
diagnosis_Vanilla_file<-file_names[[4]][str_detect(file_names[[4]], paste(strings[1], names[1], sep='_'))]

# Read the corresponding data frames
diagnosis_CLR <- data.table::fread(diagnosis_CLR_file, header = TRUE)
diagnosis_Vanilla <- data.table::fread(diagnosis_Vanilla_file, header = TRUE)

# Subset to ID and Qval
df_CLR<-diagnosis_CLR[, c('#', 'prevalence', 'minQ')]
colnames(df_CLR)<-c('ID', 'prevalence', 'qval_1')
df_Vanilla<-diagnosis_Vanilla[, c('#', 'prevalence', 'minQ')]
colnames(df_Vanilla)<-c('ID', 'prevalence', 'qval_2')

# Merge
df_diagnosis<-merge(df_CLR, df_Vanilla, c('ID', 'prevalence'))

# Create Table
nonVanilla_Significant<-df_diagnosis$ID[df_diagnosis$qval_1<FDR_Threshold]
nonVanilla_nonSignificant<-df_diagnosis$ID[df_diagnosis$qval_1>=FDR_Threshold]
Vanilla_Significant<-df_diagnosis$ID[df_diagnosis$qval_2<FDR_Threshold]
Vanilla_nonSignificant<-df_diagnosis$ID[df_diagnosis$qval_2>=FDR_Threshold]


# Contingency
A<-length(intersect(Vanilla_Significant, nonVanilla_Significant))
B<-length(intersect(Vanilla_Significant, nonVanilla_nonSignificant))
C<-length(intersect(Vanilla_nonSignificant, nonVanilla_Significant))
D<-length(intersect(Vanilla_nonSignificant, nonVanilla_nonSignificant))

# McNemar's Exact Test - Calculate P value
Mat<-matrix(c(A, B, C, D),nrow=2,ncol=2)
rownames(Mat)<-c('Significant', 'Not Significant')
colnames(Mat)<-c('Significant', 'Not Significant')
mcnemarTest<-ifelse(is.finite(round(mcnemar.test(Mat)$p.value,4)), round(mcnemar.test(Mat)$p.value,4), 1)
fisherTest<-round(fisher.test(Mat)$estimate,2)


# Create Data for Plot
Mat0<-as.data.frame(Mat)
Mat0<-rownames_to_column(Mat0, 'ID')
Mat_Final<-melt(Mat0)
Mat_Final$variable<-as.factor(Mat_Final$variable)
Mat_Final$ID<-as.factor(Mat_Final$ID)
Mat_Final$percent <- Mat_Final$value / sum(Mat_Final$value)
Mat_Final$label<-ifelse(Mat_Final$percent<1e-04, paste0(comma(Mat_Final$value), ' (', round(100 * Mat_Final$percent, 4), '%)'), paste0(comma(Mat_Final$value), ' (', round(100 * Mat_Final$percent, 2), '%)'))

# Contingency Plot (McNemar)
p_Diagnosis<-ggplot(Mat_Final, aes(variable, ID, fill=percent, label=label)) +
    geom_tile(colour = "black") +
    geom_text(aes(label = label), colour = "white", size=10) +
    scale_fill_gradient(low = "#6bb6ff", high = "#006ad1",  name="Percent", breaks = 0.1*0:9, labels = percent(0.1*0:9)) +
    theme_classic() +
    labs(title = "") +
    labs(subtitle = paste("Diagnosis",
                          ", McNemar P = ",  mcnemarTest,
                          ", OR = ",   round(fisherTest),
                          sep="")) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(plot.subtitle = element_text(hjust = 0.5)) +
    theme(axis.ticks = element_blank(),
          text = element_text(size = 30),
          axis.line.x = element_blank(),
          axis.line.y = element_blank(),
          legend.background = element_blank(),
          axis.text=element_text(size=30),
          legend.text = element_text( size = 8 ),
          axis.title = element_text( size = 30 )) +
    xlab('AST') + ylab('CLR') + guides(fill=FALSE)

#############
# Dysbiosis #
#############

Dysbiosis_CLR_file<-file_names[[3]][str_detect(file_names[[3]], paste(strings[2], names[1], sep='_'))]
Dysbiosis_Vanilla_file<-file_names[[4]][str_detect(file_names[[4]], paste(strings[2], names[1], sep='_'))]

# Read the corresponding data frames
Dysbiosis_CLR <- data.table::fread(Dysbiosis_CLR_file, header = TRUE)
Dysbiosis_Vanilla <- data.table::fread(Dysbiosis_Vanilla_file, header = TRUE)

# Subset to ID and Qval
df_CLR<-Dysbiosis_CLR[, c('#', 'prevalence', 'minQ')]
colnames(df_CLR)<-c('ID', 'prevalence', 'qval_1')
df_Vanilla<-Dysbiosis_Vanilla[, c('#', 'prevalence', 'minQ')]
colnames(df_Vanilla)<-c('ID', 'prevalence', 'qval_2')

# Merge
df_Dysbiosis<-merge(df_CLR, df_Vanilla, c('ID', 'prevalence'))

# Create Table
nonVanilla_Significant<-df_Dysbiosis$ID[df_Dysbiosis$qval_1<FDR_Threshold]
nonVanilla_nonSignificant<-df_Dysbiosis$ID[df_Dysbiosis$qval_1>=FDR_Threshold]
Vanilla_Significant<-df_Dysbiosis$ID[df_Dysbiosis$qval_2<FDR_Threshold]
Vanilla_nonSignificant<-df_Dysbiosis$ID[df_Dysbiosis$qval_2>=FDR_Threshold]


# Contingency
A<-length(intersect(Vanilla_Significant, nonVanilla_Significant))
B<-length(intersect(Vanilla_Significant, nonVanilla_nonSignificant))
C<-length(intersect(Vanilla_nonSignificant, nonVanilla_Significant))
D<-length(intersect(Vanilla_nonSignificant, nonVanilla_nonSignificant))

# McNemar's Exact Test - Calculate P value
Mat<-matrix(c(A, B, C, D),nrow=2,ncol=2)
rownames(Mat)<-c('Significant', 'Not Significant')
colnames(Mat)<-c('Significant', 'Not Significant')
mcnemarTest<-ifelse(is.finite(round(mcnemar.test(Mat)$p.value,4)), round(mcnemar.test(Mat)$p.value,4), 1)
fisherTest<-round(fisher.test(Mat)$estimate,2)


# Create Data for Plot
Mat0<-as.data.frame(Mat)
Mat0<-rownames_to_column(Mat0, 'ID')
Mat_Final<-melt(Mat0)
Mat_Final$variable<-as.factor(Mat_Final$variable)
Mat_Final$ID<-as.factor(Mat_Final$ID)
Mat_Final$percent <- Mat_Final$value / sum(Mat_Final$value)
Mat_Final$label<-ifelse(Mat_Final$percent<1e-04, paste0(comma(Mat_Final$value), ' (', round(100 * Mat_Final$percent, 4), '%)'), paste0(comma(Mat_Final$value), ' (', round(100 * Mat_Final$percent, 2), '%)'))

# Contingency Plot (McNemar)
p_Dysbiosis<-ggplot(Mat_Final, aes(variable, ID, fill=percent, label=label)) +
    geom_tile(colour = "black") +
    geom_text(aes(label = label), colour = "white", size=10) +
    scale_fill_gradient(low = "#6bb6ff", high = "#006ad1",  name="Percent", breaks = 0.1*0:9, labels = percent(0.1*0:9)) +
    theme_classic() +
    labs(title = "") +
    labs(subtitle = paste("Dysbiosis",
                          ", McNemar P = ",  mcnemarTest,
                          ", OR = ",  round(fisherTest),
                          sep="")) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(plot.subtitle = element_text(hjust = 0.5)) +
    theme(axis.ticks = element_blank(),
          text = element_text(size = 30),
          axis.line.x = element_blank(),
          axis.line.y = element_blank(),
          legend.background = element_blank(),
          axis.text=element_text(size=30),
          legend.text = element_text( size = 8 ),
          axis.title = element_text( size = 30 )) +
    xlab('AST') + ylab('CLR') + guides(fill=FALSE)



# Combine Figures
pdf("././Benchmarking/CLR_Response_Figure.pdf", height=8, width=22)
library(cowplot)
p<-plot_grid(p_Diagnosis, p_Dysbiosis, ncol = 2)
print(p)
dev.off()


