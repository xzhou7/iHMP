#### Remove All the Variables from the Workspace ####
rm(list = ls())

#### Setting Directory ####
setwd('/Users/hmallick/Dropbox (Huttenhower Lab)/HMP2/analysis/differential_abundance/Vanilla')

#### Load libraries ####
library(data.table)
library(openxlsx)
library(stringr)


# Extract file names
subdir1<-'./diagnosis/'
subdir2<-'./active_vs_inactive/'
foo <-list.files(subdir1, pattern="*.tsv")
files1<-paste(subdir1, foo, sep='')
bar <-list.files(subdir2, pattern="*.tsv")
files2<-paste(subdir2, bar, sep='')
files0<-append(files1, files2)

# Re-order files to impose systematic ordering
files<-c(files0[str_detect(files0, '/diagnosis/IBD_relations_bugs')],
         files0[str_detect(files0, 'diagnosis/IBD_relations_metabolites')],
         files0[str_detect(files0, 'diagnosis/IBD_relations_viruses')],
         files0[str_detect(files0, 'diagnosis/IBD_relations_bugkoRNADNAs')],
         files0[str_detect(files0, 'diagnosis/IBD_relations_pwyDNAs')],
         files0[str_detect(files0, 'diagnosis/IBD_relations_pwyRNAs')],
         files0[str_detect(files0, 'diagnosis/IBD_relations_pwyRNADNAs')],
         files0[str_detect(files0, 'diagnosis/IBD_relations_ecDNAs')],
         files0[str_detect(files0, 'diagnosis/IBD_relations_ecRNAs')],
         files0[str_detect(files0, 'diagnosis/IBD_relations_koDNAs')],
         files0[str_detect(files0, 'diagnosis/IBD_relations_koRNAs')],
         files0[str_detect(files0, 'diagnosis/IBD_relations_proteins')],
         files0[str_detect(files0, 'diagnosis/IBD_relations_proteinECs')],
         files0[str_detect(files0, 'diagnosis/IBD_relations_proteinKOs')],
         files0[str_detect(files0, 'Active_IBD_relations_bugs')],
         files0[str_detect(files0, 'Active_IBD_relations_metabolites')],
         files0[str_detect(files0, 'Active_IBD_relations_viruses')],
         files0[str_detect(files0, 'Active_IBD_relations_bugkoRNADNAs')],
         files0[str_detect(files0, 'Active_IBD_relations_pwyDNAs')],
         files0[str_detect(files0, 'Active_IBD_relations_pwyRNAs')],
         files0[str_detect(files0, 'Active_IBD_relations_pwyRNADNAs')],
         files0[str_detect(files0, 'Active_IBD_relations_ecDNAs')],
         files0[str_detect(files0, 'Active_IBD_relations_ecRNAs')],
         files0[str_detect(files0, 'Active_IBD_relations_koDNAs')],
         files0[str_detect(files0, 'Active_IBD_relations_koRNAs')],
         files0[str_detect(files0, 'Active_IBD_relations_proteins')],
         files0[str_detect(files0, 'Active_IBD_relations_proteinECs')],
         files0[str_detect(files0, 'Active_IBD_relations_proteinKOs')])

# Remove Redundant Files
files<-files[!str_detect(files, 'metabolites_Sena')]
files<-files[!str_detect(files, 'viruses_metaphlan')]


# Read files
list_of_datasets<-list()
for (i in 1:length(files)) {
    dat <- data.table::fread(files[i], header = TRUE)
    data<-dat[complete.cases(dat[,-c(1:2)]),]
    # Rename column names
    if (sum(colnames(data) %in% '#') == 1)     temp<-data.table::setnames(data, old = '#', new = 'Feature')
    if (sum(colnames(data) %in% 'prevalence') == 1)     temp<-data.table::setnames(data, old = 'prevalence', new = 'Prevalence')
    if (sum(colnames(data) %in% 'coefCD') == 1)     temp<-data.table::setnames(data, old = 'coefCD', new = 'Coefficient (CD)')
    if (sum(colnames(data) %in% 'coefUC') == 1)     temp<-data.table::setnames(data, old = 'coefUC', new = 'Coefficient (UC)')
    if (sum(colnames(data) %in% 'tvalCD') == 1)     temp<-data.table::setnames(data, old = 'tvalCD', new = 't Statistic (CD)')
    if (sum(colnames(data) %in% 'tvalUC') == 1)     temp<-data.table::setnames(data, old = 'tvalUC', new = 't Statistic (UC)')
    if (sum(colnames(data) %in% 'pvalCD') == 1)     temp<-data.table::setnames(data, old = 'pvalCD', new = 'P Value (CD)')
    if (sum(colnames(data) %in% 'pvalUC') == 1)     temp<-data.table::setnames(data, old = 'pvalUC', new = 'P Value (UC)')
    if (sum(colnames(data) %in% 'qvalCD') == 1)     temp<-data.table::setnames(data, old = 'qvalCD', new = 'Q Value (CD)')
    if (sum(colnames(data) %in% 'qvalUC') == 1)     temp<-data.table::setnames(data, old = 'qvalUC', new = 'Q Value (UC)')
    if (sum(colnames(data) %in% 'minQ') == 1)     temp<-data.table::setnames(data, old = 'minQ', new = 'Minimum Q Value')
    if (sum(colnames(data) %in% 'zvalCD') == 1)     temp<-data.table::setnames(data, old = 'zvalCD', new = 'Z Statistic (CD)')
    if (sum(colnames(data) %in% 'zvalUC') == 1)     temp<-data.table::setnames(data, old = 'zvalUC', new = 'Z Statistic (UC)')
    if (sum(colnames(data) %in% 'coefActiveCD') == 1)     temp<-data.table::setnames(data, old = 'coefActiveCD', new = 'Dysbiosis Coefficient (CD)')
    if (sum(colnames(data) %in% 'coefActiveUC') == 1)     temp<-data.table::setnames(data, old = 'coefActiveUC', new = 'Dysbiosis Coefficient (UC)')
    if (sum(colnames(data) %in% 'coefActivenonIBD') == 1) temp<-data.table::setnames(data, old = 'coefActivenonIBD', new = 'Dysbiosis Coefficient (non-IBD)')
    if (sum(colnames(data) %in% 'tvalActiveCD') == 1)     temp<-data.table::setnames(data, old = 'tvalActiveCD', new = 'Dysbiosis t Statistic (CD)')
    if (sum(colnames(data) %in% 'tvalActiveUC') == 1)     temp<-data.table::setnames(data, old = 'tvalActiveUC', new = 'Dysbiosis t Statistic (UC)')
    if (sum(colnames(data) %in% 'tvalActivenonIBD') == 1)     temp<-data.table::setnames(data, old = 'tvalActivenonIBD', new = 'Dysbiosis t Statistic (non-IBD)')
    if (sum(colnames(data) %in% 'pvalActiveCD') == 1)     temp<-data.table::setnames(data, old = 'pvalActiveCD', new = 'Dysbiosis P Value (CD)')
    if (sum(colnames(data) %in% 'pvalActiveUC') == 1)     temp<-data.table::setnames(data, old = 'pvalActiveUC', new = 'Dysbiosis P Value (UC)')
    if (sum(colnames(data) %in% 'pvalActivenonIBD') == 1)     temp<-data.table::setnames(data, old = 'pvalActivenonIBD', new = 'Dysbiosis P Value (non-IBD)')
    if (sum(colnames(data) %in% 'qvalActiveCD') == 1)     temp<-data.table::setnames(data, old = 'qvalActiveCD', new = 'Dysbiosis Q Value (CD)')
    if (sum(colnames(data) %in% 'qvalActiveUC') == 1)     temp<-data.table::setnames(data, old = 'qvalActiveUC', new = 'Dysbiosis Q Value (UC)')
    if (sum(colnames(data) %in% 'qvalActivenonIBD') == 1)     temp<-data.table::setnames(data, old = 'qvalActivenonIBD', new = 'Dysbiosis Q Value (non-IBD)')
    if (sum(colnames(data) %in% 'zvalActiveCD') == 1)     temp<-data.table::setnames(data, old = 'zvalActiveCD', new = 'Dysbiosis Z Statistic (CD)')
    if (sum(colnames(data) %in% 'zvalActiveUC') == 1)     temp<-data.table::setnames(data, old = 'zvalActiveUC', new = 'Dysbiosis Z Statistic (UC)')
    if (sum(colnames(data) %in% 'zvalActivenonIBD') == 1)     temp<-data.table::setnames(data, old = 'zvalActivenonIBD', new = 'Dysbiosis Z Statistic (non-IBD)')
    temp$Feature<-gsub("_", " ", temp$Feature) # Nicefy feature names
    list_of_datasets[[i]]<- temp

}
names(list_of_datasets)<-paste('S', 1:length(list_of_datasets), sep='')
style <- createStyle(textDecoration = "bold", halign = "center")
openxlsx::write.xlsx(list_of_datasets,
                     file = "/Users/hmallick/Dropbox (Huttenhower Lab)/HMP2/manuscript/tables/Supplementary_Tables_DA.xlsx",
                     headerStyle=style, colWidths='auto')







