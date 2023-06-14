
source("./common/load_bugs.r")
source("./common/load_ecs.r")
source("./common/load_metabolites.r")
source("./common/load_fecalcal.r")
source("./common/load_serology.r")
source("./common/load_proteins.r")
source("./common/load_biopsies.r")
source("./common/disease_activity.r")

fecalcal.pcl$ns
# 652
length(unique(fecalcal.pcl$meta$subject))
# 98

hbi.pcl <- ec.dna.unstrat.pcl %>%
    pcl.filter.s(!is.na(hbi))
hbi.pcl$ns
# 680
length(unique(hbi.pcl$meta$subject))
# 65

sccai.pcl <- ec.dna.unstrat.pcl %>%
    pcl.filter.s(!is.na(sccai))
sccai.pcl$ns
# 429
length(unique(sccai.pcl$meta$subject))
# 38

ec.dna.unstrat.pcl$ns
# 1595
length(unique(ec.dna.unstrat.pcl$meta$subject))
# 130

ec.rna.unstrat.pcl$ns
# 818
length(unique(ec.rna.unstrat.pcl$meta$subject))
# 109
sum(!is.na(merge_disease_activity(ec.rna.unstrat.pcl)$meta$active))
# 785
length(unique(ec.rna.unstrat.pcl$meta$subject[!is.na(merge_disease_activity(ec.rna.unstrat.pcl)$meta$active)]))
# 109

metabolites.pcl$ns
# 546
length(unique(metabolites.pcl$meta$subject))
# 106
sum(!is.na(merge_disease_activity(metabolites.pcl)$meta$active))
# 461
length(unique(metabolites.pcl$meta$subject[!is.na(merge_disease_activity(metabolites.pcl)$meta$active)]))
# 106

# Adjacent pairs of MBX samples
sum(table(metabolites.pcl$meta$subject)-1)

serology.pcl$ns
# 211
length(unique(serology.pcl$meta$subject))
# 68
sum(!is.na(merge_disease_activity(serology.pcl, lenience=2)$meta$active))
# 146
length(unique(serology.pcl$meta$subject[!is.na(merge_disease_activity(serology.pcl, lenience=2)$meta$active)]))
# 61

biopsy_htx.counts.pcl$ns
# 249
length(unique(biopsy_htx.counts.pcl$meta$subject))
# 91

proteins.pcl$ns
# 449
length(unique(proteins.pcl$meta$subject))
# 89

viruses.pcl$ns
# 619
length(unique(viruses.pcl$meta$subject))
# 105
sum(!is.na(merge_disease_activity(viruses.pcl)$meta$active))
# 592
length(unique(viruses.pcl$meta$subject[!is.na(merge_disease_activity(viruses.pcl)$meta$active)]))
# 104

ileum <- !is.na(biopsy_htx.counts.pcl$meta$biopsy_location) & biopsy_htx.counts.pcl$meta$biopsy_location=="Ileum"
sum(ileum)
# 88
length(unique(biopsy_htx.counts.pcl$meta$subject[ileum]))
# 84

ileum_inf <- !is.na(biopsy_htx.counts.pcl$meta$biopsy_location) & biopsy_htx.counts.pcl$meta$biopsy_location=="Ileum"
sum(ileum)
# 88
length(unique(biopsy_htx.counts.pcl$meta$subject[ileum]))
# 84

hmp1ii.pcl <- pcl.read(file.path(HMP2_data, "hmp1-ii", "hmp1-II_metaphlan2-mtd-qcd.tsv"))
hmp1ii.pcl.stool <- hmp1ii.pcl %>% pcl.filter.s(STSite=="Stool")
hmp1ii.pcl.stool$ns
# 553
length(unique(hmp1ii.pcl.stool$meta$RANDSID))
# 249
