#DATA2
#Xin Zhou
#https://benjjneb.github.io/dada2/tutorial.html
#updated 2021-02-17

library("dada2")

path1 <- "~/Desktop/UBIOMEData_Ori/01_uBiome_DataSharing_Sept-2017 /fastq_files/"
path2 <- "~/Desktop/UBIOMEData_Ori/02_uBiome_DataSharing-4-Jan-2018/fastq_files/"
path3 <- "~/Desktop/UBIOMEData_Ori/03_uBiome_DataSharing-27-Mar-2018/fastq_files/"
path4 <- "~/Desktop/UBIOMEData_Ori/04_uBiome-DataSharing-18-07-2018/fastq_files/"
path5 <- "~/Desktop/UBIOMEData_Ori/05_uBiome-DataSharing-29-Nov-2018/fastq_files/"
path <- "~/Desktop/UBIOME Filtered"

list.files(path1)
fnFs1 <- sort(list.files(path1, pattern="R1_001.fastq", full.names = TRUE))
fnFs2 <- sort(list.files(path2, pattern="R1_001.fastq", full.names = TRUE))
fnFs3 <- sort(list.files(path3, pattern="R1_001.fastq", full.names = TRUE))
fnFs4 <- sort(list.files(path4, pattern="R1_001.fastq", full.names = TRUE))
fnFs5 <- sort(list.files(path5, pattern="R1_001.fastq", full.names = TRUE))
fnFs <-  c(fnFs1, fnFs2, fnFs3, fnFs4, fnFs5)

sample.names <- sapply(strsplit(basename(fnFs), "_R1_"), `[`, 1)

plotQualityProfile(fnFs[1:4])

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))

#start 11:43 am, ends 13:17 pm 
out <- filterAndTrim(fnFs, filtFs, maxN=0, maxEE=2, truncQ=2, rm.phix=TRUE,compress=TRUE, multithread=T)
#head(out)

errF <- learnErrors(filtFs, multithread=TRUE)

errorrates_F_plot <- plotErrors(errF, nominalQ=TRUE)
errorrates_F_plot

derepFs <- derepFastq(filtFs, verbose=TRUE)

dadaFs <- dada(derepFs, err=errF, multithread=TRUE)

################using only forward reads######################################
seqtab_F <- makeSequenceTable(dadaFs)
seqtab.nochim_F <- removeBimeraDenovo(seqtab_F, method="consensus", multithread=TRUE, verbose=TRUE)

uniquesToFasta(unqs=seqtab.nochim_F, fout="~/Desktop/UBIOME_DADA2_result_V18/representOTU.fa")

taxa_F <- assignTaxonomy(seqtab.nochim_F, "~/Desktop/uBIOME_DADA2/RDP for DADA2_V18/rdp_train_set_18.fa.gz", multithread=TRUE, verbose = T)
taxa_F <- addSpecies(taxa_F, "~/Desktop/uBIOME_DADA2/RDP for DADA2_V18/rdp_species_assignment_18.fa.gz",verbose = T)

taxa.print_F <- taxa_F
rownames(taxa.print_F) <- NULL

from <- colnames(seqtab.nochim_F)
to <- paste("OTU", c(1:dim(seqtab.nochim_F)[2]), sep="_")
map = setNames(to, from)
colnames(seqtab.nochim_F) <- map[colnames(seqtab.nochim_F)]
row.names(taxa_F) <- map[row.names(taxa_F)]
row.names(seqtab.nochim_F) <- paste("UBMkid",sapply(strsplit(row.names(seqtab.nochim_F), "_"), "[[", 1), sep="")

write.csv(taxa_F, "~/Desktop/SkinOraltax.table_F_20210217.csv")
write.csv(seqtab.nochim_F, "~/Desktop/SkinOralotu.table_F_20210217.csv")
