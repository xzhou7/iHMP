library(dada2) 
library(ShortRead)

packageVersion("dada2")

############################################################
seq.ns.st <- readRDS("~/Desktop/HMP Data/seqtab01.rds")

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seq.ns.st, method="consensus", multithread=TRUE, verbose=TRUE)
#Identified 10930 bimeras out of 72540 input sequences.

#write representative sequencing
uniquesToFasta(unqs=seqtab.nochim, fout="~/Desktop/HMP_DADA2_result/foo.fa")

dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seq.ns.st)


#0.9806065

#you should expect everything is above or equal to 400 basepair
table(nchar(colnames(seqtab.nochim)) >= 400)

taxa <- assignTaxonomy(colnames(seqtab.nochim), "~/Desktop/uBIOME_DADA2/RDP for DADA2_V18/rdp_train_set_18.fa.gz", multithread=T,verbose = T)
taxa <- addSpecies(taxa, "~/Desktop/uBIOME_DADA2/RDP for DADA2_V18/rdp_species_assignment_18.fa.gz", verbose = T)
taxa.print <- taxa
rownames(taxa.print) <- NULL

from <- colnames(seqtab.nochim)
to <- paste("ASV", c(1:dim(seqtab.nochim)[2]), sep="")
map = setNames(to, from)
colnames(seqtab.nochim) <- map[colnames(seqtab.nochim)]
row.names(taxa) <- map[row.names(taxa)]

write.csv(taxa, "~/Desktop/HMP_DADA2_result/tax.table_20200207.csv")
write.csv(seqtab.nochim, "~/Desktop/HMP_DADA2_result/otu.table_20200207.csv")


