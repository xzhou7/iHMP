#reference:
#https://zenodo.org/record/4310151#.YB7gx2RKjfs

#RDP taxonomic training data formatted for DADA2 (RDP trainset 18/release 11.5)
#Author: Callahan, Benjamin

## The RDP trainset data was downloaded from: https://sourceforge.net/projects/rdp-classifier/files/RDP_Classifier_TrainingData/
path <- "~/Desktop/RDP/RDPClassifier_16S_trainsetNo18_rawtrainingdata"
dada2:::makeTaxonomyFasta_RDP(file.path(path, "trainset18_062020.fa"), 
                              file.path(path, "trainset18_db_taxid.txt"), 
                              "~/tax/rdp_train_set_18.fa.gz")
## Download "ten_16s.100.fa" from Robert Edgar's taxonomy testing page: https://drive5.com/taxxi/doc/fasta_index.html
dada2:::tax.check("~/tax/rdp_train_set_18.fa.gz", "~/Desktop/ten_16s.100.fa")

## This function creates the dada2 assignSpecies fasta file for the RDP from the RDP's _Bacteria_unaligned.fa file.
dada2:::makeSpeciesFasta_RDP("~/Desktop/RDP/current_Bacteria_unaligned.fa", "~/tax/rdp_species_assignment_18.fa.gz")
dada2:::tax.check("~/tax/rdp_species_assignment_18.fa.gz", "~/Desktop/ten_16s.100.fa", mode="species")