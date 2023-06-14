
source("./common/merge_metadata.r")

test_matching <- function(pclpath, type, use_ssc=F) {
    if (is.list(pclpath)) {
        pclnames <- rownames(pclpath$meta)
    } else {
        if (endsWith(pclpath, ".gz")) {
            f <- gzfile(pclpath, "r")
        } else {
            f <- file(pclpath, "r")
        }
        hline <- readLines(f, 1)
        pclnames <- strsplit(hline, "\t")[[1]]
        close(f)

        pclnames <- pclnames[!startsWith(pclnames, "#") & (pclnames!="")]
    }

    if (use_ssc) {
        mtnames <- as.character(hmp2_sample_metadata$Site.Sub.Coll[hmp2_sample_metadata$data_type == type])
    } else {
        mtnames <- as.character(hmp2_sample_metadata$External.ID[hmp2_sample_metadata$data_type == type])
    }

    if (!all(mtnames %in% pclnames)) {
        mismatch <- mtnames[!(mtnames %in% pclnames)]
        cat(sprintf("These %d samples are in the %s metadata but not in the pcl:\n", length(mismatch), type))
        print(mismatch)
    }
    if (any(duplicated(mtnames))) {
        dups <- mtnames[duplicated(mtnames)]
        cat(sprintf("These %d samples duplicated in the %s metadata:\n", length(dups), type))
        print(dups)
    }
    if (!all(pclnames %in% mtnames)) {
        mismatch <- pclnames[!(pclnames %in% mtnames)]
        cat(sprintf("These %d samples are in the pcl but not in the %s metadata:\n", length(mismatch), type))
        print(mismatch)
    }
    if (any(duplicated(pclnames))) {
        dups <- pclnames[duplicated(pclnames)]
        cat(sprintf("These %d samples duplicated in the %s pcl:\n", length(dups), type))
        print(dups)
    }

    cat(sprintf("Sample count for %s: %d\n", type, sum(pclnames %in% mtnames)))
}

source("./common/load_metabolites.r")

test_matching(file.path(HMP2_data,"mgx","taxonomic_profiles.pcl.tsv"), "metagenomics")
test_matching(file.path(HMP2_data,"mtx","pathabundance_relab.pcl.tsv.gz"), "metatranscriptomics")
test_matching(metabolites.pcl, "metabolomics")
test_matching(file.path(HMP2_data,"mpx","hmp2_MPX_not_normalized_ecs.names.pcl.tsv"), "proteomics")
test_matching(file.path(HMP2_data,"mvx","HMP2.Virome.VirMAP.pcl.tsv"), "viromics")
test_matching(file.path(HMP2_data,"16s","biopsy_16s_otus_Rfriendly.pcl.tsv"), "biopsy_16S")
test_matching(file.path(HMP2_data,"htx","host_tx_counts.pcl.tsv"), "host_transcriptomics")
test_matching(file.path(HMP2_data,"serology","hmp2_serology_Compiled_ELISA_Data.pcl.tsv"), "serology")





