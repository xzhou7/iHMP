#DATA2 HMP ST+NS local
#Xin Zhou
#https://benjjneb.github.io/dada2/tutorial.html
#https://wemp.app/posts/e49b0d5b-a7c9-43d9-86a3-c735f73b8eb2
#updated Feb. 07 2020


library("dada2")

path <- "~/Desktop/HMP Data/HMP_ST_NS"
pathout<- "~/Desktop/HMPfiltered"
list.files(path)
fastqfile <- sort(list.files(path, pattern = "clean.dehost.fastq", full.names = T))
JAXID <- lapply(strsplit(basename(fastqfile), ".clean"), `[`, 1)

#plotQualityProfile(fastqfile)

fastq.flitered <- file.path(pathout, "filtered", paste0(JAXID, ".flitered.clean.fastq.gz"))

out <- filterAndTrim(fastqfile, fastq.flitered, maxN = 5, truncQ=2, maxEE=2, truncLen=400, rm.phix=TRUE,compress=TRUE, multithread=T)

#################################################################
#change fliter path
#################################################################
filtpath <- "~/Desktop/HMPfiltered/filtered"
list.files(filtpath)
filts <- list.files(filtpath, pattern="fastq.gz", full.names=TRUE)
sample.names <- sapply(strsplit(basename(filts), ".clean"), `[`, 1)
names(filts) <- sample.names

# Learn error rates
set.seed(100)
err <- learnErrors(filts, nbases = 1e8, multithread=TRUE, randomize=TRUE)
plotErrors(err, nominalQ=TRUE)

# Infer sequence variants (~90 minutes for ~2000 samples)
dds <- vector("list", length(sample.names))
names(dds) <- sample.names
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derep <- derepFastq(filts[[sam]])
  dds[[sam]] <- dada(derep, err=err, multithread=TRUE)
}

# Construct sequence table and write to disk
seqtab <- makeSequenceTable(dds)
saveRDS(seqtab, "~/Desktop/HMP Data/seqtab01.rds")

seqtab_old <- readRDS("~/Desktop/HMP Data/seqtab_V16.rds")

