# Denoising-Illumina-reads
Denoising paired-end reads using DADA2 Workflow in R
The script used in R to denoise my sequences is:

## Libraries:

source("https://bioconductor.org/biocLite.R")
biocLite("biomformat")
library(biomformat)

source("https://bioconductor.org/biocLite.R")
biocLite("dada2")

source("https://bioconductor.org/biocLite.R")
biocLite("DECIPHER")

library(dada2)
library(DECIPHER)
library(phangorn)
library(ggplot2)
library(phyloseq)

##Specify path:
SEQS <- c("Tutorial DADA2/")
SEQS
list.files(SEQS)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(SEQS, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(SEQS, pattern="_R2_001.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#specify the full path to the fnFs and fnRs:
fnFs<-file.path(SEQS, fnFs)
fnRs<- file.path(SEQS, fnRs)

# plot quality profile:

plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(SEQS, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(SEQS, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(200,140),maxN=0, maxEE=c(2,4), truncQ=2, rm.phix=TRUE,
compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE
head(out)

## Error rates:

errF <- learnErrors(filtFs, multithread=FALSE)
errR <- learnErrors(filtRs, multithread=FALSE)

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

## Dereplication step:

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

## DADA Algorithm:

dadaFs <- dada(derepFs, err=errF, multithread=FALSE)
dadaRs <- dada(derepRs, err=errR, multithread=FALSE)
dadaFs[[1]]
dadaRs[[1]]
##merge paired reads:

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE, minOverlap = 10, maxMismatch = 1)

# Inspect the merger data.frame from the first sample

head(mergers[[1]])

##construct sequence table:

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

##remove chimeras:

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=FALSE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))

# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

##asign taxonomy:

taxa <- assignTaxonomy(seqtab.nochim, "Tutorial DADA2/gg_13_8_train_set_97.fa.gz", multithread=FALSE)
taxa <- addSpecies(taxa, "Tutorial DADA2/rdp_species_assignment_14.fa.gz")
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
