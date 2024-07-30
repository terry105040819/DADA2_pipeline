# DADA2 installation ------------------------------------------------------


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("dada2", version = "3.19")

library(dada2); packageVersion("dada2")



# Illumina Miseq data example ---------------------------------------------

path <- "~/Desktop/MiSeq_SOP/"
list.files(path)

#Sort Illumina Forford and Reverse sequence
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

#Grep sample name form fastq file
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#Visulization Q score plot
plotQualityProfile(fnFs)

## Trim and Filtering Forward and Revserse sequence form Illumina

# Create QC filtered output folder
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names


out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)

# Error rate 
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)


dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errF, multithread=TRUE)


dadaFs[[2]]


vec <- 1:10
vec2 <- 2:11
vec <- paste(vec,vec2)
