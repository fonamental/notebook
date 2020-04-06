# DADA2 pipeline

## Load packages

library(stringr)
library(tidyverse)
library(dada2); packageVersion("dada2")

## Load Rscript args

args <- commandArgs(trailingOnly = TRUE)

## Variables to change for EACH RUN

path<- args[1] # where the fastq files are located
name.run <- args[2] # name of the run
output <-args[3] # results/diagnostic output directory
trunclen <- as.integer(strsplit(args[4], ",")[[1]]) # trimming positions for forward and reverse reads
maxee <- as.integer(strsplit(args[5], ",")[[1]]) # max error for forward and reverse reads
pattern <- args[6] # pattern for data files we want to work with

## Taking filenames

fnFs <- sort(list.files(path, pattern=paste0(pattern,'_R1.fastq'), full.names = TRUE))
fnRs <- sort(list.files(path, pattern=paste0(pattern,'_R2.fastq'), full.names = TRUE))

head(fnFs)
head(fnRs)

sample.names <- sapply(strsplit(fnFs, "-"), `[`, 2) # adapt "-" for your delimiter

## Filtering and trimming

dir.create(file.path(path, "filtered"), showWarnings = FALSE) # create a 'filtered' directory

filt_path <- file.path(path, "filtered") # path to 'filtered' directory
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=trunclen,
  maxN=0, maxEE=maxee, truncQ=2, rm.phix=TRUE,
  compress=TRUE, multithread=TRUE)

## Learning error rates

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

### Plotting error rates

err.plotf <- plotErrors(errF, nominalQ=TRUE)
ggsave(str_c(output,"errors_",name.run,"_fwd.pdf"),plot=err.plotf)
err.plotr <- plotErrors(errR, nominalQ=TRUE)
ggsave(str_c(output,"errors_",name.run,"_rev.pdf"),plot=err.plotr)

## Dereplication

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

names(derepFs) <- sample.names # name the derep-class objects by the sample names
names(derepRs) <- sample.names

## Sample inference

dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

## Merge paired reads

mergers <- mergePairs(dadaFs, derepFs,
  dadaRs, derepRs,
  verbose=TRUE)

## Construct sequence table

seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))

saveRDS(seqtab, str_c(output,name.run,"seqtab.rds"))

## Remove chimeras

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", 
  multithread=TRUE, verbose=TRUE)

dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

saveRDS(seqtab.nochim, str_c(output,name.run,"seqtab_nochim.rds"))

## Track reads through the pipeline

getN <- function(x) sum(getUniques(x))

track <- cbind(out, sapply(dadaFs, getN),
  sapply(mergers, getN),
  rowSums(seqtab),
  rowSums(seqtab.nochim))

colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names

class(track)

saveRDS(track, str_c(output,name.run,"track_analysis.rds"))
