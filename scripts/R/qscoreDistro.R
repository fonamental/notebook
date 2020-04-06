# DADA2 pre-processing. Q score distribution

## Load packages

library(stringr)
library(tidyverse)
library(dada2); packageVersion("dada2")

## Load Rscript args

args <- commandArgs(trailingOnly = TRUE)

## Variables to change for EACH RUN

path<- args[1] # where the fastq files are located
output <- args[2] # where we will dump the results//diagnostic
run.name <- args[3] # Run name. A directory for the values will be created
pattern <- args[4] # pattern for data files we want to work with

## Taking filenames

fnFs <- sort(list.files(path, pattern = paste0(pattern,'_R1.fastq')))
fnRs <- sort(list.files(path, pattern= paste0(pattern,'_R2.fastq')))

sample.names <- sapply(strsplit(fnFs, "-"), `[`, 2) 

fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)

### Since a qprofile dir can be not present we check and create it

dir.create(file.path(output, "qprofile"), showWarnings = FALSE)
dir.create(file.path(output, "qprofile",run.name), showWarnings = FALSE)

out.diag <- file.path(output, "qprofile", run.name)

## Qscore distribution

for(l in 1:12){
  qF <- plotQualityProfile(fnFs[l]) +
    ggtitle("Forward reads")
  ggsave(plot=qF, path= out.diag,
         device="pdf",
         filename = paste0(sample.names[l],"_F.pdf"))
  qR <- plotQualityProfile(fnRs[l]) +
    ggtitle("Reverse reads")
  ggsave(plot=qR, path= out.diag,
         device="pdf",
         filename = paste0(sample.names[l],"_R.pdf"))}
