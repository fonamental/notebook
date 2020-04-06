#!/bin/bash

#$ -S /bin/bash
#$ -o analysis/logs/dada2.$JOB_ID
#$ -j y
#$ -cwd
#$ -V
#$ -P protist
#$ -N dada2
#$ -q all.q,share.q

module load gcc/4.9.0
module load Rstats/R-3.4.1

Rscript scripts/dada2.R \
data/surface \
surface \
analysis/dada2_surface/ \
240,170 \
2,2 \
MSTAReuk

Rscript scripts/dada2.R \
data/vertical_profiles/DNA \
verticalDNA \
analysis/dada2_verticalDNA/ \
240,170 \
2,2 \
trimmedFilt200

Rscript scripts/dada2.R \
data/vertical_profiles/RNA \
verticalRNA \
analysis/dada2_verticalRNA/ \
240,170 \
2,2 \
trimmedFilt200

