#!/bin/bash

#$ -S /bin/bash
#$ -o analysis/logs/qscore.$JOB_ID
#$ -j y
#$ -cwd
#$ -V
#$ -P protist
#$ -N dada2qscore
#$ -q all.q,share.q

module load gcc/4.9.0
module load Rstats/R-3.4.1

Rscript scripts/qscoreDistro.R data/surface analysis surface MSTAReuk

Rscript scripts/qscoreDistro.R data/vertical_profiles/DNA analysis verticalDNA trimmedFilt200

Rscript scripts/qscoreDistro.R data/vertical_profiles/RNA analysis verticalRNA trimmedFilt200
