#!/bin/bash

#######################
#######SGE BLOCK ######
#######################
#$ -cwd  
#$ -j y  
###$ -pe qiime 24
#$ -o ../
#$ -V
#$ -S /bin/bash
#$ -q share.q,all.q
###$ -l hostname=compute-0-2
#$ -N SAGs_pipeline
#######################

#$ -M ramiro.logares@icm.csic.es   # email to submit job status, remove one # to activate
#$ -m abe

#######################


echo "##########################################################################################"
echo "########################## PRE-ANOTTATION PIPELINE #######################################"
echo "######################   Ramiro Logares | V1 (07/1/2017)   ###############################"
echo "##########################################################################################"
echo "Note: this pipeline replaces sags_pipeline"

echo "##########################################################################################"


echo "Steps included in this pre-annotation pipeline"

echo "# 1) Remove contigs <1kb"
echo "# 2) Run Quast"
echo "# 3) repeat masking (RepeatMasker) # mask repetitive and low complexity genomic regions"
echo "# 4) tRNA detection EUKS (tRNAscan-SE-1.3.1)"
echo "# 5) rRNA detection EUKS  (RNAmmer)"
echo "# 6) rfam_scan"
echo "# 7) Cegma"
echo "# 8) Busco"
 
echo "##########################################################################################"

mv pre_gene_pred_pipeline_v1.6.sh ../


module load BLAST/2.2.28+

##Programs ##

PATH=$PATH:/home/rlogares/SOFTWARE/repeatmasker/RepeatMasker
PATH=$PATH:/home/rlogares/apps/quast-3.2/quast-master
PATH=$PATH:/home/rlogares/scripts
##PATH=$PATH:/home/rlogares/apps/GAEMR/GAEMR-1.0.1/bin
PATH=$PATH:/home/rlogares/SOFTWARE/rnamer-1.2.src:/home/rlogares/apps:/opt/perl/bin:/home/rlogares/perl5/perlbrew/bin:/home/rlogares/SOFTWARE/tRNAscan/tRNAscan-SE-1.3.1
export PATH

## export PYTHONPATH=/home/rlogares/apps/GAEMR/GAEMR-1.0.1:$PYTHONPATH




## Variables ###


## This should be run in a directory, which contains subdirectories, each containing a "contigs.fasta" file

echo "Starting pre-gene prediction pipeline"
echo ""
echo "Remove contigs <1kb"
echo ""
date

for i in $(ls -d *); do cd $i; remove_short_seqs.pl 1000  *  > contigs_min1kb.fasta ;cd ..; done


echo "Generating subdirectories: masked_contigs RNAmmer tRNAscan cmscan  augustus quast repeatmasker cegma "


for i in $(ls -d *); do cd $i; mkdir masked_contigs RNAmmer tRNAscan cmscan  augustus quast repeatmasker cegma ;cd ..; done


echo " Running Quast: get general estimates of assemblies/bins"
date

for i in $(ls -d *); do cd $i; quast.py  --min-contig 1000  -o quast/quast contigs_min1kb.fasta ;cd ..; done



## Gaemr

## GAEMR.py -c $contigs -o gaemr/GAEMR



echo " Running Repeatmasker"
date

for i in $(ls -d *); do cd $i; RepeatMasker  -gff -dir repeatmasker/  contigs_min1kb.fasta ;cd ..; done

date



echo "Running tRNAscan"
date

for i in $(ls -d *); do cd $i; perlbrew exec "tRNAscan-SE -o tRNAscan/contigsMaskedtRNAscanSE repeatmasker/contigs_min1kb.fasta.masked" ;cd ..; done
date

echo "Converting tRNAscan output to GFF"

for i in $(ls -d *); do cd $i/tRNAscan; for line in $(ls contigsMaskedtRNAscanSE); do perl < $line -lne '($contig,$trnanum,$st,$en, $type, $codon, $intron_st, $intron_en, $score)=split /\s+/;next unless $trnanum =~ /^\d+$/;print "$contig\ttRNAscan-SE-1.3\ttRNA\t" . ($st < $en ? "$st\t$en\t$score\t+\t" : "$en\t$st\t$score\t-\t") . ".\tID=$contig\_tRNA$trnanum;Name=tRNA_$type\_$codon" ' > $line.gff; done ;cd ../../; done


## BASAL command for reference
## for line in $(ls contigsMaskedtRNAscanSE); do perl < $line -lne '($contig,$trnanum,$st,$en, $type, $codon, $intron_st, $intron_en, $score)=split /\s+/;next unless $trnanum =~ /^\d+$/;print "$contig\ttRNAscan-SE-1.3\ttRNA\t" . ($st < $en ? "$st\t$en\t$score\t+\t" : "$en\t$st\t$score\t-\t") . ".\tID=$contig\_tRNA$trnanum;Name=tRNA_$type\_$codon" ' > $line.gff; done


echo "Running RNAmmer"
date

for i in $(ls -d *); do cd $i/RNAmmer; perlbrew exec "rnammer -S euk -m lsu,ssu,tsu -gff contig.rnammer -h contig.hmm.report  < ../contigs_min1kb.fasta" ;cd ../../; done  ## run with perl-5.20 installed in home

for i in $(ls -d *); do cd $i/RNAmmer; rnammer2gff.pl contig.rnammer ;cd ../../; done 

date


echo "Masking"
date

for i in $(ls -d *); do cd $i; maskFastaFromBed -fi repeatmasker/contigs_min1kb.fasta.masked -bed tRNAscan/contigsMaskedtRNAscanSE.gff -fo masked_contigs/contigs.masked_repeatmasker_trnascan_min1kb.fna ;cd ..; done

echo "Masked file: masked_contigs/contigs.masked_repeatmasker_trnascan_min1kb.fna"


echo "Running Cegma"
date

for i in $(ls -d *); do cd $i; perlbrew exec "cegma -g masked_contigs/contigs.masked_repeatmasker_trnascan_min1kb.fna -o cegma/cegma_min1kb" ;cd ..; done


echo "Running BUSCO"
echo "----------------------------------------------------------"
echo ""

python3=/opt/python/bin/python3.2 
busco=/home/rlogares/apps/busco/BUSCO_v1.1b1/BUSCO_v1.1b1.py

## Reference datasets

bacteria_busco=/home/rlogares/apps/busco/datasets/bacteria
fungi_busco=/home/rlogares/apps/busco/datasets/fungi
metazoa_busco=/home/rlogares/apps/busco/datasets/metazoa
eukaryota_busco=/home/rlogares/apps/busco/datasets/eukaryota


echo "BUSCO Genome assembly assessment"


for i in $(ls -d *); do cd $i; $python3 $busco -o BUSCO_genome_euk -in masked_contigs/contigs.masked_repeatmasker_trnascan_min1kb.fna -l $eukaryota_busco -m genome ;cd ..; done

### $python3 $busco -o test_busco_genome_euk -in $contigs -l $eukaryota_busco -m genome # -c 12


echo "End of pipeline"
echo ""  
echo "Files should be ready for Augustus | prior trainning with cegma output is needed if not available"
echo ""
echo "R. Logares; 21/02/2016 | Modified 13/07/2017"













