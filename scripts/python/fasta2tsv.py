#! /usr/bin/env/ python
"""
# fasta to tsv converter
# usage: fasta_tsv.py input.fasta > output.txt
# output: dump to standard stdout
# Dev: __author__ = 'aung'
# Date: 2013 03 16
"""
import sys
from Bio import SeqIO
fasta_file = sys.argv[1]
for seq_record in SeqIO.parse(fasta_file, "fasta"):
    sys.stdout.write(str(seq_record.id)+'\t'+str(seq_record.seq)+'\n')