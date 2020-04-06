#!/usr/bin/perl

use strict;
use warnings;

my $usage = 'USAGE: perl phylomatic.pl fasta_file.fas';

die $usage unless @ARGV;

my $fasta_file;

while (my $fasta_file = shift @ARGV) {
	system "usearch -sortbylength $fasta_file -fastaout $fasta_file.sorted.fas -minseqlength 10";
	print "sorting $fasta_file\n";
	system "usearch -cluster_smallmem $fasta_file.sorted.fas -id 0.97 -centroids $fasta_file.clustered.fas -uc $fasta_file.clusters.uc";
	print "clustering $fasta_file.sorted.fas\n";
	system "cat $fasta_file.clustered.fas $fasta_file.frame.fas > $fasta_file.all.fas";
	system "mafft --thread 4 --reorder --auto $fasta_file.all.fas > $fasta_file.all.pir";
	system "trimal -in $fasta_file.all.pir -out $fasta_file.all.out.pir -gt 0.3 -st 0.001";
	system "Fasttree -nt $fasta_file.all.out.pir > $fasta_file.all.tre";
	}

exit;
