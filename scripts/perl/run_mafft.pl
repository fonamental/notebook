#!/usr/bin/perl

use strict;
use warnings;

my $usage = '

Run Mafft in a batch mode on a bunch of fasta files

USAGE: perl run_mafft.pl *.fasta

';

die $usage unless @ARGV;

while (my $fasta_file = shift @ARGV) {
	system "mafft $fasta_file > $fasta_file.out";
	print "aligning $fasta_file...\n";
}

exit;
