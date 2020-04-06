#!/usr/bin/perl

use strict;
use warnings;
use Bio::SeqIO;

#USAGE: perl fasta2excel.pl fasta_file > output

my $file = shift;
my $in  = Bio::SeqIO->new(-file => "$file", -format => 'Fasta');
while (my $seq = $in->next_seq()) {
	my $id = $seq->id;
	my $sequence = $seq->seq();
	print "\>$id\t$sequence\n";
}

exit;