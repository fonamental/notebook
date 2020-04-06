#!/usr/bin/perl

use strict;
use warnings;
use Bio::SeqIO;

###USAGE: perl usearch_parser.pl file.fasta file.uc

my $fasta = $ARGV[0];
my $clusters = $ARGV[1];
open IN, "<$clusters";

my %cluster_size;
while (my $line = <IN>) {
	chomp $line;
	if ($line =~ /^C/) {
		my @columns = split ('\t', $line);
		my $size = $columns[2];
		my $seq_id = $columns[8];
		$cluster_size{$seq_id} = $size;
	}
}

my $in  = Bio::SeqIO->new(-file => "$fasta", -format => 'Fasta');
open OUT, ">$fasta".".modif";

while (my $seq = $in->next_seq()) {
	my $id = $seq->id();
	my $sequence = $seq->seq;
	if (exists ($cluster_size{$id})) {
		my $new_id = $id;
		my $value = $cluster_size{$id};
		$new_id =~ s/(\S+)/$1_$value/;
		print OUT "\>$new_id\n$sequence\n";
		print "$new_id\n";
	}
}