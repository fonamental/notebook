#!/usr/bin/perl

##USAGE: perl FigTree2seq.pl nexus_file original_fasta_file
##After colouring sequences in FigTree, this script removes the corresponding coloured sequences from a fasta file

use strict;
use warnings;
use Bio::SeqIO;

my $nexus = $ARGV[0];
chomp $nexus;
my $short_nexus = $nexus;
$short_nexus =~ s/\.\w+$//;
my $fasta = $ARGV[1];

open IN, "<$nexus";
my $in  = Bio::SeqIO->new(-file => "$fasta", -format => 'Fasta');
my $out = Bio::SeqIO->new(-file => ">$short_nexus".'_only.fasta', -format => 'Fasta');

my $all_lines = do {local $/; <IN>};
my ($a, $b) = split ('taxlabels', $all_lines);
my ($c, $d) = split ("\n\;", $b);
$c =~ s/\n//g;
$c =~ s/\'//g;
my @ids = split ("\t", $c);

my %hashids;
foreach my $id (@ids) {
	if ($id =~ /(\S+)\[\&\!color\=\#\-(\d+)/) {
		next;
	}
	else {
		$hashids{$id} = 1;
	}
}

while (my $seq = $in->next_seq()) {
	my $id = $seq->id();
	if ($hashids{$id}) {
		print "$id\n";
    	$out->write_seq($seq); 
	}
}