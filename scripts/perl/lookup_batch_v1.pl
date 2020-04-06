#!usr/local/bin/perl

### This script indexes a fasta database and generates fasta files from files containing a list of ids.

### USAGE: perl lookup_batch_v1.pl fastafile hitsnames
### fastafile is the name of the indexed file (without .idx).
### hitsnames is a file containing the .hits file names (on name per line)

use strict;
use warnings;
$ENV{BIOPERL_INDEX_TYPE} = "SDBM_File";
use Bio::Index::Fasta;
use Bio::SeqIO;
 
my $fasta = shift or die "no fasta file provided\n";
my $hit_names = shift or die "no hitsnames file provided\n";
my %hashids;

#make index of the fasta db, if already made skip it
my $idx = Bio::Index::Fasta->new(-filename => $fasta . ".idx", -write_flag => 1);
$idx->make_index($fasta);

#read in and go through the file with the .hits names
open IN1, "<$hit_names" or die "cannot open $hit_names";
while (my $hit = <IN1>) {
	chomp $hit;
	my $out = Bio::SeqIO->new (-format => 'Fasta',
                               -file => ">>$hit.fas");
	#read in the .hits files
	open IN2, "<$hit" or die "Can't open $hit : $!";
	print "$hit\n---\n";
	while (my $line = <IN2>) {
		chomp $line;
		#if line QUERY_NAME is found -> do nothing
		if ($line =~/^QUERY_NAME: /) {
			next;
		}
		my ($key, $evalue) = split "\t", $line;
		print "$key\t$evalue\n";
		#fetch and write to file the fasta sequences
		my $seq_obj = $idx->fetch($key);
		$out->write_seq($seq_obj);
	}
	close (IN2);
	print "\n\n";
}
close (IN1);
exit;