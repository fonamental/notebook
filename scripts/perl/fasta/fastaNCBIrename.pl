#!/usr/bin/perl

use strict;
use warnings;

my $usage = '

NCBIfasta_rename - read a fasta file with one or more sequences retrieved from GenBank and
                   reformat the fasta headers to Genus_species@id
                   
USAGE: perl NCBIfasta_rename.pl file.fasta

';

die $usage unless @ARGV;

while (my $file = shift @ARGV) {
	open OUT, ">$file.rename";
	open (IN, "<$file") or die ("cannot open $file");
	while (my $line = <IN>) {
		chomp $line;
		my @lines = split ('\n', $line);
		foreach my $line (@lines) {
			$line =~ s/^\>gi\|(\d{1,})\|\S+\s.*\[(\S+)\s(\w+)\.{0,1}\s{0,1}\S{0,}\s{0,1}\S{0,}\s{0,1}\S{0,}\s{0,1}\S{0,}\s{0,1}\S{0,}\s{0,1}\S{0,}\]$/>$2_$3\@$1/;
			$line =~ s/^\>gi\|(\d{1,})\|pdb\|.*$/>$1/;
			print OUT "$line\n";
		}
	}
	close OUT;
}
exit;