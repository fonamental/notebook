#!/usr/bin/perl 

use strict;
use warnings;

use Bio::SeqIO;

my $usage = '

USAGE: perl RNA2DNA.pl *.fasta

';

die $usage unless @ARGV;

while (my $file = shift @ARGV) {

	my $in = Bio::SeqIO->new(-file => $file, -format => 'fasta');
	my $out = Bio::SeqIO->new(-file => ">$file".".DNA", -format => 'fasta');
	
	while (my $seq_obj = $in->next_seq) {
		my $tmp = $seq_obj->seq;
		$tmp =~ s/T/U/g;
#		print "$tmp\n";
		$seq_obj->seq($tmp);
		$out->write_seq($seq_obj);
	}
}

exit;