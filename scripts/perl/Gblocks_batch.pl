#!/usr/bin/perl

use strict;
use warnings;

my $usage = '

Run Gblocks in a batch mode on a bunch of fasta files

USAGE: perl Gblocks_batch.pl *.linsi

';

die $usage unless @ARGV;

while (my $file = shift @ARGV) {
	
	#count the number of sequences in $file
	my $count = `grep -c \">\" $file`; 
	chomp $count;
	print "\nNumber of sequences in $file: $count\n";
	
	### set up the Gblocks arguments ###
	
	# Minimum Number Of Sequences For A Conserved Position (default: 50%+1)
	my $b1 = int($count/2)+1;
	# Minimum Number Of Sequences For A Flank Position (default: 85% if number of sequence)
	my $b2 = int($count/2)+1;
	# Maximum Number Of Contiguous Nonconserved Positions (default: 8)
	my $b3 = "12";
	# Minimum Length Of A Block (default: 10)
	my $b4 = "4";
	# Allowed Gap Positions (default: n)
	my $b5 = "h";
	# Generic File Extension
	my $e = ".gb";
	
	print "Gblocks $file -b1=$b1 -b2=$b2 -b3=$b3 -b4=$b4 -b5=$b5 -e=$e\n";
	system ("Gblocks $file -b1=$b1 -b2=$b2 -b3=$b3 -b4=$b4 -b5=$b5 -p=s -e=$e");
}

print "\n";

exit 0;
