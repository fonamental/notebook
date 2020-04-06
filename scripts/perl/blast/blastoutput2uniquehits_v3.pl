#!/usr/bin/perl

use strict;

##############################################################
##
## Script to remove redundance based on the hit accession of a blast output
##
## version3 prints out the queries with no hit, version2 keeps only the queries with hits
##
## USAGE:perl blastoutput2uniquehits_v3.pl try.blastout 
##
## use redirection '>' to print the stdout in a file
##
## example: perl blastoutput2uniquehits_v1.pl try.blastout > myhit_table
##
## dated: April 19, 2010


my $infile=$ARGV[0];
chomp $infile;

## open blast output file
open(RD1,$infile) or die "cannot open input file named: $infile\n";

## Take out the first line from blast output file
my $first_line=<RD1>;

## split it to check the correct label-data format
my @lable=split('\t',$first_line);

## Uncomment if you want to print the label 
#print "$lable[0] $lable[1] $lable[2] $lable[3]\n";

## Declare empty array to parse column which need to be printed and to create a uniq list of gene=*
my @BlastOutput=();
my @list=();

## open the file line by line
while(my $line=<RD1>)
{
	chomp $line;
	my @line=split('\t',$line); ## Split the tab separated columns of blast output file
	#print "$line[0]--$line[1]--$line[2]--$line[3]--$line[4]--$line[5]--$line[6]--$line[7]--$line[8]--$line[9]--\n"; ## check the label-data
	my $new_line=$line[0].'-'.$line[1].'-'.$line[3].'-'.$line[5]; ## add new columns if needed -->  here
	push(@BlastOutput,$new_line);
	push(@list,$line[5]); # @list contains the hit_accession
	if ($line[3] eq 'No hits found') {
		print "$line[0]\n";
	}
	else {
		next;
	}
}
close RD1;

## create uniq list of hit accession 
my %hashTemp = map { $_ => 1 } @list; # map -> to lower case
my @list_out = sort keys %hashTemp; ##uniq list goes into @list_out
#print "list_out: @list_out\n\n";
#print scalar(@BlastOutput),"\n";

## loop all the blast output on the uniq list and print the sequence with largest length
for(my $i=0;$i<=$#list_out;$i++)
{
	my @compare=();
	my @rest=();
	#print "$list_out[$i]\n";
	for(my $j=0;$j<=$#BlastOutput;$j++)
	{
		my @l_list=split('-',$BlastOutput[$j]);
		if($l_list[3] ne '') #if hit_accession is not empty
		{	
			if($l_list[3] eq $list_out[$i])
			{
#				print "$list_out[$i] == $l_list[2]\n";
				push(@rest,$l_list[0]);
#				print "@rest\n";
				push(@compare,$l_list[1]); # @compare contains the query_length
#				print "@compare\n";
			}
		}
	}
	my $largest=check_largest(@compare); ## check_largest function call.
	my $index = indexArray("$largest",@compare);
#	print "$largest\t$index\n";
	print "$rest[$index]\n"; ## print the results
}

######################## sub section ################

## sub to find out the largest number from all the number passed through an array
sub check_largest
{
	my @file_data=@_;
	my @sorted=sort{$a<=>$b}@file_data;
	return $sorted[-1];
}


## sub to find the index of an array, given a keyword (here $largest)
sub indexArray($@)
{
 my $s=shift;
 $_ eq $s && return @_ while $_=pop;
 -1;
}
