#!/usr/bin/perl

### This script performs parsing of a blast search, and generates files containing the hit names
### (one file per query). The user provides a number of hits and an e-value threshold.

### USAGE: perl parse2hits.pl file.blastout /path/to/output/folder -nbhits 5 -evalue 1e-25

### file.blastout contains the blast output


use strict;
use warnings;

use Bio::SearchIO;

my $blast_output = $ARGV[0];
my $path = $ARGV[1];
my $nb_hits;
my $evalue_threshold;

#check for the correct argument -nbhits on the command line
if ($ARGV[2] eq '-nbhits') {
	$nb_hits = $ARGV[3];
}
else {
	print "ERROR: unknown $ARGV[2] argument\n";
	exit (-1);
	}

#check for the correct argument -evalue on the command line
if ($ARGV[4] eq '-evalue') {
	$evalue_threshold = $ARGV[5];
}
else {
	print "ERROR: unknown $ARGV[4] argument\n";
	exit (-1);
	}

my $report_obj = Bio::SearchIO -> new (-format => 'blast', 
                     			  	   -file   => $blast_output);
	
#parse the blast output
while (my $result_obj = $report_obj -> next_result) {
	my %hits = ();
	my $query_name = $result_obj -> query_name;
   	
	#do not create a file if no hit was found
	if ($result_obj -> num_hits == 0) {
		next;
  	}
		
	else {
		my $hit_counter = 0;
		while (my $hit_obj = $result_obj -> next_hit and ($hit_counter < $nb_hits)) {
			$hit_counter = $hit_counter+1;
			my $hit_name = $hit_obj -> name;
			my $evalue = $hit_obj -> significance;				
			if ($evalue <= $evalue_threshold) {
				$hits{$hit_name} = $evalue #store hit names and evalue in %hits
			}	
		}			
		if (keys %hits) { #check if %hits is not empty, i.e. has been created
			open (OUT, ">$path$query_name".".hits");
			foreach my $keys (sort {$hits{$a} <=> $hits{$b}} keys(%hits)) {	#sort %hits by evalue			
				print OUT "$keys\t$hits{$keys}\n"; #print to file hit names and evalue
			}			
		}
		
	}
}
print "\n\n";

exit;