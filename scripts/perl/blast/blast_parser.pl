#!/usr/bin/perl

### This script performs parsing of a blast output.

### USAGE: perl parsing_blast.pl blastoutput
### The user is prompted to give the number of hits and hsps to be diplayed after parsing

use strict;
use warnings;

use Bio::SearchIO;

my $blast_output = shift or die "no blast output provided";

print "How many hits: \n";
my $nbhits = <STDIN>;
print "How many hsps per hit: \n";
my $nbhsps = <STDIN>;
print "Provide an evalue threshold:\n";
my $evalue_threshold = <STDIN>;

#read in the output file
open IN, "<$blast_output" or die "cannot open $blast_output";

my $report_object = Bio::SearchIO -> new (-format => 'blast', 
                          				  -file   => $blast_output);
                          				  	  
close (IN);

open (OUT, ">$blast_output"."_parsed.txt");
	
print OUT "query_name\tquery_length\thit_name\thit_length\thit_accession\thit_description\tscore\tevalue\tidentities\n";

#parse the blast output
while (my $result_object = $report_object -> next_result) {
	my $query_name = $result_object -> query_name;
	my $query_length = $result_object -> query_length;
	#my $query_description = $result_object -> query_description;
	
	#output "no hits found" if that's the case
	if ($result_object -> num_hits == 0) {
	print OUT "$query_name\t$query_length\tNo hits found\n";
   	} 
   	
   	else {
		
		my $hit_counter = 0;
		while (my $hit_object = $result_object -> next_hit and ($hit_counter < $nbhits)) {
		$hit_counter = $hit_counter+1;
			my $hit_name = $hit_object -> name;
			my $hit_length = $hit_object -> length;
			my $hit_accession = $hit_object -> accession;
			my $hit_description = $hit_object -> description;
			my $score = $hit_object -> raw_score;
			my $evalue = $hit_object -> significance;
			my $hsp_counter = 0;
			while (my $hsp_object = $hit_object -> next_hsp and ($hsp_counter < $nbhsps)) {
				$hsp_counter = $hsp_counter+1;
				my $identities = $hsp_object -> frac_identical;
				my $identities_formatted = sprintf ("%.2f", $identities);
				if (($evalue < $evalue_threshold) and ($hit_counter == 1)) {
					print OUT "$query_name\t$query_length\t$hit_name\t$hit_length\t$hit_accession\t$hit_description\t$score\t$evalue\t$identities_formatted\n";
				}
				elsif (($evalue < $evalue_threshold) and ($hit_counter > 1)) {
					print OUT "\t\t$hit_name\t$hit_length\t$hit_accession\t$hit_description\t$score\t$evalue\t$identities_formatted\n";
				}
				elsif ($evalue >= $evalue_threshold) {
					last;
				}
			}
		}
	}
}
close (OUT);

exit;