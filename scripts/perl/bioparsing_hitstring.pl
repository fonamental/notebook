#!/usr/bin/perl

use Bio::SearchIO;

use strict;
use warnings;

#parsing a blast output to get tab-formated infos

my $blast_report = new Bio::SearchIO('-format' => 'blast',
				     '-file' => $ARGV[0]);

#WARNING: set how many hits and hsps you want !

my $nbhits=100;
my $nbhsps=100;

print "Query name\tHit accession\tE-value\tScore\tQuery coverage\thsp length\tIdentities\taligned string\n";

while (my $result = $blast_report->next_result) {
	#output "no hits found" if that's the case
	if ($result -> num_hits == 0) {
		print $result->query_name(), "\t";
#		print $result->query_length(), "\t";
		print "No hits\n";
   	}    	
   	else {
    	my $counter=0;
    	while (my $hit = $result->next_hit()) {
			if ($counter < $nbhits) {
				while (my $hsp = $hit->next_hsp()) {
					if ($counter < $nbhsps) {
						print $result->query_name(), "\t";
#			    		print $result->query_length(), "\t";
#						print $hit->name(), "\t";
						print $hit->accession, "\t";
#						print $hit->description, "\t";
	    				print $hsp->evalue(),"\t";
	    				print $hit->raw_score(),"\t";
	    				print $hit->frac_aligned_query(), "\t";
	    				my $hsp_length = $hsp->hsp_length();
	    				print "$hsp_length\t";
						my $identities = $hsp->frac_identical;
						my $identities_formatted = sprintf ("%.2f", $identities);
						print "$identities_formatted\t";
						print $hsp->hit_string(), "\n";
		    			$counter+=1;
		    		}
				}
    		}
		}
 	}
 }
print "\n";
exit;
