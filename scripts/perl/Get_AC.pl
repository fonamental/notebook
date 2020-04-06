#!/usr/bin/perl

use Bio::SearchIO;

use strict;
use warnings;

#USAGE: perl Get_AC.pl blastout 80 [identity cutoff]

###Parse text blast output and get the accession number corresponding to sequences
###above a similarity threshold

my $blast_report = new Bio::SearchIO('-format' => 'blast', '-file' => $ARGV[0]);
my $identity_cutoff = $ARGV[1];
my %hashAC;

open OUT, ">AC.txt";

while (my $result = $blast_report->next_result) {
    while (my $hit = $result->next_hit()) {
		while (my $hsp = $hit->next_hsp()) {
			my $AC = $hit->accession;
			my $identity = $hsp->percent_identity;
			unless (exists $hashAC{$AC}) {
				if ($identity >= $identity_cutoff) {
					print $result->query_name(), "\t";
					print "$AC\t";
					print OUT "$AC\n";
					print "$identity\n";
				}
				else {print "$identity is below threshold\n";}
			}
			$hashAC{$AC} = 1;
		}
	}
}

exit;
