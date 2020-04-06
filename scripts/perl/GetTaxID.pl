#!/usr/bin/perl

###USAGE: perl GetTaxID.pl blastout_file [database] (the sequence db (nucleotide or protein) the GI is for)

###Requires an internet connection to access the taxonomy information

use Bio::SearchIO;
use Bio::DB::Taxonomy;

my $blast_report = new Bio::SearchIO('-format' => 'blastxml', # modify this line according to your blast output format (blasttable if tabular; blast if default, etc)
				     				 '-file' => $ARGV[0]);
				     				 
my $db = new Bio::DB::Taxonomy(-source => 'entrez', -verbose => $verbose);

open (OUT, ">taxIDs.txt");

#WARNING: set how many hits and hsps you want !

my $nbhits=1;
my $nbhsps=1;

while (my $result = $blast_report->next_result) {
	#output "no hits found" if that's the case
	if ($result -> num_hits == 0) {
		print OUT $result->query_name(), "\t";
		print $result->query_name(), "\t";
#		print $result->query_length(), "\t";
		print OUT "No hits\n";
		print "No hits\n";
   	}    	
   	else {
    	my $counter=0;
    	while (my $hit = $result->next_hit()) {
			if ($counter < $nbhits) {
				while (my $hsp = $hit->next_hsp()) {
					if ($counter < $nbhsps) {
						print OUT $result->query_name(), "\t";
						print $result->query_name(), "\t";
						my $hitname = $hit->name();
						if ($hitname =~ /gi\|(\d+)\|/) {
							my $gi = $1;
#							print "$gi\t";
							my @nodes= $db->get_Taxonomy_Node(-gi => $gi, -db => $dbname);
							my @array;
							foreach my $node (@nodes) {
								print OUT $node->scientific_name;
								print $node->scientific_name;
								my $parent = $db->get_Taxonomy_Node($node->parent_id);
#								print $parent->node_name, "\t";
								while (defined $parent && $parent->node_name ne 'root') { 									
									if ($parent->rank eq 'superkingdom') {
										push @array, $parent->node_name;
									}
									elsif ($parent->rank eq 'kingdom') {
										push  @array, $parent->node_name;
									}
									elsif ($parent->rank eq 'phylum') {
										push @array, $parent->node_name;
									}							
									$parent = $db->get_Taxonomy_Node($parent->parent_id);
								}
							}
							my $size = @array;
							if ($size == 1) {
								print OUT "\t\t\t$array[0]\t";
								print "\t\t\t$array[0]\t";
							}
							elsif ($size == 2) {
								print OUT "\t\t$array[0]\t$array[1]\t";
								print "\t\t$array[0]\t$array[1]\t";
							}
							elsif ($size == 3) {
								print OUT "\t$array[0]\t$array[1]\t$array[2]\t";
								print "\t$array[0]\t$array[1]\t$array[2]\t";
							}
						}
					print OUT $hsp->evalue(),"\n";
	    			print $hsp->evalue(),"\n";
		    		$counter+=1;	
					}
	    		}
			}
   		}
	}
 }

print "\n";
exit;