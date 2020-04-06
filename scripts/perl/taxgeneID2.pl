#!/usr/bin/perl -w
use strict;
use warnings;
use Bio::SeqIO;
use Bio::DB::GenBank;
use Bio::DB::Query::GenBank;
use Bio::DB::Taxonomy;
use Bio::Tree::Tree;
use Data::Dumper;

my $infile = shift;
my $outfile = shift;
open (IN, '<', $infile) or die "Can't open $infile $!\n";

open (OUT, '>', $outfile) or die "cant' open $outfile $!\n";

my $seq_obj;

while (my $line = <IN>) {
    chomp $line;
    my ($specie,$gene) = split ("\t", $line);
    #print "input ($specie,$gene)\n";
    my $query = "txid".$specie."[ORGANISM] AND ".$gene."[GENE]";
    my $query_obj = Bio::DB::Query::GenBank->new(-db=>'nucleotide',-query=>$query);
    my $tax_obj = Bio::DB::Taxonomy->new(-source=>'entrez',);
    my $node = $tax_obj->get_Taxonomy_Node(-taxonid => $specie);
    my $tree_functions = Bio::Tree::Tree->new();
    my $lineage = $tree_functions->get_lineage_string($node);
#print Dumper $query_obj;
   
    my $gb_obj = Bio::DB::GenBank->new;
 
    my $stream_obj = $gb_obj->get_Stream_by_query($query_obj);
 
    while ($seq_obj = $stream_obj->next_seq) {
	foreach my $feature ($seq_obj->top_SeqFeatures) {
	    if ($feature->primary_tag eq 'CDS') {
		my $CDSeq =$feature->seq;
		print OUT ">", $seq_obj->display_id,"_",$lineage,"_",$gene,"_",$seq_obj->length,"\n", $CDSeq->seq, "\n";		
	    }
	    else {
		print OUT ">", $seq_obj->display_id,"_",$lineage,"_",$gene,"_",$seq_obj->length,"\n", $seq_obj->seq, "\n";
	    }
	}
#	print OUT $CDSeq->seq,"\n";
    #do something with the sequence object
#	print OUT ">", $seq_obj->display_id,"_",$lineage,"_",$gene,"_",$seq_obj->length,"\n", $seq_obj->seq, "\n";
    }
}
