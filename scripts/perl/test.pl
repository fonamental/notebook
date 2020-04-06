#!/usr/bin/perl
use strict;
use warnings;
use Bio::DB::EUtilities;
use Data::Dumper;
use Bio::SeqIO;
use Bio::DB::GenBank;
use Bio::Tree::Tree;
use Bio::DB::Taxonomy;

## infile
# txid\tgene_name
my $infile = shift;

open( IN,  '<', $infile ) or die "Can't open $infile\n";
#open( TAXA,  '<', 'taxa_ids.txt' );
#open( GENES, '<', 'gene_names.txt' );
#chomp( my @taxa_ids = (<TAXA>) );
#chomp( my @gene_names = (<GENES>) );

my $db = Bio::DB::GenBank->new();
my $tax_obj = Bio::DB::Taxonomy->new(-source=>'entrez',);
open( my $out, '>', "${infile}.SB.fa" ) || die "Can't open file:$!";	#outfile
while (my $line = <IN>){
  chomp $line;
  my ($taxa_id,$gene_name) = split /\t/ , $line;
    my $search_tag = '[gene]';
    if ($gene_name =~ /RNA/i){
       $search_tag = '';
    }
    my $search = "txid$taxa_id\[organism\] AND $gene_name$search_tag";   #database search input file
    print "SEARCH: $search\n";
    my $eutil = Bio::DB::EUtilities->new(		#object listing DB search parameters 
      -eutil => 'esearch',
      -term   => $search,
      -db     => 'nuccore',
      -retmax => 5,
      -rettype    => 'gb',
      -usehistory => 'y',
      -email      => 'josephorkin@gmail.com'
    );  
    $eutil->get_Response( -file => 'test2.query' ); 	 #test object
    $gene_name =~ s/ NOT .+$//;				
    my @gis = $eutil->get_ids();			
    if (@gis){						#output to STDOUT about search success
      print "SEARCH RESULTS: @gis\n";
    }else{
      print "NOT FOUND: $search\n";
    }
    EACH_GI:
    foreach my $gi (@gis) {				#for each located ID, record following info
      my $seq_obj = $db->get_Seq_by_gi( $gi );		

      my $acc  = $seq_obj->id;
      my $desc = $seq_obj->desc;
      my $seq  = $seq_obj->seq;
      my $binomial = $seq_obj->species->binomial;
      
      my $node = $tax_obj->get_Taxonomy_Node(-gi => $gi);
      my $tree_functions = Bio::Tree::Tree->new();
      my $lineage = $tree_functions->get_lineage_string($node);
	  print "FOUND: $gene_name gi|$gi|$acc  $binomial [$lineage]\n";
      my $print_name = "$acc|$gene_name|$binomial";
	  $print_name =~ s/\s+/_/g;
	  print $out ">$print_name gi|$gi|$acc $gene_name $binomial [$lineage]\n$seq\n";
	  }
}