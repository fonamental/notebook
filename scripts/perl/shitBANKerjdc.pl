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
    foreach my $gi (@gis) {				#for each located ID, record follwing info
      my $seq_obj = $db->get_Seq_by_gi( $gi );		

      my $acc   = $seq_obj->id;
      my $desc = $seq_obj->desc;
      my $seq  = $seq_obj->seq;
      my $binomial = $seq_obj->species->binomial;
      #my $sub_species = $seq_obj->species->sub_species;
    

      my $node = $tax_obj->get_Taxonomy_Node(-gi => $gi);
      my $tree_functions = Bio::Tree::Tree->new();
      my $lineage = $tree_functions->get_lineage_string($node);




      #if (!defined $sub_species){			#extra loop to add subspecies name space
       #  $sub_species =  "";
      #}else{
       #  $sub_species = " $sub_species";
      #}
      foreach my $feature_obj ( $seq_obj->get_SeqFeatures ) {	#genbank search and retrieval of feature object
        #my $primary_tag = $feature_obj->primary_tag;
        my $tag= "$gene_name";
        #if ($primary_tag eq "CDS"){				#different loop for gene or rRNA product 
         # $tag = "gene";
        #}elsif ($primary_tag eq "rRNA"){
         # $tag = "product";
        #}else{
         # next;
        #}
        #print "$tag";
        my ( $start, $end ) = ( $feature_obj->start, $feature_obj->end ); 	#genbank search for start/end of fasta file
        my $range      = $start . ".." . $end;
        if ($feature_obj->has_tag($tag)){
          my @tag_values = $feature_obj->get_tag_values($tag);
          my $tag_values = join ",", @tag_values;
          next unless $tag_values =~ /$gene_name/i;
          my $len = $end - $start + 1;
       	  my $subseq = $seq_obj->subseq($start, $end);				#download subsequence of total fasta file
          my $Ns = $subseq =~ tr/Nn/Nn/;
          $subseq =~ s/(.{60})/$1\n/g;
          print "FOUND: $gene_name gi|$gi|$acc ($range) LEN:${len}bp Ns:$Ns $binomial [$lineage]\n";
        	my $print_name = "$acc|$gene_name|$binomial";
	  	$print_name =~ s/\s+/_/g;
	print $out ">$print_name ($range) LEN:${len}bp Ns:$Ns gi|$gi|$acc $gene_name $binomial [$lineage]\n$subseq\n";
          #last EACH_GI;							
          }
        else {
          next;
        }
        #last;# last if you want only first hit
      }
     }
}
