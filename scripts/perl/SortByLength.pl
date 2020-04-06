#!/usr/bin/perl

use strict;
use warnings;

use Bio::SeqIO;

my $input_file = shift;

my $seq_in  = Bio::SeqIO->new( -format => 'fasta', -file => $input_file);

    # loads the whole file into memory - be careful
    # if this is a big file, then this script will
    # use a lot of memory

my $seq;
my @seq_array;
while( $seq = $seq_in->next_seq() ) {
	push(@seq_array,$seq);
}

    # now do something with these. First sort by length,
    # find the average and median lengths and print them out

@seq_array = sort { $b->length <=> $a->length } @seq_array;

my $seq_out = Bio::SeqIO->new(-file => ">$input_file.sort", -format => 'fasta');

# write each entry in the input file to the output file
foreach my $sequences (@seq_array) {
	$seq_out->write_seq($sequences);
}

exit;