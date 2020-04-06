#!/bin/env perl
while (<>) {
  chomp;
  if (/^>(.+)/) {
    my $new = $1;
    if ($seq) {
      process_sequence($seq, $id, @headers);
      $seq = '';
    }
    ($id, @headers) = split /\s*\|\s*/, $new;
  } else {
    $seq .= $_;
  }
}
process_sequence($seq, $id, @headers) if $seq;
#...
sub process_sequence {
  my($seq, $id, @headers) = @_;

  #do something useful with all the info, like storing in a database
}