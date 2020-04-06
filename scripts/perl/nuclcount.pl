#!/usr/bin/perl

# Time-stamp: <2004-06-05 21:20:24 kasper>

# TODO: make basebias work.

use warnings;
use strict;
use Bio::Tools::SeqStats;
use Bio::SeqIO;
use Bio::Seq;
use Getopt::Long;
use constant USAGE =><<END;

SYNOPSIS:

 nuclcount.pl [OPTIONS] [infile [outfile]]

DESCRIPTION:

This script counts n-mers in a set of sequences. NOTE: Does not handle
Ns in the entries.

OPTIONS:

       --help
             Prints this help.
       --order <int>
             Order of nucleotide dependency. Eg. order 2 counts
             3-mers. Default is 0.
       --format <format>
             Specifies format. FASTA is default.
       --normalise
             Return frequencies not counts.
       --biasnormalise
             Normalises to take case of base composition bias.
       --relative
             Returns the frequencies relative to the occurrence of the
             reverese complementary word.
       --sortby word|count
             What to sort output by. word is default.


EXAMPLES:

 cat file.fa | nuclcount.pl --order 3 --normalise > file.hist

 nuclcount.pl --format GenBank --normalise file.gb file.tbl

 nuclcount.pl --sortby count --order 4 file.fa file.tbl

AUTHOR:

Kasper Munch

COPYRIGHT:

This program is free software. You may copy and redistribute it under
the same terms as Perl itself.

END

my $help = 0;
my $normalise = 0;
my $order = 0;
my $basebias = 0;
my $relative = 0;
my $format = 'FASTA';
my $sortby = 'word';

GetOptions(
           "help" => \$help,
           "normalise" => \$normalise,
           "basebias" => \$basebias,
           "relative" => \$relative,
           "order=i" => \$order,
           "format=s" => \$format,
           "sortby=s" => \$sortby,
          ) or die USAGE;

$help and die USAGE;

@ARGV = ('-') unless @ARGV;
my $input = shift @ARGV;
@ARGV = ('>&STDOUT') unless @ARGV;
my $output = shift @ARGV;
my $in = Bio::SeqIO->newFh(-file => $input);

$basebias && $order == 0 and die "No point in base composition normaliseing single bases\n";

$sortby =~ /(word)|(count)/ or die "Sort criterion '$sortby' not recognised.\n";

# Do the counting:
my $wordcount = {};
my $basecount = {};
my $total = 0;
my $disregarded = 0;
while (<$in>) {
  my $seq = $_->seq();
  $seq = uc $seq;
  $total += length $seq;
  my $id = $_->id();
  my $count = nucleotide_count($seq, $order);
  while (my ($nuc, $count) = each %$count) {
    $$wordcount{$nuc} += $count;
  }
  # See how many words that has Ns in them:
  while (my ($k, $v) = each %$wordcount) {
    if ($k =~ /N/) {
      $disregarded += $v;
      delete $$wordcount{$k};
    }
  }
  if ($basebias) {
    # Count bases for bias purposes:
    my $count = nucleotide_count($seq, 0);
    while (my ($nuc, $count) = each %$count) {
      $$basecount{$nuc} += $count;
    }
    $$basecount{'N'} ||= 0;
  }
}

# Turn into array:
my @stats;
while (my ($nuc, $count) = each %$wordcount) {
  push @stats, { word => $nuc, count => $count };
}

# Take care of bias if asked for.
if ($basebias) {

  die "Option --basebias  not working just now.\n";

  $$basecount{A} /= $total - $$basecount{N};
  $$basecount{T} /= $total - $$basecount{N};
  $$basecount{C} /= $total - $$basecount{N};
  $$basecount{G} /= $total - $$basecount{N};
  for (my $i=0; $i<@stats; $i++) {
    $stats[$i]{count} /= $total - $order - $disregarded;
    my $word = $stats[$i]{word};
    $word =~ s/(.)/$stats[$i]{count} \/= $$basecount{$1}/eg;
  }
  # Soft limit to make the counts sum to the original total:
  my $newtotal;
  $newtotal += $stats[$_]{count} for (0..$#stats);
  @stats = map { $_->{count} /= $newtotal/$total ; $_ } @stats;
}

if ($relative) {
  for my $i (0..$#stats) {
    $stats[$i]{count} /= $$wordcount{revcompl($stats[$i]{word})};
  }
} elsif ($normalise) {
  @stats = map { $_->{count} /= $total-$order-$disregarded ; $_ } @stats;
}

# Sort:
@stats = sort { $a->{$sortby} cmp $b->{$sortby} } @stats;

for (my $i=0; $i<@stats; $i++) {
  printf "%s\t%.3f\n", $stats[$i]{word}, $stats[$i]{count};
}

sub revcompl {
    my $seq = shift;
    $seq =~ tr/ATGC/TACG/;
    $seq = reverse(split(//, $seq));
    return $seq;
}

sub nucleotide_count {
  my $str = shift;
  my $order = shift;
  my $hash = {};
  if ($order == 0) {
    my $q = substr($str,0,1,"");
    $$hash{$q}+=1;
    $str =~ s/(.)/ $$hash{$1}+=1; $q=$1 /ge;
  } else {
    my $q = substr($str,0,$order,"");
    $str =~ s/(.)/ $$hash{$q.$1}++; $q = substr($q,1).$1 /ge;
  }
  return $hash;
}


=head1 SYNOPSIS:

 nuclcount.pl [OPTIONS] [infile [outfile]]

=head1 DESCRIPTION:

This script counts n-mers in a set of sequences. NOTE: Does not handle
Ns in the entries.

=head1 OPTIONS:

=over 4

=item --help

Prints this help.

=item --order <int>

Order of nucleotide dependency. Eg. order 2 counts
3-mers. Default is 0.

=item --format <format>

Specifies format. FASTA is default.

=item --normalise

Return frequencies not counts.

=item --adjustcounts

Normalises with the base composition.

=back

=head1 EXAMPLES:

 cat file.fa | nuclcount.pl --order 3 --normalise > file.hist

 nuclcount.pl --format GenBank --normalise file.gb file.tbl

 nuclcount.pl --adjustcounts --order 4 file.fa file.tbl

=head1 AUTHOR:

Kasper Munch

=head1 COPYRIGHT:

This program is free software. You may copy and redistribute it under
the same terms as Perl itself.


=cut
