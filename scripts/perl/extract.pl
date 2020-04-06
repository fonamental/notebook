#!/usr/bin/perl

use strict;
use warnings;

my $list = $ARGV[0];
my $table = $ARGV[1];

open $list, "< $list" or die "could not open $list\n";
my $keyRef;
while (<$list>) {
   chomp;
   $keyRef->{$_} = 1;
}
close $list;

open $table, "$table" or die "could not open $list\n";
while (<$table>) {
    chomp;
    my ($testKey, $label, $count) = split("\t", $_);
    if (defined $keyRef->{$testKey}) {
        print STDOUT "$_\n";
    }
}
close $table;