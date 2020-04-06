#!usr/bin/perl 
use warning;
use sctric;

my $sam_file = shift;
   open (my $in_file, '<', $sam_file) 
   or die "can't open file $sam_file $!\n";

my $db_fasta = shift;


my %hash;
while (my $line = <$in_file>){
   chomp $line;
   #my ($key, $other, $value) = split 
   my @line = split /\t/, $line;
   $hash{$line[2]}++;
}
print "gene\tcount\tdesc\n";
foreach my $key (sort { $hash{$b} <=> $hash{$a} }  keys %hash){
   my $value = $hash{$key};

   my $header = `grep "$key" $db_fasta`;
   $header =~ s/^.+Ns:\d+ //;
  
   print "$key\t$value\t$header\n";
}


