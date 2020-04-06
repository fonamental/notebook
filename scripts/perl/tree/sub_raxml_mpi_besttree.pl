#!/usr/bin/perl

use strict;
use warnings;

my $file = shift;
chomp $file;
my $short_file = $file;
$short_file =~ s/\.\w+$//;

my $range = 100000000000;
my $minimum = 100000;
my $p = int(rand($range)) + $minimum;

open OUT, ">sub_RAXMLbesttree"."_"."$short_file".".sh";

print OUT "#!/bin/bash\n#PBS -N raxml-mpi\n#PBS -A hii-894-aa\n#PBS -l nodes=16:ppn=8\n#PBS -l walltime=24:00:00\n cd \"\${PBS_O_WORKDIR}\"\n\n";
print OUT "module load compilers/intel/12.0.4\n";
print OUT "module load mpi/openmpi/1.4.5_intel\n\n\n";
print OUT "mpirun -n 128 /home/fburki/bin/raxmlHPC-MPI-SSE3 -e 1.0 -m GTRCATI -f d -N 100 -p $p -s $file -n $short_file\n";
print "$short_file \.\.\. DONE\n";