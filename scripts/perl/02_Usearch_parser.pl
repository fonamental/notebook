#!/usr/bin/perl -w
#Jordi Paps, Oxford 2013
use strict;

#Program to parse Usearch files and produce a tabular output, using the following format for seq name:
#Acc_numb_length_U|Eukaryota|Opisthokonta|Metazoa|Phylum|Class|Order|Genus|Genus_species
#AY698072.1.1860_U|Eukaryota|Opisthokonta|Metazoa|Mollusca|Gastropoda|Vetigastropoda|Haliotis|Haliotis diversicolor

#Program works with command-line arguments, it needs a GenBank file and optional value for minimum length of sequence
my($USAGE) = "\nProgram to parse Usearch files and produce a tabular output, recovering Accession number, Taxonomy and sequence in a fashion similar to PR2; it generates a tabular file.\n
WARNING: the program will append new results to old files, you should delete or move outputs from former runs
USAGE:\n$0 Usearch_output\n\n";
unless ($ARGV[0]) {print $USAGE; exit;}

#Check if files exist
my $FASTA_file = $ARGV[0];

#Open files and check them, count records
my @FASTA_file = open_file ($FASTA_file);

#Define some variables
my @new_file;
my $new_line = '';
my $new_seq = '';

#Parsing Usearch-> CSV
foreach my $line (@FASTA_file){
    if ($line =~ /^#/ || $line =~ /^\s*$/) { next;              #Skip comments and blank lines
    } elsif ($line =~ /^>/) {                                   #FASTA headers
	unless ($new_seq eq ''){                                 #"Unless" will only apply for first seq, that will be empty
	    $new_line .= ",$new_seq\n";
	    push (@new_file, $new_line);
	}
	$new_seq = '';                                          #Clear scalar for next iteration
	$new_line = '';                                          #Clear scalar for next iteration
	$line =~ s/\n$//; $line =~ s/\r$//; $line =~ s/\s/\_/;  #Clean seq name
	$line =~ s/\>//g;
	$line =~ s/\|/,/g;	
	$line =~ s/\_/,/g;
	$new_line .= $line;
    } else {
	$line =~ s/\s//g; $line =~ s/\n$//; $line =~ s/\r$//;   #Remove internal seq breakers
	$new_seq .= "$line";
	}    
}
push (@new_file, $new_line);                                     #Just for the last seq

#Save new FASTA file
my $outfile_fasta = "$ARGV[0]\_Usearch.csv";
unless ( open(FILE_FASTA, ">$outfile_fasta") ) {print "Cannot open file \"$outfile_fasta\" to write to!!\n\n"; exit;}
print FILE_FASTA @new_file;
close( FILE_FASTA );

# End of program
print "
│▒│ /▒/
│▒│/▒/
│▒ /▒/─┬─┐
│▒│▒|▒│▒│
┌┴─┴─┐-┘─┘
│▒┌──┘▒▒▒│
└┐▒▒▒▒▒▒┌┘
└┐▒▒▒▒┌┘
";
print "\n\nFile $ARGV[0] processed!\nGuifré Torruella is lame\n";
print "\nend\n";

################################################################################
#                                Subroutines                                   #
################################################################################

#OPEN_FILE: subroutine to open a file
sub open_file {
    my ($file) = @_;
    my $file_checks = file_checks ($file);
    unless (open (FILE, $file)) {print "Can't open file $file\n" && die;}
    my @file = <FILE>;
    close FILE;
    #print "\nShowing first 10 lines of file $file:\n", @file [0..9],"\n";
    return @file;
}

#FILE_CHECKS: subroutine to check that files exist, are "flat" and are not empty
sub file_checks {
    my ($filename) = @_;
    unless (-e $filename) {print "File $filename does NOT exists\n"; exit;}
    unless (-f $filename) {print "File $filename is a NOT a flat file (it has some kinf of formating)\n"; exit;}
    unless (-s $filename) {print "File $filename is empty\n"; exit;} 
}
