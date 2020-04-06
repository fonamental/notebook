#!/usr/bin/perl

### This script takes a set of parsed blast output files and generate the CREATE and INSERT
### statement needed to populate a sqlite table.

### USAGE: perl CREATE_INSERT_sqlite_table_v3 specielist parsedout
### specieslist is a file containing the species used in the blast searches (one per line). This will be used for the column names in the table.
### parsedout is a file listing the files containing the parsed blast outputs.

use strict;
use warnings;

use Data::Dumper;

#declare the variables
my @species_names = ();
my %blast_parsed = (); #this is a hash of references to anonymous arrays

#SPECIFY FILES ON THE COMMAND LINE: 
#$species_names_file is a file listing the species names used in the blast searches (one name per line)
#$parseout_file is a file listing the parsed-file names (one name per line)
my $species_names_file = shift;
my $parseout_file = shift;

print "\nProvide the name of the table as it will appear in the database:\n";
my $table_name = <STDIN>;
chomp $table_name;

#if more than one db were used in blast searches, the table needs to be larger
print "\nHow many databases did you use to perform the blast searches?\n";
my $number_of_databases = <STDIN>;
chomp $number_of_databases;

#read in the file with the species names, for the column labels in the table
open (IN1, "<$species_names_file") or die ("cannot open $species_names_file");
while ($species_names_file = <IN1>) {
	chomp $species_names_file;
	my @species_names = push @species_names, $species_names_file;
}
close (IN1);

my $size = @species_names;
print "\nFOR INFO: size of the species_names array: $size\n\n";

#check if the number of databases corresponds to the number of species names provided
if ($number_of_databases != $size) {
	print "WRONG: provide the same number of species names as the number of databases used in the blast searches\n";
	exit(1);
} else {};

#create the CREATE statement and print it to a file
open (OUT, ">CREATE_INSERT_table_$table_name.sql");
print OUT "DROP TABLE IF EXISTS $table_name;\n";
print OUT "CREATE TABLE $table_name (\n";
print OUT "\tquery_name        TEXT PRIMARY KEY,\n";
print OUT "\tquery_len         INTEGER NOT NULL,\n";
print OUT "\thit_desc          TEXT NOT NULL,\n";
print OUT "\thit_len           INTEGER,\n";
print OUT "\tscore             INTEGER,\n";
print OUT "\tevalue            REAL,\n";
print OUT "\tidentities        REAL,\n";
print OUT "\tspecies           TEXT NOT NULL\n";
print OUT ");\n\n";

#read in the parseout file (file with all the parsed-file names)
open (IN2, "<$parseout_file") or die ("cannot open $parseout_file");
my $j = 0;
while (my $data_file = <IN2>) {
	chomp $data_file;
	
    # read each file of (parsed) BLAST data and store it in %blast_parsed
	open IN3, "<$data_file" or die "can't open $data_file";	
	while (my $line = <IN3>) {
		chomp $line;
		my ($query_name, $rest) = split "\t", $line, 2;
		my @rest_split = split "\t", $rest;
		$blast_parsed{$query_name} = [@rest_split];
	}
	foreach my $blast_parsed_keys (sort keys %blast_parsed) { #$blast_parsed_keys are the query_names		
		my $species = $species_names[$j];
    	# single-quote all the values
    	my @data = @{$blast_parsed{$blast_parsed_keys}}; #dereference the reference to the anonymous array
    	for (my $i=0; $i < @data; $i++) {
	    	$data[$i] = q{'} . $data[$i] . q{'};
		}	
		# join data into a comma-separated string
		my $data_string = join(',', @data); # join data into a comma-separated string
		print OUT "INSERT INTO \"$table_name\" VALUES (\'$blast_parsed_keys\',$data_string,\'$species\');\n";
	}
	$j = $j + 1;
	close (IN3);
}
close (IN2);

exit;