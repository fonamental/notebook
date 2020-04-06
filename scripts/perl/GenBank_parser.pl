#!/usr/bin/perl -w
#Jordi Paps, Oxford 2013
use strict;
no warnings qw{qw};
no warnings ('uninitialized', 'substr');

#Program to parse GenBank files and produce a tabular output, using the following format for seq name:
#Acc_numb_length_U|Eukaryota|Opisthokonta|Metazoa|Phylum|Class|Order|Genus|Genus_species
#AY698072.1.1860_U|Eukaryota|Opisthokonta|Metazoa|Mollusca|Gastropoda|Vetigastropoda|Haliotis|Haliotis diversicolor

#Program works with command-line arguments, it needs a GenBank file and optional value for minimum length of sequence
my($USAGE) = "\nProgram to parse GenBank file, recovering Accession number, Taxonomy and sequence in a fashion similar to PR2; it generates a tabular file. Optionally, it allows to keep sequences of an specific length.\n
WARNING: the program will append new results to old files, you should delete or move outputs from former runs
USAGE:\n$0 GenBank_filename minimum_length \n\n";
unless ($ARGV[0]) {print $USAGE; exit;}

#Check if files exist
my $GB_file = $ARGV[0];
print "genbank source file: $GB_file\n";
my $min_length = $ARGV[1];
print "minimum sequence length: $min_length\n";
unless ( defined $min_length ) {
    $min_length = 0;	#If minimum length of the seq is not defined by the user, then set it to 0 so all seqs are kept
}

#Open files and check them, count records
my @GB_file = open_file ($GB_file);
my $total_records = grep /\/\//, @GB_file;
print "\nWARNING!!!: the program will append new results to old files, you should delete or move outputs from former runs\n";
print "\nNumber of GB records: $total_records\n"; 

#Define some variables for parsing GB-> PR2
my $record_counter = 0;
my $percentage = 0;
my $percentage_old = 0;

my $new_FASTA_header;
my $accession;
my $SSU_check = 0;
my $from = 0;
my $to = 0;
my $length = 0;
my $extracted_seq;
my $in_taxonomy = 0;
my $clade;
my $taxonomy;			#NCBI Taxonomy
my @taxonomy;			#NCBI Taxonomy splitted
#my @parsed_taxonomy = '';	#Taxonomy PR2 style
#my $crap_in_the_middle = 'U|Eukaryota|Opisthokonta|Metazoa';
my $newtaxastring;
my $genus;
my $species;
my $in_features = 0;
my $error = 0;
my @positions = 0;
my $RNA_flag = 0;
my $check_label = 0;
my $get_seq = 0;

my $in_sequence = 0;	
my $sequence;

#Parsing GB-> PR2
print "Percentage done:";
foreach my $line (@GB_file) {
    chomp $line;
    
    #Despite being the first part of the loop, this is the last block executed for each GenBank record
    #Here we save the results of the parsing for a given record, and clean all the variables for the next record.
    if( $line =~ /^\/\// ) { 	# If $line is end-of-record line...
	
	$record_counter++;
	$percentage = int (($record_counter / $total_records) * 100);
	if ($percentage != $percentage_old) {
	print "\n$percentage\%";	    
	}
	$percentage_old = $percentage;

	#Parse taxonomy and produce new FASTA name
	@taxonomy = split (';', $taxonomy);	#Process the taxonomy classification
	#EVAN'S NOTE, SEPT 2016: the old way works only for opisthokonts. we need non opisthokonts too (plants, stramenophiles, etc)
	#do {$clade = shift @taxonomy;
	#    unless ( defined $clade ) {
	#	$clade = "NONE";
	#    }
	#    unless (($clade =~ /Eukaryota/i) || ($clade =~ /Metazoa/i) || ($clade =~ /Ecdysozoa/i) || ($clade =~ /Scalidophora/i) || ($clade =~ /Lophotrochozoa/i) || ($clade =~ /Deuterostomia/i) || ($clade =~ /Xenacoelomorpha/i)) {
	#	push (@parsed_taxonomy, $clade);
	#    } 
	#} until (scalar @parsed_taxonomy == 5);
	
	#SEPT 2016: just keeping the entire taxonomy now
	$newtaxastring = join "|", @taxonomy;
	
	#Dump the sequence, clean it, extract the positions that belong to 18S
	$sequence =~ s/[\s0-9]//g;	#Remove spaces and numbers
	my $length = ( $to - $from ) + 1;
	my $new_FASTA_header = "$accession\_$length\_$newtaxastring\|$genus\|$species";
	
	if ($SSU_check == 1 && ($length > $min_length)) {
	    
	    my $extracted_seq = substr($sequence, $from-1, $length);
	    ##EXCEL database for each phylum
	    #my $outfile_csv1 = "Phylum_$parsed_taxonomy[1]\.csv";
	    #unless ( open( FILE_CSV1, ">>$outfile_csv1") ) {print "Cannot open file \"$outfile_csv1\" to write to!!\n\n"; exit;}
	    #print FILE_CSV1 "$accession,$length,$newtaxastring,$genus,$species,$extracted_seq\n";
	    #close( FILE_CSV1 );
	    
	    #EXCEL database for ALL phyla
	    unless ( open( ALL_PHYLA_CSV, ">>${min_length}.csv") ) {print "Cannot open file \"$min_length.csv\" to write to!!\n\n"; exit;}
	    print ALL_PHYLA_CSV "$accession,$length,$newtaxastring,$genus,$species,$extracted_seq\n";
	    close( ALL_PHYLA_CSV );
	    
	    #FASTA file for each phylum	    
	    my $outfile_fasta = "${min_length}\.fasta";
	    unless ( open(FILE_FASTA, ">>$outfile_fasta") ) {print "Cannot open file \"$outfile_fasta\" to write to!!\n\n"; exit;}
	    print FILE_FASTA "\>$new_FASTA_header\n";
	    print FILE_FASTA "$extracted_seq\n";
	    close( FILE_FASTA );
	    
	} elsif ( $length < $min_length ) {
	    #unless ( open( TOO_SHORT, ">>${min_length}_TOO_SHORT.log") ) {print "Cannot open file \"TOO_SHORT.log\" to write to!!\n\n"; exit;}
	    #print TOO_SHORT "$accession\t$length\n";
	    #close( TOO_SHORT );
	} elsif ( $error == 1 ||undef $to || undef $from || ($to - $from) <= 0 ) {
	    #unless ( open( TO_CHECK, ">>${min_length}_TO_CHECK.log") ) {print "Cannot open file \"TO_CHECK.log\" to write to!!\n\n"; exit;}
	    #print TO_CHECK "$accession\t$length\n";
	    #close( TO_CHECK );
	}
	
	#Save of previous record finished, clean all variables before new record
	$accession = '';
	$SSU_check = 0;
	$length = 0;
	$in_taxonomy= 0;
	$clade = '';
	$taxonomy = '';		#NCBI Taxonomy
	@taxonomy = '';		#NCBI Taxonosplitted
	#@parsed_taxonomy = '';	#Taxonomy PR2 style
	$genus = '';
	$species = '';
	$in_sequence = 0;	
	$sequence = '';
	next; 			# then break out of the FOREACH loop
    
    #Here starts the actual parsing of the GenBank record, starting on the first line
    #We parse using as a guide which terms are found at the beginning of the line and flags
    
    #The next four ELSIF blocks only collect info about the record (Accession, Taxonomy,...)
    } elsif ( $line =~ /^ACCESSION/ ) {		#Parse Accession code
	my @accession_line = split (' ', $line);
	$accession = $accession_line [1];
    } elsif ( ($line =~ /^REFERENCE/) || ($line =~ /^COMMENT/)) {
	$in_taxonomy = 0;	#Flag to indicate taxonomy classification has finished
    } elsif ( $line =~ /^  ORGANISM/ ) {	#Parse Genus and species
	my @species = split (' ', $line);
	$genus = $species [1];
	$species = $species [1].'_'.$species [2];
	$in_taxonomy = 1;	#Flag to indicate next line contains the taxonomy classification
    } elsif ($in_taxonomy) {
	$line =~ s/\s//g;	# Remove spaces from taxonomy classification
	$line =~ s/\.//g;	# Remove dot from taxonomy classification
	$taxonomy .= $line;
    } elsif ( $line =~ /^ORIGIN/) {
	$in_features = 0;	#Flag to indicate "features" are finished
	$in_sequence = 1;	#Flag to indicate we are in the sequence lines
	
    #Next 3 ELSIF blocks check if the line belongs to the FEATURES section of the record,
    #extract the numbers of the start and end positions of the seq, and only keeps the numbers
    #if they belong to an 18S seq
    } elsif ( $line =~ /^FEATURES/) {
	$in_features = 1;	#Flag to indicate "features" are starting
    } elsif ( ($in_features == 1) && (($line =~ /     rRNA            /) || ($line =~ /     misc_RNA        /)) ) {
	#print $line;
	##SEPT 2016, adding a thing to remove a repeat of the accession (usually ACCESSION.1) in the length line
	$line =~ s/$accession\.\d/ /g; #remove any mention of the accession
	$line =~ s/\D/ /g;	#Remove non-numerical characters
	$line =~ s/ +/ /g;	#Remove spaces
	@positions = split (" ", "$line");
	#print STDERR "$positions[0]\t$positions[1]\n"; #SANITY CHECK
	if (scalar @positions < 2) { #there's only one coordinate
		$error = 1; 
		print STDERR "there's only one coordinate for this accession: $accession\n";
	}
	if (scalar @positions > 2) { #there's a join	
		print STDERR "there's a join statement in these coordinates from accession $accession: $line\n";
		print STDERR "taking the first and last base coordinates from the join\n";
		my $last = pop @positions; #just take the first and last elements of the array, there may be a 'join' or something similar
		my $first = shift @positions;
		@positions=(); #clear the array
		$positions[0]=$first; #add the first and last coordinates, this is the sequence stretch to keep
		$positions[1]=$last;
	}
	else {
		#do nothing, @positions has two elements already
	}
	$RNA_flag = 1;		#Flag to indicate seq is RNA
    #Next ELSIF block cheks if former positions belong to 18S, if so keeps them
    } elsif ( (($in_features == 1) && ($RNA_flag == 1)) ) {
	$check_label = check_label ($line);
	if ($check_label == 1) {
	    $SSU_check = 1;	#Flag to indicate seq is SSU
	    $from = $positions[0];
	    $to = $positions[1];
		my $templength = $to - $from;
		#TODO: figure out why length reporting works perfectly here but not in the first block of code
		if ($templength < $min_length) { #SEPT 2016 I'm moving this too short reporting to this section, there's something buggy about what's happening with it up in the first block of code, but i know it works here, and length reporting also works here, so screw it
			print STDERR "18s sequence in $accession from positions $from to $to, length $templength\n";
			unless ( open( TOO_SHORT, ">>${min_length}_TOO_SHORT.log") ) {print "Cannot open file \"TOO_SHORT.log\" to write to!!\n\n"; exit;}
			print TOO_SHORT "$accession\t$templength\n";
			close TOO_SHORT;
		}
		$templength=0;
	}
	$RNA_flag = 0;		#Flag to indicate we are done collecting RNA
    } elsif ($in_sequence == 1) { 	#If we know we're in a sequence...
        $sequence .= $line; 		#...keep loading the $line into $sequence
    }
}

# End of program
print "\n\nTo remove results type:
rm *.csv *.fasta *.log\n";
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

#CHECK_LABEL: subroutine that checks if the labels /product or /note are related to 18S
#The labels have been manually curated from all the labels available in the Genbank records that
#result from searching "txid33208[Organism:exp] (18S OR SSU)"
#Just prove how heterogeneous are gene annotations... and how many typos people do!!!
sub check_label {
    my ($line) = @_;
    my @SSU_labels = (
'                     /gene="18 rRNA"',
'                     /gene="18S"',
'                     /gene="18S rDNA"',
'                     /gene="18S ribosomal RNA"',
'                     /gene="18s rRNA"',
'                     /gene="18 S rRNA"',
'                     /gene="18S rRNA"',
'                     /gene="18s rRNA gene"',
'                     /gene="18S rRNA gene"',
'                     /gene="rRNA"',		#This one may be too generic...
'                     /gene="ssu"',
'                     /gene="SSU"',
'                     /gene="ssu rRNA"',
'                     /gene="SSU rRNA"',
'                     /note="18S rDNA"',
'                     /note="18S ribosomal rDNA"',
'                     /note="18S ribosomal RNA precursor"',
'                     /note="18S ribosomal RNA, putative"',
'                     /note="18S RNA"',
'                     /note="18S rRNA gene"',
'                     /note="18S rRNA; putative"',
'                     /note="3\'end of 18S rRNA"',
'                     /note="containing SSU rDNA, ITS1, 5.8S rDNA, ITS2 and LSU',
'                     /note="contains 18S ribosomal RNA"',
'                     /note="contains 18S ribosomal RNA, and internal',
'                     /note="contains 18S ribosomal RNA and internal transcribed',
'                     /note="contains 18S ribosomal RNA, internal transcribed',
'                     /note="contains 18S ribosomal RNA, ITS1, 5.8S ribosomal',
'                     /note="contains 18S rRNA and ITS1"',
'                     /note="contains 18S rRNA, ITS1"',
'                     /note="contains 18S rRNA, ITS1, 5.8S rRNA"',
'                     /note="contains 18S rRNA, ITS1, 5.8S rRNA and ITS2"',
'                     /note="contains 18S rRNA, ITS1, 5.8S rRNA, ITS2, 28S rRNA"',
'                     /note="contains 18S rRNA, ITS1, 5.8S rRNA, ITS2 and 28S',
'                     /note="contains 18S rRNA, ITS1 and 5.8S rRNA"',
'                     /note="contains partial 18S ribosomal RNA, internal',
'                     /note="may contain 18S ribosomal RNA and internal',
'                     /note="may contain 18S ribosomal RNA, internal transcribed',
'                     /note="obtained using ribosomal RNA primers; may contain',
'                     /note="putative 18S ribosomal RNA"',
'                     /note="ribosomal RNA"',
'                     /note="rRNA primary transcript"',
'                     /note="rRNA primary transcript \(put.\); putative"',
'                     /note="sequence contains 18S rRNA gene"',
'                     /note="sequence contains 18S rRNA gene and ITS1"',
'                     /note="sequence contains 18S rRNA gene, ITS1"',
'                     /note="sequence contains 18S rRNA gene, ITS1, 5.8S rRNA',
'                     /note="sequence contains 18S rRNA gene, ITS1 and 5.8S rRNA',
'                     /note="sequence contains 18S rRNA gene \(partial\), ITS1,',
'                     /note="sequence contains partial 18S rRNA gene and ITS1"',
'                     /note="sequence contains partial 18S rRNA gene, ITS1, 5.8S',
'                     /note="similar to 18S ribosomal RNA"',
'                     /note="small subunit ribosomal RNA"',
'                     /note="SSU"',
'                     /product="18S large subunit ribosomal RNA"',
'                     /product="18S ribosomal RNA"',
'                     /product="18S Ribosomal RNA"',
'                     /product="18S ribosomal RNA gene"',
'                     /product="18S ribosomal RNA type I"',
'                     /product="18S ribosomal RNA type II"',
'                     /product="18S ribosonal RNA"',
'                     /product="18S ribsomal RNA"',
'                     /product="18S small subunit ribosomal"',
'                     /product="18S small subunit ribosomal RNA"',
'                     /product="18S type II ribosomal RNA"',
'                     /product="18S type I ribosomal RNA"',
'                     /product="3\'end of 18S ribosomal RNA"',
'                     /product="contains 18S ribosomal RNA and internal',
'                     /product="contains 18S ribosomal RNA, internal transcribed',
'                     /product="put. 18S ribosomal RNA"',
'                     /product="putative 18S ribosomal RNA"',
'                     /product="putative 28S ribosomal RNA"',
'                     /product="rRNA-18S_rRNA-related"',
'                     /product="small subunit 18S ribosomal RNA"',
'                     /product="small subunit ribosomal RNA"',
'                     /product="small subunit ribosomal RNA gene"',
'                     /product="s-rRNA"',
'                     /product="SSU ribosomal RNA"',
'                     /product="type I 18S ribosomal RNA"',
'                     /product="type II 18S ribosomal RNA"'
		    );
    my $check = 0;
    foreach my $label (@SSU_labels) {
	if ($line =~ m/$label/) {$check = 1}
    }
    return $check;
}
