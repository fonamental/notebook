#!/opt/local/bin/perl -w

#Quick and dirty perl script to get certain metainformation from genebank accession numbers.

# Mallo, D. 2012

use strict;
use warnings;
use Data::Dumper;
use Bio::DB::GenBank;

#Objects
my $gb=Bio::DB::GenBank->new(-complexity=>1);
my $seq_ob;
my $file;
my $filename;
my @ans;
my ($feat,$val);
my $outfile;

#Numerical values
my $max_it=10;
my $i=0;
my $downloaded=0;

#Strings
my ($cln,$org,$mol,$cou,$iso,$env,$seq);


$gb->retrieval_type('io_string'); #############VERY IMPORTANT!!!!: The default value 'pipeline' uses fork and pipes and these generate huge problems being inside the eval!!!!!

if (scalar(@ARGV)!=1 || !(-f $ARGV[0]))
{
	die "Error, you should run this script as: $0 accession_number_list_file\n";
}

$filename=$ARGV[0];

open ($file,$filename);
@ans=<$file>;
close $file;

open ($outfile,">$filename\_out\.csv");

foreach my $an (@ans)
{
	$i=0;
	$downloaded=0;
	chomp $an;
	print "\nAccesion number to get: $an\n";
	
	#Genbank search (NCBI)
	while($downloaded==0)
	{
		eval
		{
			$seq_ob = $gb->get_Seq_by_acc($an);
		};
		if (!($@)) 
		{
			print "Downloaded correctly\n";
      		$downloaded=1;
   		}
		else
		{
			$seq_ob=undef;
			if ($i==$max_it)
			{
				die "Max number of download attempts reached, net error\n";
			}
			print "Retrying the download, attempt $i of $max_it\n";
			++$i;
			sleep 2;
		}
		
	}
	#print "Toda la secuencia\n";
	#print Dumper $seq_ob;
	($cln,$org,$mol,$cou,$iso,$env,$seq)=(undef,undef,undef,undef,undef,undef,undef,undef);

	FEAT: for $feat($seq_ob->get_SeqFeatures) 
	{
		#print "Nueva feature\n";
		#print Dumper $feat;
		if($feat->has_tag("clone"))  
    	{
        	for $val($feat->get_tag_values("clone")) 
        	{
        		$cln.="$val";
        	}
        }
    	if($feat->has_tag("organism"))  
    	{
        	for $val($feat->get_tag_values("organism")) 
        	{
        		$org.="$val";
        	}
        }

        if($feat->has_tag("mol_type"))  
    	{
        	for my $val($feat->get_tag_values("mol_type")) 
        	{
        		$mol.="$val";
        	}
        }

        if($feat->has_tag("country"))  
    	{
        	for my $val($feat->get_tag_values("country")) 
        	{
        		$cou.="$val";
        	}
        }

        if($feat->has_tag("isolation_source"))  
    	{
        	for my $val($feat->get_tag_values("isolation_source")) 
        	{
        		$iso.="$val";
        	}
        }

        if($feat->has_tag("environmental_sample"))  
    	{
        	for my $val($feat->get_tag_values("environmental_sample")) 
        	{
        		$env.="1";
        	}
        }


        $seq=$seq_ob->seq;
        
    }
    
	#Determine sequence length, start and end nucelotides for folding
    my $length = length($seq);
    
    !$an and $an='n/a';
    !$cln and $cln='n/a';
    !$org and $org='n/a';
    !$mol and $mol='n/a';
    !$cou and $cou='n/a';
    !$iso and $iso='n/a';
    !$env and $env='0';
    !$seq and $seq='n/a';

    print $outfile "$an;$length;$cln;$org;$env;$mol;$cou;$iso;$seq\n";
    print "$an;$length;$cln;$org;$env;$mol;$cou;$iso\n";
	#exit; #DEBUG
}
close $outfile;
exit;
# 
# ACCESSION
# SOURCE
# 	ORGANISM 
# FEATURES
# /organism
# /mol_type
# /country
# /isolation_source
# /environmental_sample
# sequence

################################################################################
#                                Subroutines                                   #
################################################################################

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
