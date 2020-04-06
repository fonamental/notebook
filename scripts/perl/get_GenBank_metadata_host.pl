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
my ($hst,$org,$mol,$cou,$iso,$env);


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
	($hst)=(undef);

	FEAT: for $feat($seq_ob->get_SeqFeatures) 
	{
		#print "Nueva feature\n";
		#print Dumper $feat;
		if($feat->has_tag("host"))  
    	{
        	for $val($feat->get_tag_values("host")) 
        	{
        		$hst.="$val";
        	}
        }
    	
  
    }
    
    !$an and $an='n/a';
    !$hst and $hst='n/a';


    print $outfile "$an;$hst\n";
    print "$an;$hst\n";
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
