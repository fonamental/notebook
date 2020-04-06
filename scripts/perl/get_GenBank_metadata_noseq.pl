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
my ($cln,$org,$mol,$cou,$iso,$env);


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
	($cln,$org,$mol,$cou,$iso,$env)=(undef,undef,undef,undef,undef,undef,undef);

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
  
    }
    
    !$an and $an='n/a';
    !$cln and $cln='n/a';
    !$org and $org='n/a';
    !$mol and $mol='n/a';
    !$cou and $cou='n/a';
    !$iso and $iso='n/a';
    !$env and $env='0';

    print $outfile "$an;$cln;$org;$env;$mol;$cou;$iso\n";
    print "$an;$cln;$org;$env;$mol;$cou;$iso\n";
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
