#!/usr/bin/env perl
# protcalc.pl - run w/o any arguments for help (or see below)
# needs a fasta-like file with >name and then the AA seqs
# calculates absorbance factor, pI, etc, and molec wt for AA seq 
# v. 1.5 searches for a particular subsequence
# v. 1.4 which prints a list of the sequences with AAs exceeding a certain threshold 
# Change the values at line 24 and 25 to match the AA or %. 
# v. 1.3 prints percent of each AA.
# 
# usage: printpct.pl FILENAME.fta [printpct( 0/1)?] [motif (e.g. YG)] 
# 1.41 - 2009 - Steven Haddock : haddock [at] mbari [dot] org

$printpct=0;

if (!@ARGV) {
	printintro();
exit;
} else	{
	$searchfor="emptystring";
	$file=$ARGV[0];	
	(-e $file) || die  "Can't open $file: $!\n Run w/o arguments to see help...\n";
	if ($ARGV[1]){ 
		$printpct=$ARGV[1];
		}
	if ($ARGV[2]){ 
		$searchfor=$ARGV[2];
		}
}

# *** CHANGE HERE for flagging specific AA compositions
$AAinput = 'p'.'Y'; # Change the 2nd variable to the name of the AA you want to monitor 
$specialPCT=5.0;    # Change to the threshold percent

$specialmatch=sprintf("\nSequences with %s greater than %.1f percent: \n", $AAinput,$specialPCT);

$totalsite=0; $numotu=0;
$fasta=0;
$ticker=0; 
$foundspecial=0;

# set whether to print a table to AA percentages

$name="";
$foundnames="";
$foundone=0;
$nothing="";$noname="";
$makelut=1;
$macfile=0;
@names=();
@seqs=();
$len=0.4;    
$maxlenname=0;

open(FILE, "< $file") || die  "Can't open $file: $!\n";


while(<FILE>)	{		

	if (!$macfile && /\r/){ # if you find mac cr in line
		print "Macfile detected.\n";
		$/="\r";  # set indicator to mac CR
		$macfile=1;
		seek FILE, 0,0;
		next;
	}
	#			next if /^\s*$/;
	#			next if /^\s*#/;  # comment line skip
	if (/^>/) {
		die "First name super long -- Mac format not detected?" if (length($_)>600);
		#			$name = $_;
		if (/^>DL;/ || /^>P1;\ */){
			print "NBRF format detected.\n" unless ($nbrf || $textout);
			$nbrf=1;
			($name, $nothing) = ($_ =~ /^>\w\w;\ *(.*)(\r|\n)/ ) ;
			#				print "$name\n";
			$noname=<FILE>; #grab line after next
		}
		else { #normal fasta
			($name) = ($_ =~  /^\>(.*)$/);
		}
		#				print "N$name\n";
		$name =~  s/\t/_/g ;
		$name =~  s/(\r|\n)//g ;
		push (@names,$name);
		#				print "name: $name";
		$numotu++;
		if (length($name)>$maxlenname) {
			$maxlenname=length($name);  
#			 print $name.": ".$maxlenname."\n";
		}
	}	
	#need a provision for blank lines before first name			
	elsif ($name ne "") {
		chomp;
		$seq = $_;
		$seq =~ s/\s//g;
		if ($nbrf){
			$seq =~ s/\*//g;
		}
		$seqs{$name} .= $seq;   
	}
} # while file w/in fasta	
close FILE;

	# END BORROWED FROM ENCODENAME
     
$i=1;

foreach $name (@names){
	$seqname=$name;
	$rev = $seqs{$name};
	$rev =~ tr/[a-z]/[A-Z]/;
	$rev =~ s/[^A-Z]//g;

	$_= $rev;                    
	$pad_len = $maxlenname-length($seqname);   
	
	# removed padding for names so don't have to fix 
	# $padname = " " x ($pad_len) . $seqname;
	$padname =  $seqname;
	#sprintf("%-${pad_len}s", $seqname);
 
	$len = length($rev);
	if ($len==0) { 
	print "Length error in seq $seqname";
	$len=0.5;
	};
	
	# Count the AAs
	$nA = s/A/A/g ; #	 71.08	
	$nC = s/C/C/g ; #	 103.14	
	$nD = s/D/D/g ; #	 115.09	
	$nE = s/E/E/g ; #	 129.12	
	$nF = s/F/F/g ; #	 147.18	
	$nG = s/G/G/g ; #	 57.05	
	$nH = s/H/H/g ; #	 137.14	
	$nI = s/I/I/g ; #	 113.16	
	$nK = s/K/K/g ; #	 128.17	
	$nL = s/L/L/g ; #	 113.16	
	$nM = s/M/M/g ; #	 131.19	
	$nN = s/N/N/g ; #	 114.10	
	$nP = s/P/P/g ; #	 97.12	
	$nQ = s/Q/Q/g ; #	 128.13	
	$nR = s/R/R/g ; #	 156.19	
	$nS = s/S/S/g ; #	 87.08	
	$nT = s/T/T/g ; #	 101.11	
	$nV = s/V/V/g ; #	 99.13	
	$nW = s/W/W/g ; #	 186.21	
	$nY = s/Y/Y/g ; #	 163.18        

	$MW = ($nA * 71.07) + ($nR * 156.18) + ($nN * 114.08) + ($nD * 115.08) + ($nC * 103.10) + ($nQ * 128.13) + ($nE * 129.11) + ($nG * 57.05) + ($nH * 137.14) + ($nI * 113.15) + ($nL * 113.15) + ($nK * 128.17) + ($nM * 131.19) + ($nF * 147.17) + ($nP * 97.11) + ($nS * 87.07) + ($nT * 101.10) + ($nW * 186.20) + ($nY * 163.17) + ($nV * 99.13) + 18.02 ;

	#pi
	$pi= ($nA * 6.00 + $nC * 5.02 + $nD * 2.77 + $nE * 3.22 + $nF * 5.48 + $nG * 5.97 + $nH * 7.47 + $nI * 5.94 + $nK * 9.59 + $nL * 5.98 + $nM * 5.74 + $nN * 5.41 + $nP * 6.30 + $nQ * 5.65 + $nR * 11.15+ $nS * 5.68 + $nT * 5.64 + $nV * 5.96 + $nW * 5.89 + $nY * 5.66 )/$len;
# pi from expasy...
	$pix =($nA * 5.57 + $nC * 5.52 + $nD * 4.30 + $nE * 4.60 + $nF * 5.52 + $nG * 5.52 + $nH * 6.74 + $nI * 5.52 + $nK * 8.75 + $nL * 5.52 + $nM * 5.28 + $nN * 5.52 + $nP * 5.96 + $nQ * 5.52 + $nR * 9.75 + $nS * 5.24 + $nT * 5.19 + $nV * 5.49 + $nW * 5.52 + $nY * 5.52)/$len;
	$lenw=($nS + $nT + $nD + $nE + $nR + $nK + $nH + $nY );
	if ($lenw==0) {
		print "Length error in seq $seqname";
		$lenw=0.7;
	};
	
	if (/$searchfor/) {
		if (!$foundone) {
			$foundnames="";
			$foundone=1;
			}
		$foundnames.= $name."\n" ; 
	}
	
	$piw=($nS * 13 + $nT * 13 + $nD * 3.9 + $nE * 4.1 + $nR * 12.5 + $nK * 10.8 + $nH * 6.0 + $nY * 10.1)/$lenw;
    

	# WIKI Values
	$vNH2=8.2 ;                #NH2
	$vCO= 3.65;               #COOH                         
	$vnC= 8.18;       #C    
	$vnD= 3.9 ;        #D   
	$vnE= 4.07;       #$vE    
	$vnH= 6.04;  	#H          
	$vnK= 10.54;      #K     
	$vnR= 12.48;      #R      
	$vnY= 10.46;      #Y  

	$piwi=picalc();


	# EMBOSS		 
	$vNH2= 8.6 ;	$vCO= 3.6 ;	$vnC= 8.5 ;	$vnD= 3.9 ;
	$vnE= 4.1 ;	$vnH= 6.5 ;	$vnK= 10.8; $vnR= 12.5 ;	$vnY= 10.1 ;
	
	$piem=picalc();
	
	# Solomon 
	$vNH2=9.6 ;	$vCO= 2.4 ;	$vnC= 8.3 ;	$vnD= 3.9 ;
	$vnE= 4.3 ;	$vnH= 6.0 ;	$vnK= 10.5 ;	$vnR= 12.5 ;	$vnY= 10.1 ;
	
	$piso=picalc();


# Find percentages
	$pA=$nA/$len * 100;
	$pC=$nC/$len * 100;
	$pD=$nD/$len * 100;
	$pE=$nE/$len * 100;
	$pF=$nF/$len * 100;
	$pG=$nG/$len * 100;
	$pH=$nH/$len * 100;
	$pI=$nI/$len * 100;
	$pK=$nK/$len * 100;
	$pL=$nL/$len * 100;
	$pM=$nM/$len * 100;
	$pN=$nN/$len * 100;
	$pP=$nP/$len * 100;
	$pQ=$nQ/$len * 100;
	$pR=$nR/$len * 100;
	$pS=$nS/$len * 100;
	$pT=$nT/$len * 100;
	$pV=$nV/$len * 100;
	$pW=$nW/$len * 100;
	$pY=$nY/$len * 100;



	# Alternate formula !!
	$ext = (($nW * 5500) +( $nY * 1490) );
	
	#	Note: Cystine is the amino acid formed when of a pair of cysteine molecules are joined by a disulfide bond.
	# don't use C's because can't determine....
#	 $ext = (($nW * 5500) +( $nY * 1490) + ($nC * 125));
	 $abs = $ext / $MW;
	# Molar extinction coefficient
#	$abs = (($nW * 5690) +( $nY * 1280) + ($nC * 120)) / $MW;
	$charge = (($nR + $nK) - ($nE + $nD));

	#                                    
	# print "\n   Original AA: \n    $rev\n";
	if ($i==1){   
		$pad_len = $maxlenname-length("Sequence");   
		if ($pad_len>0) {
			$padder= " " x ($pad_len);
		}   
		else {
			$padder = "";
		}
		#removing padder so it is imported more easily
		print       "Name\tLength\tCharge\tpiwi\tpiem\tpiso\tAbs\tExtcoef\tMolecWt\n";  
		# $tablenum = "Amino Acid occurrence by count\n".$padder  ."Sequence\tA\tC\tD\tE\tF\tG\tH\tI\tK\tL\tM\tN\tP\tQ\tR\tS\tT\tV\tW\tY\n"; 
		# $tablepct = "Amino Acid occurrence by percent\n".$padder."Sequence\tA\tC\tD\tE\tF\tG\tH\tI\tK\tL\tM\tN\tP\tQ\tR\tS\tT\tV\tW\tY\n";
		$tablenum = "Amino Acid occurrence by count\n"  ."Sequence\tA\tC\tD\tE\tF\tG\tH\tI\tK\tL\tM\tN\tP\tQ\tR\tS\tT\tV\tW\tY\n"; 
		$tablepct = "Amino Acid occurrence by percent\n"."Sequence\tA\tC\tD\tE\tF\tG\tH\tI\tK\tL\tM\tN\tP\tQ\tR\tS\tT\tV\tW\tY\n";
   	} 
	printf "%s\t%d\t%d\t%.2f\t%.2f\t%.2f\t%.3f\t%d\t%.2f\n", $padname,$len,$charge,$piwi,$piem,$piso,$abs,$ext,$MW ;
	if ($printpct==1){
		$tablenum .= sprintf("%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
	 		,$padname,$nA, $nC, $nD, $nE, $nF, $nG, $nH, $nI, $nK, $nL, $nM, $nN, $nP, $nQ, $nR, $nS, $nT, $nV, $nW, $nY);
 		$tablepct .= sprintf("%s\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\n",
	 		,$padname,             $pA, $pC,  $pD,  $pE,  $pF,  $pG,  $pH,  $pI,  $pK,  $pL,  $pM,  $pN,  $pP,  $pQ,  $pR,  $pS,  $pT,  $pV,  $pW,  $pY);
	}
	$i++;    
	# print "AA",$$AAinput,"\n";
	
	if ($$AAinput>$specialPCT) {
		$specialmatch.=sprintf("%s : %.1f\n",$padname,$$AAinput ) ; 
		$foundspecial=1;
	}
	
	
} # end foreach sequence

if ($foundspecial) {
	print $specialmatch;
}

if ($printpct==1) {
	print "\n\n";
	print $tablenum;
	print "\n\n";
    print $tablepct;
}

if ($foundone){
	print "\nFound ". $searchfor . " in:\n". $foundnames ."\n";
}

sub picalc {
	# iterative pi calculation modified from http://isoelectric.ovh.org/files/practise-isoelectric-point.html#mozTocId496531
	# needs values for $vCO, etc provided above as globals
	my $keeptrying=1;
	my $piter=6.5;
	my $piterprev = 0.0;         #of finding the solution
	my $piternext = 14.0;        #0-14 is possible $pH range
	my $X = 0.0;
	my $E = 0.01;               # epsilon means precision [pI = $pH ± E]
	my $temp = 0.0;

	while ($keeptrying) {

	    $qN1= -1/ (1+10**($vCO-$piter));               #COOH                         
	    $qN2=-$nD/(1+10**($vnD-$piter));        #D   
	    $qN3=-$nE/(1+10**($vnE-$piter));       #$E    
	    $qN4=-$nC/(1+10**($vnC-$piter));       #C    
	    $qN5=-$nY/(1+10**($vnY-$piter));      #Y  
	    $qP1= $nH/(1+10**($piter-$vnH));  #H          
	    $qP2=   1/(1+10**($piter-$vNH2));                #NH2
	    $qP3= $nK/(1+10**($piter-$vnK));      #K     
	    $qP4= $nR/(1+10**($piter-$vnR));      #R      

	    $NQ=$qN1+$qN2+$qN3+$qN4+$qN5+$qP1+$qP2+$qP3+$qP4;

	   if ($pH>=14.0)
	   {   # you should never see this message
	   	$piter=99; # error
		$keeptrying=0;
	   }     
		else{
	 	  if($NQ<0)              #we are out of range, thus the new $piter value must be smaller    
		    {                   
		       $temp = $piter;
		       $piter = $piter-(($piter-$piterprev)/2);
		       $piternext = $temp;
		   }
		   else                  #we used to small $piter value, so we have to increase it
		   {                     
		       $temp = $piter;
		       $piter = $piter + (($piternext-$piter)/2);
		       $piterprev = $temp;

		   }

		   if (($piter-$piterprev<$E)&&($piternext-$piter<$E)) { # terminal condition, finding isoelectric point with given precision
		       $keeptrying=0; # found it
			}
	    }
	}
	return $piter;
	# Amino acid	NH2	COOH	C	D	E	H	K	R	Y
	# 	EMBOSS		 8.6	3.6	8.5	3.9	4.1	6.5	10.8	12.5	10.1
	# 	DTASelect	 8.0	3.1	8.5	4.4	4.4	6.5	10.0	12.0	10.0
	# 	Solomon		 9.6	2.4	8.3	3.9	4.3	6.0	10.5	12.5	10.1
	# 	Sillero		 8.2	3.2	9.0	4.0	4.5	6.4	10.4	12.0	10.0
	# 	Rodwell		 8.0	3.1	8.33	3.68	4.25	6.0	11.5	11.5	10.07
	# 	Patrickios	11.2	4.2	-	4.2	4.2	-	11.2	11.2	-
	# 	Wikipedia	 8.2	3.65	8.18	3.9	4.07	6.04	10.54	12.48	10.46  
		
}
sub printintro {
		# body...
print <<END_OF_INTRO; 
######################################################
protcalcfta.pl -- performs calculations on protein sequences

Takes an input file of fasta or nbrf AA sequences

v1.3 - 2008 Steven Haddock : haddock [at] mbari [dot] org

1. Calculates molecular weight, length, charge ,pI, ext. coeff.
2. Estimates absorbance coefficients for protein concentration
3. Estimates protein concentration, if given A260 and A280
   where AminoAcidSequence is one line (without spaces)
   and A280 and A260 are for one-cm pathlength.  
4. Prints percent composition of each AA
5. Searches for sequences containing a particular motif
   (Type the search sequence after the file name)

Notes: 
* The absorbance coeff is calculated based on numbers of 
   W,Y & C using:
   A280(1 mg/ml) = (5500 x nW + 1490 x nY + 125 x nC) / MW
The A280 value, divided by this coeff, will estimate concentration.

* If given just A280, normalizes by the absorbance coeff
  and calculates a "Pure Protein" index:
      Conc. (mg/ml) = (A280/abscoeff)

* If given A260, uses a different formula which corrects for DNA:
       Conc. (mg/ml) = (1.55 x A280/abscoeff) - 0.76 x A260)  

* The pi values are determined a new way using the "Wikipedia", EMBOSS, 
  and Solomon values as starting points.

# v. 1.5 searches for a particular subsequence
# v. 1.4 which prints a list of the sequences with AAs exceeding a certain threshold 
# Edit line starting with AAINPUT to indicate the desired AA or %. 

Usage:
          printpct.pl FILENAME.fta [printpct( 0/1)?] [motif (e.g. YG)] 
-----------------------------------------------------------
END_OF_INTRO
}

