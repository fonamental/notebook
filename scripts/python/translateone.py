#! /usr/bin/env python

	# version 1.79 - Forked from original program to return 1 w/o name of frame
	# Version 1.78 - Return only one ORF and option for fragment only
	# Version 1.77 - Return only --- if length is zero
	# Version 1.76 - Remove gaps before translating
	# Version 1.75 - Return an ordered list of sequences
	# Version 1.7  - Return multiple matches in tie among sequence lengths
	# Version 1.61 - Take user input in addition to file name
	# Version 1.5  - Split on X as well as * for size tests
	# Version 1.41 - Minor bug fix
	# Version 1.4  - Added checkfile option for faster connection
	# Version 1.3  - Added min length of ORF to retain. 
	#                Default is 1. Use 0 to print all
	# Version 1.2  - Reports more of the ambiguity info
	# Version 1.1  - Requires the mybio module for some functions
	# 			   www.mbari.org/staff/haddock/scripts/


usage="""
TRANSLATEDNA.PY - version 1.79

This program takes a fasta file of nucleotide sequences,
does a six-frame translation, and tries to determine the 
longest open reading frame from each sequence. 

If it finds one, it prints that out in fasta format. 

PRINT BEST FRAME 0/[1] 
To print all the ORFs and not just the best, add an argument 
of 0 after the file name. (Adding 1 will print best match.)

PRINT ORF FRAGMENT ONLY [0]/1 (second parameter)
Will just print the longest piece of uninterrupted ORF, 
not the rest of translation. Only works if PrintBestFrame is not zero

To save to a file, put this after the command name:
 	> outputname.txt 

You can also input a single sequence (w/o spaces or returns) 
instead of a file name and it will translate that to the screen

Bug reports  - beroe [at] mac {dot} com

Requires: mybio.py mini-library at www.mbari.org/~haddock/scripts/

Usage: 
	translatedna_one.py filename.fasta [best frame only] [fragment of ORF only] 
	

"""

import re
import sys
import os
# needed for maketrans -- not sure if this is obsolete
import string
import mybio as mb
standardWithAmbig = {
	'AAA':'K', 'AAC':'N', 'AAG':'K', 'AAR':'K', 'AAT':'N', 'AAY':'N', 'ACA':'T', 'ACB':'T', 
	'ACC':'T', 'ACD':'T', 'ACG':'T', 'ACH':'T', 'ACK':'T', 'ACM':'T', 'ACN':'T', 'ACR':'T', 
	'ACS':'T', 'ACT':'T', 'ACV':'T', 'ACW':'T', 'ACY':'T', 'AGA':'R', 'AGC':'S', 'AGG':'R', 
	'AGR':'R', 'AGT':'S', 'AGY':'S', 'ATA':'I', 'ATC':'I', 'ATG':'M', 'ATH':'I', 'ATM':'I', 
	'ATT':'I', 'ATW':'I', 'ATY':'I', 'CAA':'Q', 'CAC':'H', 'CAG':'Q', 'CAR':'Q', 'CAT':'H', 
	'CAY':'H', 'CCA':'P', 'CCB':'P', 'CCC':'P', 'CCD':'P', 'CCG':'P', 'CCH':'P', 'CCK':'P', 
	'CCM':'P', 'CCN':'P', 'CCR':'P', 'CCS':'P', 'CCT':'P', 'CCV':'P', 'CCW':'P', 'CCY':'P', 
	'CGA':'R', 'CGB':'R', 'CGC':'R', 'CGD':'R', 'CGG':'R', 'CGH':'R', 'CGK':'R', 'CGM':'R', 
	'CGN':'R', 'CGR':'R', 'CGS':'R', 'CGT':'R', 'CGV':'R', 'CGW':'R', 'CGY':'R', 'CTA':'L', 
	'CTB':'L', 'CTC':'L', 'CTD':'L', 'CTG':'L', 'CTH':'L', 'CTK':'L', 'CTM':'L', 'CTN':'L', 
	'CTR':'L', 'CTS':'L', 'CTT':'L', 'CTV':'L', 'CTW':'L', 'CTY':'L', 'GAA':'E', 'GAC':'D', 
	'GAG':'E', 'GAR':'E', 'GAT':'D', 'GAY':'D', 'GCA':'A', 'GCB':'A', 'GCC':'A', 'GCD':'A', 
	'GCG':'A', 'GCH':'A', 'GCK':'A', 'GCM':'A', 'GCN':'A', 'GCR':'A', 'GCS':'A', 'GCT':'A', 
	'GCV':'A', 'GCW':'A', 'GCY':'A', 'GGA':'G', 'GGB':'G', 'GGC':'G', 'GGD':'G', 'GGG':'G', 
	'GGH':'G', 'GGK':'G', 'GGM':'G', 'GGN':'G', 'GGR':'G', 'GGS':'G', 'GGT':'G', 'GGV':'G', 
	'GGW':'G', 'GGY':'G', 'GTA':'V', 'GTB':'V', 'GTC':'V', 'GTD':'V', 'GTG':'V', 'GTH':'V', 
	'GTK':'V', 'GTM':'V', 'GTN':'V', 'GTR':'V', 'GTS':'V', 'GTT':'V', 'GTV':'V', 'GTW':'V', 
	'GTY':'V', 'MGA':'R', 'MGG':'R', 'MGR':'R', 'TAA':'*', 'TAC':'Y', 'TAG':'*', 'TAR':'*', 
	'TAT':'Y', 'TAY':'Y', 'TCA':'S', 'TCB':'S', 'TCC':'S', 'TCD':'S', 'TCG':'S', 'TCH':'S', 
	'TCK':'S', 'TCM':'S', 'TCN':'S', 'TCR':'S', 'TCS':'S', 'TCT':'S', 'TCV':'S', 'TCW':'S', 
	'TCY':'S', 'TGA':'*', 'TGC':'C', 'TGG':'W', 'TGT':'C', 'TGY':'C', 'TRA':'*', 'TTA':'L', 
	'TTC':'F', 'TTG':'L', 'TTR':'L', 'TTT':'F', 'TTY':'F', 'YTA':'L', 'YTG':'L', 'YTR':'L',
	'---': '-', '...': '-', '~~~': '-'
}

FrameName = ['f1','f2','f3','r1','r2','r3']

def revcomp(dna):
	#reverse complement of a DNA sequence
	comp = dna.translate(string.maketrans("ATGCatgcRYMKrymkHBVDhbvd", "TACGTACGYRKMYRKMDVBHDVBH"))

	lcomp = list(comp)
	lcomp.reverse()
	return "".join(lcomp)

def dna_translate(cdna, code = standardWithAmbig):
	# translate a cDNA sequence to a protein
	allprot=[]
	allprot.append("".join([ code.get(cdna[i:i+3], "X") for i in xrange(0,len(cdna),3) ]) )
	allprot.append("".join([ code.get(cdna[i:i+3], "X") for i in xrange(1,len(cdna),3) ]) )
	allprot.append("".join([ code.get(cdna[i:i+3], "X") for i in xrange(2,len(cdna),3) ]) )
	revcdna = revcomp(cdna) 
	allprot.append("".join([ code.get(revcdna[i:i+3], "X") for i in xrange(0,len(revcdna),3) ]) )
	allprot.append("".join([ code.get(revcdna[i:i+3], "X") for i in xrange(1,len(revcdna),3) ]) )
	allprot.append("".join([ code.get(revcdna[i:i+3], "X") for i in xrange(2,len(revcdna),3) ]) )
	return allprot

def findlongest(sequencelist,fragonly=False):
	# finds the longest continuous ORF
	
	lenscores=[]
	fragscores=[]
	bestorflist=[]
	for seq in sequencelist:
		orfs = seq.rstrip().replace('X','*').split('*')
		# print orfs
		lenlist= [len(frag) for frag in orfs]
		# print >> sys.stderr, "LENGTHS:", lenlist
		lenscores.append(max(lenlist))
		if fragonly:
			bestorflist.append( orfs[lenlist.index(max(lenlist))] )
		fragscores.append(len(lenlist)-lenlist.count(0))
	# print >>sys.stderr, "FRAGSCORE:", fragscores
	# print >>sys.stderr, "BESTORFLI:", bestorflist
	# print >>sys.stderr, "LENSCORES:", lenscores

	maxlenlist = [i for i,x in enumerate(lenscores) if x==max(lenscores)]
	
#	maxlenind=lenscores.index(max(lenscores)) # max of the stored fragment lengths
	# print maxlenind
	minfraglist = [i for i,x in enumerate(fragscores) if x==min(fragscores)]

#	minfragind=fragscores.index(min(fragscores)) # min of the stored number of fragments
#	print maxlenlist,minfraglist
	# print lenscores, maxlenind
	# print fragscores, minfragind
	# make a master list with all the "winners"

	# Getting rid of minimum number of fragments requirement...
	# masterlist= list(set(minfraglist + maxlenlist))
	masterlist= list(set(maxlenlist))
	
	# print "masterlist",masterlist
	if (len(masterlist)==1):
		if fragonly:
			return [ [bestorflist[ masterlist[0] ]] , masterlist[0] ]
		else:
			return [ [sequencelist[ masterlist[0] ]] , masterlist[0] ]
		
	else:
		if (max(fragscores) > 1):
			# print "Greater XX", minfraglist
			if fragonly:
				return [[bestorflist[x] for x in masterlist],masterlist]# return a count of 1 if just one each mismatch, 
			else:
				return [[sequencelist[x] for x in masterlist],masterlist]
		else: # zero length
			# print "Less XX ", minfraglist
			return [["---"], [0] ]
			# EMD FUNCTION


##########################
# START MAIN PROGRAM
skipcount=0
checkfile=0
printbest=1
blockwidth=71
noread=False
sequences={}
ordernames=[]
orfonly=0
	
if len(sys.argv)<2:
	sys.stderr.write(usage)

else:
	# second parameter for printing all seqs
	# use zero to print all. otherwise it will be a minimum size
	# default is 1 (no size limit)
	if len(sys.argv)>2:
		printbest=int(sys.argv[2])
	if len(sys.argv)>3:
		orfonly=int(sys.argv[3])
		
	inputfilename = sys.argv[1]
	print >> sys.stderr, "FILE:",inputfilename
	
	if not(os.path.exists(inputfilename)):
		sequences[0]=inputfilename.replace('\n','').replace('-','').replace(' ','').replace('\t','') # make a fake name
		noread=True  # Take the input from the command line instead of the file
		sys.stderr.write("\n# No file found, translating input as sequence:\n\n")
		ordernames=['userinput']
	
	else:
		# This function is in the mybio module
		ordernames,sequences = mb.readfasta_list(inputfilename) #
		# print >> sys.stderr, ordernames

#	for name in sequences.keys():
	for ind,name in enumerate(ordernames):
		# print >> sys.stderr, "IND, NAME:", ind, name
		sequence=sequences[ind]
		# print >> sys.stderr, "SEQ: ", sequences
		newseq = re.sub(r'[\n\s]','',sequence).upper()
		newseq = newseq.replace('-','')

		# print newseq

		# returns all 6-frames of translation as a list
		seqlist=dna_translate(newseq)
		
		# print all six frames
		if printbest == 0:

			for frame,myseq in enumerate(seqlist):
				print ">"+name[:20]+ "_"+ FrameName[frame]
				if len(myseq)==0:
					print "----"
				else:
					mb.printblock(myseq)
	
		else:
			# finds the lonest sequence and what frame it's in
			[longorf,frame] = findlongest(seqlist,orfonly)
			# print "longorf,frame",longorf,frame
			if len(longorf) > 1:
				print >> sys.stderr, "## Ambiguous best hits for %s" % (name)
			# 	print >> sys.stderr, longorf
			# 	for orfind  in range(len(longorf)):
			# 		print ">" + name
			# 		mb.printblock(longorf[orfind],blockwidth)
			# 	
			# else:
			print ">" + name
			mb.printblock(longorf[0],blockwidth)
			
			
		# print seqlist[0].split('*')
	
