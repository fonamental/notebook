#! /usr/bin/env python

# Get rid of or keep sequences that match certain criteria

usage="""
FINDMOTIF.PY

print out fasta of only seqs whose names contain 
(or don't contain) a sub-string

To save to a file, put this after the command name:
 	> outputname.txt 

<REQUIRES mybio.py from http://www.mbari.org/~haddock/scripts/>


versions 
1.1 - you can enter a regexp as the string?
1.0 - first version 

Usage: 
	findmotif.py filename.fasta removestring
	findmotif.py filename.fasta -k keepstring
	
"""
import sys
import re
import mybio as mb
import time

blockwidth=71 # optional argument for mb.blockprint
slicecount=0
removing=True
	
if len(sys.argv)<3:
	print usage	

else:
	# second parameter for printing all seqs
	inputfilename = sys.argv[1]

	if len(sys.argv)>3:		# if more than two arguments
		searchstring=sys.argv[3] # no matter what search for arg 3
		if "-k" in sys.argv[2]: # check for -k flag
			removing=False

	else:
		searchstring  = sys.argv[2]
	print >> sys.stderr, "Reading sequences... ", time.asctime()
	# This function is in the mybio module
	sequences = mb.readfasta_nocheck(inputfilename)

	"""		Using the XOR function!
		
		Truth Table:
		-------------------------------------------------
		  string in name?  removing?   print seq
					1		1			0			
					0		1			1
					1		0			1
					0		0			0
		-------------------------------------------------
		You can get this using the XOR function ^
"""
	print >> sys.stderr, "Searching... ", time.asctime()
	for name in sequences.keys():
		sequence=sequences[name]

		if removing ^ ( not re.search(searchstring,sequence) == None ):
#		if removing ^ ( searchstring.upper() in sequence ):
			# These two lines print the sequence
			print ">" + name 
			# can add a blockwidth parameter, but defaults to 71
			mb.printblock(sequence)
			slicecount+=1
			
	# Using a python equivalent of a ternary operator!!		
	print >> sys.stderr, ["Keeping","Removing"][removing] , [slicecount,len(sequences.keys())-slicecount][removing], "records out of",len(sequences.keys())
	print >> sys.stderr, time.asctime()
