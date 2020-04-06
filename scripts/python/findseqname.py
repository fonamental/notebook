#! /usr/bin/env python

# Get rid of or keep sequences that match certain criteria

usage="""
FINDSEQNAME.PY

print out fasta of only seqs whose name contains (or doesn't 
contain) a specified sub-string

To save to a file, put this after the command name:
 	> outputname.txt 

<REQUIRES mybio.py from http://www.mbari.org/~haddock/scripts/>


version 1.0 - first version 

Usage: 
	findseqname.py filename.fasta removestring
	findseqname.py filename.fasta -k keepstring
	
"""
import sys
import re
import mybio as mb

blockwidth=71 # optional argument for mb.blockprint
slicecount=0
removing=True
	
if len(sys.argv)<3:
	print usage	

else:
	# second parameter for printing all seqs
	inputfilename = sys.argv[1]

	if len(sys.argv)>3:		# if more than two arguments
		searchstring=sys.argv[3].upper() # no matter what search for arg 3
		if "-k" in sys.argv[2]: # check for -k flag
			removing=False

	else:
		searchstring  = sys.argv[2].upper()

	# This function is in the mybio module
	sequences = mb.readfasta_nocheck(inputfilename)

	"""		Using the XOR function!
		
		Truth Table:
		-------------------------------------------------
		  string in name?  removing?   print seq
					1		1			0			
					0		1			1
					1		0			0
					0		0			1
		-------------------------------------------------
		You can get this using the XOR function ^
"""

	for name in sequences.keys():
		sequence=sequences[name]

		if removing ^ ( searchstring.upper() in name.upper() ):
			# These two lines print the sequence
			print ">" + name 
			# can add a blockwidth parameter, but defaults to 71
			mb.printblock(sequence)
			slicecount+=1
			
	# Using a python equivalent of a ternary operator!!		
	# print ["Keeping","Removing"][removing] , [slicecount,len(sequences.keys())-slicecount][removing], "records out of",len(sequences.keys())

