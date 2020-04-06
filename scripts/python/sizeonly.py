#! /usr/bin/env python

# version 1.1 - cleaned up a bit, used an ordered list
# version 1.0 - first version 


usage="""
SIZEONLY.PY
  Print out the sizes of the sequences in a file
  Warning: currently prints them out alphabetically 
  instead of in the original sequence...
  
  Requires the mybio.py library from 
  http://www.mbari.org/~haddock/scripts/ 

Usage: 
	sizeonly.py Tr100_filename.fasta > lengths.txt
	
"""
import sys
import mybio as mb

if len(sys.argv)<2:
	print usage	

else:
	
	inputfilename = sys.argv[1]
	sequences,names = mb.readfasta_nocheck(inputfilename,False,False,False,True) 
	# don't append sequences, nor remove dupes, not a quality file but do return an ordered list
	
	sys.stderr.write("\nFound %d sequences\n" % len(names))
	for sname in names:
		print len(sequences[sname])



