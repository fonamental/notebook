#! /usr/bin/env python


usage="""
blast.py

run blastx against swissprot and print out a table of results
table fields separated by tabs:
Short match name, e-value, alignment length, query, 
  url for uniprot record, url for fasta of match, matchstring

Requires biopython and blastall installed on your system

-t option after file name searches trembl database

This could take a long time to finish.

steps to process a file of contigs:

wcd -c CT1_all.contigs > CT1_clusters.txt

uniquecluster.py CT1_all.contigs CT1_clusters.txt > CT1_extracted.fta

blasttable.py CT1_extracted.fta > CT1_blasttable.txt

Get a sorted list of unique blasttable entries using:
sort VM1_blasttable.txt | cut -f 1 | uniq -c | sort -n -r | more

To run "headless", use something like
nohup blasttable.py CT1_all.contigs > CT1_blasttable.txt 2>&1 &

This could easily be modified to run with blastx, or to use a different 
 program depending on a run-time parameter.

ver. 2.2 - Accounted for non-matching sequences
ver. 2.1 - Fixed ?error about /tmp folder. Added flag to search trembl database
ver. 2   - rewritten to work with blastplus (blast+)
ver. 1   - 12 Jun 2009

Usage: 
  blasttable.py sequencelist.fta > resultstable.txt
  # to search trembl database (if so-named), append -t
  blasttable.py sequencelist.fta -t > resultstable.txt
"""

import sys
import time
import re

MaxNum=2

from Bio.Blast.Applications import NcbiblastpCommandline
if len(sys.argv)<2:
	print usage	
	# blast_file = "subset.fasta"
else: 

	blast_file = sys.argv[1].strip()
# 	blast_file = "FPsubexcerpt.fta"


	blast_db   = "swissprot"
	if (len(sys.argv) > 2) and (sys.argv[2].strip()=='-t'):
		blast_db   = "trembl"

	# Blastp_Command = NcbiblastpCommandline(query=blast_file, 
	# db=blast_db, evalue=0.001,outfmt=5, num_threads=2)
	# 
	
	ProgramName="blastp"
	
	Blastp_Command = NcbiblastpCommandline(cmd= ProgramName,query=blast_file, 
	db=blast_db, evalue=0.1,outfmt=5, num_threads=3, out="/tmp/blastout.xml", max_target_seqs = MaxNum,\
		num_alignments=MaxNum,num_descriptions=MaxNum)

	startclock=time.asctime()
	starttime= time.time()
	print >> sys.stderr, '# Started BLAST on file',blast_file,'at:  ' + startclock
	print >> sys.stderr, 'COMMAND:', Blastp_Command

	# stdout, stderr = Blastp_Command()
	# This actually takes all the time - searches all the records
	# To get incremental search/parse, you would need to read and 
	# blast each sequence from the file separately.
	
	result_handle = Blastp_Command()
	
	
	from Bio.Blast import NCBIXML
	result_handle = open("/tmp/blastout.xml")

	## print raw info about seq
	debug =  False

	# E_VALUE_THRESH = 10
	recnum = 0
	nomatch = 0
	if result_handle:
		# most normal swissprot record tags look like this
		# this will probably have to be modified for other DB
		titleRE=re.compile("RecName: Full=([^;>]+)")
		# gi|132745|sp|P02385|RL4B_XENLA 60S ribosomal prote
		# from start to first space
		titleBackup=re.compile("^[^ ]+ (.*)")
	
	
		records = NCBIXML.parse(result_handle)
		print >> sys.stderr,"%s\t%s\t%s\t%s\t%s\t%s" % ("Resulttitle", "Evalue", "Align_length","Query","UniprotAccession","Matched") 
		
		for record in records:
			recnum+=1
			#print a progress dot for every 10 records, + every 100
			# if not(recnum%10):
# 				print >> sys.stderr, ".",
# 
# 			if not(recnum%100):
# 				print >> sys.stderr, "+", recnum
# 			if not(recnum%1000):
# 				print >> sys.stderr, "------------------- + ", time.asctime()
			
			# Try because if there is no match, just skip the record
			try:
				d=record.descriptions[0]
				matchT = titleRE.search(d.title)
				if matchT:
					shorttitle= matchT.group(1)[:70]
				else:
					matchB=titleBackup.search(d.title)
					if matchB:
						shorttitle=matchB.group(1)[:70]
					else:	
						shorttitle=d.title[:70]
					
				al=record.alignments[0]
	
		
				# THE REAL ONE LINER
				print "%s\t%s\t%s\t%s\t" % (shorttitle, d.e, al.hsps[0].align_length,record.query) + \
				   "http://www.uniprot.org/uniprot/"+d.accession + "\t" + \
				   al.hsps[0].match.replace(' ','-')
				# This is all only printed out if you are in "debug" mode
				if debug:
					print >> sys.stderr, "##### ",record.query
					print >> sys.stderr, "## q letters ",record.query_letters
					print >> sys.stderr, "## q id ",record.query_id
					print >> sys.stderr, "# DESCR:", d
					# print >> sys.stderr, "<RECORD"
					print >> sys.stderr, "http://www.uniprot.org/uniprot/"+d.accession+".fasta"
					print >> sys.stderr, "expect:",d.e
					print >> sys.stderr, "## ShortTitle: ", shorttitle
					print >> sys.stderr, "len(rec.aln):",len(record.alignments)
					print >> sys.stderr, "al.len: ", al.length, al.hsps[0].align_length
					print >> sys.stderr, "al.hsp.matchstr: ", al.hsps[0].match.replace(' ','-')
				# for al in record.alignments[:3]:
				# 	# for hsp in alignment.hsps[0]:
				# 	hsp = al.hsps[0]
				# 	# print '****Alignment****'
				# 	print "###########",  record.query
				# 	print al.title[:50], al.length, 'bp', len(al.hsps), 'HSPs' 
				# 	print al.title, al.length, hsp.expect
			except:
				sys.stderr.write("###  No match for %s\n" % record.query)
				nomatch+=1
		result_handle.close()
			
	else:
		print >> sys.stderr, "ERROR - no handle"
		
	
	print >> sys.stderr, '# Started BLAST:  ', startclock
	print >> sys.stderr, '# Finished BLAST: ', time.asctime()
	print >> sys.stderr, "# Processed %d records in %.1f minutes. %d had no match" % (recnum, (time.time()-starttime)/60, nomatch)
"""
	###########################################################
	# contents of a record in records:

['__doc__', '__init__', '__module__', 'alignments', 'application', 'blast_cutoff', 
	'database', 'database_length', 'database_letters', 'database_name', 
	'database_sequences', 'date', 'descriptions', 'dropoff_1st_pass', 
	'effective_database_length', 'effective_hsp_length', 'effective_query_length', 
	'effective_search_space', 'effective_search_space_used', 'expect', 'filter', 
	'frameshift', 'gap_penalties', 'gap_trigger', 'gap_x_dropoff', 'gap_x_dropoff_final', 
	'gapped', 'hsps_gapped', 'hsps_no_gap', 'hsps_prelim_gapped', 'hsps_prelim_gapped_attemped', 
	'ka_params', 'ka_params_gap', 'matrix', 'multiple_alignment', 'num_good_extends', 'num_hits', 
	'num_letters_in_database', 'num_seqs_better_e', 'num_sequences', 'num_sequences_in_database', 
	'posted_date', 'query', 'query_id', 'query_length', 'query_letters', 'reference', 'sc_match', 
	'sc_mismatch', 'threshold', 'version', 'window_size']

	# contents of a description in a record
	dir(r.descriptions[0])
	['__doc__', '__init__', '__module__', '__str__', 'accession', 'bits', 'e', 
	'num_alignments', 'score', 'title']

	# contents of an alignment in a record
	dir(r.alignments[0])
	['__doc__', '__init__', '__module__', '__str__', 'accession', 'hit_def', 
	'hit_id', 'hsps', 'length', 'title']

	dir(r.alignments[0].hsps[0])
	['__doc__', '__init__', '__module__', '__str__', 'align_length', 'bits', 'expect', 'frame', 
	'gaps', 'identities', 'match', 'num_alignments', 'positives', 'query', 'query_end', 
	'query_start', 'sbjct', 'sbjct_end', 'sbjct_start', 'score', 'strand']


	"""
