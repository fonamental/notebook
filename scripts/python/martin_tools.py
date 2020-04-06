import os
import sys
import glob

'''makes dictionary of fasta file'''
def read_fasta(filename):
  infile = open(filename)
  line = infile.read()
  infile.close()
  seqs = line.split('>')[1:]

  seq_d = {}
  for seq in seqs:
    seq_d[seq.split('\n')[0]] = ''.join(seq.split('\n')[1:])
  
  return seq_d

'''NOT SURE'''
def cut_fasta(infile, outfile, cutoff):
  d = read_fasta(infile)
  out = open(outfile,'w')
  
  for i in d:
    if len(d[i]) >= int(cutoff):
      out.write('>%s\n%s\n' % (i, d[i]))

  out.close()

'''converts fasta file into one per line fasta file'''
def one_lining(infile, outfile):
  d = read_fasta(infile)
  out = open(outfile,'w')
  for i in d:
    out.write('>%s\n%s\n' % (i, d[i]))
  out.close()
  
'''makes blast hit dictionary out of xml formtated blast output'''
def make_blast_dict(file):
  records = NCBIXML.parse(open(file))
  records_d = {}
  for record in records:
    records_d[record.query] = record
  return records_d
  
'''makes dictionary of hit_ids from xml formated blast output 70% lenght'''
def make_hitID_dict(file):
  records = NCBIXML.parse(open(file))
  records_d = {}
  for record in records:
    list = []
    try:
      for i in record.alignments:
        #print 'Query= ', float(len(record.query.replace('-','')))/float(record.query_letters)
        #print 'SBJCT= ', float(len(i.hsps[0].sbjct.replace('-','')))/float(i.length)
        if float(len(i.hsps[0].query.replace('-','')))/float(record.query_letters) > 0.50 and float(len(i.hsps[0].sbjct.replace('-','')))/float(i.length) > 0.50:
          list.append(i.title.split()[0])
    except IndexError:
      pass
    records_d[record.query] = list
  return records_d
  
'''makes dictionary of hit_ids from txt formated blast output'''
def make_hitID_dict_tectblast(file, evalue):
  infile = open(file)
  line = infile.readlines()
  infile.close()
  records = line.split('Query= ')[1:]
  records_d = {}
  for rec in records:
    if rec.count('***** No hits found *****') != 0:
      records_d[rec.split('\n')[0]] = 'no_hit'
    else:
      hits = rec.split('(Bits)  Value\n')[1]
      hit = hits.split('\n')[0]
      records_d[rec.split('\n')[0]] = hit.strip()
  print records_d
  return records_d
  
  
  
  

