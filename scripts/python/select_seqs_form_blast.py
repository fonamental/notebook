import os
import sys
import glob

infile = open(sys.argv[1])
line = infile.read()
infile.close()

print 'fasta read'

seqs = line[1:].split('\n>')

d = {}

for seq in seqs:
	d[seq.split()[0]] = ''.join(seq.split('\n')[1:])

print 'd done'

infile = open(sys.argv[2])
lines = infile.readlines()
infile.close()

out = open('Selected_seqs.fas','w')

for i in lines:
	out.write('>%s\n%s\n' % (i.strip(), d[i.strip()]))

out.close()





