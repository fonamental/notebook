import glob
import os
import sys


seqfile = sys.argv[1]
number_seqs = sys.argv[2]


def split_seqfile(file):
  n = int(number_seqs) * 2
  print n
  os.system('split -l %s %s' % (n, file))


def make_sub_file(c, query, db, outname, outfmt, evalue, max_target_seqs, num_threads):
  out = open('sub_%s.sh' % (c),'w')
  out.write(main)
  out.write('blastn -task blastn -query %s -db /scratch/hii-894-aa/DBS_nr/env_nt -out nt_%s -outfmt %s -evalue %s -max_target_seqs %s -num_threads %s\n' % (query, outname, outfmt, evalue, max_target_seqs, num_threads))
  out.close()


split_seqfile(seqfile)

files = [fname for fname in glob.glob('x*')]

infile = open('main1')
main = infile.read()
infile.close()

for index,file in enumerate(files):
  make_sub_file(index, file, '/scratch/hmr-143-aa/DBS/', '%s_%s' % (file, index), '6', '10', '100', '8')
  os.system('msub sub_%s.sh' % (index))
