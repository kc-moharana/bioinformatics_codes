'''
 Filter large FASTA sequence file by given ID list. 
''' 

import sys,os, textwrap
import Bio.SeqIO as BP

try:
 id_list = sys.argv[1]
 fasta_inp = sys.argv[2]
except IndexError:
 print >> sys.stderr, "Usage: python %s ID.txt seq.fasta"%sys.argv[0]
 sys.exit(1)

if not os.path.isfile(fasta_inp+'.idx'):
 fasta_dict = BP.index_db(fasta_inp+'.idx',fasta_inp, 'fasta')
 fasta_dict.close()
fasta_dict = BP.index_db(fasta_inp+'.idx')
with open(id_list) as F:
 for id in F:
  id = id.rstrip()
  print ">%s"%id
  print "\n".join(textwrap.wrap(str(fasta_dict[id].seq), 50))
fasta_dict.close()