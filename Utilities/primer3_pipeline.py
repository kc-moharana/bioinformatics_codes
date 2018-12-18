'''
#Pipeline for designing primers for qPCR 
 - batch properssing
 - similar to primer-BLAST, but cann't be use non-model organism
 - customized gene annotations (gff3), genomes, transcripts
 - pick primers using primer3
 - runs blast against genome and transcriptome seq.
 - finally reports primers as .tsv.
 - plot primer locations to PDF.
 
# dependancies:
 - Python libs: Biopython, gffutils, dna_features_viewer
 - external programs: blastn, primer3
'''
import sys,os,re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import DNAAlphabet
#from Bio.SeqUtils import GC
#from Bio.SeqUtils import MeltingTemp as mt
import gffutils
import textwrap
import Bio.Emboss.Primer3 as p3
import dna_features_viewer as dna_plot

try:
 gff3 = sys.argv[1]
 genome = sys.argv[2]
 transcripts_seq = sys.argv[3]
 gene_list = sys.argv[4]
 opt_primer_len = 20  # max, min = +-2
 min_five_overlap = 7
 min_three_overlap=4
 min_product_size = 150 #bp
 max_product_size = 250 #bp
 min_GC_perc = 55
 min_Tm = 55
 min_distance_3_primer = 5 ## higher value ensures unique primers[default -1]
 Blastn_extra_params = '-num_threads 10 -word_size 4 -perc_identity 100 -qcov_hsp_perc 100'
 generate_primer3_formatted_output = 1  #[1/0]
 temp_file_loc = 'temp'

except IndexError:
 print >> sys.stderr,"+%s+"%('#'*38)
 print >> sys.stderr, "||%36s||"%('PIPELINE FOR batch Primer designing')
 print >> sys.stderr,"+%s+"%('#'*38)
 print >> sys.stderr,""
 print >> sys.stderr, "Usage:python %s <gene_annotations.gff3> <masked_genome.fa> <transcripts.fa> <selected_gene_list.txt>\n\tRequires blastn,primer3_core on $PATH "%(sys.argv[0])
 print >> sys.stderr,'See example data for format.'
 print >> sys.stderr, ""
 sys.exit(1)


######################################################
#Function to display all params
# - returns: 
######################################################
def print_params(opt_primer_len, min_five_overlap, min_three_overlap, min_product_size, max_product_size,min_GC_perc, min_Tm,min_distance_3_primer,Blastn_extra_params):
 os.system('clear')
 print >> sys.stderr, "%60s"%('#'*60)
 print >> sys.stderr, "# %58s#"%('IMPORTANT PARAMs')
 print >> sys.stderr, "%60s"%('#'*60)
 print >> sys.stderr, "%50s:%d"%('OPTIMAL PRIMER LEN'.upper(),opt_primer_len)
 print >> sys.stderr, "%50s:%d"%('MIN. 5`-overlap'.upper(),min_five_overlap)
 print >> sys.stderr, "%50s:%d"%('MIN. 3`-overlap'.upper(),min_three_overlap)
 print >> sys.stderr, "%50s:%d"%('MIN. product size (bp)'.upper(),min_product_size)
 print >> sys.stderr, "%50s:%d"%('MAX. product size (bp)'.upper(),max_product_size)
 print >> sys.stderr, "%50s:%d"%('MIN. GC (%)'.upper(),min_GC_perc)
 print >> sys.stderr, "%50s:%d"%('MIN. Tm (%)',min_GC_perc)
 print >> sys.stderr, "%50s:%d"%('MIN. 3`-primer disntace',min_distance_3_primer)
 print >> sys.stderr, "%60s"%('#'*60)
 confirm = raw_input('Proceed?[y/n]: ').rstrip()
 if not confirm.upper() == 'Y':
  print >> sys.stderr, "Aborting.."
  sys.exit(1)
 return(1)


######################################################
#Function for parsing primer3 default output file
# - returns: a dict() of all tags
######################################################
def parse_primer3_detailed_output(primer3_output):
 primer3_output_dict=dict()
 f = open(primer3_output, 'r')
 for l in f:
  l = re.split(r'=',l.rstrip())
  primer3_output_dict[l[0]]=l[1]
 f.close()
 return(primer3_output_dict)

######################################################
#Function for parsse tabular blast output
# - returns: a dict() of quryies as key, num of hits as value
######################################################
def parse_blastn_tab_output(blastn_output):   
 blastn_output_dict=dict()
 f = open(blastn_output, 'r')
 for l in f:
  l = re.split("\t",l.rstrip())
  if not l[0] in blastn_output_dict:
   blastn_output_dict[l[0]]=0
  blastn_output_dict[l[0]]+=1
 f.close()
 return(blastn_output_dict)

######################################################
#Function for plotting gene str and primer loc
# - 
######################################################
def draw_primers(exon_array,gene_id):
 features =list()
 st,end,ex_count = 0,0,0
 mRNA_len = 0
 for ex in exon_array:
  ex_count,st,end = ex_count+1,end+1,end+((ex.end-ex.start)+1)  
  mRNA_len = end
  strand=+1
  color="#ffd700" 
  if ex.strand == '-':
   strand=-1
   color='#cffccc'
  features.append(dna_plot.GraphicFeature(start=st, end=end, strand=strand, color=color,label='exon_'+str(ex_count)))
  #if ex_count != len(exon_array):
   #features.append(dna_plot.GraphicFeature(start=end-1, end=end, strand=strand, color='black',label='intron_'+str(ex_count)))
 if os.path.isfile(os.path.join(temp_file_loc, gene_id+'.primers.t.blastn_output.tsv')):
  t_blastn = open(os.path.join(temp_file_loc, gene_id+'.primers.t.blastn_output.tsv'), 'r')
  for l in t_blastn:
   l = re.split(r"\t",l.rstrip())
   strand=+1
   st, end = int(l[8]), int(l[9])
   if end < st:
    st, end = int(l[9]), int(l[8])
    #strand=-1
   features.append(dna_plot.GraphicFeature(start=int(l[8]), end=int(l[9]), strand=strand, color='purple',label=l[0]))
  t_blastn.close()
 else:
  print >> sys.stderr, os.path.join(temp_file_loc,gene_id+'.primers.t.blastn_output.tsv')+' Not found !!'
  
 record = dna_plot.GraphicRecord(sequence_length=mRNA_len, features=features)
 ax, _ = record.plot(figure_width=12)
 ax.figure.savefig(gene_id+'primers_plot.pdf') #bbox_inches='tight'
 print>> sys.stderr, "\tlength %s %d"%(gene_id, mRNA_len)
 return(1)

######################################################
#Main function 
# 
######################################################
def main():
 chk_blastn = not(os.system('blastn -version'))
 chk_primer3 = int(os.system('primer3_core -help')) #  65280
 chk_genome_index = os.path.isfile(genome+'.idx')
 chk_gff_index = os.path.isfile(gff3+'.db.idx')
 if chk_primer3 !=  65280:
  print >> sys.stderr, "primer3_core not found on $PATH!!"
  sys.exit(1)
 if not chk_blastn :
  print >> sys.stderr, "blastn not found on $PATH!!"
  sys.exit(1)
 if not chk_genome_index:
  print >> sys.stderr, "Creating Genome index"
  idx = SeqIO.index_db(genome+'.idx', genome, "fasta")
  idx.close()
  print >> sys.stderr, " --DONE--"
 genome_idx = SeqIO.index_db(genome+'.idx')
 #print genome_idx['Chr01'].seq[1:10]
 if not chk_gff_index:
  print >> sys.stderr, "Creating GFF index.."
  gffutils.create_db(gff3, gff3+'.db.idx')
  print >> sys.stderr, " --DONE--"
 gff_db = gffutils.FeatureDB(gff3+'.db.idx')
 if not os.path.isdir(temp_file_loc):
  os.makedirs(temp_file_loc) 
 if not os.path.isfile(genome+'.nin'):
  cmd = "makeblastdb -in %s -dbtype nucl"%(genome)
  print >> sys.stderr, "Creating genome blast db..\n"+cmd
  os.system(cmd)
 if not os.path.isfile(transcripts_seq+'.nin'):
  cmd = "makeblastdb -in %s -dbtype nucl"%(transcripts_seq)
  print >> sys.stderr, "Creating transcript seq blast db..\n"+cmd
  os.system(cmd)
 print_params(opt_primer_len, min_five_overlap, min_three_overlap, min_product_size, max_product_size,min_GC_perc, min_Tm,min_distance_3_primer,Blastn_extra_params)
 
 target_mRNAs = open(gene_list,'r')
 for l in target_mRNAs:
  l = l.rstrip()
  mRNA = gff_db[l]
  print>> sys.stderr, "Processing: %s %s %d bp"%(mRNA.id, mRNA.strand, (mRNA.end-mRNA.start)+1)
  mRNA_seq = ''
  primer_3_seq = ''
  exon_array = list()
  exon_junctions_list = list()	# list of exon-exon junctions
  exons = sorted([f.id for f in gff_db.children(l, featuretype='exon')])
  if exons[0] !=l+'.exon.1':
   print >> sys.stderr, "Exon ID format mismatched!!. Expecting %s. Found %s"%(l+'.exon.1',exons[0])
  #else:
  # print >> sys.stderr, 'Correct exon format!!'  
  
  for e in xrange(1,len(exons)+1):
   # For each exon segment
   ex = gff_db[l+'.exon.'+str(e)]
   exon_array.append(ex)  ## keepin order constant
   #print ex.id
   s = genome_idx[ex.seqid].seq[ex.start-1:ex.end] # as one 1-indexed gff3, python uses 0-index base
   if mRNA.strand =='-':
    s = s.reverse_complement()
   mRNA_seq +=str(s)
   primer_3_seq += str(s)
   if e!=len(exons):
    primer_3_seq += '-'
    exon_junctions_list.append(len(mRNA_seq))
  ## Store as Bio.Seq() object 
  mRNA_seq =Seq(mRNA_seq,DNAAlphabet())  
  
  ## Generate Primer3_input file
  output = open(os.path.join(temp_file_loc,mRNA.id+'.primer3_input.txt'), 'w')
  input_data= ['SEQUENCE_ID='+mRNA.id,
  'SEQUENCE_TEMPLATE='+str(mRNA_seq),'PRIMER_TASK=pick_pcr_primers',
  'PRIMER_MIN_3_PRIME_OVERLAP_OF_JUNCTION='+str(min_three_overlap),
  'PRIMER_MIN_5_PRIME_OVERLAP_OF_JUNCTION='+str(min_five_overlap),
  'SEQUENCE_OVERLAP_JUNCTION_LIST='+" ".join([str(j) for j in exon_junctions_list]),
  'PRIMER_MIN_THREE_PRIME_DISTANCE='+str(min_distance_3_primer),
  'PRIMER_OPT_SIZE='+str(opt_primer_len),
  'PRIMER_MIN_SIZE='+str(opt_primer_len-2), 
  'PRIMER_MAX_SIZE='+str(opt_primer_len+2),
  'PRIMER_MIN_GC='+str(min_GC_perc),
  'PRIMER_MAX_NS_ACCEPTED=1',
  'PRIMER_PRODUCT_SIZE_RANGE='+str(min_product_size)+'-'+str(min_product_size),
  'P3_FILE_FLAG=0','PRIMER_EXPLAIN_FLAG=1',
  'PRIMER_THERMODYNAMIC_PARAMETERS_PATH=/media/winterfell/kanhu/SOFTWARES/primer3-2.3.7/src/primer3_config/',
  '=']
  print >> output, "\n".join(input_data)  
  output.close()  
  
  print >> sys.stderr, "## RUNNING Primer3 ##"
  if generate_primer3_formatted_output:
   cmd ="primer3_core -format_output -output=%s  < %s"%(os.path.join(temp_file_loc,mRNA.id+'.primer3_Formatted_output.txt'), os.path.join(temp_file_loc, mRNA.id+'.primer3_input.txt'))
   print >> sys.stderr, "\t ",cmd
   os.system(cmd)
  cmd ="primer3_core -output=%s  < %s"%(os.path.join(temp_file_loc,mRNA.id+'.primer3_detailed_output.txt'), os.path.join(temp_file_loc, mRNA.id+'.primer3_input.txt'))
  print >> sys.stderr, "\t ",cmd
  os.system(cmd)  
  
  ## Parse primer3 default output
  pri3_results = parse_primer3_detailed_output(os.path.join(temp_file_loc, mRNA.id+'.primer3_detailed_output.txt'))
  if (not 'PRIMER_PAIR_NUM_RETURNED' in pri3_results) or (int(pri3_results['PRIMER_PAIR_NUM_RETURNED']) == 0 ):
   print >> sys.stderr, "\t No primer pairs found for mRNA: %s "%mRNA.id
   continue
  ### Generate fasta files of primers
  Fas_output = open(os.path.join(temp_file_loc,mRNA.id+'.primer3_output.fas'), 'w')
  for i in xrange(int(pri3_results['PRIMER_PAIR_NUM_RETURNED'])) :
   print >> Fas_output, ">%s\n%s"%('PRIMER_LEFT_'+str(i),pri3_results['PRIMER_LEFT_'+str(i)+'_SEQUENCE'])
   print >> Fas_output, ">%s\n%s"%('PRIMER_RIGHT_'+str(i),pri3_results['PRIMER_RIGHT_'+str(i)+'_SEQUENCE'])
  Fas_output.close()
  print >> sys.stderr, "## RUNNING Blastn Vs genome ##" 
  cmd ="blastn -db %s -query %s -outfmt 6 -out %s %s"%(genome,os.path.join(temp_file_loc,mRNA.id+'.primer3_output.fas'), os.path.join(temp_file_loc, mRNA.id+'.primers.g.blastn_output.tsv'),Blastn_extra_params)
  print >> sys.stderr, "\t ",cmd
  os.system(cmd)
  genome_blastn_out_dict = parse_blastn_tab_output(os.path.join(temp_file_loc, mRNA.id+'.primers.g.blastn_output.tsv'))
  print >> sys.stderr, "## RUNNING Blastn Vs Transcripts ##" 
  cmd ="blastn -db %s -query %s -outfmt 6 -out %s %s"%(transcripts_seq,os.path.join(temp_file_loc, mRNA.id+'.primer3_output.fas'), os.path.join(temp_file_loc, mRNA.id+'.primers.t.blastn_output.tsv'),Blastn_extra_params)
  print >> sys.stderr, "\t ",cmd
  os.system(cmd)
  transcript_blastn_out_dict = parse_blastn_tab_output(os.path.join(temp_file_loc, mRNA.id+'.primers.t.blastn_output.tsv'))

  ## Final output  
  output = open(mRNA.id+'.primer3_output.tsv', 'w')
  print >> output, "mRNA\tprimer_no\tLEFT_PRIMER\tLEFT_GC_PERCENT\tLEFT_TM\tLEFT_HAIRPIN_TH\tLEFT_END_STABILITY\tRIGHT_PRIMER\tRIGHT_GC_PERCENT\tRIGHT_TM\tRIGHT_HAIRPIN_TH\tRIGHT_END_STABILITY\tLEFT_Genome_BLASTN_HITS\tRIGHT_Genome_BLASTN_HITS\tLEFT_transcrpt_BLASTN_HITS\tRIGHT_transcript_BLASTN_HITS"
  for i in xrange(int(pri3_results['PRIMER_PAIR_NUM_RETURNED'])) :
   if not 'PRIMER_LEFT_'+str(i) in genome_blastn_out_dict:
    genome_blastn_out_dict['PRIMER_LEFT_'+str(i)] =0
   if not 'PRIMER_RIGHT_'+str(i) in genome_blastn_out_dict:
    genome_blastn_out_dict['PRIMER_RIGHT_'+str(i)] =0
   if not 'PRIMER_LEFT_'+str(i) in transcript_blastn_out_dict:
    transcript_blastn_out_dict['PRIMER_LEFT_'+str(i)] =0
   if not 'PRIMER_RIGHT_'+str(i) in transcript_blastn_out_dict:
    transcript_blastn_out_dict['PRIMER_RIGHT_'+str(i)] =0
    
   out = [ mRNA.id, str(i) , pri3_results['PRIMER_LEFT_'+str(i)+'_SEQUENCE'], pri3_results['PRIMER_LEFT_'+str(i)+'_GC_PERCENT'], pri3_results['PRIMER_LEFT_'+str(i)+'_TM'], pri3_results['PRIMER_LEFT_'+str(i)+'_HAIRPIN_TH'], pri3_results['PRIMER_LEFT_'+str(i)+'_END_STABILITY'], pri3_results['PRIMER_RIGHT_'+str(i)+'_SEQUENCE'], pri3_results['PRIMER_RIGHT_'+str(i)+'_GC_PERCENT'], pri3_results['PRIMER_RIGHT_'+str(i)+'_TM'], pri3_results['PRIMER_RIGHT_'+str(i)+'_HAIRPIN_TH'], pri3_results['PRIMER_RIGHT_'+str(i)+'_END_STABILITY'],str(genome_blastn_out_dict['PRIMER_LEFT_'+str(i)]),str(genome_blastn_out_dict['PRIMER_RIGHT_'+str(i)]),str(transcript_blastn_out_dict['PRIMER_LEFT_'+str(i)]),str(transcript_blastn_out_dict['PRIMER_RIGHT_'+str(i)]) ]
   print >> output, "\t".join(out)   
  output.close()
  print >> sys.stderr,"## Ploting "
  draw_primers(exon_array, mRNA.id)
 target_mRNAs.close()
 genome_idx.close()
 sys.exit(0)

if __name__ == "__main__":
 main()
 
