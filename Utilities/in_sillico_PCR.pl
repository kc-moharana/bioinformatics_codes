#Wrapper for Bioperl  in-sillico PCR

## Inputs
# - mRNA seq or Genome seq
# - Primer pair list 

use strict;
use Bio::PrimarySeq;
use Bio::Tools::AmpliconSearch;
use Bio::SeqIO;
use Bio::Tools::SeqStats;

#######################################
# tab separated primer pair list
# - return hash of primer pairs
######################################
sub process_primer_pair{
 my $primer_f = shift;
 my %primer_hash=();
 open(F, "<",$primer_f) or die "$! $primer_f\n";
 my $count=0;
 while(<F>){
  chomp;
  next if length($_) <3;
  my ($fwd,$rev) =split /\t/,$_;
  $primer_hash{++$count}->{-fwd} = Bio::PrimarySeq->new(-seq =>$fwd);
  $primer_hash{$count}->{-rev} = Bio::PrimarySeq->new(-seq =>$rev);
 }
 close(F);
 return(\%primer_hash);
}
######################################

if($#ARGV<2){
 print STDERR "Insuffificent input!!\nUsage: perl $0 <transcript.fa> <primer_pair.txt> <output.txt>\n\t<transcript.fa>:FASTA\n\t<primer_pair.txt>:FWD,REV primer pair tab sep per line\n";
 exit(0)
}

my $primers = process_primer_pair($ARGV[1]);
my $mRNAs = Bio::SeqIO->new(-file => $ARGV[0], -format => 'fasta',); 

#foreach my $i (keys %{$primers}){
 #print "$i\t".$primers->{$i}->{-fwd}->seq()."\t".$primers->{$i}->{-rev}->seq()."\n";
#}

open(OUT, ">$ARGV[2]") or die "$! $ARGV[2]\n";
my @header= ('#_PRIMER','LEFT_PRIMER_SEQ','RIGHT_PRIMER_SEQ','AMPLICON_SIZE(bp)','TEMPLATE','AMPLICON_START','AMPLICON_END','AMPLICON_SEQ','MOL_Wt._RANGE(KDa.)');
print OUT join("\t",@header)."\n";
while (my $template = $mRNAs->next_seq){
 #print "mRNA:".$template->id()."\n";
  
 foreach my $i (keys %{$primers}){
  my $search = Bio::Tools::AmpliconSearch->new(
   -template   => $template,
   -fwd_primer =>$primers->{$i}->{-fwd},
   -rev_primer =>$primers->{$i}->{-rev},
  );
 
  while (my $amplicon = $search->next_amplicon) {
   #print "Found amplicon at position ".$amplicon->start.'..'.$amplicon->end.":\n";
   #print $amplicon->seq->seq."\n\n";
   my $seq_stats  =  Bio::Tools::SeqStats->new(-seq => $amplicon->seq);
   my $mol_weight = $seq_stats->get_mol_wt(); ##Dalton
   my @output = ($i,$primers->{$i}->{-fwd}->seq(),$primers->{$i}->{-rev}->seq(),length($amplicon->seq->seq),$template->id(),$amplicon->start,$amplicon->end,$amplicon->seq->seq,sprintf("%0.2f",$$mol_weight[0]*0.001).'-'.sprintf("%0.2f",$$mol_weight[1]*0.001));
   print OUT join("\t",@output)."\n";
  } 
 }
} 
close(OUT)
