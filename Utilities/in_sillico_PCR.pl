#Wrapper for Bioperl  in-sillico PCR

## Inputs
# - mRNA seq or Genome seq
# - Primer pair list 

use strict;
use Bio::PrimarySeq;
use Bio::Tools::AmpliconSearch;
use Bio::SeqIO;
use Bio::Tools::SeqStats;
use GD::Simple;
use GD::Graph::colour;

#######################################
# tab separated primer pair list, first two columns used
# - return hash of primer pairs
######################################
sub process_primer_pair{
 my $primer_f = shift;
 my %primer_hash=();
 open(F, "<",$primer_f) or die "$! $primer_f\n";
 my $count=0;
 while(<F>){
  chomp;
  next if length($_) <3 or /^#/;
  my @l= split /\t/,$_;
  my ($fwd,$rev) = ($l[0],$l[1]); ## only first two columns 
  $primer_hash{++$count}->{-fwd} = Bio::PrimarySeq->new(-seq =>$fwd);
  $primer_hash{$count}->{-rev} = Bio::PrimarySeq->new(-seq =>$rev);
 }
 close(F);
 return(\%primer_hash);
}


###################################
# Rescale between a, b
###################################
sub _scale_data{
 my ($X,$min_X,$max_X,$a,$b)=@_;
 my $rescaled_X = ((($b-$a)*($X-$min_X))/($max_X-$min_X))+$a;
 return($rescaled_X);
}

###################################
# Mimics Gel-electrophretic bands
# - PCR_output_file
# - primer_hash
# - maximum amplicon_size (bp): max 3000bp
# experiment with ($width_factor,$height_factor) to change image size : scaling factors
###################################
sub draw_gel_bands{
 my $PCR_output_file = shift;
 my $primer_hash =shift;
 my $max_amplicon_size = shift || 500; #bp
 die "maximum amplicon_size  should be less than 3000bp!!!\n" if $max_amplicon_size>3000;
 my @names = GD::Simple->color_names;
 #print "@names\n";
 ## Colors
 my $gel_bg_col='darkgray';
 my $gel_band_col ='white';
 my $primer_label_col='blue';
 my $ladder_label_col='blue';
 
 ## Image height and width 
 my @ladders=(10,20,30,50,100,150,300,500,800,1000,1200,1500,1800,2000,2200,2500,2800,3000);
 my @ladder_bands =();
 foreach my $l(sort({$a<=>$b}@ladders)){
  if($l>$max_amplicon_size){
   last;
  }
  push(@ladder_bands, $l);
 }
 my @wells =('L');
 foreach my $p(sort({$a<=>$b}keys %{$primer_hash})){ push(@wells,'P_'.$p)}
 my ($width_factor,$height_factor) = (800,450);
 my ($img_width,$img_height) = (length(@wells)*$width_factor,length(@ladder_bands)*$height_factor);

 ## Gel image
 my $gel = GD::Simple->new($img_width+10,$img_height+10);
 $gel->bgcolor($gel_bg_col);
 $gel->fgcolor('white');
 $gel->rectangle(10, 10, $img_width,$img_height);

 $gel->fontsize(5);
 $gel->fgcolor($primer_label_col);
 foreach my $w (1..scalar(@wells)){
  print "$w\t_scale_data($w ,1,scalar(@wells),100/5,$img_width-(100/5)-10)),(10)\n";
  $gel->moveTo(_scale_data($w ,1,scalar(@wells),$width_factor/5,$img_width-($width_factor/5)-10),10);
  $gel->string($wells[$w-1]);
 }
 foreach my $l(sort({$b<=>$a}@ladder_bands)){
  my($x1,$y1)=(_scale_data(1,1,scalar(@wells),$width_factor/5,$img_width-($width_factor/5)-10), _scale_data($l,$max_amplicon_size,10,10,$img_height-10) ); ## y1 in reverse order
  $gel->bgcolor($gel_band_col );
  $gel->fgcolor(undef);
  $gel->rectangle($x1-$width_factor*0.05,$y1,$x1+$width_factor*0.05,$y1+5);
  $gel->fontsize(2);
  $gel->fgcolor($primer_label_col);
  $gel->moveTo(30,$y1+10);
  $gel->string($l.'bp');   
 }
 open(PCR, "$PCR_output_file") or die "$! $PCR_output_file \n";
 while(<PCR>){
  chomp;
  next if /^#/;
  my @l=split /\t/,$_;
  next if $#l<3;
  my($x1,$y1)=(_scale_data(int($l[0])+1,1,scalar(@wells),$width_factor/5,$img_width-($width_factor/5)-10), _scale_data(int($l[3]),$max_amplicon_size,10,10,$img_height-10) ); ## y1 in reverse order
  $gel->bgcolor($gel_band_col);
  $gel->fgcolor($gel_band_col);
  $gel->rectangle($x1-$width_factor*0.05,$y1,$x1+$width_factor*0.05,$y1+5);
  $gel->fontsize(0.5);
  $gel->fgcolor($primer_label_col);
  $gel->moveTo($x1-$width_factor*0.1,$y1+15);
  $gel->string($l[4].":".$l[5].'-'.$l[6]);  
 }
 close(PCR); 
 open(IMG, "> $PCR_output_file\_PCR_gel.png") or die "$! $PCR_output_file\_PCR_gel.png\n";
 binmode IMG;
 print IMG $gel->png;
 close IMG;
}


#################################### MAIN SCRIPT ######################################
if($#ARGV<2){
 print STDERR "Insufficient input!!\nUsage: perl $0 <transcript.fa> <primer_pair.txt> <output.txt>\n\t<transcript.fa>:FASTA\n\t<primer_pair.txt>:FWD,REV primer pair tab sep per line\n";
 exit(0)
}

my $primers = process_primer_pair($ARGV[1]);
my $mRNAs = Bio::SeqIO->new(-file => $ARGV[0], -format => 'fasta',); 

#foreach my $i (keys %{$primers}){
 #print "$i\t".$primers->{$i}->{-fwd}->seq()."\t".$primers->{$i}->{-rev}->seq()."\n";
#}
my $output_file = $ARGV[2];
open(OUT, "> $output_file") or die "$! $output_file\n";
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
close(OUT);

## Draw Gel image
draw_gel_bands($output_file,$primers,500);
