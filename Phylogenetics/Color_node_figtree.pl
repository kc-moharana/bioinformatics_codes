##Script to crete color scheme for phylogenetic tree in Figtree software
## @Input1: Fasta file used for MSA with IDs informat "species|seqID" format. e.g: Ath|AT4G13040 
## @INPUT2: Color code file with format "Species\tHEXCODE".One code per line and separated by tabs.  e.g Ath  #578e44
## Generate HEXCODE : http://tools.medialab.sciences-po.fr/iwanthue/

use strict;
use File::Basename;


if($#ARGV<1){ die "perl $0 <*.fasta>  <Color_code_file>\n\n";} ## avoid / at he end of folder name 
my $fasta = $ARGV[0] || die 'Fasta file mising !!';
my $color_codes = $ARGV[1] || die 'Color code file missing';

my %style;
open (F, $color_codes ) or die "$! $color_codes ";
while(<F>) {chomp; my($or,$c)=split /\s+/,$_; $style{$or}=$c; }
close F;

my $base = basename($fasta);

open (O,">$base\.style") or die "$! $base\.style";
# my $seq_count = `grep -c ">" $fasta`; chomp($seq_count); ## Linux dependancy 
my $seq_count = seq_count($fasta);
print O "
begin taxa;
        dimensions ntax=$seq_count;
        taxlabels\n";
open (I,$fasta) or die "$! $fasta";
while(<I>){

 if(/^>(\S+)/){
  my($org,$seq)=split /\|/,$1;
   print O "\t\t'$1'[&!color=$style{$org}]\n";
 }
}
 print O ";
end;";
close O; close I;

print STDERR "DONE\nNote: While importing the color scheme file   FigTree(v1.4.2) gives a warning ..
        'Missing cahracter or Data block';
  ";
 #Count no of sequences per Fasta file 
 sub seq_count{
  my $f = shift;
  open (F, $f) or die "$! $f";
  $/ = ">";
  my @no_seq = <F>;
  close F;
  $/ = "\n";
  return sprintf "%d", scalar @no_seq;
 }
