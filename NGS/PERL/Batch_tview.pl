
##Script to export read-alignment (samtools tview) using location
# location in this format: chr:start-end
# requires samootls in $PATH;

use strict;


if($#ARGV <2) {die "Batch tview export.\nUsage: perl $0 BAM_file REference_file location_file\n"}
my $bam=shift;
my $ref=shift;
my $locations=shift;


open(F,"<$locations") or die "$!locations\n";

while(<F>){
chomp;
	my ($chr,$t,$pos)=split /\s+/,$_;	##using second pos
	my $cmd="samtools tview -d T -p $chr\:$pos-5 $bam $ref > $chr\_$pos\_tview.txt";
	system($cmd);
}
close F;
