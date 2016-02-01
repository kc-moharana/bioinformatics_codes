
## SCRIPT : EXTRACT FLNAKING SEQ USING NYT POSITION
##kcm.eid@gmail.com

=input snp pos file
#ChrNo.	SNP_pos	Reference_ntd
Chr1	769857	C
Chr1	1218603	T
Chr1	6882679	T
Chr1	6897061	G
=cut




use strict;

my $pos_file=$ARGV[0];
my $fasta_file=$ARGV[1];
my $flank_len=$ARGV[2] || 150;


if(!$pos_file || !$fasta_file){ die "ERROR: insuffucient input files.\nUsage: perl $0 pos_file fasta_file [flank_length] > output.txt \npos_file: tab delimited file (chr_no,snp_pos,ref_ntd; comment header with #)[required]\nFasta_file: FASTA seq [required]\nflank_length: length of flanking sequence from side of reference position [default: 150, optional] \n\n";}

my %fasta;

open(F,"<$fasta_file") or die "$! $fasta_file\n";
my $id=0;
while(<F>){
	chomp;
	if(/^>(\S+)/){$id=$1; $fasta{$id}=""; next;}
	else{~s/\s*//g;
	$fasta{$id}.=uc $_;
	}
}
close F;

open(F,"<$pos_file") or die "$! $pos_file\n";
print STDOUT "#ID\tPosition\tRef_NTD\tLeft_flanking_seq\tRight_flanking_seq\n";
while(<F>){
	chomp;
	next if /^#/;
	my ($ch,$pos,$ntd)=split /\s+/,$_;
	if(!exists($fasta{$ch})){ warn"Warning : ID ($ch) not found in fasta file; $ch:$pos ($ntd) SKIPPED.....\n"; next;}
	my $seq=$fasta{$ch};
	my $ntd_f=substr ($seq,$pos-1,1);
	if($ntd_f ne uc $ntd) { warn"$ch : Ref nucleotide ($ntd) not matching ($ntd_f)in fasta file; $ch:$pos ($ntd) SKIPPED.....\n"; next;}
	my $left_flank = substr ($seq,$pos-1-$flank_len,$flank_len);
	my $right_flank = substr ($seq,$pos,$flank_len);

	print STDOUT "$ch\t$pos\t$ntd\t$left_flank\t$right_flank\n";	
}
close F;





