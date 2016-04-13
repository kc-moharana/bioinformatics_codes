


##calculate shanone entropy from a MSA file (extension .clus in clustalW format) for a specific column number
##cautions: nucleotide and amino acids, position number
##input column number and MSA file
# Requires Bio::Tools::dpAlign, in addition to Bio-perl; Only on Unix OS;

##KCM.EID[at]GMAIL.COM

use strict;
use Bio::SeqIO;
use Bio::DB::Fasta;
use Bio::AlignIO;
use Bio::SimpleAlign; 
use Bio::Matrix::IO;
use Bio::Tools::dpAlign;
use Bio::Seq;
use Bio::SeqIO;

my $msa_in=shift;	## MSA file; in fasta
my $pos= shift;		## ref seq;pos
my $type = shift;	## nucl|prot
die "calculate shanone entropy and ntd bases in a site from a MSA file (extension .clus in clustalW format) for a specific column number[101]\nUsage: perl $0 <MSA_in_fasta> <ref_pos> <nucl|prot> > output_file\n" if $#ARGV<2;


my $k=$msa_in;		## refrence seq;
$k=~s/\.clus//g;


warn "calculating for column : $pos\n";
my $column_data; ## hold the coulmn residues
my $in= Bio::AlignIO->new(-file   => $msa_in,-format => "clustalw" );			##see Bio::SimpleAlign - 
my $aln = $in->next_aln();

my $pos_aln=$aln->column_from_residue_number($k, $pos);		## obtain exact position of ref sequence
my %alleles;
foreach my $seq ($aln->each_seq) {				
				my $res = $seq->subseq($pos_aln,$pos_aln); $alleles{$res}++; $column_data.=$res.','; 
				}
#print "orig: $column_data\n";
my $Sj_score=calc_shanon_entropy($column_data);

my $ntds=join (',',keys %alleles );

printf "%s\t%d\t%s\t%f\t%s\n", $k,$pos,$column_data,$Sj_score,$ntds;


sub calc_shanon_entropy
{
	my $col=shift;				## enter a column like e.g C,C,G,C,D,C,G,T
	chop($col);
	$col=~s/-[,]*//g;
	my @residues=split /\,/,$col;
	#foreach  my $i(0..scalar @residues){ splice(@residues,$i,1) if $residues[$i]=~m/\-/g;}	##reomving gaps
	#print "@residues|\n";
	my @amino_acid;
	if($type eq 'nucl'){@amino_acid=('A','T','G','C');}
	elsif($type eq 'prot'){@amino_acid=('A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V','B','Z','X');}
	else{ die "ERROR: Sequence compostion unknown\n";}

	my %aa_freq_hash;
	foreach my $aa(@amino_acid){$aa_freq_hash{$aa}=0;  }
	foreach my $aa(@residues){$aa_freq_hash{$aa}++;  }
	foreach my $aa(keys %aa_freq_hash ){ $aa_freq_hash{$aa}=$aa_freq_hash{$aa}/scalar@residues;    }		##prob calculated
	my $S_j=0;
	foreach my $aa(@amino_acid){
		$S_j+= $aa_freq_hash{$aa} * log2($aa_freq_hash{$aa});
	  }
	$S_j*=(1/log2(4));
#	printf "\nShannon's Normalized Entropy: %0.2f\n",$S_j;
	return $S_j; 		## higher the value highly conserved
}


sub log2 {
	        my $n = shift;
		return 1 if $n==0;
	        return log($n)/log(2);
	    }
	
	

