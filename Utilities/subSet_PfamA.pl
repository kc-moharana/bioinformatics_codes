
#Date: Qua 22 Jun 2016 10:34:34 BRT 
###kcm.eid[at]gmail.com
#Extract given list of Pfam HMM models from PFAM_A hmm file;




use strict;

if($#ARGV<1){die "####Extract given list of Pfam HMM models from PFAM_A.hmm file###
Usage: perl $0 <PFAM_A.hmm> <file_with_ACC_1_per_line>   >   OUTPUT.hmm \n"}
my $PFAM_A=$ARGV[0] || die "HMM file missing!!!";
my $ACC_list=$ARGV[1] || die "list of ACC missing !!!";

my %acc;
open (F,"$ACC_list") or die"";
while(<F>){chomp; ~s/\s+//g; $acc{$_}=1}
close F;
my $f=0;
my $data;

open (F,"$PFAM_A") or die"";
while(<F>){
	chomp;

	if(/\/\//){
		$data .= "$_\n";
		print $data if $f;
		undef $data; undef $f;
	}
	elsif(/^ACC\s+(\S+)\.\d*/){
		
		$f=(exists $acc{$1}?1:0);
		$data .= "$_\n";
	}
	else{
		$data .= "$_\n";
	}

	

}
close F;

warn "\n\nNo. of models exported:".scalar keys %acc;
