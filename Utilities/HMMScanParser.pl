# Script to parse HMMSCAN output and filter on user defined e-val thresholds, list of genes, list of domains.


#KCM.EID[at]GMAIL.COM

##Features
#1. Filter using e-value (full-seq, i_E-valDomain,c_E-value)
#2. Filter using bit-score (domain bit score only )
#3. Filter gene list
#4. Filter Query list
#Yet to add
#5. Create graphical representation of each gene

use strict;



=input-1
#                                                                            --- full sequence --- -------------- this domain -------------   hmm coord   ali coord   env coord
# target name        accession   tlen query name           accession   qlen   E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias  from    to  from    to  from    to  acc description of target
#------------------- ---------- ----- -------------------- ---------- ----- --------- ------ ----- --- --- --------- --------- ------ ----- ----- ----- ----- ----- ----- ----- ---- ---------------------
Vang0032ss00130.1    -            219 1-cysPrx_C           PF10417.6     40   2.8e-11   43.6   0.1   1   1   1.8e-15   4.9e-11   42.8   0.1     1    40   161   200   161   200 0.93 -
Vang0045ss01500.2    -            301 12TM_1               PF09847.6    449      0.15   11.7  18.2   1   1     7e-06      0.19   11.4  18.2    12   134   130   250   119   279 0.83 -
Vang0063ss00300.1    -            259 14-3-3               PF00244.17   222  4.5e-106  353.9   7.3   1   1  2.7e-109  5.1e-106  353.7   7.3     2   222    14   235    13   235 0.99 -
Vang0017s00590.1     -            260 14-3-3               PF00244.17   222  1.5e-105  352.2   1.9   1   1    9e-109  1.7e-105  352.0   1.9     1   222    14   238    14   238 0.99 -
Vang01g10280.1       -            262 14-3-3               PF00244.17   222  7.9e-105  349.8   0.9   1   1  4.7e-108    9e-105  349.6   0.9     1   222    14   238    14   238 0.99 -
Vang04g00960.1       -            257 14-3-3               PF00244.17   222  5.2e-104  347.1   2.6   1   1  3.1e-107    6e-104  346.9   2.6     1   222    12   236    12   236 0.99 -
Vang0184s00490.1     -            280 14-3-3               PF00244.17   222   1.1e-99  333.1  11.8   1   1    3e-100   5.7e-97  324.1  11.8     2   222    14   253    13   253 0.93 -

=cut


=input-2

=cut






=data_hast DS
%data{
	-gene=>{
		-acc={
			-#_of=>[],		
			-#_of=>[],
		}	
	}

}
=cut





#############################################################################






if($#ARGV<1){die "Perl $0 [OPTIONS] -i Hmmscan_domOutput\n
#Script filters a hmmscan output (strictly domtable output) on several parameters
 -l <file> : list_of_dom
 -g <file> : list_of_genes
 -e <real_number> : e_value [default: 10]
 -b <real_number> : bit_score [default: 0]

"; }
my %opts=@ARGV;

my $e_val = $opts{-e} || 10;
my $bit_score = $opts{-b} || 0;
my $out_file = $opts{-o} || 'STDOUT';

if($out_file ne 'STDOUT'){
	open (OUT, ">$out_file" ) or die "OUTPUT file (-o)  $out_file $!\n";
	select OUT;
}


my(%data, @gene_list,@acc_list, @excluded_domain);
open(F, $opts{-i}) or die "$! $opts{-i}\n";

if($opts{-l}){
	open(L, $opts{-l}) or die "$! $opts{-l}\n";
	@acc_list=<L>;
	close L;
}

if($opts{-g}){
	open(G, $opts{-g}) or die "$! $opts{-g}\n";
	@gene_list=<G>;
	close G;
}

my(%acc_list,%gene_list);
while(<F>)
{
	chomp; 
	if (/^#/){next; }
	my @l = split /\s+/,$_;
	
	

	if (($l[6] < $e_val && ($l[11] < $e_val && $l[12] < $e_val)) && $l[13] > $bit_score ){
		$acc_list{$l[3]}=1 if !$opts{-l};
		$gene_list{$l[0]}=1 if !$opts{-g};		
		#$data{$l[0]}->{$l[3]}->{"$l[9]_$l[10]"}=$_;
		$data{$l[0]}->{$l[3]}->{$l[9]}=\@l;
	}
	else {
		push @excluded_domain,$_;
	}	
}
close F; 
@gene_list =keys %gene_list if ! $opts{-g};
@acc_list= keys %acc_list if ! $opts{-l};



print STDERR "Summary: \nTotal genes: ".($#gene_list+1)."\n (-)Excluded domains(e-val > $e_val and Bit-score < $bit_score): ".($#excluded_domain+1)."\nTotal Pfam domains: ".($#acc_list+1)."\n";

my($g,$a,$c);
foreach $g (@gene_list)
{
	$g =~s/\n//g;
	foreach $a(@acc_list)
	{
		$a =~s/\n//g;
		while (my ($c, $val) = each %{$data{$g}->{$a}})
		{
			$"="\t";
			
			print "@{$data{$g}->{$a}->{$c}}\n";

		}
	}
delete $data{$g};
}

close OUT if($out_file ne 'STDOUT');













