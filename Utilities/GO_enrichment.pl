##This script will help to perform Gene Ontology Enrichment analysis;
##Provides  Pearson's Chi-squared test with Yates' continuity correction and p-value
##


=theroy
## Gene Ontology Theory
A gene or gene product can be  included in three independent broad-categories  
> biological processes (BP) : Biological objective of the gene
> molecular process (MP) : biochemical activity of the gene
> cellular component (CC) : the plcae in the cell where the gene product is active

Each gene can have more then one GO term. Each GO term is given a unique 7 digit integer identifier (e.g GO:0002345).  There are databases where we can get GO id for each gene of a model species (e.g E.coli for bacteria). Visit : http://geneontology.org/page/download-ontology/ . I will prefere the MuSQL dump files. Text formats are also available.

### How a gene is assigned with GO id
	For each species there are several consortiums for data processing and assigning GO category. [more information to be added]
	They gather evidences from several sources to assign functional term to a gene. There are eleven types of evidences for asigning GO term to a gene.
	[i] IMP: Inferred from mutant phenotype
	[ii] IGI: Inferred from geneic interaction
	[iii] IPI: Inferred from physical interaction
	[iv] ISS: Inferred from sequence similarity
	[v] IDA: Inferred from direct assay
	[vi] IEP: Inferred from expression pattern
	[vii] IEA: Inferred from electronic annotation (lowest qality evidence)
	[viii] TAS: Tracabe author statement
	[ix] NAS: Non-tracable statements
	[x] ND: No biological data avaibale
	[xi] IC: Inferred by curator
	
###GO data representation (computer science part)
Each GO term is a node of Directed Acyclic Graph (DAG). DAG is similar to a regular graph (graph theory;node--edge-node), but here each node (or gene) is connected to multiple parents and the edges are directed. The directionality reprepresents the fact that a child-node (or gene) is either 'part of' or 'instance of' their parents. 
	    [g1][g2]
	     /\  /\
	 [gp][g3][g5] 

### GO Enrichment 	statistics
Enrichment analysis requires two datasets.
[i] A background set of genes  (usually all genes of an organism)
[ii] A test set of genes (genes of interest e.g: differentially expressed genes// hub gene list//genes specific mutation// any sub set of background genes)

Lets assume that we have N number of background genes whose GO-ids are known and here are C number of unique GO-ids. For simplisity I will consider only one GO-category named F. Now we can devide N genes in to two categories: i) Genes in F ; ii) genes NOT in F.

On mapping test-set genes list with background-set gene list; we can estimate how many genes from test-gene are in F category and how many are NOT in F-category.

This resembles with the problem similar to a urnfilled with N-number of balls with two colors namely red(F-category) and green (NOT-F-category). 
Say the number of red balls is M and hence number of green balls is N - M [N minus M].

Let's asume that we picked K number of balls randomly (random sampling: non-replacement type) and K may caontain 'x' number of red balls. 
So we want to check a null hypothesis that the number of red balls 'x' out of K sambling is by chance.  

So we will use Chi-square test for this null hypothesis test and p-value calculation. 
First of all form a 2X2 contigency matrix
				Background_gene		Test_gene	| Total
	In_category_F		a[n11]			b [n12]		| N1.
	NOT_in_category_F	c [n21]			d [n22]		| N2.
					-------------------------------------------
				Total	a+c [N.1]	b+d [N.2]	| N.. = (a+b+c+d)	
				
	Chi-square val = (N..*(|n12n22 - n12n21| - N../2)**2) / (N1.N2.N.1N.2)
	(as sample size is small we have to substract N../2 ; knownas Yates-correction; Fisher exact test can also be used instead)
	
	This Chi-square value has be to be compared with critical values obtained from Chi-squaredistributions with degree of freedom [(2-1)*(2-1) =1]; as there are 2 categories;
	
	P-value: calculated using chisqrprob() from Statistics::Distributions; Just like looking in the chi-square table
	


Chi-square PERL codes implemented from :
1. https://edwards.sdsu.edu/research/calculating-chi-squared-wiht-perl/
2. Statistics::Distributions
3. https://zhouhufeng.wordpress.com/2013/07/31/chi-square-fishers-exact-and-mcnemars-test-using-r/

=cut



=input file formats
##background_gene_file (two columns separated by a tab)
#Gene	GO:id
Gene 1	category_2
Gene 2	category_2
Gene 3	category_6
Gene 4	category_2
Gene 5	category_4
Gene 6	category_5
Gene 7	category_5
Gene 8	category_1
Gene 9	category_4
....	.....

##test_gene_file (list of genes; gene names must match with background gene name)
#
Gene 22
Gene 90
Gene 331
Gene 154
Gene 119
Gene 311
Gene 623
Gene 993
Gene 657
Gene 778
Gene 821
Gene 334
Gene 255
Gene 398
Gene 432
Gene 946
...
...


=cut

use strict;
use GD::Graph::hbars;
use constant PI => 3.1415926536;
use constant SIGNIFICANT => 5; 		# number of significant digits to be returned


if($#ARGV<1){ die "#Script for Gene Ontology Enrichment analysis\nUsage: perl $0 <background_gene_with_GO_id.tsv> <test_gene_list>\n"}
my $background_gene_file=$ARGV[0];
my $test_gene_file=$ARGV[1];

my (%all_gene_go, @test_genes, %categories,%test);
my ($a,$b,$c,$d,$df)=(0,0,0,0,1);

open (B,"<$background_gene_file") or die"$! $background_gene_file\n";
open (T,"<$test_gene_file") or die"$! $test_gene_file\n";
while(<B>)
{
 next if(/^#/);
 chomp;
 my @l=split /\t/,$_;
 $all_gene_go{$l[0]}=$l[1];
 if(!exists $categories{$l[1]} ){$categories{$l[1]}={_a=>0,_b=>0,_c=>0,_d=>0};}
 $categories{$l[1]}->{_a}++;
}

while(<T>)
{
	 next if(/^#/);
	 chomp;
	 push @test_genes,$_; 
	 $categories{$all_gene_go{$_}}->{_b}++; 
	 $test{$all_gene_go{$_}}=[] if !$test{$all_gene_go{$_}};
	 push @{$test{$all_gene_go{$_}}},$_;
}

foreach (keys %categories)
{
	$categories{$_}->{_c} = scalar (keys %all_gene_go) - $categories{$_}->{_a};
	$categories{$_}->{_d} =  scalar @test_genes - $categories{$_}->{_b} ;
}

#print "Contigency table:\n";
#print "GO_cat\tBCK_Genes\tGenes\n";

open (O, ">$test_gene_file\_enrich.txt") or die "$! $test_gene_file\_enrich.txt";
print O "#Category\tGenes\ttest_genes_in_category\tBackground_gene_in_category\tP-value\n";
my (@data,@categories, @freq, @p_vals);			## For GD plot
foreach (sort keys %categories)
{
	next if $categories{$_}->{_b} <1;			## if <1 dont print;

	#print "$_\t$categories{$_}->{_a}\t$categories{$_}->{_b}\nNOT_$_\t$categories{$_}->{_c}\t$categories{$_}->{_d}\n";
	
	####Regular Chi-square calculation
	#my $chi_sq = chi_squared($categories{$_}->{_a},$categories{$_}->{_b},$categories{$_}->{_c},$categories{$_}->{_d});
	#my $chi_p_val = chisqrprob($df,$chi_sq);
	#printf "Chi-square value:%4.3f\nP-value: %4.3f\n",$chi_sq,$chi_p_val;
	
	#Chi-square calculation with Yates-correction
	my $chi_sq_y = chi_squared_y($categories{$_}->{_a},$categories{$_}->{_b},$categories{$_}->{_c},$categories{$_}->{_d});
	my $chi_p_val_y = chisqrprob($df,$chi_sq_y);
	#printf "Chi-square value (Y):%4.3f\nP-value: %4.3f\n\n",$chi_sq_y,$chi_p_val_y,$p_val;
	
	push @categories,$_;
	push @freq,$categories{$_}->{_b};
	push @p_vals,sprintf("%4.3f",$chi_p_val_y) ;
	print O "$_\t".join(";",@{$test{$_}});
	print O "\t$categories{$_}->{_b}\t$categories{$_}->{_a}\t$chi_p_val_y\n";
		
}
close B;
close T;
close O;

@data =(\@categories,\@freq);

my $graph_chart = GD::Graph::hbars->new(500,600);
my $data_p = GD::Graph::Data->new([\@categories,\@p_vals] );

$graph_chart-> set (
	x_label=> 'Categories',
	y_label=> "Frequency \n(number at the top of bar is p-value)",
	title=> 'GO Enrichment',
	bar_spacing=>1,
	show_values => $data_p,
	correct_width=>1,
	bgclr=>"white",
	accentclr=>"white",
	labelclr=>"black",
	axislabelclr=>"black",
	borderclrs=>"black",
	fgclr =>"black",
	boxclr=>"white",
	transparent=>0,	
	
);

my $gd =$graph_chart-> plot(\@data) or die $graph_chart-> error;
my $export="$test_gene_file\_enrich.gif";
open(IMG, ">$export") or die $!;
binmode IMG;
print IMG $gd->gif;
close IMG;




##################################################################################
#FUNCTIONS

##Args: a,b,c,d
##Return: Chi-square value

sub chi_squared {
     my ($a,$b,$c,$d) = @_;
     return 0 if($b+$d == 0);
     my $n= $a + $b + $c + $d;
     return (($n*($a*$d - $b*$c)**2) / (($a + $b)*($c + $d)*($a + $c)*($b + $d)));
}

##Args: a,b,c,d
##Return: Chi-square value (Yates-correction applied)
sub chi_squared_y {
     my ($a,$b,$c,$d) = @_;
     return 0 if($b+$d == 0);
     my $n= $a + $b + $c + $d;
	return (($n*(abs($a*$d - $b*$c)-($n/2))**2) / (($a + $b)*($c + $d)*($a + $c)*($b + $d)));
}


sub factorial {
   my $n = shift;
   my $prod = 1;
   $prod *= $n-- while $n > 0;
   return $prod;
}


##Args: degree of freedom, Chi-square value
##Return: P-value as float
sub chisqrprob {
	my ($n,$x) = @_;
	my $p;

	if ($x <= 0) {
		$p = 1;
	} elsif ($n > 100) {
		$p = _subuprob((($x / $n) ** (1/3)
				- (1 - 2/9/$n)) / sqrt(2/9/$n));
	} elsif ($x > 400) {
		$p = 0;
	} else {   
		my ($a, $i, $i1);
		if (($n % 2) != 0) {
			$p = 2 * _subuprob(sqrt($x));
			$a = sqrt(2/PI) * exp(-$x/2) / sqrt($x);
			$i1 = 1;
		} else {
			$p = $a = exp(-$x/2);
			$i1 = 2;
		}

		for ($i = $i1; $i <= ($n-2); $i += 2) {
			$a *= $x / $i;
			$p += $a;
		}
	}
	return $p;
}
##MEta function
sub _subuprob {
	my ($x) = @_;
	my $p = 0; # if ($absx > 100)
	my $absx = abs($x);

	if ($absx < 1.9) {
		$p = (1 +
			$absx * (.049867347
			  + $absx * (.0211410061
			  	+ $absx * (.0032776263
				  + $absx * (.0000380036
					+ $absx * (.0000488906
					  + $absx * .000005383)))))) ** -16/2;
	} elsif ($absx <= 100) {
		for (my $i = 18; $i >= 1; $i--) {
			$p = $i / ($absx + $p);
		}
		$p = exp(-.5 * $absx * $absx) 
			/ sqrt(2 * PI) / ($absx + $p);
	}

	$p = 1 - $p if ($x<0);
	return $p;
}


