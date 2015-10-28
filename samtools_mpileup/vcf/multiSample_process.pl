
############################################################################################
#THIS SCRIPT WILL PROCESS MULTI-SAMPLE VCF FILE GENERATED FROM samtools mpileup PROGRAM.
#samtools VERSION 1.0+htslib-1.0
#MUST BE CALLED WITH '-t DP4' PARAMETERS
#	ID=DP4,Description="Number of high-quality ref-forward , ref-reverse, alt-forward and alt-reverse bases">
#
#
#
#
#
#
#
#
#
#
#chr1	14907	.	A	G	5.19219	.	DP=4;SGB=-0.516033;RPB=1;MQB=1;BQB=1;MQ0F=0;ICB=1;HOB=0.5;AC=1;AN=2;DP4=1,0,1,0;MQ=50	GT:PL:DP4	0/1:34,0,34:1,0,1,0	./.:0,0,0:0,0,0,0

=input
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##samtoolsVersion=1.0+htslib-1.0
##samtoolsCommand=samtools mpileup -I -g --min-BQ 25 -t DP4 -f /home/BMIC/RESOURCES/human_genome_sequence_hg19/hg19.fa -o mpileup_Liver_Kidney_uniq_map_hs.bcf Liver_merged_accpeted_RG_corrected.sort.bam Kidney_merged_accepted_hit_RG_corrted.sort.bam
##reference=file:///home/BMIC/RESOURCES/human_genome_sequence_hg19/hg19.fa
##contig=<ID=chr1,length=249250621>
##contig=<ID=chrUn_gl000249,length=38502>
##contig=<ID=chrX,length=155270560>
##contig=<ID=chrY,length=59373566>
##ALT=<ID=X,Description="Represents allele(s) other than observed.">
##INFO=<ID=INDEL,Number=0,Type=Flag,Description="Indicates that the variant is an INDEL.">
##INFO=<ID=IDV,Number=1,Type=Integer,Description="Maximum number of reads supporting an indel">
##INFO=<ID=IMF,Number=1,Type=Float,Description="Maximum fraction of reads supporting an indel">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw read depth">
##INFO=<ID=VDB,Number=1,Type=Float,Description="Variant Distance Bias for filtering splice-site artefacts in RNA-seq data (bigger is better)",Version=3>
##INFO=<ID=RPB,Number=1,Type=Float,Description="Mann-Whitney U test of Read Position Bias (bigger is better)">
##INFO=<ID=MQB,Number=1,Type=Float,Description="Mann-Whitney U test of Mapping Quality Bias (bigger is better)">
##INFO=<ID=BQB,Number=1,Type=Float,Description="Mann-Whitney U test of Base Quality Bias (bigger is better)">
##INFO=<ID=MQSB,Number=1,Type=Float,Description="Mann-Whitney U test of Mapping Quality vs Strand Bias (bigger is better)">
##INFO=<ID=SGB,Number=1,Type=Float,Description="Segregation based metric.">
##INFO=<ID=MQ0F,Number=1,Type=Float,Description="Fraction of MQ0 reads (smaller is better)">
##INFO=<ID=I16,Number=16,Type=Float,Description="Auxiliary tag used for calling, see description of bcf_callret1_t in bam2bcf.h">
##INFO=<ID=QS,Number=R,Type=Float,Description="Auxiliary tag used for calling">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="List of Phred-scaled genotype likelihoods">
##FORMAT=<ID=DP4,Number=4,Type=Integer,Description="Number of high-quality ref-fwd, ref-reverse, alt-fwd and alt-reverse bases">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##INFO=<ID=ICB,Number=1,Type=Float,Description="Inbreeding Coefficient Binomial test (bigger is better)">
##INFO=<ID=HOB,Number=1,Type=Float,Description="Bias in the number of HOMs number (smaller is better)">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes for each ALT allele, in the same order as listed">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##INFO=<ID=DP4,Number=4,Type=Integer,Description="Number of high-quality ref-forward , ref-reverse, alt-forward and alt-reverse bases">
##INFO=<ID=MQ,Number=1,Type=Integer,Description="Average mapping quality">
##bcftools_callVersion=1.0+htslib-1.0
##bcftools_callCommand=call -v -m mpileup_Liver_Kidney_uniq_map_hs.bcf
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	LIVER	KIDNEY
chr1	14907	.	A	G	5.19219	.	DP=4;SGB=-0.516033;RPB=1;MQB=1;BQB=1;MQ0F=0;ICB=1;HOB=0.5;AC=1;AN=2;DP4=1,0,1,0;MQ=50	GT:PL:DP4	0/1:34,0,34:1,0,1,0	./.:0,0,0:0,0,0,0
chr1	14930	.	A	G	14.9926	.	DP=1;SGB=-0.157211;MQ0F=0;AC=2;AN=2;DP4=0,0,1,0;MQ=50	GT:PL:DP4	1/1:40,3,0:0,0,1,0	./.:0,0,0:0,0,0,0
chr1	14937	.	T	C	14.9926	.	DP=1;SGB=-0.157211;MQ0F=0;AC=2;AN=2;DP4=0,0,1,0;MQ=50	GT:PL:DP4	1/1:40,3,0:0,0,1,0	./.:0,0,0:0,0,0,0
chr1	14975	.	C	T	7.47927	.	DP=3;SGB=-0.516033;RPB=1;MQB=1;MQSB=1;BQB=1;MQ0F=0;ICB=0.5;HOB=0.5;AC=2;AN=4;DP4=0,1,1,0;MQ=50	GT:PL:DP4	0/1:0,3,40:0,1,0,0	0/1:40,3,0:0,0,1,0
chr1	15053	.	C	T	6.48096	.	DP=3;SGB=-0.516033;RPB=1;MQB=1;MQSB=1;BQB=1;MQ0F=0;ICB=0.3;HOB=0.125;AC=1;AN=4;DP4=0,2,1,0;MQ=50	GT:PL:DP4	0/1:40,3,0:0,0,1,0	0/0:0,6,73:0,2,0,0
chr1	15094	.	A	G	97	.	DP=9;VDB=0.0949329;SGB=0.908071;MQ0F=0;AC=2;AN=2;DP4=0,0,0,4;MQ=50	GT:PL:DP4	1/1:122,12,0:0,0,0,4	./.:0,0,0:0,0,0,0
chr1	15150	.	C	T	44.8707	.	DP=10;VDB=0.21346;SGB=0.620439;RPB=0.150134;MQB=0.900802;MQSB=0.900802;BQB=0.900802;MQ0F=0;ICB=0.3;HOB=0.125;AC=1;AN=4;DP4=3,2,0,3;MQ=50	GT:PL:DP4	0/0:0,3,40:1,0,0,0	0/1:78,0,119:2,2,0,3
chr1	36067864	.	C	A,T	15.003	.	DP=33;VDB=0.0949329;SGB=-3.1648;RPB=0.131281;MQB=1;MQSB=1;BQB=0.615516;MQ0F=0;ICB=0.5;HOB=0.5;AC=1,1;AN=4;DP4=4,13,2,2;MQ=50	GT:PL:DP4	0/1:44,0,181,68,187,239:1,7,0,2	0/2:40,67,255,0,231,225:3,6,2,0
chr1	201342509	.	A	G,C	18.2189	.	DP=2;VDB=0.06;SGB=0.0985265;MQ0F=0;AC=1,1;AN=2;DP4=0,0,2,0;MQ=50	GT:PL:DP4	1/2:71,37,34,37,0,34:0,0,2,0	./.:0,0,0,0,0,0:0,0,0,0


=cut


=Updates


=cut


=output
V_ID	CHRM	POS-1	POS	REF	ALT	SNP_ID	ALLELE_COUNT	LIVER_A	LIVER_T	LIVER_G	LIVER_C	LIVER_DP	LIVER_TYPE	KIDNEY_A	KIDNEY_T	KIDNEY_G	KIDNEY_C	KIDNEY_DP	KIDNEY_TYPE
chr1:14907	chr1	14906	14907	A	G	.	2	1	0	1	0	2	Polymorphic	0	0	0	0	0	NO_COV
chr1:14930	chr1	14929	14930	A	G	.	2	0	0	1	0	1	ALT	0	0	0	0	0	NO_COV
chr1:14937	chr1	14936	14937	T	C	.	2	0	0	0	1	1	ALT	0	0	0	0	0	NO_COV
chr1:14975	chr1	14974	14975	C	T	.	2	0	0	0	1	1	REF	0	1	0	0	1	ALT
chr1:15053	chr1	15052	15053	C	T	.	2	0	1	0	0	1	ALT	0	0	0	0	2	REF
chr1:15094	chr1	15093	15094	A	G	.	2	0	0	4	0	4	ALT	0	0	0	0	0	NO_COV
chr1:15150	chr1	15149	15150	C	T	.	2	0	0	0	0	1	REF	0	3	0	4	7	Polymorphic
chr1:36067864	chr1	36067863	36067864	C	A,T	.	3	2	0	0	8	10	Polymorphic	0	2	0	9	11	Polymorphic
chr1:201342509	chr1	201342508	201342509	A	G,C	.	3	0	0	0	2	2	ALT	0	0	0	0	0	NO_COV

=cut





use strict;
my $vcf_inp=shift || die "No input vcf found!!!\nUsage: perl $0 <in.vcf> Exit..\n";
open(V, "<$vcf_inp") or die "$! $vcf_inp\n";

my ($vcf_flg,$samtools_v_flg,$gt_format_flg)=(0,0,0);
my $number_of_samples=0; my %gt;
my @sample_names;

while(<V>){
	
	if(/^#/){	
		if(/^##samtoolsVersion=(\S+)/){if ($1 eq '1.0+htslib-1.0'){$samtools_v_flg=1; } else{ die "samtools version is not '1.0+htslib-1.0'\n exit.." }}
		if(/^#CHROM	./g){ $vcf_flg=1; my @l=split /\s+/,$_;  $number_of_samples=scalar@l-9; print STDERR "$number_of_samples samples detected\n";
		print "V_ID\tCHRM\tPOS-1\tPOS\tREF\tALT\tSNP_ID\tALLELE_COUNT";	
		foreach my $i(9..$#l){ print "\t$l[$i]_A\t$l[$i]_T\t$l[$i]_G\t$l[$i]_C\t$l[$i]_ALT_count\t$l[$i]_DP\t$l[$i]_TYPE"; push @sample_names,$l[$i];  } print "\n";
		die "unknown samtools verion!!! You must have deleted the line ##samtoolsVersion in VCF header\nExit...\n" if !$samtools_v_flg }
		next;
	}
	else{
	die "Error in vcf file, No header and no samtools verion line in header .. Exit\n" if (!$vcf_flg or !$samtools_v_flg);
	my @r=split /\s+/,$_;
	
	if ($r[7] =~m/^INDEL;/) {warn "$r[0]:$r[1] INDEL entry found. SKIPPED...\n"; next; }
	if ($r[8] ne 'GT:PL:DP4') { die "no DP4 per sample detected, call mpileup with -t DP4 option (samtools VERSION 1.0+htslib-1.0)\n "}
	my $p_1=$r[1]-1;
	print "$r[0]:$r[1]\t$r[0]\t$p_1\t$r[1]\t$r[3]\t$r[4]\t$r[2]";
	
	#$gt{0}= uc ($r[3]); 
	my @alleles=split /,/,$r[4];
	#foreach (1..scalar @alts){$gt[$_]=$alts[$_-1]};
	unshift @alleles, uc($r[3]);
	print "\t".scalar @alleles;
	
	foreach my $i(9..$#r)
	{
	 my (%atgc_counts,$tot_DP,$typ,$alt_count);
	 
		if ($r[$i]=~m/([\.\d*])\/([\.\d*])\:\S+\:(\d+),(\d+),(\d+),(\d+)/){ 
		if($1 eq $2){
			$atgc_counts{$alleles[$1]}=$3+$4+$5+$6; 
		}
		else{
			$atgc_counts{$alleles[$1]}=$3+$4; 
			$atgc_counts{$alleles[$2]}=$5+$6; 
		}
		
		$tot_DP=$3+$4+$5+$6; $typ=((($3+$4)>=1 && ($5+$6)>=1)?"Polymorphic": (($3+$4)>=1 && ($5+$6)==0?"REF":(($3+$4)==0 && ($5+$6)>=1?"ALT":"NO_COV")));
		$alt_count=$5+$6;
		}
		printf "\t%d\t%d\t%d\t%d\t%d\t%d\t%s",$atgc_counts{A},$atgc_counts{T},$atgc_counts{G},$atgc_counts{C},$alt_count,$tot_DP,$typ;
		
	}
	print "\n";
	}
}
close V;



=further pipeline
E:\labWork_SGPGI\FINAL_RNA_EDITING_ARRANGED\3.Editings_in_human\1.mpileup\SNP_call_31-dec-2014>perl batchUCSC.pl -d hg19 -p "snp138:::" snp_extract_test.txt > snp_extract_test.out.txt

##download gne info in bed fromat from UCSC
CHROM	START	END	GENE_NAME	X	STRAND	
chr1	11873	14409	uc001aaa.3	0	+	11873	11873	0	3	354,109,1189,	0,739,1347,
chr1	11873	14409	uc010nxr.1	0	+	11873	11873	0	3	354,52,1189,	0,772,1347,
chr1	11873	14409	uc010nxq.1	0	+	12189	13639	0	3	354,127,1007,	0,721,1529,
chr1	14361	16765	uc009vis.3	0	-	14361	14361	0	4	468,69,147,159,	0,608,1434,2245,
chr1	16857	17751	uc009vjc.1	0	-	16857	16857	0	2	198,519,	0,375,

##create a database named geneInfo.db
perl -lane "print \"$F[0]|$F[1]|$F[2]|$F[3]|$F[5]\n\"" > geneInfo.inp.tmp
#write sql script to extract required positions and its related gene and strand

sqlite3.exe geneInfo.sqlite < import.sql

Where the content of import.sql is:

CREATE TABLE geneInfo (CHROM VARCHAR(7), START INT(12), END INT(12), GENE_NAME VARCHAR(30), STRAND VARCHAR(1));
.import output.csv geneInfo

E:\labWork_SGPGI\FINAL_RNA_EDITING_ARRANGED\3.Editings_in_human\1.mpileup\SNP_call_31-dec-2014>c:\sqlite3.exe -header -csv geneInfo.db "SELECT * FROM  GTF WHERE CHROM=chr1 and (START>=12344 and END<=12556);"












=cut







