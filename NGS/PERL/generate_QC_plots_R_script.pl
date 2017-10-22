
## This will auto-generate a "R script", that can be used to generate pretty box-plots; 
#	the output need to be executed in R to obtain plot

## input files must be from ´fastx_quality_stats´ program
# Command for fastx_quality_stats: fastx_quality_stats [-h] [-i INFILE] [-o OUTFILE]

=input file format
#File_name	label
ERR032207_1.fastq.txt	A1
ERR032207_2.fastq.txt	A2

=cut

use strict;

if($#ARGV<3)
{
	die "

Insufficient inputs!!!
Script to generate a R-script for pretty fastq quality plot.

Usage: perl $0 -f Input_list -t <pe|se>
	-f <file> : a tab file containing two columns (Fastx_quality_stat_output_file_names and labeles; eg. 
	  #File_name	label		<-- comment line
	  ERR032207_1.fastq.txt	A1
	  ERR032207_2.fastq.txt	A2

	-t <pe|se> : write 'pe' for paired end, 'se' for singel-end fastq files; don't mix;read length must be uniform
* Default output to Screen; use > to save to file;
#KCM.EID[at]GMAIL.com
";
}

my %opts = @ARGV;
my $inp_list = $opts{-f};
my $fq_lib = $opts{-t};
my ($rows,$cols)=(0,0);	## for PE data 2 columns, label as A1,A2
my $read_len=0;


open (I, "$inp_list") or die "$inp_list $!\n";
my @fq_files=<I>;

close I;

my ($temp,$temp2)=split /\t/,$fq_files[-1];
open (F, "$temp") or die "$temp $!\n";
my @a=<F>; $read_len=$#a;
close F;

if($fq_lib eq 'pe'){ ($rows,$cols) =((scalar @fq_files)/2,2)}
elsif($fq_lib eq 'se'){ ($rows,$cols) =(scalar @fq_files,1)}
else {die "Unknown -t option; use se or pe\n";}


print 'par(mfrow =  c('.$rows.','.$cols.'), mar=c(1,0.4,0.2,0.2), oma=c(1.5,2,1,1), mgp = c(1,0.5,0))

cex_axis=1.2
cex_lab=1.5
#dsizeGrWindow(10, 9)
#pdf(file = "QC_plot.pdf", width = 10, height = 9)
';


my $i=0;

foreach(@fq_files)
{
chomp;
next if (/^#/);
my @l=split /\t/,$_;
#xlim=c(0,'.$read_len.')
#Nucleotide positions in read
#Quality score
if($i%2==0){
print 	"#plot for $l[0]
		".'data1 <- read.table("'.$l[0].'", header=1)
		summarydata1<-list(stats=rbind(data1$lW,data1$Q1,data1$med,data1$Q3,data1$rW), n=data1$count, names=data1$column)
		bxp(summarydata1,boxwex=0.5,whisklty="solid", medlty="solid", medlwd=1, pars=list( ylim=c(0,40),xaxt="n", xlab="",ylab="",axes=F) )
		axis(1,pos=0, at = seq(10, '.$read_len.', by = 10), las=1,tck=-0.03,cex.axis=cex_axis, cex.lab=cex_lab)
		axis(1,pos=0, at = seq(5, '.$read_len.', by = 5),labels=character(length(seq(5, '.$read_len.', by = 5))), las=1,tck=-0.02)	
		axis(1,pos=0, at = seq(1, '.$read_len.', by = 1),labels=character(length(seq(1, '.$read_len.', by = 1))), las=1,tck=-0.01)			
		axis(2,pos=0,at = seq(0,40, by = 10),las=1,tck=-0.02,cex.axis=cex_axis, cex.lab=cex_lab)
		axis(2,pos=0,at = seq(0,40, by = 5),las=1,tck=-0.01 , labels=character(length(seq(0,40, by = 5))),cex.axis=cex_axis, cex.lab=cex_lab)
		rect(0,0,'.$read_len.'+1,40,col = NULL, border ="black")
		lines(cbind(data1$mean),pch = 22, col = "red", lwd = 1)
		text(2,3,"'.$l[1].'",cex=1.5)
';
	}
else{
	print "#plot for $l[0]
	".'data <- read.table("'.$l[0].'", header=1)
	summarydata<-list(stats=rbind(data$lW,data$Q1,data$med,data$Q3,data$rW), n=data$count, names=data$column)
	bxp(summarydata,boxwex=0.5,whisklty="solid", medlty="solid", medlwd=1, pars=list(ylim=c(0,40),xaxt="n",xlab="",ylab="",cex.axis=cex_axis, cex.lab=cex_lab,axes=F))
	axis(1,pos=0, at = seq(10, '.$read_len.', by = 10), las=1,tck=-0.03,cex.axis=cex_axis, cex.lab=cex_lab)
	axis(1,pos=0, at = seq(5, '.$read_len.', by = 5),labels=character(length(seq(5, '.$read_len.', by = 5))), las=1,tck=-0.02)	
		axis(1,pos=0, at = seq(1, '.$read_len.', by = 1),labels=character(length(seq(1, '.$read_len.', by = 1))), las=1,tck=-0.01)
	axis(2,pos=0,at = seq(0,40, by = 10),las=1,tck=-0.02 ,cex.axis=cex_axis, cex.lab=cex_lab)
	axis(2,pos=0,at = seq(0,40, by = 5),las=1,tck=-0.01 , labels=character(length(seq(0,40, by = 5))),cex.axis=cex_axis, cex.lab=cex_lab)
	rect(0,0,'.$read_len.'+1,40,col = NULL, border ="black")
	lines(cbind(data$mean),pch = 22, col = "red", lwd = 1)
	text(2,3,"'.$l[1].'",cex=1.5)
	';


	}
$i++;

}
#text(2,3,"Nucleotide positions in read",cex=1.5)

print 'mtext("Nucleotide position in read",1, adj=0.5, line=0.2, outer=TRUE)
mtext("Quality score",2, adj=0.6, line=0.1, outer=TRUE)


';




print STDERR '

DONE; Run the output in R-prompt using source command
';
