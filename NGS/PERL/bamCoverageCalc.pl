## Script to find Genome wide coverage from multiple BAM files in a directory
## dependency on:  samtools depth
#Date:29/12/12
# last update: 31/dec/12
## sorting flag added; threading added;
=chromosome_size_file example : 
	Chr1	43270923
        Chr2	35937250
        Chr3	36413819
        Chr4	35502694
        Chr5	29958434
        Chr6	31248787
        Chr7	29697621
=cut

use strict;
use threads;


if($#ARGV<3)
{
print "\nUsage:perl bamCoverageCalc.pl [OPTIONS...] -i <Directory_containing_BAMfiles> -s <chromosome_size_file> -sort <1|0> \n
 Directory_containing_BAMfiles should NOT contain /.
Options:
letter following a tag indicates datatype. N:Integer;F=Float;S=String;B=boolean; e.g: -r N, input for -r should be an integer.
    paramaters:-
	
	-i S: Directory Name (full path) containing BAMfiles [required]
	-s S: File name, containing chromosome sizes in clumns;
		 Col.1=Chr1  col.2=Chro.Size.[required]
	-r N: Read Depth of SNP base	[Default: N= 1]
	-sort B: sort the bam file first, 1=YES, 0=NO [Default 0]
	";

exit(1);
}


my($day, $month, $year)=(localtime)[3,4,5];
$year+=1900;
$month+=1;
my $date="$day/$month/$year";

my %argvHash=@ARGV;

			
my $BAMdir=$argvHash{'-i'};			##required
my $refSeqFile=$argvHash{'-s'};		##required
my $readDepth=$argvHash{'-r'}||1; 	# read Depth to consider[default 1]
my $toSort=$argvHash{'-sort'}||0;	##sort falg


$BAMdir=~s/\///g;
die "Directories $BAMdir Not found !!!\n" if !(-d $BAMdir);
chdir ($BAMdir);
my @bamFiles=<*.bam>;		
die"No BAM file found!!!!files must have .bam extension" if !@bamFiles;
chdir ('..');


my @threads;
foreach my $file(@bamFiles)
	{
	my $thr=threads->new(\&process_bam,$file,$BAMdir,$refSeqFile,$toSort  );
	push @threads,$thr;
    	}

$_->join for @threads;



sub process_bam
{
my $file=shift;
my $BAMdir=shift;
my $size_FH=shift;
my $sortFlg=shift;
my %size;				## all chrom size
my %printData;
open(size_FH,"< $size_FH") or die"cannt  read $size_FH";
while(<size_FH>){chomp; my($chrName,$chrLen)=split/\s+/,$_; next if !$chrName; $size{$chrName}=$chrLen; $printData{$chrName}=[0,0,1000,0,0,0]; }

close size_FH;
$file="$BAMdir/$file";
	if($sortFlg){
		my $cmd1="samtools sort $file $file\_sort"; print "$cmd1\n";
		system"$cmd1";
		my $file="$file\_sort.bam";
	}
	my $cmd2="samtools depth $file > $file\.depth"; print "$cmd2\n";
	system"$cmd2";
	
	warn "Read Depth: $readDepth\n";
	open(DEP,"< $file\.depth") or die "cannt open $file\.depth";
	open(COV,"> $file\.coverage") or die "cannt write $file\.coverage";
	print COV "$date\nRead Depth considered : $readDepth\n$cmd2\n\#ChrmoNo\tDepthCoverage(%)\tmin.Depth-Max.Depth\tCovered by ReadCoverage(%)\n";
	while (<DEP>)
		{
	    	my @a = split/\s+/;
	    	print "V@a$size{$a[0]}\n" if !$size{$a[0]};
	        #$readCover++ if $a[2]>=$readDepth;
		$printData{$a[0]}->[0]+= $a[2];	##num
		$printData{$a[0]}->[1]++;	##len
		$printData{$a[0]}->[2]=$a[2] if $printData{$a[0]}->[2] > $a[2];	##min
		$printData{$a[0]}->[3]=$a[2] if $printData{$a[0]}->[3] < $a[2];	##max
		
		$printData{$a[0]}->[4]++ if  $a[2]>=$readDepth;			#=($readCover/$size{$chrNo})*100;
		$printData{$a[0]}->[5]=($printData{$a[0]}->[4]/$size{$a[0]})*100;	#
		
		}
	close DEP;
	my ($genoCov,$eachCov,$totalSize)=(0,0,0);
	foreach my $key(sort keys %printData)
		{
		next if !$printData{$key}->[0];
		printf COV "%s\t%.1f\t%d-%d\t%\d\t%.2f\n", $key, ($printData{$key}->[0]/$printData{$key}->[1]), $printData{$key}->[2], $printData{$key}->[3],$printData{$key}->[4], $printData{$key}->[5];
		#print "\n$key\t$printData{$key}->[0]\t$printData{$key}->[1]\t$printData{$key}->[2]\t$printData{$key}->[4]!\n";
		$eachCov+=$printData{$key}->[4];
		$totalSize+=$size{$key};
		#$printData{$key}=[0,0,1000,0,0,0];
		
		}
	$genoCov=($eachCov/$totalSize)*100;
	printf COV "Genome Wide Coverage: %.2f\n",$genoCov;
  	close COV;
	print "$file processed\n\n";
}




