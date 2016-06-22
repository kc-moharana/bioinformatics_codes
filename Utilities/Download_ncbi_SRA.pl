


#Pipeline: To download and process_fastq_files
## Provide a list of files; (tab separated if downlod separate conditions)
#last update: Qua 22 Jun 2016 10:35:23 BRT 

use strict;


my $list=$ARGV[0] || die "********Pipeline: To download and process_fastq_files***********
Usage: perl $0 <FILE>
   FILE : tab delimted file with SRA acc in the frist col.
          Optionally may contain # comment line, sample descriptions 
\n ERROR\n";
my $OUT_SEP='N';
my %sra; my %cat;


#####Binary file Settings#############
my $fastq_dump_bin ='~/SOFTWARES/sratoolkit.2.5.7-ubuntu64/bin/fastq-dump';
my $fastQC_bin ='~/SOFTWARES/FastQC/fastqc';

unless(-e $fastq_dump_bin or -e $fastQC_bin)
{
	print STDERR "\nERROR:Binaries for fastq-dump or fastqc Not found!! \n";
	exit(1);
}
######################################

open(S, $list) or die "$! $list\n";
while(<S>){ next if /^#/ or !(/^(SRR|ERR|DRR)\d+/);chomp; my @r=split /\t/,$_;my ($id,$c)=($r[0],$r[1]); $c=($c?$c:'.');$sra{$id}=$c; $cat{$c}++ if $c ne '.';}
close S;
print STDERR "\n\nInput details\n--\n";
printf STDERR "\tNumber of valid IDs found: %d\n", scalar keys %sra;
printf STDERR "\tNumber of categories found: %d\n", scalar keys %cat;


a:
if((scalar keys %cat)>0){
	printf STDERR "Do you want to put FASTQ files into separate foldes?(Y/N) ";
	chomp($OUT_SEP=uc <STDIN>);
	if($OUT_SEP ne 'Y' and $OUT_SEP ne 'N'){ warn "$OUT_SEP is an invalid choice\n"; goto  a;}	
}


print STDERR "Downloading SRA files... \n";
foreach (sort keys %sra)
{	
	my $dir=$_; #SRR1174205
	$dir =~s/(^(SRR|ERR|DRR)\d{3})\d+/$1/g;	## 
	my $url="";
	if(/^SRR/){$url ="ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/$dir/$_/$_\.sra";}
	elsif(/^ERR/){$url ="ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/ERR/$dir/$_/$_\.sra";}
	elsif(/^DRR/){$url ="ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/ERR/$dir/$_/$_\.sra";}
	else {print STDERR "\tUnknown File format $_ (use SRR* or ERR* or DRR* formats)\tSkip...\n"; next; }
	print STDERR "\t$_\t";
	system ("wget -c $url");
	if(-e "$_\.sra") {print STDERR "\t..done..\n"; sleep(2);}
	else{print STDERR "\t..ERROR..\n";}
}

print STDERR "\nExtracting Reads from SRA files... 

Warning!!!!
fastq-dump creates huge garbage in \$HOME. Check .ncbi folder in your home regularly. 
\n";

foreach (sort keys %sra)
{	
	print STDERR "\t$_\.sra\t";
	system ("$fastq_dump_bin --split-3 $_");
	if(-e "$_\.fastq" || (-e "$_\_1.fastq" and -e "$_\_2.fastq") ) {print STDERR "\t..done..\n"; sleep(5);}
	else{print STDERR "\t..ERROR..\n";}
	sleep(5);
}

print STDERR "\nCompressing FastaQ files \n";

my @fq=<*.fastq>;
foreach (sort @fq)
{	
	print STDERR "\t$_\t";
	system ("gzip -c $_ > $_\.gz ");
	if(-e "$_\.fastq.gz" || (-e "$_\_1.fastq.gz" and -e "$_\_2.fastq.gz") ) {print STDERR "\t..done..\n";}
	else{print STDERR "\t..ERROR..\n";}
	sleep(2);
}

print STDERR "\nFastQC report... \n";
foreach (sort @fq)
{
	print STDERR "\t$_\.gz\t";
	system ("perl $fastQC_bin $_\.gz");
	~s/\..fastq\.gz//g;
	if(-e "$_\_fastqc.html" || (-e "$_\_1_fastqc.html" and -e "$_\_2_fastqc.html") ) {print STDERR "\t..done..\n";}
	else{print STDERR "\t..ERROR..\n";}
}

if($OUT_SEP eq 'Y'){
print STDERR "\nOrganising files... \n";
	foreach (sort keys %sra)
	{
		my $dir =$sra{$_};
		print STDERR "\tmove $_\.gz to $dir\t";

		mkdir $dir unless (-d "$dir");
		system ("mv $_\.fastq.gz $_\_fastqc.* $dir");
		if(-e "$dir/$_\.fastq.gz" || (-e "$dir/$_\_1.fastq.gz" and -e "$dir/$_\_2.fastq.gz") ) {print STDERR "\t..done..\n";}
		else{print STDERR "\t..ERROR..\n";}
	}

}
print STDERR "\n___FINISH___\nPlease delete .sra and .fastq files to free disk space\n\n";


























