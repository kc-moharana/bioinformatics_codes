

##This script willhelp to compare and contrast 2 Tabular data.
## last update: Qua 22 Jun 2016 10:35:13 BRT 
##Kcm.eid[at]gmail.com



use strict;
use File::Basename;
my($f1,$c1,$f2,$c2,$sim)=@ARGV;
 
if ($#ARGV<3){
print STDERR "
#####Compare two tabular files w.r.t specific coulmns (expects unique entries) #######
Usage: perl $0 File1 <col_no> File2 <col_no> [S]

	col_no	: The first column is 1 and so on 
	[S] 	: can be 1. [default:0],if 1 prints statistics, no output generated

Outputs 3 files:
	Intersection	: contains matched lines from both files
	Unique-A	: Unique lines in File1
	Unique-B	: Unique lines in File2

Useful for Venn-diagrams.
";

exit(1);
}
## Convert to 0-indexed 
--$c1; --$c2;
open F1, "$f1" or die "$! $f1";
open F2, "$f2" or die "$! $f2";
my $f1_b =basename($f1);
my $f2_b =basename($f2);

## Store common, unique entries 
my (%common, %file1_uniq, %file2_uniq);

while(<F1>){
	next if /^#/;
	chomp;
	my  @r=split /\t/,$_;
	$file1_uniq{$r[$c1]}=\@r;	
}

while(<F2>){
	next if /^#|^\s*$/;
	chomp;
	my  @r=split /\t/,$_;
	if (exists $file1_uniq{$r[$c2]}){$common{$r[$c2]}=[]; push ( @{$common{$r[$c2]}},@{$file1_uniq{$r[$c2]}},@r);  }
	delete $file1_uniq{$r[$c2]} if (exists $file1_uniq{$r[$c2]});
	$file2_uniq{$r[$c2]}=\@r if !(exists $common{$r[$c2]});
}
close F1; close F2;

## Default print comparative file
unless($sim){
	open (C,">$f1_b\_$f2_b\_COMMON.txt") or die "$f1_b\_$f2_b\_COMMON.txt $!";
	open (U1,">$f1_b\_UNIQ.txt") or die "$f1_b\_UNIQ.txt $!";
	open (U2,">$f2_b\_UNIQ.txt") or die "$f2_b\_UNIQ.txt $!";
	$"="\t";
	foreach my $c(sort keys %common)
	{
		print C "@{$common{$c}}\n";
	}
	close C;
	foreach my $c(sort keys %file1_uniq)
	{
		print U1 "@{$file1_uniq{$c}}\n";
	}
	close U1;

	foreach my $c(sort keys %file2_uniq)
	{
		print U2 "@{$file2_uniq{$c}}\n";
	}
	close U2;
}

printf STDERR "Uniq value in $f1_b: %d\nUniq value in $f2_b: %d\nCommon (file contails entries 2nd file): %d\n", scalar(keys %file1_uniq),scalar(keys %file2_uniq),scalar(keys %common);





