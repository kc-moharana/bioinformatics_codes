##contigs.fa at command line

my ($len,$total,$readCount)=(0,0,0);
my @x;
while(<>)
{
	if(/^[\>\@]/){if($len>0){$total+=$len;push@x,$len;};$len=0;$readCount++; }
	else{s/\s//g;$len+=length($_);}
}
print "Total Count: $readCount\n" ;
print "Total Size: ".$total/1000000;
print "mb\nAvg. Length: ".$total/$readCount;
if ($len>0){$total+=$len;push @x,$len;}
@x=sort{$b<=>$a}@x;
my ($count,$half)=(0,0);
for (my $j=0;$j<@x;$j++){$count+=$x[$j];
if(($count>=$total/2)&&($half==0)){print "\nN50: $x[$j]\n";$half=$x[$j]}
elsif($count>=$total*0.9){print "N90: $x[$j]\n";exit;}}


