##
# This script will use a list of SRX or eqivalent experiment IDs
# - output : tab delimited sample details 

##

use strict;
use LWP::Simple;
use XML::Hash;
use Data::Dumper;

# https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=biosample&id=9684
#to get id of SRX
#https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term=SRX1853183[accn]
## use id to fetch 
# https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=sra&id=18611
#
#https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id=SRX017401

my $xml_converter = XML::Hash->new();



my $db = 'sra';
my @srx_list = ();					#('SRX017401','SRX017406','SRX1853183');
chomp(@srx_list=<STDIN>);
#my @srx_list=('SRX1483308');

my @id_list=();
print "SRX_ID\tSRR_id\tLIBRARY_STRATEGY\tLIBRARY_SOURCE\tLIBRARY_SELECTION\tSAMPLE_DESCRIPTOR\tLIBRARY_DESCRIPTOR\tSAMPLE_DESCRIPTION\tSAMPLE_TITLE\tSAMPLE_ATTRIBUTES\tLIBRARY_LAYOUT\tPLATFORM\tSTUDY_TITLE\tSTUDY_ABSTRACT\tPUBMED\n";

foreach my $i(0..$#srx_list){
	print STDERR "fetching $srx_list[$i] ...\t";
	my $output = get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id='.$srx_list[$i]);
#	warn "$output\n";
	my $l='';
	if($output){
	 eval {
    		$l= $xml_converter->fromXMLStringtoHash($output); 
	 };
	 if (my $e = $@) {
    	 	print STDERR "Something went wrong while processing : $srx_list[$i] ... \n$e\n";
    	 	next;
	 }
	 
	 #print $l->{'EXPERIMENT_PACKAGE_SET'}->{'IDENTIFIERS'}->{'PRIMARY_ID'};
	 #print keys %{$l->{'EXPERIMENT_PACKAGE_SET'}->{'EXPERIMENT_PACKAGE'}};
#	 print join("\t",keys %{$l->{'EXPERIMENT_PACKAGE_SET'}});
	 #my $srx_id =  $l->{'EXPERIMENT_PACKAGE_SET'}->{'EXPERIMENT_PACKAGE'}->{'EXPERIMENT'}->{'IDENTIFIERS'}->{'PRIMARY_ID'}->{'text'};
	 #my $srr_id =  $l->{'EXPERIMENT_PACKAGE_SET'}->{'EXPERIMENT_PACKAGE'}->{'EXPERIMENT'}->{'STUDY_REF'}->{'IDENTIFIERS'}->{'PRIMARY_ID'}->{'text'}
	 my $SAMPLE_DESCRIPTOR= '';
	 $SAMPLE_DESCRIPTOR= $l->{'EXPERIMENT_PACKAGE_SET'}->{'EXPERIMENT_PACKAGE'}->{'EXPERIMENT'}->{'DESIGN'}->{'SAMPLE_DESCRIPTOR'}->{'accession'};
	 $SAMPLE_DESCRIPTOR=~s/[^[:ascii:]]+//g;
	 
	 
	 my $LIBRARY_DESCRIPTOR =  $l->{'EXPERIMENT_PACKAGE_SET'}->{'EXPERIMENT_PACKAGE'}->{'EXPERIMENT'}->{'DESIGN'}->{'LIBRARY_DESCRIPTOR'}->{'LIBRARY_NAME'}->{'text'};
	 $LIBRARY_DESCRIPTOR =~s/[^[:ascii:]]+//g;
	  
	 
	 my $LIBRARY_STRATEGY =  $l->{'EXPERIMENT_PACKAGE_SET'}->{'EXPERIMENT_PACKAGE'}->{'EXPERIMENT'}->{'DESIGN'}->{'LIBRARY_DESCRIPTOR'}->{'LIBRARY_STRATEGY'}->{'text'};
	 $LIBRARY_STRATEGY =~s/[^[:ascii:]]+//g; #RNAseq
	 my $LIBRARY_SOURCE =  $l->{'EXPERIMENT_PACKAGE_SET'}->{'EXPERIMENT_PACKAGE'}->{'EXPERIMENT'}->{'DESIGN'}->{'LIBRARY_DESCRIPTOR'}->{'LIBRARY_SOURCE'}->{'text'};
	 $LIBRARY_SOURCE =~s/[^[:ascii:]]+//g; #TRANSCRIPTOMIC
	 my $LIBRARY_SELECTION =  $l->{'EXPERIMENT_PACKAGE_SET'}->{'EXPERIMENT_PACKAGE'}->{'EXPERIMENT'}->{'DESIGN'}->{'LIBRARY_DESCRIPTOR'}->{'LIBRARY_SELECTION'}->{'text'};
	 $LIBRARY_SELECTION =~s/[^[:ascii:]]+//g; #cDNA
	 
	 
	 my $LIBRARY_LAYOUT = join("", keys %{$l->{'EXPERIMENT_PACKAGE_SET'}->{'EXPERIMENT_PACKAGE'}->{'EXPERIMENT'}->{'DESIGN'}->{'LIBRARY_DESCRIPTOR'}->{'LIBRARY_LAYOUT'}});
	 #my $READ_TYPE = $l->{'EXPERIMENT_PACKAGE_SET'}->{'EXPERIMENT_PACKAGE'}->{'EXPERIMENT'}->{'DESIGN'}->{'SPOT_DESCRIPTOR'}->{'SPOT_DECODE_SPEC'}->{'READ_SPEC'}->{'READ_TYPE'}->{'text'};
	 $LIBRARY_LAYOUT=~s/[^[:ascii:]]+//g;
	 
	 
	 my $PLATFORM = '';
	 foreach my $k(keys %{$l->{'EXPERIMENT_PACKAGE_SET'}->{'EXPERIMENT_PACKAGE'}->{'EXPERIMENT'}->{'PLATFORM'}}){
	  $PLATFORM= $k.'_'.$l->{'EXPERIMENT_PACKAGE_SET'}->{'EXPERIMENT_PACKAGE'}->{'EXPERIMENT'}->{'PLATFORM'}->{$k}->{'INSTRUMENT_MODEL'}->{'text'};
	 }
	 $PLATFORM=~s/[^[:ascii:]]+//g;
	 my $STUDY_TITLE =  $l->{'EXPERIMENT_PACKAGE_SET'}->{'EXPERIMENT_PACKAGE'}->{'STUDY'}->{'DESCRIPTOR'}->{'STUDY_TITLE'}->{'text'};
	 $STUDY_TITLE =~s/[^[:ascii:]]+//g;
	 my $STUDY_ABSTRACT = $l->{'EXPERIMENT_PACKAGE_SET'}->{'EXPERIMENT_PACKAGE'}->{'STUDY'}->{'DESCRIPTOR'}->{'STUDY_ABSTRACT'}->{'text'};
	 my $PUBMED = "";
	if(ref $l->{'EXPERIMENT_PACKAGE_SET'}->{'EXPERIMENT_PACKAGE'}->{'STUDY'}->{'STUDY_LINKS'}->{'STUDY_LINK'} eq 'ARRAY'){
	 foreach my$k(@{$l->{'EXPERIMENT_PACKAGE_SET'}->{'EXPERIMENT_PACKAGE'}->{'STUDY'}->{'STUDY_LINKS'}->{'STUDY_LINK'}}){
	  $PUBMED.=','.$k->{'XREF_LINK'}->{'ID'}->{'text'} if $k->{'XREF_LINK'}->{'DB'}->{'text'} eq 'pubmed';
	  }	  
	 }
	else{
	 $PUBMED= $l->{'EXPERIMENT_PACKAGE_SET'}->{'EXPERIMENT_PACKAGE'}->{'STUDY'}->{'STUDY_LINKS'}->{'STUDY_LINK'}->{'XREF_LINK'}->{'ID'}->{'text'} if $l->{'EXPERIMENT_PACKAGE_SET'}->{'EXPERIMENT_PACKAGE'}->{'STUDY'}->{'STUDY_LINKS'}->{'STUDY_LINK'}->{'XREF_LINK'}->{'DB'}->{'text'} eq 'pubmed';
	 }
	 $PUBMED=~s/[^[:ascii:]]+//g;
	 my $SAMPLE_TITLE = $l->{'EXPERIMENT_PACKAGE_SET'}->{'EXPERIMENT_PACKAGE'}->{'SAMPLE'}->{'TITLE'}->{'text'};
	 $SAMPLE_TITLE=~s/[^[:ascii:]]+//g;
 	 my $SAMPLE_DESCRIPTION = $l->{'EXPERIMENT_PACKAGE_SET'}->{'EXPERIMENT_PACKAGE'}->{'SAMPLE'}->{'DESCRIPTION'}->{'text'};
 	 $SAMPLE_DESCRIPTION=~s/[^[:ascii:]]+//g;
 	 my $SAMPLE_ATTRIBUTES='';
 	 #print "\n=== ".$l->{'EXPERIMENT_PACKAGE_SET'}->{'EXPERIMENT_PACKAGE'}->{'SAMPLE'};
 	 if(exists ${$l->{'EXPERIMENT_PACKAGE_SET'}->{'EXPERIMENT_PACKAGE'}->{'SAMPLE'}}{'SAMPLE_ATTRIBUTES'}){
 	  if(ref $l->{'EXPERIMENT_PACKAGE_SET'}->{'EXPERIMENT_PACKAGE'}->{'SAMPLE'}->{'SAMPLE_ATTRIBUTES'}->{'SAMPLE_ATTRIBUTE'} eq 'ARRAY'){
 	   foreach my $k(@{$l->{'EXPERIMENT_PACKAGE_SET'}->{'EXPERIMENT_PACKAGE'}->{'SAMPLE'}->{'SAMPLE_ATTRIBUTES'}->{'SAMPLE_ATTRIBUTE'}}){
 	   ;
 	    $SAMPLE_ATTRIBUTES.=$k->{'TAG'}->{'text'}.':'.$k->{'VALUE'}->{'text'}.';';
 	    }
 	  }
 	  else{
 	   $SAMPLE_ATTRIBUTES=$l->{'EXPERIMENT_PACKAGE_SET'}->{'EXPERIMENT_PACKAGE'}->{'SAMPLE'}->{'SAMPLE_ATTRIBUTES'}->{'SAMPLE_ATTRIBUTE'}->{'TAG'}->{'text'}.':'.$l->{'EXPERIMENT_PACKAGE_SET'}->{'EXPERIMENT_PACKAGE'}->{'SAMPLE'}->{'SAMPLE_ATTRIBUTES'}->{'SAMPLE_ATTRIBUTE'}->{'VALUE'}->{'text'}.';';
 	  }
 	 }
 	 $SAMPLE_ATTRIBUTES=~s/[^[:ascii:]]+//g;
 	 
 	# my $TAXON_ID=
 	 
 	 my $SRR_id='';
 	 if( ref $l->{'EXPERIMENT_PACKAGE_SET'}->{'EXPERIMENT_PACKAGE'}->{'RUN_SET'}->{'RUN'} eq 'ARRAY'){ 	 
 	  foreach my $k(@{$l->{'EXPERIMENT_PACKAGE_SET'}->{'EXPERIMENT_PACKAGE'}->{'RUN_SET'}->{'RUN'}}){
  	   $SRR_id =$k->{'accession'};
 	   print "$srx_list[$i]\t$SRR_id\t$LIBRARY_STRATEGY\t$LIBRARY_SOURCE\t$LIBRARY_SELECTION\t$SAMPLE_DESCRIPTOR\t$LIBRARY_DESCRIPTOR\t$SAMPLE_DESCRIPTION\t$SAMPLE_TITLE\t$SAMPLE_ATTRIBUTES\t$LIBRARY_LAYOUT\t$PLATFORM\t$STUDY_TITLE\t$STUDY_ABSTRACT\t$PUBMED\n";
 	   }
 	  }
 	 else{
 	  $SRR_id =$l->{'EXPERIMENT_PACKAGE_SET'}->{'EXPERIMENT_PACKAGE'}->{'RUN_SET'}->{'RUN'}->{'accession'};
 	  print "$srx_list[$i]\t$SRR_id\t$LIBRARY_STRATEGY\t$LIBRARY_SOURCE\t$LIBRARY_SELECTION\t$SAMPLE_DESCRIPTOR\t$LIBRARY_DESCRIPTOR\t$SAMPLE_DESCRIPTION\t$SAMPLE_TITLE\t$SAMPLE_ATTRIBUTES\t$LIBRARY_LAYOUT\t$PLATFORM\t$STUDY_TITLE\t$STUDY_ABSTRACT\t$PUBMED\n";
 	  }
 	   $SRR_id=~s/[^[:ascii:]]+//g; 	 
	}
	else{print STDERR "$output ... No ID found $srx_list[$i]... SKIPPED.";}
	#my $xml_hash = $xml_converter->fromXMLStringtoHash($output);
	#print Dumper( $xml_hash);
	print STDERR "-DONE-\n";
	sleep(1) if $i%5;
}




=l
#foreach my $i(0..$#srx_list){
#	my $output = get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term='.$srx_list[$i].'[accn]');
	#110\n2647205\n       SRX1853183[accn]    accn    1    N      GROUP  SRX1853183[accn]
#	if($output){my$l=$xml_converter->fromXMLStringtoHash($output); $id_list[$i]=$l->{'eSearchResult'}->{'IdList'}->{'Id'}->{'text'};}
#	else{$id_list[$i]=0; warn "$output ... No ID found $srx_list[$i]...\n exit. \n\n\n"}
#}
#my $id_list=join(",",@id_list);
#$id_list =~s/,0,/,/g; ## replcae the missing ids. careful !!!
#print "$id_list\n";

my $base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
#my $url = $base . "esummary.fcgi?db=$db&id=$id_list";

#post the esummary URL
#my $output = get $url;
#unless(length($output)>50){die "$output ... Fetched nothing ...\n exit. \n\n\n"}

#open(F,">test_srx.xml") or die "$!";
#print F $output;
#close F; 

my $output='';
open(F,"test_srx.xml") or die "$!";
$output.=$_ while(<F>);
close F; 

my $xml_hash = $xml_converter->fromXMLStringtoHash($output);
print Dumper( $xml_hash);


=cut












# Download protein records corresponding to a list of GI numbers.
#SAMN00009684	SRX017401	SRR037374
#SAMN05258394	SRX1853183	SRR3680368
#SAMN05258394	SRX1853183	SRR3680369
#SAMN05258394	SRX1853183	SRR3680370


=biosample
my $db = 'biosample';
my @id_list = ('SAMN00009684','SAMN05258394');
foreach my $i(0..$#id_list){$id_list[$i]=~s/SAMN[0]*//;$id_list[$i]=int($id_list[$i])}
my $id_list=join(",",@id_list);

my $base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
my $url = $base . "esummary.fcgi?db=$db&id=$id_list";

#post the esummary URL
#my $output = get $url;
#unless(length($output)>50){die "$output ... Fetched nothing ...\n exit. \n\n\n"}


#open(F,">test.xml") or die "$!";
#print F $output;
#close F; 

my $output='';
open(F,"test.xml") or die "$!";
$output.=$_ while(<F>);
close F; 
# Convertion from a XML String to a Hash
my $xml_hash = $xml_converter->fromXMLStringtoHash($output);

#print Dumper( $xml_hash);
my $dat=$xml_hash->{'eSummaryResult'}->{'DocumentSummarySet'}->{'DocumentSummary'};

                                                      
foreach my $sum(@{$dat}){
#	print Dumper($sum);
	my $bios_smp_id= $sum->{'Accession'}->{'text'};
	my $samp_dat_hash=$xml_converter->fromXMLStringtoHash($sum->{'SampleData'}->{'text'});
#	print Dumper $samp_dat_hash;
	my $title = $samp_dat_hash->{'BioSample'}->{'Description'}->{'Title'}->{'text'};
	print "$bios_smp_id\t$title\n";
	
}
=cut



=test
#assemble the epost URL
$base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
$url = $base . "epost.fcgi?db=$db&id=$id_list";

#post the epost URL
$output = get($url);

#parse WebEnv and QueryKey
$web = $1 if ($output =~ /<WebEnv>(\S+)<\/WebEnv>/);
$key = $1 if ($output =~ /<QueryKey>(\d+)<\/QueryKey>/);

### include this code for EPost-ESummary
#assemble the esummary URL
$url = $base . "esummary.fcgi?db=$db&query_key=$key&WebEnv=$web";

#post the esummary URL
$docsums = get($url);
print "$docsums";

### include this code for EPost-EFetch
#assemble the efetch URL
$url = $base . "efetch.fcgi?db=$db&query_key=$key&WebEnv=$web";
$url .= "&rettype=fasta&retmode=text";

#post the efetch URL
$data = get($url);
print "$data";
=cut
