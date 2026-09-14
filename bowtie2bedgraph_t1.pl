#!/usr/bin/perl

use warnings;								#generate warnings
use strict 'vars';							#generate warnings when undefined variables are referenced
use threads;								#allow use of threads
use Getopt::Std;							#use get options package
use DBI;									#use Mysql integration

my @cols;												#initialize variables					
my (%plus_table,%minus_table,%plus_bin_table,%minus_bin_table,%merge_table,%binned_merge_table,%option,%length_table);
my ($key,$plus_thread,$minus_thread,$plus_bin_thread,$minus_bin_thread,$merge_thread,$binned_merge_thread,$temp);
my ($normalized,$binned,$merged,$binned_merged,$bin_size,$shift_val,$trim,$list,$shifted,$swap);
my $norm_value=0;
my $name;
my $dbh;
my $rh;
my @result;

getopts('o:b:s:t:l:r:x',\%option) || print_usage();				#extract options from the command line argument
									#allow -o, -b, and -s, if other arguments are passed, run print_usage subroutine
if(exists($option{o}))							#if the -o option is passed
	{
	if($option{o}=~/[^nbmv]/)					#check for characters other than n, b, v, or m in its argument
		{
		print_usage("o");					#if other characters are found, run the print_usage subroutine
		}
	$normalized=($option{o}=~/n/i) || 0;				#set flag values for each optional file time, 
	$binned=($option{o}=~/b/i) || 0;
	$merged=($option{o}=~/m/i) || 0;
	$binned_merged=($option{o}=~/v/i) || 0;
	}

print_usage("bz") if (exists($option{b}) && $option{b}<2);		#if shift value is less than 1 or bin size is less than 2, run the print_usage subroutine
print_usage("sz") if (exists($option{s}) && $option{s}<1);
print_usage("b") if ((!$binned && !$binned_merged) && exists($option{b}));			#if bin size is specified and binned output is not requested, run print_usage 
print_usage("s") if ((!$merged && !$binned_merged) && exists($option{s}));			#if shift value is specified and merged output is not requested, run print_usage

$bin_size=$option{b} || 25;								#set bin size to specified value, otherwise set it to 25
$shift_val=$option{s} || 75;							#if shift value is not specified, set to 75
$trim=$option{t} || 1;									#if trim value is not specified, set $trim to 0
$swap=$option{x} || 0;									#if -x option is not specified, set $swap to 0, otherwise 1

if(exists($option{r}))								
	{
	$list=$option{r};
	print_usage("r") if (!($merged || $binned || $binned_merged));										#if merged or binned output is not requested, and a reference genome is specified, run print_usage
	$dbh=DBI->connect("DBI:mysql:database=".$list.";host=genome-mysql.cse.ucsc.edu","genome");	#connect the the user-specified UCSC database
	print_usage("db") if $dbh->{Active}==0;
	$rh=$dbh->prepare("select chrom, size from chromInfo");								#prepare and execute a request for a list of chromosomes and their lengths
	$rh->execute();						
	while(@result=$rh->fetchrow_array())												#store each row of the list in the length_table hash
		{
		$length_table{$result[0]}=$result[1];
		}
	print_usage("db") if keys(%length_table)==0;
	}
elsif(exists($option{l}))
	{
	$list=$option{l};
	print_usage("l") if (!($merged || $binned));				#if merged or binned output is not requested, and a list file name is specified, run print_usage
	open(LISTFILE,$list) || print_usage("lf");	#open chromosome list file, run print_usage subroutine if file cannot be opened
	while(<LISTFILE>)							#split each line by whitespace and store in the length_table hash, with chromosome name as key and length as value
		{
		/(\S+)\s+(\S+)/;
		$length_table{$1}=$2;
		}
	close(LISTFILE);
	}
else
	{
	$list="";									#if list file name is not specified, assign blank value
	}									
					
open(INFILE,"$ARGV[0]") || print_usage("f");		#open input file, run print_usage subroutine if file cannot be opened

$name=$ARGV[1] || print_usage("t");			#check for output file previx in argument, run print_usage subroutine if not specified
$name=~s/\/?(\S*\/)*//;						#strip path from output file prefix

while(<INFILE>)								#for each line of the input file
	{
	$norm_value++;							#increment overall hit count
	@cols=/^[^\t]+\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t/;			#extract the second through fourth fields of the current line
	$cols[1]="chrM" if $cols[1]=~m/mito/i;				#
	if($cols[1]!~m/chr/i)						#if the chromosome field does not begin with 'chr', add it
		{
		$cols[1]="chr".$cols[1];
		}
	if($cols[0] eq '+')						#if hit is to forward strand
		{
		$cols[2]+=(1-$trim);					#adjust position from continuum starting at 0 to continuum starting at 1 and adjust for trimming
		$key="$cols[1]\t$cols[2]";				#concatenate the chromosome and location, and use it as a hash key, allows quick locating of previously hit locations
		$plus_table{$key}++;					#increase hit count for that chromosome and location
		if($binned)						#if binned output is requested
			{
			$temp=1 if ($temp=int($cols[2]/$bin_size)*$bin_size)<1;		#determine alignment's bin, use 1 if calculated value is zero
			$key="$cols[1]\t$temp";						#concatenate the chromosome and bin start location, and use it as a hash key
			$plus_bin_table{$key}++;					#increase bin hit count
			}
		$shifted=$cols[2]+$shift_val;
		}
	elsif($cols[0] eq '-')						#if hit is to reverse strand
		{
		$cols[2]+=(length($cols[3])+$trim); 			#adjust position by adding the length of the aligned read, now specifies 5' location, also add trim value
		$key="$cols[1]\t$cols[2]";						#perform same tasks as above
		$minus_table{$key}++;
		if($binned)
			{
			$temp=1 if ($temp=int($cols[2]/$bin_size)*$bin_size)<1;
			$key="$cols[1]\t$temp";
			$minus_bin_table{$key}++;
			}
		$shifted=$cols[2]-$shift_val;
		}
	if($merged)
		{
		$temp=$shifted;
		$temp=1 if $temp<1;
		if(exists($length_table{$cols[1]}))						
			{
			$temp=$length_table{$cols[1]} if $temp>$length_table{$cols[1]};				
			}
		elsif($list ne "")																
			{
			die "Error: chromosome $cols[1] could not be found in the chromosome lengths list\n";
			}
		$key="$cols[1]\t$temp";				#concatenate the chromosome and hit location, and use it as a hash key
		$merge_table{$key}++;					#increase the location hit count, do this for both forward and reverse strand hits
		}
	if($binned_merged)
		{
		$temp=$shifted;
		$temp=1 if ($temp=int($temp/$bin_size)*$bin_size)<1;
		if(exists($length_table{$cols[1]}))						
			{
			$temp=int($length_table{$cols[1]}/$bin_size)*$bin_size if $temp>$length_table{$cols[1]};				
			}
		elsif($list ne "")																
			{
			die "Error: chromosome $cols[1] could not be found in the chromosome lengths list\n";
			}
		$key="$cols[1]\t$temp";					#concatenate the chromosome and bin start location, and use it as a hash key				
		$binned_merge_table{$key}++;			#increase bin hit count, do this for both forward and reverse strand hits
		}
	}
close(INFILE);								#close input file

$norm_value/=1000000;							#divide overall hit count by one million, hit values will be normalized to millions of overall hits

if($swap==0) {
	$plus_thread=threads->create('sort_tables',\%plus_table,"forward");			#create separate threads for sorting the forward and reverse hits
	$minus_thread=threads->create('sort_tables',\%minus_table,"reverse");			#pass a reference to the appropriate hash along with a string to append to the output file
	$plus_bin_thread=threads->create('sort_tables_alt',\%plus_bin_table,"forward_binned") if $binned;
	$minus_bin_thread=threads->create('sort_tables_alt',\%minus_bin_table,"reverse_binned") if $binned;
}
else {
	$plus_thread=threads->create('sort_tables',\%minus_table,"forward");			#create separate threads for sorting the forward and reverse hits
	$minus_thread=threads->create('sort_tables',\%plus_table,"reverse");			#pass a reference to the appropriate hash along with a string to append to the output file
	$plus_bin_thread=threads->create('sort_tables_alt',\%minus_bin_table,"forward_binned") if $binned;
	$minus_bin_thread=threads->create('sort_tables_alt',\%plus_bin_table,"reverse_binned") if $binned;
}
$merge_thread=threads->create('sort_tables_alt',\%merge_table,"merged") if $merged;
$binned_merge_thread=threads->create('sort_tables_alt',\%binned_merge_table,"binned_merged") if $binned_merged;	

$plus_thread->join();							#join threads
$minus_thread->join();
$plus_bin_thread->join() if $binned;
$minus_bin_thread->join() if $binned;
$merge_thread->join() if $merged;
$binned_merge_thread->join() if $binned_merged;

sub sort_tables								#subroutine used by plus and minus threads for sorting
	{
	my $ref=shift;							#retrieve arguments passed to subroutine					
	my $strand=shift;
	my @result;
	my $value;
	open(OUTFILE,">$ARGV[1]_${strand}\.bedgraph") || die "Could not create output file: $ARGV[1]_${strand}\.bedgraph\n";			#create output file for standard hit count
	@result=sort 																#retrieve the hash keys, split them into chromosome and location components,
		{																#then sort first by chromosome name alphabetically, next by location numerically
		$a->[0] cmp $b->[0] || $a->[1] <=> $b->[1]											#store sorted keys in result array
		} map
			{
			[/(^.*)\t(.*$)/]
			} keys(%$ref);
	print OUTFILE "track type=bedGraph name=${name}_$strand description=${name}_$strand visibility=full color=179,27,27 altColor=179,27,27 priority=20\n\n";				#print track information to standard output file
	if($normalized)																						#if normalized output is requested
		{
		open(NORMFILE,">$ARGV[1]_norm1m_${strand}\.bedgraph") || die "Could not create output file: $ARGV[1]_norm1m_${strand}\.bedgraph\n";						#create normalized output file
		print NORMFILE "track type=bedGraph name=${name}_norm1m_${strand} description=${name}_norm1m_${strand} visibility=full color=179,27,27 altColor=179,27,27 priority=20\n\n";	#print track information to normalized file
		foreach (@result)						
			{							 
			$value=$ref->{"$_->[0]"."\t"."$_->[1]"};		#for each key in the result array, retrieve hit value by rejoining the key components
			print OUTFILE "$_->[0]\t$_->[1]\t".($_->[1]+1)."\t$value\n";	#print chromosome, location, location+1 (0-based half-open end), and hit value to output file
			$value/=$norm_value;					#divide hit value by normalization constant
			print NORMFILE "$_->[0]\t$_->[1]\t".($_->[1]+1)."\t$value\n";	#print chromosome, location, location+1 (0-based half-open end), and normalized hit value to output file
			}
		close(NORMFILE);						#close the normalized file
		}
	else									#if normalized output is not requested, only write to the standard output file
		{
		foreach (@result)						
			{							 
			$value=$ref->{"$_->[0]"."\t"."$_->[1]"};		#for each key in the result array, retrieve hit value by rejoining the key components
			print OUTFILE "$_->[0]\t$_->[1]\t".($_->[1]+1)."\t$value\n";	#print chromosome, location, location+1 (0-based half-open end), and hit value to output file
			}
		}
	close(OUTFILE);							#close the standard output file
	return;
	}
	
sub sort_tables_alt							#subroutine used by bin and merge threads for sorting
	{
	my $ref=shift;							#retrieve arguments passed to subroutine					
	my $strand=shift;
	my @result;
	my ($value,$end);
	open(OUTFILE,">$ARGV[1]_${strand}\.bedgraph") || die "Could not create output file: $ARGV[1]_${strand}\.bedgraph\n";			#create binned or merged output file
	@result=sort 																#retrieve the hash keys, split them into chromosome and location components,
		{																#then sort first by chromosome name alphabetically, next by location numerically
		$a->[0] cmp $b->[0] || $a->[1] <=> $b->[1]											#store sorted keys in result array
		} map
			{
			[/(^.*)\t(.*$)/]
			} keys(%$ref);
	print OUTFILE "track type=bedGraph name=${name}_$strand description=${name}_$strand visibility=full color=179,27,27 altColor=179,27,27 priority=20\n\n";				#print track information to file
	foreach (@result)						
		{							 
		$value=$ref->{"$_->[0]"."\t"."$_->[1]"};						#for each key in the result array, retrieve hit value by rejoining the key components
		if($strand=~/bin/)														#if outfile is binned
			{
			$end=($_->[1]==1 ? ($bin_size-1) : $_->[1]+($bin_size-1));			#calculate end position based on bin size
			}
		else																	#otherwise use start value + 1 as end value (0-based half-open end)
			{
			$end=$_->[1]+1;
			}
		if($list ne "")
			{
			die "Error: chromosome $_->[0] could not be found in the chromosome lengths list\n" if !exists($length_table{$_->[0]});
			$end=$length_table{$_->[0]} if $end>$length_table{$_->[0]};			#if end value exceeds the current chromosome length and a length list was provided, set it's value to the chromosome's end
			}
		print OUTFILE "$_->[0]\t$_->[1]\t$end\t$value\n";					#print chromosome, location, location, and hit value to output file
		}
	close(OUTFILE);							#close output files
	return;
	}
	
sub print_usage								#subroutine called when an error in command line input is detected
	{
	my $arg=shift;							#retrieve argument
	if($arg eq "b")							#if -b option is present and binned output is not requested
		{
		print "Error: binned output files must be requested (-o b or -o v) for -b option to be accepted\n";
		}
	elsif($arg eq "s")						#if -s option is present and merged output is not requested
		{
		print "Error: merged output file must be requested (-o m) for -s option to be accepted\n";
		}
	elsif($arg eq "o")						#if letters other than n, b, v, or m are included in the -o option's argument
		{
		$option{o}=~s/n|b|m|v//g;					#remove valid letters so invalid letters only may be printed
		print "Error: invalid output file(s) requested: $option{o}\n";
		}
	elsif($arg eq "bz")						#if bin size specified is less than 1
		{
		print "Error: bin size must be greater than 1\n";
		}
	elsif($arg eq "sz")						#if shift value specified is less than 1
		{
		print "Error: shift value must be greater than 0\n";
		}
	elsif($arg eq "f")						#if input file name is invalid
		{
		print "Error: could not open input file \"$ARGV[0]\"\n"
		}
	elsif($arg eq "t")						#if track name is not specified
		{
		print "Error: output file prefix must be specified\n"
		}							
	elsif($arg eq "l")
		{
		print "Error: binned or merged output must be requested for -l option to be accepted\n";
		}
	elsif($arg eq "r")
		{
		print "Error: binned or merged output must be requested for -r option to be accepted\n";
		}
	elsif($arg eq "db")
		{
		print "Error: could not fetch list of chromosome lengths for reference genome \"$list\"\n";
		}
	elsif($arg eq "lf")
		{
		print "Error: could not open chromosome list file \"$list\"\n";
		}							#print small description of proper command use and exit
	die "Usage: bowtie2bedgraph.pl [options] [input file] [output file prefix]\n\t-o [nbmv]\tspecify additional files to be generated: n=normalized,\n\t\t\tb=binned, m=merged, v=binned and merged\n\t-b [number>=2]\tspecify bin size, requires -o b or -o v\n\t-s [number>=1]\tspecify number of bps to shift ChIP-seq hits prior to\n\t\t\tmerging, requires -o m or -o v\n\t-l [file name]\tprovide list of chromosome names and lengths to prevent\n\t\t\thits from being shifted beyond the end of the\n\t\t\tchromosome, requires -o b, -o m, or -o v\n\t-r [genome]\trather than specifying a chromosome length list,\n\t\t\tspecify a UCSC reference genome identifier (e.g. mm9),\n\t\t\tto have chromosome lengths fetched automatically,\n\t\t\trequires -o b, -o m, or -o v\n\t-t [number>=1]\tspecify number of nt trimmed prior to alignment\n\t-x\t\trequest stranded output files be swapped,\n\t\t\tforward->reverse and reverse->forward\n";
	}#!/usr/bin/perl

use warnings;								#generate warnings
use strict 'vars';							#generate warnings when undefined variables are referenced
use threads;								#allow use of threads
use Getopt::Std;							#use get options package
use DBI;									#use Mysql integration

my @cols;												#initialize variables					
my (%plus_table,%minus_table,%plus_bin_table,%minus_bin_table,%merge_table,%binned_merge_table,%option,%length_table);
my ($key,$plus_thread,$minus_thread,$plus_bin_thread,$minus_bin_thread,$merge_thread,$binned_merge_thread,$temp);
my ($normalized,$binned,$merged,$binned_merged,$bin_size,$shift_val,$trim,$list,$shifted,$swap);
my $norm_value=0;
my $name;
my $dbh;
my $rh;
my @result;

getopts('o:b:s:t:l:r:x',\%option) || print_usage();				#extract options from the command line argument
									#allow -o, -b, and -s, if other arguments are passed, run print_usage subroutine
if(exists($option{o}))							#if the -o option is passed
	{
	if($option{o}=~/[^nbmv]/)					#check for characters other than n, b, v, or m in its argument
		{
		print_usage("o");					#if other characters are found, run the print_usage subroutine
		}
	$normalized=($option{o}=~/n/i) || 0;				#set flag values for each optional file time, 
	$binned=($option{o}=~/b/i) || 0;
	$merged=($option{o}=~/m/i) || 0;
	$binned_merged=($option{o}=~/v/i) || 0;
	}

print_usage("bz") if (exists($option{b}) && $option{b}<2);		#if shift value is less than 1 or bin size is less than 2, run the print_usage subroutine
print_usage("sz") if (exists($option{s}) && $option{s}<1);
print_usage("b") if ((!$binned && !$binned_merged) && exists($option{b}));			#if bin size is specified and binned output is not requested, run print_usage 
print_usage("s") if ((!$merged && !$binned_merged) && exists($option{s}));			#if shift value is specified and merged output is not requested, run print_usage

$bin_size=$option{b} || 25;								#set bin size to specified value, otherwise set it to 25
$shift_val=$option{s} || 75;							#if shift value is not specified, set to 75
$trim=$option{t} || 1;									#if trim value is not specified, set $trim to 0
$swap=$option{x} || 0;									#if -x option is not specified, set $swap to 0, otherwise 1

if(exists($option{r}))								
	{
	$list=$option{r};
	print_usage("r") if (!($merged || $binned || $binned_merged));										#if merged or binned output is not requested, and a reference genome is specified, run print_usage
	$dbh=DBI->connect("DBI:mysql:database=".$list.";host=genome-mysql.cse.ucsc.edu","genome");	#connect the the user-specified UCSC database
	print_usage("db") if $dbh->{Active}==0;
	$rh=$dbh->prepare("select chrom, size from chromInfo");								#prepare and execute a request for a list of chromosomes and their lengths
	$rh->execute();						
	while(@result=$rh->fetchrow_array())												#store each row of the list in the length_table hash
		{
		$length_table{$result[0]}=$result[1];
		}
	print_usage("db") if keys(%length_table)==0;
	}
elsif(exists($option{l}))
	{
	$list=$option{l};
	print_usage("l") if (!($merged || $binned));				#if merged or binned output is not requested, and a list file name is specified, run print_usage
	open(LISTFILE,$list) || print_usage("lf");	#open chromosome list file, run print_usage subroutine if file cannot be opened
	while(<LISTFILE>)							#split each line by whitespace and store in the length_table hash, with chromosome name as key and length as value
		{
		/(\S+)\s+(\S+)/;
		$length_table{$1}=$2;
		}
	close(LISTFILE);
	}
else
	{
	$list="";									#if list file name is not specified, assign blank value
	}									
					
open(INFILE,"$ARGV[0]") || print_usage("f");		#open input file, run print_usage subroutine if file cannot be opened

$name=$ARGV[1] || print_usage("t");			#check for output file previx in argument, run print_usage subroutine if not specified
$name=~s/\/?(\S*\/)*//;						#strip path from output file prefix

while(<INFILE>)								#for each line of the input file
	{
	$norm_value++;							#increment overall hit count
	@cols=/^[^\t]+\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t/;			#extract the second through fourth fields of the current line
	$cols[1]="chrM" if $cols[1]=~m/mito/i;				#
	if($cols[1]!~m/chr/i)						#if the chromosome field does not begin with 'chr', add it
		{
		$cols[1]="chr".$cols[1];
		}
	if($cols[0] eq '+')						#if hit is to forward strand
		{
		$cols[2]+=(1-$trim);					#adjust position from continuum starting at 0 to continuum starting at 1 and adjust for trimming
		$key="$cols[1]\t$cols[2]";				#concatenate the chromosome and location, and use it as a hash key, allows quick locating of previously hit locations
		$plus_table{$key}++;					#increase hit count for that chromosome and location
		if($binned)						#if binned output is requested
			{
			$temp=1 if ($temp=int($cols[2]/$bin_size)*$bin_size)<1;		#determine alignment's bin, use 1 if calculated value is zero
			$key="$cols[1]\t$temp";						#concatenate the chromosome and bin start location, and use it as a hash key
			$plus_bin_table{$key}++;					#increase bin hit count
			}
		$shifted=$cols[2]+$shift_val;
		}
	elsif($cols[0] eq '-')						#if hit is to reverse strand
		{
		$cols[2]+=(length($cols[3])+$trim); 			#adjust position by adding the length of the aligned read, now specifies 5' location, also add trim value
		$key="$cols[1]\t$cols[2]";						#perform same tasks as above
		$minus_table{$key}++;
		if($binned)
			{
			$temp=1 if ($temp=int($cols[2]/$bin_size)*$bin_size)<1;
			$key="$cols[1]\t$temp";
			$minus_bin_table{$key}++;
			}
		$shifted=$cols[2]-$shift_val;
		}
	if($merged)
		{
		$temp=$shifted;
		$temp=1 if $temp<1;
		if(exists($length_table{$cols[1]}))						
			{
			$temp=$length_table{$cols[1]} if $temp>$length_table{$cols[1]};				
			}
		elsif($list ne "")																
			{
			die "Error: chromosome $cols[1] could not be found in the chromosome lengths list\n";
			}
		$key="$cols[1]\t$temp";				#concatenate the chromosome and hit location, and use it as a hash key
		$merge_table{$key}++;					#increase the location hit count, do this for both forward and reverse strand hits
		}
	if($binned_merged)
		{
		$temp=$shifted;
		$temp=1 if ($temp=int($temp/$bin_size)*$bin_size)<1;
		if(exists($length_table{$cols[1]}))						
			{
			$temp=int($length_table{$cols[1]}/$bin_size)*$bin_size if $temp>$length_table{$cols[1]};				
			}
		elsif($list ne "")																
			{
			die "Error: chromosome $cols[1] could not be found in the chromosome lengths list\n";
			}
		$key="$cols[1]\t$temp";					#concatenate the chromosome and bin start location, and use it as a hash key				
		$binned_merge_table{$key}++;			#increase bin hit count, do this for both forward and reverse strand hits
		}
	}
close(INFILE);								#close input file

$norm_value/=1000000;							#divide overall hit count by one million, hit values will be normalized to millions of overall hits

if($swap==0) {
	$plus_thread=threads->create('sort_tables',\%plus_table,"forward");			#create separate threads for sorting the forward and reverse hits
	$minus_thread=threads->create('sort_tables',\%minus_table,"reverse");			#pass a reference to the appropriate hash along with a string to append to the output file
	$plus_bin_thread=threads->create('sort_tables_alt',\%plus_bin_table,"forward_binned") if $binned;
	$minus_bin_thread=threads->create('sort_tables_alt',\%minus_bin_table,"reverse_binned") if $binned;
}
else {
	$plus_thread=threads->create('sort_tables',\%minus_table,"forward");			#create separate threads for sorting the forward and reverse hits
	$minus_thread=threads->create('sort_tables',\%plus_table,"reverse");			#pass a reference to the appropriate hash along with a string to append to the output file
	$plus_bin_thread=threads->create('sort_tables_alt',\%minus_bin_table,"forward_binned") if $binned;
	$minus_bin_thread=threads->create('sort_tables_alt',\%plus_bin_table,"reverse_binned") if $binned;
}
$merge_thread=threads->create('sort_tables_alt',\%merge_table,"merged") if $merged;
$binned_merge_thread=threads->create('sort_tables_alt',\%binned_merge_table,"binned_merged") if $binned_merged;	

$plus_thread->join();							#join threads
$minus_thread->join();
$plus_bin_thread->join() if $binned;
$minus_bin_thread->join() if $binned;
$merge_thread->join() if $merged;
$binned_merge_thread->join() if $binned_merged;

sub sort_tables								#subroutine used by plus and minus threads for sorting
	{
	my $ref=shift;							#retrieve arguments passed to subroutine					
	my $strand=shift;
	my @result;
	my $value;
	open(OUTFILE,">$ARGV[1]_${strand}\.bedgraph") || die "Could not create output file: $ARGV[1]_${strand}\.bedgraph\n";			#create output file for standard hit count
	@result=sort 																#retrieve the hash keys, split them into chromosome and location components,
		{																#then sort first by chromosome name alphabetically, next by location numerically
		$a->[0] cmp $b->[0] || $a->[1] <=> $b->[1]											#store sorted keys in result array
		} map
			{
			[/(^.*)\t(.*$)/]
			} keys(%$ref);
	print OUTFILE "track type=bedGraph name=${name}_$strand description=${name}_$strand visibility=full color=179,27,27 altColor=179,27,27 priority=20\n\n";				#print track information to standard output file
	if($normalized)																						#if normalized output is requested
		{
		open(NORMFILE,">$ARGV[1]_norm1m_${strand}\.bedgraph") || die "Could not create output file: $ARGV[1]_norm1m_${strand}\.bedgraph\n";						#create normalized output file
		print NORMFILE "track type=bedGraph name=${name}_norm1m_${strand} description=${name}_norm1m_${strand} visibility=full color=179,27,27 altColor=179,27,27 priority=20\n\n";	#print track information to normalized file
		foreach (@result)						
			{							 
			$value=$ref->{"$_->[0]"."\t"."$_->[1]"};		#for each key in the result array, retrieve hit value by rejoining the key components
			print OUTFILE "$_->[0]\t$_->[1]\t".($_->[1]+1)."\t$value\n";	#print chromosome, location, location+1 (0-based half-open end), and hit value to output file
			$value/=$norm_value;					#divide hit value by normalization constant
			print NORMFILE "$_->[0]\t$_->[1]\t".($_->[1]+1)."\t$value\n";	#print chromosome, location, location+1 (0-based half-open end), and normalized hit value to output file
			}
		close(NORMFILE);						#close the normalized file
		}
	else									#if normalized output is not requested, only write to the standard output file
		{
		foreach (@result)						
			{							 
			$value=$ref->{"$_->[0]"."\t"."$_->[1]"};		#for each key in the result array, retrieve hit value by rejoining the key components
			print OUTFILE "$_->[0]\t$_->[1]\t".($_->[1]+1)."\t$value\n";	#print chromosome, location, location+1 (0-based half-open end), and hit value to output file
			}
		}
	close(OUTFILE);							#close the standard output file
	return;
	}
	
sub sort_tables_alt							#subroutine used by bin and merge threads for sorting
	{
	my $ref=shift;							#retrieve arguments passed to subroutine					
	my $strand=shift;
	my @result;
	my ($value,$end);
	open(OUTFILE,">$ARGV[1]_${strand}\.bedgraph") || die "Could not create output file: $ARGV[1]_${strand}\.bedgraph\n";			#create binned or merged output file
	@result=sort 																#retrieve the hash keys, split them into chromosome and location components,
		{																#then sort first by chromosome name alphabetically, next by location numerically
		$a->[0] cmp $b->[0] || $a->[1] <=> $b->[1]											#store sorted keys in result array
		} map
			{
			[/(^.*)\t(.*$)/]
			} keys(%$ref);
	print OUTFILE "track type=bedGraph name=${name}_$strand description=${name}_$strand visibility=full color=179,27,27 altColor=179,27,27 priority=20\n\n";				#print track information to file
	foreach (@result)						
		{							 
		$value=$ref->{"$_->[0]"."\t"."$_->[1]"};						#for each key in the result array, retrieve hit value by rejoining the key components
		if($strand=~/bin/)														#if outfile is binned
			{
			$end=($_->[1]==1 ? ($bin_size-1) : $_->[1]+($bin_size-1));			#calculate end position based on bin size
			}
		else																	#otherwise use start value + 1 as end value (0-based half-open end)
			{
			$end=$_->[1]+1;
			}
		if($list ne "")
			{
			die "Error: chromosome $_->[0] could not be found in the chromosome lengths list\n" if !exists($length_table{$_->[0]});
			$end=$length_table{$_->[0]} if $end>$length_table{$_->[0]};			#if end value exceeds the current chromosome length and a length list was provided, set it's value to the chromosome's end
			}
		print OUTFILE "$_->[0]\t$_->[1]\t$end\t$value\n";					#print chromosome, location, location, and hit value to output file
		}
	close(OUTFILE);							#close output files
	return;
	}
	
sub print_usage								#subroutine called when an error in command line input is detected
	{
	my $arg=shift;							#retrieve argument
	if($arg eq "b")							#if -b option is present and binned output is not requested
		{
		print "Error: binned output files must be requested (-o b or -o v) for -b option to be accepted\n";
		}
	elsif($arg eq "s")						#if -s option is present and merged output is not requested
		{
		print "Error: merged output file must be requested (-o m) for -s option to be accepted\n";
		}
	elsif($arg eq "o")						#if letters other than n, b, v, or m are included in the -o option's argument
		{
		$option{o}=~s/n|b|m|v//g;					#remove valid letters so invalid letters only may be printed
		print "Error: invalid output file(s) requested: $option{o}\n";
		}
	elsif($arg eq "bz")						#if bin size specified is less than 1
		{
		print "Error: bin size must be greater than 1\n";
		}
	elsif($arg eq "sz")						#if shift value specified is less than 1
		{
		print "Error: shift value must be greater than 0\n";
		}
	elsif($arg eq "f")						#if input file name is invalid
		{
		print "Error: could not open input file \"$ARGV[0]\"\n"
		}
	elsif($arg eq "t")						#if track name is not specified
		{
		print "Error: output file prefix must be specified\n"
		}							
	elsif($arg eq "l")
		{
		print "Error: binned or merged output must be requested for -l option to be accepted\n";
		}
	elsif($arg eq "r")
		{
		print "Error: binned or merged output must be requested for -r option to be accepted\n";
		}
	elsif($arg eq "db")
		{
		print "Error: could not fetch list of chromosome lengths for reference genome \"$list\"\n";
		}
	elsif($arg eq "lf")
		{
		print "Error: could not open chromosome list file \"$list\"\n";
		}							#print small description of proper command use and exit
	die "Usage: bowtie2bedgraph.pl [options] [input file] [output file prefix]\n\t-o [nbmv]\tspecify additional files to be generated: n=normalized,\n\t\t\tb=binned, m=merged, v=binned and merged\n\t-b [number>=2]\tspecify bin size, requires -o b or -o v\n\t-s [number>=1]\tspecify number of bps to shift ChIP-seq hits prior to\n\t\t\tmerging, requires -o m or -o v\n\t-l [file name]\tprovide list of chromosome names and lengths to prevent\n\t\t\thits from being shifted beyond the end of the\n\t\t\tchromosome, requires -o b, -o m, or -o v\n\t-r [genome]\trather than specifying a chromosome length list,\n\t\t\tspecify a UCSC reference genome identifier (e.g. mm9),\n\t\t\tto have chromosome lengths fetched automatically,\n\t\t\trequires -o b, -o m, or -o v\n\t-t [number>=1]\tspecify number of nt trimmed prior to alignment\n\t-x\t\trequest stranded output files be swapped,\n\t\t\tforward->reverse and reverse->forward\n";
	}