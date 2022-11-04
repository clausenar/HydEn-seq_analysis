#!/usr/bin/perl -w

use strict;
use File::Basename;

####### Programs and diretories #########

####### Programs and directories #########
my $dirname = dirname(__FILE__); #get the script directory 

my $bowtie2bedgraph                   = glob($dirname."/bin/bowtie2bedgraph-t1.pl");
my $bedGraphToBigWig                  = glob($dirname."/bin/bedGraphToBigWig-linux");
my $rectify_trimmed_pairs             = glob($dirname."/bin/rectify_trimmed_pairs.pl");
my $bin                               = glob($dirname."/bin");
my $a1                                = glob($dirname."/adapters/forAdapters.tsv");
my $a2                                = glob($dirname."/adapters/revAdapters.tsv");
my $organisms                         = glob($dirname."/organisms/org.txt");
my $oligos                            = glob($dirname."/oligo");

# can maybe move this information to the repo
my $chrM_hg38                         = "/data3/liam/Homo_sapiens/UCSC/chrM/ChromInfo.txt";

######## Available data ################

my $adapter1        = `cut -f2 $a1 | sed 's/^/-a /' |tr '\012' ' '`;
my $adapter2        = `cut -f2 $a2 | sed 's/^/-a /' |tr '\012' ' '`;

my $orgs            = "";
my %index           = ();
my %chrSize         = ();

$index{"oligo"}     = $oligos;
$chrSize{"chrM"}    = $chrM_hg38;
my $map             = "uniq|multi";

open(ORG,$organisms); #reading std files
while(<ORG>){
	chomp($_);
	my @a           = split("\t",$_);
	$index{$a[0]}   = $a[1];
	$chrSize{$a[0]} = $a[2];
	$orgs          .= $a[0]."|";
}
close(ORG);

if(scalar(@ARGV) != 4){print STDERR "\n ***ERROR: Too many/few arguments\n". &help();exit;}


#### Arguments ############

my $org       = $ARGV[0];
my $mapping   = $ARGV[1];
my $FastqDir  = $ARGV[2];
my $OutDir    = $ARGV[3];

my $m         = ""; if($mapping eq "uniq"){$m = "-m1"}
my $algSuffix = "multi"; if($mapping eq "uniq"){$algSuffix = "uniq"}

my $tmp       = $OutDir."/tmp";
my $Log       = $OutDir."/Log";
my $Stats     = $OutDir."/Stats";
my $bg        = $OutDir."/bg";
my $bw        = $OutDir."/bw";
my $bwChrM    = $OutDir."/bwChrM";
my $bgChrM    = $OutDir."/bgChrM";

#### checking input variables ####

&checkArguments($org,$orgs);
&checkArguments($mapping,$map);
&checkFiles($a1);
&checkFiles($a2);
&checkFiles($organisms);
&checkFiles("$oligos/oligo.1.ebwt");

if(!-d $FastqDir){print STDERR "\n*** ERROR: $FastqDir does not exist\n".&help(); exit}    ## hacer una funcion para que cheque si existe o no


################ pipeline ##########

### creating directories if not created

`mkdir -p $OutDir $tmp $Log $Stats $bg $bw $bwChrM $bgChrM`;

### Looping over each fastq file

opendir(UN,$FastqDir) or die print "cannot open FastqDir:$!\n";
while (my $file = readdir(UN)) {
	next if($file !~ m/R1_001.fastq.gz$/);
	chomp($file);
    
	my $sample    = $file; $sample =~ s/_R1_001.fastq.gz//;   # names shouldn't have a _
	my $sampleOUT = $file; $sampleOUT =~ s/_.*//;   # names shouldn't have a _
    print $sampleOUT;
	my $log       = "$OutDir/Log/$sampleOUT"."_"."$algSuffix"."_"."$org.log";
	my $stats     = "$OutDir/Stats/$sampleOUT"."_"."$algSuffix"."_"."$org.stats.tsv";
    my $index_report     = "$OutDir/Stats/$sampleOUT"."_"."$algSuffix"."_"."$org.index_report.txt";

	open(LOG,">$log") or die print "cannot open $log:$!\n";

	#    print STATS "$sample\t";

	print STDERR "\n\n****** Processing: $sample\n\n";

	print LOG "#-----------------------------------------------------------\n".
	"#CMD: $org $mapping $FastqDir $OutDir on $sample\n".
	"#-----------------------------------------------------------\n\n";
    
	print LOG `date`;
	print LOG "#-----------------------------------------------------------\n\n";

	print LOG `\tbowtie-inspect -s $index{"oligo"}/oligo > $index_report\n\n`;
	`bowtie-inspect -s $index{"oligo"}/oligo > $index_report`;

	# removing adapter 
	print LOG "## remove adapter\n\n";
        
	###-a file:$dirname/adapters/forAdapters.fa  # an option for specifying a fasta file of sequences to remove

	# removing adapter 
	print LOG "## remove adapter\n\n";
        
	###-a file:$dirname/adapters/forAdapters.fa
	print LOG  "\tcutadapt $adapter1 -f fastq --match-read-wildcards --quiet -m 15 -q 10 $FastqDir/".$sample."_R1_001.fastq.gz > $tmp/$sample.R1.cutadapt\n\n";
	`cutadapt $adapter1 -f fastq --match-read-wildcards --quiet -m 15 -q 10 $FastqDir/"$sample"_R1_001.fastq.gz > $tmp/$sample.R1.cutadapt`;

	print LOG "\tcutadapt $adapter2 -f fastq --match-read-wildcards --quiet -m 15 -q 10 $FastqDir/".$sample."_R2_001.fastq.gz > $tmp/$sample.R2.cutadapt\n\n";
	`cutadapt $adapter2 -f fastq --match-read-wildcards --quiet -m 15 -q 10 $FastqDir/"$sample"_R2_001.fastq.gz > $tmp/$sample.R2.cutadapt`;

	# rectifying pairs
	print LOG "## rectifying pairs\n\n";
	print LOG "\t$rectify_trimmed_pairs $tmp/$sample.R1.cutadapt $tmp/$sample.R2.cutadapt\n\n";
	`$rectify_trimmed_pairs $tmp/$sample.R1.cutadapt $tmp/$sample.R2.cutadapt`;

	# mapping towards all oligos 
	print LOG "## mapping towards all oligos\n\n";

	print LOG "\tbowtie $m -v2 --max $tmp/$sample.R1.cutadapt.paired.oligo_max --un $tmp/$sample.R1.cutadapt.paired.oligo_unhit $index{\"oligo\"}/oligo $tmp/$sample.R1.cutadapt.paired.fastq $tmp/$sample.R1.cutadapt.paired.oligo_out 2>> $stats.1.tmp\n\n";
	`bowtie $m -v2 --max $tmp/$sample.R1.cutadapt.paired.oligo_max --un $tmp/$sample.R1.cutadapt.paired.oligo_unhit $index{"oligo"}/oligo $tmp/$sample.R1.cutadapt.paired.fastq $tmp/$sample.R1.cutadapt.paired.oligo_out 2>> $stats.1.tmp`;

	print LOG "\tbowtie $m -v2 --max $tmp/$sample.R1.cutadapt.paired.oligo_max --un $tmp/$sample.R1.cutadapt.paired.oligo_unhit $index{\"oligo\"}/oligo $tmp/$sample.R1.cutadapt.paired.fastq $tmp/$sample.R1.cutadapt.paired.oligo_out 2>> $stats.1.tmp\n\n";
	`bowtie $m -v2 --max $tmp/$sample.R1.cutadapt.paired.oligo_max --un $tmp/$sample.R1.cutadapt.paired.oligo_unhit $index{"oligo"}/oligo $tmp/$sample.R1.cutadapt.paired.fastq $tmp/$sample.R1.cutadapt.paired.oligo_out 2>> $stats.1.tmp`;

	# extracting read2 from unhit
	print LOG "extracting read2 from unhit\n\n";
	print LOG "\t&matchUnhit($tmp/$sample.R1.cutadapt.paired.oligo_unhit, $tmp/$sample.R2.cutadapt.paired.fastq, $tmp/$sample.R2.cutadapt.paired.unhit)\n\n";
	&matchUnhit("$tmp/$sample.R1.cutadapt.paired.oligo_unhit", "$tmp/$sample.R2.cutadapt.paired.fastq", "$tmp/$sample.R2.cutadapt.paired.unhit");

	# alignment towards organism    
	print LOG "## alignment towards organism\n\n";
	print LOG "\tbowtie $m -v2 -p8 -X2000 --un $tmp/$sample.unhit $index{$org} -1 $tmp/$sample.R1.cutadapt.paired.oligo_unhit -2 $tmp/$sample.R2.cutadapt.paired.unhit $tmp/$sample.pair 2>> $stats.2.tmp\n\n";
	`bowtie $m -v2 -p8 -X2000 --un $tmp/$sample.unhit  $index{$org} -1 $tmp/$sample.R1.cutadapt.paired.oligo_unhit -2 $tmp/$sample.R2.cutadapt.paired.unhit $tmp/$sample.pair 2>> $stats.2.tmp`;

	# extract mate1
	print LOG "## extract mate1\n\n";
	print LOG "\t&mate1($tmp/$sample.pair,$tmp/$sample.pair.mate1)\n\n";
	&mate1("$tmp/$sample.pair","$tmp/$sample.pair.mate1");

	# realign unhit to organism
	print LOG "## realign unhit to organism\n\n";
	print LOG "\tbowtie $m -v2 -p8 $index{$org} $tmp/".$sample."_1.unhit $tmp/".$sample."_1.unhit.out 2>>  $stats.3.tmp \n\n";
	`bowtie $m -v2 -p8 $index{$org} $tmp/"$sample"_1.unhit $tmp/"$sample"_1.unhit.out 2>>  $stats.3.tmp`; 

	# concatenating mate1 and realigned unhit
	print LOG "## concatenating mate1 and realigned unhit\n\n";
	print LOG "\tcat $tmp/$sample.pair.mate1 $tmp/".$sample."_1.unhit.out > $tmp/$sample.pair.mate1_1.unhit.out\n\n";
	`cat $tmp/$sample.pair.mate1 $tmp/"$sample"_1.unhit.out > $tmp/$sample.pair.mate1_1.unhit.out`;

	# generating bedgraphs y bigwig
	print LOG "generation bedgraphs and bigwig\n\n";
	print LOG "\t$bowtie2bedgraph $tmp/$sample.pair.mate1_1.unhit.out $bg/$algSuffix.$org.$sampleOUT\n\n";
	`$bowtie2bedgraph $tmp/$sample.pair.mate1_1.unhit.out $bg/$algSuffix.$org.$sampleOUT`;

	print LOG "\t&UCSC($bg/".$algSuffix.".".$org.".".$sampleOUT."_forward.bedgraph,$org)\n\n".
	"\t&UCSC($bg/".$algSuffix.".".$org.".".$sampleOUT."_reverse.bedgraph,$org)\n\n";
	&UCSC($bg."/".$algSuffix.".".$org.".".$sampleOUT."_forward.bedgraph",$org);
	&UCSC($bg."/".$algSuffix.".".$org.".".$sampleOUT."_reverse.bedgraph",$org);

	print LOG "\t$bedGraphToBigWig $bg/$algSuffix.$org.$sampleOUT"."_"."forward.fix.bedgraph $chrSize{$org} $bw/$algSuffix.$org.$sampleOUT"."_"."forward.fix.bw\n\n";
	`$bedGraphToBigWig $bg/$algSuffix.$org.$sampleOUT"_"forward.fix.bedgraph $chrSize{$org} $bw/$algSuffix.$org.$sampleOUT"_"forward.fix.bw`;
    
	print LOG "\t$bedGraphToBigWig $bg/$algSuffix.$org.$sampleOUT"."_"."reverse.fix.bedgraph $chrSize{$org} $bw/$algSuffix.$org.$sampleOUT"."_"."reverse.fix.bw\n\n";
	`$bedGraphToBigWig $bg/$algSuffix.$org.$sampleOUT"_"reverse.fix.bedgraph $chrSize{$org} $bw/$algSuffix.$org.$sampleOUT"_"reverse.fix.bw`;

	# generating bedgraphs and bigwig for human chrM
	print LOG "##generating bedgraphs and bigwig for human chrM\n\n";
	if($org eq "hg38"){
		print LOG "\tgrep \"chrM\" $bg/$algSuffix.$org.$sampleOUT"."_"."forward.bedgraph > $bgChrM/$algSuffix.$org.$sampleOUT.chrM_forward.bedgraph\n\n".
		"\tgrep \"chrM\" $bg/$algSuffix.$org.$sampleOUT"."_"."reverse.bedgraph > $bgChrM/$algSuffix.$org.$sampleOUT.chrM_reverse.bedgraph\n\n".
		"\t&UCSC($bgChrM/$algSuffix.$org.$sampleOUT.chrM_forward.bedgraph,$org)\n\n".
		"\t&UCSC($bgChrM/$algSuffix.$org.$sampleOUT.chrM_reverse.bedgraph,$org)\n\n".
		"\t$bedGraphToBigWig $bgChrM/$algSuffix.$org.$sampleOUT.chrM_forward.fix.bedgraph $chrSize{chrM} $bwChrM/$algSuffix.$org.$sampleOUT.chrM_forward.fix.bw\n\n".
		"\t$bedGraphToBigWig $bgChrM/$algSuffix.$org.$sampleOUT.chrM_reverse.fix.bedgraph $chrSize{chrM} $bwChrM/$algSuffix.$org.$sampleOUT.chrM_reverse.fix.bw\n\n";

		`grep "chrM" $bg/$algSuffix.$org.$sampleOUT"_"forward.bedgraph > $bgChrM/$algSuffix.$org.$sampleOUT.chrM_forward.bedgraph`; 
		`grep "chrM" $bg/$algSuffix.$org.$sampleOUT"_"reverse.bedgraph > $bgChrM/$algSuffix.$org.$sampleOUT.chrM_reverse.bedgraph`;
		&UCSC("$bgChrM/$algSuffix.$org.$sampleOUT.chrM_forward.bedgraph",$org);
		&UCSC("$bgChrM/$algSuffix.$org.$sampleOUT.chrM_reverse.bedgraph",$org);
		`$bedGraphToBigWig $bgChrM/$algSuffix.$org.$sampleOUT.chrM_forward.fix.bedgraph $chrSize{chrM} $bwChrM/$algSuffix.$org.$sampleOUT.chrM_forward.fix.bw`;
		`$bedGraphToBigWig $bgChrM/$algSuffix.$org.$sampleOUT.chrM_reverse.fix.bedgraph $chrSize{chrM} $bwChrM/$algSuffix.$org.$sampleOUT.chrM_reverse.fix.bw`;
	}

	print STDERR "Process $sampleOUT"."_"."$algSuffix"."_"."$org . . . . Done!\n\n";

	# adding info to the report

	print STDERR "Compiling report . . .\n\n";

	&report($sample, "$OutDir/Stats/$sampleOUT"."_"."$algSuffix"."_"."$org.stats.tsv");    

	print STDERR "Report printed at $OutDir/Stats/$sampleOUT"."_"."$algSuffix"."_"."$org.stats.tsv \n\n";
	print LOG "Report printed at $OutDir/Stats/$sampleOUT"."_"."$algSuffix"."_"."$org.stats.tsv \n\n";
	
	print LOG "#-----------------------------------------------------------\n";
	print LOG `date`;
	print LOG "#-----------------------------------------------------------\n";

	`rm -rv $tmp/$sample*`;
	#`rm -rv $stats*.tmp`;

}
close(LOG);
`cat $Stats/*.tsv >> $Stats/run_final_log.tsv`;
`rm -rv $Stats/*.tmp`;
`rm -rv $tmp`;


################  SUBS ##########################

sub report{
    my ($sample,$stats) = @_;
 
    open(REPORT, ">$stats") or die print "cannot open $stats:$!\n";
    
    my @header          = ("Library","Alignment","organism","Raw_Read_Count","Post-Trim_Read_Pairs","%Retained","Pairs_Aligned_to_Oligos","%Aligned_to_Oligos","Uniquely_Aligned_Pairs","%Uniquely_Aligned_as_Pairs","Pairs_Mapped_to_Multiple_Locations","%Mapped_to_Multiple_Locations_as_Pair","Uniquely_Mapped_Single_End","%Mapped_Uniquely_as_Single_End","Single_End_Mapped_to_Multiple_Locations","Mapped_to_Multiple_Locations_as_Single_End","Unmapped_Pairs","%Unmapped");
    foreach (@header){print REPORT "$_\t";} print REPORT "\n";

    my $raw             = (`gunzip -c $FastqDir/$sample"_R1_001.fastq.gz" | wc -l`/4); chomp($raw);
    my $trim            = (`wc -l $tmp/$sample.R1.cutadapt |cut -d ' ' -f1` / 4); chomp($trim); 
    my $retained        = sprintf ("%.2f%s",(100*$trim/$raw),"%");
  
    my $oligoMulti      = (`grep "reads with at least one reported alignment" $stats.1.tmp | sed 's/.*: //' |sed 's/ .*//'`); chomp($oligoMulti);
    my $oligoOne        = (`grep "reads with alignments suppressed due to" $stats.1.tmp | sed 's/.*: //' |sed 's/ .*//'`); chomp($oligoOne);
    my $oligoAln        = chomp($oligoOne)  + chomp($oligoMulti);
    my $oligoAlnRate    =  sprintf("%.2f%s",(100*$oligoAln/$trim),"%"); 

    my $pairedMulti     = (`grep "reads with at least one reported alignment" $stats.2.tmp | sed 's/.*: //' |sed 's/ .*//'`); chomp($pairedMulti); 
    my $pairedMultiRate = (`grep "reads with at least one reported alignment" $stats.2.tmp | sed 's/.*(//' |sed 's/).*//'`); chomp($pairedMultiRate);
    my $pairedOne       = (`grep "reads with alignments suppressed due to" $stats.2.tmp | sed 's/.*: //' |sed 's/ .*//'`); chomp($pairedOne); 
    my $pairedOneRate   = (`grep "reads with alignments suppressed due to" $stats.2.tmp | sed 's/.*(//' |sed 's/).*//'`); chomp($pairedOneRate); 

    my $singleMulti     = (`grep "reads with at least one reported alignment" $stats.3.tmp | sed 's/.*: //' |sed 's/ .*//'`); chomp($singleMulti); 
    my $singleMultiRate = (`grep "reads with at least one reported alignment" $stats.3.tmp | sed 's/.*(//' |sed 's/).*//'`); chomp($singleMultiRate);
    my $singleOne       = (`grep "reads with alignments suppressed due to" $stats.3.tmp | sed 's/.*: //' |sed 's/ .*//'`); chomp($singleOne); 
    my $singleOneRate   = (`grep "reads with alignments suppressed due to" $stats.3.tmp | sed 's/.*(//' |sed 's/).*//'`); chomp($singleOneRate); 

    my $allAln          = $oligoAln + $pairedOne + $pairedMulti + $singleOne + $singleMulti;
    my $allNotAln       = $trim - $allAln; 
    my $allNotAlnRate   = sprintf ("%.2f%s",(100*$allNotAln/$raw),"%");
	
	print REPORT "$sample\t$algSuffix\t$org\t$raw\t$trim\t$retained\t$oligoAln\t$oligoAlnRate\t$pairedMulti\t$pairedMultiRate\t$pairedOne\t$pairedOneRate\t$singleMulti\t$singleMultiRate\t$singleOne\t$singleOneRate\t$allNotAln\t$allNotAlnRate\n";
	close(REPORT);
}


sub UCSC{
	my %table = ();
	my @files = @_;
	my %chr   = ();

	open(IN, "$chrSize{$files[1]}") or die print "cannot open $chrSize{$files[1]}:$!\n";;
	while(<IN>){ 
		chomp; 
		my @G         = split(/\t/,$_); 
		$table{$G[0]} = $G[1];
		$chr{$G[0]}   = 1;	
	} 
	close(IN); 

	my $out = $files[0]; $out =~ s/bedgraph/fix.bedgraph/;

	open(OUT,">$out") or die print "cannot open $out!\n";; 
	open(IN,"$files[0]") or die print "cannot open $files[0]:$!\n";; 
	while(<IN>){ 
		chomp; 
		if($_ =~ m/^chr/) {
			my @F = split(/\t/,$_);
			if($chr{$F[0]} && $F[1]>0 && $F[2]<=$table{$F[0]} ){
				$F[1]--;
				print OUT "$F[0]\t$F[1]\t$F[2]\t$F[3]\n";
			}
		} 
	}
	close(IN);close(OUT);
}


sub mate1{
	my @files = @_;
	
	open(OUT,">$files[1]") or die print "cannot open $files[1]:$!\n";; 
	open(IN,"$files[0]") or die print "cannot open $files[0]:$!\n";; 
	while(<IN>){ 
		chomp; 
		if($_  =~ m/\/1\t/) { 
			print OUT "$_\n"; 
		} 
	}
	close(IN);
	close(OUT);
}


sub matchUnhit{
	my @files = @_;
	my %table = ();
	my $toss  = "";

	open(IN,"$files[0]") or die print "cannot open $files[0]:$!\n";; 
	while(<IN>) { 
		chomp; 
		$table{$_} = 1; 
		$toss      = <IN>; 
		$toss      = <IN>; 
		$toss      = <IN>; 
	}
	close(IN); 
    
	$toss="";

	open(OUT,">$files[2]") or die print "cannot open $files[2]:$!\n";; 
	open(IN,"$files[1]") or die print "cannot open $files[1]:$!\n";; 
	while(<IN>){ 
		chomp; 
		if(exists($table{$_})){ 
			print OUT "$_\n"; 
			my $toss = <IN>; print OUT $toss; 
			$toss = <IN>; print OUT  $toss; 
			$toss = <IN>; print OUT $toss; 
		} 
		else{ 
			$toss = <IN>; 
			$toss = <IN>; 
			$toss = <IN>; 
		} 
	} 
	close(IN);
	close(OUT);
}


sub adapters{
	my $file = $_[0];
	my %hash = ();
	
	open(IN,"$file") or die print "cannot open $file:$!\n";;
	while(<IN>){
		chomp($_);
		my @a = split("\t",$_);
		$hash{$a[0]} = $a[1];
	}
	close(IN);
	return(%hash);
}


sub checkArguments{
	my $list  = $_[1];
	my $check = $_[0];
	my @list  = split("\\|",$list);
	my $exist = 0;
    
	foreach(@list){if($_ eq $check){$exist=1;}}
		if($exist == 0){print STDERR "\n*** ERROR: $check is not a valid option\n".&help(); exit}
}


sub checkFiles{
	my $file = $_[0];
	if(!-e $file ){print STDERR "\n*** ERROR: $file does not exist\n".&help(); exit}

}

sub help {
	my $help = "\n\t usage: HydEn.pl [org] [alnReport] [FastqDir] [OutDir]

\t         org                $orgs
\t         alnReport          $map
\t         FastqDir           absolut path to the Fastq files
\t         OutDir             absolute path for the results


\tForward Adapters at $a1
\tReverse Adapters at $a2
\tOrganisms are at $bin/org.txt [org Bowtie_index chrSizes]
\tOligo Bowtie index at /$bin/oligo
\tchrM size at /home/anders/Index/chrM.hg38.chrom.sizes


\t         Run one level up from Fastq folder. HydEn will process all FASTQ files within the Fastq folder.\n\n";
	return $help;
}

=FINAL:

for i in *R1_001.fastq; do cutadapt -f fastq --match-read-wildcards -m 15 -q 10 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG $i > $i.R1.cutadapt;done
 
for i in *R2_001.fastq; do cutadapt -f fastq --match-read-wildcards -m 15 -q 10 -a AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT $i > $i.R2.cutadapt; done
 
Remove mates of reads trimmed shorter than 15nt (script is attached):

perl ~/Tools/rectify_trimmed_pairs R1.cutadapt R2.cutadapt

IÕm using the following to align to the oligos:

for i in *R1_001.paired.fastq; do bowtie -m1 -v2 --max $i.max_oligo --un $i.unhit_oligo ~/Index/oligo $i $i.out.oligo;done

You could use -m1 or -k1 for this, as the important thing is to have an unhit file that is free of reads alignable to any oligo.Ê Using -m1 would be a problem if I wasnÕt also redirecting the multi-mapped reads to the max file.

I then extract the end 2 reads to match those in the unhit file using the following:

perl -e 'open(IN,"mtDNA-KOH-150512_S7_L000_R1_001.paired.fastq.unhit_oligo"); while(<IN>) { chomp; $table{$_}=1; $toss=<IN>; $toss=<IN>; $toss=<IN>; } close(IN); open(IN," mtDNA-KOH-150512_S7_L000_R2_001.paired.fastq"); while(<IN>) { chomp; if(exists($table{$_})) { print "$_\n"; $toss=<IN>; print $toss; $toss=<IN>; print $toss; $toss=<IN>; print $toss; } else { $toss=<IN>; $toss=<IN>; $toss=<IN>; } } close(IN);' > mtDNA-KOH-150512_S7_L000_R2_001.paired.fastq_unhit

Perform paired end alignment, multimapped:

bowtie -v2 -p8 -X2000 --un unhit_KOH-150512 ~/Index/GCA_000001405.15_GRCh38_no_alt_analysis_set -1 mtDNA-KOH-150512_S7_L000_R1_001.paired.fastq.unhit_oligo -2 mtDNA-KOH-150512_S7_L000_R2_001.paired.fastq_unhit mtDNA-KOH-150512_S7_pair

Extract mate 1:

perl -ane 'print $_ if $F[0]=/\/1/;' mtDNA-KOH-150512_S7_pair > mtDNA-KOH-150512_S7_pair.mate1

Realign unhit read1:

bowtie -v2 -p8 ~/Index/GCA_000001405.15_GRCh38_no_alt_analysis_set unhit_KOH-150512_1 out-unhit_KOH-150512_1

Concatenate mate 1 and realigned unhit:

cat mtDNA-KOH-150512_S7_pair.mate1 out-unhit_KOH-150512_1 > KOHl-150512-cat

Make bedgraph:

perl ~/Tools/bowtie2bedgraph-t1.pl KCl-150512-cat t1.KCl-150512-cat

for i in *bedgraph; do grep "chrM" $i > "chrM".$i; done

Prior to converting the files, I create a ÒfixedÓ version of the bedgraphs to conform to UCSC standards (column 2 is a 0-based coordinate, column 3 is 1-based) and remove any entries that were shifted beyond the end of a chromosome:
Ê
for f in *.bedgraph; do (perl -ane 'BEGIN { open(IN, "L03.genome"); while(<IN>) { chomp; @G=split(/\t/,$_); $table{$G[0]}=$G[1]; } close(IN); } next if $_=~/track|^$/; $F[1]--; next if $F[1]<0; next if $F[2]>$table{$F[0]}; print "$F[0]\t$F[1]\t$F[2]\t$F[3]\n";' $f > ${f}.fix &); done
make executable:
chmod +x bedGraphToBigWig-linux

for i in *forward.bedgraph.fix; do ~/Index/bedGraphToBigWig-linux $i ~/Index/chrM.hg38.chrom.sizes $i.bw
