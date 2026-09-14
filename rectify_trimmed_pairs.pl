#!/usr/bin/perl

use warnings;
use strict 'vars';

my $outname1;
my $outname2;
my $current1;
my $current2;
my $max;
my $id;
my $seq;
my $qual;
my $plus;
my @F;
my %table;

if(scalar(@ARGV)<2) {
	die "Usage: rectify_trimmed_pairs.pl [End 1 FASTQ] [End 2 FASTQ] [maximum length (optional)]\n";
}

if($ARGV[0]=~/(.*)\.fastq/) {
	$outname1="$1.paired.fastq";
}
else {
	$outname1="$ARGV[0].paired.fastq";
}

if($ARGV[1]=~/(.*)\.fastq/) {
	$outname2="$1.paired.fastq";
}
else {
	$outname2="$ARGV[1].paired.fastq";
}

$max=$ARGV[2] || 0;

open(ONE,$ARGV[0]) || die "Error: could not open \"$ARGV[0]\"\n";
open(TWO,$ARGV[1]) || die "Error: could not open \"$ARGV[1]\"\n";
open(OUT1,">$outname1") || die "Error: could not create output file \"$outname1\"\n";
open(OUT2,">$outname2") || die "Error: could not create output file \"$outname2\"\n";

while(<ONE>) {
	$id=$_;
	$seq=<ONE>;
	$plus=<ONE>;
	$qual=<ONE>;
	@F=split(/\s+/,$id);
	$F[0]=~s/\/1$//;
	if(length($seq)>($max+1) && $max!=0) {
		$seq=substr($seq,0,$max);
		$seq.="\n";
		$qual=substr($qual,0,$max);
		$qual.="\n";
	}
	else {
		$seq=substr($seq,0,-2);
		$seq.="\n";
		$qual=substr($qual,0,-2);
		$qual.="\n";
	}
	$table{$F[0]}=[$seq,$qual];
}
close(ONE);

while(<TWO>) {
	$id=$_;
	$seq=<TWO>;
	$plus=<TWO>;
	$qual=<TWO>;
	if(length($seq)>($max+1) && $max!=0) {
		$seq=substr($seq,0,$max);
		$seq.="\n";
		$qual=substr($qual,0,$max);
		$qual.="\n";
	}
	else {
		$seq=substr($seq,0,-2);
		$seq.="\n";
		$qual=substr($qual,0,-2);
		$qual.="\n";
	}
	@F=split(/\s+/,$id);
	$F[0]=~s/\/2$//;
	if(exists($table{$F[0]})) {
		print OUT1 "$F[0]\n".$table{$F[0]}->[0]."+\n".$table{$F[0]}->[1];
		print OUT2 "$F[0]\n".$seq."+\n".$qual;
	}
}
close(TWO);

close(OUT1);
close(OUT2);