open(IN,"R1.cutadapt.paired.oligo_unhit");
while(<IN>) { chomp; $table{$_}=1; 
$toss=<IN>;
$toss=<IN>;
$toss=<IN>;
} close(IN);


open(IN," SRR1609188_R1_001001.out.paired.fastq");
while(<IN>) { chomp;
if(exists($table{$_}))
{
  print "";
$toss=<IN>;
print $toss;
$toss=<IN>;
print $toss;
$toss=<IN>;
print $toss;
}
else
{$toss=<IN>;
$toss=<IN>;
$toss=<IN>;
} } close(IN);
' > unhit_oligo-KOH_5
