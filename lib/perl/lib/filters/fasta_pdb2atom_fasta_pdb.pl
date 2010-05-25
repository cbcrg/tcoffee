#!/usr/bin/env perl

$fasta=$ARGV[0];

open (F, $fasta);
$/='>';
$n=0;
<F>;
while ( <F>)
  {
    
    $entry=">$_";
    chop $entry;
    $entry=~/(\>.*)/g;$header=$1;
    $header=~/>(\S+)/;$name=$1;
    $name=~s/[_]//g;
    
    $seq=`t_coffee -other_pg extract_from_pdb -netfile $name -seq -seq_field ATOM | grep -v '>'`;
    print "$header\n$seq\n";
    
    if ( $n==5){die;}
    $n++;
  }


