#!/usr/bin/env perl
#This script reads an aln in fasta format and modifies names so thathere are no duplicates
#if two sequences have the same name, the second one is renamed name_1 and so on
use Env qw(HOST);
use Env qw(HOME);
use Env qw(USER);
#Specific Files and Default Parameters
# VersionTag  1.73
#mode=1 Compares lib1 with lib2
#mode=2 Compares lib1 with elib2
#mode=3 Outputs  elib2 for the sequences in elib1
$lib1=$ARGV[0];
$lib2=$ARGV[1];
$filter1=-1;
$filter2=-1;
foreach $a (@ARGV){$cl.=" $a";}
if (($cl=~/-mode (\S+)/)){$mode=$1;}
else {$mode=1;}
if (($cl=~/-filter1 (\d+)/)){$filter1=$1;}
else {$filter1=0;}

if (($cl=~/-filter2 (\d+)/)){$filter2=$1;}
else {$filter2=0;}

%lib1=read_tc_lib ($lib1, $filter1);
%lib2=read_tc_lib ($lib2, $filter2);

$n1=$lib1{size};
$n2=$lib2{size};
$n=0;







if ( $mode==1)
  {
    my ($n,$m);
    foreach $i (keys(%lib2))
      {
	
	if (!($i=~/##/)){next;}
	foreach $j(keys(%{$lib2{$i}}))
	  {
	    if (!($j=~/##/)){next;}
	    if ($lib2{$i}{$j}>$filter2)
	      {
		if ($lib1{$i}{$j}){$tp++;}
		else {$fp++;}
	      }
	    elsif ($lib2{$i}{$j}<$filter2)
	      {
		if ($lib1{$i}{$j}){$fn++;}
		else {$tn++;}
	      }
	  }
      }
    
    $sn=$tp/($tp+$fn);
    $sp=$tn/($tn+$fp);
    $co=$tp/($n1);
    print  "TP: $tp FN: $fn\n";
    printf "SN %.2f SP %.2f CO: %.2f", $sn, $sp, $co;
    exit (0);
  }
elsif ($mode="auc")
  {
    my ($i, $j, %h, @k, @tpr, @fpr);
    foreach $i (keys(%lib2))
      {
	if (!($i=~/##/)){next;}
	foreach $j(keys(%{$lib2{$i}}))
	  {
	    $h{$nh}{score}= $lib2{$i}{$j};
	    $h{$nh}{stat }=($lib1{$i}{$j})?1:0;
	    $nh++;
	    if ($lib1{$i}{$j}){$pp++;}
	    else {$pn++;}
	  }
      }
    @k=(sort {$h{$a}{score} <=> $h{$b}{score}} keys %h);
    
    
    foreach $k (@k)
      {
	if ($h{$k}{stat}==0){$fp++;}
	else {$tp++;}
	

	$tpr=$tp/$lib2{size};
	$fpr=$fp/$lib2{size};
	#printf "%.2f %.2f\n", $tpr, $fpr;
	$y1=$tpr;
	$x1=$fpr;
	$height=(($y2 + $y1)/2);
	$base = abs($x2 - $x1);
	$AUC += $base*$height;
	$y2 = $y1;
	$x2 = $x1;
      }
    print "AUC $AUC\n";
    exit (0);
  }
elsif ( $mode ==2)
  {
    foreach $i (%lib1)
      {
	if (!($i=~/##/)){next;}
	foreach $j(%{$lib1{$i}})
	  {
	    if (!($j=~/##/)){next;}
	    @l1=keys (%{$lib2{$i}});
	    @l2=keys (%{$lib2{$j}});
	    $s=0;
	    foreach $m1 (@l1)
	      {
		foreach $m2 (@l2)
		  {

		    $s+=($m1 eq $m2)?1:0;
		  }
	      }
	    $n+=($s>0)?1:0;
	    
	  }
      }
  }
elsif ($mode ==3)
  {
    $seq1=$lib1{index}{1};
    $seq2=$lib1{index}{2};
    $len1=$lib1{len}{$seq1};
    $len2=$lib1{len}{$seq2};
    print "! TC_LIB_FORMAT_01\n";
    print "2\n";
    print "$seq1 $len1 $lib1{seq}{$seq1}\n";
    print "$seq2 $len2 $lib1{seq}{$seq2}\n";
    print "#1 2\n";
    
    for ($i=1; $i<=$len1; $i++)
      {
	$k1="##_$seq1\_$i";
	for ($j=1; $j<=$len2; $j++)
	  {
	    $k2="##_$seq2\_$j";
	  
	    @l1=keys (%{$lib2{$k1}});
	    @l2=keys (%{$lib2{$k2}});
	  
	    $s=0;
	    foreach $m1 (@l1)
	      {
		foreach $m2 (@l2)
		  {	
		    $s+=($m1 eq $m2)?1:0;
		  }
	      }
	    if ($s) {print "\t$i $j $s\n";}
	  }
      }
    print "! SEQ_1_TO_N\n";
    die;
  }
$n/=2;$n1/=2;$n2/=2;

$sn=$n/$n1;
$sp=$n/$n2;

printf "\nSIZE: $lib1 $n1\n";
printf "\nSIZE: $lib2 $n2\n";

printf "COMP: SN %.2f SP %.2f", $sn, $sp;






sub read_tc_lib 
  {
    my ($f,$filter)=@_;
    my ($line, $nseq,$cl, $k, %index, @l, %lib, $seq1, $seq2);
    open (F, "$f");

    while (<F>)
      {
	$line=$_;
	
	if ( /^!/){;}
	elsif (/^\d+\n/){;}
	elsif (!$cl && $line=~/^(\S+) (\d+) (\S+)/)
	  {
	 
	    $nseq++;
	    $lib{name}{$1}=$nseq;
	    $lib{index}{$nseq}=$1;
	    $lib{len}{$1}=$2;
	    $lib{seq}{$1}=$3;
	  }
	elsif (/^#/)
	  {
	    
	    $cl=1;
	    ($i, $j)=($line=~/\#\s*(\d+)\s+(\d+)/);
	    $seq1=$lib{index}{$i};
	    $seq2=$lib{index}{$j};
	    
	  }
	elsif ( @l=($line=~/(\d+)/g))
	  {
	    $k1="##_$seq1\_$l[0]";
	    $k2="##_$seq2\_$l[1]";
	    $lib{$k1}{$k2}=$lib{$k2}{$k1}=$l[2];
	    $lib{size}+=2;
	  }
      }
    
  return %lib;
  }

