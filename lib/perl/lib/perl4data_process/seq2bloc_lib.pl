#!/usr/bin/env perl
#This script reads an aln in fasta format and modifies names so thathere are no duplicates
#if two sequences have the same name, the second one is renamed name_1 and so on
use Env qw(HOST);
use Env qw(HOME);
use Env qw(USER);
#Specific Files and Default Parameters
# VersionTag  1.73
# CL: 
#-lib: output lib, 
#-stat: output stat
#-deg : use degenerated alphabet

foreach $a (@ARGV){$cl.=" $a";}

%hseq=read_fasta_seq($ARGV[0]);
$ktup=$ARGV[1];

foreach $s (keys %hseq)
  {
    $seq=$hseq{$s}{seq};
   
    if ( ($cl=~/-deg/)){$seq=simplify_seq($seq)};
       
    $l=(length ($seq)-$ktup)+1;
    for ( $a=0; $a<$l; $a++)
      {
	$m=substr ($seq, $a, $ktup);
	
	$lu{$m}.="sep__sep$a\sep2__sep2$s";
	$count{$m}{N}++;
	if (!$count{$m}{$s})
	  {
	    $count{$m}{UN}++;
	    $count{$m}{$s}=1;
	  }
	else 
	  {
	    $count {$m}{UN}--;
	  }
      }
  }


foreach $k (keys (%lu))
  {
    
    $l=$lu{$k};
  
    @list=split (/sep__sep/, $l);shift (@list);
    $n=$#list+1;
    
    for ($i=0; $i<$n-1; $i++)
      {
	($p_i, $seq_i)=split /sep2__sep2/, $list[$i];
	
	for ($j=$i+1; $j<$n; $j++)
	  {
	    ($p_j, $seq_j)=split /sep2__sep2/, $list[$j];
	    
	    for ($k=0; $k<$ktup; $k++)
	      {
		$pos1=$p_i+$k+1;
		$pos2=$p_j+$k+1;
		if ($mat{$seq_i}{$pos1}{$seq_j}{$pos2} ||$mat{$seq_j}{$pos2}{$seq_i}{$pos1}  ){next;}
		else
		  {
		    $mat{$seq_i}{$pos1}{$seq_j}{$pos2}=$mat{$seq_j}{$pos2}{$seq_i}{$pos1}=1;
		    $ml{$seq_i}{$seq_j}.="\t$pos1 $pos2 10 0 0\n";
		    $count {$seq_i}{$seq_j}++;
		    $count {$seq_j}{$seq_i}++;
		    $count {$seq_i}{TOT}++;
		    $count {$seq_j}{TOT}++;
		    $tot++;
		  }
	    }
	  }
      }
  }

$output_lib=0;
@sl=keys(%hseq);
$nseq=$#sl+1;
if ( $cl=~/lib/)
  {
    print "! TC_LIB_FORMAT_01\n";
    print "$nseq\n";
    for ($a=0; $a<$nseq; $a++)
      {
	$seq=$hseq{$sl[$a]}{seq};
	$len=length($seq);
	print "$sl[$a] $len $seq\n";
      }
    
    for ($a=0; $a<$nseq-1; $a++)
      {
	for ($b=$a+1; $b<$nseq; $b++)
	  {
	    $i=$a+1;$j=$b+1;
	    print "#$i $j\n";
	    $seq_i=$sl[$a];
	    $seq_j=$sl[$b];
	    print "$ml{$seq_i}{$seq_j}";
	    
	  }
      }
    print "! CPU 0\n";
    print "! SEQ_1_TO_N\n";
  }

if ($cl =~/stat/)
  {
    foreach $s (keys(%hseq))
      {
	print "\n$s $count{$s}{TOT}";
      }
  }

sub read_fasta_seq 
  {
    my $f=$_[0];
    my %hseq;
    my (@seq, @com, @name);
    my ($a, $s,$nseq);

    open (F, $f);
    while (<F>)
      {
	$s.=$_;
      }
    close (F);
    
    
    @name=($s=~/>(.*).*\n[^>]*/g);
    
    @seq =($s=~/>.*.*\n([^>]*)/g);
    @com =($s=~/>.*(.*)\n([^>]*)/g);


    $nseq=$#name+1;
  
    for ($a=0; $a<$nseq; $a++)
      {
	my $n=$name[$a];
	my $s;
	$hseq{$n}{name}=$n;
	$s=$seq[$a];$s=~s/\s//g;
	
	$hseq{$n}{seq}=$s;
	$hseq{$n}{com}=$com[$a];
      }
    return %hseq;
  }

sub simplify_seq
  {
    my $s=$_[0];
    
    $s=~s/[avilmAVILM]/a/g;
    $s=~s/[dekrDEKR]/d/g;
    $s=~s/[wfyWFY]/w/g;
    return $s;
  }
	
