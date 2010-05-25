#!/usr/bin/env perl
use Getopt::Long;
$io="fasta_compound";
$T=50;
$minI=-1;
$maxI=10000;
$ncomp=0;
$dup="yes";
GetOptions("-in=s" => \$in, "-io=s" => \$io, "-T=s" => \$T, "-max=s" => \$maxI,"-min=s" => \$minI, "-aln=s" => \$aln, "-dup=s" => \$dup  );

if ($aln)
  {
    open (F, "$aln");
    %hseq=read_fasta_seq ($aln);
  }

open (F, "$in");
while (<F>)
  {
    $l=$_;
    
    if ( /compound/)
      {

	$l=~s/,/ /g;
	@values=($l=~/(\w+)/g);
	shift (@values);
	$ngene=$#values+1;
	for ($a=0; $a<$ngene; $a++)
	  {
	    $gene{$a}=$values[$a];
	    $gene{$values[$a]}=$a;
	  }
      }
    else
      {
	$l=~/\s*(\w.*\w)\s*,/;
	$c=$1;
	$c=~s/\s/_/g;
	
	$comp{$ncomp}=$c;
	$comp{$c}=$ncomp;
	
	$l=~s/.*,//;
	@values=($l=~/(\S+)/g);
	foreach ($b=0; $b<$ngene; $b++)
	  {
	    $h{$gene{$b}}{$comp{$ncomp}}=$values[$b];
	  }
	$ncomp++;
      }
  }
close (F);



if(!%hseq)
  {
    $nseq=$ngene;
    %seq=%gene;
  }
else
  {
    my ($p,$s);
    $nseq=0;
    foreach $s (keys(%hseq))
      {
	$p=$hseq{$s}{order};
	$seq{$p}=$s;
	$nseq++;
      }
  }
	  

#Count #I
for ($a=0; $a<$ncomp; $a++)
      {
	for ($b=0; $b<$ngene; $b++)
	  {
	    $h{$gene{$b}}{$comp{$a}}-=$h{$gene{$b}}{$comp{$ncomp-1}};
	    $v=$h{$gene{$b}}{$comp{$a}};
	    $v2=($v>$T)?'I':'O';
	    if ($v2 eq "I"){$nI{gene}{$gene{$b}}++;$nI{comp}{$comp{$a}}++;}
	  }
      }
if ($io eq "fasta_seq")
  {
   
    for ($a=0; $a<$nseq; $a++)
      {
	print "\n>$seq{$a}\n";
	for ($b=0; $b<$ncomp-1; $b++)
	  {
	    if ($nI{comp}{$comp{$b}}<=$minI || $nI{comp}{$comp{$b}}>=$maxI ){next;}
	    #$t=$h{$seq{$a}}{$comp{$ncomp-1}}*$T;
	    $t=$T;
	    $v=$h{$seq{$a}}{$comp{$b}};
	    $v2=($v>$T)?'I':'O';
	    print $v2; 
	  }
      }
  }
elsif ( $io eq "fasta_compound")
    {
      for ($a=0; $a<$ncomp-1; $a++)
	{
	  if ($nI{comp}{$comp{$a}}>$minI && $nI{comp}{$comp{$a}}<$maxI )
	    {
	      $string="";
	      for ($b=0; $b<$nseq; $b++)
		{
		  $t=$T;
		  $v=$h{$seq{$b}}{$comp{$a}};
		  $v2=($v>$T)?'I':'O';
		  $string.="$v2";
		}
	      if ( $duplicate {$string} && $dup eq "no"){;}
	      else
		{
		  print "\n>$comp{$a} $nI{comp}{$comp{$a}}\n";
		  print "$string";
		  $duplicate{$string}=1;
		}

	    }
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
	$hseq{$n}{order}=$a;
	
	$s=$seq[$a];$s=~s/\s//g;
	
	$hseq{$n}{seq}=$s;
	$hseq{$n}{com}=$com[$a];
      }
    
    return %hseq;
  }
