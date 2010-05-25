#!/usr/bin/env perl
#Version 1.01 (25/02/03)
$upper_threshold=9;
$lower_threshold=0;
$char='-';
$output="fasta_aln";
$keep=0;
$output="triangle";
foreach ($np=0; $np<=$#ARGV; $np++)
    {
    $value=$ARGV[$np];

    if ($value eq "-al1")
      {
      $alphabet1= $ARGV[++$np];
      }
    elsif ($value eq "-al2")
      {
        $alphabet2= $ARGV[++$np];
      }
    elsif ($value eq "-mat1")
      {
        $mat1= $ARGV[++$np];
      }
   
    elsif ($value eq "-output")
      {
        $output= $ARGV[++$np];
      }
     elsif ($value eq "-outformat")
      {
        $outformat= $ARGV[++$np];
      }
    else
      {
	printf "convert_substitution_matrix -al1 alphabet1 -al2 -mat1 -mat2 -output square|triangle -outformat code|text";
      }
  }



@al1=($alphabet1=~/(.{1})/g);
@al2=($alphabet2=~/(.{1})/g);
open (F, $mat1);
$nline=0;
while (<F>)
  {
    $line=$_;
  
    @values=($line=~/([-+0123456789]+)[\s,.}]+/g);
    print "$#values\n";
    for ( $a=0; $a<=$#values; $a++)
      {
	$full_matrix{$al1[$nline]}{$al1[$a]}=$full_matrix{$al1[$a]}{$al1[$nline]}=$values[$a];
	      }
   $nline++;


  }


for ( $a=0; $a<=$#al2; $a++)
  {
    $max_b=($output eq "triangle")?$a:$#al2;
    for ( $b=0; $b<=$max_b; $b++)
      {
	
	$c=$full_matrix{$al2[$a]}{$al2[$b]};
	printf " %3d",$c;
	if ($b!=$max_b){print",";}
      }
    if ( $a==$#al2){print "};";}
    else {print ",";}
    print "\n";    
  }
   

    
    
	 
	 
