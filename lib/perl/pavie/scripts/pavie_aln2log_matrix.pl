#!/usr/bin/env perl
#This script reads alignments and outputs a log odd matrix 

use Env qw(HOST);
use Env qw(HOME);
use Env qw(USER);
#Specific Files and Default Parameters
# VersionTag  0.01

#Input:
#Input the alignments and count the substitutions/Consevations
#Format:
#name: sequence
#name: sequence
#Empty line
#Note: Cna handle multiple alignments
#Output:
#alphabet
#symbol1 symbol2 log odd
#...

while ( <>)  
  {
    if (/al_seqEnd/){$count=1;}
    elsif (/al_seqStart/){$n=0;}
    
    elsif ( /al_seqStart/){;}
    elsif ( /al_seq/)
      {
	/[\S]*\s(.*)/;
	$sl[$n]=$1;

	$seq=$sl[$n]; 
	
	@l=($seq=~/([^-])/g);
	$tot_symbol+=$#l+1;
	$n++;
	$count=0;
      }
    else 
      {
	next;
      }


#1 Measure The % id to weight aln contribution
    if ( $count)
      {
	$id=0;$tot=0;
	for ($a=0;$a<$n-1; $a++)
	  {
	    @seq1=($sl[$a  ]=~/(.)/g);
	    for ($b=$a+1; $b<$n; $b++)
	      {
		@seq2=($sl[$b  ]=~/(.)/g);
		for  ($c=0;$c<=$#seq2;$c++)
		  {
		    $s1=$seq1[$c];$s2=$seq2[$c];
		    $id+=($s1 eq $s2 && $s1 ne '-')?1:0;
		    $tot+=($s1 ne '-' && $s2 ne '-')?1:0;
		  }
	      }
	  }
	if ($tot==0){next;}

	$w=int(($id*100)/$tot);
	$tot_symbol*=$w;
	
#2 Count the vrious pairs in the alignment
	
	for ($a=0; $a<$n-1; $a++)
	  {
	    @seq1=($sl[$a  ]=~/(.)/g);
	    for ($b=$a+1; $b<$n; $b++)
	      {
		@seq2=($sl[$b  ]=~/(.)/g);
		for  ($c=0;$c<=$#seq2;$c++)
		  {
		    $s1=$seq1[$c];$s2=$seq2[$c];
		    if ( $s1 ne '-'){$count_mat{'scount'}{$s1}+=$w;}
		    if ( $s2 ne '-'){$count_mat{'scount'}{$s2}+=$w;}
		    
		    if ( $s1 ne '-' && $s2 ne '-')
		      {$npairs+=$w;
		       $count_mat{$s1}{$s2}+=$w;
		     }
		  }
	      }
	  }
      }
  }


#3 compute and Output the matrix
@symbol_list=keys %{$count_mat{'scount'}};
print "Symbol_list: @symbol_list\n";

$tot2=$tot_symbol*$tot_symbol;
for ($a=0; $a<=$#symbol_list; $a++)
  {
    for ($b=$a; $b<=$#symbol_list; $b++)
      {
	$s1=$symbol_list[$a];
	$s2=$symbol_list[$b];
#Philipp Distances:
	
	$r1=($count_mat{$s1}{$s2}+1)/($count_mat{$s1}{$s1}+1);
	$r2=($count_mat{$s2}{$s1}+1)/($count_mat{$s2}{$s2}+1);
	$value=(log ($r1)+log($r2))/2;
	printf "%s %s %8.3f\n",$s1,$s2,$value;
	$gep+=($s1 ne $s2)?$value:0;
	$ngep+=($s1 ne $s2)?1:0;
      }
  }
printf "GAP %.2f",($gep/$ngep);
