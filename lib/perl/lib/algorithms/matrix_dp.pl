#!/usr/bin/env perl

use Env qw(HOST);
use Env qw(HOME);
use Env qw(USER);
#Specific Files and Default Parameters
# VersionTag  1.00


@seq1=read_fasta($ARGV[0]);
@seq2=read_fasta($ARGV[1]);




for ( $a=0; $a< $l1[0]; $a++)
  {
    for ( $b=0; $b< $l1[1]; $b++)
      {
	for ( $c=0; $c< $l2[0]; $c++)
	  {
	    for ( $d=0; $d< $l2[1]; $d++)
	      {
		matrix_dp (%mat1,$a, $b, %mat2, $c, $d);
	      }
	  }
      }
  }

twomat2onemat
  {
    
    my %mat=@_[0];
    
    
    for ($i=0,$ni=$StartX;$ni<$EndX; $i++, $ni++)
      {
	for ($j=0,$nj=$StartY;$nj<$EndY; $j++, $nj++)
	  {
	    $new_mat{$i}{$j}{0}=$mat{0}{$ni}{$nj};
	    $new_mat{$i}{$j}{1}=$mat{1}{$ni-$mat0{SizeX}}{$nj--$mat0{SizeY}};
	  }
      }
    
		


while (<>)
  {
    if ( /^>(\S+)/)
      {;}
    else 
      {
	@list=/(\S)/g;

	foreach $s (@list)
	  {

	    $freq{$s}++;
	    $tot++;
	  }
      }
  }

@symbol_list=keys %freq;
foreach $s (@symbol_list)
  {
    printf "%s\t%.2f\n", $s, $freq{$s}/$tot;
  }

      
