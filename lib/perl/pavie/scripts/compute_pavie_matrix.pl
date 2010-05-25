#!/usr/bin/env perl
#Takes a list of file where each file contains one sequence in fasta format
#Makes all the pairwise alignments
#outputs a matrix

@al=@ARGV;
for ( $a=0; $a<=$#al; $a++)
  {
    if ( $al[$a] eq "-seq")
      {
	$seq=$al[++$a];
      }
    elsif ($al[$a] eq "-in_mat")
      {
	$in_mat=$al[++$a];
      }
    elsif ($al[$a] eq "-out_mat")
      {
	$out_mat=$al[++$a];
      }
  }

if ( !$out_mat)
  {
    print " You must give a name for the output matrix (-out_mat )\n";
    die;
  }
if ( !$seq)
  {
    print " You must provide sequences (-seq )\n";
    die;
  }

open (F, $seq);
while ( <F>)
  {
    /(.*)/;
    @list=(@list, $1);
  }
close (F);


$tmp_file="Pavie_tmp_file$$.tmp";
$max=(($#list+1)*($#list+1)-($#list+1))/2;
for ($a=0; $a<=$#list-1; $a++)
 {
   for ($b=$a+1; $b<=$#list; $b++, $c++)
     {
       $completion=($c*100)/$max;
       print "\rAlign $a with $b [$completion % Complete]";
       
       `align_nw4pavie.pl $list[$a] $list[$b] $in_mat >>$tmp_file`;
     }
 }
open (F, $tmp_file);
while (<F>){print "$_";}



`pavie_aln2log_matrix.pl $tmp_file >$out_mat `;
#unlink $tmp_file;




