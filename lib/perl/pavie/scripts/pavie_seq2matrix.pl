#!/usr/bin/env perl
#Takes a list of file where each file contains one sequence in fasta format
#outputs a matrix
@al=@ARGV;
$delta_min=1;
for ( $a=0; $a<=$#al; $a++)
  {
    if ( $al[$a] eq "-seq")
      {
	$seq=$al[++$a];
      }

    elsif ($al[$a] eq "-out_mat")
      {
	$out_mat=$al[++$a];
      }
    elsif ($al[$a] eq "-delta_min")
      {
	$delta_min=$al[++$a];
      }
  }
if ( !$out_mat)
  {
    print " You must give a name for the output matrix (-out_mat )\n";
    die;
  }
$n=0;
$diff=$delta_min+1;
$new_matrix=$out_mat.$n;
while ( $diff>$delta_min)
  {
    #1
    if ($old_matrix){$in_matrix="-in_mat $old_matrix";}
    else {$in_matrix="";}
    `compute_pavie_matrix.pl -seq $seq $in_matrix -out_mat $new_matrix`;
    if ($old_matrix)
      {

	$diff=`compare_pavie_matrix.pl $old_matrix $new_matrix`;
	print "\nNew_matrix: $new_matrix Diff $diff\n";
      }
    else
      {
	print "\nFirst Matrix: $new_matrix\n";
      }
    $n++;
    $old_matrix=$new_matrix;
    $new_matrix=$out_mat.$n;
  }
