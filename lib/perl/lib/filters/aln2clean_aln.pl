#!/usr/bin/env perl 

$upper=9;
$lower=0;
 foreach ($np=0; $np<=$n_para; $np++)
    {
    $value=$ARGV[$np];

    if ($value eq "-infile")
      {
      $infile= $ARGV[++$np];
      }
    elsif ($value eq "-lib")
      {
        $lib= $ARGV[++$np];
      }
    elsif ($value eq "-upper_threshold")
      {
        $upper_threshold= $ARGV[++$np];
      }
     elsif ($value eq "-lower_threshold")
      {
        $lower_threshold= $ARGV[++$np];
      }
    elsif ($value eq "-scorefile")
      {
	$scorefile= $ARGV[++$np];
      }
    else
      {
	print "$value is not a valid argument";
	exit (1);
      }

    $infile=~/([^.]*)\.[^.]*/;
    $suffix=$1;
    
    if ( ! -e $lib  && ! -e $scorefile)
      {
	if (!$lib){$lib="$suffix.tc_lib";}
	
	`t_coffee S$infile -out_lib=$lib -convert -outfile=no -quiet stdout`;
       }

    if (! -e $scorefile)
      {
	if (!$scorefile){$scorefile="$suffix.score_ascii";}
	`t_coffee -infile $infile  -in L$lib -output score_ascii -score -outfile= $scorefile`;
      }
    
    $conversion="[$upper_threshold lower_threshold]";
    `seq_reformat -in $infile -struc_in $scorefile -struc_in_f number_aln -action +convert $conversion '#x -output fasta_aln'`;
   
      
   
    
    
	 
	 
