#!/usr/bin/env perl
use File::Copy;
use Env qw(HOST);
use Env qw(HOME);
use Env qw(USER);
$x_field=0;
$y_field=1;
$interval=0;
$file="stdin";
$print_avg=1;
$print_sd=0;
$print_sum=0;
$print_n=0;
foreach $value ( @ARGV)
    {
	if ($value ne $ARGV[$np]) 
	    {
	    ;
	    }
	elsif($value eq "-print_all")
	    {
	    $print_sd=$print_avg=$print_n=$print_sum=1;
	    $np++;
	    }
	elsif($value eq "-print_sum")
	    {
	    $print_sum=1;
	    $print_avg=0;
	    $np++;
	    }
	elsif($value eq "-print_n")
	    {
	    $print_n=1;
	    $print_avg=0;
	    $np++;
	    }
	elsif($value eq "-print_avg")
	    {
	    $print_avg=1;
	    $print_avg=0;
	    $np++;
	    }
	elsif($value eq "-print_sd")
	    {
	    $print_sd=1;
	    $print_avg=0;
	    $np++;
	    }
	elsif($value eq "-header")
	    {
	    $header=1;
	    $np++;
	    }
	elsif ($value eq "-interval")
	    {
	    $interval= $ARGV[++$np];
	    $np++;
    	    }
	elsif ($value eq "-x_field")
	    {
	    $x_field= $ARGV[++$np]-1;
	    $np++;
    	    }
	elsif ($value eq "-y_field")
	    {
	      
	    while ($ARGV[$np+1] && !($ARGV[$np+1]=~/\-/))
	      {
		$y_field[$nyf++]=$ARGV[++$np]-1;
		$y_field_set=1;
	      }

	    $np++;
    	    }
	elsif ($value eq "-file")
	    {
	    $file= $ARGV[++$np];
	    $file_set=1;
	    $np++;
    	    }       
	elsif ( $value eq "-spec")
	  {
	    $do_specificity=1;
	  }
	elsif ( $value eq "h" ||  $value eq "-h" || $value eq "-H" || $value eq "-help" || $value eq "help")
	  {
	    print STDOUT "data_analyse: Analyse and discretization of data\n";
	    print STDOUT "       -file:    <file containing the data to analyze>,.<def=STDIN>\n";
	    print STDOUT "       -x_field: <field containing the X>,...............<Def=0>\n";
	    print STDOUT "       -y_field: <field containing the Y>,...............<Def=1>\n";
	    print STDOUT "       -interval:<Interval size on the Y>,...............<Def=0>\n";
	    print STDOUT "       -interval:<0:only one interval>\n";
	    print STDOUT "       -header  : print column header \n";
	    print STDOUT "       -spec    : read TP TN FP FN \n";
	    exit (0);
	  }
	elsif ($value=~/-/)
	  {
	    print "$value is not a valid FLAG[FATAL]\n";
	    exit (0);
	   } 
	elsif ($list eq "") 
	    {
	    $file=$ARGV[$np];
	    $np++;
	    }
	
	
      }





if ($file eq "stdin")
	{
	$remove_file=1;
	$file="tmp$$";
	open (F, ">$file");
	while (<STDIN>)
		{
		print F $_;
		}
	close (F);
	 
	;}


open(F,$file);

while (<F>)
          {
	    $line=$_;
	    @list=/(\S+)/g;
	
	    if ( $do_specificity)
	      {
		$tp+=$list[0];
		$tn+=$list[1];
		$fp+=$list[2];
		$fn+=$list[3];
	      }
	    if ($interval==0){$bin=0;}
	    else{$bin=$list[$x_field]/$interval;}
	
	    for ( $a=0; $a<$nyf; $a++)
	      {
		$sum{$a}{$bin}+=$list[$y_field[$a]];
		$sum2{$a}{$bin}+=$list[$y_field[$a]]*$list[$y_field[$a]];
		$n{$a}{$bin}++;
	      }
	    if ( $bin>$max_bin){$max_bin=$bin;}
	  }


if ( !$do_specificity)
  {
    
    if ( $header)
      {
	for ( $a=0; $a<=$max_bin; $a++)
	  {
	    $i=$interval*$a;
	    printf "%10.2f ", $i;
	  }
	 printf "\n";
      }
    for ( $b=0; $b<$nyf; $b++)
      {
	for ( $a=0; $a<=$max_bin; $a++)
      	  {
	    $i=$interval*$a;
	    if ( $n{$b}{$a}==0)
	      {
	      $avg=0;
	      $sd=0;
	      }
	    else
	      {
		$avg=$sum{$b}{$a}/$n{$b}{$a};
		$sd=sqrt($sum2{$b}{$a}*$n{$b}{$a}-$sum{$b}{$a}*$sum{$b}{$a})/($n{$b}{$a}*$n{$b}{$a});
		
		if ($print_n) {printf "%10.2f ", $n{$b}{$a};}
		if ($print_sum){printf "%10.2f ", $sum{$b}{$a};}
		if ($print_avg){printf "%10.2f ", $avg}
		if ($print_sd) {printf "%10.2f ", $sd;}
	      }
	  }
	printf ("\n");
      }
  }
else
  {
    $positive=$tp+$fn;
    $negative=$tn+$fp;

    $sn=$tp/($tp+$fn);
    $sp=$tn/($tn+$fp);
    $sp2=$tp/($tp+$fp);
    $smc=($tp+$fn)/($tp+$fn+$fp+$tn);
    printf  ("|Sn %10.2f |Sp %10.2f |Sp2=%10.2f |Pos %10d |Neg %10d\n", $sn, $sp, $sp2, $positive, $negative);
  }


if ( $remove_file){unlink $file;}

