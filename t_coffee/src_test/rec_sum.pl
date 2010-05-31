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
	elsif($value eq "-sd")
	    {
	    $print_sd=1;
	    $print_avg=0;
	    $np++;
	    }
	elsif($value eq "-h")
	    {
	    $header=1;
	    $np++;
	    }
	elsif ($value eq "-i")
	    {
	    $interval= $ARGV[++$np];
	    $np++;
    	    }
	elsif ($value eq "-r")
	    {
	    $min= $ARGV[++$np];
	    $max= $ARGV[++$np];
	    $np++;
    	    }
	
	elsif ($value eq "-x")
	    {
	    $x_field= $ARGV[++$np]-1;
	    $np++;
    	    }
	elsif ($value eq "-y")
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
	elsif ( $value eq "h" ||  $value eq "-h" || $value eq "-H" || $value eq "-help" || $value eq "help")
	  {
	    print STDOUT "data_analyse: Analyse and discretization of data\n";
	    print STDOUT "       -file:    <file containing the data to analyze>,.<def=STDIN>\n";
	    print STDOUT "       -x: <field containing the X>,...............<Def=0>\n";
	    print STDOUT "       -y: <field containing the Y>,...............<Def=1>\n";
	    print STDOUT "       -i:<Interval size on the X>,...............<Def=0>\n";
	    print STDOUT "       -i:<0:only one interval>\n";
	    print STDOUT "       -r:<Range of the X>\n";
	    print STDOUT "       -sd: print standard deviation on the Y";
	    print STDOUT "       -h  : print column header \n";
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

if ($interval)
  {
    $interval_size=($max-$min)/$interval;
  }
while (<F>)
  {
    $line=$_;
    if (!/\S/){next;}
    @list=($line=~/(\S+)/g);
    
    if ($interval==0){$bin=0;}
    else{$bin=int (($list[$x_field]-$min)/($interval_size));}

    
    if ($bin && $bin==$interval){$bin--;}
    for ( $a=0; $a<$nyf; $a++)
      {
	$sum{$a}{$bin}+=$list[$y_field[$a]];
	$sum2{$a}{$bin}+=$list[$y_field[$a]]*$list[$y_field[$a]];
	$n{$a}{$bin}++;
      }
  }

if (!$interval){$interval=1;}
for ( $a=0; $a<$interval; $a++)
  {
    printf ( "%3d %3d ", $interval_size*$a, $interval_size*($a+1));
    for ( $b=0; $b<$nyf; $b++)	
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
	  }
	if ($print_n) {printf "%10.4f ", $n{$b}{$a};}
	if ($print_sum){printf "%10.4f ", $sum{$b}{$a};}
	if ($print_avg){printf "%10.4f ", $avg}
	if ($print_sd) {printf "%10.4f ", $sd;}
      }
    printf ("\n");
  }


if ( $remove_file){unlink $file;}

