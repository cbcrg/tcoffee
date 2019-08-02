#!/usr/bin/env perl
use File::Copy;
use Env qw(HOST);
use Env qw(HOME);
use Env qw(USER);
$x_field=0;
$y_field=1;
$y_field_set=1;
$nyf=1;

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
	elsif($value eq "-s")
	     {
	       $step=$ARGV[++$np];
	       $np++;
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
	    $nyf=0;  
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
	    print STDOUT "       -s:<Step on the  X, 0 means non sliding bins>\n";
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



if ($interval && $step)
  {
    my $nl;
    open(F,$file);
    while (<F>)
      {
	$line=$_;
	
	if (!/\S/){next;}
	@list=($line=~/(\S+)/g);
	$val{$nl}{x}=$list[$x_field];
	$val{$nl}{y}=$list[$y_field[0]];
	$nl++
      }
    close (F);
    
    for (my $a=$min; $a<($max+$interval); $a+=$step)
      {
	my ($avgx, $avgy, $cn);
	
	my $rmin=$a-$interval;
	my $rmax=$a;
	$cn=0;
	for (my $b=0; $b<$nl; $b++)
	  {
	    my $x=$val{$b}{x};
	    my $y=$val{$b}{y};
	    if ($x<=$rmax && $x>=$rmin)
	      {
		$avgx+=$x;
		$avgy+=$y;
		$cn++;
		$tcn++;
		$val{$b}{used}=1;
	      }
	  }
	if ($cn)
	  {
	    $avgx/=$cn;
	    $avgy/=$cn;
	  }
	printf "%.3f %.3f %.3f\n", $avgx, $avgy, $avgx-$avgy;
      }
    for (my $a=0; $a<$nl; $a++)
      {
	if ( !$val{$a}{used})
	  {
	    print "---$val{$a}{x}; $val{$a}{y}\n";
	  }
      }
  }
else
  {
    if ($interval && $max)
      {
	$interval_size=($max-$min)/$interval;
      }
    elsif ($interval)
      {
	open(F,$file);  
	my $set_max=0;
	my $set_min=0;
	while (<F>)
	  {
	    my $v=$_;
	    chomp($v);
	    print "--$v--";
	    
	    if ($v<$min ||!$set_min){$set_min=1;$min=$v;}
	    if ($v>$max ||!$set_max){$set_max=1;$max=$v;}
	  }
	close (F);
	print "$min $max uuuu";
	$interval_size=($max-$min)/$interval;
      }
    open(F,$file);  
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
	printf ( "%4d %4d ", $interval_size*$a, $interval_size*($a+1));
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
	    if ($print_n) {printf "%15.4f ", $n{$b}{$a};}
	    if ($print_sum){printf "%15.4f ", $sum{$b}{$a};}
	    if ($print_avg){printf "%15.4f ", $avg}
	    if ($print_sd) {printf "%15.4f ", $sd;}
	  }
	printf ("\n");
      }
  }

if ( $remove_file){unlink $file;}

