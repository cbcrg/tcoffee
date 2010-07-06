#!/usr/bin/env perl

use HTTP::Date;
use strict;
use FileHandle;
srand();
my $index;
my $SHIFT=4;
my $LOG_ZERO=-999999999;
my $LOG_UNDERFLOW=0.0000000000000000000001;
my $HEADER;
my $d={};
my $M={};
my $A={};
my $file=shift (@ARGV);
$d=parse_data ($file);
foreach my $c (split (/\-+/,join(" ", @ARGV)))
  {
    $A=string2hash ($c,$A);
    if ($c=~/filter/)
      {
	($d,$A)=filter_data ($d,$A);
      }
    elsif ($c=~/seq2model/)
      {
	($d,$A)=seq2model ($d,$A);
      }
    elsif ($c=~/filter/)
      {

	($d,$A)=seq2data ($d,$A);
      } 
    elsif ($c=~/filter/)
      {
	($d,$A)=seq2data ($d,$A);
      } 
    elsif ($c=~/decode/)
      {
	($d,$A)=decode ($d,$A);
      }
    elsif ($c=~/outdata/)
      {
	display_data ($A->{modelR}, $A->{outdata});
      }
    elsif ($c=~/outmodel/)
      {
	display_model ($A->{dataR}, $A->{outmodel});
      }
    elsif ($c=~/test/)
      {
	test_bw_trainning();
      }
  }

die;
#test_bw_trainning ();

my $d={};
$d=&parse_data ($ARGV[0], $d,0,0.01);

$d=filter_data ($d, "Nature", "char", "keep", "food");
#%data=filter_data (\%data, "Cage", "int", "keep", 1);
$d=data2index ($d);
$d=data2bin ($d, "Value", 10, 0.02);
my ($P, $M)=multi_baum_welch ($d,5,10,10,10);
display_decode ($d);
display_model ($M, "MODEL");die;




my @int=data2intervals ($d,"full");
for ( $a=1; $a<=12; $a++){data2string ($d, $a);}


for (my $a=0; $a<$#int; $a++)
  {
    my $t=interval2count($d,$int[$a], $int[$a+1]);
    print "$int[$a] $int[$a+1] => $t\n";
  }
data2log_odd($d,"all","sc_", $int[$a], $int[$a+1]);
data2log_odd($d,"all","cd_", $int[$a], $int[$a+1]);

foreach my $diet (("sc_","cd_"))
  {
    for (my $c=1; $c<=12; $c++)
      {
	for (my $a=0; $a<$#int; $a++)
	  {
	    
	    data2log_odd($d, $c,$diet, $int[$a], $int[$a+1]);
	  }
      }
  }
print "\n\n\n";

sub display_data
  {
    my $d=shift;
    my $file=shift;
    my $F= new FileHandle;
    
    if (!$file){$F=0;}
    else {open ($F, ">$file");}
    
    print $F "$HEADER";
    foreach my $c (sort ({$a<=>$b}keys(%$d)))
      {
	foreach my $i (sort {$a<=>$b}keys (%{$d->{$c}}))
	  {
	    if ($F){print $F "#d;CAGE;$c;";}
	    else {print "#d;CAGE;$c;";}
	    foreach my $k (sort (keys (%{$d->{$c}{$i}})))
	      {
		
		if ($F){print $F "$k;${$d->{$c}{$i}}";}
		else {print $F "$k;${$d->{$c}{$i}}";}
	      }
	    if ($F){print "\n";}
	    else {print $F "\n";}
	  }
      }
    if ($file){close ($F);}
  }

sub parse_data
  {
    my $file=shift;
    my $data={};
    
    my $F=new FileHandle;
    my $lineN;
    open($F, $file);
    
    while (<$F>)
      {
	my $line=$_;
	my $LineN++;
	
	if ( $line=~/#d/)
	  {
	    my $L={};
	    chomp $line;
	    my @v=split (/;/,$line);
	    shift @v; #get rid of the line header
	    while (@v)
	      {
		my $key=shift @v;
		my $value= shift @v;
		$L->{$key}=$value;
		
	      }
	    $L->{LineN}=$LineN;
	    
	    if ($L->{Type})
	      {
		my $c=$L->{CAGE};
		my $ch=$L->{Channel};
		my $t=$L->{StartT};
	
		foreach my $k (%$L)
		  {
		    $data->{$c}{$t}{$k}=$L->{$k};
		  }
	      }
	  }
	else {$HEADER-=$line;}
      }
    #reformat/annotate fields fields
    $data=&data2overlap($data);
    $data=&channel2Nature($data);
    $data=&channel2correct_channel ($data);
    #$data=&filter_overlap ($data,1);
    return $data;
  }

  #  #Filter data
#    $data=&filter_data ($data,"Value","float","rm",-9999999,$T);
#    $data=&filter_overlap ($data, $dup);
    
    
    
    
#    return $data;
#  }

sub channel2correct_channel
    {
      #THis function corrects all sorts of labelling errors made by the acquisition equipment
      my $d=shift;
      my ($tot, $n);

      foreach my $c (sort(keys (%$d)))
	{
	  foreach my $t (sort(keys (%{$d->{$c}})))
	    {
	      my $Name=$d->{$c}{$t}{Name};
	      my $Channel=$d->{$c}{$t}{Channel};
	      my $Caption=$d->{$c}{$t}{Caption};
	      
	      $Channel=~/.*(\d)/;
	      my $i=$1;
	      
	      #This is meant to correct the wrong labelling of the intakes: 1 and 2 are always frink, 3 and 4 are always food
	      if ($i==1){$d->{$c}{$t}{Caption}="Drink 1";}
	      elsif ($i==2){$d->{$c}{$t}{Caption}="Drink 2";}
	      elsif ($i==3){$d->{$c}{$t}{Caption}="Food 1";}
	      elsif ($i==4){$d->{$c}{$t}{Caption}="Food 2";}
	      else
		{
		  print "\n*** ERROR: unknown index for the Intake\n";
		}
	      if ($Caption ne $d->{$c}{$t}{Caption}){$tot++;}
	      $n++;
	    }
	}
      print "\nclean_data: relabled $tot values out of $n\n";
      return $d;
    }

sub channel2Nature
      {
	# This function creates a label describing the precise content of each intake
	my $d=shift;
		
	foreach my $c (sort(keys (%$d)))
	  {
	    foreach my $t (sort(keys (%{$d->{$c}})))
	      {
		my $Name=$d->{$c}{$t}{Name};
		my $Channel=$d->{$c}{$t}{Channel};
		my $Caption=$d->{$c}{$t}{Caption};
		
		$Channel=~/.*(\d)/;
		my $i=$1;
		my $Nature="";
		
		
		$d->{$c}{$t}{SlotI}=$i;
		
		if ($Caption=~/Food/){$Nature="food";}
		else {$Nature="drink";}
		
		$Name=lc ($Name);
		$Name=~s/\s//g;
		
		
		
		if ($Nature eq "food")
		  {
		    # print "--$Name--\n";
		    if ($Name eq "sc"){$Nature.="_sc";}
		    elsif ($Name =~/cd/ && $Nature eq "food")
		      {
			if    ($i==1 && (($Name =~/slota/) ||($Name =~/ina/) )){$Nature.="_cd";}
			elsif ($i==2 && (($Name =~/slotb/) ||($Name =~/inb/) )){$Nature.="_cd";}
			elsif ($i==3 && (($Name =~/slotc/) ||($Name =~/inc/) )){$Nature.="_cd";}
			elsif ($i==4 && (($Name =~/slotd/) ||($Name =~/ind/) )){$Nature.="_cd";}
			else  {$Nature.="_sc";}
		    }
		    else
		      {
			print "ERROR: $Name\n";
		      }
		  }
		if    ($Name =~/sc/){$Nature="sc_".$Nature;}
		elsif ($Name =~/cd/){$Nature="cd_".$Nature;}
		
		$d->{$c}{$t}{Nature}=$Nature;
		
	      }
	  }
	return $d;
      }
sub filter_data
   {
     my $d=shift;
     my $A=shift;
     
     my $field=$A->{field};
     my $min=$A->{min};
     my $max=$A->{max};
     
     
     my $tot=0;
     my $n=0;
     
     #foreach my $c ( keys (%$A)){print "$c : $A->{$c}\n";}die;
     
     foreach my $c (sort(keys (%$d)))
       {
	 foreach my $t (sort(keys (%{$d->{$c}})))
	   {
	     $n++;
	     my $mark=0;
	     
	     if (defined ($A->{matches}) && $c=~$A->{matches}){$mark=1;}
	     elsif (defined ($A->{equals}) && $c eq $A->{equals}){$mark=1;}
	     elsif ( defined ($A->{min})   && $d->{$c}{$t}{$field}>=$min && $d->{$c}{$t}{$field}<=$max ){$mark=1;}
	     
	     if    ($mark==1 && $action eq "rm"){delete($d->{$c}{$t}); $tot++;}
	     elsif ($mark==0 && $action eq "keep"){delete($d->{$c}{$t}); $tot++;}
	   }
       }
     print "\nFiltering: Removed $tot values out of $n (Min=$min Max=$max)\n";
     return $d;
    }
sub filter_overlap
    {
      
      my $d=shift;
      my $T=shift;
      my $tot=0;
      my $n=0;
      if ($T eq "no"){return $d;}
      foreach my $c (sort(keys (%$d)))
	{
	  foreach my $t (sort(keys (%{$d->{$c}})))
	    {
	      $n++;
	      if ( $d->{$c}{$t}{Collision} && $d->{$c}{$t}{Collision}>$T){delete ($d->{$c}{$t});$tot++;}
	    }
	}
      print "\nCollisions: Removed $tot values out of $n (T: $T)\n";
      return $d;
    }
    
sub display_data 
    {
      my $d=shift;
      
      foreach my $c (sort(keys (%$d)))
	{
	  foreach my $t (sort(keys (%{$d->{$c}})))
	    {
	      
	      print "$c - $t - $d->{$c}{$t}{Channel} -  $d->{$c}{$t}{Value}\n";
	    }
	}
    }
    
sub parse_header
  {
    my $file=shift;
    my $data=shift;
    my $F=new FileHandle;
    
    open($F, $file);
    
    while (<$F>) 
      {
	my $line=$_;
	if ( $line=~/#h/)
	  {
	    chomp $line;
	    my @k=split ($line, /;/);
	    
	    $data->{$k[1]}{$k[2]}{$k[3]}{$k[4]}=$k[5];
	  }
	elsif ($line=~/#d/)
	  {
	    last;
	  }
      }
    close ($F);
    return $data;
  }
sub data2overlap
    {
      my $d=shift;
      my $print=shift;
      my  $nc;
      my ($n,$tot);
      
      foreach my $c (sort(keys (%$d)))
	{
	  my $pStartT=-1;
	  my $pEndT=-1;
	  my $pChannel=-1;
	  my $pStartL;
	  my $pEndL;
	  my $pFile;
	  my $pValue;
	  my $pc;
	  my $pt;
	  
	  if ($print)
	    {
	      print "\nCHECK CAGE $c\n";
	    }
	  foreach my $t (sort(keys (%{$d->{$c}})))
	    {
	      $n++;
	      my $StartT=$d->{$c}{$t}{StartT};
	      my $EndT=$d->{$c}{$t}{EndT};
	      my $Channel=$d->{$c}{$t}{Channel};
	      my $StartL=$d->{$c}{$t}{StartL};
	      my $EndL=$d->{$c}{$t}{EndL};
	      my $File=$d->{$c}{$t}{File};
	      my $Value=$d->{$c}{$t}{Value};
	      if ($pStartT!=-1)
		{
		  if (($StartT<$pEndT))
		    {
		      
		      my $delta=$pEndT-$StartT;
		      my $v1=$delta/($EndT-$StartT);
		      my $v2=$delta/($pEndT-$pStartT);
		      my $v=($v1<$v2)?$v2:$v1;
		      
		      $d->{$c}{$t}{Collision}=$v1;
		      $d->{$pc}{$pt}{Collision}=$v2;
		      
		      
		      
		      if ($print)
			{
			  print "***** ERROR: OVERLAP: CAGE $c --- $delta :\n";
			  print "\t\tC: $pChannel [$pStartT -- $pEndT] [$pStartL -- $pEndL] VALUE: $pValue File: $pFile\n";
			  print "\t\tC: $Channel [$StartT -- $EndT] [$StartL -- $EndL] File: VALUE: $Value$File\n";
			  }
		      $tot++;
		    }
		}
	      $pc=$c;
	      $pt=$t;
	      $pStartT=$StartT;
	      $pEndT=$EndT;
	      $pStartL=$StartL;
	      $pEndL=$EndL;
	      $pFile=$File;
	      $pChannel=$Channel;
	      $pValue=$Value;
	      
	    }
	}
      print "\nOverlap: $tot values out of $n\n";
      return $d;
    }
sub data2string 
    {
      my $d=shift;
      my $Cage=shift;
      my  $nc;
      
      foreach my $c (sort(keys (%$d)))
	{
	  if ($Cage!=$c){next;}
	  foreach my $t (sort(keys (%{$d->{$c}})))
	    {
	      printf "$c:   %12s => %5.2f\n",$d->{$c}{$t}{Nature}, $d->{$c}{$t}{Value};
	    }
	}
    }
sub data2log_odd 
    {
      my $d=shift;
      my $Cage=shift;
      my $Nature=shift;
      my $start=shift;
      my $end=shift;
      
      my  $nc;
      my $M={};
      my $chc={};
      
            
      foreach my $c (sort(keys (%$d)))
	{
	  if (!($Cage eq "all" || $Cage==$c)){next;}
	  my ($ch,$pEndT);
	  foreach my $t (sort(keys (%{$d->{$c}})))
	    {
	      if ($start && ($t<$start || $t>$end)){;}
	      else
		{
		  my $cch=$d->{$c}{$t}{Nature};
		  if ($cch=~/$Nature/)
		    {
		      if ($ch)
			{
			  $M->{$ch}{$cch}{Count}{tot}++;
			  $M->{$ch}{$cch}{Interval}{tot}+=$d->{$c}{$t}{StartT}-$pEndT;
			}
		      if (1==2 && $pEndT>$d->{$c}{$t}{StartT})
			{
			  print "\n****** ERROR: $c => Start= $d->{$c}{$t}{StartT} pEnd: $pEndT\n";
			}

		      $chc->{$cch}{Count}{tot}++;
		      $chc->{$cch}{Value}{tot}+=$d->{$c}{$t}{Value};
		      $chc->{$cch}{Duration}{tot}+=$d->{$c}{$t}{Duration};
		      if ($pEndT){$chc->{$cch}{Interval}{tot}+=$d->{$c}{$t}{StartT}-$pEndT;}
		      $ch=$cch;
		      $pEndT=$d->{$c}{$t}{EndT};
		      $nc++;
		    }
		}
	    }
	}
      if (!$nc){return;}
	  
      for my $c(keys(%$chc))
	{
	  for my $x (keys (%{$chc->{$c}}))
	   {
	     $chc->{$c}{$x}{avg}=$chc->{$c}{$x}{tot}/$nc;
	   }
       }

      print "CAGE: $Cage DIET: $Nature START: $start END: $end Nitervals: $nc\n";
      foreach my $c(sort(keys(%$chc)))
	    {
	      printf "\t%12s Freq: %6.2f Dur: %6d Value: %6.2f Int: %6d N: %6d\n", $c,$chc->{$c}{Count}{avg},$chc->{$c}{Duration}{avg}, $chc->{$c}{Value}{avg},$chc->{$c}{Interval}{avg},$chc->{$c}{Count}{tot} ;
	    }
      
      foreach my $c1 (sort(keys (%$chc)))
	{
	  foreach my $c2 (sort(keys (%$chc)))
	    {
	      my $avgInt=$M->{$c1}{$c2}{Interval}{tot}/$M->{$c1}{$c2}{Count}{tot};
	      my $v1=$M->{$c1}{$c2}{Count}{tot}/($nc-1);
	      my $v2=$chc->{$c1}{Count}{avg}*$chc->{$c2}{Count}{avg};
	      my $logv=($v2>0 && $v1>0)?mylog($v1/$v2):0;
	      printf "\t\t%12s ---> %12s  %6.2f Interval: %6d Count: %6d\n", $c1, $c2 ,$logv, $avgInt,$M->{$c1}{$c2}{Count}{tot};
	    }
	}
    }


sub interval2count
    {
      my $d=shift;
      my $s=shift;
      my $e=shift;
      my ($tot);
      
      foreach my $c (sort(keys (%$d)))
	{
	  foreach my $t (sort(keys (%{$d->{$c}})))
	    {
	      if ($t>$s && $t<$e){$tot++;}
	    }
	}
      return $tot;
    }

sub data2intervals
    {
      my $d=shift;
      my $n=shift;
      
      my ($a, $start, $end, @list, $delta);
      
      $start=$end=-1;
      
      foreach my $c (sort(keys (%$d)))
	{
	  foreach my $t (sort(keys (%{$d->{$c}})))
	    {

	      my $StartT=$d->{$c}{$t}{StartT};
	      my $EndT=$d->{$c}{$t}{EndT};
	      
	      if ($start==-1 || $start>$StartT){$start=$StartT;}
	      if ($end==-1    || $end<$EndT){$end=$EndT;}
	    }
	}
      
      if ($n eq "hour"){$delta=3600;}
      elsif ($n eq "day"){$delta=3600*24;}
      elsif ($n eq "week"){$delta=3600*24*7;}
      elsif ($n eq "month"){$delta=3600*24*31;}
      elsif ($n eq "year"){$delta=3600*24*365;}
      elsif ($n eq "full"){return ($start, $end);}
      else {$delta=($end-$start)/$n;}
      
      for ($a=$start; $a<$end; $a+=$delta)
	{
	  push (@list, $a);
	}
      push (@list, $a);
      

      
      return @list;
    }

sub data2index
  {
    my $d=shift;
    my $nd={};
    
    #indexes data so that it can be dealt with by viterbi trainning

    foreach my $c (keys(%$d))
      {
	my $i;
	foreach my $t (sort (keys(%{$d->{$c}})))
	  {
	    ++$i;
	    foreach my $k(keys(%{$d->{$c}{$t}}))
	      {
		$nd->{$c}{$i}{$k}=$d->{$c}{$t}{$k};
	      }
	    delete ($d->{$c}{$t});
	  }
      }
    return $nd;
  }
	    
sub data2bin 
  {
    my $S=shift;
    my $field=shift;
    my $nbin=shift;
    my $delta=shift;

    

    my $count={};
    
    if (!$delta)
      {
	my($min,$max)=datafield2minmax($S,$field);
	my $delta=($max-$min)/$nbin;
	print "$max - $min\n";
      }
    else
      {
	$delta=0.02;
      }
    
    foreach my $c (keys(%$S))
      {
	foreach my $i (keys (%{$S->{$c}}))
	  {
	    my $v=$S->{$c}{$i}{$field};
	    if ( $v<0){;}
	    else {
	      my $bin=int($S->{$c}{$i}{$field}/$delta)+1;
	      if ($bin>$nbin){$bin=$nbin;}
	      elsif ($bin<-$nbin){$bin=$nbin*-1;}
	      $S->{$c}{$i}{bin}=$bin;
	    }
	  }
      }
    return $S;
  }

sub datafield2minmax
  {
    my $d=shift;
    my $field=shift;
    my $set;
    my ($max, $min);
    
    print "Field=$field";
    
    foreach my $c (keys(%$d))
      {
	
	foreach my $i (keys %{$d->{$c}})
	  {
	    my $v=$d->{$c}{$i}{$field};
	    print "$c -- $i -- $field -- $v\n";
	    if (!$set)
	      {
		$max=$min=$v;
		$set=1;
	      }
	    else
	      {
		$max=($max<$v)?$v:$max;
		$min=($min>$v)?$v:$min;
	      }
	  }
      }
    return ($min,$max);
  }
	
	

				  
#####################################################3
#emits the model for the occasionally dishonnest casino
sub get_file 
  {
    my $A=shift;
    my $t=shift;
    my $r={};
    my $test;
    my $ref=$t."R";
    
    if (!$A->{$t} && ! $A->{$ref} )
      {
	die;
      }
    elsif ($A->{$ref})
      {
	my $b=$A->{$ref};
	$r=$b;
      }
    else
      {
	if ($t =~/model/){$r=parse_model ($A->{$t});}
	else {$r=parse_data ($A->{$t});}
	$A->{$ref}=$r;
      }
    return $r;
  }
    
sub seq2model 
  {
    my $A=shift;
    my $M={};
    
    my $d=get_file($A, "data");
    
    if (!$A->{ntries}){$A->{ntries}=5;}
    if (!$A->{nit}){$A->{nit}=1000;}
    if (!$A->{nbin}){$A->{nbin}=10;}
    if (!$A->{delta}){$A->{delta}=0.02;}
    if (!$A->{field}){$A->{field}="Value";}
    
    if (!$A->{nstate}){$A->{nstates}=2;}
    if (!$A->{nemit}){$A->{nemit}=$A->{nbin};}
    
    $d=data2bin ($d, $A->{"Value"}, $A->{nbin},$A->{delta});
    my @topo;
    for (my $a=0; $a<$A->{nstate}; $a++){push (@topo, $A->{nemit});}
				       
    my ($P, $M)=multi_baum_welch($d,$A->{ntries}, $A->{nit}, @topo);
    $A->{modelP}=$P;
    $A->{modelR}=$M;
  }
sub decode
  {
    my $A=shift;
    
    my $M=get_file($A, "model");
    my $d=get_file($A, "data");
    
    ($d,$A->{modelP})=viterbiL ($M, $d);
    
    $d=posteriorL($M, $d);
  }
  
  
sub test_bw_trainning()
  {
    my ($RP,$P);
    my $M={};
    my $score;
    my $RM=ODHC2model();
    my $RM2;
    
    my $nseeds=1;
    my $nseq=1;
    my $len=200;
    
    my $S=model2sequence ($RM,1,1000);
   
    
    $RM=model2modelL($RM);
    $RM2=sequence2model ($S);
    my $RP=seq2probaL($RM, $S);
    
    ($P,$M)=multi_baum_welch ($S,5,100, 6, 6);
    
    
    
    ($S,$score)=viterbiL ($M, $S);
    $S=posteriorL($M, $S);
    
    display_decode ($S);
    $M  =modelL2model($M);
    $RM =modelL2model($RM);
    $RM2=modelL2model($RM2);
    
    
    display_decode ($S);
    
    display_model ($RM);
    display_model ($RM2);
    
    display_model ($M);
    die;
  }

sub ODHC2model
  {
    my $M={};

    $M->{'ST::START'}{'ST::fair'  }=0.5;
    $M->{'ST::START'}{'ST::unfair'}=0.5;

    $M->{'ST::fair'} {'ST::fair'  }=0.95;
    $M->{'ST::fair'} {'ST::unfair'}=0.05;
    
    $M->{'ST::unfair'} {'ST::fair'  }=0.1;
    $M->{'ST::unfair'} {'ST::unfair'}=0.9;
    
    $M->{'ST::fair'} {'ST::END'  }=0.5;
    $M->{'ST::unfair'} {'ST::END'}=0.5;
    
    $M->{'ST::fair'} {'1'}=1/6;
    $M->{'ST::fair'} {'2'}=1/6;
    $M->{'ST::fair'} {'3'}=1/6;
    $M->{'ST::fair'} {'4'}=1/6;
    $M->{'ST::fair'} {'5'}=1/6;
    $M->{'ST::fair'} {'6'}=1/6;
    
    $M->{'ST::unfair'} {'1'}=1/10;
    $M->{'ST::unfair'} {'2'}=1/10;
    $M->{'ST::unfair'} {'3'}=1/10;
    $M->{'ST::unfair'} {'4'}=1/10;
    $M->{'ST::unfair'} {'5'}=1/10;
    $M->{'ST::unfair'} {'6'}=5/10;
    
    return $M;
  }

sub model2sequence
  {
    my $M=shift;
    my $n=shift; #number of strings
    my $l=shift; #string length
    my $S={};
        
    my $state="ST::START";
    my $symbol;
    for (my $j=0; $j<$n; $j++)
      {
	for (my $i=1; $i<=$l; $i++)
	  {
	    ($state, $symbol)=model2emit($M, $state);
	    
	    $S->{$j}{$i}{bin}=$symbol;
	    $S->{$j}{$i}{RST}=$state;
	    
	  }
      }
    return $S;
  }
  
sub model2emit
  {
    my $M=shift;
    my $start=shift;
    my ($state, $bin);
    my $r_state=rand(1);
    my $r_bin=rand(1);
    my $p=0;
    
    foreach my $k (keys(%$M))
      {
	$p+=$M->{$start}{$k};
	if ( $r_state<=$p){$state=$k;last;}
      }
    if ($state eq "ST::END" || $state eq "ST::START"){return model2emit ($M, $start);}
    $p=0;
    foreach my $bin (keys (%{$M->{$state}}))
      {
	if ( $bin =~/ST::/){;}
	else
	  {
	    $p+=$M->{$state}{$bin};
	    if ($r_bin<=$p)
	      {
		return ($state, $bin);
	      }
	  }
      }
    
    return ($state, $bin);
  }
sub modelL2model
  {
     my $M=shift;
     my $tag=shift;
         
     
     foreach my $k(keys(%$M))
       {
	 foreach my $l (keys(%$M))
	   {
	     $M->{$k}{$l}=(!$M->{$k}{$l} ||$M->{$k}{$l}==$LOG_ZERO)?0:exp($M->{$k}{$l});
	   }
	}
     
      foreach my $k(keys(%$M))
	{
	  foreach my $l (keys(%{$M->{$k}}))
	    {
	      if (!($l=~/ST::/))
		{
		  $M->{$k}{$l}=(!$M->{$k}{$l} ||$M->{$k}{$l}==$LOG_ZERO)?0:exp($M->{$k}{$l});
		}
	    }
	}
      return $M;
    }
sub model2modelL
  {
     my $M=shift;
     my $tag=shift;
     
     foreach my $k(keys(%$M))
       {
	 foreach my $l (keys(%$M))
	   {
	     $M->{$k}{$l}=(!$M->{$k}{$l} ||$M->{$k}{$l}<0.00000001)?$LOG_ZERO:mylog($M->{$k}{$l});
	   }
	}
      
      foreach my $k(keys(%$M))
	{
	  foreach my $l (keys(%{$M->{$k}}))
	    {
	      
	      if (!($l=~/ST::/))
		{
		  $M->{$k}{$l}=(!$M->{$k}{$l} ||$M->{$k}{$l}<0.00000001)?$LOG_ZERO:mylog($M->{$k}{$l});
		}
	    }
	}
     return $M;
    }
sub display_hash
    {
      my $dp=shift;
      my $h=shift;
      my @kl=@_;
      
      
      foreach my $k (keys(%$h))
	{
	  print "**** $k\n";
	  if ($dp==1)
	    {
	      foreach my $l (@kl)
		{
		  print "-- $l";
		}
	      print "-- $k -- $h->{$k}\n";
	    }
	  else
	    {
	      @kl=push (@kl, $k);
	      display_hash ($dp-1, $h, @kl);
	    }
	}
    }
sub display_decode 
    {
      my $S=shift;
      my $tag=shift;
      my $T={};
      
      foreach my $j (sort {$a<=>$b}keys (%$S))
	{
	  my $L=keys(%{$S->{$j}});
	  for (my $i=1; $i<=$L; $i++)
	    {
	      my $bin   =$S->{$j}{$i}{bin};
	      my $state =$S->{$j}{$i}{viterbi};
	      my $rstate=$S->{$j}{$i}{RST};
	      $T->{$rstate}{$state}++;
	      print "$tag :: $j :: $i :: $bin == VITERBI: $state";
	      if ($rstate){print "== REAL: $rstate";}
	      if ( exists ($S->{$j}{$i}{bpost}))
		   {
		     printf "== POSTERIOR: $S->{$j}{$i}{bpost}{k}";
		     
		   }
	      print "\n";
	    }
	}
      
      foreach my $k (keys(%$T))
	{
	  foreach my $l (keys (%{$T->{$k}}))
	    {
	      print "TOT: $k -- $l => $T->{$k}{$l}\n";
	    }
	}
    }
#################################################
sub viterbi_trainningL
  {
    my $M=shift;
    my $S=shift;
    
    my $A={};
    my $E={};
    my $ns;
    my $score;
    
    ($S,$score)=viterbiL($M,$S);
    
    foreach my $j (keys(%$S))
      {
	my $L=keys(%{$S->{$j}});
	for (my $i=1; $i<=$L; $i++)
	  {

	    my $l=$S->{$j}{$i}{viterbi};
	    if ($i>1)
	      {
		my $k=$S->{$j}{$i-1}{viterbi};
		$A->{$k}{$l}++;
	      }
	    $E->{$l}{$S->{$j}{$i}{bin}}++;
	  }
      }
    
    foreach my $k (keys (%$M))
      {
	foreach my $l (keys(%{$M->{$k}}))
	  {
	    if (($l=~/ST::/)){$A->{$k}{$l}+=1;}
	    else {$E->{$k}{$l}+=1;}
	  }
      }
    #update A/Model
    foreach my $k (keys (%$M))
      {
	my $num;
	foreach my $lp (keys(%$M)){$num+=$A->{$k}{$lp};}
	foreach my $l (keys(%$M)){$M->{$k}{$l}=$A->{$k}{$l}/$num;}
      }
    
    # update E/model 
    foreach my $k (keys (%$M))
      {
	my $num;
	foreach my $lp (keys(%{$M->{$k}}))
	  {
	    if (!($lp =~/ST::/)){$num+=$E->{$k}{$lp};}
	  }
 	foreach my $l (keys(%{$M->{$k}}))
	  {
	    if (!($l =~/ST::/)){$M->{$k}{$l}=$E->{$k}{$l}/$num;}
	  }
      }
    
    return model2modelL($M);
  }


sub multi_baum_welch
  {
    my $S=shift;
    my $multi=shift;
    my $nit=shift;
   

    my @topology=@_;
    
    my $best_score;
    my $BM;
    my $maxidle=5;
    my $score;
    
    for (my $i=0;$i<$multi; $i++)
      {
	my ($idle, $PP);
	print "---- $i -----\n";
	my $P;
	my $M=topology2modelL ("rand",@topology);
	my $cont=1;
	for (my $it=0; $it<$nit && $idle<$maxidle && $cont; $it++)
	  {
	    ($P,$M)=baum_welchL ($S,$M);
	    
	    if ($PP)
	      {
		my $delta=$P-$PP;
		if ($delta<-1){$cont=0;}
		elsif ($delta<0.001){$idle++;}
		else {$idle=0;}
	      }
	    int($P);
	    print "\t$it ==> $P [$best_score --] [$idle]\n";
	    $PP=$P;
	  }
	my $score=seq2probaL($M, $S);
	if (!$BM || $score>$best_score)
	  {
	    $BM=$M;
	    $best_score=$score;
	  }
      }
    
    
    ($S, $score)=viterbiL($BM, $S);
    $S=posteriorL($BM, $S);
    
    return ($best_score, $BM);
  }
sub posteriorL
    {
      my $M=shift;
      my $S=shift;
   
   foreach my $j (keys (%$S))
      {
	
	my $L=keys (%{$S->{$j}});
	my $F={};
	my $B={};

	my ($P,$B)=backwardL($M,$S->{$j});#log_space
	my ($P,$F)=forwardL ($M, $S->{$j});
	
	for (my $i=1; $i<=$L; $i++)
	  {
	    my $symbol=$S->{$j}{$i}{'bin'};
	    my $bpost_score;
	    my $bpost_k;
	    foreach my $k (keys (%$M))
	      {
		if (!(exists ($M->{$k}{$symbol}))){next;}
		my $p=log_divide (log_multiply($F->{$i}{$k},$B->{$i}{$k}),$P);
		$S->{$j}{$i}{post}{$k}=$p;
		if (!$bpost_score || $p>$bpost_score)
		  {
		    $bpost_score=$p;
		    $bpost_k=$k;
		  }
	      }
	    $S->{$j}{$i}{bpost}{k}=$bpost_k;
	    $S->{$j}{$i}{bpost}{score}=$bpost_score;
	  }
      }
   return $S;
 }

sub baum_welchL
  {

    my $S=shift;
    my $M=shift;
    my $A={};
    my $E={};
    my $P;

    
    foreach my $j (keys (%$S))
      {
	
	my $L=keys (%{$S->{$j}});
	my $F={};
	my $B={};
	
	my ($P,$B)=backwardL($M, $S->{$j});#log_space
	my ($P,$F)=forwardL ($M, $S->{$j});
	
	
	
	#Update A
	foreach my $k (keys (%$M))
	  {
	    foreach my $l(keys(%$M))
	      {
		$A->{$j}{$k}{$l}=$LOG_ZERO;
		if (!$A->{$k}{$l}){$A->{$k}{$l}=$LOG_ZERO;}
		
		for (my $i=1; $i<$L; $i++)
		  {
		    my $symbol=$S->{$j}{$i+1}{'bin'};
		    if (!(exists ($M->{$l}{$symbol}))){next;}
		    
		    my $fo=$F->{$i}{$k};#log space
		    my $ba=$B->{$i+1}{$l};#log_space
		    my $tr=$M->{$k}{$l}; #log_space
		    my $em=$M->{$l}{$symbol};;#log_space
		    $A->{$j}{$k}{$l}=log_add ($A->{$j}{$k}{$l},log_multiply($fo,$tr,$em,$ba));
		  }
		$A->{$j}{$k}{$l}=log_divide($A->{$j}{$k}{$l},$P);
		$A->{$k}{$l}=log_add ($A->{$k}{$l},$A->{$j}{$k}{$l});
	      }
	  }
	
	#update Emissions
	foreach my $k (keys (%$M))
	  {
	    foreach my $b (keys (%{$M->{$k}}))
	      {
		if ($b=~/ST::/){next;}
		$E->{$j}{$k}{$b}=$LOG_ZERO;
		if ( ! exists ($E->{$k}{$b})){$E->{$k}{$b}=$LOG_ZERO;}
		
		for (my $i=1; $i<=$L; $i++)
		  {
		    if ($S->{$j}{$i}{bin} eq $b)
		      {
			my $p=$E->{$j}{$k}{$b};
			my $q=log_multiply ($F->{$i}{$k},$B->{$i}{$k});
			$E->{$j}{$k}{$b}=log_add($p,$q);
		      }
		  }
		
		$E->{$j}{$k}{$b}=log_divide($E->{$j}{$k}{$b},$P);
		$E->{$k}{$b}=log_add($E->{$k}{$b},$E->{$j}{$k}{$b});
	      }
	  }
      }
    
    
    foreach my $k (keys(%$M))
      {
	foreach my $l (keys (%{$M->{$k}}))
	  {
	    if (($l=~/ST::/ )){$A->{$k}{$l}=exp($A->{$k}{$l});}
	    else {$E->{$k}{$l}=exp($E->{$k}{$l});}
	  }
      }
    
     #add pseudo-counts
    foreach my $k (keys (%$M))
      {
	foreach my $l (keys(%{$M->{$k}}))
	  {
	    if (($l=~/ST::/)){$A->{$k}{$l}+=1;}
	    else {$E->{$k}{$l}+=1;}
	  }
      }
    #update A/Model
    foreach my $k (keys (%$M))
      {
	my $num;
	foreach my $lp (keys(%$M)){$num+=$A->{$k}{$lp};}
	foreach my $l (keys(%$M)){$M->{$k}{$l}=$A->{$k}{$l}/$num;}
      }
    
    # update E/model 
    foreach my $k (keys (%$M))
      {
	my $num;
	foreach my $lp (keys(%{$M->{$k}}))
	  {
	    if (!($lp =~/ST::/)){$num+=$E->{$k}{$lp};}
	  }
 	foreach my $l (keys(%{$M->{$k}}))
	  {
	    if (!($l =~/ST::/)){$M->{$k}{$l}=$E->{$k}{$l}/$num;}
	  }
      }
    $M=model2modelL($M);
    $P=seq2probaL($M,$S);
    return ($P,$M);
  }
    
sub seq2probaL
  {
    my $M=shift;
    my $S=shift;
    
    my $TP;

    for my $j (keys(%$S))
      {
	my ($P, $f)=forwardL($M, $S->{$j});
	$TP+=$P;
      }
    return $TP;
  }
sub forwardL 
    {
      my $M=shift;
      my $S=shift;
      my $f={};
      my $P;
      my $L=keys(%$S);

      foreach my $k (keys(%$M)){$f->{0}{$k}=$LOG_ZERO;}
      $f->{0}{'ST::START'}=0;
     
      for (my $i=1; $i<=$L; $i++)
	{
	  foreach my $l (keys(%$M))
	    {

	      $f->{$i}{$l}=$LOG_ZERO;
	      my $emit=(!exists($M->{$l}{$S->{$i}{bin}}))?$LOG_ZERO:$M->{$l}{$S->{$i}{bin}};
	      
	      foreach my $k (keys(%$M))
		{
		  $f->{$i}{$l}=log_add($f->{$i}{$l}, log_multiply ($f->{$i-1}{$k},$M->{$k}{$l}));
		  my $v1=myexp($f->{$i}{$l});
		  my $v2=myexp($f->{$i-1}{$k});
		  my $v3=myexp($M->{$k}{$l});
		  #print "\t----L: V1: $v1 V2: $v2 V3: $v3\n";
		}
	      my $v1=myexp($f->{$i}{$l});
	      my $v2=myexp($emit);
	      
	      $f->{$i}{$l}=log_multiply($f->{$i}{$l},$emit);
	      my $v3=myexp($f->{$i}{$l});
	      #print "----L: $i V1: $v1 V2: $v2 V3: $v3\n";
	      
	    }
	}

      
      $P=$LOG_ZERO;
      foreach my $k (keys (%$M))
	{
	  $P=log_add ($P, log_multiply ($f->{$L}{$k},$M->{$k}{'ST::END'}));
	}
      
      return ($P,$f);
    }

sub backwardL 
    {
      my $M=shift;
      my $S=shift; 
      
      my $B={};
     
      my $P;
      my $L=keys (%$S);
      

      foreach my $k (keys(%$M)){$B->{$L}{$k}=$M->{$k}{'ST::END'};}

      for (my $i=$L-1; $i>=1; $i--)
	{
	  foreach my $k (keys(%$M))
	    {
	      $B->{$i}{$k}=$LOG_ZERO;
	      
	      foreach my $l (keys(%$M))
		{
		  if (!exists($M->{$l}{$S->{$i+1}{bin}})){next;}
		  my $x=$M->{$k}{$l};
		  my $y=$M->{$l}{$S->{$i+1}{bin}};
		  my $z=$B->{$i+1}{$l};
		  my $p=$B->{$i}{$k};
		  my $q=$x+$y+$z;
		  
		  $B->{$i}{$k}=log_add ( $B->{$i}{$k}, log_multiply($x,$y,$z));
		  
		}
	    }
	}
      return (0,$B);
    }
sub viterbiL 
  {
    my $M=shift;
    my $S= shift;
    my ($max_k, $ptr_k);
    
    foreach my $j (keys(%$S))
      {
	my $L=keys(%{$S->{$j}});
	my $PTR={};
	my $V={};
	my ($path, $ppath);
	
	foreach my $k (keys(%$M)){$V->{0}{$k}=$LOG_ZERO;}
	$V->{0}{'ST::START'}=0;
	
	for (my $i=1; $i<=$L; $i++)
	  {
	    my $symbol=$S->{$j}{$i}{bin};
	    
	    foreach my $l (keys (%$M))
	      {
		if ( !exists ($M->{$l}{$symbol}))
		  {
		    $V->{$i}{$l}=$LOG_ZERO;
		    $PTR->{$i}{$l}=$LOG_ZERO;
		  }
		else
		  {
		    $max_k=$LOG_ZERO;
		    $ptr_k="";
		    foreach my $k (keys (%$M))
		      {
			my $v=log_multiply($V->{$i-1}{$k},$M->{$k}{$l});
			my $v1=myexp($V->{$i-1}{$k});
			my $v2=myexp($M->{$k}{$l});
			my $v3=myexp($v);
			#print "\t$k $l---> V1: $v1 * V2: $v2 = V3: $v3\n";
			if ($v>$max_k || $max_k==$LOG_ZERO)
			  {
			    $max_k=$v;
			    $ptr_k=$k;
			  }
		      }
		    $V->{$i}{$l}=log_multiply($M->{$l}{$S->{$j}{$i}{bin}},$max_k);
		    $PTR->{$i}{$l}=$ptr_k;
		    my $v=exp($max_k);
		  }
	      }
	  }
	
	$max_k=$LOG_ZERO;
	$ptr_k="";
	foreach my $k (keys (%$M))
	  {

	    my $vv=log_multiply($V->{$L}{$k},$M->{$k}{'ST::END'});
	    if ($vv>$max_k  || $max_k==$LOG_ZERO)
	      {
		$max_k=$vv;
		$ptr_k=$k;
	      }
	  }
	for (my $i=$L; $i>=1; $i--)
	  {
	    $S->{$j}{$i}{viterbi}=$ptr_k;
	    $ptr_k=$PTR->{$i}{$ptr_k};
	  }
	
      }
    return ($S,$max_k);
  }

sub topology2modelL
  {
    my ($model,@topo)=@_;
    my $M=topology2model ($model,@topo);
    return model2modelL($M);
  }
sub topology2model
  {
    my ($model,@states)=@_;
    my $SL={};
    my $M={};
   
    for (my $a=1; $a<=($#states+1); $a++)
      {
	$SL->{"ST::$a"}=$states[$a-1];
      }
    $SL->{"ST::START"}=0;
    $SL->{"ST::END"}=0;
    
    foreach my $st1 (keys (%$SL))
      {
	my $tot;

	#set emmissions
	for (my $a=1; $a<=$SL->{$st1}; $a++)
	  {
	    $tot+=$M->{$st1}{$a}=($model eq "rand")?rand (1000):100;
	  }
	for (my $a=1; $a<=$SL->{$st1}; $a++)
	  {
	    $M->{$st1}{$a}/=$tot;
	  }


	#Set Transitions
	$tot=0;
	foreach my $st2 (keys (%$SL))
	  {
	    $tot+=$M->{$st1}{$st2}=($model eq "rand")?rand (1000):100;
	  }
	foreach my $st2 (keys (%$SL))
	  {
	    $M->{$st1}{$st2}/=$tot;
	  }
      }
    display_model ($M);
    return $M;
  }
sub display_path 
    {
      my $M=shift;
      my $T=shift;
      my $tag=shift;
      
      my $L=keys (%$T);
      print "L=$L $tag\n";
      for ( my $i=1; $i<=$L; $i++)
	{
	  my $v;
	  foreach my $k (keys (%$M))
	    {
	      if ($tag eq "exp"){$v=exp($T->{$i}{$k});}
	      else {$v=$T->{$i}{$k};}
	      print "-- $i $k ==> $v ($tag)\n";
	    }
	}
    }

sub display_model
  {
    my $M=shift;
    my $fname=shift;
    my $F= new FileHandle;
    
    if ($fname){open ($F,">$fname");}
    else {$F=0;}
    print "#### MODEL: $fname\n";
    print "#### STATES\n";
    foreach my $k(keys(%$M))
      {
	foreach my $l (keys(%$M))
	  {
	    if ($fname){printf "$k;$l;%7.5f\n",$M->{$k}{$l};}
	    else {printf "$k;$l;%7.5f\n",$M->{$k}{$l};}
	    
	  }
      }
    print "#### EMISSIONS\n";
    foreach my $k(keys(%$M))
      {
	foreach my $l (keys(%{$M->{$k}}))
	  {
	    if (!($l=~/ST::/))
	      {
		
		if ($fname){printf $F "$k;$l;%7.5f\n",$M->{$k}{$l};}
		else {printf "$k;$l;%7.5f\n",$M->{$k}{$l};}
	      }
	    
	  }
      }
    print "#### END\n";
    if ($fname) 
      {
	print STDERR "---- Dumped Model in file $fname\n";
	close ($F);
      }
    
    return;
  }
sub sequence2model
      {
	my $S=shift;
	my $A={};
	my $E={};
	my $M={};
	
	foreach my $j (keys(%$S))
	  {
	    my $L=keys (%{$S->{$j}});
	    for (my $i=2; $i<=$L; $i++)
	      {
		my $s=$S->{$j}{$i}{bin};
		my $cstate=$S->{$j}{$i}{viterbi};
		my $pstate=$S->{$j}{$i-1}{viterbi};
		$A->{$pstate}{$cstate}++;
		$E->{$cstate}{$s}++;
	      }
	  }
	foreach my $k (keys (%$A))
	  {
	    my $tot;

	    $tot=0;
	    foreach my $l (keys(%{$A->{$k}}))
	      {
		$tot+=$A->{$k}{$l};
	      }
	    foreach my $l  (keys(%{$A->{$k}}))
	      {
		$M->{$k}{$l}=$A->{$k}{$l}/$tot;
	      }
	    
	    $tot=0;
	    foreach my $l (keys(%{$E->{$k}}))
	      {
		$tot+=$E->{$k}{$l};
	      }
	    foreach my $l (keys(%{$E->{$k}}))
	      {
		$M->{$k}{$l}=$E->{$k}{$l}/$tot;
	      }
	  }


	foreach my $k (keys(%$M))
	  {
	    $M->{$k}{"ST::END"}=1;
	    $M->{"ST::START"}{$k}=1;
	  }
	
	return model2modelL($M);
      }
	
sub log_divide 
  {
    my $x=shift;
    my $y=shift;
    
    if ($x==$LOG_ZERO || $y==$LOG_ZERO || $y==0){return $LOG_ZERO;}
    return $x-$y;
  }
sub log_multiply
  {
    my @l=@_;
    my $r;
    
    foreach my $v (@l)
      {
	if ($v==$LOG_ZERO){return $LOG_ZERO;}
	$r+=$v;
      }
    return $r;
  }
sub log_add 
  {
    my ($x,$y)=@_;
    
    
    if ($x==$LOG_ZERO){return $y;}
    elsif ($y==$LOG_ZERO){return $x;}
    elsif ($x>=$y)
      {
	return $x+log(1+exp($y-$x));
	$x=(($x==$LOG_ZERO) || ($y-$x)>=$LOG_UNDERFLOW)?$y:mylog($y-$x)+$x;
      }
    else{return log_add ($y,$x);}
  }

sub mylog
  {
   my $x=shift;
   if ( $x<$LOG_UNDERFLOW){return $LOG_ZERO;}
   else
     {
       return log($x);
       
      }
   if ($x <= 1.00) 
     {return ((-0.009350833524763 * $x + 0.130659527668286) * $x + 0.498799810682272) * $x + 0.693203116424741;}
   if ($x <= 2.50)
     {return ((-0.014532321752540 * $x + 0.139942324101744) * $x + 0.495635523139337) * $x + 0.692140569840976;}
   if ($x <= 4.50) 
     {return ((-0.004605031767994 * $x + 0.063427417320019) * $x + 0.695956496475118) * $x + 0.514272634594009;}
  
  return ((-0.000458661602210 * $x + 0.009695946122598) * $x + 0.930734667215156) * $x + 0.168037164329057;
} 
sub expLookup
  {
    my $x=shift;
    
   if ($x > -2)
     {
       if ($x > -0.5)
	 {
	   if ($x > 0){return exp($x);}
	   return (((0.03254409303190190000*$x + 0.16280432765779600000)*$x + 0.49929760485974900000)*$x + 0.99995149601363700000)*$x + 0.99999925508501600000;
	 }
       if ($x > -1){return (((0.01973899026052090000*$x + 0.13822379685007000000)*$x + 0.48056651562365000000)*$x + 0.99326940370383500000)*$x + 0.99906756856399500000;}
       return (((0.00940528203591384000*$x + 0.09414963667859410000)*$x + 0.40825793595877300000)*$x + 0.93933625499130400000)*$x + 0.98369508190545300000;
     }
   if ($x > -8)
     {
       if ($x > -4){return (((0.00217245711583303000*$x + 0.03484829428350620000)*$x + 0.22118199801337800000)*$x + 0.67049462206469500000)*$x + 0.83556950223398500000;}
       return (((0.00012398771025456900*$x + 0.00349155785951272000)*$x + 0.03727721426017900000)*$x + 0.17974997741536900000)*$x + 0.33249299994217400000;
     }
   if ($x > -16){return (((0.00000051741713416603*$x + 0.00002721456879608080)*$x + 0.00053418601865636800)*$x + 0.00464101989351936000)*$x + 0.01507447981459420000;}
   return 0;
 }
sub myexp
    {
      my $x=shift;
      if ( $x==$LOG_ZERO){return 0;}
      return exp($x);
    }
sub myexp2
  {
    return shift;
  }
sub array2hash

  {
    my $arrayR=shift;
    my $A=shift;
    my ($v, $k);
    my @array=@$arrayR;
    
    while (($k=shift(@array)))
      {
	my $v=shift (@array);
	$k=~s/-//g;
	$A->{$k}=$v;
      }
    return $A;
  }

sub string2hash 
  {
    my $s=shift;
    my $h=shift;

    my @l=split (/\s+/, $s);
    shift @l;
    return array2hash (\@l, $h);
  }
