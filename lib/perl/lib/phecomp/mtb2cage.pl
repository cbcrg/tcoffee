#!/usr/bin/env perl

use HTTP::Date;
use strict;
use FileHandle;
my $index;
my $SHIFT=4;


&mtb2intervals ("Intake 1;Intake 2;Intake 3;Intake 4", @ARGV);
	      
#process_mbt (@ARGV);

sub process_mbt
  {
    my $infile=shift;
    
    my ($infile)=@ARGV;
    my %data;


    %data=&mtb2parse ($infile);
    %data=&merge_fields (\%data, "Food", "Intake 3", "Intake 4");
    %data=&intake2intervals (\%data, "Food");
    
    
    
    print "START: $data{$infile}{'HEADER'}{'EHEADER'}{'StartTime'} END: $data{'HEADER'}{'EndTime'}";
    die;
  }
  
sub mtb2intervals
    {
      my $channels=shift;
      my @files=@_;
      my %data;
      my @ch;
      my $ncages;


      @ch=split (/;/,$channels);
      foreach my $f (@files)
	{
	  my (@ci,%pv,%pe,%pt, $count, $count2);
	  my $F = new FileHandle;
	  
	  %data=&mtb2header($f, \%data);
	  my @chL=keys (%{$data{$f}{"[Intake Channels]"}});
	  my $stime=$data{$f}{"HEADER"}{'EHEADER'}{"StartStamp"};
	  $ncages=$data{$f}{'HEADER'}{'EHEADER'}{'Ncages'};
	  &display_header  (\%data);
	  $F=set_handle($f,"DATASECTION");

	  while (<$F>)
	    {
	      my $line=$_;
	      
	      my %dataline=mtb2parse_line ($f,$line,\%data,\@chL);
	      if ($count==1000)
		{
		  print STDERR "------>$f: $count2\n";
		  $count=0;
		}
	      $count++; $count2++;
	      
	      for (my $c=1; $c<=$ncages; $c++)
		{
		  foreach my $ch (@ch)
		    {
		      
		      my $v=$dataline{$c}{$ch}{Value};
		      my $t=$dataline{$c}{$ch}{Time};
		      my $l=$t-$stime;
		      my $e=$dataline{$c}{$ch}{Event};
		      my $deltaE =$e-$pe{$c}{$ch};
		      my $deltaV=$v-$pv{$c}{$ch};
		      my $deltaT=$t-$pt{$c}{$ch};
		      my $ci;

		      $ci=$data{'INTERVALS'}{$c}{$ch}{0};
			  
		      $data{'INTERVALS'}{$c}{$ch}{$ci}{EndT}=$t;
		      $data{'INTERVALS'}{$c}{$ch}{$ci}{Value}=$deltaV;
		     # print "[** $e $deltaE]";
		      
		      if (!$data{'INTERVALS'}{$c}{$ch}{0})
			{
			  $data{'INTERVALS'}{$c}{$ch}{0}=0;
			  $ci=++$data{'INTERVALS'}{$c}{$ch}{0};
			  $data{'INTERVALS'}{$c}{$ch}{$ci}{StartT}=$t;
			  $data{'INTERVALS'}{$c}{$ch}{$ci}{Type}=$e;
			  $data{'INTERVALS'}{$c}{$ch}{$ci}{Index}=$ci;
	
			  $pv{$c}{$ch}=$v;
			  $pt{$c}{$ch}=$t;
			}
		      elsif ($deltaE)
			{
			  
			  
			  $ci=++$data{'INTERVALS'}{$c}{$ch}{0};
			  
			  $data{'INTERVALS'}{$c}{$ch}{$ci}{StartT}=$t;
			  $data{'INTERVALS'}{$c}{$ch}{$ci}{Type}=$e;
			  $data{'INTERVALS'}{$c}{$ch}{$ci}{Index}=$ci;
			  $pv{$c}{$ch}=$v;
			  $pt{$c}{$ch}=$t;
			}
		      $pe{$c}{$ch}=$e;
		    }
		}
	    }
	}
      
      for (my $c=1; $c<=$ncages; $c++)
	{
	  foreach my $ch (@ch)
	    {
	      
	      my $ci=$data{'INTERVALS'}{$c}{$ch}{0};
	      for (my $i=1; $i<=$ci; $i++)
		{
		  my $start=$data{'INTERVALS'}{$c}{$ch}{$i}{StartT};
		  my $end=$data{'INTERVALS'}{$c}{$ch}{$i}{EndT};
		  my $v=$data{'INTERVALS'}{$c}{$ch}{$i}{Value};
		  my $T=$data{'INTERVALS'}{$c}{$ch}{$i}{Type};
		  my $duration=$end -$start;
		  
		  printf "CAGE ; %2d ; CH; %s ; INDEX; %5d ; START ; %6d ; END ; %6d ; DURATION ; %5d ; VALUE; %5.2f ; TYPE ; %d\n",$c, $ch, $i, $start, $end, $duration, $v, $T;
		}
	    }
	}
      die;
      
      return %data;
    }

sub mtb2parse_line 
    {
      my $file=shift;
      my $line=shift;
      my $dataR=shift;
      my $channelsR=shift;
     
      
      my %data=%$dataR;
      my @channels=@$channelsR;
      my @list;
      my ($Ncages, $stime,$time, $id, %dataline);
      
      $id="[Intake Channels]";
      $stime=$data{$file}{"HEADER"}{'EHEADER'}{"StartStamp"};
      
      $Ncages=$data{$file}{"HEADER"}{'EHEADER'}{"Ncages"};
      

     
      
      if ($line)
	{
	  my ( $i, $cage, $t, $lineN);
	  $line =~ s/[=]/ /g;
	  $line =~ s/,/\./g;
	  
	  @list=split ( /\s+/, $line);
	  $time=$lineN=shift (@list);
	  $time=$lineN+$stime;
	  for (my $a=0; $a<$SHIFT; $a++){shift(@list);}
	  
	  foreach my $ch (@channels)
	    {
	      
	      my $chi=$data{$file}{$id}{$ch}{Index}-1;
	      for (my $cage=1; $cage<=$Ncages;$cage++)
		{
		  my $event;
		  my $i=$chi*$Ncages+($cage-1);
		  my $v=$list[$i];
		  if ($v=~/\*/){$event=1;$v=~s/\*//;}
		  else {$event=0;}
		  
		  $dataline{$cage}{$ch}{'Value'}=$list[$i];
		  $dataline{$cage}{$ch}{'Event'}=$event;
		  my $t=$dataline{$cage}{$ch}{'Time'}=$time;
		}
	    }
	 
	}
      
      
      return %dataline;
    }
sub mtb2header 
    {
    my $file=shift;
    my $dataR=shift;
    my $F = new FileHandle;
    my %header=%$dataR;
    
    my ($id,@id_list,$index, $v, $time);
    my ($Nintake, $Ncages);
     
    
   
    $index=&set_channel (\%header,$file,"[Intake Channels]",$index,"activity",("Type=output", "Captions=activity", "unit=undef"));
   
    $index=&set_channel (\%header,$file,"[Intake Channels]",$index,"rearing" ,("Type=output", "Captions=rearing",  "unit=undef"));
   
    open ($F, "$file");
    while (<$F>)
      {
	my $line=$_;
	if ($line=~/(\[.*\])/)
	  {
	    $id=$1;
	    push (@id_list, $id);
	  }
	else
	  {
	    if ($id=~/ANIMALS DATA/)
	      {
		if (($line=~/(\D+)(\d+)=([\w .\/:]*)/))
		  {
		    $header{$file}{$id}{$2}{$1}=$3;
		  }
	      }
	    elsif ($id=~/Intake/)
	      {
		my ($channel, $name, $val);
		if (($line=~/(Intake\s+\d+)\s+(.+)=([\w .\/:]*)/))
		  {
		    $channel=$1;
		    $name=$2;
		    $val=$3;
		    
	
		    if (!$header{$file}{$id}{$channel})
		      {
			$index++;
			$header{$file}{$id}{$channel}{Type} ="output";
			$header{$file}{$id}{$channel}{Index}= $index;
			
		      }
		    $header{$file}{$id}{$channel}{$name}=$val;
		  }
	      }
	    elsif ($id=~/DATASECTION/){last;}
	    else 
	      {
		if (($line=~/(.*)=([\w .\/:]*)/))
		  {
		    $header{$file}{$id}{$id}{$1}=$2;
		  }
	      }
	  }
      }
    close ($F);
    $header{$file}{'HEADER'}{'EHEADER'}{'Ncages'}=$Ncages=&header2value("Number of cages", \%header, $file);
    $time=$header{$file}{"HEADER"}{'EHEADER'}{'StartStamp'}=str2time(&header2value("Date and time", \%header, $file));
   
    return %header;
  }
    

sub intake2intervals (\%data, "Food")
{
  my $file=shift;
  my $dataR= shift;
  my $ch =shift;
  my %data =%$dataR;
  my $start_t=$data{$file}{'HEADER'}{'EHEADER'}{'StartTime'};
  my $end_t  =$data{$file}{'HEADER'}{'EHEADER'}{'EndTime'};
  my $ncages =$data{$file}{'HEADER'}{'EHEADER'}{'Ncages'};
  
  for (my $c=1; $c<=$ncages; $c++)
    {
      my ($ci,$pv);
      for (my $t=$start_t; $t<$end_t; $t++)
	{
	  my $v=$data{$file}{"DATASECTION"}{$c}{$ch}{$t};
	  my $lineN=$t-$data{$file}{"HEADER"}{'EHEADER'}{"StartStamp"};
	  
	  if ($v!=$pv || !$ci)
	    {
	      $ci++;
	      my $delta=$pv-$v;
	      
	      if ($delta>=0.02)
		{
		  printf "\n\t---- Big   Delta: Line: %d Cage: %d Delta: %.2f [$pv :: $v]\n",$lineN , $c, $delta;
		}
	      elsif ($delta>0)
		{
		 printf "\n\t---- Small Delta: Line: %d Cage: %d Delta: %.2f [$pv :: $v]\n",$lineN , $c, $delta;
		}
	      
	      $data{'INTERVALS'}{$c}{$ch}{$ci}{'start'}=$t;
	      $data{'INTERVALS'}{$c}{$ch}{$ci}{'delta'}=$delta;
	      $data{'INTERVALS'}{$c}{$ch}{0}++;
	      }
	  $data{'INTERVALS'}{$c}{$ch}{$ci}{'end'}=$t;
	  $pv=$v;
	}
      print STDERR "CAGE: $c N_intervals: $data{'INTERVALS'}{$c}{$ch}{0}\n";
      die;
    }

  for (my $a=1; $a<=$data{'INTERVALS'}{1}{$ch}{0}; $a++)
    {
      my $d= $data{'INTERVALS'}{1}{$ch}{$a}{'end'}-$data{'INTERVALS'}{1}{$ch}{$a}{'start'};
      printf "I:%5d %5d =>%.2f\n", $a, $d, $data{'INTERVALS'}{1}{$ch}{$a}{'delta'};
    }
  return %data;
}

 
sub merge_fields
  {
    #merge the channels in CH list and creat a new channel

    my $dataR=shift;
    my $new_ch=shift;
    my @ch_list=@_;
    my ($c, $t,$ncages, $start_t, $end_t, $ch, $time, $x);
    my %data;

    %data=%$dataR;
    
   
    $start_t=$data{'HEADER'}{'EHEADER'}{'StartTime'};
    $end_t=$data{'HEADER'}{'EHEADER'}{'EndTime'};
    $ncages=$data{'HEADER'}{'EHEADER'}{'Ncages'};
    print STDERR "MERRGE CHANNELS";

   
    
    
    for ($c=1; $c<=$ncages; $c++)
      {
	
	for ($t=$start_t; $t<$end_t; $t++)
	  {
	    my $v;
	    foreach $ch (@ch_list)
	      {

		$v+=$data{"DATASECTION"}{$c}{$ch}{$t};
	      }
	    $x=$data{"DATASECTION"}{$c}{$new_ch}{$t}=$v;
	  }
	print "$c: $x\n";
      }
    return %data;
  }
    

sub mtb2parse
  {
    my $count;
    my $infile= shift;
    my $F = new FileHandle;
    my %header;
    
    my ($id,@id_list,@output, $index, $v, $time);
    my (%header,$Nintake, $Ncages);
    
    open ($F, "$infile");
   
    %header=&set_channel (\%header,$infile,"[Intake Channels]","activity",("Type=output", "Captions=activity", "unit=undef"));
    push (@output, "activity");
    %header=&set_channel (\%header,$infile,"[Intake Channels]","rearing" ,("Type=output", "Captions=rearing",  "unit=undef"));
    push (@output, "rearing");
    
    
    while (<$F>)
      {
	my $line=$_;
	if ($line=~/(\[.*\])/)
	  {
	    $id=$1;
	    
	    push (@id_list, $id);
	    print "$id---";
	    
	  }
	else
	  {
	    if ($id=~/ANIMALS DATA/)
	      {
		if (($line=~/(\D+)(\d+)=([\w .\/:]*)/))
		  {
		    $header{$infile}{$id}{$2}{$1}=$3;
		  }
	      }
	    elsif ($id=~/Intake/)
	      {
		my ($channel, $name, $val);
		if (($line=~/(Intake\s+\d+)\s+(.+)=([\w .\/:]*)/))
		  {
		    $channel=$1;
		    $name=$2;
		    $val=$3;
		    
		    if (!$header{$infile}{$id}{$channel})
		      {
			$index++;
			$header{$infile}{$id}{$channel}{Type} ="output";
			$header{$infile}{$id}{$channel}{Index}= $index;
			push (@output, "$channel");
		      }
		    $header{$infile}{$id}{$channel}{$name}=$val;
		  }
	      }
	    elsif ($id=~/DATASECTION/){last;}
	    else 
	      {
		if (($line=~/(.*)=([\w .\/:]*)/))
		  {
		    $header{$infile}{$id}{$id}{$1}=$2;
		  }
	      }
	  }
	if ($id=~/DATASECTION/){last;}
      }
    
    &display_header(\%header);
    
    $header{$infile}{'HEADER'}{'EHEADER'}{'Ncages'}=$Ncages=&header2value("Number of cages", \%header);
    $header{$infile}{"HEADER"}{'EHEADER'}{"StartStamp"}=str2time(&header2value("Date and time", \%header));
        
    while (<$F>)
      {
	my $line=$_;
	my (@list, $time, $i, $cage, $t, $lineN);
	$line =~ s/[*=]/ /g;
	$line =~ s/,/\./g;
	
	
	
	@list=split ( /\s+/, $line);
	$lineN=$t=$time=shift (@list);
	
	$time+=$header{$infile}{"HEADER"}{"StartStamp"};
	
	if (!$header{$infile}{"HEADER"}{'EHEADER'}{"StartTime"})
	  {$header{$infile}{"HEADER"}{'EHEADER'}{"StartTime"}=$time;}
	
	$header{$infile}{"HEADER"}{"EndTime"}=$time;
	
	for (my $a=0; $a<$SHIFT; $a++){shift(@list);}
	
	if ( $count == 10000){print STDERR "\n\t$time\t";$count=0;}
	for (my $o=0; $o<=$#output;$o++)
	  {
	    for ($cage=1; $cage<=$Ncages; $cage++)
	      {
		$header{$infile}{'DATASECTION'}{$cage}{$output[$o]}{$time}=shift (@list);
		#$header{'DATASECTION'}{$cage}{$output[$o]}{$time}{'line'}=$lineN;
	      }
	  }
      $count++;
      }
    close (F);
    return %header;
  }

sub display_header 
{
  my ($hr)=shift;
  my $f= shift;
  my %h = %$hr;

  foreach my $k1 (keys (%{$h{$f}}))
    {
      if ($k1 eq "INTERVALS"){next;}
      foreach my $k2 (keys (%{$h{$f}{$k1}} ))
	{
	  foreach my $k3 (keys (%{$h{$f}{$k1}{$k2}}))
	    {
	      print "#h $f $k1 $k2 $k3 $h{$f}{$k1}{$k2}{$k3}\n";
	    }
	}
    }
}
      
sub header2value
{
  my ($name,$hr, $file)=@_;
  my %h = %$hr;
  my ($k1,$k2,$k3, @fl);
  
  if ($file){@fl=($file);}
  else {@fl=keys(%h);}
	
  foreach my $f (@fl)
    {
      foreach $k1 (keys (%{$h{$f}}))
	{
	  if ($k1 eq "DATAVALUE"){next;}
	  foreach $k2 (keys (%{$h{$f}{$k1}} ))
	    {
	      foreach $k3 (keys (%{$h{$f}{$k1}{$k2}}))
		{
		  if ($k3=~/$name/)
		    {return $h{$f}{$k1}{$k2}{$k3};}
		}
	    }
	}
    }
}
sub set_channel 
  {
    
    my $hin=shift;
    my $file=shift;
   
    my $id=shift;
    my $index=shift;
    my $channel=shift;
    my @values=@_;
    my ($v);
    my %h=%$hin;
    $index++;
   
    $h{$file}{$id}{$channel}{Index}=$index;
   
    foreach $v (@values)
      {
	my ($l, $r);
	($l, $r)=split (/=/, $v);
	$h{$file}{$id}{$channel}{$l}=$r;
      }
    return $index;
  }


sub set_handle
    {
      my $f=shift;
      my $exp= shift;
      my $F=new FileHandle;
      my $LineN
      open($F, $f);
      while (<$F>)
	{
	  $LineN++;
	  if (/$exp/){return ($F,$LineN);}
	}
      close ($F);
      return 0;
    }
