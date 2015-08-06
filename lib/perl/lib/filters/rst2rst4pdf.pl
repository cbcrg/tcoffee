#!/usr/bin/env perl
use strict;
if (!$ARGV[0])
  {
    print "\nusage: rst2rst4pdf <input file> <max line length> > stdou\n";
    die;
  }
my $temp="temp$$";
if (-e $temp){unlink $temp;}

my $l;
open (F, $ARGV[0]);
while(<F>){$l.=$_;}close F;
$l=~s/\\\n //g;

open (T, ">$temp");
print T "$l";
close (T);

open (F,$temp);
if (!$ARGV[1])
  {
    while (<F>){print "$_";}
  }
else
  {
    my $current;
    my $in;
    
    while (<F>)
      {
	my $l=$_;
	if ($l=~/^::/)
	  {
	    $current=$l;
	    $in=1;
	  }
	elsif ($in==1 && $l eq "\n")
	  {
	    $current.=$l;
	    $in=2;
	  }
	elsif ($in==2 && $l eq "\n")
	  {
	    $current.=$l;
	    wrap_box ($current, $ARGV[1]);
	    $in=0;
	  }
	elsif ($in==2)
	  {
	    $current.=$l;
	  }
	else
	  {
	    print "$l";
	  }
      }
  }
close (F);

sub wrap_box
  {
    my ($string,$size)=@_;
    my $cc=0;
    my $k=0;
    my $ret=0;
    
    
    my @cl=split (//, $string);

    
    for (my $a=0; $a<=$#cl; $a++)
      {
	
	my $char=$cl[$a];
	if ($char eq "\n"){$ret++;}
	else {$ret=0;}
	
	if ($char eq "\n")
	  {
	    print "$char";
	    $cc=0;
	  }
	elsif ($cc==$size)
	  {
	    print "\\\n $char";
	    $cc=1;
	  }
	else
	  {
	    print "$char";
	    $cc++;
	  }	
      }
    return;
  }







