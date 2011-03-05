#!/usr/bin/env perl
use HTTP::Date;
use strict;
use FileHandle;
my $p;

my $F=new FileHandle;
open ($F, shift(@ARGV));

while (<$F>)
  {
    my $l=$_;
    chomp($l);
    if ( $l=~/#d/)
      {
	my $d={};
	my @v=($l=~/([^;]+)/g);
	shift @v;
	my $exp=shift (@v);
	my $record=shift(@v);
	for (my $a=0; $a<=$#v; $a+=2)
	  {
	    $d->{$v[$a]}=$v[$a+1];
	  }
	foreach my $f (@ARGV)
	  {
	    printf "%s %15s ", $f, $d->{$f};
	  }
	print "\n";
      }
  }	
close ($F);
