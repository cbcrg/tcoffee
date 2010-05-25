#!/usr/bin/env perl

if ( $#ARGV!=2)
  {
    print " Retuns all the sequences at list <Thrshold> identical to <seq>";
    print " <sim file> <seq of interest> <threshold>\n";
    print " sim file:\n<TOP|BOT> <seqname1> <seqname2> <percentid>\n";
    die;
  }

$tmp_file="tmp_file$$.tmp";

`grep $ARGV[1] $ARGV[0] >$tmp_file`;
  
  
open (F, $tmp_file);

while ( <F>)
  {
    @l=($_=~/(\b[\S]+\b)/g);
    if (($l[0] eq "TOP" || $l[0] eq "BOT") && $l[1] eq $ARGV[1]){$seq{$l[2]}=$l[3];}    
  }

sub byid {$seq{$b} <=> $seq{$a}}

print ">$ARGV[1] (100, $ARGV[1])\n\n";

for  $s (sort byid keys %seq)
      {
	if ($seq{$s}>$ARGV[2]){ print ">$s ($seq{$s}, $ARGV[1])\n\n";}
      }

unlink $tmp_file;
