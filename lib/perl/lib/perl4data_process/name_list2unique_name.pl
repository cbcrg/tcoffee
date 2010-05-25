#!/usr/bin/env perl

if ( $#ARGV==100)
  {
    print " Retuns a list of names that occur at least one time in an unordered list";
    print " <file>";
    die;
  }
    
while ( <>)
  {
    @l=($_=~/(\b[\S]+\b)/g);
    foreach $name (@l)
      {
	if ( $list{$name}==1){;}
	else
	 {
	   print "$name\n";
	   $list{$name}=1;
	 }
     }
  }
    
