#!/usr/bin/env perl
#This script reads an aln in fasta format and modifies names so thathere are no duplicates
#if two sequences have the same name, the second one is renamed name_1 and so on
use Env qw(HOST);
use Env qw(HOME);
use Env qw(USER);
#Specific Files and Default Parameters
while (<>)
  {
    if ( /^>(\S+)/)
      {
	if ($list{$1})
	  {
	    print ">$1_$list{$1}\n";
	    $list{$1}++;
	  }
	else
	  {
	    print $_;
	    $list{$1}=1;
	  }
      }
    else
      {
	print $_;
      }
  }
      
