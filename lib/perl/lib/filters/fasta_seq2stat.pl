#!/usr/bin/env perl

use Env qw(HOST);
use Env qw(HOME);
use Env qw(USER);
#Specific Files and Default Parameters
# VersionTag  1.00

while (<>)
  {
    if ( /^>(\S+)/)
      {;}
    else 
      {
	@list=/(\S)/g;

	foreach $s (@list)
	  {

	    $freq{$s}++;
	    $tot++;
	  }
      }
  }

@symbol_list=keys %freq;
foreach $s (@symbol_list)
  {
    printf "%s\t%.2f\n", $s, $freq{$s}/$tot;
  }

      
