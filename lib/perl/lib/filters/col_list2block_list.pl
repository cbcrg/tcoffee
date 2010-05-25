#!/usr/bin/env perl

open (F, $ARGV[0]);
while (<F>)
	{
	/(.*)/;
	$n1=$1;
	$n2=$n1+1;
	if ($n1){print "cons $n1 $n2\n";}
	}
	 
