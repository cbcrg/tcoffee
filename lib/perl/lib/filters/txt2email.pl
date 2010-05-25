#!/usr/bin/env perl

while (<>)
	{
	/(\S+\@\S+)/;
	$m=$1;
	if ($m=~/[^\w,;.-@\:]/ || $m=~/[\>\<]/ || $m=~/\d\d\d/ || $m=~/AA/ || $m=~/hotmail/ || $m=~/yahoo/ || $m=~/wanado/ || $m=~/gmail/)
	  {next;}
	$m=~/(\S+)@/;
	$name=$1;
	$n=lc $name;
	
	if( !$mail{$n})
	  {
	    $mail{$n}=1;
	    print "\n$m; ";
	  }
      }
