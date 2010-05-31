#!/usr/bin/env perl
use Env qw(HOST);
use Env qw(HOME);
use Env qw(USER);


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
      

